#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Tyagi (1968) nearest
neigbor heuristic. It extends the basic procedure in nearest_neighbor.py with
an intra route improvement scheme.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has moderate dependencies: a TSP solver (the built in local
search solver can be used), and numpy and scipy for reading and preparing the
problem instance.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from sys import stderr

import numpy as np
from math import ceil
from logging import log, DEBUG

from classic_heuristics.nearest_neighbor import nearest_neighbor_init

try:
    # Intead of the outdated TSP heuristic described in Tyagi (1968) just use LKH.
    from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as solve_tsp

from util import sol2routes, routes2sol, totald, objf

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

TYAGI_SEED_METHOD = "farthest"

def _calculate_penalty(sol_or_route, D, d, C):
    penalty = 0
    routes = sol2routes(sol_or_route) 
    for r in routes:
        penalty += C-totald(r,d)*objf(r,D)
    return penalty 

def _calculate_route_demand_statistics(sol, d, K=None):
    if K is None:
        K = sol.count(0)-1
    route_demands = np.zeros((K))
    route_demand = 0.0
    route_index = 0
    for n in sol:
        if n==0 and route_demand!=0.0:
            route_demands[route_index] = route_demand
            route_demand = 0.0
            route_index+=1
        else:
            route_demand+=d[n]
    # bigger is better, therefore the std is negative (smaller std better)
    #  in case there is a tie, the one with bigger average wins (that is, 
    #  the loads are both balanced, but the one with larger loads wins).
    return -np.std(route_demands), np.average(route_demands)

def _push_if_popped_to_peek_queue(pq, n, queue_key):
    """ Ensures the element n will be popped again, that is inserts it back 
        to the PeekQueue (see nn.py). This does nothing if n is not yet popped.
        """
    pushed = False
    if pq.posleft>=0:
        if pq.l[pq.posleft][1]<=queue_key:
            # push back (replace previous)
            pq.l[pq.posleft] = (n, queue_key) 
            pq.posleft -= 1
            pushed = True
            #if __debug__: log(DEBUG, "REMOVEME: PUSH TO LEFT")
    if not pushed and pq.posright<-1:
        if pq.l[pq.posright][1]>=queue_key:
            # push back (replace previous)
            pq.l[pq.posright] = (n, queue_key) 
            pq.posright += 1
            #if __debug__: log(DEBUG, "REMOVEME: PUSH TO RIGHT")
    
    return
    
def route_end_interchange(route, route_d, route_l, D, d, C, L, 
                     served, node_nearest_neighbors,
                     try_to_interchange_all=True):
    """ This is the proto-route (group) optimization scheme proposed in 
    (Tyagi 1968). It concerns trying to replace the first or last node
    in the proto-route with the node (next_nearer) that would've been added
    next if the constraint had not been violated. 
    
    if try_to_interchange_all is set to be True (default False), a variant 
    where any of the proto-route nodes can be replaced with the next_nearer.
    """
    
    ## Get the node that is "the next nearer point to the last included point"
    #  and assyme Tyagi meant of the unrouted nodes.
    
    first_added = route[1]
    last_added = route[-2]
    next_nearer = None
    while True:
        next_nearer = node_nearest_neighbors[last_added].peekleft()[0]
        if served[next_nearer]:
            node_nearest_neighbors[last_added].popleft()
        else:
            break

    ## Try to make the interchanges, and keep the best improving.
    
    if try_to_interchange_all:
        # as Tyagi 1968 seems to have done it if inferred from the example
        interchange_candidates = zip( range(1,len(route)-1), route[1:-1] )
    else:
        # as descibed in Tyagi 1968
        # try to iterchange the next_nearer with last / first
        interchange_candidates = [(-2,last_added), (1, first_added)]
    
    to_add_back_idx = None
    best_d_delta = 0.0
    best_l_delta = 0.0
    for to_remove_idx, to_remove_node in interchange_candidates:
        if C:
            d_delta = -d[to_remove_node]+d[next_nearer]
        if L:
            insert_after = route[to_remove_idx-1]
            insert_before = route[to_remove_idx+1]
            
            l_delta = -D[insert_after, to_remove_node]\
                      -D[to_remove_node,insert_before]\
                      +D[insert_after, next_nearer] \
                      +D[next_nearer,insert_before]

        is_L_feasible = not L or route_l+l_delta-S_EPS<=L
        is_C_feasible = not C or route_d+d_delta-C_EPS<=C
        if is_L_feasible and is_C_feasible and (
              (C and d_delta>best_d_delta) or
              (not C and L and l_delta<best_l_delta)):
            # insert last_added back!
            to_add_back_idx = to_remove_idx
    
    ## Found an improving interchange, now push the node to be removed back to 
    # the unrouted nodes. 
    
    if to_add_back_idx is not None:
        to_add_back = route[to_add_back_idx]
        served[next_nearer] = True
        for i in range(len(D)):
            if (i==0 or not served[i]) and i!=to_add_back:
                nnli = node_nearest_neighbors[i]
                _push_if_popped_to_peek_queue(nnli,to_add_back, D[i,to_add_back])                                             
                
        served[to_add_back] = False
        
        new_route = route[:to_add_back_idx]+[next_nearer]+route[to_add_back_idx+1:]
        if __debug__:
            log(DEBUG-1, "Interchange n%d to n%d"%(to_add_back, next_nearer))
            log(DEBUG, "Improved route : %s"%str(new_route))
        return new_route
    
    ## No changes, keep the original route 
    
    return route
    

def tyagi_init(D, d, C, L=None,
               only_large_and_close_customer_single_routes = False,
               try_interchange_with_all_group_customers = True,
               select_grouping = "min_penalty" ):
    """ This is an implementation of the Tyagi (1968) nearest neighbor 
    heurstic. In the first phase the method distributes customers into groups
    using bulding a nearest neighbor chains leaving from the depot. The group
    is deemed complete when a constraint is violated. Then an interchange 
    heuristic tries to remove one node from the group so that the one left out 
    can be fitted. This interchange is done if the demand of the route 
    increases due to the swap. If a group with only one customer is generated, 
    all of the candidates that are close to the depot and have large demand are 
    considered to be served with a route with a single customer. Furthermore 
    the entire procedure is repeated for all single customer route candidates 
    and the one that has "maximum use of the capacity" is used (see option 
    select_grouping). Then a TSP procedure is used to route the customer
    groups. This this implementation uses 3-Opt local search to do the routing 
    instead of the greedy TSP heuristic described in Tyagi (1968).
    
    The procedure above is used when the function is called with parameters:
     * only_large_and_close_customer_single_routes = True,
     * try_interchange_with_all_group_customers = False,
     * select_grouping = "max_demand"
    The "max_demand" selects the grouping where the non-single routes carry
    most capacity.
      
    However, these settings do not greate the grouping nor solution presented
    as an illustrative example in Tyagi (1968). It seems he uses different
    undocumented variant of his heuristic to solve the 12 customer problem
    from Dantzig and Ramser (1958). Therefore, following changes and options 
    were added, and can be used, to replicate the result:
     * only_large_and_close_customer_single_routes = False, which allows any 
        customer be served with a single customer route. 
     * try_interchange_with_all_group_customers = True, which tries to do the
        interchange with all of the nodes in a group (proto-route), instead of
        trying just the first and last customers.
     * select_grouping = "balanced", which selects the most balanced
        alternative among the grouping candidates with single customer route.
        Here, only the non-single customer routes are considered.
    Alternatively, another option for selecting grouping (which more closely
    follows the idea of Tyagi) can be used:
     * select_grouping = "min_penalty", which selects the option where 
        the routes with lot of free capacity are shortest (individual edges 
        or the used vehicle capacity when traversing them are not considered).
        This is calculated on the optimized routes.
        
    Note that if C is not set the select_grouping and its functionality is not
     used.
    
    With these options the grouping of customers given in illustrative
    example is replicated and as is the solution quality. Unfortunately it 
    is impossible to say if the modifications are the ones Tyagi (1968) used.
        
    Tyagi, M.S. (1968), "A Practical Method for Truck Dispatching Problem",
        J. Operations Research Society of Japan, 10, 76-92.
    Dantzig, G. B., & Ramser, J. H. (1959), "The truck dispatching problem",
        Management science, 6(1), 80-91.
    Helsgaun, K. (2006), "An Effective Implementation of K-opt Moves for the 
      Lin-Kernighan TSP Heuristic." DATALOGISKE SKRIFTER, No. 109, 2006. 
      Roskilde University.
    Helsgaun, K. (2009), "General k-opt submoves for the Lin-Kernighan TSP 
      heuristic." Mathematical Programming Computation, 1(2), 119-163.    
    """
    
    
    ## 0. Determine the number of vehicles
    may_have_single_node_route = False
    if C:
        tot_d = sum(d)
        K = ceil(tot_d/float(C))
        may_have_single_node_route = min(d[1:]) > tot_d-C*(K-1)
    
    ## 1. "Grouping the delivery points" (nodes) to routes with nearest nbour.
    interchange_variant = lambda r, r_d, r_l, D, d, C, L, svd, nnn : \
        route_end_interchange(r, r_d, r_l, D, d, C, L, svd, nnn, 
                              try_interchange_with_all_group_customers)
    solution = \
        nearest_neighbor_init(D, d, C, L, 
                              initialize_routes_with=TYAGI_SEED_METHOD,
                              emerging_route_count=1, add_only_to_end=True,
                              route_improvement_callback=interchange_variant)
    
    routes = sol2routes(solution)

    if not C:
        best_C_utilization = None
    elif select_grouping == "balanced" :
        # the best will be one that will have balanced routes
        best_C_utilization = \
            _calculate_route_demand_statistics(solution, d, K=len(routes))
    elif select_grouping == "max_demand" :
        # the best will be the one that is best to top up the vehicle capacity
        best_C_utilization = sum(d[n] for n in solution)/float(len(routes)*C)
    elif select_grouping == "min_penalty" :
        # the best will be the one with least distance traveled with unused capacity
        best_C_utilization = _calculate_penalty(solution, D, d, C)
    else:
        raise ValueError('"Invalid select_grouping parameter "%s"'%select_grouping)
        
    if __debug__:
        if C:
            log(DEBUG-1, "Initial Grouping = %s\n"%(solution))
          
    
    ## Tyagi heuristic has a special case for the single node routes
    try:
        single_nodes = [r[1] for r in routes if len(r)==3]
        if C and (may_have_single_node_route or len(single_nodes)>0):
            N = len(D)
            best_K = solution.count(0)-1
            
            if only_large_and_close_customer_single_routes:
                # According to Tyagi1967, only the large demand nodes close to the
                #  depot are considered ...
                max_dist_from_depot = max(D[0,:])
                candidates = [n for n in range(1,N) \
                          if (d[n]>0.5*C and D[0,n]<0.5*max_dist_from_depot)]
            else:
                # ... however, using this criteria does not allow replicating the
                #  results (the illustrative example) of the paper. Therefore,
                #  there is an option to try with complete set of nodes.
                candidates = list(range(1,len(D)))
                
            for single_route_candidate in candidates:
                if __debug__:
                    log(DEBUG-2, "Try to find grouping with n%d forbidden"%(single_route_candidate))
                cndt_sol = nearest_neighbor_init(D, d, C, L, 
                               initialize_routes_with=TYAGI_SEED_METHOD,
                               emerging_route_count=1,
                               forbidden_nodes=[single_route_candidate],
                               route_improvement_callback=interchange_variant)
                               
                cndt_K = cndt_sol.count(0)-1+1 #the single_route_candidate_node
                if cndt_K>best_K:
                    if __debug__:
                        log(DEBUG-1, "No valid grouping found\n")
                
                    continue
                
                if select_grouping == "balanced" :
                    cndt_C_utilization = \
                        _calculate_route_demand_statistics(cndt_sol, d, K=cndt_K-1)
                elif select_grouping == "max_demand" :
                    cndt_C_utilization = sum(d[n] for n in cndt_sol)/float(cndt_K*C)
                elif select_grouping == "min_penalty" :
                    cndt_C_utilization = _calculate_penalty(cndt_sol, D, d, C)
    
                if __debug__:
                    log(DEBUG-1, "Grouping = %s (utilization %s)\n"%(
                        cndt_sol+[single_route_candidate,0],str(cndt_C_utilization)))
                 
                # select the one with maximum possible use of carrier capacity
                if cndt_C_utilization>best_C_utilization:
                    best_C_utilization = cndt_C_utilization
                    solution = cndt_sol+[single_route_candidate, 0]
                    routes = sol2routes(solution)
            
        ## 2. "Finding the optimal tours"
        #instead of the heuristic of Tyagi, just solve it with a TSP solver
        if __debug__: 
            log(DEBUG-1, "Post-optimize solution %s (%.2f)"%(solution, objf(solution,D)))
        for i in range(len(routes)):
            routes[i],route_f = solve_tsp(D, routes[i][:-1])
            
            if __debug__:
                log(DEBUG-2, "Got TSP solution %s (%.2f)"%(str(routes[i]), route_f))
    except KeyboardInterrupt: #or SIGINT
        raise KeyboardInterrupt(solution)
    
    return routes2sol(routes)
                               

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ty_algorithm():
    algo_name = "Ty68-NN"
    algo_desc = "Tyagi (1968) Nearest Neighbor construction heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if minimize_K:
            raise NotImplementedError("Nearest neighbor algorithm does not support minimizing the number of vehicles")
        return tyagi_init(D, d, C, L )
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ty_algorithm())
