#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementation of the Clarke and Wright (1964)
savings functions for parallel (as in multiple route) savings heuristic. The
provided algorithm is built to be extendable and many other savings variants
implemented in VeRyPy use the generic parallel savings heuristic implemented 
here.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has minimal dependencies: only numpy and scipy
are needed for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from logging import log, DEBUG
from util import TW, routes2sol, objf
from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def clarke_wright_savings_function(D, ctrs):
    N = len(D)
    n = N-1
    savings = [None]*int((n*n-n)/2)
    idx = 0
    for i in range(1,N):
        for j in range(i+1,N):
            s = D[i,0]+D[0,j]-D[i,j]
            savings[idx] = (s,-D[i,j],i,j)
            idx+=1
    savings.sort(reverse=True)
    return savings 

def parallel_savings_init(D, d, ctrs, minimize_K=False,
                          savings_callback=clarke_wright_savings_function):
    """
    Implementation of the basic savings algorithm / construction heuristic for
    capaciated vehicle routing problems with symmetric distances (see, e.g.
    Clarke-Wright (1964)). This is the parallel route version, aka. best
    feasible merge version, that builds all of the routes in the solution in
    parallel making always the best possible merge (according to varied savings
    criteria, see below).
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is a dict containing the constraints. Usually some combination of:
        - "C" is the capacity constraint limit for the identical vehicles.
        - "L" is the optional constraint for the maximum route length/duration.
        - "TWs" a list of timewindow pairs for each customer (0 being the depot)
    
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring a merge that would increase the total distance.
       WARNING: This only works when the solution from the savings algorithm is
       final. With postoptimimization this non-improving move might have still
       led to improved solution.
   
    * optional savings_callback is a function of the signature:
        sorted([(s_11,x_11,i_1,j_1)...(s_ij,x_ij,i,j)...(s_nn,x_nn,n,n) ]) =
            savings_callback(D)
      where the returned (sorted!) list contains savings (that is, how much 
       solution cost approximately improves if a route merge with an edge
       (i,j) is made). This should be calculated for each i \in {1..n},
       j \in {i+1..n}, where n is the number of customers. The x is a secondary
       sorting criterion but otherwise ignored by the savings heuristic.
      The default is to use the Clarke Wright savings criterion.
        
    See clarke_wright_savings.py, gaskell_savings.py, yellow_savings.py etc.
    to find specific savings variants. They all use this implementation to do 
    the basic savings procedure and they differ only by the savings
    calculation. There is also the sequental_savings.py, which builds the 
    routes one by one.
    
    Clarke, G. and Wright, J. (1964). Scheduling of vehicles from a central
     depot to a number of delivery points. Operations Research, 12, 568-81.
    """
    N = len(D)
    ignore_negative_savings = not minimize_K

    ## 1. make route for each customer
    routes = [[i] for i in range(1,N)]
    route_demands = d[1:] if ('C' in ctrs) else [0]*N
    if ('TWs' in ctrs):
        wait_times = [max(0, ctrs['TWs'][i][TW.OPEN]-D[0,i]) for i in range(0,N)]
        route_costs = [D[0,i]+wait_times[i]+D[i,0] for i in range(1,N)]
    elif ('L' in ctrs):
        route_costs = [D[0,i]+D[i,0] for i in range(1,N)]
    
    try:
        ## 2. compute initial savings 
        savings = savings_callback(D, ctrs)
        
        # zero based node indexing!
        endnode_to_route = [0]+list(range(0,N-1))
        
        ## 3. merge
        # Get potential merges best savings first (second element is secondary
        #  sorting criterion, and it it ignored)
        for best_saving, _, i, j in savings:
            if __debug__:
                log(DEBUG-1, "Popped savings s_{%d,%d}=%.2f" % (i,j,best_saving))
                
            if ignore_negative_savings:
                cw_saving = D[i,0]+D[0,j]-D[i,j]
                if cw_saving<0.0:
                    break
            
            # "left" and "right" are names by convention, no real orderings
            left_route = endnode_to_route[i]
            right_route = endnode_to_route[j]
            
            # the node is already an internal part of a longer segment
            if ((left_route is None) or
                (right_route is None) or
                (left_route==right_route)):
                continue
            
            if __debug__:
                log(DEBUG-1, "Route #%d : %s"%
                             (left_route, str(routes[left_route])))
                log(DEBUG-1, "Route #%d : %s"%
                             (right_route, str(routes[right_route])))
                
            # check capacity constraint validity
            if ('C' in ctrs):
                merged_demand = route_demands[left_route]+route_demands[right_route]
                if merged_demand-C_EPS > ctrs['C']:
                    if __debug__:
                        log(DEBUG-1, "Reject merge due to "+
                            "capacity constraint violation")
                    continue
            # check time window constraits
            if ('TWs' in ctrs):
                arrival_cost_at_j = route_costs[left_route]-D[0,i]+D[i,j]
                arrives_late = arrival_cost_at_j-S_EPS>ctrs['TWs'][j][TW.CLOSE]

                if arrives_late:
                    if __debug__:
                        log(DEBUG-1, "Reject merge due to it being outside "+
                            "of allowed time windows.")
                    continue
                
                j_tw_open = ctrs['TWs'][j][TW.OPEN]
                merged_wait_time_at_j = max(0, j_tw_open-arrival_cost_at_j)

                # Recalculate wait times and route duration (cost) and check for
                #  TW constraints for the right hand side of the merge.
                updated_wait_times = {j:merged_wait_time_at_j}
                arrival_cost_at_k = arrival_cost_at_j+merged_wait_time_at_j
                prev_k = j
                for k in routes[right_route]:
                    arrival_cost_at_k+=D[prev_k,k]
                    arrives_late = arrival_cost_at_k-S_EPS>ctrs['TWs'][k][TW.CLOSE]

                    if arrives_late:
                        if __debug__:
                            log(DEBUG-1, "Reject merge due to it causing being "
                                " late of the later time windows.")
                        continue
                    
                    k_tw_open = ctrs['TWs'][k][TW.OPEN]
                    updated_wait_times[k]=max(0, k_tw_open-arrival_cost_at_k)
                    arrival_cost_at_k+=updated_wait_times[k]
                    prev_k = k
                    

                merged_cost = arrival_cost_at_j+merged_wait_time_at_j\
                              +arrival_cost_at_k+D[prev_k, 0]
                
                
            # if there are route cost constraint, check its validity        
            if ('L' in ctrs):
                merged_cost = route_costs[left_route]-D[0,i]+\
                                route_costs[right_route]-D[0,j]+\
                                D[i,j]
                if merged_cost-S_EPS > ctrs['L']:
                    if __debug__:
                        log(DEBUG-1, "Reject merge due to "+
                            "maximum route length constraint violation")
                    continue

            
            
            # update bookkeeping only on the recieving (left) route
            #  right hand route is no more
            if ('C' in ctrs): route_demands[left_route] = merged_demand
            if ('TWs' in ctrs):
                for k, uwt in updated_wait_times.items():
                    wait_times[k] = uwt
                route_costs[left_route] = merged_cost
            elif ('L' in ctrs):
                route_costs[left_route] = merged_cost
            
            # TODO: VRPTW. Uhh. No. What to do here???
            # merging is done based on the joined endpoints, reverse the 
            #  merged routes as necessary
            if (not 'TWs' in ctrs):
                if routes[left_route][0]==i:
                    routes[left_route].reverse()
                if routes[right_route][-1]==j:
                    routes[right_route].reverse()
                    
            # the nodes that become midroute points cannot be merged
            if len(routes[left_route])>1:
                endnode_to_route[ routes[left_route][-1] ] = None
            if len(routes[right_route])>1:
                endnode_to_route[ routes[right_route][0] ] = None
            
            # all future references to right_route are to merged route
            endnode_to_route[ routes[right_route][-1] ] = left_route
            
            # merge with list concatenation
            routes[left_route].extend( routes[right_route] )
            routes[right_route] = None
            
            if __debug__:
                dbg_sol = routes2sol(routes)
                log(DEBUG-1, "Merged, resulting solution is %s (%.2f)"%
                             (str(dbg_sol), objf(dbg_sol,D)))
    except KeyboardInterrupt: # or SIGINT
        interrupted_sol = routes2sol(routes)
        raise KeyboardInterrupt(interrupted_sol)
        
    return routes2sol(routes)


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ps_algorithm():
    algo_name = "CW64-PS"
    algo_desc = "Clarke & Wright (1964) parallel savings algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        # Convert legacy call to new VRPTW supported init
        ctrs = {}
        if C: ctrs['C']=C
        if L: ctrs['L']=L
        return parallel_savings_init(D,d,ctrs,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ps_algorithm())
