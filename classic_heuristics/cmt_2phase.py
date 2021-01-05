#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Christofides, Mingozzi
and Toth (1979) two phase heuristic. Also a deterministic variant is included.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: orderedset, a TSP solver, 
and numpy and scipy for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import numpy as np
from collections import deque, namedtuple 
from logging import log, DEBUG

# for the original stochastic version (needs to be enabled separately)
from random import shuffle
from sys import stderr

# One of the few non standard-lib additions, used in algorithm's second 
#  phase to keep track of the non-routed nodes. Install it e.g. with:
#
# $ pip install orderedset
#
from orderedset import OrderedSet

try:
    ## For tiny instances you might want to get the optimal solution
    #from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp
    ## This is the default TSP solver. Used mainly if L constraint is set.
    from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp
    ## For largest instances we have reserved an option to rely rely on a CUSTOM acotsp 
    #from tsp_solvers.tsp_solver_acotsp import solve_tsp_acotsp as solve_tsp
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as solve_tsp

from util import is_better_sol, routes2sol, without_empty_routes, objf
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS


__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


## Helper functions for different initialization schemes ##

def _first_seed(D, d, unrouted):
    return unrouted[0]

def _farthest_seed(D, d, unrouted):
    largest_D_0l_unrouted_idx = int(np.argmax( D[[0], unrouted] ))
    return unrouted[largest_D_0l_unrouted_idx]

def _closest_seed(D, d, unrouted):
    largest_D_0l_unrouted_idx = int(np.argmin( D[[0], unrouted] ))
    return unrouted[largest_D_0l_unrouted_idx]

def _biggest_seed(D, d, unrouted):
    npd = np.array(d)
    largest_d_unrouted_idx = int(np.argmax(npd[unrouted]))
    return unrouted[largest_d_unrouted_idx]

def _phase_one(lambda_multiplier, D,d,C,L, seed_f, rr):
    """ This imlements the fist phase of the algorithm. Sequentally
     add nodes to an emerging node. Different seed node selection 
     functions (above) can be used. """
    
    route_seeds = []
    N = len(D)
    sol = []
    total_cost = 0.0
    
    customer_nodes = list(range(1,N))
    if rr is not None:
        shuffle(customer_nodes)
        rr-=1
    unrouted = OrderedSet(customer_nodes)
    
    if __debug__:
        log(DEBUG, "## Sequential route bulding phase ##")
    
    route_idx = 0
    
    try:
        while unrouted:
            ## Step 1: choose unrouted i_k to act as a route seed point
            
            route_seed_k = seed_f(D, d, unrouted)
            unrouted.remove(route_seed_k)
            
            route_seeds.append(route_seed_k)
            route_demand = d[route_seed_k] if C else 0
            route = [0, route_seed_k]
            route_cost = D[0,route_seed_k]+D[route_seed_k,0]
            route_l_updated = True
            route_idx += 1
            
            if __debug__:
                log(DEBUG, "Initialize route #%d with n%d"%
                           (route_idx, route_seed_k))
            
            ## Step 2: Compute savings
            
            s_vals = (D[[0],unrouted]+lambda_multiplier*
                      D[unrouted,[route_seed_k]]).tolist()
            savings = list(zip(s_vals, unrouted))
            savings.sort()
            for best_saving, i in savings:
                ## Step 3: insert until feasibility is broken
            
                if __debug__:
                    log(DEBUG, "Check feasibility of inserting n%d with savings=%.2f"%
                               (i,best_saving))
        
                #TODO: To improve performance keep track of minimal unrouted d and
                # break the savings loop if we see that there are no feasible
                # insertions to be made.
                
                if C and route_demand+d[i]-C_EPS>C:
                    # capacity constraint violated, route complete
                    if __debug__:
                        log(DEBUG,"Insertion would break C constraint, skip.")
                    continue
    
                #TODO: To improve performance, keep track of minimal possible 
                # cost increase to the route. This would involve updating e.g.
                # Held-Karp (1970) lower bound for the route.
                
                # Use upper bound estimate to save some computations. We know 
                #  that the maximum route cost constraint cannot be violated.
                UB_route_cost = route_cost-D[route[-1],0]+D[route[-1],i]+D[i,0]            
                if L and UB_route_cost-S_EPS>L:
                    new_route, new_route_cost = solve_tsp(D, route+[i])
                    if __debug__:
                        log(DEBUG-1,"Got TSP solution %s (%.2f)"%
                            (list(new_route),new_route_cost))
             
                    if new_route_cost-S_EPS>L:
                        if __debug__:
                            log(DEBUG, "Insertion would break L constraint, skip.")
                        continue
                    route_cost = new_route_cost
                    route_l_updated = True
                else:
                    route_l_updated = False
                    route_cost = UB_route_cost
                    new_route = route+[i,0]
                    
                # accept including node i
                route=new_route[:-1]
                if C: route_demand+=d[i]
                unrouted.remove(i)
    
                if __debug__:
                    log(DEBUG, "Inserted n%d to create a route %s (%.2f)."%
                        (i, route, route_cost))
        
            # if L is not set, optimize TSP after the route is full
            if not route_l_updated:
                new_route, route_cost = solve_tsp(D, route)
                route = new_route[:-1]
                if __debug__:
                    if not L:
                        log(DEBUG-1,"Got TSP solution %s (%.2f)" %
                            (str(route+[0]), route_cost))
            if __debug__:
                log(DEBUG, "Route %s (%.2f) complete.\n"%
                          (str(route+[0]),route_cost))
     
            total_cost += route_cost
            sol+=route
    except KeyboardInterrupt:  #or SIGINT
        interrupted_sol = sol+routes2sol([n] for n in unrouted if n not in sol)
        raise KeyboardInterrupt(interrupted_sol)
        
    sol+=[0]
    
    if __debug__:
        log(DEBUG, "Phase 1 solution %s (%.2f) complete.\n"%(str(sol),total_cost))
        log(DEBUG-1, "Pass on route seeds %s to Phase 2.\n"%str(route_seeds))            
    
    return route_seeds, sol, total_cost, rr

RouteState = namedtuple('RouteState', 'route demand cost updated')
def _routestates2solution(routes, D):
    sol = []
    total_cost = 0.0
    for route,_,route_cost,route_l_updated in routes:
        # if L is not set, optimize TSP *in the end*
        if not route_l_updated:
            new_route, route_cost = solve_tsp(D, route)
            route = new_route[:-1]
            if __debug__:
                log(DEBUG-1, "Got TSP solution %s (%.2f)" % (str(new_route), route_cost))
        total_cost += route_cost
        sol+=route
    sol+=[0]
    return sol, total_cost
    
def _phase_two(mu_multiplier,route_seeds, D,d,C,L, rr,
               choose_most_associated_route = True,
               repeated_association_with_n_routes=1):
    ## Step 0: reuse seed nodes from phase 1 to act as route seed points
    N = len(D)
    K = len(route_seeds)

    customer_nodes = list(range(1,N))
    if rr is not None: 
        #->stochastic version, resolve the ties randomly
        shuffle(customer_nodes)
    unrouted_nodes = OrderedSet(customer_nodes)
    unrouted_nodes.difference_update( route_seeds )
    
    # routes are stored in dict with a key of route seed,
    #  the route, demand, cost, and if it is updated are all stored
    routes = [ RouteState( [0, rs],            #initial route
                            d[rs] if C else 0, #initial demand 
                            D[0,rs]+D[rs,0],   #initial cost
                            True) for rs in route_seeds ]
    
    insertion_infeasible = [[False]*N for rs in route_seeds ] 
            
    if __debug__:
        log(DEBUG, "## Parallel route building phase with seeds %s ##" %
                  str(list(route_seeds)))
    
    ## Step 1.1: vectorized calculation of eps
    # TODO: this also calculates eps for depot and seeds. Omitting those would 
    #  be possible and save few CPU cycles, but it would make indexing more
    #  complex and because accuracy>simplicity>speed, it is the way it is.
        
    eps = ( np.tile(D[[0],:],(K, 1)) 
            +mu_multiplier*D[:,route_seeds].transpose()
            -np.tile(D[0,route_seeds],(N,1)).transpose() )
    
    associate_to_nth_best_route = 1
    insertions_made = False
    route_seed_idxs = []
    insertions_made = False
    first_try = True
    
    try:
        while unrouted_nodes:                
            ## Main while loop bookkeeping
            if not route_seed_idxs:
                idxs = list(range(K))
                if rr is not None: #->stocastic version, construct routes in random order
                    shuffle(idxs)
                route_seed_idxs = deque(idxs)
                
                if not first_try:
                    # The CMT1979 exits when all routes have been tried
                    if repeated_association_with_n_routes is None:
                        break            
                
                    if not insertions_made:
                        associate_to_nth_best_route+=1 # for the next round
                    if associate_to_nth_best_route>\
                       repeated_association_with_n_routes:
                        break
                
                first_try = False
                insertions_made = False
                                
            ## Step 2.1: Choose a (any) route to add customers to.        
        
            # some nodes may have been routed, update these
            eps_unrouted = eps[:,unrouted_nodes]
            
            ## Step 1.2: Associate each node to a route
            # note: the assignments cannot be calculated beforehand, as we do 
            #  not know which customers will be (were?) "left over" in the 
            #  previous route building steps 3. 
            
            if associate_to_nth_best_route==1 or len(route_seed_idxs)==1:
                r_stars_unrouted = np.argmin(eps_unrouted[route_seed_idxs,:], axis=0)    
            else:
                ## note: an extension for the deterministic variant,
                # get smallest AND 2. smallest at the same time using argpartition
                #top = np.argsort(eps, axis=0)[:associate_with_n_routes, :]
                #r_stars = [ top[i,:] for i in range(associate_with_n_routes) ]
                if len(route_seed_idxs)<associate_to_nth_best_route:
                    route_seed_idxs = []
                    continue
                    
                #take_nth = min(len(route_seed_idxs), associate_to_nth_best_route)-1
                take_nth = associate_to_nth_best_route-1
                reorder_rows_per_col_idxs = np.argpartition(eps_unrouted[route_seed_idxs,:],
                                              take_nth, axis=0)
                nth_best, unrouted_node_order = np.where(
                        reorder_rows_per_col_idxs==take_nth)
                r_stars_unrouted = nth_best[ np.argsort(unrouted_node_order) ]
            
            if choose_most_associated_route:
                unique, counts = np.unique(r_stars_unrouted, return_counts=True)
                seed_idx_idx = unique[np.argmax(counts)]
                route_seed_idx = route_seed_idxs[seed_idx_idx]
                route_seed_idxs.remove(route_seed_idx)
                associated_cols = list(np.where(r_stars_unrouted==seed_idx_idx)[0])
            else:
                route_seed_idx = route_seed_idxs.popleft()
                associated_cols = list(np.where(r_stars_unrouted==0)[0])
                
            route,route_demand,route_cost,route_l_updated = routes[route_seed_idx]
                
            ## Step 2.2: Vectorized calculation of sigma score for the customers
            #             associated to the chosen route.
            eps_bar = eps_unrouted[route_seed_idx,associated_cols]
            
            # NOTE: CMT 1979 does not specify what happens if S is empty, we assume
            #  we need (and can) omit the calculation of eps_prime in this case.
            
            brdcast_rs_idxs = [[rsi] for rsi in route_seed_idxs]
            if route_seed_idxs:
                eps_prime = np.min(eps_unrouted[brdcast_rs_idxs, associated_cols],
                                   axis=0)
                sigmas = eps_prime-eps_bar
            else:
                # last route, try to add rest of the nodes
                eps_prime = None
                sigmas = -eps_bar
                
            col_to_node = [unrouted_nodes[int(c)] for c in associated_cols]
            sigma_ls = list(zip(sigmas.tolist(), col_to_node))
            sigma_ls.sort(reverse=True)
            
            if __debug__:
                log(DEBUG, "Assigning associated nodes %s to a route %s (seed n%d)"%
                          (str(col_to_node), str(route+[0]),route_seeds[route_seed_idx]))
            
            ## Step 3: insert feasible customers from the biggest sigma first
            for sigma, l_star in sigma_ls:
                if __debug__:
                    log(DEBUG-1, "Check feasibility of inserting "+\
                        "n%d with sigma=%.2f"%(l_star,sigma))
                
                if C and route_demand+d[l_star]-C_EPS>C:
                    if __debug__:
                        log(DEBUG-1, "Insertion would break C constraint.")
                    continue
                
                # use cached L feasibility check
                if L and insertion_infeasible[route_seed_idx][l_star]:
                    continue
                
                # Do not run TSP algorithm after every insertion, instead calculate
                #  a simple a upper bound for the route_cost and use that.
                UB_route_cost = (route_cost-
                                   D[route[-1],0]+
                                   D[route[-1],l_star]+D[l_star,0])
                if L and UB_route_cost-S_EPS>L:
                       
                    # check the real TSP cost
                    new_route, new_route_cost = solve_tsp(D, route+[l_star])
    
                    if __debug__:
                        log(DEBUG-1, "Got TSP solution %s (%.2f)" %
                            (str(new_route), new_route_cost, ))
                        
                    if new_route_cost-S_EPS>L:
                        if __debug__:
                            log(DEBUG-1,"DEBUG: Insertion would break L constraint.")
                        insertion_infeasible[route_seed_idx][l_star] = True
                        continue
                    
                    route_cost = new_route_cost
                    route=new_route[:-1]
                    route_l_updated = True
                else:
                    route_l_updated = False
                    route_cost = UB_route_cost
                    route = route+[l_star]
                
                if C: route_demand+=d[l_star]
                unrouted_nodes.remove(l_star)
                insertions_made = True
                
                if __debug__:
                    log(DEBUG, "Inserted n%d to create a route %s."%(l_star, route))
        
        
            # All feasible insertions of the associated customers is done, record
            #  the modified route.
            if insertions_made:
                routes[route_seed_idx] =  RouteState( route,       #updated route
                                                  route_demand,    #updated demand
                                                  route_cost,      #updated cost 
                                                  route_l_updated) #cost state
    except KeyboardInterrupt: #or SIGINT
        rs_sol, _ = _routestates2solution(routes, D)
        interrupted_sol = rs_sol[:-1]+routes2sol([n] for n in unrouted_nodes
                                                if n not in rs_sol)
        raise KeyboardInterrupt(interrupted_sol)
    
    ## Step 4: Redo step 1 or construct the solution and exit 
    if len(unrouted_nodes)>0:
        if __debug__:
            log(DEBUG, "Phase 2 failed to create feasbile solution with %d routes."%K)
            log(DEBUG-1, "Nodes %s remain unrouted."%str(list(unrouted_nodes)))
        return 0, None, None, rr
    else:
        sol, total_cost = _routestates2solution(routes, D)
        if __debug__:
            log(DEBUG, "Phase 2 solution %s (%.2f) complete."%(str(sol),total_cost))        
        return K, sol, total_cost, rr
    

def cmt_2phase_init(D, d, C, L=None, minimize_K=False,
                    lambda_multiplier=2.0, mu_multiplier=1.0,
                    phase1_seed_selection_method = "farthest",
                    phase2_choose_most_associated_route = True,
                    phase2_repeated_association_with_n_routes = 1,
                    number_of_randomized_retries = None):
    
    """ Implementation of the Christofides, Mingozzi & Toth (1979) two phase
    heuristic. In the first phase a customer is selected to act as a seed node 
    and initialize a route. Then, a savings criteria parametrized with
    lambda_multiplier is used to determine which customers to insert.
    Insertions are done until a constraint is violated and then a new seed is
    selected and the insertions continue. This is repeated until no unrouted
    customers remain or we run out of route seeds. Finally, the routes are made
    r-optimal with 3-opt. 
    
    The seed customers are carried over the the second phase of the algorithm.
    Here, each customer is associated to a seed customer based on a second 
    savings criteria parametrized with mu_multiplier. Also the next closest 
    seed customer has an effect to the score used when associating the nodes.
    Then, a route is built around each seed customer with the nodes associated
    to that route taking care not to violate feasibility of the route. Finally,
    if a feasible solution was generated, the routes from the second phase 
    are made r-optimal with 3-opt. 
    
    A better of the solutions from the first and second phases is selected 
    and returned.
    
    Note that the default parameters are for a deterministic variant of the 
    stochastic algorithm described in (Christofides et al 1979).
    
    Basic parameters:
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route cost/length/duration.
    
    Objective parameter:  
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes and the current best is always replaced
       with a solution with smaller K. If set to False (default) the algorithm 
       optimizes only for the mimimum solution/routing cost. 
    
    Route shape parameters:
    * lambda_multiplier   specifies how closely the customer is associated to 
                           the emerging route seed customer in the first phase.
    * mu_multiplier       specifies how closely the customer is associated to
                           route seed customers in the second phase.

    The implementation includes some improvements to the CMT (1979) algorithm
    to improve the chance of second phase producing feasible solutions:
        
    * phase1_seed_selection_method
                          instead of selecting a seed customer for emerging 
                           route at random in the first phase, select the 
                           "farthest" or "closest" to the depot or the one with 
                           the "biggest" demand. Can also be "first", which 
                           will be random if randomized_retries is set.
                           
    * phase2_choose_most_associated_route   
                           instead of building the routes in random order in 
                           phase 2, start from the route with most associated
                           customers. If set to False implements the original
                           behaviour of (CMT 1979).
    
    * phase2_repeated_association_with_n_routes           
                           if set to None, the original behaviour of (CMT 1979)
                           is used. That is, terminate phase 2 without an
                           feasible solution if the first route building pass
                           over the route seed customers leaves unrouted
                           customers when S=0. If this is set to 1, the
                           procedure is repeated until a) all customers are
                           routed or b) no feasible insertions can be made.
                           If this parameter is set to be >1, also insertion of
                           2. best alternatives to associate with a seed
                           customers are tried. Can also be "K". Then the
                           number of routes generated in the first phase is
                           used as the value of this parameter.
    * number_of_randomized_retries  
                           If None the algorithm is deterministic. If set to an 
                           integer value. The first phase can generate this 
                           many seed customer configurations to second phase 
                           in case second phase is unable to produce feasible
                           solutions.
    """
    
    if phase1_seed_selection_method=="first":
        seed_f = _first_seed
    elif phase1_seed_selection_method=="farthest":
        seed_f = _farthest_seed
    elif phase1_seed_selection_method=="closest":
        seed_f = _closest_seed
    elif phase1_seed_selection_method=="biggest":
        seed_f = _biggest_seed
    
    rr = number_of_randomized_retries 

    best_sol = None
    best_f = None
    best_K = None
    interrupted = False
    
    while (rr is None) or (rr>0):
        
        phase1_sol, phase1_f, phase1_K = None, float("inf"), float("inf")
        phase2_sol, phase2_f, phase2_K = None, float("inf"), float("inf")
        
        try:
            phase1_seeds, phase1_sol, phase1_f, rr = \
                _phase_one(lambda_multiplier,D,d,C,L, seed_f, rr)
            phase1_K = len(phase1_seeds)
            
            # extension to CMT, option to associate customers multiple times 
            #  (to other routes, starting from the route with minimal eps).
            associate_routes = phase2_repeated_association_with_n_routes
            if phase2_repeated_association_with_n_routes=="K":
                associate_routes = phase1_K
                
            phase2_K, phase2_sol, phase2_f, rr = \
                _phase_two(mu_multiplier,phase1_seeds,D,d,C,L, rr,
                    phase2_choose_most_associated_route, associate_routes)
        
        except KeyboardInterrupt as e: #or SIGINT
            # Phase 1 OR phase 2 was interrupted. 
            if len(e.args)>0 and type(e.args[0]) is list:
                if phase1_sol is None:
                    phase1_sol = without_empty_routes(e.args[0])
                    phase1_f = objf(phase1_sol)
                    phase1_K = phase1_sol.count(0)-1
                    
                elif phase2_sol is None:
                    phase2_sol = without_empty_routes(e.args[0])
                    phase2_f = objf(phase2_sol)
                    phase2_K = phase2_sol.count(0)-1
            interrupted = True
        
        # Pick the better out of the two
        p1_better_than_p2 = is_better_sol(phase2_f, phase2_K,
                                          phase1_f, phase1_K, minimize_K)
        p1_best_so_far    = is_better_sol(best_f, best_K,
                                          phase1_f, phase1_K, minimize_K)
        p2_best_so_far    = is_better_sol(best_f, best_K,
                                          phase2_f, phase2_K, minimize_K)
        if p1_better_than_p2 and p1_best_so_far:
                best_sol = phase1_sol
                best_f = phase1_f
                best_K = phase1_K
        if not p1_better_than_p2 and p2_best_so_far:
                best_sol = phase2_sol
                best_f = phase2_f
                best_K = phase2_K
        
        if interrupted:
            # pass on the current best solution
            raise KeyboardInterrupt(best_sol)
        
        # deterministic version, no retries
        # stochastic version terminates as soon as phase2 succeeds
        if (rr is None) or (phase2_sol is not None):
            break

        
    return best_sol


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_cmt2p_algorithm():
    algo_name = "CMT79-2P"
    algo_desc = "Christofides, Mingozzi & Toth (1979) two phase heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return cmt_2phase_init(D, d, C, L, minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_cmt2p_algorithm())
