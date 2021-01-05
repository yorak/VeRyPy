#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Stewart & Golden (1984)
3-opt* heuristic with Lagrangean relaxation.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: a TSP solver (the built-in one 
can be used) and the internal 3-opt* implementation. Also, numpy is needed for
by the algorithm for few convenience functions, as well as scipy for reading 
and preparing the problem instance.
"""
###############################################################################

import numpy as np
from logging import log, DEBUG
from random import shuffle
from functools import partial
from sys import stderr

try:
    # The original algorithm of Stewart & Golden (1984) was stochastic and used restarts
    #  with multiple random TSP solutions. For a deterministic variant, we use LKH to get
    #  a single high quality TSP solution.
    from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as solve_tsp
    

from cvrp_ops import fast_constraint_check, calculate_objective, normalize_solution
from local_search.solution_operators import do_3optstar_move
from local_search import LSOPT
from util import without_empty_routes, sol2routes, routes2sol
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def _make_checker_function(lambda1, lambda2):
    return 

def _check_lr3opt_move(D, C, L, removed_weights, best_delta,
                     edges, end_p, end_n, cum_d, cum_l,
                     ldepot_12, ldepot_34, lambdas):
    """ The Lagrangian constrant violation penalty calculation.
    
    The penalties for each of the 14 possible 3-opt* moves are calculated here.
    They are lifted as the edges are removed and recalculated when the new
    edges involved in the move are added back. This also updates the delta.

    .. note:: a good performance optimization would be to do this in parts when
    i,j,k are changed but that would mean duplication of the do_3optstar_move
    code, which was not preferable."""
    
    ## Calculate the penalties that are lifted
    removed_penalties = 0
    if C:
        route_d = cum_d[0]+cum_d[1]
        if ldepot_12:
            if route_d-C_EPS>C:
                removed_penalties+=(route_d-C)*lambdas[0]
            route_d = cum_d[2]+cum_d[3]
        else:
            route_d+=cum_d[3]
        if ldepot_34:
            if route_d-C_EPS>C:
                removed_penalties+=(route_d-C)*lambdas[0]
            route_d = cum_d[4]+cum_d[5]
        else:
            route_d += cum_d[5]
        if route_d-C_EPS>C:
            removed_penalties+=(route_d-C)*lambdas[0]
        
    if L:
        route_l = cum_l[0]+cum_l[1]+D[end_n[0],end_n[1]]
        if ldepot_12:
            if route_l-S_EPS>L:
                removed_penalties+=(route_l-L)*lambdas[1]
            route_l = cum_l[2]+cum_l[3]+D[end_n[2],end_n[3]]
        else:
            route_l+=cum_l[3]+D[end_n[2],end_n[3]]
        if ldepot_34:
            if route_l-S_EPS>L:
                removed_penalties+=(route_l-L)*lambdas[1]
            route_l = cum_l[4]+cum_l[5]+D[end_n[4],end_n[5]]
        else:
            route_l += cum_l[5]+D[end_n[4],end_n[5]]
        if route_l-S_EPS>L:
            removed_penalties+=(route_l-L)*lambdas[1]
    
    ## Calculate the delta (still missing the move induced penalties)
    added_weights = D[end_n[edges[0][0]], end_n[edges[0][1]]]+\
                    D[end_n[edges[1][0]], end_n[edges[1][1]]]+\
                    D[end_n[edges[2][0]], end_n[edges[2][1]]]
    delta = added_weights-removed_weights-removed_penalties
    
    ## Calculate the move induced Lagrangian relaxation penalties
    # Note that the move has to have some potential to be an improvement.
    if (delta+S_EPS<best_delta):
        
        prev_edge = None
        route_d = 0
        route_l = 0
        
        # Assume that the edges are in right order and to the right direction
        #  so that the final solution can be formed by applying them
        #  consecutively.
        for curr_edge in edges:
            # Check if there is a visit to the depot on the previous segment, 
            #  that is, between previous END and the current edge START nodes.
            if (curr_edge[0]==0 or prev_edge[1]==5 or
                (ldepot_12 is not None and (
                 end_p[prev_edge[1]]<=ldepot_12<=end_p[curr_edge[0]] or
                 end_p[prev_edge[1]]>=ldepot_12>=end_p[curr_edge[0]])) or 
                ( ldepot_34 is not None and (      
                end_p[prev_edge[1]]<=ldepot_34<=end_p[curr_edge[0]] or
                end_p[prev_edge[1]]>=ldepot_34>=end_p[curr_edge[0]]))):
                
                if C:
                    if route_d-C_EPS>C:
                        delta += (route_d-C)*lambdas[0]
                    route_d = cum_d[curr_edge[0]] 
                if L:
                    if route_l-S_EPS>L:
                        delta += (route_l-L)*lambdas[1]
                    route_l = cum_l[curr_edge[0]]
                
                # abort poor moves early
                if (delta+S_EPS>=best_delta):
                    return None
            if C:
                route_d += cum_d[curr_edge[1]]
            if L:
                e_wt = D[end_n[curr_edge[0]],end_n[curr_edge[1]]]
                route_l+=e_wt+cum_l[curr_edge[1]]                   
            prev_edge = curr_edge
        
        # The last edge has been connected, time to do check if the last formed
        #  route induces any penalties.
        if C and route_d-C_EPS>C:
            delta += (route_d-C)*lambdas[0]
        if L and route_l-S_EPS>L:
            delta += (route_l-L)*lambdas[1]
        if (delta+S_EPS<best_delta):
            return delta
        
    # was not an improving move
    return None

def _init_with_random(D,d,C,L):
    customers = list(range(1,len(D)))
    shuffle(customers)
    random_sol = [0]+customers+[0]
    return random_sol, calculate_objective(random_sol, D)

def _init_with_tsp(D,d,C,L):
    route_tsp_sol, route_f = solve_tsp(D, list(range(0,len(D))))
    return route_tsp_sol+[0], route_f

def _force_feasible(sol, D, d, C, L):
    # Force an incomplete solution feasible
    feasible_routes = []
    routes = sol2routes(sol)
    for r in routes:
        feasibler,totald,totall = [0], 0.0, 0.0
        prevn = 0
        for n in r[1:]:
            C_violated = (C and totald+d[n]-C_EPS > C)
            L_violated = (L and totall+D[prevn,n]+D[n,0]-S_EPS > L)
            if (n==0 or C_violated or L_violated) and len(feasibler)>1:
                feasibler.append(0)
                feasible_routes.append(feasibler)
                feasibler,totald,totall = [0], 0.0, 0.0
                prevn = 0
            if C:
                totald+=d[n]
            if L:
                totall+=D[prevn,n]
            if n!=0:
                feasibler.append(n)
            prevn = n
    return routes2sol(feasible_routes)

def _get_max(D, with_sol):
    aidx = np.array( np.unique( with_sol ) )
    return np.max( D[aidx][:, aidx] )

def lr3opt_init(D, d, C, L,
                initial_lambda1_C=None, initial_lambda1_L=None,
                initialization_algorithm=_init_with_tsp,
                postoptimize_with_3optstar=True,
                max_concecutive_lamba_incs=None):
    """ An implementation of the Stewart & Golden [1]_ 3-opt* heuristic
    with Lagrangean relaxation.
    
    The algorithm starts from a solution that can be either feasible or
    infeasible and uses local search to move towards better and feasible 
    solutions. More specifically, it works by replacing the constraint checks 
    of the 3-opt* with a penalty that depends on how much the constraint was 
    violated. The 3-opt* that operates on the entire solution, that is, checks 
    for both intra and inter route moves on one pass, was used. The penalties
    are iteratively doubled as the local search progresses and it is assumed
    that this eventually forces the solutions to feasible region.

    .. [1] Stewart, W. R. and Golden, B. L. (1984). A lagrangean relaxation 
           heuristic for vehicle routing. European Journal of Operational
           Research, 15(1):84–88.
    
    Parameters
    ----------
    D : numpy.ndarray
        is the full 2D distance matrix.
    d : list
        is a list of demands. d[0] should be 0.0 as it is the depot.
    C : float
        is the capacity constraint limit for the identical vehicles.
    L : float
        is the optional constraint for the maximum route length/duration/cost.
    
    initial_lambda1_C : float
        is the initial Langrange multiplier value for the capacity constraint C.
        If left empty (None) the formula ``l1_C=average(d)/(20*max(D))`` is used.
        The alternative value suggested by Stewart & Golden (1984) was 0.05.
    initial_lambda1_L : float
        is the initial Langrange multiplier value for the maximum route cost/
        duration/length constraint. If left empty (None) the formula
        ``l1_L=average(distance to nearest neighbor)/(10*max(D))`` is used.
    initialization_algorithm (function): is a function that retuns a TSP or VRP 
        solution and its objective function value. The default is to use LKH TSP 
        solution, but the function _init_with_random can be used to replicate the
        results of Stewart & Golden (1984) where a random solution is used.
    
    Returns
    -------
    list
        The solution as a list of node indices to visit.
    
    .. todo:: due to how the algorithm works, introducing minimize_K would require
    balancing between penalizing constraint violations and penalizing new 
    routes with an additional multipiler. This was not implemented.
    """
    
    
    sol = None
    try:
        ## STEP 1: Generate an initial solution
        sol, initial_f = initialization_algorithm(D,d,C,L)
        
        max_D = None
        lambdas = [initial_lambda1_C,initial_lambda1_L]
        if C and lambdas[0]==None:
            max_D = _get_max(D, sol)
            lambdas[0] = np.average(d)/(20*max_D)
        if L and lambdas[1]==None:
            # Stewart & Golden (1984) did not propose an extension for the maximum
            #  route duration/length/cost constraint, but here we have something 
            #  similar to L than they used for C constraint relaxation.
            max_D = _get_max(D, sol) if (max_D is None) else max_D
            closest_neighbor_D = D.copy()
            np.fill_diagonal(closest_neighbor_D, max_D)
            lambdas[1] = np.average(closest_neighbor_D.min(axis=0))/(10*max_D)

        if __debug__:
            log(DEBUG, "Start from initial solution %s (%.2f), and with l1=%.2f, l2=%.2f"%
                (sol, calculate_objective(sol, D),
                 (0 if lambdas[0] is None else lambdas[0]),
                 (0 if lambdas[1] is None else lambdas[1])))
       
        checker_function = partial(_check_lr3opt_move, lambdas=lambdas)
        
        # STEP 2: Solve the relaxed problem using 3-opt*
        c_lambda_incs = 0
        while True:
            # Make sure there is an empty route (for giving the 3-opt* procedure
            #  the option of adding vehicles)
            while not ( sol[-1]==0 and sol[-2]==0 ):
                sol+=[0]
            
            if __debug__:
                log(DEBUG-2, "Finding a LR3OPT move for %s (%.2f)"%
                    (sol, calculate_objective(sol, D)))
            new_sol, delta = do_3optstar_move(sol, D, d, C, L,
                         strategy=LSOPT.FIRST_ACCEPT, 
                         move_checker = checker_function)
            
            # local optima reached, tighten the relaxation
            # TODO: it should not happen that the sol==new_sol. However it happens and as a quickfix check for it.
            if delta is None or sol==new_sol:
                # return the first feasible solution (note: does not check for covering)
                if fast_constraint_check(sol,D,d,C,L):
                    if __debug__:
                        log(DEBUG, "Reached feasible solution %s (%.2f)"%
                            (sol, calculate_objective(sol,D)))
                    while postoptimize_with_3optstar:
                        opt_sol, delta = do_3optstar_move(sol, D, d, C, L,
                                         strategy=LSOPT.FIRST_ACCEPT)
                        if delta is None:
                            return normalize_solution(sol) # remove any [0,0]'s
                        else:
                            sol = opt_sol
                            #print("REMOVEME improved with post-optimization 3-opt*")
                            log(DEBUG, "Found improving 3-opt* move leading to %s (%.2f)"%
                                (sol, calculate_objective(sol,D)))                
                               
                    return normalize_solution(sol) # remove any [0,0]'s
                else:
                    # STEP 3: Update lambdas
                    lambda_at_inf = False
                    if lambdas[0] is not None:
                        lambdas[0] = lambdas[0]*2
                        lambda_at_inf = lambdas[0]==float('inf')
                    if lambdas[1] is not None:
                        lambdas[1] = lambdas[1]*2
                        lambda_at_inf = lambda_at_inf or lambdas[0]==float('inf')
                    if __debug__:
                        log(DEBUG-1, "No improving moves left, increasing lambda to l1=%.2f, l2=%.2f"%
                            ((0 if lambdas[0] is None else lambdas[0]),
                             (0 if lambdas[1] is None else lambdas[1])))
                    #print("No improving moves left, increasing lambda to l1=%.2f, l2=%.2f"%
                    #        ((0 if lambdas[0] is None else lambdas[0]),
                    #         (0 if lambdas[1] is None else lambdas[1])))
                        
                    #TODO: if penalty >> cost, break (stuck on a infeasible region)
                    # how much bigger can be determined by finding such a
                    # pathological problem instance?

                    # safeguard for getting stuck
                    c_lambda_incs+=1
                    #print("REMOVEME: c_lambda_incs", c_lambda_incs)
                    if lambda_at_inf or (max_concecutive_lamba_incs is not None and
                                         c_lambda_incs > max_concecutive_lamba_incs):
                        return _force_feasible(sol, D, d, C, L)

            else:
                if __debug__:
                    log(DEBUG, "Found improving LR3OPT move leading to %s (%.2f)"%
                        (new_sol, calculate_objective(new_sol,D)))
                    log(DEBUG-2, "However, routes %s remain infeasible."%
                        [r for r in sol2routes(new_sol) if not fast_constraint_check(r,D,d,C,L)])
                        
                    
                sol = new_sol
                c_lambda_incs = 0



    except KeyboardInterrupt: # or SIGINT
        # Pass on the current solution forced feasbile by splitting routes
        #  according to the constraints.
        raise KeyboardInterrupt(_force_feasible(sol, D, d, C, L))

    return without_empty_routes(sol)
                
                
# Wrapper for the command line user interface (CLI)
def get_lr3opt_algorithm():
    algo_name = "SG84-LR3OPT"
    algo_desc = "Stewart & Golden (1984) Lagrangian relaxed 3-opt* heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if minimize_K:
            raise NotImplementedError("LR3OPT does not support minimizing the number of vehicles")
        return lr3opt_init(D, d, C, L)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_lr3opt_algorithm())
