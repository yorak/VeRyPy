#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Desrochers and Verhoog
(1989) matchings based savings algorithm. In the algoritm a weighted maximum
matching problem is solved repeadedly to find the best feasible merge.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has extensive dependencies: MIP solver Gurobi, a TSP
solver (the built in local search solver can be used), and numpy and scipy for
reading and preparing the problem instance.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range
from sys import stderr

from logging import log, DEBUG
from signal import signal, SIGINT, default_int_handler

# this uses Gurobi to solve the maximum matching problem
from gurobipy import Model, GRB, GurobiError                         

try:
    # The algoritm uses Gurobi to solve the TSPs of the maximum matching problem.
    from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as default_solve_tsp
    ## For very large instances one can optionally use faster LKH or acotsp local search.
    #from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as default_solve_tsp  
    #from tsp_solvers.tsp_solver_acotsp import solve_tsp_acotsp as default_solve_tsp
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as default_solve_tsp

from util import objf

from config import MAX_MIP_SOLVER_RUNTIME, MIP_SOLVER_THREADS
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"
    
def _mmp_add_cnts_sum(m, x_ij, x_ij_keys, w_ij, n):
    ## constraints
    for k in range(0,n):
        active_vars = x_ij.sum(k,'*')
        if active_vars.size()>0:
            m.addConstr( active_vars == 1, "m_%d"%k )

def _mmp_solve(w1_ij, x_ij_keys, n, w2_ij = None):
    """A helper function that solves a weighted maximum matching problem.
    """
    
    m = Model("MBSA")
    
    if __debug__:
        log(DEBUG,"")
        log(DEBUG,"Solving a weighted maximum matching problem with "+
                  "%d savings weights." % len(w1_ij))

    # build model
    x_ij = m.addVars(x_ij_keys, obj=w1_ij, vtype=GRB.BINARY, name='x')    
    _mmp_add_cnts_sum(m, x_ij, x_ij_keys, w1_ij, n)
    m._vars = x_ij
    m.modelSense = GRB.MAXIMIZE
    m.update()

    # disable output
    m.setParam('OutputFlag', 0)    
    m.setParam('TimeLimit', MAX_MIP_SOLVER_RUNTIME)
    m.setParam('Threads', MIP_SOLVER_THREADS)
    #m.write("out.lp")
    m.optimize()

    # restore SIGINT callback handler which is changed by gurobipy
    signal(SIGINT, default_int_handler)

    if __debug__:
        log(DEBUG-1, "Gurobi runtime = %.2f"%m.Runtime)
    
    if m.Status == GRB.OPTIMAL:
        if w2_ij==None:
            max_wt, max_merge = max( (w1_ij[k], x_ij_keys[k])
                                      for k, v in enumerate(m.X) if v )
        else:
            max_wt, _, max_merge = max( (w1_ij[k], w2_ij[k], x_ij_keys[k])
                                      for k, v in enumerate(m.X) if v )
            
        return max_wt, max_merge[0], max_merge[1]
    elif m.Status == GRB.TIME_LIMIT:
        raise GurobiError(10023, "Gurobi timeout reached when attempting to solve GAP")
    elif m.Status == GRB.INTERRUPTED:
        raise KeyboardInterrupt()
    return None

def _get_tsp_sol(nodes, cache, D, solve_tsp):
    if not nodes in cache:
        cache[nodes] = solve_tsp(D, list(nodes))
    return cache[nodes]

def _get_demand(nodes, cache, d):
    if not nodes in cache:
        cache[nodes] = sum( d[n] for n in nodes )
    return cache[nodes]

def _calculate_secondary_criteria(w_ij, x_ij_keys, n):
    # calcuate potential improvements for each route
    total_route_wt = [0]*n
    for k, wt in enumerate(w_ij):
        i,j = x_ij_keys[k]
        total_route_wt[i] += wt
        total_route_wt[j] += wt
    
    secondary_criteria = [-(total_route_wt[x_ij_keys[k][0]]+
                            total_route_wt[x_ij_keys[k][1]]-2*wt)
                            for k, wt in enumerate(w_ij)]
    return secondary_criteria


def _calculate_savings(rs1,rs2,demand_cache,tsp_cache,W,D,d,C,L,solve_tsp):    
    # Check that it is C feasible
    if C:
        s1s2_d = _get_demand(rs1,demand_cache,d)+_get_demand(rs2,demand_cache,d)
        if s1s2_d-C_EPS>C:
            return W # capacity constraint violated
    
    # Check that it is L feasible    
    s1s2 = rs1.union(rs2)
    s1s2_len = _get_tsp_sol(s1s2, tsp_cache, D, solve_tsp)[1]
    if L and s1s2_len-S_EPS>L:
        return W # max. tour len. constraint violated
    
    # C and L feasible, pass on the savings value
    s = _get_tsp_sol(rs1, tsp_cache, D, solve_tsp)[1]\
        +_get_tsp_sol(rs2, tsp_cache, D, solve_tsp)[1]\
        -s1s2_len
    return s

def _geedy_merge(D, d, C, L, W, savings, route_sets, tsp_cache, demand_cache):
    """ This is not the part of the main algorithm proper but a greedy fallback 
    the procedure relies on in case of an interrupt. It simply uses the 
    existing  savings list and caches to create solution with minimal amount of
    computation (it is an O(n) algorithm at this point).
    """
    savings.sort(reverse=True)
    
    try:    
        
        # Do a first pass where the routes in the tsp_cache are used
        first_pass_joined_routes = set()
        for saving, route_ids in savings:
            if saving==W:
                break # from now one all are infeasible
            
            route_i, route_j = route_ids
            if route_i in first_pass_joined_routes or\
               route_j in first_pass_joined_routes:
                continue # either route is already merged
            
            nodes_i = route_sets[route_i]
            nodes_j = route_sets[route_j]
            nodes_ij = nodes_i.union(nodes_j)
            
            if nodes_ij not in tsp_cache:
                greedy_routing = [0]+list(nodes_ij)+[0]
                greedy_f = objf(greedy_routing, D)
                if not L or greedy_f-S_EPS<L:
                    # remember for the next stage
                    tsp_cache[nodes_ij] = (greedy_routing, greedy_f)
                continue # all in tsp_cache are feasible
            
            # record and do the merge
            first_pass_joined_routes.add(route_i)
            first_pass_joined_routes.add(route_j)
            route_sets[route_i] = nodes_ij
            route_sets[route_j] = route_i # mark as a reference to the another
            
        # the second stage just smash routes together (in all possible 4-ways)
        for saving, route_ids in savings:
            # follow the merge path (if any)
            route_i, route_j = route_ids
            while isinstance(route_sets[route_i], int):
                route_i = route_sets[route_i]
            while isinstance(route_sets[route_j], int):
                route_j = route_sets[route_j]
            if route_i==route_j: continue #already merged
        
            # get nodes, TSP chain, demand, and lenght for both routes
            nodes_i = route_sets[route_i]
            nodes_j = route_sets[route_j]
            chain_i = tsp_cache[nodes_i][0][1:-1]
            l_i = tsp_cache[nodes_i][1]
            d_i = _get_demand(nodes_i, demand_cache, d)        
            chain_j = tsp_cache[nodes_j][0][1:-1]
            l_j = tsp_cache[nodes_j][1]
            d_j = _get_demand(nodes_j, demand_cache, d)
            
            # check constraints
            if C and d_i+d_j-C_EPS>C: continue
            na, nb, nc, nd = chain_i[0], chain_i[-1], chain_j[0], chain_j[-1]
            alt_costs = [l_i+l_j+D[nb,nc]-D[nb,0]-D[nc,0], #ab-cd
                         l_i+l_j+D[na,nc]-D[na,0]-D[nc,0], #ba-cd
                         l_i+l_j+D[nb,nd]-D[nb,0]-D[nd,0], #ab-dc
                         l_i+l_j+D[na,nd]-D[na,0]-D[nd,0]] #cd-ab
            min_cost = min(alt_costs)
            min_alt = alt_costs.index(min_cost)                                
            if L and min_cost-S_EPS>L: continue
        
            # find the best merge among the 4-way merges and apply it
            if min_alt==0: ij_route = [0]+chain_i+chain_j+[0] #ab-cd
            if min_alt==1: ij_route = [0]+chain_i[::-1]+chain_j+[0] #ba-cd
            if min_alt==2: ij_route = [0]+chain_i+chain_j[::-1]+[0] #ab-dc
            if min_alt==3: ij_route = [0]+chain_j+chain_i+[0] #cd-ab
            nodes_ij = nodes_i.union(nodes_j)
            route_sets[route_i] = nodes_ij
            route_sets[route_j] = route_i # mark as a reference to the another
            tsp_cache[nodes_ij] = (ij_route, min_cost)
    except KeyboardInterrupt: #or SIGINT
        pass
    
    greedy_solution = []
    for nodes in route_sets:
        if isinstance(nodes, int): continue
        greedy_solution +=tsp_cache[nodes][0][:-1]
    greedy_solution.append(0)
    return greedy_solution 
        
def mbsa_init(D, d, C, L, minimize_K=False, W=0.0, solve_tsp=default_solve_tsp,
              primary_criteria_callback=_calculate_savings,
              secondary_criteria_callback=_calculate_secondary_criteria):
    """ An implementation of Desrochers & Verhoog (1989) Matching Based Savings 
    Algortihm. It uses Gurobi to solve the maxumum matching problem (MMP) and
    the TSP. The routes are merged according to the MMP until no valid merges
    remain.
    
    The parameters for this implementation are:
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.
    
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring negative savings values.
      
    * W is the savings cost assigned to infeasible merges. One may need to 
        set this to a large negative value if minimize_K is true.
    
    * solve_tsp which TSP solver is used when calculating the savings values
    * primary_criteria_callback and secondary_criteria_callback can be changed
       from the defaults to use different savings criteria. See the reference
       implementations for details on the fucntion signatures.

    Desrochers, M. and Verhoog, T. (1989). G-89-04 : A matching based
    savings algorithm for the vehicle routing problem. Technical report,
    GERAD, Montreal, Canada.
    """
    
    
    N = len(D)
    tsp_cache = {}
    demand_cache = {}
    savings = []
    
    ignore_negative_savings = not minimize_K
    
    ## Step 0: initalization
    # serve each customer with a single route
    route_sets = [frozenset([0,i]) for i in range(1,N)]

    try:    
        # calculate initial savings or merging *routes* i and j
        savings = [ (primary_criteria_callback(rs1,rs2,demand_cache,tsp_cache,
                                        W,D,d,C,L,solve_tsp), (i,j))
                    for i, rs1 in enumerate(route_sets)
                    for j, rs2 in enumerate(route_sets) if i<j ]
        
        while True:
            ## Step 1: Evaluate the weights and solve the weighted matching prolbem
            
            w_ij,x_ij_keys = zip(*savings)
            valid_matchings = any( w for w in w_ij if w!=W )
            if valid_matchings==0:
                break
            
            s_ij = None
            if secondary_criteria_callback:
                s_ij = secondary_criteria_callback(w_ij, x_ij_keys, len(route_sets))
            best_matching = _mmp_solve(w_ij,x_ij_keys,len(route_sets), s_ij)   
            if best_matching is None:
                break # no valid matchings found
     
            best_savings, i_prime,j_prime = best_matching
            
            if best_savings==W:
                break
            
            if ignore_negative_savings and best_savings<0:
                continue # do not allow merges that would make the solution worse
                         # (even if it would mean fewer routes)
    
            if __debug__:
                log(DEBUG,"Best matching joins routes %s (idx:%d) and %s (idx:%d)."%(
                        tsp_cache[route_sets[i_prime]][0],i_prime,
                        tsp_cache[route_sets[j_prime]][0],j_prime))
            
            ## Step 2: merge the route combination
            
            # merge routes i and j
            new_i_prime_set=route_sets[i_prime].union(route_sets[j_prime])
            
            ## Step 3: Update savings
            
            ## Filter out all edges adjacent to i' and j'
            # select only the savings that *do not* contain routes i and j
            #  and while at it, update the indexing
            savings = [ (sa,(si if si<j_prime else si-1,
                             sj if sj<j_prime else sj-1))
            
                        for sa,(si,sj) in savings
                        if (si!=i_prime and sj!=i_prime and
                            si!=j_prime and sj!=j_prime) ]
            
            del route_sets[j_prime]
            if i_prime>j_prime: i_prime=i_prime-1
            route_sets[i_prime]=new_i_prime_set
            
            # Evaluate the savings associated with the new route 
            rs1 = route_sets[i_prime]
            merged_route_savings = [ (primary_criteria_callback(rs1,rs2,
                                                         demand_cache,tsp_cache,W,
                                                         D,d,C,L,solve_tsp),(i_prime,k))
                                     for k, rs2 in enumerate(route_sets) if k!=i_prime]
            savings.extend(merged_route_savings)
            
            if __debug__:
                dbg_r, dbg_w = _get_tsp_sol(rs1, tsp_cache, D, solve_tsp)
                log(DEBUG-1,"Route#%d after merge %s (%.2f)"%
                    (i_prime, str(dbg_r), objf(dbg_r, D)))

    except KeyboardInterrupt: #or SIGINT
        interrupted_sol = _geedy_merge(D, d, C, L, W, savings, route_sets,
                                       tsp_cache, demand_cache)
        raise KeyboardInterrupt(interrupted_sol)
        
    # the optimized TSP tours for the routes are cached, reuse
    final_routes = [_get_tsp_sol(rs, tsp_cache, D, solve_tsp)[0] for rs in route_sets]
    sol = [0]+[n for route in final_routes for n in route[1:]]
    
    return sol 
        
# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_mm_algorithm():
    algo_name = "DV89-MM"
    algo_desc = "Desrochers and Verhoog (1989) maximum matching based "+\
                "savings algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return mbsa_init(D, d, C, L, minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_mm_algorithm())
