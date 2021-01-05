#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the route-first, cluster-
second heuristic of Newton Thomas (1969) and Beasley (1983).

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has moderate dependencies: a TSP solver. LKH and the built-
in local search solvers are used, but if needed, only the built-in solver will
suffice just fine. Furhtermore, numpy and scipy for reading and preparing the 
problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range
from sys import stderr

import numpy as np
from logging import log, DEBUG

try:
    # The original implementation of Beasley (1983) generated multiple 2-optimal initial
    #  solutions. We just generate one that is very, very good with state-of-the-art TSP
    #  solver LKH.
    from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp_intial
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as solve_tsp_intial
# Always use built-in TSP solver to guarantee 2 and 3-optimality of single routes.
from tsp_solvers.tsp_solver_ropt import solve_tsp_2opt, solve_tsp_3opt

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS


__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

def _partition_to_routes_with_cost_matrix(D, d, C, L, giant_tour,
                                          route_cost_constant,
                                          partition_tsp_opt_algo,
                                          route_tsp_opt_algo):
    """ Partitions the giant tour solution optimally to routes using
    Floyd-Warshall algorithm.
    """
    
    N = len(D) 
    gtl = len(giant_tour)
    if giant_tour[0]==giant_tour[-1]:
        gtl = gtl-1
    
    # these are the dist and next tables of Floyd-Warshall algorithm
    sequence_D = np.full((N,N), np.inf)
    next_split = np.zeros((N,N), dtype='int')
    interrupted = False
    
    try:
        ## Build the feasible sequences of nodes on the giant tour
        u_idx = 0
        while True:
            cum_demand = 0.0
            
            v_idx = u_idx+1
            if v_idx==gtl:
                v_idx = 0
                
            # remember prev route to make local search initial solution better 
            #  (and thus, convergence faster)
            route_nodes = [0, 0]
            while True:
                v = giant_tour[v_idx]
                
                # C constraint would be violated, increment u_idx
                if (C and cum_demand+d[v]-C_EPS>C):
                    break
                
                 # only start checking if the TSP giant tour would break the constraint
                route_nodes[-1] = v
                
                # recycle the tsp result, to save few local search cycles
                route_nodes, route_cost = partition_tsp_opt_algo(D, route_nodes)
                
                # L constraint would be violated, increment u_idx
                if L and route_cost-S_EPS>L:
                    break
                if C: cum_demand+=d[v]
                
                # "cost of serving [u+1...v] **after** u"
                sequence_D[u_idx,v_idx] = route_cost+route_cost_constant
                next_split[u_idx,v_idx] = v_idx
                
                # leave out the last node of the giant tour (same as the first)
                v_idx+=1
                if v_idx==gtl:
                    v_idx = 0
                if v_idx==u_idx:
                    break
            
            # leave out the last node of the giant tour (same as the first)
            u_idx+=1
            if u_idx==gtl:
                break
    except KeyboardInterrupt:
        interrupted = True
        
    # Solve using Floyd's algorithm
    #  note that we cannot use floyd_warshall from scipy.sparse.csgraph, as 
    #  this initializes the i==i to 0, which is exactly the distance we are 
    #  interested in!
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if sequence_D[i,j] > sequence_D[i,k]+sequence_D[k,j]:
                    sequence_D[i,j] = sequence_D[i,k]+sequence_D[k,j]
                    next_split[i,j] = next_split[i,k]
    
    # build a solution based on the Floyd's algorithm *next* matrix for the 
    #  least cost tour, 
    best_start_idx = np.argmin(sequence_D.diagonal())
    sol = [0]
    sol_K = 0
    sol_f = 0.0
    
    u_idx = best_start_idx 
    v_idx = next_split[u_idx, best_start_idx]
    while True:
        
        # Build the route in same way as with the shortcut matrix i.e. adding 
        # cusotmers to the existing TSP solution and re-optimizing it after 
        # each step with partition_tsp_opt_algo. Also do it here to avoid a
        # rare issue where e.g. 3-optimal first accept result may be worse than 
        # the 2-optimal one for the shortcut matrix which may cause L constr.
        # violation.
        route = [0,0]
        while len(route)==2 or u_idx!=v_idx:
            if u_idx+1 == gtl:
                u_idx = -1
            route[-1] = giant_tour[u_idx+1]
            if not interrupted:
                route, _ = partition_tsp_opt_algo(D, route)
            else:
                route.append(0)
            u_idx+=1
        
        if not interrupted:
            route, route_cost = route_tsp_opt_algo(D, route )
        
        # record the route
        sol.extend(route[1:])
        sol_K += 1
        sol_f += route_cost
        
        # take next route
        u_idx = v_idx
        v_idx = next_split[u_idx, best_start_idx]
        if (u_idx==best_start_idx):
            break
    
    if __debug__: 
        if not interrupted:
            log(DEBUG, "Partitioning result was %s, K = %d, f = %.2f"%
                       (str(sol), sol_K, sol_f))
    
    if interrupted:
        raise KeyboardInterrupt(sol)
    
    return sol

def route_first_cluster_second_init(D, d, C, L=None, minimize_K=False,
                tsp_gen_algo=solve_tsp_intial,
                partition_tsp_opt_algo=solve_tsp_2opt,
                route_tsp_opt_algo=solve_tsp_3opt):
    
    """ Implementation of the Route First Cluster Second heuristic for solving
    symmetric CVRPs. The idea has been presented e.g. in Beasley (1983).
    However this implementation is deterministic and does not rely on
    generating several random giant tours.
    
    The algorithm works as follows. First, a the problem is solved as a TSP.
    Then the giant tour TSP solution is partitioned into paths of consecutive
    nodes that satisfy problem constraints e.g. the capacity constraint. From
    the all valid parititions the one with shortest total travel cost is
    selected as the optimized solution. The method based on Floyd-Warshall
    algorithm like in from Beasley (1983) is used to find a shortest path of 
    "shortcuts" (i.e. VRP routes) through this giant tour.
    
    The algorithm expects following parameters:
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/cost/duration.
    * partitioning_method sets how partitions are made. I
       
    Furthermore, the algorithm can be exdended by using a different giant tour
    generation algorithm tsp_gen_algo and a cluster routing algorihm
    partition_opt_algo. Both have a signature of algo(D, node_idxs), where D is the 
    (full) distance matrix and node_idxs a list of node indices to route.
    Giving a random giant tour initialization with 2-opt improvement step as
    tsp_gen_algo and 3-opt improvement as tsp_opt_algo the results of 
    Beasley (1983) can be replicated.

    As in Beasley (1983), we do not include the depot when solving the TSP.
    Also, in this implementation relies on external solver solving the TSP,
    which is by default the estate-of-the-art LKH (Helsgaun 2006 & 2009). 
    
    Beasley, J. E. (1983). "Route first—cluster second methods for vehicle
      routing". Omega, 11(4), 403-408.
    Helsgaun, K. (2006), "An Effective Implementation of K-opt Moves for the 
      Lin-Kernighan TSP Heuristic." DATALOGISKE SKRIFTER, No. 109, 2006. 
      Roskilde University.
    Helsgaun, K. (2009), "General k-opt submoves for the Lin-Kernighan TSP 
      heuristic." Mathematical Programming Computation, 1(2), 119-163.
    """
    N = len(D)
    if N<=2:
        return list(range(N))+[0]
    
    giant_tour_sol, giant_tour_l  = tsp_gen_algo(D, list(range(1,N)))
        
    if __debug__:
        log(DEBUG, "TSP tour solution %s (%.2f)"%(str(giant_tour_sol),giant_tour_l))
        
    route_cost_constant = giant_tour_l if minimize_K else 0.0    
    return _partition_to_routes_with_cost_matrix(D, d, C, L,
                                giant_tour_sol, route_cost_constant, 
                                partition_tsp_opt_algo, route_tsp_opt_algo)


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_rfcs_algorithm():
    algo_name = "Be83-RFCS"
    algo_desc = "Route-first-cluster-second heuristic of Beasley (1983)"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return route_first_cluster_second_init(D, d, C, L, minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_rfcs_algorithm())
        
        
    

    
    
        
