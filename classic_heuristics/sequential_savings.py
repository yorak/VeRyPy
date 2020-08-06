#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the (parametrized)
sequential savings heuristic of Webb (1964). In contrast to the more
popular parallel savings heuristic the sequential version build the routes one
at the time.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has minimal dependencies: only numpy and scipy are
required to be installed for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from logging import log, DEBUG
from util import objf, routes2sol

from collections import deque

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

def clarke_wright_savings(unrouted, i, D):
    savings = [(D[i,0]+D[0,j]-D[i,j],-D[i,j],i,j) for j in unrouted if i!=j]
    savings.sort()
    return savings 

def sequential_savings_init(D, d, C, L=None, minimize_K=False,
                          initialize_routes_with = "closest",
                          savings_callback=clarke_wright_savings):
    """
    Implementation of the Webb (1964) sequential savings algorithm /
    construciton heuristic for capaciated vehicle routing problems with
    symmetric distances.
    
    This is the sequential route version, which builds the solution one route 
    at the time by making always the best possible merge for that active route.
    
    It has been implemented as separate function from the parallel version, 
    because their working principle and datastrucutres differ significantly.
    
    The parameters are mostly the same as in parallel_savings.py:
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/cost/duration.
    
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring negative savings values.
    
    * optional savings_callback is a function of the signature
        sorted([(s_i1,x_i1,i,j_1)...(s_ij,x_ij,i,j)...(s_in,x_in,i,n) ]) =
            savings_callback(unrouted, i, D)
      where the returned (sorted!) list contains savings s_ij (that is, how
       much route cost approximately decreases if a merge with an edge (i,j)
       is done). The x_ij is a secondary sorting criterion and it is ignored 
       by the savings heuristic. Note that the signature is different to the 
       parallel version  as the savings calculation is done only for the
       currently unrouted  customers and for customer i. 
      The default is a the Clarke Wright savings criterion.
       
	Webb, M. (1964). A study in transport routing. Glass Technology, 5:178181
    """
    
    N = len(D)
    ignore_negative_savings = not minimize_K
    
    unrouted =  set(range(1,N))
    
    ## 1. Generate a list of seed nodes for emerging route inititialization 
    
    if initialize_routes_with=="farthest" or initialize_routes_with=="closest":
        # build a ordered priority queue of potential route initialization nodes
        seed_customers = list(range(1,N))
        seed_customers.sort(reverse = initialize_routes_with=="closest",
                            key = lambda i: (D[0][i],i) )
    elif initialize_routes_with=="first":
        seed_customers = list(range(1,N))
    elif initialize_routes_with=="last":
        seed_customers = list(range(N-1,-1,-1))
    elif initialize_routes_with=="savings":
        # Calculate the savings A i,j, sort them, and generate a j list
        seed_customers = sum(
            (clarke_wright_savings(unrouted, i, D) for i in range(1,N)), [])
        seed_customers.sort()
        _,_,_,seed_customers = zip(*seed_customers)
        seed_customers = list(seed_customers)
    elif callable(initialize_routes_with):
        seed_customers = initialize_routes_with(D)
    else:
        raise ValueError("No such route initialization method")
    
    ## 2. Initialize a single emerging route and make all feasible merges on it
    
    solution = [0]
    savings = None
    emerging_route_nodes = None
    try:
        while unrouted:
            # start a new route
            if not savings:
                while True:
                    seed = seed_customers.pop()
                    if seed in unrouted:
                        break
                    
                emerging_route_nodes = deque([seed])
                unrouted.remove(seed)
                
                route_d = d[seed] if C else 0.0
                route_l = D[0,seed]+D[seed,0] if L else 0.0
                
                dbg_route = [0]+list(emerging_route_nodes)+[0]
                log(DEBUG-1, "Initialize a new route as %s (%.2f)"%
                             (str(dbg_route), objf(dbg_route,D)) )
            
                savings = savings_callback(unrouted, seed, D)
                
            while len(savings)>0:
                # Note: i is the one to merge with, and j is the merge_candidate
                best_saving, _, i, j = savings.pop()
                if __debug__:
                    log(DEBUG-1, "Popped savings s_{%d,%d}=%.2f" % (i,j,best_saving))
                
                if ignore_negative_savings:
                    cw_saving = D[i,0]+D[0,j]-D[i,j]
                    if cw_saving<0.0:
                        savings = [] # force route change
                        break
                    
                if not j in unrouted:
                    continue
     
                if C and route_d+d[j]-C_EPS > C:
                    if __debug__:
                        log(DEBUG-1, "Reject merge due to C constraint violation")
                    continue #next savings
                    
                # it is still a valid merge?
                do_left_merge = emerging_route_nodes[0] == i
                do_right_merge = emerging_route_nodes[-1] == i and\
                                 len(emerging_route_nodes)>1
                if not (do_left_merge or do_right_merge):
                    if __debug__:
                        log(DEBUG-1,("Reject merge because the %d"%i)+
                                      "is no longer next to the depot.")
                    continue #next savings
                   
                if L:
                    # symmetric distances allow left and right merge use same check
                    l_delta = D[0,j]+D[i,j]-D[0,i]
                    if route_l+l_delta-S_EPS>L:
                        if __debug__:
                            log(DEBUG-1, "Reject merge due to L constraint violation")
                        continue #next savings
                    # it is OK, update
                    route_l+=l_delta
                if C: route_d+=d[j]
                
                if do_left_merge:
                    emerging_route_nodes.appendleft(j)
                if do_right_merge:
                    emerging_route_nodes.append(j)
                unrouted.remove(j)
                
                # update the savings list
                savings += savings_callback(unrouted, j, D)
                savings.sort()
                
                if __debug__:
                    dbg_route = [0]+list(emerging_route_nodes)+[0]
                    log(DEBUG-1, "Merged, resulting route is %s (%.2f)"%
                                 (str(dbg_route), objf(dbg_route,D)) )
            
            if __debug__:
                dbg_route = [0]+list(emerging_route_nodes)+[0]
                log(DEBUG, "Route %s (%.2f) COMPLETE"%
                             (str(dbg_route), objf(dbg_route,D)) )
            
            # All savings merges tested, complete the route
            emerging_route_nodes.append(0)
            solution+=emerging_route_nodes
            emerging_route_nodes = None
            
    except KeyboardInterrupt:
        interrupted_solution = solution
        if emerging_route_nodes:
            emerging_route_nodes.append(0)
            interrupted_solution+=emerging_route_nodes
        interrupted_solution+=routes2sol([n] for n in unrouted)[1:]
        raise KeyboardInterrupt(interrupted_solution)
                 
    return solution


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ss_algorithm(lambda_multiplier='auto'):
    algo_name = "We64-SS"
    algo_desc = "Webb (1964) sequential savings algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return sequential_savings_init(D,d,C,L,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ss_algorithm())
