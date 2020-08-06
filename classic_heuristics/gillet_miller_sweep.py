#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Gillet & Miller (1974)
Sweep heuristic.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs.It has moderate dependencies: the sweep procedure from
sweep.py, the OrderedSet module, and a TSP solver (the built-in one based on 
local search can be used). As usual, numpy and scipy are required for reading
and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import numpy as np
from logging import log, DEBUG

# One of the few non standard-lib additions, used in algorithm's second 
#  phase to keep track of the non-routed nodes. Install it e.g. with:
#
# $ pip install orderedset
#
from orderedset import OrderedSet

from routedata import RouteData
from util import produce_nn_list
from classic_heuristics.sweep import _step,sweep_init,BEST_ALTERNATIVE

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

#from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp
from tsp_solvers.tsp_solver_ropt import solve_tsp_3opt as solve_tsp

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"



def _shortest_path_through_nodes(D, from_node, to_node, through_nodes):
    # make sure to_node is the last by setting the distance between the end
    # points to 0.0. Also, instead of copying (potentially large) D, just.
    #  switch out the values in the D for the TSP algorithm.
    if not through_nodes:
        return [from_node,to_node], D[from_node, to_node]

    old_path_end_distance = D[from_node, to_node]
    D[from_node, to_node] = 0.0
    D[to_node, from_node] = 0.0
    route,cost = solve_tsp(D, [from_node]+list(through_nodes)+[to_node])
    # restore original weights
    D[from_node, to_node] = old_path_end_distance
    D[to_node, from_node] = old_path_end_distance
    
    # contingency plan for the special case where there are nodes that overlap 
    #  the depot. Check for it and fix it with a simple swap if needed.
    if (route):
        if route[-2]!=to_node:
            last_node_idx = route.index(to_node)
            route[last_node_idx]=route[-2]
            route[-2]=to_node     
            if __debug__:
                log(DEBUG-1, "WARNING: the chain does not end to node "+
                           "(%d). Swapped two nodes to get route %s"%
                           (to_node,route))
                
    return route[:-1], cost
    
def _pack_datastructures_callback(D,d,C,L,sweep_phi_rho_nodes):
    """ This callback is called in the beginning of the the sweep(...) 
    the heuristic (in our case Gillett and Miller improvement heuristic)
    can then pick, choose and preprocess the data. Again, in this case 
    a nearest neighbor list in constructed only once, and the same data 
    structure can then be used in all route improvement calls.
    """
    pos_to_node = dict((i,int(n)) for i,n in enumerate(sweep_phi_rho_nodes[2]))
    node_to_pos = dict((int(n),i) for i,n in enumerate(sweep_phi_rho_nodes[2]))
    avr = np.average( sweep_phi_rho_nodes[1] )
    return (len(D), D,produce_nn_list(D),d,C,L,node_to_pos,pos_to_node, avr)
    
def _improvement_callback(route_data, callback_datastructures,
                          sweep_rhos, sweep_phis,
                          sweep_J_pos, step_inc, routed):
    """ This callback implements the Gillett and Miller (1974) improvement
    heuristic. It involves trying to remove a node and insertion of several 
    candidates. Therefore, the algorithm involves Steps 8-15 in the appedix
    of Gillett and Miller (1974). """

    # unpack callback data structures (packed in pack_datastructures)
    N, D, NN_D, d, C, L, node_to_pos, pos_to_node, avr = callback_datastructures
    
    # do not try to improve, if the J node is already routed (already Swept
    #  full circle)
    node_J =  pos_to_node[sweep_J_pos]
    if routed[node_J]:
        # nothing to do, no nodes added, no nodes removed, route is complete 
        return route_data, [], [], True
    
    # unpack route information, (route, route_cost, route_demand)    
    D1_route, D1, D1_demand, D1_nodes = route_data
    
    
    # G&M Step 8. This messy looking line vectorizes the minimization of
    #  R(K(I))+An(K(I))*AVR 
    #
    # What makes it a little more messier, is that the angle is pointing 
    #  "the right way" depending on the cw/ccw direction (encoded in step_inc)
    # +1 is there to convert indexing as there is no depot in node_rho_phis
    route_K_nodes = list(D1_nodes[1:])
    route_K_positions = [node_to_pos[n] for n in route_K_nodes]   
    route_rhos = sweep_rhos[route_K_positions]
    route_phis = sweep_phis[route_K_positions]
    
    rem_choose_function = route_rhos+route_phis*avr
    to_remove_node_KII = route_K_nodes[np.argmin(rem_choose_function)]
    
    if __debug__:
        log(DEBUG-2, "G&M improvement phase for route %s (%.2f). Trying to replace KII=n%d."%
            (str(D1_route), D1, to_remove_node_KII))
        log(DEBUG-3, "This is due to R(K(I))+An(K(I)*AVR = %s"%
            str(zip(route_K_nodes,route_phis,list(rem_choose_function))))    
    # take the node J-1 (almost always the node last added on the route)
    sweep_prev_of_J_pos = _step(sweep_J_pos, -step_inc, N-2)
    prev_node_J =  pos_to_node[sweep_prev_of_J_pos]
    
    
    # Get the insertion candidates
    try:
        candidate_node_JJX = next( (node_idx for node_idx,_ in NN_D[prev_node_J]
            if not routed[node_idx]) )
    except StopIteration:
        if __debug__:
            log(DEBUG-2,"G&M Step 9, not enough unrouted nodes left for JJX.")
            log(DEBUG-2,"-> EXIT with no changes")
        return route_data, [], [], False    
    try:
        candidate_node_JII = next( (node_idx for node_idx,_ in NN_D[candidate_node_JJX]
            if (not routed[node_idx] and node_idx!=candidate_node_JJX)) )
    except StopIteration:
        candidate_node_JII = None    

    # construct and route to get the modified route cost D2
    D2_route_nodes = OrderedSet(D1_nodes)
    D2_route_nodes.remove(to_remove_node_KII)
    D2_route_nodes.add(candidate_node_JJX)
    D2_route,D2 = solve_tsp(D, list(D2_route_nodes))
    D2_demand = D1_demand-d[to_remove_node_KII]+d[candidate_node_JJX] if C else 0  
    
    ## G&M Step 9    
    
    if not ((L and D2-S_EPS<L) and (C and D2_demand-C_EPS<=C)):
        if __debug__:
            log(DEBUG-2,"G&M Step 9, rejecting replacement of KII=n%d with JJX=n%d"%
                  (to_remove_node_KII, candidate_node_JJX))
            log(DEBUG-3," which would have formed a route %s (%.2f)"%(str(D2_route), D2))
            if (C and D2_demand-C_EPS>C):
                log(DEBUG-3," violating the capacity constraint")
            else:
                log(DEBUG-3," violating the maximum route cost constraint")
            log(DEBUG-2, " -> EXIT with no changes")
            
        # go to G&M Step 10 -> 
        # no changes, no skipping, route complete
        return route_data, [], [], True
            
    ## G&M Step 11
    D3_nodes = OrderedSet() # the min. dist. from 0 through J,J+1...J+4 to J+5
    D4_nodes = OrderedSet() # the min. dist. /w JJX excluded, KII included
    D6_nodes = OrderedSet() # the min. dist. /w JJX and JII excl., KII incl.
    JJX_in_chain = False
    JII_in_chain = False
            
    # step back so that the first node to lookahead is J
    lookahead_pos = sweep_prev_of_J_pos
    for i in range(5):
        lookahead_pos = _step(lookahead_pos, step_inc, N-2)
        lookahead_node = pos_to_node[lookahead_pos]
        if routed[lookahead_node]:
            continue
        
        D3_nodes.add(lookahead_node)
        if lookahead_node==candidate_node_JJX:
            # inject KII instead of JJX
            D4_nodes.add(to_remove_node_KII)
            D6_nodes.add(to_remove_node_KII)
            JJX_in_chain = True
        elif lookahead_node==candidate_node_JII:            
            D4_nodes.add(lookahead_node)
            JII_in_chain = True
        else:
            D4_nodes.add(lookahead_node)
            D6_nodes.add(lookahead_node) 
    
    # if JJX was not in the sequence J, J+1, ... J+5
    if not JJX_in_chain:
        if __debug__:
            log(DEBUG-2, "G&M Step 11, JJX=n%d not in K(J)..K(J+4)"%
                  candidate_node_JJX)
            log(DEBUG-3, " which consists of nodes %s"%str(list(D3_nodes)))
            log(DEBUG-2, "-> EXIT with no changes")
        # go to G&M Step 10 -> 
        # no changes, no skipping, route complete
        return route_data, [], [], True

    # The chain *end point* J+5
    last_chain_pos = _step(lookahead_pos, step_inc, N-2)
    last_chain_node = pos_to_node[last_chain_pos]
    
    if routed[last_chain_node]:
        last_chain_node = 0
    
    # D3 -> EVALUATE the MINIMUM distance from 0 through J,J+1...J+4 to J+5
    _, D3 = _shortest_path_through_nodes(D, 0, last_chain_node, D3_nodes)
    # D4 -> DETERMINE the MINIMUM distance with JJX excluded, KII included 
    _, D4 = _shortest_path_through_nodes(D, 0, last_chain_node, D4_nodes)

    if not (D1+D3 < D2+D4):
        ## G&M Step 12     
        if __debug__:
            log(DEBUG-2, "G&M Step 12, accept an improving move where "+
               "KII=n%d is removed and JJX=n%d is added" %
               (to_remove_node_KII, candidate_node_JJX))
            log(DEBUG-3, " which forms a route %s (%.2f)"%(str(D2_route),D2))
            log(DEBUG-2, " -> EXIT and continue adding nodes")
            
        ignored_nodes = [to_remove_node_KII]
        if candidate_node_JJX!=node_J:
            ignored_nodes += [node_J]
   
        # go to G&M Step 4 -> 
        # route changed, KII removed and skip current node J, not complete
        return RouteData(D2_route, D2, D2_demand, D2_route_nodes),\
               [candidate_node_JJX], ignored_nodes, False
    
    else:
        ## G&M Step 13
        
        # JII and JJX (checked earlier) should be in K(J)...K(J+4) to continue
        if not JII_in_chain:
            if __debug__:
                if candidate_node_JII is None:
                    log(DEBUG-2, "G&M Step 13, no unrouted nodes left for JII.")
                else:
                    log(DEBUG-2, "G&M Step 13, JII=n%d not in K(J)..K(J+4)"%candidate_node_JII)
                    log(DEBUG-3, " which consists of nodes %s"%str(list(D3_nodes)))
                log(DEBUG-2, "-> EXIT with no changes")
            # go to G&M Step 10 -> no changes, no skipping, route complete
            return route_data, [], [], True
            
        # construct and route to get the modified route cost D2
        D5_route_nodes = D2_route_nodes
        D5_route_nodes.add(candidate_node_JII)
        D5_route,D5 = solve_tsp(D, list(D5_route_nodes))
        D5_demand = D2_demand+d[candidate_node_JII] if C else 0  
        if not ((L and D5-S_EPS<L) and (C and D5_demand-C_EPS<=C)):
            if __debug__:
                log(DEBUG-2, "G&M Step 13, rejecting replacement of KII=n%d with JJX=n%d and JII=n%d"%
                    (to_remove_node_KII, candidate_node_JJX, candidate_node_JII))
                log(DEBUG-3, "  which would have formed a route %s (%.2f)"%(str(D5_route),D5))
                if D5_demand-C_EPS>C:
                    log(DEBUG-3," violating the capacity constraint")
                else:
                    log(DEBUG-3," violating the maximum route cost constraint")
                log(DEBUG-2, "-> EXIT with no changes")
            # go to G&M Step 10 -> no changes, no skipping, route complete 
            return route_data, [], [], True
        
        ## G&M Step 14
        # D6 -> DETERMINE the MINIMUM distance with JJX and JII excluded and
        #  KII ncluded 
        _, D6 = _shortest_path_through_nodes(D, 0, last_chain_node, D6_nodes)
        
        if D1+D3<D5+D6:
            if __debug__:
                log(DEBUG-2, "G&M Step 14, rejecting replacement of KII=n%d with JJX=n%d and JII=n%d"%
                      (to_remove_node_KII, candidate_node_JJX, candidate_node_JII))
                log(DEBUG-3, " which would have formed a route %s (%.2f)"%(str(D5_route),D5))
                log(DEBUG-2,"-> EXIT with no changes")
            # go to G&M Step 10 -> no changes, no skipping, route complete 
            return route_data, [], [], True
    
        ## G&M Step 15  
        if __debug__:
            log(DEBUG-2, "G&M Step 15, accept improving move where "+
                "KII=n%d is removed and JJX=n%d and JII=n%d are added"%
                (to_remove_node_KII, candidate_node_JJX, candidate_node_JII))
            log(DEBUG-3, " which forms a route %s (%.2f)"%(str(D2_route),D2))
            log(DEBUG-2," -> EXIT and continue adding nodes")
  
        ignored_nodes = [to_remove_node_KII]
        if candidate_node_JJX!=node_J and candidate_node_JII!=node_J:
            ignored_nodes += [node_J]
            
        # go to G&M Step 4 -> 
        # route changed, KII removed and skip current node J, not complete
        return RouteData(D5_route, D5, D5_demand, D5_route_nodes),\
               [candidate_node_JJX, candidate_node_JII],\
               ignored_nodes, False

def gillet_miller_init(points, D, d, C, L=None, minimize_K=False, 
                       direction="both", seed_node=BEST_ALTERNATIVE):
    """ This is an implementation of the Giller and Miller (1974) Sweep
    algorithm. The basic scheme is similar to Sweep, but there is an
    online intra-route improvement heuristic, that looks ahead on the sweep and
    checks if including a node from there is expected to improve the resulting
    solution.
    
    * pts, D, d, C, L define the problem.
    
    * direction can be "cw" or "ccw" for clockwise or counter-clockwise sweep
                 or "both" (default) to try and do it in both ways
    * seed_node can be used to set the first customer of the sweep. It also 
                 accepts the same search modes as the sweep.py
                 - BEST_ALTERNATIVE (default), that takes the best result of 
                    all possible sweep start poisitions
                 - CLOSEST_TO_DEPOT, that starts the sweep from the customer
                    closest to the depot
                 - SMALLEST_ANGLE selects the (somewhat arbitary) customer
                    that has the smallest polar coordinate phi.
    """
    if not points:
        raise ValueError("The algorithm requires 2D coordinates for the points")
    
    return sweep_init(points, D, d, C, L, minimize_K,
               direction, seed_node, routing_algo=solve_tsp,
               prepare_callback_datastructures=_pack_datastructures_callback,
               intra_route_improvement=_improvement_callback)

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_gm_algorithm():
    algo_name = "GM74-SwRI"
    algo_desc = "Gillett & Miller (1974) Sweep algorithm with emering "+\
                "route improvement"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if single:
            direction="ccw"
            seed_node=1 #first node    
        else:
            direction="both"
            seed_node=BEST_ALTERNATIVE
            
        return gillet_miller_init(points, D, d, C, L, minimize_K,
                                  direction, seed_node)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_gm_algorithm())
