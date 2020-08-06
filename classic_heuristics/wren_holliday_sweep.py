#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Wren & Holliday (1972)
Sweep algorithm with a local search post-optimization phase.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs.It has moderate dependencies: the sweep procedure from
sweep.py, the OrderedSet module, and a TSP solver (the built-in one based on 
local search can be used). As usual, numpy and scipy are required for reading
and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

import numpy as np
from math import pi
from itertools import permutations
from logging import log, DEBUG

from classic_heuristics.sweep import sweep_init, get_sweep_from_polar_coordinates,\
                  cart2pol, BEST_ALTERNATIVE
from routedata import RouteData
         
from local_search import LSOPT               
from local_search.intra_route_operators import do_2opt_move, do_relocate_move
from local_search.inter_route_operators import do_1point_move, do_chain_move,\
                                               do_insert_move,do_redistribute_move

from util import objf, without_empty_routes, is_better_sol

from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


#shorthand
FAS = LSOPT.FIRST_ACCEPT

# defines (do not change)
BEST_OF_FOUR = -99
LEAST_DENSE = -98
LEAST_DEMAND = -97

#~consts for readabiliy
NEG_RANGE = 0
POS_RANGE = 1
RANGE_MIN = 0
RANGE_MAX = 1

def _get_route_phi_range(route_phis):
    # to take into account the discontinuity of the angles, 
    #  the ranges are stored separately for both halves 
    
    if len(route_phis)==0:
        return None
        
    # negative, positive
    half_ranges = [[0,-pi],[pi,0]]
    
    sign = None
    prev_phi = route_phis[0]
    for phi in route_phis:
        new_sign = phi/abs(phi) if phi!=0 else sign
        
        # route edge crosses the valid range of [-pi,pi] or origo
        if (sign is not None) and (new_sign != sign):
            # wraps
            if abs(phi-prev_phi)>pi:
               half_ranges[NEG_RANGE][RANGE_MIN]=-pi
               half_ranges[POS_RANGE][RANGE_MAX]=pi
            # jumps over origo
            else:
               half_ranges[NEG_RANGE][RANGE_MAX]=0
               half_ranges[POS_RANGE][RANGE_MIN]=0
        if phi==0 or phi<0:
            if half_ranges[NEG_RANGE][RANGE_MAX]!=0:
                half_ranges[NEG_RANGE][RANGE_MAX] = \
                    max(half_ranges[NEG_RANGE][RANGE_MAX],phi)
            if half_ranges[NEG_RANGE][RANGE_MIN]!=-pi:
                half_ranges[NEG_RANGE][RANGE_MIN] = \
                    min(half_ranges[NEG_RANGE][RANGE_MIN],phi)
        if phi==0 or phi>0:
            if half_ranges[POS_RANGE][RANGE_MAX]!=pi:
                half_ranges[POS_RANGE][RANGE_MAX] = \
                    max(half_ranges[POS_RANGE][RANGE_MAX],phi)
            if half_ranges[POS_RANGE][RANGE_MIN]!=0:
                half_ranges[POS_RANGE][RANGE_MIN] = \
                    min(half_ranges[POS_RANGE][RANGE_MIN],phi)
        if phi!=0:
            sign = new_sign
        prev_phi = phi
            
    return half_ranges

def _remove_empty_in_place(routes):
    routes[:] = [r for r in routes if not r.is_empty()]

def inspect_heuristic(routes,D,C,d,L):
    improvement_found = False
    for rd in routes:
        if rd.is_empty():
            continue
        
        delta = 0.0
        while delta is not None:
            new_route, delta = do_2opt_move(rd.route, D, FAS)
            if delta is not None:
                rd.route = new_route
                rd.cost+=delta
                improvement_found = True
    return improvement_found 
                
def single_heuristic(routes,D,C,d,L):
    # TODO: this does not exactly correspond to the Wren and Holliday
    #   as the intra and inter route operations are different local 
    #   search moves, and, therefore, all options for a node are 
    #   not considered in the same order.
    #  Therefore, it is a possible source of replication problems.

    improvement_found = False
    for rd1_idx, rd1 in enumerate(routes):
        if rd1.is_empty():
            continue
        
        route_improvement_found = True
        while route_improvement_found:
            route_improvement_found = False
            
            # move on same route
            new_route, delta = do_relocate_move(rd1.route, D, FAS)
            if delta is not None:
                rd1.route = new_route
                rd1.cost += delta
                route_improvement_found = True
                improvement_found = True
                
            # move between routes
            for rd2_idx, rd2 in enumerate(routes):
                if rd1_idx==rd2_idx:
                    continue
                new_rd1, new_rd2, delta =\
                    do_1point_move(rd1, rd2, D, d, C, L, FAS)
                if delta is not None:
                    routes[rd1_idx] = new_rd1
                    routes[rd2_idx] = new_rd2       
                    rd1 = new_rd1
                    route_improvement_found = True
                    improvement_found = True
               
    return improvement_found
    

def pair_heuristic(routes,D,C,d,L):
    #todo: modify so that rd3==rd1 is possible (node swap)
    # do_chain_move is very expansive, use demands to shrink the search space
    np_d = np.array(d)
    min_d = np.min(np_d)
    
    if C:
        smallest_demands = {}
        for rd in routes:
            if rd.is_empty():
                smallest_demands[rd] = 0
            else:
                smallest_demands[rd] = np.min(np_d[rd.route[1:-1]])    
        
    improvement_found = False
    for rd1_idx,rd2_idx,rd3_idx in permutations(range(len(routes)),3):
        rd1 = routes[rd1_idx]
        rd2 = routes[rd2_idx]
        rd3 = routes[rd3_idx]
        
        if rd1.is_empty() or rd2.is_empty():
            continue
        
        if C:
            # there has to be room for smallest demand AND the smallest demand 
            #  in route2. 
            r3_room_left = C+C_EPS-smallest_demands[rd3]
            if (r3_room_left < min_d or
                r3_room_left < smallest_demands[rd2]):
                continue
            
        # the actual chain / "pair" operation
        new_rd1,new_rd2,new_rd3,delta = do_chain_move(rd1,rd2,rd3, D,d,C,L,FAS)

        if delta is not None:
            routes[rd1_idx] = new_rd1
            routes[rd2_idx] = new_rd2
            routes[rd3_idx] = new_rd3
            
            if C:
                smallest_demands[new_rd1] = 0 if len(new_rd1.route)==2 else\
                                            np.min(np_d[new_rd1.route[1:-1]])  
                smallest_demands[new_rd2] = 0 if len(new_rd2.route)==2 else\
                                            np.min(np_d[new_rd2.route[1:-1]])  
                smallest_demands[new_rd3] = 0 if len(new_rd3.route)==2 else\
                                            np.min(np_d[new_rd3.route[1:-1]])
            improvement_found = True
       
    return improvement_found

def delete_heuristic(routes, D,C,d,L):
    """
    If several routes can be deleted, apply the operation that improves the
    solution most (or worsens it it least).
    """
    
    improvement_found = False    
    best_insertion_delta = float("inf")
    best_insertion_routes = None

    for rd1_idx, rd1 in enumerate(routes):
        if rd1.is_empty():
            continue
        
        # do not allow empty routes
        other_routes = [r for r in routes if (not r.is_empty() and r!=rd1)]
        
        result = do_redistribute_move(rd1, other_routes, D,d,C,L,
                               best_delta = best_insertion_delta,
                               recombination_level=0)
        
        delta = result[-1]
        if (delta is not None) and (delta<best_insertion_delta):
            improvement_found = True
            best_insertion_delta = result[-1]
            best_insertion_routes = result[:-1]
    
    if improvement_found:
        routes[:] = best_insertion_routes
    return improvement_found

def complain_heuristic(routes,unrouted_nodes, #input / outputs
                       D,C,d,L):
    """ This procedure tries to add unrouted nodes that were e.g. omitted from
    the initial routes on the existing routes by removing one node from them
    and reinserting it to another route.
    
    It is possible that not all nodes can be routed. Then a list of nodes that
    could not be inserted, or that were removed to make a better insertions, is
    returned. Empty list means that unrouted_nodes were successfully inserted.
    """
    was_assigned = set()
    was_omitted = set()
    improvement_found = False
    
    for node in unrouted_nodes:
        best_replace_delta = None
        best_replace_move = None
        node_routed = False
        for rd1_idx in range(len(routes)):
            # unpack route, current cost, and current demand
            rd1 = routes[rd1_idx]
            if rd1.is_empty():
                continue
            route1, r1l, r1d, _ = rd1
                  
            for i in range(1,len(route1)-1):
                replace_after = route1[i-1]
                was_replaced = route1[i]
                replace_before = route1[i+1]
                
                
                replace_d = 0
                if C:
                    replace_d = r1d-d[was_replaced]+d[node]
                    if replace_d-C_EPS>C:
                        continue
                
                replace_delta = +D[replace_after, node]\
                        +D[node, replace_before]\
                        -D[replace_after, was_replaced ]\
                        -D[was_replaced , replace_before]
                if L and r1l+replace_delta-S_EPS>L:
                    continue
                
                if (best_replace_delta is None) or\
                   (replace_delta<best_replace_delta):
                    best_replace_delta = replace_delta
                    best_replace_move = (rd1_idx, i)
                      
                # the unrouted insertion is feasible, now fit the removed back
                for rd2_idx in range(len(routes)):
                    if rd1_idx==rd2_idx:
                        continue
                    rd2 = routes[rd2_idx]
                                        
                    # the insertion position does not matter, but it must fit
                    #  and be C and L feasible
                    _, new_rd2, insert_delta = do_insert_move(
                                                  was_replaced,rd2, D,d,C,L, 
                                                  LSOPT.FIRST_ACCEPT)
                    
                    if insert_delta is not None:
                        # operation successful, accept the move
                        routes[rd2_idx] = new_rd2
                        routes[rd1_idx] = RouteData(
                            route1[:i]+[node]+route1[i+1:],
                            r1l+replace_delta,
                            replace_d)
                        node_routed = True
                        improvement_found = True
                        was_assigned.add(node)
                        break
                
                if node_routed:
                    break
            if node_routed:
                break
                
        
        if not node_routed:
            if best_replace_delta is not None:
                (best_rd_idx, i) = best_replace_move
                # place that of highest priority on the route
                best_route, best_rl, best_rd, _ = routes[best_rd_idx]
                to_remove = best_route[i]
                
                # Wren and Holliday propose that priority is given to the customer
                #  requesting the greater load.
                if C and d[node]>d[to_remove]:
                    routes[best_rd_idx] = RouteData(
                        best_route[:i]+[node]+best_route[i+1:],
                        best_rl+best_replace_delta,
                        best_rd-d[to_remove]+d[node])
                    
                    was_assigned.add(node)
                    was_omitted.add(to_remove)
                    improvement_found = True      
              
    # update the situation of unrouted_nodes               
    unrouted_nodes -= was_assigned
    unrouted_nodes |= was_omitted
    return improvement_found
        
def combine_heuristic(routes, D,C,d,L):
    raise NotImplementedError("Was not implemented in Wren and Holliday 1972")
    
def disentangle_heuristic(routes,sweep,node_phis,D,C,d,L):
    # find a overlapping / "tanged" pair
    improvement_found = False
    entangled = True

    route_phi_ranges = []
    for rd in routes:
        route = rd.route
        if rd.is_empty():
            route_phi_ranges.append(None)
        else:
            r_phis = node_phis[route[1:-1]]
            r_phi_range = _get_route_phi_range(r_phis)
            route_phi_ranges.append(r_phi_range)
    
    while (entangled):
        entangled = False
        for rd1_idx,rd2_idx in permutations(range(len(routes)),2):  
            # unpack route, current cost, and current demand
            route1, r1l, r1d, _ = routes[rd1_idx]
            route2, r2l, r2d, _ = routes[rd2_idx]
            r1_phi_range = route_phi_ranges[rd1_idx]
            r2_phi_range = route_phi_ranges[rd2_idx]
            
            # the route may have become empty
            if (r1_phi_range is None) or (r2_phi_range is None):
                continue

            # ranges overlap            
            overlap = ((r1_phi_range[0][0] < r2_phi_range[0][1] and
                        r2_phi_range[0][0] < r1_phi_range[0][1]) or\
                       (r1_phi_range[1][0] < r2_phi_range[1][1] and\
                        r2_phi_range[1][0] < r1_phi_range[1][1]))
            if overlap:
                et_nodes = route1[:-1]+route2[1:-1]
                et_nodes.sort()
                
                # Check for backwards compatibility: numpy.isin was introduced
                #  in 1.13.0. If it not availabe (as is the case e.g. in Ubuntu
                #  16.04 LTS), use the equivalent list comprehension instead.
                if hasattr(np, 'isin'):
                    sweep_mask = np.isin( sweep[2], et_nodes, assume_unique=True )
                else:
                    sweep_mask = np.array([item in et_nodes for item in sweep[2]])
                    
                et_sweep = sweep[:,sweep_mask]
                
                
                over_wrap_phi = (et_sweep[0][0]+2*pi)-et_sweep[0][-1]
                et_seed_node = [ np.argmax( np.ediff1d(et_sweep[0], to_end=over_wrap_phi )) ]
                
                reconstructed_et_sol = sweep_init( et_sweep, D, d, C, L,
                    seed_node=et_seed_node, direction="cw", routing_algo=None)
                    
                # translate nodes back to master problem node indices
                #reconstructed_et_sol = [et_nodes[n] for n in reconstructed_et_sol]
                   
                # construct an array of RouteData objects for heuristics
                et_routes = RouteData.from_solution(reconstructed_et_sol, D, d)
                omitted = set()
                for extra_route in et_routes[2:]:
                    omitted|=set(extra_route.route[1:-1])
                et_routes = et_routes[:2]
                
                if __debug__:
                    log(DEBUG-1, "DISENTANGLE omitted %s"%str(list(omitted)))
                    _log_after_ls_op("DISENTANGLE/SWEEP", True, et_routes, D)
                    
                refined = True
                while refined:
                    refined = False
                    refined|=inspect_heuristic(et_routes,D,C,d,L)
                    if __debug__:
                        _log_after_ls_op("DISENTANGLE/INSPECT", refined, et_routes, D)
                    refined|=single_heuristic(et_routes,D,C,d,L)
                    if __debug__:
                        _log_after_ls_op("DISENTANGLE/SINGLE", refined, et_routes, D)
                    # there were omitted nodes, try to complain them in
                    if omitted:
                        refined|=complain_heuristic(et_routes,omitted,D,C,d,L)
                        if __debug__:
                            _log_after_ls_op("DISENTANGLE/COMPLAIN", refined, et_routes, D)
                    
                    
                # all nodes routed and disentangling improves -> accept this
                if len(et_routes)==1:
                    # keep an empty route (for now) to avoid indexing errors
                    et_routes.append( RouteData([0,0], 0.0, 0.0) )
                    disentange_improves = True
                else:
                    disentange_improves = et_routes[0].cost+et_routes[1].cost+S_EPS < r1l+r2l
                    
                if len(omitted)==0 and disentange_improves:
                    routes[rd1_idx] = et_routes[0]
                    routes[rd2_idx] = et_routes[1]
                    
                    # The routes were updated, update also the ranges
                    r1_phis = node_phis[et_routes[0].route[1:-1]]
                    r2_phis = node_phis[et_routes[1].route[1:-1]]
                    r1_phi_range = _get_route_phi_range(r1_phis)
                    r2_phi_range = _get_route_phi_range(r2_phis)
                    route_phi_ranges[rd1_idx] = r1_phi_range
                    route_phi_ranges[rd2_idx] = r2_phi_range
                
                if __debug__:
                    if len(omitted)>0:
                        log(DEBUG-1, "DISENTANGLE rejected ( due to omitted %s)"%
                                     str(omitted))
                    elif disentange_improves:
                        log(DEBUG-1, "DISENTANGLE rejected ( due to not improving %.2f vs. %.2f )"%
                                     (et_routes[0].cost+et_routes[1].cost, r1l+r2l))
                    else:
                        log(DEBUG-1, "DISENTANGLE accepted ( due to inserting omitted and improving %.2f vs. %.2f )"%
                                     (et_routes[0].cost+et_routes[1].cost, r1l+r2l))
                        
    return improvement_found

def _log_after_ls_op(name, changed, routes, D):
    log(DEBUG-1, "After %s (%s):"%(name, str(changed)))
    sol_f = 0.0
    for ri, rd in enumerate(routes):
        log(DEBUG-1, "Route #%d : %s"%(ri, str(list(rd.route))))
        sol_f += objf(rd.route, D)
    log(DEBUG-1, "Cost %.2f\n"%sol_f)
    
def wren_holliday_init(points, D, d, C, L=None, minimize_K=False,
                       seed_node=BEST_OF_FOUR, direction='both',
                       full_convergence = True):
    """ This implements the Wren and Holliday improvement heuristic. The
    initial solution is generated using the generic sweep procedure of
    `sweep.py`, and the improvement procedure works as specified in 
    Wren & Holliday (1972) Fig 1 (Flowchart of program). Basically, the 
    refining processes, also known as local search improvement /
    post-optimization phase, repeadedly applies local search operators to 
    improve the initial solutions and returns the best one.
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route cost/duration/length.
    
    * seed_node sets how many different seed nodes are tried for the Sweep
       based initial solution generation. If LEAST_DENSE, the sweep is started
       from the direction from the depot that has lowest customer density. The
       same applies to BEST_OF_FOUR (default), but also 3 other directions
       spaced by ~90deg are considered. Also a complete (but costly) generation
       of all possible initial solutions with BEST_ALTERNATIVE can be used. 
    * direction can be 'ccw' for counter-clockwise Sweep, 'cw' (default) for
       clockwise or 'both' for trying both.
    * full_convergence determines if the improvement is stopped after the
       "delete" operation fails the first time (False) or if the local search
       continues until no operation is capable of finding an improving more 
       (True, default).
       
    Returns the solution.
    
    Wren, A., & Holliday, A. (1972). Computer scheduling of vehicles from one
    or more depots to a number of delivery points. Journal of the Operational
    Research Society, 23(3), 333-344.
    """                      
    
    if not points:
        raise ValueError("The algorithm requires 2D coordinates for the points")
    N = len(D)
    
    # Calculate the sweep coordinates
    # this is 99% same as _get_sweep..., but we need to get also
    #  the node_phis for later use, so duplicated some code here.
    np_pts = points if isinstance(points, np.ndarray) else np.asarray(points)
    depot_x, depot_y = points[0]
    node_rhos, node_phis = cart2pol(np_pts[:,0]-depot_x, np_pts[:,1]-depot_y)
    sweep = get_sweep_from_polar_coordinates(node_rhos,node_phis)
    sweep_phis = sweep[0]        
    
    directions = ['cw', 'ccw'] if direction=='both' else [direction]
    sweeps = []
    
    for cur_dir in directions:
        if seed_node==BEST_OF_FOUR or seed_node==LEAST_DENSE:
            # Wren & Holliday method of selecting starting customers for the
            #  Sweep initialization from the "least dense direction".
            # Turn the polar coordinates so that the smallest angles are to the
            #  direction of the smallest density, weighted by their demands.
            o = np.array([depot_x, depot_y])
            if d:
                weighted_xy_from_origin = np.multiply(
                    np_pts[1:]-o,
                    np.transpose(np.array([d[1:],d[1:]])))
                avgx, avgy = np.average(weighted_xy_from_origin, axis = 0)
            else:
                avgx, avgy = np.average(np_pts[1:]-o, axis = 0)
            avg_rho, avg_phi = cart2pol(np.array([avgx]),
                                        np.array([avgy]))
            # to range [-pi, pi]
            if avg_phi>0:
                least_dense_phi = avg_phi[0]-pi
            else:
                least_dense_phi = avg_phi[0]+pi
            
            # evenly spaced by pi/2
            for i in range(4):
                angle_tgt = least_dense_phi+i*pi/2
                if angle_tgt>pi:
                    angle_tgt-=pi*2
                    
                # take the first satisfying the condition
                start_from_here = np.argmax(sweep_phis>angle_tgt)
                start_node = start_from_here-1 if cur_dir=="cw" \
                                                 else start_from_here
                if start_node==-1:
                   start_node = N-2
                sweeps.append( (start_node, cur_dir) )
                
                if seed_node==LEAST_DENSE:
                    break # do not take the ones at 90 deg intervals
            
        elif seed_node==BEST_ALTERNATIVE:
            nodes = list(range(N-1))
            sweeps.extend( zip( nodes, [cur_dir]*len(nodes) ) )
        elif type(seed_node)==int:
            sweeps.append( [seed_node, cur_dir] )
        
    ## PHASE 1 : "Generate ... initial solutions and choose the best"    
    initial_sols  = []
    
    try:
        for start_node, cur_dir in sweeps:
            isol = sweep_init(sweep, D, d, C, L, seed_node=[start_node],
                              direction=cur_dir, routing_algo=None)
            initial_sols.append( isol )
            
    except KeyboardInterrupt as e: # or SIGINT
        # if interrupted on initial sol gen, return the best of those
        if len(e.args)>0 and type(e.args[0]) is list:
            initial_sols.append( e.args[0] )
        if not initial_sols:
            raise e
        else:
            best_isol, best_if, best_iK = None, float('inf'), float('inf')
            for isol in initial_sols:
                isol_f = objf(isol, D)
                isol_K = isol.count(0)-1
                if is_better_sol(best_if, best_iK, isol_f, isol_K, minimize_K):
                    best_isol = isol
                    best_if = isol_f
                    best_iK = isol_K
            raise KeyboardInterrupt(best_isol)
    
    
    best_sol, best_f, best_K = None, float('inf'), float('inf')
    interrupted = False
    for sol in initial_sols:
        
        # Construct an array of RouteData objects for local search improvement
        #  heuristics to use
        routes = RouteData.from_solution(sol, D, d)

        if __debug__: 
            log(DEBUG, "Improving solution %s (%.2f)"%
                       (sol, objf(sol,D)))
            _log_after_ls_op("SWEEP(S)", False, routes, D)
                
        ## PHASE 2 : Improvement phase (see Wren & Holliday 1974, Figure 1)
        
        # Variables storing the state
        converging = False
        deleted = True
        omitted_nodes = set()
        prev_Q_point_sol_f = None    
        prev_iteration_sol_f = None
        changed = False
        
        try:
            while True:
                _remove_empty_in_place(routes)
                if not minimize_K:
                    # +1 empty route can allow local search to find an improvement
                    routes.append( RouteData() )            
                            
                changed = False
                
                ## "INSPECT, SINGLE" ##
                # run 2opt on each route to remove any crossing edges
                inspect_improved = inspect_heuristic(routes,D,C,d,L)
                if inspect_improved and minimize_K:
                    _remove_empty_in_place(routes)
                changed |= inspect_improved
                if __debug__: _log_after_ls_op("INSPECT", changed, routes, D)
                
                # move a node to a better position on the route or other routes
                single_improved = single_heuristic(routes,D,C,d,L)
                if single_improved and minimize_K:
                    _remove_empty_in_place(routes)
                changed |= single_improved
                if __debug__: _log_after_ls_op("SINGLE", changed, routes, D)
                
                ## "Are customers omitted?" ##
                omitted_were_assigned = False
                if omitted_nodes:
                    inserted_nodes = set()
                    ## "Take omitted ... in order and try to fit into existing routes" ##
                    for node in sorted(list(omitted_nodes)):
                        for rdi, rd in enumerate(routes):
                            _, new_rd, delta = do_insert_move(node, rd, D,d,C,L,
                                                              LSOPT.BEST_ACCEPT)
                            if delta is not None:
                                routes[rdi] = new_rd
                                inserted_nodes.add(node)
                                omitted_were_assigned = True
                    omitted_nodes-=inserted_nodes
                
                    if omitted_were_assigned and minimize_K:
                        _remove_empty_in_place(routes)
                    changed |= omitted_were_assigned
                    if __debug__: _log_after_ls_op("INSERT", changed, routes, D)
                
                ## "Are customers omitted still?" ##
                if omitted_nodes:
                    omitted_were_assigned |= complain_heuristic(routes,omitted_nodes, D,C,d,L)
                    if omitted_were_assigned and minimize_K:
                        _remove_empty_in_place(routes)
                    changed |= omitted_were_assigned                
                    if __debug__: _log_after_ls_op("COMPLAIN", changed, routes, D)
                    
                sol_f = 0
                for rd in routes:
                    sol_f+=rd.cost
                    
                ## Q-point : "Has distance been reduced by more that 5% OR
                # has a previously omitted customer been assigned?" ##
                if (prev_Q_point_sol_f is None) or\
                   (sol_f<prev_Q_point_sol_f*0.95) or \
                   omitted_were_assigned:
                    prev_Q_point_sol_f = sol_f
                    converging = False
                    continue
                else:
                    prev_Q_point_sol_f = sol_f
                    converging = True
                
                ## "Is problem small?" -> PAIR ##
                if len(D)<=80:
                    pair_improved = pair_heuristic(routes,D,C,d,L)
                    if pair_improved and minimize_K:
                        _remove_empty_in_place(routes)
                    changed |= pair_improved 
                    if __debug__: _log_after_ls_op("PAIR", changed, routes, D)
                    
                ## "Is deleted true?" ##
                if deleted:
                    # "DELETE" -> "Was delete succesful?" ##
                    deleted = delete_heuristic(routes, D,C,d,L)
                    if deleted and minimize_K:
                        _remove_empty_in_place(routes)
                    changed |= deleted
                    if __debug__: _log_after_ls_op("DELETE", changed, routes, D)
                    
                ## DISENTANGLE ##
                disentangle_improved = disentangle_heuristic(routes,sweep,node_phis,D,C,d,L)
                if disentangle_improved and minimize_K:
                    _remove_empty_in_place(routes)
                changed |= disentangle_improved
                if __debug__: _log_after_ls_op("DISENTANGLE", changed, routes, D)
                
                     
                ## "Has situation changed in interation?" ##
                solution_changed_between_iterations = True
                if prev_iteration_sol_f:
                    if prev_iteration_sol_f == sol_f:
                        solution_changed_between_iterations = False
                prev_iteration_sol_f = sol_f
                
                if converging and (
                  (full_convergence and not changed) or
                  (not full_convergence and not deleted) or 
                  (not solution_changed_between_iterations)):
                    ## STOP ##
                    break
                
        except KeyboardInterrupt:
            interrupted = True
            
        # return the optimized solutios
        sol = [0]+[n for r in routes for n in r.route[1:]]
        # LS may cause empty routes
        sol = without_empty_routes(sol)
        sol_K = sol.count(0)-1
        sol_f = objf(sol, D)
        
        if __debug__:
            log(DEBUG, "Improved solution %s (%.2f)"% (sol,sol_f))
        
        if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
            best_sol = sol
            best_f = sol_f
            best_K = sol_K
            
        if interrupted:
            raise KeyboardInterrupt(best_sol)
            
    return best_sol


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_wh_algorithm():
    algo_name = "WH72-SwLS"
    algo_desc = "Wren and Holliday (1972) Sweep heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        seed_node = 1 if single else BEST_OF_FOUR
        return wren_holliday_init(points, D, d, C, L, minimize_K, seed_node)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_wh_algorithm())
