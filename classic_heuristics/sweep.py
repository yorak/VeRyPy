#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the general Sweep approach
of Gillett and Miller (1974) and Wren & Holliday (1972). This basic 
implementation does not include improvement heuristics, but is written in a
way that it can be extended with online and post-optimization procedures.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has extensive dependencies: a TSP solver (the built in
local search solver can be used), and numpy and scipy for reading and preparing
the problem instance.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
    
from math import pi
from logging import log, DEBUG

import numpy as np

#todo: could use OrderedDict, but it is slow in Python2.7
#  the ordered property is used by wren_holliday (the routes are built in the
#  order the nodes are added during the sweep).
from orderedset import OrderedSet

from util import objf, without_empty_routes, is_better_sol
from routedata import RouteData

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS


__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


CLOSEST_TO_DEPOT = 0
SMALLEST_ANGLE = -1
BEST_ALTERNATIVE = -2

#~defines
PHI = 0
RHO = 1
NODE = 2

def cart2pol(x, y):
    """ helper for that converts x,y to polar coords.
    """
    # numpy sqrt and atan2 are robust and fast
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)

def bisect_angle(angle1, angle2, ratio=0.5, direction=1):
    """ A helper function that gives an angle between angle1 and angle2
    (in rads) based on the bisect ratio. Default gives the angle that bisects
    the sector in the middle.
    """
    while direction*angle2<direction*angle1:
        angle2+=direction*2*pi            
    return angle1-(angle1-angle2)*ratio

def _step(current, inc, max_val):
    current+=inc
    if current>max_val:
        current = 0
    if current<0:
        # reverse direction
        current = max_val
    return current

def _find_earliest_previous(current, previous, inc, max_val):
    min_d = float('inf')
    min_p = None
    for p in previous:
        if inc==1 and p>=current:
            d = p-current
        elif inc==1 and p<current:    
            d = p+max_val-current
        elif inc==-1 and p<=current:
            d = current-p
        elif inc==-1 and p>current:    
            d = current+max_val-p
        if d < min_d:
            min_d = d
            min_p = p
    return min_p

def _ensure_L_feasible(D, d, C, L, current_route,
                       L_feasible_sweep_pos, step_inc, max_sweep_idx,
                       routed, sweep_pos_to_node_idx, routing_algo):
    """ This makes sure the current route in satisfies the
        route length/cost/duration constraint (L). We follow Gillett & Miller 
        (1974) where nodes are removed from the route set until the route
        cost satisfies the L constraint.
        
        Potentially modifies current_route and routed"""
    pos_to_remove = L_feasible_sweep_pos
    while True:
        if current_route.cost-S_EPS<=L:
            break
        else:
            pos_to_remove = _step(pos_to_remove, -step_inc, max_sweep_idx)   
            node_to_remove = sweep_pos_to_node_idx(pos_to_remove)
            if node_to_remove in current_route.node_set:
                current_route.node_set.remove(node_to_remove)       
                routed[node_to_remove] = False
                tsp_sol, tsp_cost = routing_algo(D, list(current_route.node_set))
                # update current route
                current_route.route = tsp_sol
                current_route.cost = tsp_cost
                if C:
                    current_route.demand-=d[node_to_remove]
                
                # also the sweep goes back one step
                L_feasible_sweep_pos = pos_to_remove   
                
                if __debug__:
                    log(DEBUG-2, "L constraint was violated, removed"+
                      " n%d from the route set"%(node_to_remove))
                    log(DEBUG-3, "Got TSP solution %s (%.2f)" %
                                 (str(tsp_sol),tsp_cost)) 
            
    return L_feasible_sweep_pos

def do_one_sweep(N, D, d, C, L, routing_algo,
                  sweep, start, step_inc,
                  generate_alternative_first_routes = False,
                  intra_route_callback=None, inter_route_callback=None,
                  callback_data=None):
    """ This function does one full circle starting from start index of the 
    sweep and proceeds to add the nodes one by one in the direction indicated 
    by the step_inc parameter. A new route is started if either of the
    constraints constraints C and L (L can be None) are violated.

    If generate_alternative_first_routes is set, every possible route until 
    C OR L constraints are generated and then the sweep terminates.
    This is used by some constuction heuristics to build potential routes.
    
    The generated routes are TSP optimized using routing_algo and may be 
    further optimized with intra and inter route  improvement callbacks.""" 
    
    # make sure all before the start node have larger angles than the start node
    # and that they are of range [0-2*pi]
    if intra_route_callback or inter_route_callback:
        sweep_phis = np.copy( sweep[0] )
        #sweep_phis = -step_inc/abs(step_inc)*sweep_phis
        sweep_phis-=(pi+sweep_phis[start]) # the first is at -pi
        sweep_phis[start-step_inc::-step_inc]+=2*pi
        sweep_phis = -1*sweep_phis
    #node_idx_to_sweep_pos = lambda n: int(sweep[3,n-1])
    sweep_pos_to_node_idx = lambda idx: int(sweep[2,idx])
    max_sweep_idx = len(sweep[2])-1
    
    if __debug__:
        int_sweep = list(sweep[2].astype(int))
        log(DEBUG-1, "Sweep node order %s\n"%str(int_sweep[start:]+int_sweep[:start]))

    # Routes
    routes = [] 
    routed = [False]*N
    routed[0] = True
    routed_cnt = 0
    total_to_route = len(sweep[0])
    
    # Emerging route
    current_route = RouteData([0])
    current_route.node_set = OrderedSet([0])
    current_route_cost_upper_bound = 0
    route_complete = False
    
    # Backlog
    # improvement heuristics may leave some nodes unrouted, collect them here
    #  and insert them on next (or last) routes
    blocked_nodes = set()
    
    # THE MAIN SWEEP LOOP
    # iterate until a full sweep is done and the backlog is empty    
    sweep_pos = start
    sweep_node = sweep_pos_to_node_idx(sweep_pos)  
    while True:
        if __debug__:
            if sweep_node:
                prev_pos = _step(sweep_pos, -step_inc, max_sweep_idx)
                next_pos = _step(sweep_pos, step_inc, max_sweep_idx)
                prev_ray = bisect_angle(sweep[0][prev_pos], sweep[0][sweep_pos], direction=step_inc) 
                next_ray = bisect_angle(sweep[0][sweep_pos], sweep[0][next_pos], direction=step_inc) 
                log(DEBUG-2, "Considering n%d between rays %.2f, %.2f" %
                             (sweep_node,prev_ray,next_ray))
        
        # append route to the list of routes if
        # 1. all nodes have been added to the routes OR
        # 2. the route demand would be higher than the vehicle capacity 
        #    (check possible L constraint later) OR
        # 3. if there is only L constraint, check if adding sweep node would
        #    break the maximum route distance constraint.
        if not route_complete and C:
            would_break_C_ctr = current_route.demand+d[sweep_node]-C_EPS>C
            route_complete = would_break_C_ctr            
        if not route_complete and L and current_route_cost_upper_bound>L:
            tsp_sol, tsp_cost = \
                routing_algo(D, list(current_route.node_set)+[sweep_node]) 
            current_route_cost_upper_bound = tsp_cost
            would_break_L_ctr = tsp_cost-S_EPS>L
            route_complete = would_break_L_ctr
            
        route_interesting = generate_alternative_first_routes and \
                            len(current_route.node_set)>1
        if route_complete or route_interesting:  
            tsp_sol, tsp_cost = routing_algo(D, list(current_route.node_set))
            current_route.route = tsp_sol
            current_route.cost = tsp_cost
            current_route_cost_upper_bound = current_route.cost
            
            if __debug__:
                if route_complete:
                    log(DEBUG-2, "Route %s full."%str(list(current_route.node_set)))
                log(DEBUG-3, "Got TSP solution %s (%.2f)"%(str(tsp_sol),tsp_cost))
                     
            if L:
                L_feasible_pos = _ensure_L_feasible(D, d, C, L, 
                    current_route, sweep_pos, step_inc, max_sweep_idx, routed, 
                    sweep_pos_to_node_idx, routing_algo)
                current_route_cost_upper_bound = current_route.cost
                
                if L_feasible_pos!=sweep_pos:
                    sweep_pos = L_feasible_pos  
                    if sweep_node!=None:
                        sweep_node = sweep_pos_to_node_idx(sweep_pos)
                    # route is full if L constraint was violated
                    route_complete = True
                    route_interesting = False
                
            if sweep_node!=None and intra_route_callback:                
                # Do intra route improvement. It may include new unrouted
                #  customers to the route or remove existing ones!
                current_route, included_nodes, ignored_nodes, route_complete = \
                    intra_route_callback( current_route, callback_data,
                                          sweep[1], sweep_phis,
                                          sweep_pos, step_inc, routed )
                current_route_cost_upper_bound = current_route.cost
                
                if __debug__:
                    log(DEBUG-1, "Improved the route to %s (%.2f)"%(current_route.route, current_route.cost))
                
                # Make sure all routed nodes are routed.
                for rn in included_nodes:
                    routed[rn]=True  
                # Make sure skipped and removed nodes are added in later.
                for lon in ignored_nodes:
                    routed[lon] = False
                    if __debug__:
                        log(DEBUG-2, "Block n%d for now"%lon)
                    # avoid trying to add it on next iter.
                    blocked_nodes.add( lon )
                    if lon==sweep_node:
                        # skip adding current node
                        sweep_node = None
                
             # If route is complete, store it and start a new one   
            if route_complete:
                if __debug__:
                    log(DEBUG-1, "Route completed %s\n"%str(current_route.route))
                    #log(DEBUG-1, "with backlog of %s\n"%str(node_backlog))

                # Check if we have all routed, and can exit the main sweep loop                
                route_customer_cnt = len(current_route.node_set)
                if route_customer_cnt>1:
                    routed_cnt+=route_customer_cnt-1
                    routes.append( current_route )
                    
                if routed_cnt>=total_to_route:
                    if __debug__:
                        log(DEBUG-1,"Last route completed %s\n"%str(current_route.route))
                    break # SWEEP
                
                current_route = RouteData([0])
                current_route.node_set = OrderedSet([0])
                current_route_cost_upper_bound = 0.0
                route_complete = False
                
                # New route, allow customers in backlog
                blocked_nodes = set()
                
                if generate_alternative_first_routes:
                    break
            
            # if the route demand is higher than the lower bound, record it
            elif route_interesting:
                if __debug__:
                    log(DEBUG-1, "Route recorded %s\n"%str(current_route.route))
                # store a deep copy
                routes.append( RouteData(list(current_route.route),
                                         current_route.cost,
                                         current_route.demand,
                                         OrderedSet(current_route.node_set)) )
                
        # Route improvement can route nodes, so ensure this was not routed.
        if (sweep_node is not None) and (not routed[sweep_node]):
            current_route.node_set.add(sweep_node)
            routed[sweep_node] = True
            if C:
                current_route.demand+=d[sweep_node]
            if L:
                prev_node = 0
                if len(current_route.route)>2:
                    prev_node = current_route.route[-1] \
                                if current_route.route[-1]!=0\
                                else current_route.route[-2]
                # calculate and add the delta of just appending the swept node
                ub_delta = -D[prev_node, 0]+D[prev_node, sweep_node]+D[sweep_node,0]
                current_route_cost_upper_bound+=ub_delta
            
            if __debug__:
                log(DEBUG-2, "Added n%d to the route set"%sweep_node)
        
        if __debug__:
            log(DEBUG-3, "Step to a next sweep node from position %d (n%d) with %s blocked."%
                (sweep_pos, sweep_pos_to_node_idx(sweep_pos), list(blocked_nodes)))
        start_stepping_from = sweep_pos
        while True:
            sweep_pos = _step(sweep_pos, step_inc, max_sweep_idx)
            sweep_node = sweep_pos_to_node_idx(sweep_pos)
            
            if (not routed[sweep_node]) and (sweep_node not in blocked_nodes):
                break # found an unrouted node continue with it
                
            if sweep_pos == start_stepping_from:
                # We checked, and it seems there is no unrouted non-blocked
                # nodes left -> start a new route, reset blocked and try again.
                sweep_node = None
                route_complete = True
                blocked_nodes = set()
                break
        
    # (optional) improvement phase         
    if inter_route_callback:
        routes = inter_route_callback(routes, callback_data)
        
    return routes

def get_sweep_from_polar_coordinates(rhos,phis):
    N = len(rhos)
    # stack the arrays, so that we can sort them (along different dimensions)
    customer_phirhos = np.stack( (phis[1:],rhos[1:],np.arange(1,N)) )
    sweep_node_order = np.argsort(customer_phirhos[0])
    sweep = customer_phirhos[:,sweep_node_order]
    return sweep

def get_sweep_from_cartesian_coordinates(pts):
    """Convert cartesian coordinates of the customer to polar coordinates
    centered at the depot. Returns a sweep which is a stack of 3 vectors:
        
     * phi angles for all customers 
     * rho distances all customers 
     * customer indexes (as float)
    
    Also the phi angles are returned (including the depot)
    Note that the sweep is a size of N-1, as it does not include the depot.
    """
    
    np_pts = pts if isinstance(pts, np.ndarray) else np.asarray(pts)
    depot_x, depot_y = pts[0]
    rhos, phis = cart2pol(np_pts[:,0]-depot_x, np_pts[:,1]-depot_y)
    return get_sweep_from_polar_coordinates(rhos,phis)


def sweep_init(coordinates, D, d, C, L=None, minimize_K=False,
               direction="both", seed_node=BEST_ALTERNATIVE,
               routing_algo=None, **callbacks):
    """
    This algorithm was proposed in Wren (1971) and in Wren & Holliday
    (1972). Sweep was also proposed in Gillett and Miller (1974) who
    gave the algorithm its name. The proposed variants differ in on how many
    starting locations (seed) for the sweep are considered: four in Wren &
    Holliday (1972) and all possible in both directions in Gillett and Miller
    (1974). Also, the improvement procedures differ. The version in this file
    is basebones as as it does not include any route improvement heuristics.
    For implementations of Gillett and Miller (1974) or  Wren & Holliday (1972)
    algorithms, please see their Python files (gillet_miller_sweep.py and
    wren_holliday_sweep.py).
    
    The basic principle of the Sweep algorithm is simple: The algorithm assumes
    that the distances of the CVRP are symmetric, and, furthermore, that the
    points are located on a plane. The catresian coordinates of these points
    in relation to the depot are converted to polar coordinates (rho, phi), and
    then sorted by phi. Starting from an arbitary node (in this implementation
    the default is the one closest to the depot) create a new route and add 
    next  adjecent unrouted node according to their angular coordinate. Repeat 
    as long as the capacity is not exceeded. When this happens, start a new 
    route and repeat the procedure until all nodes are routed. Finally, the 
    routes can optionally be optimized using a TSP algorithm. 
       
    Note that the algorithm gives different results depending on the direction 
    the nodes are inserted. The direction parameter can be "cw" for clockwise 
    insertion order and "ccw" for counterclockwise. As the algorithm is quite
    fast, it is recommended to run it in both directions.
    
    Please note that the actual implementation of the sweep procedure is in the
     do_one_sweep function.
    
    * coordinates can be either 
        a) a list/array of cartesian coordinates (x,y)
        b) a lists/arrays (3) of polar coodinates WITH node indexes (i.e. a 
            numpy stack of phi,rho,idx)
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/cost/duration.
    * direction is either "cw" or "ccw" depending on the direction the nodes
       are to be processed
    * seed_node is optional parameter that specifies how the first node of the 
       sweep is determned. This can be one of CLOSEST_TO_DEPOT (0),
       SMALLEST_ANGLE (-1), BEST_ALTERNATIVE (-2), which tries every possible 
       staring id, or a positive integer explicitly specifying the node id to
       start from. Also, a list of indexes can be given. These are explicit
       sweep indexes and it is adviseable to give also the sweep parameter.
    
    Wren, A. (1971), "Computers in Transport Planning and Operation", Ian 
      Allan, London.
    Wren, A., and Holliday, A. (1972), "Computer scheduling of vehicles from
      one or more depots to a number of delivery points", Operations Research
      Quarterly 23, 333-344.
    Gillett, B., and Miller, L., (1974). "A heuristic algorithm for the vehicle
      dispatch problem". Operations Research 22, 340-349.
    """
    
    N = len(D)
    if len(coordinates[0])==2:
        sweep = get_sweep_from_cartesian_coordinates(coordinates)
    elif len(coordinates)==3 and (len(coordinates[0])==len(coordinates[1])==len(coordinates[2])):
        # it is not necessarily to sweep to contain all nodes in D and d
        sweep = coordinates
    else:
        raise ValueError("The coordinates need to be (x,y) or (phi,rho,node_index,sweep_index_for_node-1). Not "+str(coordinates))
        
    
    ## specify the direction
    if direction == "ccw":
        step_incs = [1]    
    elif direction == "cw":
        step_incs = [-1]        
    elif direction == "both":
        step_incs = [1,-1]        
    else:
        raise ValueError("""Only "cw", "ccw", and "both" are valid values for the direction parameter""")

    ## specify where to start
    if seed_node==CLOSEST_TO_DEPOT:
        starts = [np.argmin(sweep[1])]
    elif seed_node==SMALLEST_ANGLE:
        starts = [0]
    elif seed_node==BEST_ALTERNATIVE:
        starts = list(range(0,N-1))
    elif type(seed_node) is int:
        # we interpret it as a node idx
        starts = [np.where(sweep[2]==abs(seed_node)%N)[0][0]]
    elif type(seed_node) is list:
        # we interpret it as a node idx
        starts = seed_node
     
    ## Make sure there is a valid route improvement method
    if routing_algo is None:
        # Default generates the route from the list of nodes in the order they
        #  were swept. Assume that depot (0) is the first of node_set.
        routing_algo = lambda D, node_set: (list(node_set)+[0],
                                            objf(list(node_set)+[0],D))
        
    ## for exteding Sweep with improvement heuristics
    callback_data = None
    intra_route_callback = None
    inter_route_callback = None
    if 'prepare_callback_datastructures' in callbacks:
        pcds_callback = callbacks['prepare_callback_datastructures']
        callback_data = pcds_callback(D,d,C,L,sweep)
    if 'intra_route_improvement' in callbacks:        
        intra_route_callback = callbacks['intra_route_improvement']
    if 'inter_route_improvement' in callbacks:        
        inter_route_callback = callbacks['inter_route_improvement']
        
    ## Do the search with the parameter specified above
    best_sol = None
    best_f = None  
    best_K = None
    
    try:
        for step_inc in step_incs:
            for start in starts:
                if __debug__:
                    log(DEBUG, "\nDo a sweep from position %d (n%d) by steps of %d"%
                                 (start,sweep[2][start],step_inc))
                
                ## This does one sweep from one start location to one direction
                routes = do_one_sweep(N, D, d, C, L, routing_algo,
                                           sweep, start, step_inc,
                                           False,
                                           intra_route_callback,
                                           inter_route_callback,
                                           callback_data)            
                    
                sol = [n for rd in routes for n in rd.route[:-1]]+[0]
                # LS of the callbacks may cause empty routes
                sol = without_empty_routes(sol)
                sol_f = objf( sol, D )   
                sol_K = sol.count(0)-1
        
                if __debug__:
                    log(DEBUG, "Previous sweep produced solution %s (%.2f)\n\n" %
                                 (str(sol),sol_f))
                    
                if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
                    best_sol = sol
                    best_f = sol_f
                    best_K = sol_K
    except KeyboardInterrupt: # or SIGINT
        raise KeyboardInterrupt(best_sol)
        
    return best_sol


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_swp_algorithm():
    algo_name = "Sweep"
    algo_desc = "Sweep algorithm without route improvement heuristics"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        seed_search = SMALLEST_ANGLE if single else BEST_ALTERNATIVE
        direction = "cw" if single else "both"
        return sweep_init(points, D, d, C, L, minimize_K,
                          direction=direction, seed_node=seed_search)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_swp_algorithm())    
    
    
