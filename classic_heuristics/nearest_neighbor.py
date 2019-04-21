#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the basic nearest neighbour
heuristic. Both parallel and sequential algorithms are provided. A description
of the algorithm can be found e.g. in van Breedam (1994), but it seems that the
 idea of adapting the approach for CVRP originates from Tyagi (1967).

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has really minimal dependencies: only numpy and scipy 
are needed for reading and preparing the problem instance.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from logging import log, DEBUG
from operator import itemgetter
from collections import deque
from util import objf, without_empty_routes, is_better_sol, routes2sol

from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


class _PeekQueue:
    """ This is a helper data strucure, that allows peeking into a list. It is
    used to calculate the capacity total with ne next nearest neighbor 
    without consuming it from the iterator.
    
    TODO: replace with deque (instead of peek, just pop and push)
    it will be faster, simpler, and reduce SLOC.
    
    TODO: While at it implement minimize_K so that multiple customers can 
    be peeked. Keep track of the smallest d (if C is set) and smallest inc 
    (min dist from 0 to j + j to k where k is closest node to j) and do not
    continue if improvement is not possible. Adds complexity!
    """
    
    def __init__(self, l):
        self.posleft = -1
        self.posright = 0
        self.l = list(l)
        
    def __len__(self):
        return len(self.l)-(self.posleft+1)+(self.posright)
    
    def __getitem__(self, idx):
        if idx>=len(self):
            raise IndexError
        return self.l[self.posleft+1+idx]
        
    def peekleft(self):
        if len(self)==0:
            raise IndexError
        return self.l[self.posleft+1]

    def popleft(self):
        if len(self)==0:
            raise IndexError
        self.posleft+=1
        return self.l[self.posleft]    
    
    def peekright(self):
        if len(self)==0:
            raise IndexError
        return self.l[self.posright-1]

    def popright(self):
        if len(self)==0:
            raise IndexError
        self.posright-=1
        return self.l[self.posright]

def get_seed_node(seed_mode, D, node_nearest_neighbors, served):
    N = len(node_nearest_neighbors)
    seed_node = None
    while True:
        if seed_mode=='farthest':
            seed_node  = node_nearest_neighbors[0].popleft()[0]
        elif seed_mode=='closest':
            seed_node  = node_nearest_neighbors[0].popright()[0]
        elif seed_mode=='nearest':
            nearest_distance = float('inf')
            for i in range(1,N):
                if  served[i]:
                    continue
                elif seed_node==None:
                    seed_node = i # in case there is no pair                                
                while not served[i] and len(node_nearest_neighbors[i])>0:
                    nearest = node_nearest_neighbors[i].peekleft()[0]
                    if served[nearest]:
                        #print("pop", nearest)
                        node_nearest_neighbors[i].popleft()
                    else:
                        break
                if D[i,nearest] < nearest_distance:
                    seed_node  = i
                    nearest_distance = D[i,nearest]
        else:
            raise ValueError("Only 'farthest', 'closest', and "+
                             "'nearest' route initialization "+
                             "methods are supported")
        if not served[seed_node]:
            break 
    return seed_node 
            
def nearest_neighbor_init(D, d, C, L=None, 
                          emerging_route_count=1,
                          initialize_routes_with="farthest",
                          add_only_to_end=False,
                          forbidden_nodes=None,
                          route_improvement_callback=None):
    """ A greedy nearest neighbor algorithm for solving symmetric CVRPs. 
    The general idea of the heuristic has been proposed e.g. in (Tyagi 1967).
    The nearest unrouted node is added to a route until a constraint (C or L)
    is violated and then a new route (taking the nearest unrouted node from 
    the depot) is initialized. This is repeated until all nodes are routed.

    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.
    
    * emerging_route_count sets how many routes are constructed in parallel
    * initialize_routes_with can be "farthest" (default) or "closest", which
       selects if the unrouted node farthest from the depot or closest to the
       depot is used to initialize the emerging route. Can also be "nearest", 
       and then the route is initialized with the two unrouted nearest 
       neighbors.
    * add_only_to_end is same to van Breedam (1994, 2002) "places to add stops"
       parameter. It sets if the customers are added to the end or at either
       end of the emering route (default).
    
    The forbidden_nodes and route_improvement_callback can be used to extend 
    nearest_neighbor_init. forbidden_nodes is a list of nodes that are not
    considered for insertion and route_improvement_callback is an route
    improvement heuristic with the signature:
    route_improvement_callback(route_data, D, d, C, L,
                               served, node_nearest_neighbors), where
    Route data is a 4-tuple (or similar), with 
      (route, route_demand, route_cost, None)
    
    TODO: implement other route seeding methods (now only "closest" and "farthest")
    TODO: find a lit. ref. for sequential and parallel
     esp. check what the parallel lit. algo does when a route is full...
     
    Tyagi, M. (1968). A practical method for the truck dispatching problem. 
    Journal of the Operations Research Society of Japan, 10:76-92
    """
    
    # build the nearest neighbor lists
    N = len(D)
    node_nearest_neighbors = [None]*N
    for i in range(N):
        node_nearest_neighbors[i] = _PeekQueue( sorted(enumerate(D[i][:]),
                                                        key=itemgetter(1)) )
        node_nearest_neighbors[i].popleft() # pop reference to self
        
    # bookkeeping on the nodes that we have already served
    served = [False]*N
    served[0]=True
    
    if forbidden_nodes:
        for fn in forbidden_nodes:
            served[fn] = True
    
    sol = [0]
    route_nodes = [None]*emerging_route_count
    prev_added = [None]*emerging_route_count
    route_demands = [0.0]*emerging_route_count
    route_costs = [0]*emerging_route_count
    
    try:
        k = 0        
        while True:

            if not route_nodes[k]:
                # Needs to initialize a new emerging route
                seed_node =  get_seed_node(initialize_routes_with, D,
                                            node_nearest_neighbors, served)
                    
                # build a route around the seed node
                route_nodes[k] = deque()
                route_nodes[k].append(seed_node)
                prev_added[k] = seed_node
                if C: route_demands[k] = d[seed_node]
                if L: route_costs[k] = 0.0
                served[seed_node]=True
                                
                if __debug__:
                    route = [0]+list(route_nodes[k])+[0]
                    log(DEBUG-1, "Initialize route #%d with n%d, resulting to %s (%.2f)"%
                        (k,seed_node, str(route), objf(route, D)))
                
            else:        
                # add as a new nonserved node after the previous one
                first_node = route_nodes[k][0] 
                last_node = route_nodes[k][-1] 
                prev_node = prev_added[k]
                
                nn_node = None
                while True:
                    nn_node = node_nearest_neighbors[prev_node].peekleft()[0]
                    if served[nn_node]:
                        #print("pop", nn_node)
                        node_nearest_neighbors[prev_node].popleft()
                    else:
                        break
                    
                if __debug__:
                    log(DEBUG-2, "Consider to extend route #%d with n%d"%(k,nn_node))
                
                delta = D[0,first_node]+D[last_node,nn_node]+D[nn_node,0]
                add_to_end = True
                if not add_only_to_end:
                    delta_start =  D[0,nn_node]+D[nn_node,first_node]+D[last_node,0]
                    if delta_start<delta:
                        delta = delta_start
                        add_to_end = False        
                    
                # check constraints
                constraints_violated = \
                  (C and route_demands[k]+d[nn_node]-C_EPS>C) or \
                  (L and route_costs[k]+delta-S_EPS>L)
                
                if constraints_violated:
                    current_route = [0]+list(route_nodes[k])+[0]
                    if route_improvement_callback:
                        current_cost = D[0,first_node]+route_costs[k]+D[last_node,0]
                        
                        # note: route_improvement_callback may make some served
                        #  and move from from served to not-served!
                        improved_route = route_improvement_callback(
                            current_route, route_demands[k], current_cost,
                            D, d, C, L, served, node_nearest_neighbors)
                        
                        if __debug__:
                            if improved_route != current_route:
                                log(DEBUG-2, "Route improved from %s (%.2f) to %s (%.2f)"%
                                    (str(current_route), objf(current_route, D),
                                     str(improved_route), objf(improved_route, D)))
                                
                        current_route = improved_route # leave out the last 0
                        
                    if __debug__:
                        log(DEBUG-1, "Constraint violated, storing route #%d as %s (%.2f)"%
                            (k,str(current_route), objf(current_route,D)))
                    
                    sol.extend( current_route[1:] )
                    route_nodes[k] = None
                    # Did not add a node (as there was no capacity left)
                    #  do not increment k, as it is still the turn of this
                    #  route to get a new node.
                        
                else:
                    if C:
                        route_demands[k] += d[nn_node]
                    if add_to_end:
                        route_nodes[k].append(nn_node)
                        if L:
                            route_costs[k] += D[last_node,nn_node]
                    else: # add_to_front
                        route_nodes[k].appendleft(nn_node)
                        if L:
                            route_costs[k] += D[nn_node,first_node]
                    prev_added[k] = nn_node
                    served[nn_node]=True
                    
                    if __debug__:
                        route = [0]+list(route_nodes[k])+[0]
                        log(DEBUG-1, "Extending route #%d with n%d, resulting to %s (%.2f)"%
                            (k, nn_node, str(route), objf(route,D)))
                    
                    # a node was added to the route. "pop" the peeked value
                    # off the iterator and continue to the next route. 
                    node_nearest_neighbors[prev_node].popleft()
                    k+=1
            
            # Parallel version builts K routes in prarallel
            if k==emerging_route_count:
                k = 0
            
    except (IndexError, KeyboardInterrupt) as e:
        interrupted = type(e) is KeyboardInterrupt
            
        #TODO: modify the code so that route_improvement_callback is called
        # and search continued if a node is removed. Remember to avoid infinite
        # improvement loop if nodes are re-marked as unrouted.
        # 
        # The algorithm ran out of nearest neighbors, which means that all
        #  nodes are served. Close the remaining routes with a return to depot.
        for rcustomers in route_nodes:
            if rcustomers:
                sol.extend(rcustomers)
                sol.append(0)

                if __debug__:
                    if not interrupted:
                        current_route = [0]+list(rcustomers)+[0]
                        croute_f = objf(current_route,D)
                        log(DEBUG-1, "Constraint violated, storing route "
                            "#%d as %s (%.2f)"%(k,str(current_route),croute_f))                    
        
        if interrupted:
            if sol:
                unroutedr = [[n] for n in range(N)
                             if (n not in sol)]
                if sol and len(sol)>1:
                    interrupted_sol = sol[:-1]+routes2sol( unroutedr )
                else:
                    interrupted_sol = routes2sol( unroutedr )
                raise KeyboardInterrupt(interrupted_sol)
                
        elif len(sol)==1:
            sol = None

    return sol     
    
 
# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_snn_algorithm():
    algo_name = "vB95-SNN"
    algo_desc = "van Breedam (1994) Sequential Nearest Neighbor construction "+\
                "heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if minimize_K:
            raise NotImplementedError("Nearest neighbor algorithm does "+
                                          " not support minimizing the number"+
                                          " of vehicles")
        return nearest_neighbor_init(D,d,C,L, emerging_route_count=1)
    return (algo_name, algo_desc, call_init)

def get_pnn_algorithm(emerging_route_count="auto"):
    algo_name = "vB95-PNN"
    algo_desc = "Parallel Nearest Neighbor construction heuristic"
    if emerging_route_count=="auto":
        def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
            if minimize_K:
                # todo: remove this when supprot (see TODO notes in algo desc)
                raise NotImplementedError("Nearest neighbor algorithm does "+
                                          " not support minimizing the number"+
                                          " of vehicles")
              
            sol_snn = nearest_neighbor_init(D, d, C, L, emerging_route_count=1)            
            if single:
                return sol_snn
            
            auto_route_count = sol_snn.count(0)-1
            
            # NN is so fast we can try with several K and take the best
            best_sol = sol_snn
            best_f = objf(sol_snn,D)
            best_K = auto_route_count
            for k in range(2,auto_route_count+1):
                sol = nearest_neighbor_init(D, d, C, L, emerging_route_count=k)
                sol = without_empty_routes(sol)
                sol_f = objf(sol,D)
                sol_K = sol.count(0)-1
                
                if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
                    best_sol = sol
                    best_f = sol_f
                    best_K = sol_K
                    
            return best_sol
    elif emerging_route_count>1:
        def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
            if minimize_K:
                # todo: remove this when supprot (see TODO notes in algo desc)
                raise NotImplementedError("Nearest neighbor algorithm does "+
                                          " not support minimizing the number"+
                                          " of vehicles")
            return nearest_neighbor_init(D, d, C, L, emerging_route_count)
    else:
        raise ValueError("Not a valid emerging_route_count value "+
                         "(%s) for parallel algorithm"%
                         str(emerging_route_count))
    return (algo_name, algo_desc, call_init)
        
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_pnn_algorithm())
