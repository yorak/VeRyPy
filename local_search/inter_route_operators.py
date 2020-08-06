# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and implements inter route (i.e. between two or more routes)
local search improvement heuristics such as one and two point move etc. All 
operators assume we start from a feasible solution. Also, all functions 
implementing the operations have the following signature:

do_Y_move(route1_data, route2_data, ... , D,C,d,L, strategy, best_delta), where

* routeZ_data is a RouteData object, the data members are
        1. route as a list of node indices. Must begin and end to depot (0).
        2. current route cost
        3. current route demand, can be None is C = None
* D is the numpy-compatibe distance matrix as with intra route heuristics.
* C is the optional capacity constraint (can be None). If set, also give d.
* d is a list of customer demands, and of the depot it is d[0]=0,
    The parameter must be set if  C is given.
* L is the maximun route length/duration/cost constraint
* strategy is either FIRST_ACCEPT (default) or BEST_ACCEPT. 
    First accept returns a modified route as soon as the first improvement is  
    encoutered. Best accept tries all possible combinations and returns the
    best one.
* best_delta is the required level of change in the in route cost. It can be
   used to  set the upper bound (requirement) for the improvement. Usually it 
   is None, which sets the level to 0.0. In that case only improving deltas are
   accepted. If set to (large) positive value, best worsening move can also be
   returned.
   
All inter route improvement operators return the new improved routes and the
improvement (delta) as (n_r+1)-tuple, where n_r is the number of routes given 
to the operator or (None,...,None) if no improvements were found."""
###############################################################################

#TODO the efficiency of these implementations could be greatly improved using
# numba, doubly linked lists for routes, and in place editing of routes.
#  - use doubly linked list (dllist) as the data structure for the routes
#  - do not copy routes, make changes in place
#  - use numba (and numba compatible dllist) or even cython 

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from itertools import product, permutations

from routedata import RouteData
from local_search import LSOPT, ROUTE_ORDER_SENSITIVE_OPERATORS
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

# Create a decorator to mark routes sensitive to the route order
def routeordersensitive(f):
    ROUTE_ORDER_SENSITIVE_OPERATORS.add(f)
    return f
    
def do_2optstar_move(route1_data, route2_data, D, d=None, 
                      C=None, L=None, # constraints
                      strategy=LSOPT.FIRST_ACCEPT,
                      best_delta = None):
    
    """ 2-opt* inter-route local search operation for the symmetric distances D
    Remove 2 edges from different routes and check if swapping the edge halves
    (in two different ways) would yield an improvement, while making sure the
    move does not violate constraints.
    """

    # use 2-opt
    if route1_data==route2_data:
        raise ValueError("Use do_2opt_move to find intra route moves")

    # make sure we have the aux data for constant time feasibility checks   
    if not isinstance(route1_data, RouteData):
        route1_data = RouteData(*route1_data)
    if not isinstance(route2_data, RouteData):
        route2_data = RouteData(*route2_data)
    if not route1_data.aux_data_updated:
        route1_data.update_auxiliary_data(D,d)
    if not route2_data.aux_data_updated:
        route2_data.update_auxiliary_data(D,d)

    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    #print("REMOVEME: 2opt* on %s %s"%(list(route1_data.route), list(route2_data.route)) )
    
    
    for i in range(0,len(route1_data.route)-1):
        for j in range(0,len(route2_data.route)-1):
            a = route1_data.route[i]
            b = route1_data.route[i+1]
            c = route2_data.route[j]
            d = route2_data.route[j+1]
            
            #print("REMOVEME: attempt remove %d-%d and %d-%d"%(a,b,c,d) )
            
            # a->c b->d
            #       __________
            #      /          \
            # 0->-a   b-<-0-<-c   d->-0
            #         \___________/  
            #
            delta = D[a,c] + D[b,d] \
                     -D[a,b]-D[c,d]
            if delta+S_EPS<best_delta:
                # is an improving move, check feasibility 
                constraint_violated = False
                r1_new_demand = None
                r2_new_demand = None
                if C:
                    r1_new_demand = route1_data.fwd_d[i]+route2_data.fwd_d[j]
                    r2_new_demand = route1_data.rwd_d[i+1]+route2_data.rwd_d[j+1]
                    if r1_new_demand-C_EPS>C or r2_new_demand-C_EPS>C:
                        constraint_violated = True
                if not constraint_violated and L and (    
                  route1_data.fwd_l[i]+D[a,c]+route2_data.fwd_l[j]-S_EPS>L or
                  route1_data.rwd_l[i+1]+D[b,d]+route2_data.rwd_l[j+1]-S_EPS>L):
                    constraint_violated = True
                
                if not constraint_violated:
                    # store segments, deltas, and demands
                    best_move = (((None,i+1,1), (j,None,-1),
                                  D[a,c]-D[a,b], r1_new_demand),
                                 ((None,i,-1), (j+1,None,1),
                                  D[b,d]-D[c,d], r2_new_demand))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # j loop
            
            # a->d c->b
            #       ______________
            #      /              \
            # 0->-a   b-<-0-<-c   d->-0
            #         \______/  
            #
            delta = D[a,d] + D[b,c] \
                     -D[a,b]-D[c,d]
            if delta+S_EPS<best_delta:
                # is an improving move, check feasibility 
                constraint_violated = False
                r1_new_demand = None
                r2_new_demand = None
                if C:
                    r1_new_demand = route1_data.fwd_d[i]+route2_data.rwd_d[j+1]
                    r2_new_demand = route1_data.rwd_d[i+1]+route2_data.fwd_d[j]
                    if r1_new_demand-C_EPS>C or r2_new_demand-C_EPS>C:
                        constraint_violated = True
                if not constraint_violated and L and (
                  route1_data.fwd_l[i]+D[a,d]+route2_data.rwd_l[j+1]-S_EPS>L or
                  route1_data.rwd_l[i+1]+D[b,c]+route2_data.fwd_l[j]-S_EPS>L):
                    constraint_violated = True
                    
                if not constraint_violated:
                    # store segments, deltas, and demands
                    best_move = (((None,i+1,1), (j+1,None,1),
                                  D[a,d]-D[a,b], r1_new_demand),
                                 ((None,i,-1), (j,None,-1),
                                  D[b,c]-D[c,d], r2_new_demand))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # j loop
                     
        if accept_move:
            break # i loop
                
    if best_move:
        # unpack the move
        ((r1_sgm1, r2_sgm1, r1_delta, r1_new_demand),
         (r1_sgm2, r2_sgm2, r2_delta, r2_new_demand)) = best_move
         
        return (
            # route 1
            RouteData(
                route1_data.route[r1_sgm1[0]:r1_sgm1[1]:r1_sgm1[2]]+\
                route2_data.route[r2_sgm1[0]:r2_sgm1[1]:r2_sgm1[2]],
                route1_data.cost+r1_delta,
                r1_new_demand),
            # route 2
            RouteData(
                route1_data.route[r1_sgm2[0]:r1_sgm2[1]:r1_sgm2[2]]+\
                route2_data.route[r2_sgm2[0]:r2_sgm2[1]:r2_sgm2[2]],
                route1_data.cost+r1_delta,
                r1_new_demand),
            # delta
            best_delta)
 
    return None, None, None

# These define all possible 3opt moves. Describes which segment ends to connect
#  This makes implementing the do_3optstar_move (below) operator simpler.

#     
#    __0  1__      route 1
#   /        \        
# n0---2  3--n0    route 2
#   \__    __/      
#      4  5        route 3  
#
MOVES_3OPTSTAR_3ROUTES = (
    #((0, 1), (2, 3), (4, 5)), 
    #2-opt
    ((0, 1), (2, 4), (3, 5)), 
    ((0, 1), (2, 5), (4, 3)), 
    ((0, 4), (2, 3), (1, 5)),
    ((0, 5), (2, 3), (4, 1)),
    ((0, 2), (1, 3), (4, 5)), 
    ((0, 3), (2, 1), (4, 5)), 
    #3-opt
    ((0, 2), (5, 3), (4, 1)),
    ((0, 2), (4, 3), (1, 5)),
    ((0, 3), (2, 5), (4, 1)),
    ((0, 3), (2, 4), (1, 5)),
    ((0, 4), (2, 1), (3, 5)),
    ((0, 4), (1, 3), (2, 5)),
    ((0, 5), (2, 1), (4, 3)),
    ((0, 5), (1, 3), (4, 2)),
)

#
#    _0  1_2  3_     route 1
#   /           \  
# n0            n0  
#   \____    ___/  
#        4  5        route 2
#
# All moves have 3 segments on the second route
#  (connected by the two first edges)
MOVES_3OPTSTAR_2ROUTES = (
    #((0, 1), (2, 3), (4, 5)), 
    
    #2-opt
    ((3, 5), (0, 1), (2, 4)), # note, route 2 first
    ((4, 3), (0, 1), (2, 5)), # note, route 2 first
    ((0, 4), (3, 2), (1, 5)), 
    ((0, 5), (4, 1), (2, 3)), 
    ((4, 5), (0, 2), (1, 3)), # note, route 2 first
    #((0, 3), (2, 1), (4, 5)), would form a loop
    #3-opt
    ((3, 5), (0, 2), (1, 4)), # note, route 2 first
    ((4, 3), (0, 2), (1, 5)), # note, route 2 first
    ((0, 3), (4, 1), (2, 5)),
    ((0, 3), (4, 2), (1, 5)),
    #((0, 4), (2, 1), (3, 5)), would form a loop
    ((0, 4), (3, 1), (2, 5)),
    #((0, 5), (2, 1), (4, 3)), would form a loop
    ((0, 5), (4, 2), (1, 3)),
)

@routeordersensitive
def do_3optstar_2route_move(route1_data, route2_data, D, demands=None, 
                      C=None, L=None, # constraints
                      strategy=LSOPT.FIRST_ACCEPT,
                      best_delta = None):

    if route1_data==route2_data:
        raise ValueError("Use do_3opt_move to find intra route moves")
        
    # make sure we have the aux data for constant time feasibility checks   
    if not isinstance(route1_data, RouteData):
        route1_data = RouteData(*route1_data)
    if not isinstance(route2_data, RouteData):
        route2_data = RouteData(*route2_data)
    if not route1_data.aux_data_updated:
        route1_data.update_auxiliary_data(D,demands)
    if not route2_data.aux_data_updated:
        route2_data.update_auxiliary_data(D,demands)
        
    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    # segment end nodes and cumulative d and l at those nodes is specified 
    # by the indices i,j,k. In the two route case this means:
    #
    #    _0/i  1_2/j  3_     route 1
    #   /               \  
    # n0                n0  
    #   \____      _____/  
    #        4/k  5        route 2
    #
    end_n = [0]*6
    cum_d = [0]*6
    cum_l = [0]*6
    
    for i in range(0,len(route1_data.route)-1):
        end_n[0] = route1_data.route[i]
        end_n[1] = route1_data.route[i+1]

        # make it so that i<j
        for j in range(i+1 ,len(route1_data.route)-1):
            
            # the edge endpoints
            end_n[2] = route2_data.route[j]
            end_n[3] = route2_data.route[j+1]
            if C:
                cum_d[0] = route1_data.fwd_d[i]
                cum_d[1] = route1_data.rwd_d[i+1]\
                          -route1_data.rwd_d[j+1]
                cum_d[2] = route1_data.fwd_d[j]\
                           -route1_data.fwd_d[i]
                cum_d[3] = route1_data.rwd_d[j+1]
            if L:
                cum_l[0] = route1_data.fwd_l[i]
                cum_l[1] = route1_data.rwd_l[i+1]\
                           -route1_data.rwd_l[j]
                cum_l[2] = route1_data.fwd_l[j]\
                           -route1_data.fwd_l[i+1]
                cum_l[3] = route1_data.rwd_l[j+1]
                
            
            for k in range(0,len(route2_data.route)-1):
                
                # the edge endpoints
                end_n[4] = route2_data.route[k]
                end_n[5] = route2_data.route[k+1]
                if C:
                    cum_d[2] = route2_data.fwd_d[k]
                    cum_d[3] = route2_data.rwd_d[k+1]
                if L:
                    cum_l[2] = route2_data.fwd_l[k]
                    cum_l[3] = route2_data.rwd_l[k+1]
                
                removed_weights = D[end_n[0],end_n[1]]+\
                                  D[end_n[2],end_n[3]]+\
                                  D[end_n[4],end_n[5]]
      
                for e1,e2,e3 in MOVES_3OPTSTAR_2ROUTES:
                    e1_wt = D[end_n[e1[0]],end_n[e1[1]]]
                    e2_wt = D[end_n[e2[0]],end_n[e2[1]]]
                    e3_wt = D[end_n[e3[0]],end_n[e3[1]]]
                
                    delta = e1_wt+e2_wt+e3_wt-removed_weights
                    
                    if ((delta+S_EPS<best_delta) and
                        (not C or (cum_d[e1[0]]+cum_d[e1[1]]-C_EPS<C and
                         cum_d[e2[0]]+cum_d[e2[1]]+cum_d[e3[1]]-C_EPS<C)) and
                        (not L or ( cum_l[e1[0]]+cum_l[e1[1]]+e1_wt-S_EPS<L and
                          cum_l[e2[0]]+cum_l[e2[1]]+cum_l[e3[1]]+e2_wt+e3_wt-S_EPS<L))):
                        best_move = ((i,j,k),(e1,e2,e3),delta)
                        best_delta = delta  
                        
                        if strategy==LSOPT.FIRST_ACCEPT:
                            accept_move = True
                            break # move loop
                if accept_move:
                    break # k loop
            if accept_move:
                break # j loop
        if accept_move:
            break # i loop
    
    
    if best_move:
        raise NotImplementedError("Not yet completely implemented")
        
        # TODO: 
        # 1. unpack the move
        # 2. construct the modified routes
        # 3. return the modified routes and delta
          
        # unpack the move
        #(best_ijk,best_edges,best_delta) = best_move
        #best_e1, best_e2, best_e3 = best_edges
    
#         new_route1 = 
#        
#            if edge[0]%2 and edge[1]%2:
#                new_route = (
#                        routes[r1_idx].route[None:r1_i :-1]+
#                        routes[r2_idx].route[r2_i+1:None:1])
#                new_cost = ( D[b, d]+
#                        routes[r1_idx].rwd_l[r1_i+1]+
#                        routes[r2_idx].rwd_l[r2_i+1])
#                new_demand = (
#                        routes[r1_idx].rwd_d[r1_i+1]+
#                        routes[r2_idx].rwd_d[r2_i+1])
#            # combine the head of route1 with the tail of route2
#            elif not edge[0]%2 and edge[1]%2:
#                new_route = (
#                        routes[r1_idx].route[None:r1_i+1:1]+
#                        routes[r2_idx].route[r2_i+1:None:1])
#                new_cost = (D[a, d],
#                        routes[r1_idx].fwd_l[r1_i]+
#                        routes[r2_idx].rwd_l[r2_i+1])
#                new_demand = (
#                        routes[r1_idx].fwd_d[r1_i]+
#                        routes[r2_idx].rwd_d[r2_i+1])
        
def do_3optstar_3route_move(route1_data, route2_data, route3_data, D, demands=None, 
                      C=None, L=None, # constraints
                      strategy=LSOPT.FIRST_ACCEPT,
                      best_delta = None):
    
    """ 3-opt* inter-route local search operation for the symmetric distances D
    Remove 3 edges from different routes and check if reconnecting them in any
    configuration would improve the solution, while also making sure the
    move does not violate constraints.
    
    For two route checks, the first and second route can be same i.e. 
    route1_data==route2_data. However, if route2_data==route3_data or 
    route1_data==route3_data a ValueError is raised.
    """

    if route1_data==route2_data or\
       route2_data==route3_data or\
       route1_data==route3_data:
        raise ValueError("Use do_3opt_move to find intra route moves")
   
   # make sure we have the aux data for constant time feasibility checks   
    if not isinstance(route1_data, RouteData):
        route1_data = RouteData(*route1_data)
    if not isinstance(route2_data, RouteData):
        route2_data = RouteData(*route2_data)
    if not isinstance(route3_data, RouteData):
        route3_data = RouteData(*route3_data)
    if not route1_data.aux_data_updated:
        route1_data.update_auxiliary_data(D,demands)
    if not route2_data.aux_data_updated:
        route2_data.update_auxiliary_data(D,demands)
    if not route3_data.aux_data_updated:
        route3_data.update_auxiliary_data(D,demands)

    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    # segment end nodes and cumulative d and l at those nodes
    # Note the indexing of segment end nodes here:
    #     
    #    __0/i  1__      route 1
    #   /          \        
    # n0---2/j  3--n0    route 2
    #   \__      __/      
    #      4/k  5        route 3  
    #
    end_n = [0]*6
    cum_d = [0]*6
    cum_l = [0]*6
    
    for i in range(0,len(route1_data.route)-1):
        end_n[0] = route1_data.route[i]
        end_n[1] = route1_data.route[i+1]

        # make it so that i<j if route1==route2
        for j in range(0 ,len(route2_data.route)-1):
            
            # the edge endpoints
            end_n[2] = route2_data.route[j]
            end_n[3] = route2_data.route[j+1]
            if C:
                cum_d[0] = route1_data.fwd_d[i]
                cum_d[1] = route1_data.rwd_d[i+1]
                cum_d[2] = route2_data.fwd_d[j]
                cum_d[3] = route2_data.rwd_d[j+1]
            if L:
                cum_l[0] = route1_data.fwd_l[i]
                cum_l[1] = route1_data.rwd_l[i+1]
                cum_l[2] = route2_data.fwd_l[j]
                cum_l[3] = route2_data.rwd_l[j+1]
            
            for k in range(0,len(route3_data.route)-1):
                
                # the edge endpoints
                end_n[4] = route3_data.route[k]
                end_n[5] = route3_data.route[k+1]
                if C:
                    cum_d[4] = route3_data.fwd_d[k]
                    cum_d[5] = route3_data.rwd_d[k+1]
                if L:
                    cum_l[4] = route3_data.fwd_l[k]
                    cum_l[5] = route3_data.rwd_l[k+1]
                
                removed_weights = D[end_n[0],end_n[1]]+\
                                  D[end_n[2],end_n[3]]+\
                                  D[end_n[4],end_n[5]]
      
                for e1,e2,e3 in MOVES_3OPTSTAR_3ROUTES:
                    e1_wt = D[end_n[e1[0]],end_n[e1[1]]]
                    e2_wt = D[end_n[e2[0]],end_n[e2[1]]]
                    e3_wt = D[end_n[e3[0]],end_n[e3[1]]]
                    
                    delta = e1_wt+e2_wt+e3_wt-removed_weights

                    if ((delta+S_EPS<best_delta) and
                       (not C or ( cum_d[e1[0]]+cum_d[e1[1]]-C_EPS<C and
                                   cum_d[e2[0]]+cum_d[e2[1]]-C_EPS<C and
                                   cum_d[e3[0]]+cum_d[e3[1]]-C_EPS<C ))
                       and
                       (not L or ( cum_l[e1[0]]+cum_l[e1[1]]+e1_wt-S_EPS<L and
                                   cum_l[e2[0]]+cum_l[e2[1]]+e2_wt-S_EPS<L and
                                   cum_l[e3[0]]+cum_l[e3[1]]+e3_wt-S_EPS<L ))):
                        best_move = ((i,j,k),(e1,e2,e3),delta)
                        best_delta = delta  
                        
                        
                        if strategy==LSOPT.FIRST_ACCEPT:
                            accept_move = True
                            break # move loop
                if accept_move:
                    break # k loop
            if accept_move:
                break # j loop
        if accept_move:
            break # i loop
                
    if best_move:
        # unpack the move
        (best_ijk,best_edges,best_delta) = best_move
    
        routes = [route1_data,route2_data,route3_data]
        ret = []
        for edge in best_edges:
            
            # Unfortunately, the splicing is a bit complex, but basic idea 
            #  is to use the move and its edges in best_edges to get the
            #  route to splice from and if the splice should be inverted.
            #
            # edge[i]//2 gets the route index of the move endpoint
            # edge[i]%2 tells if it is the head segment or tail segment
            #
            #    ___a c___    r1
            #   /         \
            # n0           n1
            #   \___   ___/
            #       b d       r2
            #                 
            
            r1_idx = edge[0]//2
            r2_idx = edge[1]//2
            r1_i = best_ijk[r1_idx]
            r2_i = best_ijk[r1_idx]
            a = routes[r1_idx].route[r1_i]
            b = routes[r1_idx].route[r1_i+1]
            c = routes[r2_idx].route[r2_i]
            d = routes[r2_idx].route[r2_i+1]
            
            # combine the tails of route1 and route2
            if edge[0]%2 and edge[1]%2:
                new_route = (
                        routes[r1_idx].route[None:r1_i :-1]+
                        routes[r2_idx].route[r2_i+1:None:1])
                new_cost = ( D[b, d]+
                        routes[r1_idx].rwd_l[r1_i+1]+
                        routes[r2_idx].rwd_l[r2_i+1])
                new_demand = (
                        routes[r1_idx].rwd_d[r1_i+1]+
                        routes[r2_idx].rwd_d[r2_i+1])
            # combine the head of route1 with the tail of route2
            elif not edge[0]%2 and edge[1]%2:
                new_route = (
                        routes[r1_idx].route[None:r1_i+1:1]+
                        routes[r2_idx].route[r2_i+1:None:1])
                new_cost = (D[a, d],
                        routes[r1_idx].fwd_l[r1_i]+
                        routes[r2_idx].rwd_l[r2_i+1])
                new_demand = (
                        routes[r1_idx].fwd_d[r1_i]+
                        routes[r2_idx].rwd_d[r2_i+1])
            # combine the heads of route1 and route2
            elif not edge[0]%2 and not edge[1]%2:
                new_route = (
                        routes[r1_idx].route[None:r1_i+1:1]+
                        routes[r2_idx].route[r2_i:None:-1])
                new_cost = (D[a, c],
                        routes[r1_idx].fwd_l[r1_i]+
                        routes[r2_idx].fwd_l[r2_i])
                new_demand = (
                        routes[r1_idx].fwd_d[r1_i]+
                        routes[r2_idx].fwd_d[r2_i])
            else:
                assert False, "there is no move to combine the tail of"+\
                              "route1 with the (reversed) head of route2"
                            
            ret.append(RouteData(new_route, new_cost, new_demand))
        ret.append(best_delta)
        return tuple(ret)
     
    return None, None, None, None


@routeordersensitive
def do_1point_move(route1_data, route2_data, D, d=None,
                      C=None, L=None, # constraints
                      strategy=LSOPT.FIRST_ACCEPT,
                      best_delta = None):
    """ Move one point from route1 to route2. Tries all possible combinations
    of moving a node from route 1 to any valid position in route 2. Sometimes
    called "relocate" (e.g., Bräysy & M. Gendreau 2005, Savelsbergh 1992), but
    here we use the name "one point move" (Groër et al 2010) to differentiate
    from the intra-route version.
    
    If an improving move was found and made, operation returns new routes in
    same format as route_data inputs, but if C and d are not given the 3.
    field of the tuple is None. If there is no improving move, returns None.

    Groër, C., Golden, B. and Wasil, E., 2010. A library of local search
     heuristics for the vehicle routing problem. Mathematical Programming
     Computation, 2(2), pp.79-101.
    Bräysy, O. & Gendreau, M. 2005. Vehicle Routing Problem, Part I: Route 
     Construction and Local Search Algorithms Transportation Science 39(1),
     pp. 104–118
    Savelsbergh, M. W. P. 1992. The vehicle routing problem with time windows:
     Minimizing route duration. J.Comput. 4 146-154
    """

    if route1_data==route2_data:
        return None,None,None
    
    # unpack route, current cost, and current demand
    route1, r1_l, r1_d, _ = route1_data
    route2, r2_l, r2_d, _ = route2_data
    
    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    for i in range(1,len(route1)-1):
        remove_after = route1[i-1]
        to_move = route1[i]
        remove_before = route1[i+1]
        
        remove_delta = D[remove_after,remove_before]\
                      -D[remove_after,to_move]\
                      -D[to_move,remove_before]
        
        # capacity constraint feasibility check 
        if C and r2_d+d[to_move]-C_EPS>C:
            continue

        for j in range(1,len(route2)):
            insert_after = route2[j-1]
            insert_before = route2[j]
            
            
            insert_delta = D[insert_after,to_move]+D[to_move,insert_before]\
                          -D[insert_after,insert_before]
            
            # route cost constraint feasibility check 
            if L and r2_l+insert_delta-S_EPS>L:
                continue
            
            delta = remove_delta+insert_delta
            if delta+S_EPS<best_delta:                    
                best_delta = delta
                best_move = (i, j, remove_delta, insert_delta)
                if strategy==LSOPT.FIRST_ACCEPT:
                    accept_move=True
                    # break j-loop
                    break
        # break i-loop
        if accept_move:
            break
        
    if best_move:
        # unpack best move
        i, j, remove_delta, insert_delta = best_move
        to_move = route1[i] 
        return (RouteData(route1[:i]+route1[i+1:], r1_l+remove_delta,
                 None if not C else r1_d-d[to_move]),
                RouteData(route2[:j]+[to_move]+route2[j:], r2_l+insert_delta,
                 None if not C else r2_d+d[to_move]),
                remove_delta+insert_delta)
                
    return None,None,None


def do_2point_move(route1_data, route2_data, D, d=None,
                   C=None,L=None, # constraints
                   strategy=LSOPT.FIRST_ACCEPT,
                   best_delta = None):                 
    """ Swap one point from route1 with one point on route2 if it improves the
    solution. This operation is sometimes referred to as "exchange" (e.g. in 
    Bräysy & M. Gendreau 2005, Savelsbergh 1992), but we use the name
    "two point move" (Groër et al 2010) to differentiate it from the inter-
    route one.
    
    If an improving move was found and made, operation returns new routes in
     same format as route_data inputs, but if C and d are not given the 3.
     field of the tuple is None. If there is no improving move, returns None.

    Groër, C., Golden, B. and Wasil, E., 2010. A library of local search
     heuristics for the vehicle routing problem. Mathematical Programming
     Computation, 2(2), pp.79-101.
    Bräysy, O. & Gendreau, M. 2005. Vehicle Routing Problem, Part I: Route 
     Construction and Local Search Algorithms Transportation Science 39(1),
     pp. 104–118
    Savelsbergh, M. W. P. 1992. The vehicle routing problem with time windows:
     Minimizing route duration. J.Comput. 4 146-154
    """
    
    # unpack route, current cost, and current demand
    route1, r1_l, r1_d, _ = route1_data
    route2, r2_l, r2_d, _ = route2_data
    
    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    for i in range(1,len(route1)-1):
        to_swap1 = route1[i]

        swap1_after = route1[i-1]
        swap1_before = route1[i+1]     


        for j in range(1,len(route2)-1):
            to_swap2 = route2[j]
            
            # capacity constraint feasibility check 
            if C and (r1_d-d[to_swap1]+d[to_swap2]-C_EPS>C or
                      r2_d-d[to_swap2]+d[to_swap1]-C_EPS>C):
                continue
            
            swap2_after = route2[j-1]
            swap2_before = route2[j+1]    
            
            route1_delta = -D[swap1_after,to_swap1]\
                           -D[to_swap1,swap1_before]\
                           +D[swap1_after,to_swap2]\
                           +D[to_swap2,swap1_before]
            route2_delta = -D[swap2_after,to_swap2]\
                           -D[to_swap2,swap2_before]\
                           +D[swap2_after,to_swap1]\
                           +D[to_swap1,swap2_before]
                           
            # route cost constraint feasibility check 
            if L and (r1_l+route1_delta>L or r2_l+route2_delta>L):
                continue
                        
            delta = route1_delta+route2_delta
            if delta+S_EPS<best_delta:                    
                best_delta = delta
                best_move = (i, j, route1_delta, route2_delta)
                if strategy==LSOPT.FIRST_ACCEPT:
                    accept_move=True
                    break # j loop
        if accept_move:
            break # i loop
        
    if best_move:
        # unpack best move
        i, j, route1_delta, route2_delta = best_move
        to_swap1 = route1[i] 
        to_swap2 = route2[j] 
        return (RouteData(route1[:i]+[to_swap2]+route1[i+1:],
                 r1_l+route1_delta,
                 None if not C else r1_d-d[to_swap1]+d[to_swap2]),
                RouteData(route2[:j]+[to_swap1]+route2[j+1:],
                 r2_l+route2_delta,
                 None if not C else r2_d-d[to_swap2]+d[to_swap1]),
                          
                route1_delta+route2_delta)
                
    return None,None,None

@routeordersensitive
def do_insert_move(unrouted, recieving_route_data, D,d=None,C=None,L=None,
                   strategy=LSOPT.FIRST_ACCEPT,
                   best_delta = None):
    """ Try to insert (unrouted) node(s) on the existing recieving route. The
    operation succeeds only if all customers can be inserted on the recieving
    route. The difference to the one point move is that there is no delta for
    removing the unrouted node from its RouteData / list.
    
    Returns a 3-tuple: An empty RouteData, A copy of the updated RouteData and
     the route cost (usually route length) delta; or (None,None,None) if the
     insertion fails.
    
    * unrouted
       Can be a single customer (int), a list of those, or a route (RouteData).
       If RouteData is given its .demand field must be set. 
    * recieving_route_data
       is the route (RouteData) that the customers to insert are inserted to.
       Its .demand field must be set. 
    * D,d,C,L
       the problem definition and constraints
    * strategy
       is the applied loca search strategy to use. Please note that if there is
       more than one customer to insert, the maximum route cost constraint L
       is set, and strategy is set to first accept, the best insertion position
       for the customers are still searched for to leave the most amount of
       room for subsequent insertions.
    """

    if type(unrouted) is int:
        unrouted_demand = d[unrouted] if C else 0
        unrouted = [unrouted]
    if isinstance(unrouted, RouteData):
        unrouted_demand = unrouted.demand if C else 0
        unrouted = unrouted.route[1:-1]
    elif C:
        unrouted_demand = sum(d[n] for n in unrouted)
    
    # make an ansatz for fitting the inserted customers into
    ansatz_route = list(recieving_route_data.route)
    ansatz_d = recieving_route_data.demand
    ansatz_l = recieving_route_data.cost
    total_delta = 0
    
    # is not possible because of constraint C
    if C and ansatz_d+unrouted_demand-C_EPS>C:
        return None, None, None
        
    for ni, node in enumerate(unrouted):
        best_insert_pos = None
        if strategy==LSOPT.BEST_ACCEPT or L:
            # need to find a place where it can be inserted (if at all)
            best_insert_delta = None
            for i in range(1, len(ansatz_route)):
                insert_after = ansatz_route[i-1]
                insert_before = ansatz_route[i]
                
                insert_delta = +D[insert_after,node]\
                               +D[node,insert_before]\
                               -D[insert_after,insert_before]
                        
                if not L or ansatz_l+insert_delta-S_EPS<=L:
                    if (best_insert_delta is None) or\
                       (insert_delta<best_insert_delta):
                        best_insert_pos = i
                        best_insert_delta = insert_delta 
                    # to avoid complexity, only allow early termination with 
                    #  FIRST_ACCEPT when customers are inserted one
                    #  by one
                    if strategy==LSOPT.FIRST_ACCEPT and ni==len(unrouted)-1:
                        break # i loop
        else:
            # we know it fits, just append it 
            best_insert_pos = -1
            insert_after = ansatz_route[-2] # the last non-depot node
            best_insert_delta = +D[insert_after,node]+D[node,0]-D[insert_after,0] 
        
        # no L valid routing
        if best_insert_pos is None:
            return None, None, None
        else:
            if C: ansatz_d+=d[node]
            ansatz_l+=best_insert_delta
            # we found a valid insertion location accept it
            ansatz_route = ansatz_route[:best_insert_pos]\
                           +[node]\
                           +ansatz_route[best_insert_pos:]
            total_delta+=best_insert_delta
            
            # must be better than the preset level (if applicable)
            if (best_delta is not None) and (total_delta+S_EPS>best_delta):
                return None, None, None
            
    return RouteData(), RouteData(ansatz_route, ansatz_l, ansatz_d), total_delta

@routeordersensitive
def do_redistribute_move(redisributed_route_data,
                         receiving_route_or_routes_data,
                         D,d=None, C=None, L=None,
                         strategy=LSOPT.FIRST_ACCEPT,
                         best_delta = None,
                         recombination_level=0):
    """
    Try to insert the nodes of the first route on the other routes if possible.
    Note that second argument receiving_route_data can also be a list of routes.
    
    The insertion order matters, so it is possible to try all permutations
    of nodes to be inserted and routes to insert to (and return the best).
    However, this will be very computationally intesive as the number of
    combinations explode. Therefore, there are different options available:
        
    * recombination_level=0, insert the nodes in the order they are on the
       route that will be destroyed and consider the recieving routes in the
       order of the "routes" parameter.
    * recombination_level=1, the nodes are inserted in all possible 
       order. The recieving routes are checked in the order of the "routes"
       parameter.
    * recombination_level=2, insert the nodes in the order they are on the
       route that will be destroyed. However, every possible order of the 
       recieving routes is checked.
    * recombination_level=3, all possible orderings of the nodes to be
       inserted AND the recieving routes is checked.       
       
    There could be yet another level (try_all_insertions_level=4), where for
    each insertion of each insertion node ordering, but it has not been
    implemented.
    """
 
    # allow only single redistributing to a single other route e.g. when using
    #  (local_search API)
    if isinstance(receiving_route_or_routes_data, RouteData):
        receiving_route_or_routes_data = [receiving_route_or_routes_data]
    
    # before iterating the permuations, check that it is at all possible to 
    #  fit all from the route to be destroyed to the other routes.
    if C:
        total_d_slack = sum(C-rrd.demand for rrd in receiving_route_or_routes_data)
        if redisributed_route_data.demand-C_EPS > total_d_slack:
            return [None]*(len(receiving_route_or_routes_data)+2)
        
    if not best_delta:
        # allows move that makes the overall solution cost worse
        best_delta = float("inf")
    redistribution_succesfull = False    
    best_insertion_routes = None
 
    # unpack route, current cost, and current demand
    route1, r1l, r1d, _ = redisributed_route_data    
    if recombination_level==0:
        insertion_node_orderings = [(route1[1:-1], receiving_route_or_routes_data)]
    elif recombination_level==1:
        insertion_node_orderings = product(
            permutations(route1[1:-1]), 
            [receiving_route_or_routes_data])
    elif recombination_level==2:
        insertion_node_orderings = product(
            [route1[1:-1]], 
            permutations(receiving_route_or_routes_data))
    elif recombination_level==3:
        insertion_node_orderings = product(
            list(permutations(route1[1:-1])), 
            list(permutations(receiving_route_or_routes_data)))
    
    #print(list(po for po in insertion_node_orderings))
    for route_customers, candidate_routes in insertion_node_orderings:
        total_delta = 0
        all_succesfully_inserted = True        
        ansatz_route_data = list(candidate_routes)
        
        for to_insert in route_customers:            
            node_succesfully_inserted = False
            for rd2_ar_idx, rd2 in enumerate(ansatz_route_data):
                _, new_rd2, delta = do_insert_move(to_insert, rd2, D,d,C,L,
                                                   strategy=strategy)
                if delta is None:
                    continue
                
                #print("REMOVEME:", "succesfully inserted %d to create %s"%(to_insert, new_rd2.route))
                ansatz_route_data[rd2_ar_idx] = new_rd2
                total_delta+=delta
                if total_delta+S_EPS<best_delta:
                    node_succesfully_inserted = True
                    break
                else:
                    continue
                
            if not node_succesfully_inserted:
                all_succesfully_inserted = False
                break
            
        if all_succesfully_inserted:
            #print("REMOVEME:", "all inserted to create %s (delta %.2f)"%([rd.route for rd in ansatz_route_data], total_delta))   
            if total_delta<best_delta:
                best_delta = total_delta
                best_insertion_routes = ansatz_route_data
                redistribution_succesfull = True
                    
    if redistribution_succesfull:
        # return an empty route and n routes the customers were distributed to
        return tuple( [RouteData()]+best_insertion_routes+[best_delta] )
    else:
        return [None]*(len(receiving_route_or_routes_data)+2)
        
@routeordersensitive
def do_chain_move(route1_data, route2_data, route3_data, D, d=None,
                      C=None, L=None, # constraints
                      strategy=LSOPT.FIRST_ACCEPT,
                      best_delta = None):  
    """ This is the "pair" operation described in Wren and Holliday (1972). 
    It involves moving a repacing a node on route 2 with a node on route 1.
    The replaced node is then inserted on route 3 (if able). Route 1!=2!=3 
    
    This is a very expensive operation, corresponding to 5-opt with chain
    length of 1, use with care.
    
    Wren, A. and Holliday, A., 1972. Computer scheduling of vehicles from one
    or more depots to a number of delivery points. Journal of the Operational
    Research Society, 23(3), pp.333-344.
    """

    if not best_delta:
        best_delta = 0
    best_move = None
    accept_move = False
    
    route1, r1_l, r1_d, _ = route1_data    
    route2, r2_l, r2_d, _ = route2_data
    route3, r3_l, r3_d, _ = route3_data
    
    for i in range(1,len(route1)-1):
        
        remove_after = route1[i-1]
        to_move = route1[i]
        remove_before = route1[i+1]
        
        remove_delta =  +D[remove_after,remove_before]\
                        -D[remove_after,to_move]\
                        -D[to_move,remove_before]
        
        for j in range(1,len(route2)-1):
            
            replace_after = route2[j-1]
            to_replace = route2[j]
            replace_before = route2[j+1]
            
            # can do capacity constraint feasibility check here
            if C and (r2_d-d[to_replace]+d[to_move]>C or
                      r3_d+d[to_replace]>C): 
                continue
                            
            replace_delta = +D[replace_after,to_move]\
                            +D[to_move,replace_before]\
                            -D[replace_after,to_replace]\
                            -D[to_replace,replace_before]
                            
            if L and r2_l+replace_delta>L:
                continue
        
            for k in range(1,len(route3)):
                
                insert_after = route3[k-1]
                insert_before = route3[k]
                
                insert_delta =  +D[insert_after, to_replace]\
                                +D[to_replace, insert_before]\
                                -D[insert_after,insert_before]
                                
                if L and (r3_l+insert_delta>L):
                    continue

                # check if this is best so far            
                delta = remove_delta+replace_delta+insert_delta
                if delta+S_EPS<best_delta:                    
                    best_delta = delta
                    best_move = (i, j, k, remove_delta, replace_delta, insert_delta)
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move=True
                        break # k loop
            if accept_move:
                break # j loop
        if accept_move:
            break # i loop
            
    if best_move:
        # unpack best move
        i, j, k, remove_delta, replace_delta, insert_delta = best_move
        to_move = route1[i] 
        to_replace = route2[j] 
        return (RouteData(route1[:i]+route1[i+1:],
                    r1_l+remove_delta,
                    None if not C else r1_d-d[to_move]),
                RouteData(route2[:j]+[to_move]+route2[j+1:],
                    r2_l+replace_delta,
                    None if not C else r2_d-d[to_replace]+d[to_move]),
                RouteData(route3[:k]+[to_replace]+route3[k:],
                    r3_l+insert_delta,
                    None if not C else r3_d+d[to_replace]),
                best_delta)
    return None,None,None,None             
