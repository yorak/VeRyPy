# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and contains an INCOMPLETE IMPLEMENTATION of the 3-opt inter
route (aka. 3-opt*) local search improvement heuristic. Thus, for a complete 
working and tested 3-opt* implementation please see _solution_operators.py_ 

This implementation is left here for some brave soul to later complete it.
The reason for a separate non-full-solution 3-opt* would be to make the 
3-opt* more performant because the more granular control of the 
`do_local_search` function in _local_search_ module allows avoiding already
tested moves. Also, it would be convenient to have a more granular control
over 3-opt*, for example for a more geographically targeted improvement
attempts.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from itertools import product, permutations

from routedata import RouteData
from local_search import LSOPT, ROUTE_ORDER_SENSITIVE_OPERATORS
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

# Create a decorator to mark routes sensitive to the route order
def routeordersensitive(f):
    ROUTE_ORDER_SENSITIVE_OPERATORS.add(f)
    return f
    

# This defines all possible 3opt moves involving two routes. The data describes
#  which segment ends to connect when attempting the move. This makes implementing
#  the do_3optstar_[2|3]route_move (below) operator simpler.
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
     
