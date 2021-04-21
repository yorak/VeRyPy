# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and contains a set of naive implementations of the popular 
local search operations.
 
Note: these are for educational use only as the objective function and feasibi-
      lity (if needed) are computed for each possible modification of the
      solution. Thus, it is quite slow. However, the implementations are more
      readable than the more optimized ones found in inter_route_operations.py
      and intra_route_operations.py. The implementations in this file are used
      to verify the correct operation of the optimized code (see.
      tests/test_local_search_parallel.py
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from itertools import groupby
from time import time

from util import objf
from local_search import LSOPT
from cvrp_ops import fast_constraint_check

from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

def _apply_as_intra_route(operation,
                     solution,D,d,C,L,
                     strategy=LSOPT.BEST_ACCEPT,
                     best_delta=None):
    """ This allows using an educational version of the move operations to 
    implement the intra route (per route) variants.
    
    Note that the operation is accepted on all routes if improving moves 
    are found. Even if the strategy is LSOPT.BEST_ACCEPT (this only selects 
    the best move *within* the route).
    """
    
    total_delta = 0
    if not best_delta:
        best_delta = 0
    route_start = 0
    route_end = 1 

    ansatz_sol = list(solution)
    while route_end<len(ansatz_sol):
        
        # find where the route ends
        while ansatz_sol[route_end]!=0:
            route_end+=1
            
        route = ansatz_sol[route_start:route_end+1]
        new_route, improvement = operation(route,D,d,C,L,
                               strategy=strategy,
                               best_delta=None)
        if improvement:
            total_delta+=improvement
            ansatz_sol[route_start:route_end+1] = new_route
            
        route_start = route_end
        route_end = route_start+1
    
    if total_delta+S_EPS<best_delta:
        return ansatz_sol, total_delta
    else:
        return None, None   

def do_naive_2opt_move(solution,D,d,C,L,
                       strategy=LSOPT.BEST_ACCEPT,
                       best_delta=None):
    """ This is an educational version of the 2-opt intra route move. The move
    tries to find a best (default) or first improving reversal of a sequence
    of customers.
    """
    
    sol_f = objf(solution, D)
    best_sol = None
    best_found = False
    if not best_delta:
        best_delta = 0
    
    route_start = 0
    route_end = 1 

    while route_end<len(solution):
        
        # find where the route ends
        while solution[route_end]!=0:
            route_end+=1
        
        # all possible ways of reversing the intra route sequence of length 2+
        for a in range(route_start, route_end):
            b = a+1
            for c in range(a+2, route_end):
                d = c+1
                
                # a->c b->d
                #       ______
                #      /      \
                # 0->-a   b-<-c   d->-0
                #         \______/  
                #
                ansatz_sol = solution[:b]+\
                             solution[c:a:-1]+\
                             solution[d:]
                             
                ansatz_sol_f = objf(ansatz_sol, D)      
                delta = ansatz_sol_f-sol_f
                if delta+S_EPS<best_delta:            
                    best_delta = delta
                    best_sol = ansatz_sol
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        best_found = (a,c)
                        break # reverse_end loop
            if best_found:
                break # reverse_start loop
        if best_found:
            break # route loop
            
        route_start = route_end
        route_end = route_start+1
        
    if best_sol:
        print("best_found", best_found)
        return best_sol, best_delta
    else:
        return None, None          
        
 
def do_naive_2optstar_move(solution,D,d,C,L,
                       strategy=LSOPT.BEST_ACCEPT,
                       best_delta=None):
    """ A naive and inefficient inter-route version of the do_naive_2opt_move.
    Used to verify the operation of inter_route_operations.do_2optstar_move.
    """
    # free d for other uses
    demands = d
    
    sol_f = objf(solution, D)
    best_sol = None
    best_found = False
    if not best_delta:
        best_delta = 0
        
    # all possible ways of reconnecting the sequences
    for a in range(0, len(solution)-3):
        b = a+1
        left_mid_depot = b if solution[b]==0 else None
        right_mid_depot = left_mid_depot
        for c in range(b, len(solution)-1):
            d = c+1
            if solution[c]==0:
                if left_mid_depot is None:
                    left_mid_depot = c
                right_mid_depot = c
            
            for alternative in [1,2]:
                if alternative == 1:
                    #print("\nalt1",a,b,c,d, "depot", mid_depot)
                    # a->c b->d
                    #       ______________
                    #      /              \
                    # 0->-a   b-<-0->-0-<-c   d->-0
                    #         \_______________/  
                    #
                    if left_mid_depot is None:
                        ansatz_sol = solution[:b]+\
                                     solution[c:a:-1]+\
                                     solution[d:]
                    else:
                         ansatz_sol = solution[:b]+\
                                 solution[c:right_mid_depot:-1]+\
                                 solution[left_mid_depot:right_mid_depot+1]+\
                                 solution[left_mid_depot-1:a:-1]+\
                                 solution[d:]
                if alternative == 2:
                    # This is not possible if there is no depot on the string
                    if left_mid_depot is None:
                        continue
                    #print("\nalt2",a,b,c,d, "depot", mid_depot)
                    # a->d c->b
                    #       __________________
                    #      /                  \
                    # 0->-a   b-<-0->-0-<-c   d->-0
                    #         \___________/  
                    #
                    
                    next_depot_from_d = d+solution[d:].index(0)
                    ansatz_sol = solution[:b]+\
                                 solution[d:next_depot_from_d]+\
                                 solution[left_mid_depot:right_mid_depot+1]+\
                                 solution[left_mid_depot-1:a:-1]+\
                                 solution[c:right_mid_depot:-1]+\
                                 solution[next_depot_from_d:]
                    
                ansatz_sol_f = objf(ansatz_sol, D)      
                delta = ansatz_sol_f-sol_f

         
                if delta+S_EPS<best_delta:
                    if fast_constraint_check(ansatz_sol,D,demands,C,L):
                        best_delta = delta
                        best_sol = ansatz_sol
                        
                        if strategy==LSOPT.FIRST_ACCEPT:
                            best_found = True
                            break # alternative loop
            if best_found:
                break # second edge loop          
        if best_found:
            break # first edge loop
        
        
    if best_sol:
        # remove 0,0´s
        best_sol = [x[0] for x in groupby(best_sol)]
        
        return best_sol, best_delta
    else:
        return None, None          
       
        
def do_naive_3optstar_move(solution,D,demands,C,L,
                       strategy=LSOPT.BEST_ACCEPT,
                       best_delta=None):
    """ A naive and inefficient inter-route version of the move.
    Used to verify the operation of inter_route_operations.do_naive_3optstar_move.
    """
    
    # The implementation is incomplete, sorry!
    raise NotImplementedError
    
    sol_f = objf(solution, D)
    best_sol = None
    best_found = False
    if not best_delta:
        best_delta = 0
        
    # all possible ways of reconnecting the sequences
    for a in range(0, len(solution)-3):
        b = a+1
        ab_left_mid_depot = b if solution[b]==0 else None
        ab_right_mid_depot = ab_left_mid_depot
        
        for c in range(b, len(solution)-1):
            d = c+1
            for e in range(d, len(solution)):
                f = e+1
                
                #       ______________
                #      /              \
                # 0->-a   b-<-0->-0-<-c   d-0->->0->-e->-f->-0
                #         \_______________/  
                #
        
                ansatz_sol = solution[:b]+\
                             solution[c:a:-1]+\
                             solution[d:]
            
                ansatz_sol = solution[:d]+\
                             solution[e:c:-1]+\
                             solution[f:]
                
                #       _________________________________
                #      /                                 \   
                # 0->-a   b-<-0->-0-<-c---d-0->->0->-e   f->-0
                #         \__________________________/  
                #
                ansatz_sol = solution[:b]+\
                             solution[f:]+\
                             solution[b:f]
                             
                ansatz_sol = solution[:d]+\
                             solution[e:c:-1]+\
                             solution[f:]
                
                #TODO: Write naive 3-opt* implementation some day
                raise NotImplementedError("The other move operators have not been implemented")
                             
            
            if best_found:
                break # second edge loop          
        if best_found:
            break # first edge loop
        
        
    if best_sol:
        # remove 0,0´s
        best_sol = [x[0] for x in groupby(best_sol)]
        
        return best_sol, best_delta
    else:
        return None, None   
           
    
    # these take one route as list of nodes, e.g. [0,1,2,3,0]
    #do_3opt_move:1,
    #do_relocate_move:1,    
    #do_exchange_move:1,  
    # these take multiple routes as RouteData (see routedata.py) objects
    #do_2optstar_move:2,
    #do_1point_move:2,
    #do_2point_move:2,
    #do_insert_move:2, 
    #do_redistribute_move:2, # could be also 3,4 etc.
    #do_chain_move:3

def do_naive_1point_move(solution,D,d,C,L,
                         strategy=LSOPT.BEST_ACCEPT,
                         best_delta=None):
    """ This is an educational version of the one point/relocate move.
    A customer is moved to the different position of the same or different
    route if it improves the solution.
    """
    
    sol_f = objf(solution, D)
    best_sol = None
    best_found = False
    if not best_delta:
        best_delta = 0
    
    for n in solution:
        if n==0:
            continue
            
        for j in range(1,len(solution)-1):
            ansatz_sol = list(solution)
            ansatz_sol.remove(n)
            ansatz_sol.insert(j,n)
       
            if fast_constraint_check(ansatz_sol,D,d,C,L):
                ansatz_sol_f = objf(ansatz_sol, D)
                delta = ansatz_sol_f-sol_f

                if delta+S_EPS<best_delta:
                    best_delta = delta
                    best_sol = ansatz_sol
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        best_found = True
                        break
        if best_found:
            break
        
    if best_sol:
        return best_sol, best_delta
    else:
        return None, None
    
    
def do_naive_relocate_move(solution,D,d,C,L,
                         strategy=LSOPT.BEST_ACCEPT,
                         best_delta=None):
    """ This is an educational version of the intra-route relocate move. A
    customer is moved to a different position on the same route if it improves
    the solution. It uses do_naive_1point_move to implement the move.
    
    As do_naive_1point_move also considers intra-route moves, this is a more
    restricted operation.
    """
    return _apply_as_intra_route(do_naive_1point_move,
                                 solution,D,d,C,L,
                                 strategy=strategy,
                                 best_delta=None)

def do_naive_2point_move(solution,D,d,C,L,
                           strategy=LSOPT.BEST_ACCEPT,
                            best_delta=None):
    """ This is an educational version of the inter route node exchange move
    (two point move). A customer is swapped with another on the same or
    different route if it improves the solution.
    """
    
    sol_f = objf(solution, D)
    best_sol = None
    best_found = False
    if not best_delta:
        best_delta = 0
    
    for i in range(len(solution)):
        n1 = solution[i]
        if n1==0:
            continue
            
        for j in range(1,len(solution)-1):
            n2 = solution[j]
            if n2==0:
                continue
        
            ansatz_sol = list(solution)
            ansatz_sol[i], ansatz_sol[j] = ansatz_sol[j], ansatz_sol[i]
            
#            print("\nmove", i, j)
#            print("orig",solution )            
#            print("check",ansatz_sol )
#            print("quality_delta", objf(ansatz_sol, D)-objf(solution, D))
#            print("feasibility",fast_constraint_check(ansatz_sol,D,d,C,L) )
            
            if fast_constraint_check(ansatz_sol,D,d,C,L):
                ansatz_sol_f = objf(ansatz_sol, D)
                delta = ansatz_sol_f-sol_f
                
                if delta+S_EPS<best_delta:
                    #print("set as the best")
                
                    best_delta = delta
                    best_sol = ansatz_sol
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        best_found = True
                        break
        if best_found:
            break
        
    if best_sol:
        return best_sol, best_delta
    else:
        return None, None

def do_naive_exchange_move(solution,D,d,C,L,
                         strategy=LSOPT.BEST_ACCEPT,
                         best_delta=None):
    """ This is an educational version of the intra-route node exchange move 
    (aka. swap). A customer is swapped with another customer on the same route 
    if it improves the solution. It uses the inter-route do_naive_2point_move
    to implement the move.
    """
    return _apply_as_intra_route(do_naive_2point_move,
                                 solution,D,d,C,L,
                                 strategy=strategy,
                                 best_delta=None)
    
def do_naive_local_search(ls_ops, sol, D, d, C, L = None,
                          operator_strategy=LSOPT.FIRST_ACCEPT,
                          max_iterations=None):
    """ Repeatedly apply naive ls_ops until no more improvements can be made.
    Optionally maximum number of iterations (of applying all ls_ops operators)
    can be given.
    """

    iteration = 0
    improved = True
    while improved:
        improved = False    
        for ls_op in ls_ops:
            start_t = time()
            new_sol, delta = ls_op(sol,D,d,C,L,operator_strategy)
            if delta is None:
                continue
            else:
                improved = True
            elapsed_t = time()-start_t 
            if __debug__:
                print("%s improved from %s (%.2f) to %s (%.2f) in %.2f s"%
                  (ls_op.__name__, sol, objf(sol, D),
                   new_sol, objf(new_sol, D), elapsed_t))
            
            sol = new_sol            
            iteration += 1
            if max_iterations and iteration==max_iterations:
                break
    return sol

        
    
