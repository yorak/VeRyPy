# -*- coding: utf-8 -*-
###############################################################################
""" This file implements a basic hill climbing (or valley seeking, as we are
minimizing) local search procedure, which applies the given intra- or inter-
route operations repeadedly until a local optima is reached.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from collections import defaultdict, Counter
from logging import log, DEBUG
from itertools import permutations
from inspect import getargspec

from routedata import RouteData
from util import objf, is_sorted
from config import COST_EPSILON as S_EPS

# enum like
class LSOPT:
    FIRST_ACCEPT = 1 # Accept the first improving move
    BEST_ACCEPT = 2 # Accept the best improving move
class ITEROPT(LSOPT):  
    ALL_ACCEPT = 3 # Improve all improving moves as they are discovered
    REPEATED_ACCEPT = 4 # Check if the same move and be made again to improve even more

# All operators that are sensitive to route order must be added here.
#  See inter_route_operations.py decorator for details.
ROUTE_ORDER_SENSITIVE_OPERATORS = set()

def do_local_search(ls_ops, sol, D, d, C, L=None,
                    operator_strategy=LSOPT.FIRST_ACCEPT,
                    iteration_strategy=ITEROPT.ALL_ACCEPT,
                    max_iterations=None):
    """ Repeatedly apply ls_ops until no more improvements can be made. The
    procedure keeps track of the changed routes and searches only combinations
    that have been changed.
    
    Optionally the operator_strategy FIRST_ACCEPT (default)/BEST_ACCEPT can be
    given as well as the maximum number of iterations (that is, how many times
    all given operations are applied until giving up on reaching local optima).
    
    The iteration_strategy has an effect on which order the operations
    are applied. The default is ALL_ACCEPT (default), but the options are:
        
     * ITEROPT.FIRST_ACCEPT accepts every improving move found by the active 
        operator, and after an improving move is found and made, it starts
        again from the first operator of the list. This is repeated until no
        operator finds any improvements.
     * ITEROPT.BEST_ACCEPT accepts the very best (single) move over all
        operators. Stops when no operator finds any improvements. This, when
        used together with operator_strategy=LSOPT.BEST_ACCEPT, will be the
        very steepest descent. However, it is also the most computationally
        heavy alternative as it has to evaluate every operator for each move.
     * ITEROPT.ALL_ACCEPT accepts every improving move of each operator it
        finds. It continues to the next operator after making that move or
        noting that there is no improvement to be made. The strategy stops when
        no operator finds any improvements.
     * ITEROPT.REPEATED_ACCEPT runs a single operator until no improving moves 
        are found for it. Then, the strategy moves on to the next operator. The
        strategy starts again from the first operator if any of the operators
        found an improvement.
        
    Note that these may freely be combined with the operator_strategy.    
    """
    
    current_sol = sol
    route_datas = RouteData.from_solution(sol, D, d)
    route_data_idxs = list(range(len(route_datas)))

    # We keep track of the operations to avoid search when it has already been
    #  unsuccesfully applied   
    at_lsop_optimal = defaultdict(set)
    customer_to_at_lsopt_optimal = defaultdict(list)
    
    iteration = 0
    improving_iteration = True
    while improving_iteration:
        improving_iteration = False
        
        best_iteration_result = None
        best_iteration_delta = None 
        best_iteration_operator = None
        
        ls_op_idx = 0
        while ls_op_idx<len(ls_ops):
            ls_op = ls_ops[ls_op_idx]
            # TODO: For full Python 3.X replace this with modern approach.
            ls_op_args = getargspec(ls_op)[0]
            route_count = ls_op_args.index('D')
            op_order_sensitive = ls_op in ROUTE_ORDER_SENSITIVE_OPERATORS
            
            op_improved = False
            
            if __debug__:
                log(DEBUG-1, "Applying %s on %s"%(ls_op.__name__, str(current_sol)))
            
            # TODO: Consider using a counter to check for this
            # check if we already reached local optima on all routes with ls_op
            #if all( (ls_op in lsop_optimal[ri]) for ri in route_data_idxs ):
            #    if __debug__:
            #        log(DEBUG-2, "All route combinations already searched for %s, skipping it."%ls_op.__name__)
            #    break
            
            best_delta = None
            best_result = None
            
            no_improving_lsop_found = set()                
            for route_indices in permutations(route_data_idxs,route_count):
                # If the order does not matter, require that the route indices
                #  are ordered from smallest to largest.
                if (not op_order_sensitive) and (not is_sorted(route_indices)):
                    continue
                
                # ls_op is already at local optima with this combination of routes
                if ls_op in at_lsop_optimal[route_indices]:
                    if __debug__:
                        log(DEBUG-2, "Route combination %s already searched for %s, skipping it."%
                            (str(route_indices), ls_op.__name__))
                    continue

                # The one route case has different call signature
                if route_count==1:
                    op_params = [route_datas[route_indices[0]].route,
                                 D, operator_strategy]
                else:
                    op_params = [route_datas[ri] for ri in route_indices]+\
                                 [D, d, C, L, operator_strategy]
                                 # Ideally, best_delta can be used as an upper
                                 # bound to avoid unnecessary result generation
                                 # and to allow early ls_op termination.
                                 # However, then we lose the ability to mark
                                 # some route combinations as ls_optimal.
                                 #+[best_delta]
                result = ls_op(*op_params)
                #print("REMOVEME:",route_datas[route_indices[0]].route, "->", result)
                
                # route was changed, record the change in route datas
                delta = result[-1]
                if delta is None:
                    no_improving_lsop_found.update((route_indices,))
                else:
                    # For route_count==1 every route contributes for the same
                    # best_delta (unless trying to find the very best *single*
                    # move!)
                    if route_count==1:
                        skip_result = (
                            (best_delta != None and delta+S_EPS>best_delta) and
                            (iteration_strategy==ITEROPT.BEST_ACCEPT) )
                        
                        if not skip_result:
                            if ((best_result is None) or
                                (iteration_strategy==ITEROPT.BEST_ACCEPT)):
                                best_result = []
                                best_delta = 0
                           
                            old_rd = route_datas[route_indices[0]]
                            new_rd = RouteData(result[0],old_rd.cost+delta,old_rd.demand)
                            best_result.append( (route_indices[0], new_rd) )
                            best_delta+=delta
                    else:
                        if (best_result is None) or (delta+S_EPS<best_delta):
                            best_result = zip(route_indices, result[:-1])
                            best_delta = delta
                    
                    # Found a first improving with this operator, move on.
                    if operator_strategy==LSOPT.FIRST_ACCEPT:
                        break # route combination loop
                
            # end route combination loop
                        
            # Mark the routes that had no potential improvements to be at
            #  local optima to avoid checking the same moves again.
            for ris in no_improving_lsop_found:
                at_lsop_optimal[ris].add(ls_op)
                for ri in ris:
                    customer_to_at_lsopt_optimal[ri].append(ris)
                
            if best_result is not None:    
                if iteration_strategy==ITEROPT.BEST_ACCEPT:
                    if (best_iteration_result is None) or \
                       (best_delta+S_EPS<best_iteration_delta):
                        best_iteration_result = best_result
                        best_iteration_delta = best_delta 
                        best_iteration_operator = ls_op.__name__
                else:
                    op_improved = True
                    improving_iteration = True
                    for ri, new_rd in best_result:
                        route_datas[ri] = new_rd
                        # The route was modified, allow other operators to 
                        #  check if it can be improved again.
                        for ris in customer_to_at_lsopt_optimal[ri]:
                            at_lsop_optimal[ris].clear()
                            
                        # Check if route is [0,0] or [0] or []
                        if len(new_rd.route)<=2:
                            # remove this route from the future search 
                            route_data_idxs.remove(ri)
                        
                    if __debug__:
                        op_improved = True
                        opt_sol = RouteData.to_solution(route_datas)
                        log(DEBUG, "Improved from %s (%.2f) to %s (%.2f) using %s"%
                                (str(current_sol),objf(current_sol,D),str(opt_sol),objf(opt_sol,D),ls_op.__name__))
                        current_sol = opt_sol
                        
                    if iteration_strategy==ITEROPT.FIRST_ACCEPT:
                        ls_op_idx = 0
                        break # the ls_op loop (start from the beginning)
                    
            if __debug__:
                 if best_result is None:
                    log(DEBUG-1, "No improving move with %s"%ls_op.__name__)
                
            if op_improved and iteration_strategy==ITEROPT.FIRST_ACCEPT:
                # after an improvement start from the first operator
                ls_op_idx = 0
            if op_improved and iteration_strategy==ITEROPT.REPEATED_ACCEPT:
                # keep repeating the operator until no improvement is found
                ls_op_idx = ls_op_idx 
            else:
                # BEST_ACCEPT and ALL_ACCEPT always move on
                ls_op_idx += 1
                
            #END OF LS_OP LOOP
        
        if (iteration_strategy==ITEROPT.BEST_ACCEPT) and\
           (best_iteration_result is not None):
            improving_iteration = True
            
            for ri, new_rd in best_iteration_result:
                route_datas[ri] = new_rd
                # The route was modified, allow other operators to 
                #  check if it can be improved again.
                for ris in customer_to_at_lsopt_optimal[ri]:
                    at_lsop_optimal[ris].clear()
                # Check if route is [0,0] or [0] or []
                if len(new_rd.route)<=2:
                    # remove this route from the future search 
                    route_data_idxs.remove(ri)

            if __debug__:
                op_improved = True
                opt_sol = RouteData.to_solution(route_datas)
                log(DEBUG, "Improved from %s (%.2f) to %s (%.2f) using %s"%
                        (str(current_sol),objf(current_sol,D),str(opt_sol),objf(opt_sol,D),best_iteration_operator))
                current_sol = opt_sol
        
        iteration+=1
        if max_iterations and iteration>=max_iterations:
            break # iteration loop 

    current_sol = RouteData.to_solution(route_datas)              
    if __debug__:
        log(DEBUG,"Repeadedly applying %s resulted in %s"%
            (",".join(ls_op.__name__ for ls_op in ls_ops),str(current_sol)))
                  
                  
    return current_sol