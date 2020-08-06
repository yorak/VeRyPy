#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Mole & Jameson (1976) 
sequential cheapest insertion heuristic.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has moderate dependencies: a TSP solver (the built in local
search solver can be used), and numpy and scipy for reading and preparing the
TSPLIB problem instance.
"""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from logging import log, DEBUG

from classic_heuristics.cheapest_insertion import parametrized_insertion_criteria,\
                               cheapest_insertion_init
from local_search import LSOPT, do_local_search
from local_search.intra_route_operators import do_2opt_move
from local_search.inter_route_operators import do_1point_move,\
                                               do_redistribute_move
from util import objf, without_empty_routes, is_better_sol
from routedata import RouteData
from config import COST_EPSILON as S_EPS


# One of the few non standard-lib additions, a doubly linked list. Used to make
#  insertions constant time operations. Install it e.g. with:
#
# $ pip install llist
#
from llist import dllist

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def _try_insert_2opt_and_update(insertion, rd, D, L, minimize_K):
    
    inserted = rd.route.insert(insertion.customer, insertion.before_node)   
    #print("inserting", insertion, "resulting to", list(rd.route))

    ## keep the route 2-optimal
    # unfortunately we cannot do local_search.py do_2opt_move on the reoute,
    #  because it is an dllist. However, as the route is kept 2-optimal, it 
    #  is probable that not all search attempts lead to updating the route,
    #  so convert and udpate only as needed.
    route_2opt_improved = list(rd.route)
    total_delta = insertion.cost_delta
    applied_2opt_moves = False
    while True:
        updated_route, delta = do_2opt_move(route_2opt_improved, D,
                                            strategy=LSOPT.BEST_ACCEPT)
        if delta is not None:
            applied_2opt_moves = True
            total_delta+=delta
            route_2opt_improved = updated_route
        else:
            break
    
    # do not accept insertions/moves that make the solution worse!
    if not minimize_K:
        # Compared to a solution where the customer is served individually
        insertion_cost = total_delta-D[0,insertion.customer]-D[insertion.customer,0]
        if insertion_cost > 0:
            if __debug__:
                log(DEBUG-2,("Rejecting insertion of n%d "%insertion.customer)+
                            ("as it is expected make solution worse."))
            rd.route.remove(inserted) # reverse the change
            return False, None
    
    if L and rd.cost+total_delta-S_EPS > L:
        # The L constraint cannot be satisfied, even with 2-opt
        rd.route.remove(inserted) # reverse the change
        return False, None
    
    rd.cost += total_delta
    rd.used_capacity += insertion.demand_delta 
        
    if applied_2opt_moves:
        # this is very ineffective, but as the dllist nodes and next/prev
        #  are read only, and the llist module does not offer a sequence 
        #  reversing functionality (!) we have no choice but scrap the entire 
        #  insertion queue and calculate the insertions all over for the 
        #  updated route. :(
        #
        #TODO: an improvement would be to use an another doubly-linked-list 
        # implementaiton that allows manipulation of the nodes / list
        del rd.potential_insertions[:]
        
        updated_dllist = dllist(route_2opt_improved)
        new_edges = []
        prev_node = updated_dllist.first
        current_node = prev_node.next
        while prev_node!=updated_dllist.last:
            new_edges.append( (prev_node, current_node) )
            prev_node = current_node
            current_node = prev_node.next

        if __debug__:
            log(DEBUG,"Applying a 2-opt moves, which transforms "+
                " %s to %s with an %.2f improvement"%
                (str(list(rd.route)),str(route_2opt_improved),
                 total_delta-insertion.cost_delta))
        
        rd.route = updated_dllist
        return True, new_edges
    
    else:
        new_edges = [ (insertion.after_node, inserted), (inserted, insertion.before_node) ]
        return True, new_edges

def _create_new_criteria_function(lm,mm):
    """ Due to how Python manages scope, the creation of closure is wrapped 
    into a function. The lambda cannot be created in the loop, because then
    all function calls will refer to the same variable!.
    see e.g. https://stackoverflow.com/questions/2295290 """
    return lambda D, i, u, j: parametrized_insertion_criteria(D, i, u, j,
                                                              lm=lm, mm=mm)
 
def _refine_solution(sol,  D, d, C, L, minimize_K):
    # refine until stuck at a local optima
    local_optima_reached = False
    while not local_optima_reached:
        sol = without_empty_routes(sol)
        if not minimize_K:
            sol.append(0) #make sure there is an empty route to move the pt to
            
        # improve with relocation and keep 2-optimal
        sol = do_local_search([do_1point_move, do_2opt_move], sol,
                              D, d, C, L, LSOPT.BEST_ACCEPT)
        
        # try to redistribute the route with smallest demand
        sol = without_empty_routes(sol)
        routes = RouteData.from_solution(sol, D, d)
        min_rd = min(routes, key=lambda rd:rd.demand)
        routes.remove(min_rd)
        
        if not minimize_K:
            routes.append(RouteData())
        
        if __debug__:
            log(DEBUG,"Applying do_redistribute_move on %s (%.2f)"%(str(sol), objf(sol,D)))
      
        redisribute_result = do_redistribute_move(min_rd, routes,
                           D,d,C,L,strategy=LSOPT.FIRST_ACCEPT,
                           #Note: Mole and Jameson do not specify exactly 
                           # how the redistribution is done (how many 
                           # different combinations are tried).
                           # Increase the recombination_level if the 
                           # for more agressive and time consuming search
                           # for redistributing the customers on other
                           # routes.
                           recombination_level=0)
        redisribute_delta = redisribute_result[-1]
        
        if (redisribute_delta is not None) and\
           (minimize_K or redisribute_delta<0.0):
               
            updated_sol = RouteData.to_solution(redisribute_result[:-1])
            if __debug__:
                log(DEBUG-1,("Improved from %s (%.2f) to %s (%.2f)"%
                    (sol, objf(sol, D), updated_sol, objf(updated_sol, D)))
                    +"using inter route heuristic do_redistribute_move\n")
            sol = updated_sol
        else:
            local_optima_reached = True
            if __debug__:
                log(DEBUG-1,"No move with do_redistribute_move\n")
    return sol
    
def mole_jameson_insertion_init(D, d, C, L=None, minimize_K=False,
                                strain_criterion='all' ):
    """ This is the implementation of Mole and Jameson (1976) cheapest
    insertion algorithm. The emerging route is first initialized according to
    which strain criterion (insertion cost calculation method) is used,
    On each step an unrouted customer for which the insertion cost is lowest
    (between any two nodes on the emerging route) is searched and the insertion
    made until no feasible insertions remain. For details see the insertion 
    implementation in cheapest_insertion.py:cheapest_insertion_init
    
    * strain_criterion can be one of 
        - 'proximity_ranking'
        - 'min_strain'
        - 'clarke_wright'
        each implementing a sligtly different insertion criteria
        - 'gaskell'
        - 'augumented_min_strain'
        try several values and
        - 'all' (default)
        tries all of the above
        
    For 'clarke_wright' and 'gaskell' and when lambda is 2.0 in
    'augumented_min_strain' the routes are initialized according to the 
    primary value of the current strain criteria. For 'proximity_ranking',
    'min_strain', and rest of 'augumented_min_strain' the emerging route
    is initialized with farthest unrouted customer.
    
    Mole, R. and Jameson, S. (1976). A sequential route-building algorithm 
      employing a generalised savings criterion. Journal of the Operational
      ResearchSociety, 27(2):503-511.
    """

    callback_configurations = []
    if strain_criterion=='proximity_ranking' or strain_criterion=='all':
        callback_configurations.append( 
            ( _create_new_criteria_function(lm=0.0,mm=0.0), "farthest") )
        
    if strain_criterion=='min_strain':
       #or strain_criterion=='all': # <- this is already  is in 'augumented_min_strain'
        callback_configurations.append(
            (_create_new_criteria_function(lm=0.0,mm=1.0), "farthest") )
        
    if strain_criterion=='clarke_wright':
       #or strain_criterion=='all': # <- this is already in 'augumented_min_strain'
        callback_configurations.append( 
            # when mu = lambda-1,  initiate route with savings criteria
            (_create_new_criteria_function(lm=2.0,mm=1.0), "strain") )
        
    if strain_criterion=='gaskell' or strain_criterion=='all':
        lambda_mults = [1.25, 1.5, 1.75, 2.0]
        for glm in lambda_mults:
            callback_configurations.append( 
                # when mu = lambda-1,  initiate route with savings criteria
                (_create_new_criteria_function(lm=glm,mm=glm-1.0), "strain") )
            
    if strain_criterion=='augumented_min_strain' or strain_criterion=='all':
        # the lm=2.0, mm=1.0 is already in 'gaskell'
        lambda_mults = [0, 0.5, 1, 1.5] if strain_criterion=='all' else\
                       [0, 0.5, 1, 1.5, 2.0]
        for alm in lambda_mults:
            callback_configurations.append( 
                ( _create_new_criteria_function(lm=alm,mm=1.0),
                   "strain" if  alm-1.0==1.0 else "farthest" ) )
            
    ## Find the best solution among the active strain criterions        
    best_sol = None
    best_f = None
    best_K = None
    interrupted = False
    for strain_function, init_method in callback_configurations:
        
        sol, sol_f, sol_K = None, float('inf'), float('inf')
        try:
            sol = cheapest_insertion_init(D, d, C, L, minimize_K=False,
                                          emerging_route_count=1,
                                          initialize_routes_with=init_method,
                                          insertion_strain_callback=strain_function,
                                          insert_callback=_try_insert_2opt_and_update)
            sol = _refine_solution(sol, D, d, C, L, minimize_K)
            # LS may make some of the routes empty
            sol = without_empty_routes(sol)
            
        except KeyboardInterrupt as e: #or SIGINT
            # some of the strain function insertion runs was interrupted
            if len(e.args)>0 and type(e.args[0]) is list:
                sol = e.args[0]
                interrupted = True
           
        if sol:
            sol_f = objf(sol, D)
            sol_K = sol.count(0)-1
        if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
            best_sol = sol
            best_f = sol_f
            best_K = sol_K
            
        if interrupted:
            raise KeyboardInterrupt(best_sol)
            
    return best_sol

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
# (takes an extra argument for enabling parallel version)
def get_mj_algorithm():
    algo_name = "MJ76-INS"
    algo_desc = "Mole & Jameson (1976) sequential cheapest insertion "+\
                "heuristic with a route improvement phase"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if single:
            return mole_jameson_insertion_init(D, d, C, L, minimize_K,
                                               strain_criterion="clarke_wright")
        else:
            return mole_jameson_insertion_init(D, d, C, L, minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_mj_algorithm())
        
