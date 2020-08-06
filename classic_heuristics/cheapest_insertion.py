#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides a basis for implementing of the sequential and 
parallel cheapest insertion heuristics (see e.g. Mole and Jameson (1976)).

The script is callable and can be used as a standalone sequental insertion
solver for TSPLIB formatted CVRPs. It has moderate dependencies as it requires
3rs party module dllist, in addition to numpy and scipy for reading and
preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from collections import namedtuple
from heapq import heappush, heappop
from logging import log, DEBUG

from util import objf, routes2sol
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

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


class _EmergingRouteData:
    """ A data structure to keep all of the data of the route that is being
    constructed in the same place. """
    def __init__(self, seed_customers, D, d):
        self.potential_insertions = []
        if seed_customers:
            lroute = [0]+seed_customers+[0]
            self.route = dllist(lroute)
            self.cost = objf(lroute, D)
            self.used_capacity = sum(d[n] for n in lroute) if d else 0
        else:
            self.route = dllist([0,0])
            self.used_capacity = 0
            self.cost = 0
            
    @staticmethod
    def export_solution(route_datas):
        sol = []
        for rd in route_datas:
            if rd and len(rd.route)>2:
                sol += list(rd.route)[1:]
        return sol

def _initialize_new_route(seed_customers, unrouted, D, d):
    chosen = None
    while True:
        _, candidates = seed_customers.pop()
        chosen = []
        for c in candidates:
            if c not in unrouted:
                chosen = []
                break
            if c not in chosen:
                chosen.append(c)
        if not chosen:
            continue        
        
        # remove duplicates
        for c in chosen:
            unrouted.remove(c)
        break
    
    if (chosen is None) and (len(unrouted)>0):
        chosen = [unrouted.pop()]
        
    return _EmergingRouteData(chosen, D, d)


Insertion = namedtuple('Insertion', ['customer', 'after_node', 'before_node',
                                     'cost_delta', 'demand_delta'])
def _new_potential_insertions(unrouted, D, d, C, L, rd, after, before,
                              insertion_strain_callback):
    for customer in unrouted:
        # do not add infeasible insertions (note that the insertions, however,
        #  can become infeasible later!)
        l_delta =  -D[after.value, before.value]\
                   +D[after.value, customer]\
                   +D[customer, before.value]
        
        d_delta = 0.0
        if C:
            d_delta = d[customer]
            if rd.used_capacity+d_delta-C_EPS>C:
                continue
        if L:
            if rd.cost+l_delta-S_EPS > L:                
                continue
        
        # Build the insertion for the heap
        insertion = Insertion(customer, after, before, l_delta, d_delta)
        strain, secondary_criterion = insertion_strain_callback(D, after.value,customer,before.value)
        l_saving =  l_delta-D[0,customer]-D[customer,0]
        
        sorting_key = (-strain, secondary_criterion, l_saving)
        #sorting_key = (-strain, l_saving, secondary_criterion)
        #sorting_key = (-strain, +D[after.value, customer]
        #                        +D[customer, before.value], -customer)
        
        heappush( rd.potential_insertions, (sorting_key, id(insertion), insertion) )
      

def _seems_valid_insertion(insertion, route_current_d, unrouted, D, d, C):  
    # The node has been already inserted somewhere
    if not insertion.customer in unrouted:
        if __debug__:
            log(DEBUG-3,"Customer n%d is already routed."%insertion.customer)
        return False
    
    # Something else has already been inserted here, ignore
    if insertion.after_node.next != insertion.before_node:
        if __debug__:
            log(DEBUG-3,"Chain n%d-n%d-n%d is no longer possible."%
                      (insertion.after_node.value, insertion.customer, insertion.before_node.value))
        return False
    
    # Also the available capacity may have changed,
    #  check if inserting would break a constraint
    if C:
        d_inc = d[insertion.customer]
        if route_current_d+d_inc-C_EPS>C:
            if __debug__:
                log(DEBUG-3,"Insertion of n%d would break C constraint."%insertion.customer)
            return False
    
    return True

def _try_insert_and_update(insertion, rd, D, L, minimize_K):
    if not minimize_K:
        # Compared to a solution where the customer is served individually
        insertion_cost = insertion.cost_delta\
                         -D[0,insertion.customer]\
                         -D[insertion.customer,0]
        if insertion_cost > 0:
            if __debug__:
                log(DEBUG-2,("Rejecting insertion of n%d "%insertion.customer)+
                            ("as it is expected make solution worse."))
            return False, None
        
    if L and rd.cost+insertion.cost_delta-S_EPS > L:
        return False, None
    else:
        inserted = rd.route.insert(insertion.customer, insertion.before_node)                
        rd.cost += insertion.cost_delta
        rd.used_capacity += insertion.demand_delta 
        new_edges = [ (insertion.after_node, inserted), (inserted, insertion.before_node) ]
        return True, new_edges
    
def parametrized_insertion_criteria(D, i, u, j, lm, mm):
    """ Mole and Jameson (1976) insertion criteria.
    
    i = insert after this customer / depot
    u = insert this customer 
    j = insert before this customer / depot
    
    * lm = lambda_multiplier (>=0) defines how strongly one prefers routes that
       visit just one customer.
    * mm = mu_multiplier (>=0)  is similar to the lambda of the savings algorithm. 
       It defines how much the edge that needs to be removed to make the  
       insertion weights when calculating the insertion cost.
    
    Note that we do this in one pass with min( -lambda*c_0_u + e(i,u,j) )
    instead of first calculating e(r,u,s) = min(e(i,u,j)) \forall i,j first and
    then computing  max( lambda*c_0_u - e(r,u,s) ) = max(sigma(i,u,s)).
    
    The reason is twofold. First, the best insertion position for u may be
    taken by other customer, which means that we should keep track on where
    to update it. Second the one pass scheme is equivalent, as the customer 
    with the best (maximal) criteria (l*) is always the best among different
    u as the c_0,u is the same \forall u therefore making the min(e(i,u,j)) 
    the decisive factor also in the one pass mode.
    """ 
    MST_c = D[i,u]+D[u,j]-mm*D[i,j]
    MSAV_c = lm*D[0,u]-MST_c
    return MSAV_c, MST_c
    
def cheapest_insertion_init(D, d, C, L=None, minimize_K=False,
        emerging_route_count=1, 
        initialize_routes_with="farthest",
        
        # these callbacks can be used to change how the algorithm works
        insertion_strain_callback=lambda D, i,u,j:\
            parametrized_insertion_criteria(D, i, u, j,1.0,0.0),
        insert_callback=_try_insert_and_update ):
    """ 
    A (Cheapest) Insertion Heuristic.
    
    The basic idea is introduced in Mole and Jameson (1976): calculate
    each cost of inserting an unrouted customer between every possible pair on
    a route. The one with cheapest insertion cost (based on the some criteria)
    is then executed. The process is repeated until all nodes have been
    inserted. The routes are initialized with a visit to the farthest / 
    losest node(s) in respect to the depot.
    
    Basic parameters:
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route cost/length/duration.
    
    Objective parameter:
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In insertion algorithms 
       this is achieved by not trying to fill the routes completely full if 
       an insertion candidate is expected to make the overall cost worse.
       
    Algorithm specific parameters:
    * emerging_route_count is the number of routes to be build in parallel, the 
       default is 1 (the sequential route building heuristic)
    * initialize_routes_with is set to "farthest" if the emerging route is 
       needed to be initialized with the node furthest from depot. Other
       alternatives are: "closest", initialize with the node closest to the
       depot; "strain" initializes the route using the savings function e.g.
       "D[i,j]-D[0,i]-D[j,0]" (the unrouted nodes of the savings are used
       to initialize the route taken). t can also be a function of the form
       `initialize_routes_with(D, unrouted, rd)`, where rd can be None (for
       the initial route).
     
     Parameters for building extensions:
     * insertion_strain_callback calculates the insertion "cost". Must return
        two values. One is the strain, second is the secondary sorting criteria
     * insert_callback can be used to change which operations are done when a 
        new customer is inserted on the route. Can be used e.g. to optimize
        the route after the insertion. Note, that the _EmergingRouteData
        fields must be updated accordingly. Returns a list of edges that are
        new and for which new insertions can be made, and for which the
        insertion strain is calculated with insertion_strain_callback. Check
        _try_insert_and_update for an example.
    
    Therefore:
      \lambda = 2, \mu = 1 | Generalized Clarke and Wright criterion
      \lambda = 0, \mu = 1 | Minimum strain criterion
      1 <= \lambda <= 2, \mu = lambda-1 | Generalized Gaskell criterion
      \lambda = 0, \mu = 0 | Proximity to two nearest neighbors
     
    See (Solomon 1987) for details on the recommendations for values of
    mu_multiplier and lambda_multiplier. The defaults mu=1.0, lambda=2.0 are
    those that Solomon reported worked best in his experiments. Solomon (1987)
    used presented the sequential version of the algorithm (K=1) and also
    initialized the routes with furthest nodes.
    
    Q: Does the heap approach really correspond to the Insertion heuristic?
    The original description the insertion candidate is found in two phases:
        1. find min(c_1)
        2. if min(c_1) is feasible, find max(c_2)
    A: This should be the case \lambda d_0c for the node to insert c, is the
    same for all potential inserts c_1(a,b) \forall adjacent a,b. As it is! Or
    at least I'm 99% certain it is. -JR
    
    
    Mole, R. and Jameson, S. (1976). A sequential route-building algorithm 
      employing a generalised savings criterion. Journal of the Operational
      ResearchSociety, 27(2):503-511.
     
    Solomon, M. M. (1987). Algorithms for the vehicle routing and scheduling
       problems with time window constraints. Operations research, 35(2).
    """
    
    complete_routes = []       
    unrouted = set(range(1,len(D)))
    
    if initialize_routes_with=="farthest" or initialize_routes_with=="closest":
        # build a ordered priority queue of potential route initialization nodes
        seed_customers = list(zip(D[0][1:], [(i,) for i in range(1,len(D))]))
        seed_customers.sort(reverse = initialize_routes_with=="closest")
    elif initialize_routes_with=="strain":
        # use the strain function to calculate the first insertion
        #  for parametrized_insertion_criteria this leads to 
        #  (\lambda-1)c_{oi}-c_{ij}+c_{ij}
        seed_customers = [
                ((insertion_strain_callback(D, 0, i, j)[0], -D[i,j]), (i,j))
                            for i in range(1,len(D)) for j in range(i,len(D))]
        seed_customers.sort(reverse=False)
        #(seed_customers)
    elif callable(initialize_routes_with):
        seed_customers = initialize_routes_with(D, unrouted, None)
    else:
        raise ValueError("Unknown route initialization method '%s'"%initialize_routes_with)
        
    route_datas = [_initialize_new_route(seed_customers, unrouted, D, d)
                   for k in range(emerging_route_count)]
    if __debug__:
        for ri, rd in enumerate(route_datas):
            log(DEBUG, "Initialized a new route #%d %s"%(ri, str(list(rd.route))))
    
    try:
        # while there are nodes to insert
        route_index = -1
        while len(unrouted)>0:
            # Get next route to insert to (esp. in parallel insertion mode)
            route_index+=1
            if route_index == len(route_datas):
                route_index = 0
                
            rd = route_datas[route_index]
            if __debug__:
                log(DEBUG, "Inserting to route #%d %s"%(route_index, str(list(rd.route))))
            
            # Generate the list of insertion candidates just-in-time
            # (this avoids adding unnecessary candidates).
            if len(rd.potential_insertions)==0:
                initial_after = rd.route.first
                while initial_after!=rd.route.last:
                    _new_potential_insertions(unrouted, D, d, C, L, rd,
                                              initial_after, initial_after.next,
                                              insertion_strain_callback)
                    initial_after = initial_after.next
                    
            while len(rd.potential_insertions)>0:
                insertion_succesfull = False
                # unpack a best candidate from the insertion queue
                _, _, insertion = heappop(rd.potential_insertions)
                
                # Check if it *seems* OK to insert
                # - The node has been already inserted somewhere
                # - Something else has already been inserted here, ignore
                # - And check if inserting would break the C constraint
                is_ok =  _seems_valid_insertion(insertion, rd.used_capacity,
                                                unrouted, D, d, C)
                if not is_ok:
                    continue
                
                if __debug__:
                    log(DEBUG-2,"Try insertion %d-%d-%d with l_delta %.2f"% (
                         insertion.after_node.value, insertion.customer,
                         insertion.before_node.value, insertion.cost_delta))
                   
                # It is all OK to insert the customer (check L constraint)?
                insertion_succesfull, new_edges = insert_callback(insertion,
                                                                  rd,D,L,
                                                                  minimize_K)
                    
                if insertion_succesfull:
                    unrouted.remove(insertion.customer)
                    for from_node, to_node in new_edges:
                        # insertion changes the route -> new insertions possible
                        _new_potential_insertions(unrouted, D, d, C, L, rd,
                                                  from_node, to_node,
                                                  insertion_strain_callback)
                
                    if __debug__:
                        log(DEBUG,"Chose to insert n%d resulting in route %s (%.2f)"%
                                   (insertion.customer, str(list(rd.route)), rd.cost))
                    break
                else:
                    if __debug__:
                        log(DEBUG-3,"Insertion of n%d would break L constraint."%
                            insertion.customer)
                    continue
                        
            # if are not able to add any more customers to the route,
            #  so start a new route
            insertions_exhausted = len(rd.potential_insertions)==0
            if insertions_exhausted and len(unrouted)>0:
                if __debug__:
                    log(DEBUG,"Route #%d finished as %s (%.2f)"%
                        (route_index, str(list(rd.route)), rd.cost))
                complete_routes.append(rd)

                # Some seed initializations rely on the state of the completed
                #  routes. Then, we can use initialize_routes_with callback 
                #  to generate more seed(s).
                if not seed_customers:
                    if callable(initialize_routes_with):
                        seed_customers = initialize_routes_with(D, unrouted,
                                                                complete_routes)
                    else:
                        raise StopIteration("Ran out of seed nodes before all"+
                                            "customers were routed")
                
                rd = _initialize_new_route(seed_customers, unrouted, D, d)
                if __debug__:
                    log(DEBUG, "Initialized a new route #%d %s"%
                               (route_index, str(list(rd.route))))
    
                route_datas[route_index] = rd            
    except KeyboardInterrupt:  #or SIGINT
        interrupted_sol = [0]+_EmergingRouteData.export_solution(route_datas)\
                          +_EmergingRouteData.export_solution(complete_routes)
        interrupted_sol.extend( routes2sol([n] for n in unrouted
                                           if n not in interrupted_sol)[1:] )
        raise KeyboardInterrupt(interrupted_sol)
                              
    return [0]+_EmergingRouteData.export_solution(route_datas)\
              +_EmergingRouteData.export_solution(complete_routes)
                


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
# (takes an extra argument for enabling parallel version)
def get_si_algorithm():
    algo_name = "vB94-SI"
    algo_desc = "Mole & Jameson (1976) sequential cheapest insertion heuristic "+\
                "without local search (van Breedam 1994, 2002)"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return cheapest_insertion_init(D, d, C, L=L, minimize_K=minimize_K,
                                       emerging_route_count=1)
    
    return (algo_name, algo_desc, call_init)

def get_pi_algorithm(emerging_route_count="auto"):
    algo_name = "vB94-PI"
    algo_desc = "van Breedam (1994, 2002) parallel insertion heuristic"
    if emerging_route_count=="auto":
        def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
            sol_ci = cheapest_insertion_init(D, d, C, L=L,
                                             minimize_K=minimize_K,
                                             emerging_route_count=1)
            one_eroute_K = sol_ci.count(0)-1
            return cheapest_insertion_init(D, d, C, L=L,
                                           minimize_K=minimize_K,
                                           emerging_route_count=one_eroute_K)
    elif emerging_route_count>1:
        def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
            return cheapest_insertion_init(D, d, C, L=L,
                                           minimize_K=minimize_K,
                                           emerging_route_count=emerging_route_count)
    else:
        raise ValueError("Not a valid parameter value, has to be 'auto' or "+
                         "an integer > 1, (was %s)"%str(emerging_route_count))
            
    return (algo_name, algo_desc, call_init)

if __name__=="__main__":
    from shared_cli import cli
    cli(*get_si_algorithm())
        
