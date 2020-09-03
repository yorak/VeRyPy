#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Foster and Ryan (1976)
Petal algorithm. A set of candidate routes are generated using the Sweep
procedure of sweep.py and the CVRP is solved as set covering problem.

The script is callable and can be used as a standalone solver for TSPLIB
formatted CVRPs. It has extensive dependencies: MIP solver Gurobi, Sweep 
procedure of sweep.py, a TSP solver (the built in local search solver can be
used), and numpy and scipy for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from signal import signal, SIGINT, default_int_handler
from collections import namedtuple
from math import ceil    
from logging import log, DEBUG

from gurobipy import Model, GRB, GurobiError, quicksum

from classic_heuristics.sweep import do_one_sweep, get_sweep_from_cartesian_coordinates
#from local_search import ITEROPT
from tsp_solvers.tsp_solver_fast import solve_tsp_fast as solve_tsp  
#from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp
from util import objf, totald, routes2sol, without_empty_routes, is_better_sol
from routedata import RouteData

from local_search.inter_route_operators import do_redistribute_move
from local_search import LSOPT

from config import MAX_MIP_SOLVER_RUNTIME, MIP_SOLVER_THREADS
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


# These define how the 'auto' setting for discard limit works. Basically if 
#  the longest route is more than the number of customers of the first number, 
#  remove at most as many customers as the second number indicates.
# The limits must be in monotonically increasing order. If the route length is
#  under every length limit the number of removed customers is not limited. 
# (see petal_init for details) 
AUTO_DISCARD_LIMITS = [(10,5), (20,3), (30,2), (50,1)]

# Usually all petals are used in the final pass. However, for large problems
#  the number of petals may become excessively large. This is the magic number
#  up to we always allow use of full (EXTENDED) set of petals. That is, finally
#  use EXTENDED set as long as the total petal count stays under this number. 
ALWAYS_USE_ALL_PETALS_LIMIT = 10000

class PTL_SET:
    RESTRICTED = 1
    REDUCED = 2
    EXTENDED = 3
       
def _regain_feasibility(to_move, r1_delta, r2_delta, rd1_aftermove, ri1, ri2,
                        discard_at_most, d_excess, l_excess,
                        route_datas, D, d, C, L):
    """ From Foster & Ryan (1976) "The result of such a (1 point) move may
    be to cause the receiving route to exceed the mileage or capacity 
    restriction and so it must be permitted to discard deliveries to regain
    feasibility. ... these ... are in turn relocated in the current solution
    schedule without causing further infeasibility and that the total effect 
    of all the relocations must be a reduction in the mileage."
    
    Thus, after the first of the stored 1 point moves causes the solution
    to become infeasible, this function tries to regain the feasibility by
    redistributing some of the customers from the recieving route to other 
    routes so that the moves still produce an improved solution. However,
    to keep the entire operation improving, the redistribution should not
    increase the cost more than delta_slack."""
    
    ## REPLACE MOVE
    # there is no room, have to replace one *OR MORE* nodes
    # We need to check removing all combinations of customers to make room and
    #  store those that still would remove the problem.
    
    route2, r2_l, r2_d, _ = route_datas[ri2]
    assert route2[0]==0 and route2[-1]==0, "All routes must start and end to the depot"
    r2_min_d = min( d[n] for n in route2[1:-1] ) if d else 0
    
    # First route is changed and reserve an empty route to receive the nodes.
    new_empty_rd = RouteData()
    routes_with_slack = [rd for ri, rd in enumerate(route_datas)
                        if (ri!=ri1 and ri!=ri2 and (not C or C-rd.demand+C_EPS>r2_min_d))]\
                        +[rd1_aftermove,new_empty_rd]
    
    # A depth first search (no recursion) for removing JUST enough customers
    if d: 
        stack  = [(i+1, [n], d[n]) for i,n in enumerate(route2[1:-1])]
    else:
        stack =  [(i+1, [n], None) for i,n in enumerate(route2[1:-1])]
        
    improving_rds = []
    while stack: 
        last_n_i, to_remove_ns, to_remove_d = stack.pop()
        
        # We have to recalculate the delta_r2 as removing nodes almost surely
        #  will change the best position to insert to.
        new_route2 = [n for n in route2 if n not in to_remove_ns]
        new_r2_delta = float('inf')
        new_r2_l = 0.0
        best_j = None
        for j in range(1,len(new_route2)):
            insert_after = new_route2[j-1]
            insert_before = new_route2[j]
            new_r2_l+=D[insert_after,insert_before]
                
            delta = +D[insert_after,to_move]\
                    +D[to_move,insert_before]\
                    -D[insert_after,insert_before]
                        
            if delta < new_r2_delta:
                new_r2_delta = delta
                best_j = j
        to_remove_l = (r2_l+r2_delta)-(new_r2_l+new_r2_delta)
        new_route2 = new_route2[:best_j]+[to_move]+new_route2[best_j:]
        
        if (not d_excess or to_remove_d+C_EPS >= d_excess) and\
           (not l_excess or to_remove_l+S_EPS >= l_excess):
            
            # Redistributing ALWAYS increases the cost at least a little,
            #  it makes no sense to even try sometimes.
            # If all of the improvements goes to inserting to_move customer,
            #  there is none left to try to redistribute.
            delta_slack = -(r1_delta+new_r2_delta)
            if delta_slack+S_EPS<0:
                continue
               
            # After removing customers in to_remove_ns route2 becomes feasible,
            #  now we have to check if redistributing the removed back to the
            #  solution is possible.
            
            to_redistribute_r = [0]+to_remove_ns+[0]
            to_redisribute_rd = RouteData( to_redistribute_r,
                                           objf(to_redistribute_r, D),
                                           totald(to_redistribute_r, d) )
            
            result = do_redistribute_move(to_redisribute_rd, routes_with_slack,
                                          D, d, C, L, strategy=LSOPT.BEST_ACCEPT,
                                          recombination_level=0,
                                          best_delta=delta_slack)
            redistribute_delta = result[-1]        
            if redistribute_delta is not None:
                    
                new_r2_d = r2_d-to_remove_d+d[to_move] if d else 0
                new_rd2 = RouteData(new_route2, new_r2_l+new_r2_delta, new_r2_d)
                improving_rds.append(new_rd2) 
                improving_rds += result[1:-1] # includes modified rd1
                
            #TODO: if one would like to explore all discard combinations, 
            # one would do branching also here or at least if the redisribute
            # fails.
            
        elif not discard_at_most or len(to_remove_ns)<discard_at_most:
            # branch out
            for candidate_j, candidate_n in enumerate(route2[last_n_i+1:-1]):
                stack.append((last_n_i+1+candidate_j,
                              to_remove_ns+[candidate_n],
                              to_remove_d+d[candidate_n] if d else None))
                
    return improving_rds
                
def _generate_solution_relaxation_petals(
        routes, discard_at_most,
        D, C=None, d=None, L=None):
    
    """ This is the route improvement, or overconstraint i&ii relaxations
    phase, where customers of a LP solution can be moved from route to another
    if the move improves the solution.
    
    Note that this differs from the do_1point_move local_search because the
    prodecure allows regaining feasibility if the move would improve the 
    solution (see _regain_feasibility for details)."""
    
    relaxed_route_datas = []
    route_datas = RouteData.from_routes(routes, D, d)
    route_indices = list(range(len(route_datas)))
    
    # Try to move a customer from route 1 to route 2
    for ri1 in route_indices:        
        route1, r1_l, r1_d, _ = route_datas[ri1]  

        i = 0
        while i<len(route1)-2:
            i+=1
            
            remove_after = route1[i-1]
            to_move = route1[i]
            remove_before = route1[i+1]
            
            r1_delta =  +D[remove_after,remove_before]\
                        -D[remove_after,to_move]\
                        -D[to_move,remove_before]
        
            for ri2 in route_indices:
                # Do not try to move nodes from ri1 to ri1!
                if ri1==ri2:
                    continue
                
                route2, r2_l, r2_d, _ = route_datas[ri2]                
                         
                # Find the best place to insert it to
                j = 0
                best_r2_delta = float('inf')
                best_j = None
                while j<len(route2)-1:
                    j+=1
                
                    insert_after = route2[j-1]
                    insert_before = route2[j]
                    
                    r2_delta = +D[insert_after,to_move]\
                                +D[to_move,insert_before]\
                                -D[insert_after,insert_before]
                                
                    if r2_delta<best_r2_delta:
                        best_r2_delta = r2_delta
                        best_j = j 
    
                        
                if r1_delta+best_r2_delta<-S_EPS:
                    
                    
                    d_excess = r2_d+d[to_move]-C \
                               if (C and r2_d+d[to_move]-C_EPS>C) \
                               else None
                    l_excess = r2_l+best_r2_delta-L \
                               if (L and r2_l+best_r2_delta-S_EPS>L) \
                               else None

                    new_route1 = route1[:i]+route1[i+1:]
                    new_rd1 = RouteData(new_route1,
                                        r1_l+r1_delta,
                                        r1_d-d[to_move] if d else 0)
                    new_route2 = route2[:best_j]+[to_move]+route2[best_j:]
                    new_rd2 = RouteData(new_route2,
                                        r2_l+best_r2_delta,
                                        r2_d+d[to_move] if d else 0)
                    
                    # If operation would break a constraint...
                    if d_excess or l_excess:
                        # ... discard and redistribute some nodes
                        mod_route_datas = _regain_feasibility( to_move,
                                r1_delta, best_r2_delta, new_rd1, ri1, ri2,
                                discard_at_most,
                                d_excess, l_excess,
                                route_datas, D, d,C, L)
                    else:
                        mod_route_datas = [new_rd1, new_rd2]
                        
                    relaxed_route_datas.extend(mod_route_datas)
                    
                        
    return relaxed_route_datas 

PetalSet = namedtuple('PetalSet', ['nodes', 'routes', 'costs'])
def _generate_overconstrained_petals(
                     restricted_route_ratio, C_t, 
                     N, points, D, d, C, L,
                     generate_resticted_ptls = True,
                     generate_reduced_ptls = True,
                     generate_extended_ptls  = True):
    
    gen_ptls = set()
    can_gen_restricted = restricted_route_ratio is not None
    
    restricted_ptls = PetalSet([],[],[]) if generate_resticted_ptls else None
    reduced_ptls = PetalSet([],[],[]) if generate_reduced_ptls else None
    extended_ptls = PetalSet([],[],[]) if generate_extended_ptls else None
    
    sweep = get_sweep_from_cartesian_coordinates(points) 
    for direction in [1]:#, -1]:
        for start in range(0,N-1):
            # Use Sweep to generate petals.
            route_datas = do_one_sweep(N, D, d, C, L, solve_tsp, 
                                       sweep, start, direction,
                                       generate_alternative_first_routes = True)
                  
            for rd in route_datas:
                rd.update_node_set()
                rd.normalize()
                tuple_r = tuple(rd.route) # convert to tuple for hashing
                if __debug__:
                    if tuple_r in gen_ptls: 
                        log(DEBUG-2, "Petal collision with %s"%str(rd.route))
                        log(DEBUG-2, "(route generated twice)")
                if tuple_r not in gen_ptls:
                    if can_gen_restricted and (
                       (C and (rd.demand+C_EPS>=restricted_route_ratio*C))):# or
                       #(L and (rd.cost+S_EPS>=restricted_route_ratio*L))):
                        if not generate_resticted_ptls:
                            continue
                        restricted_ptls.routes.append(rd.route)
                        restricted_ptls.costs.append(rd.cost)
                        restricted_ptls.nodes.append(rd.node_set)
                        if __debug__: petal_class = "resticted"
                    elif rd.demand+C_EPS>=C_t:
                        if not generate_reduced_ptls:
                            continue
                        reduced_ptls.routes.append(rd.route)
                        reduced_ptls.costs.append(rd.cost)
                        reduced_ptls.nodes.append(rd.node_set)
                        if __debug__: petal_class = "reduced"                        
                    else:
                        if not generate_extended_ptls:
                            continue
                        extended_ptls.routes.append(rd.route)
                        extended_ptls.costs.append(rd.cost)
                        extended_ptls.nodes.append(rd.node_set)                        
                        if __debug__: petal_class = "extended"
                    gen_ptls.add(tuple_r)
                    
                    if __debug__:
                        log(DEBUG-2, "Storing %s route %s (%.2f) as a petal"%
                                     (petal_class, str(rd.route), rd.cost))
                
    if __debug__:
        r_cnt = len(restricted_ptls.routes)
        e_cnt = len(reduced_ptls.routes)
        x_cnt = len(extended_ptls.routes)
        log(DEBUG, "#PETALS (RESTRICTED) = %d (cum. total %d)"%(r_cnt,r_cnt))
        log(DEBUG, "#PETALS (REDUCED) = %d (cum. total %d)"%(e_cnt,r_cnt+e_cnt))
        log(DEBUG, "#PETALS (COMPLETE) = %d (cum. total %d)"%(x_cnt,r_cnt+e_cnt+x_cnt))
        
        all_ptls = restricted_ptls.routes+reduced_ptls.routes+extended_ptls.routes
        all_costs = restricted_ptls.costs+reduced_ptls.costs+extended_ptls.costs
        for i, petal in enumerate(all_ptls):
            log(DEBUG-3, "PETAL #%d : %s (%.2f)"%(i, str(petal), all_costs[i]))
    
    return restricted_ptls, reduced_ptls, extended_ptls

def _decision_variables_to_petals(X, active_ptls, relaxed_ptls):
    nactive = len(active_ptls.routes)
    nrelaxed = len(relaxed_ptls.routes)
    routes_with_idxs = []
    
    for j, v in enumerate(X):
        if j>=nactive+nrelaxed:
            # relaxed solution has additional decision variables for choosing
            #  which constraints to relax, skip them
            break
        if v:
            if j<nactive:
                r_w_idx = (list(active_ptls.routes[j]), j)
            else:
                relaxed_j = j-nactive
                r_w_idx = (list(relaxed_ptls.routes[relaxed_j]),
                           -relaxed_j-1)
            routes_with_idxs.append(r_w_idx)
    return routes_with_idxs

def _add_set_covering_constraints(m, N, K, X_j, X_j_keys, X_j_node_sets, X_j_forbidden_combos):
    # c1, the convering constraint, each node mus be served exactly by 
    #  one route
    c1_constrs = []
    for i in range(1,N):
        active_vars_sum = quicksum([X_j[j]
                                   for j, ns in enumerate(X_j_node_sets)
                                   if i in ns])
        c1c = m.addConstr( active_vars_sum == 1, "c1_n%d"%i )
        c1_constrs.append(c1c);
        
    # c2, the forbid constraints that do not allow some petal combinations
    c2_constrs = []
    for i, fc in enumerate(X_j_forbidden_combos):
        # sum_{i \in S } x_i - sum_{j \in \not{S} } x_j < |S| 
        combo_sum = quicksum([X_j[j] if (j in fc) else -X_j[j]
                                    for j in X_j_keys])
        c2c = m.addConstr( combo_sum  <= len(fc)-1, "c2_n%d"%i )
        c2_constrs.append(c2c);
    
    # c3, number of vehicles constraint
    c3_constrs = []
    if K:
        active_vars_sum = quicksum(X_j)
        c3c = m.addConstr( active_vars_sum <= K, "c3" )
        c3_constrs.append(c3c)
        
    return c1_constrs, c2_constrs, c3_constrs
        
        
def _relax_customer_constraints_with_feasRelax(m, customer_cover_constraints, active_ptls, relaxed_ptls):
    # relax the model and relax minimal number of customer serving constraints
    pens = [1.0]*len(customer_cover_constraints)
    m.feasRelax(2, True, None, None, None, customer_cover_constraints, pens)
    # TODO: not sure if feasRelax can change Status, test it someday
    if m.Status == GRB.INTERRUPTED: 
        raise KeyboardInterrupt()
    m.optimize()

    # restore SIGINT callback handler which is changed by gurobipy
    signal(SIGINT, default_int_handler)

    status = m.Status
    if __debug__:
        log(DEBUG-2, "Relaxed problem Gurobi runtime = %.2f"%m.Runtime)
        if m.Status == GRB.OPTIMAL:
            log(DEBUG-3, "Relaxed problem Gurobi objective = %.2f"%
                         m.getObjective().getValue())

    if status == GRB.OPTIMAL:
        return _decision_variables_to_petals(m.X, active_ptls, relaxed_ptls), False
    elif status == GRB.TIME_LIMIT:
        raise GurobiError(10023, "Gurobi timeout reached when attempting to solve relaxed SCPCVRP")
    elif m.Status == GRB.INTERRUPTED:
        raise KeyboardInterrupt()
    return None,False
    
def _solve_set_covering(N, active_ptls, relaxed_ptls, forbidden_combos,
                        allow_infeasible=True, D=None, K=None):
    """A helper function that solves a set covering problem. The inputs are
    a list of node sets and corresponding costs for the set.
    --
    Fisher, M. L. and Jaikumar, R. (1981), A generalized assignment 
    heuristic for vehicle routing. Networks, 11: 109-124.
    """
    
    m = Model("SCPCVRP")
    nactive = len(active_ptls.routes)
    nrelaxed = len(relaxed_ptls.routes)
    nforbidden = len(forbidden_combos)
    
    # the order of the keys is important when we interpret the results
    X_j_keys = list(range(nactive+nrelaxed))
    # variables and the objective
    X_j_costs = active_ptls.costs+relaxed_ptls.costs
    X_j_node_sets = active_ptls.nodes+relaxed_ptls.nodes
    #update forbidden indices to match the current active petal set
    X_j_forbidden_combos = []
    for fc in forbidden_combos:
        if all( i<nactive for i in fc ):
            X_j_forbidden_combos.append( [i if i>=0 else -i-1+nactive
                                          for i in fc] ) 
        
    #print("REMOVEME: Solving with K=%d, %d node sets, and %d solutions forbidden" % (K, nactive+nrelaxed,nforbidden))
    
    if __debug__:
        log(DEBUG, "Solving over-constrained VRP as a set covering problem "+
              "with %d petals, where %d of the possible configurations are forbidden." % 
              (nactive+nrelaxed,nforbidden))
        if nforbidden>0:
            log(DEBUG-2," and with following solutions forbidden:")
            log(DEBUG-3,"(petal indices = %s)"%str(X_j_forbidden_combos))
            for fc in X_j_forbidden_combos:
                fc_sol = [0]
                for i in fc:
                    if i<nactive:
                        fc_sol.extend( active_ptls.routes[i][1:] )
                    else:
                        fc_sol.extend( relaxed_ptls.routes[i-nactive][1:] )
                fc_sol = without_empty_routes(fc_sol)
                log(DEBUG-2, "%s (%.2f)"%(fc_sol , objf(fc_sol,D)))
    
    X_j = m.addVars(X_j_keys, obj=X_j_costs , vtype=GRB.BINARY, name='x')    
    
    ## constraints
    c1_constrs, c2_constrs, c3_constrs = _add_set_covering_constraints(m, N, K,
        X_j, X_j_keys, X_j_node_sets, X_j_forbidden_combos)
    
    ## update the model and solve 
    m._vars = X_j  
    m.modelSense = GRB.MINIMIZE
    m.update()
    # disable output
    m.setParam('OutputFlag', 0)    
    m.setParam('TimeLimit', MAX_MIP_SOLVER_RUNTIME)
    m.setParam('Threads', MIP_SOLVER_THREADS)
    #m.write("petalout.lp")
    m.optimize()

    # restore SIGINT callback handler which is changed by gurobipy
    signal(SIGINT, default_int_handler)

        
    if __debug__:
        log(DEBUG-2, "Gurobi runtime = %.2f"%m.Runtime)
        if m.Status == GRB.OPTIMAL:
            log(DEBUG-3, "Gurobi objective = %.2f"%m.getObjective().getValue())
    
    if m.Status == GRB.OPTIMAL:
        return _decision_variables_to_petals(m.X, active_ptls, relaxed_ptls), True
    elif m.Status == GRB.TIME_LIMIT:
        raise GurobiError(10023, "Gurobi timeout reached when attempting to solve SCPCVRP")
    elif m.Status == GRB.INTERRUPTED:
        raise KeyboardInterrupt()
    # Sometimes the solution is infeasible, try to relax it a little.
    elif m.Status == GRB.INFEASIBLE and allow_infeasible:
        return _relax_customer_constraints_with_feasRelax(m, c1_constrs, active_ptls, relaxed_ptls)
    return None,False

def _log_debug_scp_info(ptl_set, nptl, routes_with_idxs, is_feasible, active_K):
    set_name = ""
    if ptl_set==PTL_SET.RESTRICTED:
        set_name = "RESTRICTED"
    if ptl_set==PTL_SET.REDUCED:
        set_name = "REDUCED"
    if ptl_set==PTL_SET.EXTENDED:
        set_name = "EXTENDED"
    feasible_str = "INFEASIBLE"
    if is_feasible:
        feasible_str = "FEASIBLE"
        
    log(DEBUG-2,"Solved set covering with #petals (%s set) = %d"%(set_name,nptl))
    if routes_with_idxs:
        routes, _ = zip(*routes_with_idxs)
        if active_K:
            log(DEBUG-2,"It produced %s solution with K=%d/%d"%(feasible_str, len(routes), active_K))
        else:
            log(DEBUG-2,"It produced %s solution with K=%d"%(feasible_str, len(routes)))
    else:
        log(DEBUG-2,"It produced NO solution") 

def _remove_multiserved(chosen_ptl_routes, D):
    """If the customer serve constraints are lifted, it is possible that 
    a customer is (unncessarily) served twice, fix the solution."""
    served=[[] for i in range(len(D))]
    served_too_many_times = []
    nserved = 1 #depot served always
    for ri, r in enumerate(chosen_ptl_routes):
        for ni, n in enumerate(r):
            if n!=0:
                served[n].append( (ri,ni) )
                if len(served[n])==1:
                    nserved+=1
                elif len(served[n])==2:
                    served_too_many_times.append(served[n])
                    
    for si, sl in enumerate(served_too_many_times):
        # served more than once, leave the one that makes the sol have the best
        #  total cost.
        best_delta = float('-inf')
        keep_this = None
        for ri,ni in sl:
            nfrom = chosen_ptl_routes[ri][ni-1]
            nvia = chosen_ptl_routes[ri][ni]
            nto = chosen_ptl_routes[ri][ni+1]
            delta = D[nfrom,nto]-D[nfrom,nvia]-D[nvia,nto]
            if delta>best_delta:
                best_delta = delta
                keep_this = ri,ni
        if keep_this:
            kri, kni = keep_this
            for rri,rni in sl:
                if kri==rri and kni==rni:
                    continue
                else:
                    del chosen_ptl_routes[rri][rni]
                # update indices
                for slu in served_too_many_times[si+1:]:
                    for islu, pos in enumerate(slu):
                        riu,niu = pos
                        if riu==rri and niu>rni:
                            slu[islu] = (riu,niu-1)
    return nserved==len(D)

def petal_init(points, D, d, C, L, K=None,
               minimize_K = False, 
               relaxe_SCP_solutions=True,
               required_iterations=None,
               min_iterations=None,
               restricted_route_ratio=0.75,
               allow_infeasible = True,
               can_discard_multiple_customers='auto',
               predefined_petals_generator=None):
    
    """ An implementation of Foster and Ryan (1976) Petal algorithm. The VRP
    is solved with a set covering formulation (SCP->MIP/LP). The decision 
    variables are feasible routes (petals) R_i. The three initial petal sets
    (RESTRICTED, REDUCED, and EXTENDED) are generated with a Sweep algorithm
    the RELAXED petal set is grown when improving local search moves are found.
    
    The algorithm solves SCP iteratively and the petal set size is increased
    when no improvements are found or other conditions are not fulfilled:
        RESTRICTED+RELAXED ->
        REDUCED+RESTRICTED+RELAXED ->
        EXTENDED+REDUCED+RESTRICTED+RELAXED
    
    The implementation uses Gurobi to solve the SCP. 
    
    The RESTRICTED set contains routes each with total demand:
        d_R_i > restricted_route_ratio*C
    
    The REDUCED set contains routes not in RESTRICTED but with total demand:
        d_R_i > C_l
        C_l = sum_{j}{d_j}-(K-1)C
        
    The EXTENDED set contains the ones that have not high enough total demand
    to be included in RESTRICTED or REDUCED sets. By default, the algorithm 
    terminates as soon as no improving solutions can be found, even with
    increasing the petal set size (RESTRICTED->REDUCED->EXTENDED) or loosening
    the K constraint.
    
    Args:
    * points is a list (or 2D numpy ndarray) of coordinate points.
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.

    * K can be set to require a predefined number of vehicles. However, this 
       constraint is loosened if no feasible solutions are found. If minimize_K 
       is set to True and K is None, a minimum K is found by starting from sm-
       allest possible K and increasing it until a feasible solution is found.
       Else, if K set manually, it is respected (and minimize_K has no effect).
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In practice disables
       automatic K constraint.
    
    * relaxe_SCP_solutions turns on the improvement heuristic that is similar 
       to Pair (the chain) heuristic of Wren & Holliday (1972). These improved
       routes are stored to the RELAXED petal set.
    * required_iterations can be used force the Petal algorithm to run a preset 
       number of iterations, with each  exploring one valid SCP solution. If
       this is set to -1 (default) the algorithm terminates after no SCP
       improvements seem to be possible. required_iterations can change this
       by allowing new RELAXED petals to be generated from worse SCP solutions.
       If set to None (default) the feature is disabled.
    * min_iterations works similarly to required_iterations, but it only
       specifies a lower bound for the number of iterations. The rule based
       termination can still decide to continue the search after min_iterations
       has been tried. If set to None (default) the feature is disabled.
    * can_discard_multiple_customers when building relaxed petals, one can 
       set how many customers can be redistributed. This is needed when the 
       route length grows, and the number of combinations becomes too large.
       The value can  be set to:
         False = discard at most 1
         True = discard as many as necessary
         <int> = discard at most this many
         'auto' = the number of maximum discarded customers depend on the
                  length of the route as specified by AUTO_DISCARD_LIMITS
                  e.g. at most 5 for routes with more than 10 customers, 
                       at most 3 for routes with more than 20 customers, 
                       at most 2 for routes with more than 30 customers, 
                       at most 1 for routes with more than 50 customers, 
      
    Foster, B. A. and Ryan, D. M. (1976). An integer programming approach to
    the vehicle scheduling problem. JORS, 27(2):367-384.
    
    Wren, A., & Holliday, A. (1972). Computer scheduling of vehicles from one
    or more depots to a number of delivery points. JORS, 23(3), 333-344.
    """
    
    if not points:
        raise ValueError("The algorithm requires 2D coordinates for the points")
    N = len(points)
    
    # Translate the argument to a value
    if type(can_discard_multiple_customers) is int:
        discard_at_most = can_discard_multiple_customers
    elif can_discard_multiple_customers==False:
        discard_at_most = 1
    elif can_discard_multiple_customers==True:
        discard_at_most = None # no limit
    elif can_discard_multiple_customers=='auto':
        pass # is set when _generate_solution_relaxation_petals is called
    else:
        raise ValueError("Unknown parameter value for can_discard_multiple_customers"+
                         " can be integer, True, False or 'auto'")
        
        
    
    ## 1. GENERATE PETALS
    
    # "Route capacity Lower Bound" (LB) C_t of Foster and Ryan (1976, p. 373)
    if C:
        # "smallest integer **greater** than x" (note, not equal)
        d_tot = sum(d)
        v = ceil(d_tot/float(C)+C_EPS)
        K_constraint = v if (not K and minimize_K) else K
        C_t = d_tot-(K-1)*C if K else d_tot-(v-1)*C
    else:
        K_constraint = 1 if (not K and minimize_K) else K
        C_t = 0.0
    
    #TODO: for the larger instances, get the restricted FIRST, then the
    # reduced and the extended ONLY IF ABSOLUTELY NEEDED!
    if predefined_petals_generator is None:
        restricted_ptls, reduced_ptls, extended_ptls = \
            _generate_overconstrained_petals(
                restricted_route_ratio, C_t,
                N, points, D, d, C, L)
    else:
         restricted_ptls, reduced_ptls, extended_ptls = predefined_petals_generator()

    relaxed_ptls = PetalSet([],[],[]) 
    active_ptls = None
        
    r_cnt = len(restricted_ptls.routes)
    e_cnt = len(reduced_ptls.routes)
    x_cnt = len(extended_ptls.routes)
    ptl_cnt = r_cnt+e_cnt+x_cnt
    
    ptl_set = PTL_SET.RESTRICTED
    max_ptl_set_used = PTL_SET.RESTRICTED
    seen_ptls = set()
    ptls_in_use_set = None
        
    ## 2. SOLVE SET COVERING UNTIL NO IMPROVEMENTS    
    routes_with_idxs = None    
    forbidden_combinations = []
    best_sol_feasible = False
    best_sol = None
    best_sol_f = float('inf')
    best_sol_K = len(D)
    interrupted = False
    
    #TODO: the set covering solving with reduced, then extended (if no 
    # fesible set covering problem, all petals if still no feasible sol.),
    # and on subsequent iterations relaxed tepalts, are all very similar.
    # Thus, one should use warm starting and add/modify the constraints on the 
    # fly. Now the LP/MIP is solved from the beginning on each iteration.
    # Check how this column generation can be implemented in Gurobi.

    iteration_counter = 0
    while not interrupted:
        iteration_counter+=1

        ## SET COVERING PHASE
        
        n_base_ptls = 0
        n_rlxd_ptls = len(relaxed_ptls.routes)
        routes_with_idxs = None
        sol_feasible = False
        
        try:
            if ptl_set==PTL_SET.RESTRICTED and len(restricted_ptls.nodes)>0:
                # Try set covering first with a reduced petal set, that are, by 
                #   default, at least 75% full.
                
                # Update the set of active petals
                if ptls_in_use_set!=ptl_set:
                    seen_ptls.update(set(tuple(r) for r in restricted_ptls.routes))
                    active_ptls = PetalSet(restricted_ptls.nodes,
                                           restricted_ptls.routes,
                                           restricted_ptls.costs)
                    ptls_in_use_set = PTL_SET.RESTRICTED
                
                n_base_ptls = len(active_ptls.routes)
                routes_with_idxs, sol_feasible = _solve_set_covering(N,
                    active_ptls, relaxed_ptls, forbidden_combinations,
                    allow_infeasible, D=D, K=K_constraint)
                
                if __debug__:
                    active_K = K_constraint if K_constraint else K
                    _log_debug_scp_info(ptl_set, n_base_ptls+n_rlxd_ptls,
                                        routes_with_idxs, sol_feasible, active_K)
                    
                    
            no_rsol = not routes_with_idxs
            if (no_rsol and ptl_set==PTL_SET.RESTRICTED) or ptl_set==PTL_SET.REDUCED:
                # The reduced set contains all feasible routes that hit the 
                #  constraint limit C_t = \sum(d_i)-(K-1)*C
                
                ptl_set = PTL_SET.REDUCED
                # Update the set of active petals
                if ptls_in_use_set!=ptl_set:
                    seen_ptls.update(set(tuple(r) for r in reduced_ptls.routes))
                    active_ptls = PetalSet(
                        restricted_ptls.nodes+reduced_ptls.nodes,
                        restricted_ptls.routes+reduced_ptls.routes,
                        restricted_ptls.costs+reduced_ptls.costs)
                    ptls_in_use_set = PTL_SET.REDUCED
                    
                n_base_ptls = len(active_ptls.routes)
                routes_with_idxs, sol_feasible = _solve_set_covering(N,
                    active_ptls, relaxed_ptls, forbidden_combinations,
                    allow_infeasible, D=D, K=K_constraint)
                  
                if __debug__:
                    active_K = K_constraint if K_constraint else K
                    _log_debug_scp_info(ptl_set, n_base_ptls+n_rlxd_ptls,
                                        routes_with_idxs, sol_feasible, active_K)
    
            # "If the LP defined by this reduced petal set uses more than v 
            #  vehicles we can computethe remaining petals with capacitiesless
            #  than c1 and confirm the LP solution on the extended petal set."
            #                                            - Foster & Ryan 1976
            no_esol = not routes_with_idxs
            if no_esol or (K and len(routes_with_idxs[0])>K) or ptl_set==PTL_SET.EXTENDED:
                # The complete set contains all feasible routes that were generated
                #  with the sweep procedure including single customer routes.
                
                ptl_set = PTL_SET.EXTENDED
                # Update the set of active petals
                if ptls_in_use_set!=ptl_set:
                    seen_ptls.update(set(tuple(r) for r in extended_ptls.routes))
                    active_ptls = PetalSet(
                        restricted_ptls.nodes+reduced_ptls.nodes+extended_ptls.nodes,
                        restricted_ptls.routes+reduced_ptls.routes+extended_ptls.routes,
                        restricted_ptls.costs+reduced_ptls.costs+extended_ptls.costs)
                    ptls_in_use_set = PTL_SET.EXTENDED
                    
                n_base_ptls = len(active_ptls.routes)
                routes_with_idxs, sol_feasible = _solve_set_covering(N,
                    active_ptls, relaxed_ptls, forbidden_combinations,
                    allow_infeasible, D=D, K=K_constraint)
                    
                if __debug__:
                    active_K = K_constraint if K_constraint else K
                    _log_debug_scp_info(ptl_set, n_base_ptls+n_rlxd_ptls, routes_with_idxs, sol_feasible, active_K)
                    
            max_ptl_set_used = max(max_ptl_set_used, ptl_set)
            
        except KeyboardInterrupt: #or SIGINT
            # continue and store routes_with_idxs if any
            interrupted = True

        try:
            if routes_with_idxs:
                chosen_ptl_routes, chosen_plt_indices = zip(*routes_with_idxs)
            
                #TODO: in some pathological situations (probaby due to floating point 
                # math inaccuracies) feasRelax may return infeasible solution that does 
                # not really respect the forbidden solution constraints. It seems this
                # is only in cases where there actually is no solution (e.g. due to K
                # constraint and forcing it just breaks things).
                #WARNING: This is just a quickfix, should investigate this "someday".
                if chosen_plt_indices in forbidden_combinations:
                    found_solution = False
                else:
                    found_solution = True
                    if not sol_feasible and not interrupted:
                        if __debug__:
                            log(DEBUG-2, "Check if some customers are served "+
                                         "multiple times and remove those "+
                                         "visits that are unnecessary")
                        sol_feasible = _remove_multiserved(chosen_ptl_routes, D)
                    else:
                        # Once the first feasible solution is found, we no longer 
                        #  accept infeasible ones!
                        allow_infeasible = False
            else:
                found_solution = False
            
        
            # It may be the case that constraints are so tight that we did not find 
            #  even an infeasbile solution. Do our best to find one by loosening K.
            if not found_solution:
                if ptl_set == PTL_SET.RESTRICTED:
                    ptl_set = PTL_SET.REDUCED
                elif ptl_set == PTL_SET.REDUCED:
                    ptl_set = PTL_SET.EXTENDED
                elif K_constraint:
                    ptl_set = PTL_SET.RESTRICTED
                    K_constraint = K_constraint+1
                else:
                    return best_sol
                continue
            
            petal_sol = routes2sol(chosen_ptl_routes)
            petal_sol_f = objf(petal_sol, D)
            petal_sol_K = petal_sol.count(0)-1
            if __debug__:
                feasibilitys = "feasible" if len(set(petal_sol))==len(D) else "infeasible"
                log(DEBUG, "\nGot %s Petal LP sol %s (%.2f)"%
                           (feasibilitys, str(petal_sol), petal_sol_f))
                log(DEBUG-2, "\n...with indices %s"%str(chosen_plt_indices))
            
            #print("REMOVEME: it =", iteration_counter, "f =", petal_sol_f,
            #      "k =", len(chosen_ptl_routes), "is feasible =", sol_feasible) 
        
            # "a new starting schedule is determined by banning all the routes in
            # the optimum IP solution and re-converging the IP"
            # add in the relaxed petals as negative numbers
            forbidden_combinations.append(chosen_plt_indices)
              
            # Check if we still continue: for example, check if an improvement was
            #  made and store the best  so far. The condition is a little tricky as
            #  feasible wins infeasible, always.
            
            better_or_same = is_better_sol(best_sol_f, best_sol_K, petal_sol_f,
                                           petal_sol_K, minimize_K)\
                             or\
                             (petal_sol_K==best_sol_K and petal_sol_f==best_sol_f)
            
            if (best_sol_feasible and sol_feasible and better_or_same) or\
               (not best_sol_feasible and (sol_feasible or better_or_same)):
                best_sol = petal_sol
                best_sol_f = petal_sol_f
                best_sol_K = petal_sol_K
                best_sol_feasible = sol_feasible
            elif ptl_set==PTL_SET.RESTRICTED:
                # "The region is relaxed to the complete petal set when no 
                #  further improvements can be found from the reduced set" 
                #                                       -Foster & Ryan 1976
                ptl_set=PTL_SET.REDUCED
                
            ## Try our best to find a feasible solution
            
            # If we were unable to find a feasible solution from REDUCED set,
            #  try with complete set.
            # Also, if we have used the EXTENDED set prior to relaxing K, we grow
            #  the petal set back to it's former size after checking REDUCED set.
            elif ((not best_sol_feasible or ptl_cnt<=ALWAYS_USE_ALL_PETALS_LIMIT)
                   and (ptl_set==PTL_SET.REDUCED)) or (ptl_set<max_ptl_set_used):
                if __debug__:
                    if (not best_sol_feasible):
                        log(DEBUG, "No feasible solution found, using COMPLETE petal set.")
                ptl_set=PTL_SET.EXTENDED
            elif (not best_sol_feasible) and K_constraint:
                if __debug__:
                    log(DEBUG, "No feasible solution found, loosening K constraint.")
                # As a last resort, relax K
                K_constraint = K_constraint+1
                ptl_set=PTL_SET.RESTRICTED
            elif not (min_iterations and iteration_counter<min_iterations) and\
                (required_iterations is None):
                # Finally, abort if no improvements can be found.
                break
                
            if required_iterations and iteration_counter>=required_iterations and\
                not (min_iterations and iteration_counter<min_iterations):        
                break # main iteration loop
        
            ## IMPROVEMENT PHASE 
            # aka. "relaxation ... of moving one delivery  between two routes so
            #  as to reduce the total mileage"
            
            if not interrupted and relaxe_SCP_solutions:
                if __debug__:
                    log(DEBUG, "Searching for relaxed petals between the "+
                        "%d routes of the Petal solution"%(petal_sol.count(0)-1))
                
                if can_discard_multiple_customers=='auto':
                    # Avoid combinatorial explosion and restrict how many customers
                    #  can be removed and redistributed.
                    discard_at_most = None
                    longest_route_len = max( len(r) for r in chosen_ptl_routes)-2
                    for limit, set_to in AUTO_DISCARD_LIMITS:
                        if longest_route_len > limit:
                            discard_at_most = set_to 
                        
                new_petal_candidates = _generate_solution_relaxation_petals(
                        chosen_ptl_routes, discard_at_most,
                        D, C, d, L)  
                
                #TODO: Implement the second secondary relaxation that tries to find
                # improvements where groups of concecutive customers are moved from
                # a route to another (p.381, Foster & Ryan 1976) and add an option
                # to enable it.
                
                # Store petals that are new
                for ird in new_petal_candidates:
                    ird.update_node_set()
                    ird.normalize()
                    tuple_r = tuple(ird.route)
                    if len(tuple_r)>2 and tuple_r not in seen_ptls:
                        relaxed_ptls.routes.append( ird.route )
                        relaxed_ptls.costs.append( ird.cost ) 
                        relaxed_ptls.nodes.append( ird.node_set ) 
                        if __debug__:
                            log(DEBUG-1, "Added a relaxed petal %s (%.2f)"%(str(ird.route), ird.cost))
                        seen_ptls.add(tuple_r)
                        ptl_cnt+=1
                
                # Do not store the improved solution as the best petal solution,
                #  as the next set conver solution will cover (pun intented) this.
        except KeyboardInterrupt: #or SIGINT
            interrupted = True
            break #the main loop
                
    if interrupted:
        raise KeyboardInterrupt(best_sol)
        
    return best_sol 
        
# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ptl_algorithm():
    algo_name = "FR76-1PTL"
    algo_desc = "Foster & Ryan (1976) Petal set covering algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        if single:
            return petal_init(points, D,d,C,L,
                              minimize_K=minimize_K,
                              required_iterations=1,
                              relaxe_SCP_solutions=False)
            
        else:
            return petal_init(points, D,d,C,L,
                              minimize_K=minimize_K)
                              #minimize_K=True)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ptl_algorithm())
