# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and implements heuristics that operate on the entire
solution. Also, all operators have the following signature:
    
do_Z_move(solution, D,C,d,L, strategy, best_delta), where

* solution is a giant tour encoded solution that must start and end to the
    depot (index 0) that has also visits to the depot between routes e.g.
    [0,1,2,3,0,4,5,0].
* D is the symmetric full distance matrix, given as numpy-like 2darray.
    Thus, D[i,j] gives the distance between nodes i and j. 
* C is the capacity constraint, i.e. the maximum carrying capacity of vehicles
* d is the demand of each customer (demand for depot d[0] is 0)
* L is the maximum route duration/length/cost constraint.
* strategy is either FIRST_ACCEPT (default) or BEST_ACCEPT. 
    First accept returns a modified route as soon as the first improvement is  
    encoutered. Best accept tries all possible combinations and returns the
    best one.
* best_delta is the required level of change in the in route cost. It can be
   used to set the upper bound (requirement) for the improvement. Usually it 
   is None, which sets the level to 0.0. In that case only improving deltas are
   accepted. If set to (large) positive value, best worsening move can also be
   returned.
   
All operators return the new improved solution and the improvement (delta) as a
2-tuple or (None,None) if no improvement was found."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from builtins import range

from collections import namedtuple
from local_search import LSOPT

from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

DEPOT = 0

REMOVE_ME_DEBUG = False

def _build_cumulative_lists(solution, D, d, direction):
    """ This is a helper function that is used to build auxiliary data for a
    solution. The auxiliary data has cumulative demands and route costs.
    
    * Direction should be 1 or -1 
    
    Note: depot node has always the cumulative sums set to 0.0
    (no matter the direction) """
    
    sol_cum_d = [0.0]*len(solution)
    sol_cum_l = [0.0]*len(solution)
    
    route_cum_d = 0.0
    route_cum_l = 0.0
    prev_n = 0
    start_i = 1 if direction==1 else -2
    i = start_i
    for n in solution[start_i::direction]:   
        if n==0:
            route_cum_d = 0.0
            route_cum_l = 0.0
        else:
            if d:
                route_cum_d += d[n]
            route_cum_l += D[prev_n,n]
            
        sol_cum_d[i] = route_cum_d
        sol_cum_l[i] = route_cum_l
        
        prev_n = n
        i+=direction
    
    if REMOVE_ME_DEBUG:
        print("\t".join(str(n) for n in solution))
        print("\t".join(str(v) for v in sol_cum_d))
        print("\t".join(str(v) for v in sol_cum_l))
        print()
        
    return sol_cum_d, sol_cum_l

def _check_3opt_move(D, C, L, removed_weights, best_delta,
                     edges, end_p, end_n, cum_d, cum_l,
                     ldepot_12, ldepot_34):
    """ Each move is checked in in 2 if clauses to avoid unnecessary 
      feasibility checking (the clauses "short-circuit", that is, 
      when the first clause fails, the rest are not checked).
     1. check if the move would improve
     2. if loops are formed, there is a depot in the segment
     3. capacity constraint is not violated
     4. max. route cost constraint is not violated
    """
    
    added_weights = D[end_n[edges[0][0]], end_n[edges[0][1]]]+\
                    D[end_n[edges[1][0]], end_n[edges[1][1]]]+\
                    D[end_n[edges[2][0]], end_n[edges[2][1]]]
    delta = added_weights-removed_weights

    # Would be an improvement and the move will not create 
    #  a loop route disconnected from a depot?
    if (delta+S_EPS<best_delta):
       
        # Constraint checks
        feasible = True
        # pylint: disable=unsubscriptable-object
        prev_edge = None
        route_d = 0
        route_l = 0
        
        # Assume that the edges are in right order and to the right direction
        for curr_edge in edges:
            # Check if there is a visit to the depot on the previous segment, 
            #  that is, between previous END and the current edge START nodes.
            if (curr_edge[0]==0 or prev_edge[1]==5 or
                (ldepot_12 is not None and (
                 end_p[prev_edge[1]]<=ldepot_12<=end_p[curr_edge[0]] or
                 end_p[prev_edge[1]]>=ldepot_12>=end_p[curr_edge[0]])) or 
                ( ldepot_34 is not None and (      
                end_p[prev_edge[1]]<=ldepot_34<=end_p[curr_edge[0]] or
                end_p[prev_edge[1]]>=ldepot_34>=end_p[curr_edge[0]]))):
               
                if C:
                    route_d = cum_d[curr_edge[0]] 
                if L:
                    route_l = cum_l[curr_edge[0]]
                
            if REMOVE_ME_DEBUG:
                print("edge", "abcdef"[curr_edge[0]], "abcdef"[curr_edge[1]],\
                      "positions", end_p[curr_edge[0]], end_p[curr_edge[1]],\
                      "nodes", end_n[curr_edge[0]], end_n[curr_edge[1]])
                print("PRE_CHECK d,l", route_d, route_l)
            if C:
                route_d += cum_d[curr_edge[1]]
                if REMOVE_ME_DEBUG:
                    print("POST_CHECK d", route_d)
                if route_d-C_EPS>C:
                    feasible = False
                    break # edge loop
            if L:
                e_wt = D[end_n[curr_edge[0]],end_n[curr_edge[1]]]
                route_l+=e_wt+cum_l[curr_edge[1]]
                if REMOVE_ME_DEBUG:
                    print("POST_CHECK l", route_l)
                if route_l-S_EPS>L:
                    feasible = False
                    break # edge loop
                   
            prev_edge = curr_edge
               
        # store best feasible move
        if REMOVE_ME_DEBUG:
            print("FEASIBLE\n" if feasible else "INFEASIBLE\n")
        if feasible:
            return delta
        
    # not an improving feasible move
    return None
                               
SolutionAuxiliaryData = namedtuple('SolutionAuxiliaryData',
    ['sol','fwd_d','rwd_d','fwd_l','rwd_l'])    

class MoveDefinition:
    def __init__(self, move_idx, edges, 
                 loop_12,loop_34,
                 concat_recipe):
      
        self.move_idx = move_idx
        self.edges = edges
        self.loop_12 = loop_12
        self.loop_34 = loop_34
        self.concat_recipe = concat_recipe
 
MOVES_3OPTSTAR = (
    # ->- = arbitarily many nodes and connecting edges
    # + = a recombined edge
    # (d) = the depot
    # (0-5) = a node (can be a depot)
    # (6-9/d) = a node next to a depot on a segment (if applicable)
    #
    # (d)->-(0)+(1)->-{(d)-(6)}->-{(d)-(7)}->-(2)
    #                                          +
    # (d)-<-(5)+(4)-<-{(9)-(d)}-<-{(8)-(d)}-<-(3)
    #
    # this is the initial state:
    #((0, 1), (2, 3), (4, 5)), 
    
    # 2-OPT MOVES, NO LOOPS
    
    # #a->b# c->e d->f ( 2-opt move )
    #  0->1  2->4 3->5
    #               _____
    #              /     \
    # >-a--b->-c  d-<-e  f-->
    #          \_____/  
    MoveDefinition( move_idx=0,
        edges=((0, 1), (2, 4), (3, 5)),
        loop_12=False, loop_34=False, 
        concat_recipe=((None,3,1), (4,2,-1), (5,None,1)) ),
    
   
    # a->e #d->c# b->f  (2-opt move)
    # 1->4  3->2  1->5
    #         ____________
    #        /            \
    # >--a  b-<-c--d-<-e  f-->
    #     \___________/ 
    MoveDefinition( move_idx=1,
        edges=((0, 4), (3, 2), (1, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1), (4,0,-1), (5,None,1))),
            
    # a->c b->d #e->f#  ( 2-opt move )
    # 0->2 1->3  4->5
    #        _____
    #       /     \
    # >--a  b-<-c  d->-e--f-->
    #     \____/  
    MoveDefinition( move_idx=2,
        edges=((0, 2), (1, 3), (4, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1), (2,0,-1), (3,None,1))),
    

    ## 3-OPT MOVES, NO LOOPS
            
    # a->c b->e d->f  ( 3-opt move )
    # 0->2 1->4 3->5
    #         _________
    #        /         \
    # >--a  b-<-c  d-<-e  f-->
    #     \____/    \____/
    MoveDefinition( move_idx=3,
        edges=((0, 2), (1, 4), (3, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1),(2,0,-1),(4,2,-1),(5,None,1))),
            
    # a->d e->b c->f  (3-opt move)
    # 0->3 4->1 2->5
    #         __________
    #        /          \
    # >--a  b->-c   d->-e  f-->
    #     \______\_/      /
    #            \_______/
    MoveDefinition( move_idx=4,
        edges=((0, 3), (4, 1), (2, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1),(3,5,1),(1,3,1),(5,None,1))),
            
    # a->d e->c b->f  (3-opt move)
    # 0->3 4->2 1->5
    #          ___________
    #         /   _____   \
    #        /   /     \  \  
    # >--a  b-<-c  d->-e  f-->
    #     \_______/  
    MoveDefinition( move_idx=5,
        edges=((0, 3), (4, 2), (1, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1),(3,5,1),(2,0,-1),(5,None,1))),
         
    # a->e d->b c->f  (3-opt move)
    # 0->4 3->1 2->5
    #             ________
    #            /        \
    # >--a  b->-c  d-<-e  f-->
    #     \  \____/   / 
    #     \__________/  
    MoveDefinition( move_idx=6,
        edges=((0, 4), (3, 1), (2, 5)),
        loop_12=False, loop_34=False,
        concat_recipe=((None,1,1),(4,2,-1),(1,3,1),(5,None,1))),
            
    ## 2-OPT, WITH LOOPS ON 1-2
    
    # a->d c->b #e->f#  ( 2-opt move, with a loop )
    # 0->3 4->5  2->1
    #        _______
    #       /       \
    # >--a  b->-0->-c  d->-e--f-->
    #     \___________/  
    MoveDefinition( move_idx=7,
        edges=((0, 3), (4, 5), (2, 1)),
        loop_12=True, loop_34=False,
        concat_recipe=((None,1,1), (3,None,1), (7,3,1), (1,7,1) )),
                   
                   
    # a->e d->f c->b  (3-opt move, with a loop)
    # 0->4 3->5 2->1
    #         ______     _____
    #        /      \   /     \
    # >--a  b->-0->-c  d-<-e  f-->
    #    \________________/  
    MoveDefinition( move_idx=8,
        edges=((0, 4), (3, 5), (2, 1)),
        loop_12=True, loop_34=False,
        concat_recipe=((None,1,1),(4,2,-1),(5,None,1),(7,3,1),(1,7,1))),
                   
    # TODO: This move can be evaluated twice, once for each loop option. Thus
    #  there is some overhead. Consider a faster (and more complex) way.
    #  Remember that the concat_recipe must be defined for both options.
    # a->f #c->d# e->b  (2-opt move, with a loop)
    # 0->5  2->3  4->1
    #        __________________
    #       /                  \
    # >--a  b->-0->-c--d->-0->-e  f-->
    #     \_______________________/ 
    MoveDefinition( move_idx=9,
        edges=((0, 5), (2, 3), (4, 1)),
        loop_12=True, loop_34=False,
        concat_recipe=((None,1,1), (5,None,1), (7,5,1), (1,7,1))),
                   
    # TODO: This move can be evaluated twice, once for each loop option. Thus
    #  there is some overhead. Consider a faster (and more complex) way.
    #  Remember that the concat_recipe must be defined for both options. 
    # a->f c->e d->b  (3-opt move, with a loop)
    # 0->5 2->4 3->1
    #         _________    
    #        /         \   
    # >--a  b->-0->-c  d-<-0-<-e  f-->
    #    \          \_________/  /
    #    \______________________/  
    MoveDefinition( move_idx=10,
        edges=((0, 5), (1, 3), (4, 2)),
        loop_12=True, loop_34=False,
        concat_recipe=((None,1,1),(5,None,1),
                       (7,3,1),(4,2,-1),(1,7,1))),
    
    ## 2-OPT, WITH LOOPS ON 3-4
    
    # #a->b# c->f e->d ( 2-opt move, with a loop )
    #  0->1  2->5 4->3
    #               ______
    #              /      \
    # >-a--b->-c  d->-0->-e  f-->
    #          \_____________/  
    MoveDefinition( move_idx=11,
        edges=((0, 1), (2, 5), (4, 3)),
        loop_12=False, loop_34=True, 
        concat_recipe=((None,3,1), (5,None,1), (9,5,1), (3,9,1))),
    
    # a->c b->f e->d  ( 3-opt move, with a loop )
    # 0->2 1->5 4->3
    # 
    #         _______________
    #        /               \
    # >--a  b-<-c  d->-0->e  f-->
    #     \____/    \____/
    MoveDefinition( move_idx=12,
        edges=((0, 2), (1, 5), (4, 3)),
        loop_12=False, loop_34=True,
        concat_recipe=((None,1,1),(2,0,-1),(5,None,1),(9,5,1),(3,9,1))),
            
    # TODO: This move can be evaluated twice, once for each loop option. Thus
    #  there is some overhead. Consider a faster (and more complex) way.
    #  Remember that the concat_recipe must be defined for both options.
    # a->f e->b #c->d#  (2-opt move, with a loop)
    # 0->5 4->1  2->3  
    #        __________________
    #       /                  \
    # >--a  b->-0->-c--d->-0->-e  f-->
    #     \_______________________/ 
    MoveDefinition( move_idx=13,
        edges=((0, 5), (4, 1), (2, 3)),
        loop_12=False, loop_34=True,
        concat_recipe=((None,1,1), (5,None,1), (9,5,1), (1,9,1))),
                       
    # TODO: This move can be evaluated twice, once for each loop option. Thus
    #  there is some overhead. Consider a faster (and more complex) way.
    #  Remember that the concat_recipe must be defined for both options.
    # a->f e->c b->d  (3-opt move, with a loop)
    # 0->5 4->2 1->3
    #         _________    
    #        /         \   
    # >--a  b-<-0-<-c  d->-0->-e  f-->
    #    \          \_________/  /
    #    \______________________/  
    MoveDefinition( move_idx=14,
        edges=((0, 5), (4, 2), (1, 3)),
        loop_12=False, loop_34=True,
        concat_recipe=((None,1,1),(5,None,1),
                       (9,5,1),(2,0,-1),(3,9,1))),
                   
    # a->f c->b e->d  (3-opt move, with a loop)
    # 0->5 2->1 4->3
    #         ______     ______
    #        /      \   /      \
    # >--a  b->-0->-c  d->-0->-e  f-->
    #    \_______________________/  
    MoveDefinition( move_idx=15,
        edges=((0, 5), (2, 1), (4, 3)),
        loop_12=True, loop_34=True,
        concat_recipe=((None,1,1),(5,None,1),
                       (7,3,1),(1,7,1),
                       (9,5,1),(3,9,1))),
)

# These predefined move sets are chosen depending on the loops of the i,j,k    
MOVES_3OPTSTAR_FOR_ALL = MOVES_3OPTSTAR[:7]
MOVES_3OPTSTAR_WHEN_12LOOP = MOVES_3OPTSTAR[:11]
MOVES_3OPTSTAR_WHEN_34LOOP = MOVES_3OPTSTAR[:7]+MOVES_3OPTSTAR[11:15]
MOVES_3OPTSTAR_WHEN_BOTH_LOOPS = MOVES_3OPTSTAR

def do_3optstar_move(solution, D, demands=None, 
                     C=None, L=None, # constraint
                     strategy=LSOPT.FIRST_ACCEPT, 
                     best_delta = None,
                     move_checker = _check_3opt_move,
                     return_solution_with_auxiliary_data=False):
    """ 3-opt local search operation for the symmetric distances D that 
    operates on the entire solution.
    
    Note: due to how this is implemented, it may change the route order."""

    sN = len(solution)
    best_move = None
    if best_delta is None:
        best_delta = 0
    accept_move = False
    
    # make sure we have the auxlirary data pre-computed for constant time
    #  feasibility checks
    if isinstance(solution, list):
        sol_cum_d_fwd, sol_cum_l_fwd = \
            _build_cumulative_lists(solution, D, demands, +1)
        sol_cum_d_rwd, sol_cum_l_rwd = \
            _build_cumulative_lists(solution, D, demands, -1)
        sol_data = SolutionAuxiliaryData(solution,
                    sol_cum_d_fwd, sol_cum_d_rwd,
                    sol_cum_l_fwd, sol_cum_l_rwd)
    # something else than a list? Rely on duck typing.
    else:
        sol_data = solution
        solution = sol_data.sol
    
    
    # cached data tables containing:
    # edge end positions, end nodes, and cumulative demand and length (cost)
    #  at looking at the segments from the detached edge ends.
    end_p = [0]*6
    end_n = [0]*6
    cum_d = [0]*6
    cum_l = [0]*6
    
    for i in range(0,sN-1):
        # left and rightmost visits to the depot between i and j
        ldepot_12 = None
        rdepot_12 = None
        
        # update the cached data tables 
        end_p[0] = i
        end_p[1] = i+1
        end_n[0] = solution[i]
        end_n[1] = solution[i+1]
             
        
        for j in range(i+1,sN-1):
            if solution[j]==DEPOT:
                if ldepot_12 is None:
                    ldepot_12 = j
                rdepot_12 = j
        
            # left and rightmost visits to the depot between j and k
            ldepot_34 = None
            rdepot_34 = None
            
            # update the cached data tables 
            end_p[2] = j
            end_p[3] = j+1
            end_n[2] = solution[j]
            end_n[3] = solution[j+1]
            if C:
                cum_d[0] = sol_data.fwd_d[i]
                cum_d[1] = sol_data.rwd_d[i+1]\
                 -(sol_data.rwd_d[j+1] if (ldepot_12 is None) else 0.0)
            if L:
                cum_l[0] = sol_data.fwd_l[i]                    
                cum_l[1] = sol_data.rwd_l[i+1]\
                 -(sol_data.rwd_l[j]   if (ldepot_12 is None) else 0.0) 
                 
            for k in range(j+1,sN-1):
                if solution[k]==DEPOT:
                    if ldepot_34 is None:
                        ldepot_34 = k
                    rdepot_34 = k

                # update the cached data tables 
                end_p[4] = k
                end_p[5] = k+1
                end_n[4] = solution[k]
                end_n[5] = solution[k+1]
                if C:
                    cum_d[2] = sol_data.fwd_d[j]\
                     -(sol_data.fwd_d[i]   if (ldepot_12 is None) else 0.0)
                    cum_d[3] = sol_data.rwd_d[j+1]\
                     -(sol_data.rwd_d[k+1] if (ldepot_34 is None) else 0.0)
                    cum_d[4] = sol_data.fwd_d[k]\
                     -(sol_data.fwd_d[j]   if (ldepot_34 is None) else 0.0)
                    cum_d[5] = sol_data.rwd_d[k+1]         
                if L:
                    cum_l[2] = sol_data.fwd_l[j]\
                     -(sol_data.fwd_l[i+1] if (ldepot_12 is None) else 0.0)
                    cum_l[3] = sol_data.rwd_l[j+1]\
                     -(sol_data.rwd_l[k]   if (ldepot_34 is None) else 0.0)  
                    cum_l[4] = sol_data.fwd_l[k]\
                     -(sol_data.fwd_l[j+1] if (ldepot_34 is None) else 0.0)
                    cum_l[5] = sol_data.rwd_l[k+1]                    
                
                # the cumulative costs and demands
                # 0--a->b-0?-c->d-0?-e->f--0
                #TODO: write a test for this
                # _update_segment_cumulative_costs can be found in r8092
                #cum_d, cum_l = _update_segment_cumulative_costs(
                #        i,j,k,C, L, sol_data, ldepot_12, ldepot_34)
                #
                #if REMOVE_ME_DEBUG:
                #    print("\ncum_d: "+"\t".join(str(n) for n in cum_d))
                #    print("cum_l: "+"\t".join(str(n) for n in cum_l))
                
                # After removing edges a-b, c-d, and e-f, try all of the 7
                # combinations in which the segments can be reconnected.
                removed_weights = \
                    D[sol_data.sol[end_p[0]],sol_data.sol[end_p[1]]]+\
                    D[sol_data.sol[end_p[2]],sol_data.sol[end_p[3]]]+\
                    D[sol_data.sol[end_p[4]],sol_data.sol[end_p[5]]]

                ## MOVES
                
                # The moves could be checked in a loop, however, because
                #  this is EXTREMELY HOT piece of code, only check the
                #  moves that are currently available due to how loops
                #  can be formed when edges are added back.
                
                depot_between_12 = (ldepot_12 is not None)
                depot_between_34 = (ldepot_34 is not None)
                
                if depot_between_12 and depot_between_34:
                    moves = MOVES_3OPTSTAR_WHEN_BOTH_LOOPS
                elif depot_between_12:
                    moves = MOVES_3OPTSTAR_WHEN_12LOOP
                elif depot_between_34:
                    moves = MOVES_3OPTSTAR_WHEN_34LOOP
                else:
                    moves = MOVES_3OPTSTAR_FOR_ALL
                    
                for move in moves:
                    
                    improvement_delta = move_checker(
                        D, C, L, removed_weights, best_delta,
                        move.edges, end_p, end_n, cum_d, cum_l,
                        ldepot_12, ldepot_34)
                    
                    if improvement_delta is not None:
                        best_move = (move.move_idx, improvement_delta,
                                  (i,i+1,j,j+1,k,k+1,
                                  ldepot_12+1 if depot_between_12 else None,
                                  rdepot_12+1 if depot_between_12 else None,
                                  ldepot_34+1 if depot_between_34 else None,
                                  rdepot_34+1 if depot_between_34 else None))
                        best_delta = improvement_delta
                        
                        if strategy==LSOPT.FIRST_ACCEPT:
                            accept_move = True
                            break # move loop
                if accept_move:
                    break # k-loop
            if accept_move:
                break # j-loop
        if accept_move:
            break # i-loop
        
    if best_move:
        best_move_idx, best_delta, cut_positions = best_move
        concat_recipe = MOVES_3OPTSTAR[best_move_idx].concat_recipe
        new_solution = []
        for segm in concat_recipe:
            from_pos = None if (segm[0] is None) else cut_positions[segm[0]]
            to_pos = None if (segm[1] is None) else cut_positions[segm[1]]
            direction = segm[2]
            
            #TODO: avoid reversing sections between depots!
            # (the depot locations are in cut_positions)
            
            # Avoid concatenating visits to the depot (one less route)
            if (new_solution and new_solution[-1]==DEPOT) and \
               ((from_pos is not None) and solution[from_pos]==DEPOT):
                from_pos+=direction
            
            new_solution.extend( solution[from_pos:to_pos:direction] )
            
        if return_solution_with_auxiliary_data:
            raise NotImplementedError("Changing and returning the modified aux data is not implemented, please recalculate it instead.")
        else:
            return new_solution, best_delta

    return None, None

