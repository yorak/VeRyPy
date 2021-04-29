# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and implements intra route (i.e. within one route) local 
search improvement heuristics such as 2-opt, one-point-move etc. All operators 
assume we start from a feasible solution. Also, all functions implementing the
operations have the following signature:

do_X_move(route, D, strategy, best_delta), where

* route is a list of nodes forming a route.The route must start and end to the
    depot (index 0), e.g. [0,1,2,3,0].
* D is the symmetric full distance matrix, given as numpy-like 2darray.
    Thus, D[i,j] gives the distance between nodes i and j. 
* strategy is either FIRST_ACCEPT (default) or BEST_ACCEPT. 
    First accept returns a modified route as soon as the first improvement is  
    encoutered. Best accept tries all possible combinations and returns the
    best one.
* best_delta is the required level of change in the in route cost. It can be
   used to  set the upper bound (requirement) for the improvement. Usually it 
   is None, which sets the level to 0.0. In that case only improving deltas are
   accepted. If set to (large) positive value, best worsening move can also be
   returned.
   
All intra route improvement operators return the new improved route and the
improvement (delta) as 2-tuple or (None,None) if no improvement was found."""
###############################################################################

#TODO:
# - support for asymmetric distances, the loops in 2-opt and 3-opt could keep
#    track of the mid segments cost in reverse order with neglible performance
#    impact (no added loops, just update the reverse segment cost on each iter).
# - improve performance
#  - use doubly linked list (dllist) as the data structure for the routes
#  - calculate nearest neighbour list of the route nodes and iterate trough it
#     instead of iterating over entire route when finding potential moves
#      (this is probably the single greatest pefrormance improvement that can
#       be made but requires some brain time to do right)
#  - do not copy routes, make changes in place
#  - use numba or even cython (and numba compatible dllist)

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from local_search import LSOPT
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def do_2opt_move(route, D, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """ 2-opt local search operation for the symmetric distances D
    Remove 2 edges from the route and check if swapping then endpoints
    would yield an improvement.
    """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(0,rN-1):
        for j in range(i+1,rN-1):
            a = route[i]
            b = route[i+1]
            c = route[j]
            d = route[j+1]
            
            # a->c b->d( 2-opt move )
            #        _____
            #       /     \
            # >--a  b-<-c  -->d
            #     \____/ 
            
            # For edges a-b and c-d, try if savings could be made by traveling 
            #  from a-c and b-d (effectively reversing the chain from
            #  b to c).
            delta = D[a,c] + D[b,d] \
                     -D[a,b]-D[c,d]
                     
            if delta+S_EPS<best_delta:
                best_move = (i,j)
                best_delta = delta
                if strategy==LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break # j loop
        if accept_move:
            break # i loop
                
    if best_move:
        i,j = best_move
        #print("REMOVEME:","best_move", i,j, route, best_delta)
        return route[:i+1]+route[j:i:-1]+route[j+1:], best_delta
    return None, None
    
def do_3opt_move(route, D, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """ 3-opt local search operation for the symmetric distances D """
 
    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    
    for i in range(0,rN-1):
        for j in range(i+1,rN-1):
            for k in range(j+1,rN-1):

                # the edge endpoints
                a = route[i]
                b = route[i+1]
                c = route[j]
                d = route[j+1]
                e = route[k]
                f = route[k+1]
                #print("search abcdef", a,b,c,d,e,f)
                
                # After removing edges a-b, c-d, and e-f, try all of the 7
                # combinations in which the segments can be reconnected.
                removed_weights = D[a,b] + D[c,d] + D[e,f]
                
                
                # The combinations could be iterated, but for simplicity
                #  (and speed), the loop has been unrolled below.
                
                ## 2-opt moves
                
                # #a->b# c->e d->f ( 2-opt move )
                #               _____
                #              /     \
                # >-a--b->-c  d-<-e  f-->
                #          \_____/  
                delta=(D[a,b] + D[c,e] + D[d,f])-removed_weights   
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (i+1,j+1, 1),
                                 (k,j, -1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                # a->c b->d #e->f#  ( 2-opt move )
                #        _____
                #       /     \
                # >--a  b-<-c  d->-e--f-->
                #     \____/ 
                delta==(D[a,c] + D[b,d] + D[e,f])-removed_weights
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (j,i, -1),
                                 (j+1,k+1, 1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                # a->e #d->c# b->f  (2-opt move)
                #        ____________
                #       /            \
                # >--a  b-<-c--d-<-e  f-->
                #     \___________/ 
                delta=(D[a,e] + D[d,c] + D[b,f])-removed_weights   
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (k,j, -1),
                                 (j,i, -1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                ## 3-opt moves
                
                # a->c b->e d->f  ( 3-opt move )
                #         _________
                #        /         \
                # >--a  b-<-c  d-<-e  f-->
                #     \____/    \____/
                delta=(D[a,c] + D[b,e] + D[d,f])-removed_weights;    
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (j,i, -1),
                                 (k,j, -1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                # a->d e->b c->f  (3-opt move)
                #         __________
                #        /          \
                # >--a  b->-c   d->-e  f-->
                #     \______\_/      /
                #             \______/
                delta=(D[a,d] + D[e,b] + D[c,f])-removed_weights
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (j+1,k+1, 1),
                                 (i+1,j+1, 1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                # a->d e->c b->f  (3-opt move)
                #          __________
                #         /   _____  \
                #        /   /     \  \  
                # >--a  b-<-c  d->-e  f-->
                #     \_______/ 
                delta=(D[a,d] + D[e,c] + D[b,f])-removed_weights
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (j+1,k+1, 1),
                                 (j,i, -1), (k+1,None, 1))
                    best_delta = delta
                    
                    if strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
                    
                    
                # a->e b->d c->f  (3-opt move)
                #             _______
                #            /       \
                # >--a  b->-c  d-<-e  f-->
                #     \  \____/   / 
                #      \_________/  
                delta=(D[a,e] + D[d,b] + D[c,f])-removed_weights
                if delta+S_EPS<best_delta:
                    best_move = ((None,i+1, 1), (k,j, -1),
                                 (i+1,j+1, 1), (k+1,None, 1))
                    best_delta = delta
                    
                    if best_move and strategy==LSOPT.FIRST_ACCEPT:
                        accept_move = True
                        break # k loop
            if accept_move:
                break # j loop
        if accept_move:
            break # i loop
        
    if best_move:
        sgmt1,sgmt2,sgmt3,sgmt4 = best_move
        
        # the head containing leaving the depot +
        #  the first segment, reversed or not +
        #  the second segment, reversed or not +
        #  the tail containing return to the depot
        return route[sgmt1[0]:sgmt1[1]:sgmt1[2]]+\
               route[sgmt2[0]:sgmt2[1]:sgmt2[2]]+\
               route[sgmt3[0]:sgmt3[1]:sgmt3[2]]+\
               route[sgmt4[0]:sgmt4[1]:sgmt4[2]], best_delta

    return None, None

def do_relocate_move(route, D, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """Relocate local search operation for the symmetric distances D.
    Check if a node on the route can be moved to another position on the same
    route.
    
    Please note that relocate search space is is a subset of 3-opt. 
    However, the search space of 3-opt is larger.
    """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(1,rN-1):
        for j in range(1,rN):
            if i==j or j==i-1:
                continue
            
            a = route[i-1]
            b = route[i]
            c = route[i+1]
            
            d = route[j-1]
            e = route[j]
            
            # check the no op position
            if d==b or e==b:
                continue
            
            # move b from between a and c to between d and e
            delta = -D[a,b]-D[b,c]+D[a,c]\
                     -D[d,e]+D[d,b]+D[b,e]
                                          
            if delta+S_EPS<best_delta:
                best_move = (i,j)
                best_delta = delta
                if strategy==LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break # j loop                
        if accept_move:
            break # i loop
                
    if best_move:
        i,j = best_move
        if i<j:
            return route[:i]+route[i+1:j]+[route[i]]+route[j:], best_delta
        else:
            return route[:j]+[route[i]]+route[j:i]+route[i+1:], best_delta
            
    return None, None

def do_exchange_move(route, D, strategy=LSOPT.FIRST_ACCEPT, best_delta=None):
    """(Node) exchange local search operation for the symmetric distances D.
    Checks if two nodes on the route can be swapped. 
    
    Please note that exchange search space is is a subset of 4-opt. 
    However, the search space of 4-opt is significantly larger.
    """

    rN = len(route)
    best_move = None
    if not best_delta:
        best_delta = 0
    accept_move = False
    for i in range(1,rN-1):
        for j in range(i+1,rN-1):
            if i==j:
                continue
            a = route[i-1]
            b = route[i]
            c = route[i+1]
            
            d = route[j-1]
            e = route[j]
            f = route[j+1]
            
            if c==e:
                delta = -D[a,b]-D[b,e]-D[e,f]\
                         +D[a,e]+D[e,b]+D[b,f]
            else:
                # swap b and e from between a and c to between d and f
                delta = -D[a,b]-D[b,c]+D[a,e]+D[e,c]\
                         -D[d,e]-D[e,f]+D[d,b]+D[b,f]

            #print("REMOVEME:", i,j, "(", delta, ")", "best =", best_delta)
            if delta+S_EPS<best_delta:
                best_move = (i,j)
                best_delta = delta
                
                if strategy==LSOPT.FIRST_ACCEPT:
                    accept_move = True
                    break # j loop
        if accept_move:
            break # i loop
                
    if best_move:
        i,j = best_move
        return route[:i]+[route[j]]+route[i+1:j]+[route[i]]+route[j+1:], best_delta
    return None, None
