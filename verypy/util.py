# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides shared utility functions for the algoritm
implementations."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from itertools import groupby    
from operator import itemgetter

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

  
def is_sorted(l):
    """Checks if the list is sorted """
    return all(l[i] <= l[i+1] for i in range(len(l)-1))

def produce_nn_list(D):
    """Produces a list of lists, each list has 2-tupes of node indices and 
    distances and from a node to all other nodes sorted by that distance. """
    
    n = len(D)
    # preprocess D to sorted_per_line_D
    NN_D = [None]*n
    for i in range(n):
        # sort each row
        NN_D[i] = sorted(enumerate(D[i,:]), key=itemgetter(1))
    return NN_D
    
def objf(sol, D):
    """A quick procedure for calclulating the quality of an solution (or a 
    route). Assumes that the solution (or the route) contains all visits (incl. 
    the first and the last) to the depot."""
    return sum(( D[sol[i-1],sol[i]] for i in range(1,len(sol))))

def totald(sol, d):
    """A quick procedure for calclulating the total demand of a solution
    (or a route)."""
    if not d: return 0
    return sum( d[n] for n in sol )

def first_valid(l):
    """Returns the first non-None or otherwise valid (evals to true) item from
    the list, or None if no such item exists."""
    return next((item for item in l if item), default=None)
    
def is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
    """Compares a solution against the current best and returns True if the
    solution is actually better accordint to minimize_K, which sets the primary
    optimization target (True=number of vehicles, False=total cost)."""
    
    if sol_f is None or sol_K is None:
        return False
    if best_f is None or best_K is None:
        return True
    elif minimize_K:
        return (sol_K<best_K) or (sol_K==best_K and sol_f<best_f)
    else:
        return sol_f<best_f

def without_empty_routes(sol):
    """Removes empty routes from the solution. WARNING: this also removes
    other concecutive duplicate nodes, not just 0,0!"""
    return [n[0] for n in groupby(sol)]
    
def sol2routes(sol):
    """Convert  solution to a list of routes (each a list of customers leaving 
    and returning to a depot (node 0). Removes empty routes. WARNING: this also 
    removes other concecutive duplicate nodes, not just 0,0!"""
    if not sol or len(sol)<=2: return []
    return [[0]+list(r)+[0] for x, r in groupby(sol, lambda z: z == 0) if not x]
     
def sol2edgeset(sol, symmetric=True):
    """Converts solution to a set of edges (2-tuples). If the problem is
    symmetric the tuples are directed from smaller to larger node value to 
    avoid duplicates."""
    edges = set()
    for i in range(0,len(sol)-1):
        j = i+1
        if symmetric or sol[i]<sol[j]:
            edges.add( (sol[i], sol[j]) )
        else:
            edges.add( (sol[j], sol[i]) )
    return edges

def routes2sol(routes):
    """Concatenates a list of routes to a solution. Routes may or may not have
    visits to the depot (node 0), but the procedure will make sure that 
    the solution leaves from the depot, returns to the depot, and that the 
    routes are separated by a visit to the depot."""
    if not routes:
        return None
    
    sol = [0]
    for r in routes:
        if r:
            if r[0]==0:
                sol += r[1:]
            else:
                sol += r
            if sol[-1]!=0:
                sol += [0]
    return sol


