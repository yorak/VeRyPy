#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides a brute force solver that tries all permutations
for tiny problems.

The script is callable and can be used as a standalone sequental insertion
solver for TSPLIB formatted CVRPs. It has moderate dependencies as it requires
3rs party module dllist, in addition to numpy and scipy for reading and
preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from math import ceil
from itertools import permutations, chain
from sys import float_info
from logging import log, DEBUG, WARNING

from verypy.util import objf, routes2sol
from verypy.cvrp_ops import normalize_solution, recalculate_objective
from verypy.config import COST_EPSILON as S_EPS
from verypy.config import CAPACITY_EPSILON as C_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2023, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@gmail.com"
__status__ = "Development"

MAXF = float_info.max

def _validate_solution(solution, f_cutoff, D, d, C, L=None):
    cum_c = 0
    cum_l = 0
    f = 0
    prev_n = None
    for n in solution:
        if f>f_cutoff:
            return None, MAXF
        if C:
            cum_c += d[n]
            if cum_c-C_EPS>C:
                return False, MAXF
        if prev_n is not None:
            f += D[prev_n, n]
            if L:
                cum_l += D[prev_n, n]
                if cum_l-S_EPS>L:
                    return False, MAXF

        if n == 0:
            cum_c = 0
            cum_l = 0

        prev_n = n

    return True, f

def brute_force_init(D, d, C, L=None, minimize_K=False):
    """ A brute force heuristic. Just try different permutations. """

    if (len(D)>12):
        print("ERROR: Find better use for your CPU time.")
        return None
    elif (len(D)>11):
        print("WARNING: This will take ridiculously long.")
    elif (len(D)>10):
        print("WARNING: This will take quite a while.")
    elif (len(D)>9):
        print("WARNING: This will take a while.")

    best_f = float_info.max
    best_sol = None

    try:
        min_k = 1
        if C:
            min_k = ceil(sum(d)/float(C)-C_EPS)
        for k in range(min_k,len(D)+1):
            customers_and_depot_visits = list(range(1,len(D)))+[0]*(k-1)
            if __debug__:
                log(DEBUG-2, ("Checking permutations for K=%d"%k))
            for ordering in permutations(customers_and_depot_visits):
                solution = [0]+list(ordering)+[0]
                feasible, f = _validate_solution(solution, best_f, D, d, C, L)
        
                if __debug__ and feasible:
                    log(DEBUG-2, "Found a permutation %s (%.2f vs %.2f) which was %s"%(
                        solution, recalculate_objective(solution, D), 
                        -1.0 if f is None else f,
                        "feasible" if feasible else "infeasible"))

                if feasible and f+S_EPS<best_f:
                    best_sol = solution
                    best_f = f

            if minimize_K and best_sol:
                break # look no larger Ks

    except KeyboardInterrupt:  #or SIGINT
        if not best_sol:
            pieces = [[0, n] for n in range(1,len(D))]
            best_sol = chain( *pieces )+[0]
        raise KeyboardInterrupt(best_sol)
                              
    return best_sol
                


# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_bf_algorithm():
    algo_name = "BF"
    algo_desc = "Brute force algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return brute_force_init(D, d, C, L=L, minimize_K=minimize_K)
    call_init.__doc__ = brute_force_init.__doc__
    return (algo_name, algo_desc, call_init)

if __name__=="__main__":
    from verypy.shared_cli import cli
    cli(*get_bf_algorithm())
        
