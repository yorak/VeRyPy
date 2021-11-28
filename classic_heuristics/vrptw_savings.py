#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementations of the VRPTW savings heuristic.
Note, that in contrast to other VeRyPy heuristics, this supports asymmetric
cases.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: parallel savings heuristic 
procedure from parallel_savings.py, built-in TSP solver, and numpy and scipy
for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from builtins import range

import numpy as np
from util import TW, objf, is_better_sol
from classic_heuristics.parallel_savings import parallel_savings_init

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2021, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def placeholder_vrptw_savings_f(D, ctrs):
    n = len(D)-1
    savings = []
    # Note, VRPTW is asymmetric, we need to consider full D
    for i in range(1,n+1):
        for j in range(1,n+1):
            raise NotImplemented("Please implement me to enable TW support.")
            if i==j: continue

            s_D = D[i,0]+D[0,j]-D[i,j]
            s_t = 0    # TODO: Implement some way of describing "closeness" in 
                       #  time between two points (i and j) for *s_t* using
                       #   ctrs['TWs'][i][TW.OPEN]
                       #   ctrs['TWs'][i][TW.CLOSE]
                       #   ctrs['TWs'][j][TW.OPEN]
                       #   ctrs['TWs'][j][TW.CLOSE]
                       #   ctrs['TWs'][j][TW.CAN_WAIT]
                       #   and possibly D[i,j]
                       # Note that it should be comparable to s_D so that neither
                       #  dominates!
                       #
                       # The -D[i,j] is a tiebraker in case two potential merges
                       #  have a same savings value. It also can consider TWs
                       #  if one so wishes.
            savings.append( (s_D+s_t,-D[i,j],i,j) )
    savings.sort(reverse=True)

    return savings 

def vrptw_savings_init(D,d,ctrs,minimize_K=False):
    """ Savings algorithm with VRTPW savings criteria (see, e.g. Solomon 1987 or
    Van Landeghem 1988 for inspiration).

    It uses parallel_savings.py for the savings procedure. 
    Uses a revised constraints dictionary (see below) that can contain
    timewindow information.
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is a dict containing the constraints. Usually some combination of:
       - "C" is the capacity constraint limit for the identical vehicles.
       - "L" is the optional constraint for the maximum route length/duration.
       - "TWs" a list of timewindow pairs for each customer (0 being the depot)
        
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring negative savings values.
    
    
    Solomon, M. M. (1987). Algorithms for the vehicle routing and scheduling
       problems with time window constraints. Operations research, 35(2).

    Van Landeghem, H. R. G. (1988). A bi-criteria heuristic for the vehicle
       routing problem with time windows. European Journal of Operational 
       Research, 36(2), 217-226.
    """

    try:
        sol = parallel_savings_init(D,d,ctrs,minimize_K,placeholder_vrptw_savings_f)
    except KeyboardInterrupt as e:  #or SIGINT
        # lambda or pi was interrupted
        if len(e.args)>0 and type(e.args[0]) is list:
            sol = e.args[0]

    return sol

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_gs_algorithm():
    algo_name = r"So87-PS|VRPTW"
    algo_desc = r"Parallel savings algorithm for asymmetric distances "+\
                r"with Solomon (1987) VRPTW criteria"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        # Convert legacy call to new VRPTW supported init
        ctrs = {}
        if C: ctrs['C']=C
        if L: ctrs['L']=L
        return vrptw_savings_init(D,d,ctrs,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_gs_algorithm())
