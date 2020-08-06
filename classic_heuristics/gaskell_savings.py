#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementations of the Gaskell (1967)
\pi and \lambda savings functions for parallel (as in multiple route)
savings heuristic.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: parallel savings heuristic 
procedure from parallel_savings.py, built-in TSP solver, and numpy and scipy
for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from builtins import range

import numpy as np
from util import objf, is_better_sol
from classic_heuristics.parallel_savings import parallel_savings_init

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def gaskell_lambda_savings_function(D):
    n = len(D)-1
    savings = [None]*int((n*n-n)/2)
    idx = 0
    d_avg = np.average(D[0:])
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            s_AB = D[i,0]+D[0,j]-D[i,j]
            lambda_AB = s_AB*(d_avg+abs(D[0,i]-D[0,j])-D[i,j])
            savings[idx] = (lambda_AB, -D[i,j], i,j)
            idx+=1            
    savings.sort(reverse=True)

    return savings 

def gaskell_pi_savings_function(D):
    n = len(D)-1
    savings = [None]*int((n*n-n)/2)
    idx = 0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            pi_AB = D[i,0]+D[0,j]-2*D[i,j]
            savings[idx] = (pi_AB, -D[i,j], i,j)
            idx+=1            
    savings.sort(reverse=True)    
    return savings 
    
def gaskell_savings_init(D,d,C,L, minimize_K=False, savings_method="both"):
    """ Savings algorithm with Gaskell (1967) pi and lambda savings criteria.
    Uses parallel_savings.py for the savings procedure. 
    
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.
        
    * minimize_K sets the primary optimization objective. If set to True, it is
       the minimum number of routes. If set to False (default) the algorithm 
       optimizes for the mimimum solution/routing cost. In savings algorithms 
       this is done by ignoring negative savings values.
       
    * savings_method selects the savings criteria: "lambda" or "pi". If set to
      "both" (default) the one with better results is returned.
    
    Gaskell, T. (1967). Bases for vehicle fleet scheduling. Journal of the
    Operational Research Society, 18(3):281-295.
    """

    savings_functions = []
    if savings_method=="both":
        savings_functions = [gaskell_lambda_savings_function, 
                             gaskell_pi_savings_function]
    elif savings_method=="lambda":
        savings_functions = [gaskell_lambda_savings_function]
    elif savings_method=="pi":
        savings_functions = [gaskell_pi_savings_function]
    else:
        raise ValueError("Only 'lambda', 'pi', or 'both' are supported.")
    
    best_sol = None
    best_f = None
    best_K = None
    interrupted = False
    for sav_f in savings_functions:
        sol, sol_f, sol_K = None, float('inf'), float('inf')
        try:
            sol = parallel_savings_init(D,d,C,L,minimize_K,sav_f)
        except KeyboardInterrupt as e:  #or SIGINT
            # lambda or pi was interrupted
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
            # pass on the current best_sol
            raise KeyboardInterrupt(best_sol)

    return best_sol

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_gs_algorithm():
    algo_name = r"Ga67-PS|pi+lamda"
    algo_desc = r"Parallel savings algorithm with Gaskell (1967) $\pi$ and "+\
                r"$\lambda$ criteria"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        savings_method = "both" if not single else "pi"
        return gaskell_savings_init(D,d,C,L, minimize_K, savings_method)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_gs_algorithm())
