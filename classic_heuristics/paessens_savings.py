#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementation of the Paessens (1988)
generalized savings heuristic. Also the automatic parameter search strategies
M1 and M4 are implemented.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: parallel savings heuristic 
procedure from parallel_savings.py, built-in 3-opt heuristic, and numpy and
scipy for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from builtins import range

from classic_heuristics.parallel_savings import parallel_savings_init
import numpy as np
from local_search import LSOPT, do_local_search
from local_search.intra_route_operators import do_3opt_move 
from util import objf, without_empty_routes, is_better_sol

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

M4_FINETUNE_STEP = 0.1

def _cartesian_product(nparray1, nparray2):
    return np.array(np.meshgrid(nparray1,nparray2)).T.reshape(-1,2)

def paessens_savings_function(D, g_multiplier, f_multiplier):
    """ the Paessens 1988 Savings function
    0.0<g<=3.0
    0.0<=f<=1.0
    """
    N = len(D);
    n = N-1
    savings = [None]*int((n*n-n)/2)
    idx = 0
    for i in range(1,N):
        for j in range(i+1,N):
            s = D[i,0]+D[0,j]\
                -g_multiplier*D[i,j]\
                +f_multiplier*abs(D[i,0]-D[0,j])
            savings[idx] = (s,-D[i,j],i,j)
            idx+=1
    savings.sort(reverse=True)
    return savings 

def paessens_savings_init(D,d,C,L, minimize_K=False,
                          strategy="M4", do_3opt=True):
    """
    This implements the Paesses (1988) variant of the parallel savings
     algorithm of Clarke and Wright (1964). The savings function of
     (Paesses 1988) is parametrized with multipiers g and f:
         
         S_ij  = d_0i + d_0j - g * d_ij + f * | d_0i - d_0j |
    
    If two merges have the same savings value, the one where i and j are closer
     to another takes precedence. Otherwise the impelementation details can be 
     read from parallel_savings.py as it contains the actual code implementing 
     the parallel savings procedure. The variant specific parameters are:
     
    * strategy which can be:
        - "M1" for 143 runs of the savings algorithm with all combinations of
           g = np.linspace(0.8, 2.0, num=13)
           f = np.linspace(0.0, 1.0, num=11)            
        - "M4" for 8 runs (g,f) = (1.0,0.1), (1.0,0.5), (1.4,0.0), (1.4,0.5)
           with a parameter combinations +/- 0.1 around the best of these four. 
        - or a list of (g,f) value tuples.
    * do_3opt (default True) optimize the resulting routes to 3-optimality
    
    Note: Due to the use of modern computer, and low priority in computational
     efficiency of this implementation, not all of the tecninques specified in
     "reduction of computer requirements" (Paessens 1988) were employed.
    """
    
    parameters = []
    if strategy=="M1":
        parameters.extend( _cartesian_product(np.linspace(0.8, 2.0, num=13),
                               np.linspace(0.0, 1.0, num=11)) )
    elif strategy=="M4":
        parameters.extend( [(1.0, 0.1), (1.0, 0.5), (1.4, 0.0), (1.4, 0.5)] )
    else:
        parameters.extend(strategy)
    
    best_params = None
    best_sol = None
    best_f = None
    best_K = None
    interrupted = False
    
    params_idx = 0
    while params_idx<len(parameters):
        g,f = parameters[params_idx]
        
        # Note: this is not a proper closure. Variables g and f are shared
        #  over all iterations. It is OK like this, but do not use/store the 
        #  lambda after this loop.
        gf_savings = lambda D: paessens_savings_function(D, g, f)
        
        sol, sol_f, sol_K = None, float('inf'), float('inf')
        try:
            sol = parallel_savings_init(D,d,C,L,minimize_K, gf_savings)
            if do_3opt:
                sol = do_local_search([do_3opt_move], sol, D, d, C, L,
                                      LSOPT.BEST_ACCEPT)
            # 3-opt may make some of the routes empty
            sol = without_empty_routes(sol)
        except KeyboardInterrupt as e: # or SIGINT
            # some parameter combination was interrupted
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
            best_params = (g,f)
        
        if interrupted:
            raise KeyboardInterrupt(best_sol)
        
        params_idx+=1
        # after the best of 4 for the M4 is found, check 4 more around it
        if params_idx==4 and strategy=="M4":
            g_prime, f_prime = best_params
            parameters.extend( [(g_prime-M4_FINETUNE_STEP, f_prime),
                                (g_prime+M4_FINETUNE_STEP, f_prime),
                                (g_prime, f_prime-M4_FINETUNE_STEP ),
                                (g_prime, f_prime+M4_FINETUNE_STEP )] )
    return best_sol
            

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_gps_algorithm():
    algo_name = "Pa88-PS|G2P"
    algo_desc = "Paessens (1988) parametrized parallel savings algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return paessens_savings_init(D,d,C,L,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_gps_algorithm())
