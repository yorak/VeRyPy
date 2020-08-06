#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of Holmes & Parker (1976) 
savings suppression heuristic. It runs the parallel savings heurstic multiple 
times suppressing some of the merges that have been made on the earlier
iterations.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: parallel savings heuristic 
procedure from parallel_savings.py, and numpy and scipy for reading and
preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from logging import log, DEBUG

from classic_heuristics.parallel_savings import parallel_savings_init
from util import objf, is_better_sol

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"



def supression_savings_function(D, suppressed, savings_cache=None):
    """
    This does the basic savings calculation:
        D[i,0]+D[0,j]-D[i,j]
    but with a twist: A set of *suppressed* merges is kept, and these
    (i,j) pairs are set to have savings value of 0.0 in the savings list. To
    avoid recalculation of the savings values, an empty savings_cache (list
    object) can be given. In subsequent iterations, it is then reused and
    the *topmost* suppressed values (according to the suppressed set) are set
    to savings value 0.0 and the savings list is then re-sorted."""
    
    N = len(D);
    n = N-1
    if not savings_cache:
        if type(savings_cache) is list:
            savings_cache[:] = [None]*int((n*n-n)/2)
        else:
            savings_cache = [None]*int((n*n-n)/2)
        idx = 0
        for i in range(1,N):
            for j in range(i+1,N):
                if (i,j) in suppressed:
                    savings_cache[idx] = (0,-D[i,j],i,j)
                else:
                    # Note,  Holmes & Parker (1976) do not use neg. values
                    s = max(0, D[i,0]+D[0,j]-D[i,j])
                    savings_cache[idx] = (s,-D[i,j],i,j)
                idx+=1
        savings_cache.sort(reverse=True)
    else:
        # reuse the savings from the savings cache
        modified_savings = False
        si = 0
        # remove all suppressed from the cache
        while (savings_cache[si][2], savings_cache[si][3]) in suppressed:
            # replace the savings record with one with suppressed value (0)
            savings_cache[si] = (0, savings_cache[si][1],
                         savings_cache[si][2], savings_cache[si][3],)
            modified_savings = True
            si+=1
        if modified_savings:
            savings_cache.sort(reverse=True)
            
    
    first_savings = savings_cache[0]
    suppress = (first_savings[2], first_savings[3])
    if __debug__:
        log(DEBUG, "Suppressing merge %s"%str(suppress))
    suppressed.add( suppress )
    
    
    return savings_cache 

def suppression_savings_init(D,d,C,L, minimize_K=False, Lprime="auto"):
    """
    This is the "vehicle scheduling procedure based upon savings and a solution
    pertubation scheme" of Holmes & Parker (1976). It works by suppressing a 
    certain number (at most Lprime) best available merges. Thus the parallel
    savings algorithm of Clarke & Wright (1964) is executed Lprime times.
    Please see the parallel_savings.py:parallel_savings_init for details and
    description of the parameters.
    
    * Lprime is the maximum suppression number L'. If set to "auto" (default) a
       linear prediction of a suitable value is used: L'=N/7+5, where N is the
       number of customers in the problem. Can only be int.
       
       
    Holmes, R. and Parker, R. (1976). A vehicle scheduling procedure based upon
    savings and a solution perturbation scheme. Journal of the Operational Re-
    search Society, 27(1):83–92.
    
    Clarke, G. and Wright, J. W. (1964). Scheduling of vehicles from a central
    depot to a number of delivery points. Operations research, 12(4):568-581.
    """
    N = len(D)-1
    if Lprime=="auto":
        # according to the (limited) experimental data of the Holmes & Parker
        #  (1976), for best results this can grow linearly.
        Lprime = int(N/7)+5
    Lprime = min(Lprime, int((N**2-N)/2))
    
    best_sol = None
    best_f = None
    best_K = None
    interrupted = False
    #best_sL = 1
    
    currently_suppressed_merges = set()
    savings_cache = []
    suppressed_f = lambda D:supression_savings_function(D,
                                currently_suppressed_merges, savings_cache)
    for Lcounter in range(Lprime):

        sol, sol_f, sol_K = None, float('inf'), float('inf')
        try:
            # On later invocations of parallel_savings_init, suppressed_f 
            #  generates an different set of savings values (suppressing some 
            #  the first/best of the previous iteration).
            sol = parallel_savings_init(D,d,C,L, minimize_K, suppressed_f)
            sol_f = objf(sol, D)
            sol_K = sol.count(0)-1
        except KeyboardInterrupt as e: # or SIGINT
            if len(e.args)>0 and type(e.args[0]) is list:
                sol = e.args[0]
                interrupted = True
        if sol:
            sol_f = objf(sol, D)
            sol_K = sol.count(0)-1
        if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
            if __debug__:
                log(DEBUG, "Best so far solution %s (%.2f) found with L'=%d"%
                           (sol,sol_f,Lcounter))
            best_sol = sol
            best_f = sol_f
            best_K = sol_K
            #best_sL = iterL
            
        if interrupted:
            raise KeyboardInterrupt(best_sol)
            
    return best_sol
          
# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ims_algorithm():
    algo_name = "HP76-PS|IMS"
    algo_desc = "Holmes & Parker (1976) parallel savings supression algorithm"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return suppression_savings_init(D,d,C,L,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ims_algorithm())
