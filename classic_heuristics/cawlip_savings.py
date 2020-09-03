#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides implementation of the Robbins and Turner (1979)
savings heuristic with Lin interchange (2-opt) post-optimization phase.

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has moderate dependencies: parallel savings heuristic 
procedure from parallel_savings.py, built-in 2-opt and 2-opt* heuristics,
and numpy and scipy for reading and preparing the problem instance."""
###############################################################################

from classic_heuristics.parallel_savings import parallel_savings_init
from local_search import LSOPT, do_local_search
from local_search.intra_route_operators import do_2opt_move 
#from local_search.inter_route_operators import do_2optstar_move 

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


def cawlip_savings_init(D,d,C,L,minimize_K=False):
    """
    This implements the Robbins & Turner (1979) extension to the parallel
     savings algorithm of Clarke and Wright (1964). Clarke and Wright algorithm
     and its saving function are used as is, but the savings procedure is
     followed by 2-opt* improvement phase, where all possible ways of
     connecting any two edges of the solution are searched and improving moves
     are accepted. Please see the parallel_savings.py:parallel_savings_init for
     details and description of the parameters.
    """
    
    # Note: this is not a proper closure. Variables g and f are shared
    #  over all iterations. It is OK like this, but do not use/store the 
    #  lambda after this loop.
    

    sol = parallel_savings_init(D,d,C,L,minimize_K)
    
    sol = do_local_search([do_2opt_move],#, do_2optstar_move],
                          sol, D, d, C, L, LSOPT.BEST_ACCEPT)
    
    return sol
            

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_ps2o_algorithm():
    algo_name = "RT79-CAWLIP"
    algo_desc = "Robbins and Turner (1979) CAWLIP parallel "+\
                "savings algorithm with 2-opt* improvement phase"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return cawlip_savings_init(D,d,C,L,minimize_K)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_ps2o_algorithm())
