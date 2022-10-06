# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from sys import stderr

from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from logging import log, DEBUG
from turtle import end_fill
import numpy as np
import sys
import os

from verypy.config import LKH_EXACT_DISTANCES_PRECISION_DECIMALS

import elkai


def solve_tsp_lkh(D, selected_idxs,
                  float_accuracy = LKH_EXACT_DISTANCES_PRECISION_DECIMALS,
                  num_runs = None):
    """
    A wrapper for LKH (Helsgaun 2006, 2009) TSP solver. Prepares the necessary
    problem and parameter files, calls the LKH executable, and interprets the
    results. 
    
    Helsgaun, K. (2006), "An Effective Implementation of K-opt Moves for the 
      Lin-Kernighan TSP Heuristic." DATALOGISKE SKRIFTER, No. 109, 2006. 
      Roskilde University.
    Helsgaun, K. (2009), "General k-opt submoves for the Lin-Kernighan TSP 
      heuristic." Mathematical Programming Computation, 1(2), 119-163.
      
    NOTE: Only symmetric problems are supported. If real valued D is given
    the accuracy of the optimization can be adjusted using the float_accuracy
    argument (the default accuracy is set in config.py"""
    
    if len(selected_idxs)<=3:
        p=None
        sol = list(sorted(selected_idxs))
        sol.append(sol[0])
        sol_f = 0
        for n in sol:
            if p is not None:
                sol_f+=int(D[p,n])
            p=n
        return sol, sol_f

    M = D[selected_idxs, :][:, selected_idxs]

    are_float_distances = issubclass(M.dtype.type, np.float)
    if are_float_distances:
        M = int(M * float_accuracy)
    else:
        M = M
    
    if num_runs is None:
        num_runs = 1
        
    output = elkai.solve_int_matrix(M, runs=num_runs)





    # 4. Process output
    sol = []
    obj_f = 0
    tail = []
    for l in output:
        depot_found = False

        vrp_nid = selected_idxs[l]
        if vrp_nid==0:
            tail.append(0)
            depot_found = True
        if depot_found:                    
            sol.append(vrp_nid)
        else:
            tail.append(vrp_nid)

    sol += tail
    if not depot_found:
        sol+=[sol[0]]


    # 5. Calculating the tour length
    obj_f = 0 
    for i in range(len(sol)-1):
        obj_f += D[sol[i], sol[i+1]]


    #TODO: CHECK FOR LKH FAILURE?
    #if "Press any key to continue" in stdout_data[-40:]:
    #   lkh_successfull = True




    return sol, obj_f

    
if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("lkh", solve_tsp_lkh)
