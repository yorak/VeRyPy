# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from sys import stderr

from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from logging import log, DEBUG
import numpy as np
import sys
import os

from cvrp_io import write_TSPLIB_file    

from config import LKH_EXE_PATH, LKH_EXACT_DISTANCES_PRECISION_DECIMALS

# Make sure we have access to LKH executable (download, comiple, and modify config.py).
if not os.path.isfile(LKH_EXE_PATH):
    raise ImportError("No LKH executable found at \"%s\" (the path is defined in config.py)"%LKH_EXE_PATH) 
if not os.access(LKH_EXE_PATH, os.X_OK):
    raise ImportError("LKH executable is not set executable") 

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

    # 1. Create a TSPLIB file
    are_float_distances = issubclass(D.dtype.type, np.float)
    if not are_float_distances:
        float_accuracy = None
    with NamedTemporaryFile( delete=False, suffix='.tsp') as tmpfile:
        temp_problem_file_path = tmpfile.name
    write_TSPLIB_file(temp_problem_file_path, D,
                      selected_idxs=selected_idxs, 
                      float_to_int_precision=float_accuracy)
    
    # 2. Create a parameter file for lkh.exe
    with NamedTemporaryFile( delete=False, suffix='.par') as tmpfile:
        temp_parameter_file_path = tmpfile.name
    with NamedTemporaryFile( delete=False, suffix='.tspsol') as tmpfile:
        temp_output_file_path = tmpfile.name
    with open(temp_parameter_file_path, 'w') as problem_file:
        problem_file.write("PROBLEM_FILE = %s\n" % temp_problem_file_path)
        problem_file.write("OUTPUT_TOUR_FILE = %s\n" % temp_output_file_path)
        problem_file.write("TRACE_LEVEL = 0\n")
        #problem_file.write("SEED = 42\n")
        if num_runs:
            problem_file.write("RUNS = %d\n"%num_runs)
            
    
    # 3. Call lkh
    command = [LKH_EXE_PATH, temp_parameter_file_path]
    p = Popen(command, stdout=PIPE, stdin=PIPE)
    
    # In Python3 bytes go in and come out. Special handling is required.
    #  (not pretty, but allows smoketests to pass on Py3)
    if sys.version_info[0] >= 3:
        stdout_data = p.communicate(input=b' ')[0].decode('ascii')
    else:
        stdout_data = p.communicate(input=' ')[0]

    
    if __debug__:    
        log(DEBUG-3, stdout_data)
    
    #TODO: CHECK FOR LKH FAILURE?
    #if "Press any key to continue" in stdout_data[-40:]:
    #   lkh_successfull = True
    
    # 4. Process output
    sol = []
    obj_f = 0
    tail = []
    with open(temp_output_file_path, 'r') as output_file:
        skip = True
        depot_found = False
        for l in output_file.readlines():
            l = l.strip()
            if l == "EOF" or l == "-1":
                skip = True
            elif "Length = " in l:
                obj_f = float(l.split()[-1])
            elif not skip:
                vrp_nid = selected_idxs[int(l)-1]
                if vrp_nid==0:
                    tail.append(0)
                    depot_found = True
                if depot_found:                    
                    sol.append(vrp_nid)
                else:
                    tail.append(vrp_nid)
            if l == "TOUR_SECTION":
                skip = False
        sol+=tail
        if not depot_found:
            sol+=[sol[0]]
        
    os.remove(temp_problem_file_path) 
    os.remove(temp_parameter_file_path) 
    os.remove(temp_output_file_path) 
    
    if are_float_distances:
        obj_f = obj_f/float_accuracy;
    
    return sol, obj_f
    
if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("lkh", solve_tsp_lkh)
