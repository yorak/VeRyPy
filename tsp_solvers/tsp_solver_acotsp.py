# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from sys import stderr, version_info

from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
import re
import os

import numpy as np

from cvrp_io import write_TSPLIB_file
from config import ACOTSP_EXE_PATH, ACOTSP_EXACT_DISTANCES_PRECISION_DECIMALS

# Make sure we have access to LKH executable (download, comiple, and modify config.py).
if not os.path.isfile(ACOTSP_EXE_PATH):
    raise ImportError("No acotsp executable found at \"%s\" (the path is defined in config.py)"%ACOTSP_EXE_PATH) 
if not os.access(ACOTSP_EXE_PATH, os.X_OK):
    raise ImportError("acotsp executable is not set executable") 

def solve_tsp_acotsp(D, selected_idxs,
                  float_accuracy = ACOTSP_EXACT_DISTANCES_PRECISION_DECIMALS,
                  use_ants = False):
    """
    A wrapper for ACOTSP (Stuzle 2002, Dorigo & Stuzle 2004) solver. This
    wrapper writes a temporary TSBLIB file, calls modified ACOTSP executable,
    and interprets the results. 
    
    The version of the ACOTSP solver used is modified from the v.1.03 of the
    original solver that can disable ants and only do the deterministic local 
    search. It is available at:
        
    https://github.com/juherask/ACOTSP
    
    The executable compiled from the modified code includes additional command
    line arguments for  disabling the ants by issuing command line arguments
    "-m 1 -r 1 -s 1 -n -j 0". The first arguments set only one trial for the
    deterministic version and the "-n" or "--noants" argument turns off all ant
    systems (pheromone trails etc.). The "-j 0" or "--lsrnd 0" turns off
    randomization of the local search search order. These are set if the
    function is called with the use_ants = False (default). Other modification
    includes allowing the printing of the solution to stdout when the --quiet
    flag is given. 
    
    With use_ants = False this is a very fast TSP solver as it builds solutions
    using nearest neighbors construction and then does 2-opt followed by 3-opt
    improvement steps. Using ant systems (i.e. use_ants = True) on the other
    hand is slow as it employs the stochastic metaheuristic in addition to same
    ops as without ants.
    
    Stutzle, T. (2002). ACOTSP: A software package of various ant colony
    optimization algorithms applied to the symmetric TSP.
    URL http://www.aco-metaheuristic.org/aco-code
    
    Dorigo, M. and Stuzle, T. (2004). Ant Colony Optimization, MIT Press,
    Cambridge, MA, USA.
    
    NOTE: Only symmetric problems are supported. Also, ACOTSP does not support 
    float distances, thus if a real valued D is given it is converted using the
    float_accuracy (which is a mutiplier, e.g. 1000.0, the default accuracy is
    set in config.py)"""
    
    if len(selected_idxs)<=3:
        p=None
        sol = list(sorted(selected_idxs))
        sol.append(sol[0])
        sol_f = sum(D[sol[i-1]][sol[i]] for i in range(1,len(sol)))
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
    
    # 2. Call acotsp
    if use_ants:
        # with defaults
        command = [ACOTSP_EXE_PATH, "--quiet", "-i", temp_problem_file_path]
    else:
        # use only local 2-opt+3-opt local search
        nul_f = open(os.devnull, 'w')
        command = [ACOTSP_EXE_PATH,
                   "-r", "1", # one try is enough without ants (determnistic)
                   "-m", "1", # one (deteministic) ant is still needed 
                              #  to hold the nearest neighbour + 3opt solution 
                   "-s", "1", # build only 1 route (tour)
                   "-n",      # disable all ant systemes (NEW)
                   "-j", "0", # turn of random local search order (NEW)
                   "-l", "3", # do 3-opt
                   "--quiet", # in 1.03 disables writing the files
                   "-i", temp_problem_file_path]
    p = Popen(command, stdout=PIPE, stdin=PIPE, stderr=nul_f)
    # In Python3 bytes go in and come out. Special handling is required.
    #  (not pretty, but allows smoketests to pass on Py3)
    if version_info[0] >= 3:
        stdout_data = p.communicate(input=b' ')[0].decode('ascii')
    else:
        stdout_data = p.communicate(input=' ')[0]

  
    #if __debug__:    
        #print(stdout_data)
        
    # 3. Process output
    sol = None
    sol_f = None
    
    best_re = re.compile("try [0-9]+, Best ([0-9]+(\.?[0-9]*)?),")
    best_obj = None
    for l in stdout_data.split('\n'):
        if "best solution so far is" in l[:25]:
            tour_indices = [int(n) for n in l.replace("best solution so far is","").split()]
            sol = list( selected_idxs[i] for i in tour_indices )
            first_pos = sol.index(selected_idxs[0])
            sol = sol[first_pos:]+sol[:first_pos]+[selected_idxs[0]]
            sol_f = sum(D[sol[i-1]][sol[i]] for i in range(1,len(sol)))
            #print("REMOVEME:", sol_f, best_obj)
            
        bo = best_re.search(l)
        if bo:
            best_obj = float(bo.group(1))
            if are_float_distances:
                best_obj = best_obj/float_accuracy;  
  
    os.remove(temp_problem_file_path) 
    return sol, sol_f
    
if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("acotsp", solve_tsp_acotsp)
