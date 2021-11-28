# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides a set of shared Command Line Interface (CLI)
fuctionality for the classical heuristics."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import sys
from time import time
from os import path
from glob import glob
import logging

from natsort import natsorted

import cvrp_ops
import cvrp_io
from util import objf, sol2routes, is_better_sol
from config import DEBUG_VERBOSITY as DEFAULT_DEBUG_VERBOSITY

def print_problem_information(points, D, d, C, L, service_time, tightness=None, verbosity=0):
    N=len(D)
    print("SIZE:", N)
    if C:
        if verbosity>0 and tightness:
            print("TIGHTNESS:", "%.3f"%tightness)
        print("CAPACITY:", C)
    else:
        tightness = 0
    print("DISTANCE:", L)
    print("SERVICE_TIME:", service_time)
    
    if verbosity>2:
        print("POINTS:", points)
        print("DEMANDS:", d, "\n")
    if verbosity>3:
        print("D:", D)
        
    
def print_solution_statistics(sol, D, D_cost, d, C, L=None, service_time=None,
                              verbosity=-1):
    print("\nSOLUTION:", sol)
    cover_ok,capa_ok,rlen_ok = cvrp_ops.check_solution_feasibility(
                                          sol, D_cost,d,C,L,True)
    
    if verbosity>1:
        print("ALL SERVED:", cover_ok)   
        if C:
            print("IS C FEASIBLE:", capa_ok)
        if L:
            print("IS L FEASIBLE:", rlen_ok)
    else:
        print("FEASIBLE:", cover_ok and capa_ok and rlen_ok)
    print("SOLUTION K:", sol.count(0)-1)
    
    sol_f = None if D is None else objf(sol, D) 
    sol_c = None if D_cost is None else objf(sol, D_cost) 
    if (verbosity>0 and sol_f!=sol_c) or (not sol_c):
        print("SOLUTION COST:",sol_c, "\n")
    if sol_c:
        print("SOLUTION LENGTH:",sol_f)
    
    if verbosity>1:
        routes = sol2routes(sol)
        print("ROUTES:")
        print("No.\tCost\tLength\tLoad\tRoute")
        for i, route in enumerate(routes):
            print(i+1,
                  "%.2f"%objf(route,D_cost), 
                  "%.2f"%objf(route,D),
                  sum( (d[n] for n in route )) if C else "-",
                  route, sep='\t' )
        print("Total",
              "%.2f"%objf(sol,D_cost),
              "%.2f"%objf(sol,D), sep='\t')

def read_and_solve_a_problem(problem_instance_path, with_algorithm_function,
                             minimize_K, best_of_n=1, verbosity=-1,
                             single=False, measure_time=False):
    """ Solve a problem instance with the path in problem_instance_path
    with the agorithm in <with_algorithm_function>.
    
    The <with_algorithm_function> has a signature of:
    init_f(points, D_c, d, C, L, st, wtt, verbosity, single, minimize_K)
    
    Options <verbosity>, <single> and <measure_time> may be used to adjust what
    is printed and if a restricted single iteration search (different meaning 
    for different algorithms) is made."""
    
    pfn = problem_instance_path
    N, points, dd_points, d, D, C, ewt = cvrp_io.read_TSPLIB_CVRP(pfn)
    required_K, L, st = cvrp_io.read_TSBLIB_additional_constraints(pfn)
    
    # model service time with the distance matrix
    D_c = cvrp_ops.D2D_c(D, st) if st else D
        
    if points is None:
        if dd_points is not None:
            points = dd_points
        else:
            points, ewt = cvrp_ops.generate_missing_coordinates(D)

    tightness = None
    if C and required_K:
        tightness = (sum(d)/(C*required_K))
    if verbosity>=0:
        print_problem_information(points, D, d,C,L,st,tightness,verbosity)
    
    best_sol = None
    best_f = float('inf')
    best_K = len(D)
    interrupted = False
    for repeat_n in range(best_of_n):
        
        sol, sol_f, sol_K = None, float('inf'), float('inf')
        start = time()
        try:
            sol = with_algorithm_function(points, D_c, d, C, L, st,
                                          ewt, single, minimize_K)
        except KeyboardInterrupt as e:
            print ("WARNING: Solving was interrupted, returning "+
                   "intermediate solution", file=sys.stderr)
            interrupted = True
            # if interrupted on initial sol gen, return the best of those
            if len(e.args)>0 and type(e.args[0]) is list:
                sol = e.args[0]
        elapsed = time()-start 
        
        if sol:      
            sol = cvrp_ops.normalize_solution(sol)
            sol_f = objf(sol, D_c)
            sol_K = sol.count(0)-1
            if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
                best_sol = sol
                best_f = sol_f
                best_K = sol_K
            if best_of_n>1 and verbosity>=1:
                print("SOLUTION QUALITY %d of %d: %.2f"%
                      (repeat_n+1,best_of_n, objf(best_sol, D_c)))
            if measure_time or verbosity>=1:
                print("SOLVED IN: %.2f s"%elapsed)
                
        if interrupted:
            break
            
    if verbosity>=0 and best_sol:
        n_best_sol = cvrp_ops.normalize_solution(best_sol)
        print_solution_statistics(n_best_sol, D, D_c, d, C, L, st, verbosity=verbosity)
    
    if interrupted:
        raise KeyboardInterrupt()
        
    return best_sol, objf(best_sol, D), objf(best_sol, D_c)

def get_a_problem_file_list(problem_paths):
    files_to_solve = []
    for problem_path in problem_paths:
        if path.isdir(problem_path):
            for in_fn in natsorted(glob(path.join(problem_path, "*.vrp"))):
                files_to_solve.append( in_fn )
            for in_fn in natsorted(glob(path.join(problem_path, "*.tsp"))):
                files_to_solve.append( in_fn )
            for in_fn in natsorted(glob(path.join(problem_path, "*.pickle"))):
                files_to_solve.append( in_fn )
        elif path.isfile(problem_path) and problem_path[-4:].lower()==".txt":
            with open(problem_path, 'r') as vrp_list_file:
                for line in vrp_list_file.readlines():
                    line = line.strip()
                    if path.isfile(line):
                        files_to_solve.append(line)
        elif path.isfile(problem_path) and (problem_path[-4:].lower()==".vrp" or
                                            problem_path[-4:].lower()==".tsp" or
                                            problem_path[-7:].lower()==".pickle"):
            files_to_solve.append( problem_path )
        else:
            print(problem_path, "is not a .vrp file, folder, or text file",
                  file=sys.stderr)
    return files_to_solve

def set_logger_level(level, logfile=None):
    #set the logger verbosity level
    if level>=0:            
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.DEBUG-level,
                            stream=sys.stdout)
        for lvl in range(1,10):
            logging.addLevelName(lvl, "DEBUG")

        if logfile is not None:
            fileloghandler = logging.FileHandler(logfile)
            fileloghandler.setLevel(logging.DEBUG-level)
            fileloghandler.setFormatter( logging.Formatter("%(levelname)s:%(message)s") )
            logging.getLogger('').addHandler(fileloghandler)

def enable_function_tracing():
    """ Overkill for any other situation but a very deep deep debugging.

    Shows a complete list of all functions that get called during the invocation
    of the program!

    Code is from https://stackoverflow.com/a/8315566/1788710 . One can refer to
    that for more details on this Python magic!
    """
    def tracefunc(frame, event, arg, indent=[0], blacklist=["<lambda>", "<genexpr>"]):
       con = frame.f_code.co_name
       if con in blacklist:
           return tracefunc
       if event == "call":
           indent[0] += 2
           print("-" * indent[0] + "> call function", con)
       elif event == "return":
           print("<" + "-" * indent[0], "exit function", con)
           indent[0] -= 2
       return tracefunc
    sys.setprofile(tracefunc)

def tsp_cli(tsp_f_name, tsp_f):
    # import here so that the function can be used without these dependencies
    from util import objf
    
    if len(sys.argv)==2 and path.isfile(sys.argv[1]):
        P = cvrp_io.read_TSPLIB_CVRP(sys.argv[1])
        D = P.distance_matrix
        start_t = time()
        tsp_sol, tsp_f = tsp_f(D, list(range(len(D))))
        elapsed_t = time()-start_t
        print("Solved %s with %s in %.2f s"%(path.basename(sys.argv[1]), 
                                             tsp_f_name, elapsed_t))
        tsp_o = objf(tsp_sol,D)
        print("SOLUTION:", str(tsp_sol))
        print("COST:", tsp_o)  
        assert(tsp_f==tsp_o)
    else:
        print("usage: tsp_solver_%s.py TSPLIB_file.tsp"%tsp_f_name, file=sys.stderr)
         
def cli(init_name, init_desc, init_f):
    ## Simple command line interface
    single = False # ask to run only single iteration of the algorithm
    measure_time = False
    verbosity = DEFAULT_DEBUG_VERBOSITY
    minimize_K = False
    output_logfilepath = None
    best_of_n = 1
    interrupted = False
    
    for i in range(0, len(sys.argv)-1):
        if sys.argv[i]=="-v" and sys.argv[i+1].isdigit():
            verbosity = int(sys.argv[i+1])
        if sys.argv[i]=="-n" and sys.argv[i+1].isdigit():
            best_of_n = int(sys.argv[i+1])
        if sys.argv[i]=="-1":
            single = True
        if sys.argv[i]=="-t":
            measure_time = True
        if sys.argv[i]=="-l":
            output_logfilepath = sys.argv[i+1]       
        if sys.argv[i]=="-b":
            otarget = sys.argv[i+1].lower()
            if otarget=="cost" or otarget=="c":
                minimize_K = False
            elif otarget=="vehicles" or otarget=="k":
                minimize_K = True
            else:
                print("WARNING: Ignoring unknown optimization target %s"%otarget)
                
    if verbosity>=0:  
        set_logger_level(verbosity, logfile=output_logfilepath)
        
    if sys.argv[-1].isdigit():        
        N = int(sys.argv[-1])
        problem_name = "random "+str(N)+" point problem"
        N, points, _, d, D, C,_ = cvrp_io.generate_CVRP(N, 100, 20, 5)
        d = [int(de) for de in d]
        D_c = D
        L,st = None, None
        wtt = "EXACT_2D"
        
        best_sol = None
        best_f = float('inf')
        best_K = len(D)
        for i in range(best_of_n):
            sol, sol_f, sol_K = None, float('inf'), float('inf')
            try:
                sol = init_f(points, D_c, d, C, L, st, wtt, single, minimize_K)
            except KeyboardInterrupt as e:
                print ("WARNING: Solving was interrupted, returning "+
                       "intermediate solution", file=sys.stderr)
                interrupted = True
                # if interrupted on initial sol gen, return the best of those
                if len(e.args)>0 and type(e.args[0]) is list:
                    sol = e.args[0]
            if sol:      
                sol = cvrp_ops.normalize_solution(sol)
                sol_f = objf(sol, D_c)
                sol_K = sol.count(0)-1
                
                if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
                    best_sol = sol
                    best_f = sol_f
                    best_K = sol_K
                    
            if interrupted:
                break
                
        print_solution_statistics(best_sol, D, D_c, d, C, L, st, verbosity=verbosity)
    
    problem_file_list = get_a_problem_file_list([sys.argv[-1]])
    if not problem_file_list or "-h" in sys.argv or "--help" in sys.argv:
        print ("Please give a TSPLIB file to solve with "+\
          init_name+\
          " OR give N (integer) to generate a random problem of N customers."+\
          " OR give a path to a folder with .vrp files."+\
          "\n\nOptions (before the file name):\n"+\
          "  -v <int> to set the verbosity level (default %d)\n"%DEFAULT_DEBUG_VERBOSITY+\
          "  -n <int> run the algorithm this many times and return only the best solution\n"+\
          "  -1 to run only one iteration (if applicable)\n"+\
          "  -t to print elapsed wall time\n"+\
          "  -l <file_path> to store the debug output to a file\n"+\
          "  -b <'cost'|'vehicles'> or <c|K> sets the primary optimization oBjective (default is cost)",
          file=sys.stderr)
    elif problem_file_list:
        for problem_path in problem_file_list:
            problem_name = path.basename(problem_path)
            print("Solve", problem_name ,"with", init_name)
            read_and_solve_a_problem(problem_path, init_f, minimize_K, best_of_n, 
                                     verbosity, single, measure_time)
