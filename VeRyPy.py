#!/usr/bin/env python
################################################################################
# -*- coding: utf-8 -*-
""" Provides a main callable interface to the algorithms provided by the
VeRyPy classical vehicle routing problem heuristic library. Can solve TSPLIB 
formatted CVRP problems with capacity and/or maximum route cost/length/duration
constraints.
"""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import sys
import pickle
from argparse import ArgumentParser
from os import path
from time import time
from collections import defaultdict
import numpy as np

import cvrp_io
import cvrp_ops
import shared_cli

from local_search import do_local_search
from local_search.intra_route_operators import do_2opt_move, do_3opt_move
#TODO: offer other LS ops?

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

################################################################################

#TODO: Generic TODO list for the very high level 
# - Experiment with http://numba.pydata.org/ esp. in local search

PRINT_INSTANCE_DATA = True

algo_name_aliases = {
        # savings heuristics
        "ps":"ps",   "cw64-ps":"ps",
                     "cw":"ps",
                     "parallelsavings":"ps",
                     "clarkewright":"ps",
        "gpl":"gpl", "ga67-ps":"gpl","ga67-ps|pi":"gpl","ga67-ps|lambda":"gpl",
                     "gaskellpi":"gpl", "gaskelllambda":"gpl",
        "ss":"ss",   "we64-ss":"ss", "sequentialsavings":"ss",
        
        "gps":"gps", "pa88-ps":"gps","generalizedsavings":"gps",
                     "generalizedparallelsavings":"gps",
                     "paessens":"gps",
        "ims":"ims", "hp76-ps":"ims","hp76-ps|ims":"ims",
                     "suppression":"ims",
                     "iterativemergesuppression":"ims",
                     "mergesuppressionsavings":"ims",
        "ps2o":"ps2o", "rt79-ps":"ps2o", "cawlip":"ps2o",
                     "savingswith2opt":"ps2o", 
        
        # insertion heuristics
        "si":"si",   "ci":"si",
                     "cheapestinsertion":"si",
                     "sequentialinsertion":"si",
        "mj":"mj",   "mj76-si":"si",
                     "molejameson":"mj","molejamesoninsertion":"mj",
        "pi":"pi",   "parallelinsertion":"pi", "parallelcheapestinsertion":"pi",
        
        # maximum mathcing heuristics
        "mbsa":"mbsa",   "dv89-mbsa":"mbsa", "matching":"mbsa",
                     "maximummatching":"mbsa", "mm":"mbsa",
                     "desrochersverhoog":"mbsa",
        
        # the 2-phase heuristic
        "cmt":"cmt", "cmt79-2p":"cmt","tp":"cmt","twophase":"cmt",
                     "cmt2p":"cmt","cmt2phase":"cmt","cmttwophase":"cmt",
        
        "sn":"sn", "nn":"sn","snn":"sn", # parallel nearest neighbour
        "pn":"pn", "pnn":"pn", # parallel nearest neighbour
        "ty":"ty",   "ty68-snn":"ty","tyagi":"ty", # tyagi nearest neighbour
        
        "swp":"swp", "sweep":"swp", # plain sweep
        "wh":"wh",   "wh72-swls":"wh", "whs":"wh", "whswp":"wh",
                     "wrenholliday":"wh", "wrenhollidaysweep":"wh", 
        "gm":"gm",   "gm74-swri":"gm", "gms":"gm", "gmswp":"gm",
                     "gilletmiller":"gm", "gilletmillersweep":"gm",
        
        "rfcs":"rfcs", "be83-rfcs":"rfcs","routefirstclustersecond":"rfcs",
                       "beasley":"rfcs",
                       "newtonthomas":"rfcs",
        
        # generalized assignment problem heuristic
        "gap":"gap", "fj81-gap":"gap", "fisherjaikumar":"gap",
        
        # set covering "petal" heuristic
        "ptl":"ptl", "fr76-ptl":"ptl", "petal":"ptl", "fosterryan":"ptl",
        
        # lagrangian relaxation 3-opt* heuristic
        "lr3o":"lr3o", "sg82-lr3opt":"lr3o", "lr3opt":"lr3o",
        
        # TO ENABLE THEM ALL!
        "all":"all",
        "classical":"classical", "classic":"classical",
}

def get_algorithms(names):
    has_gurobi = True
    try:
        import gurobipy
    except:
        has_gurobi = False
    
    # get algorithms
    algos = []
    for algo_name in names:
        algo_name = algo_name.lower() # ignore case
        if algo_name in algo_name_aliases:
            # translate name to standard abbreviation
            algo_name = algo_name_aliases[algo_name]

            # import only those that are needed, when needed
            if algo_name in ["ps", "all", "classical"]:
                #"ps":"Clarke & Wright (1964) parallel savings heuristic",
                from classic_heuristics.parallel_savings import get_ps_algorithm
                algos.append( ("ps",)+get_ps_algorithm() )
            if algo_name in ["ps2o", "all"]:
                #"ps2o":"Robbins and Turner (1979) CAWLIP savings algorithm",
                from classic_heuristics.cawlip_savings import get_ps2o_algorithm
                algos.append( ("ps2o",)+get_ps2o_algorithm() )
            if algo_name in ["gpl", "all", "classical"]:
                #"gpl":"Gaskell (1967) best of pi/lambda savings heuristic",
                from classic_heuristics.gaskell_savings import get_gs_algorithm
                algos.append( ("gpl",)+get_gs_algorithm() )
            if algo_name in[ "ss", "all", "classical"]:
                #"ss":"Webb (1964) sequential savings algorithm"
                from classic_heuristics.sequential_savings import get_ss_algorithm
                algos.append( ("ss",)+get_ss_algorithm(lambda_multiplier=1.0) )
            if algo_name in ["gps", "all", "classical"]:
                #"gps":"Paessens (1988) parallel savings heuristic with generalized savings function and M4 strategy",
                from classic_heuristics.paessens_savings import get_gps_algorithm
                algos.append( ("gps",)+get_gps_algorithm() )
                pass
            if algo_name in ["ims", "all", "classical"]:
                #"ims":"Holmes & Parker (1976) iterative parallel savings merge suppression heuristic",
                from classic_heuristics.suppression_savings import get_ims_algorithm 
                algos.append( ("ims",)+get_ims_algorithm() )
            if algo_name in ["si", "all"]:
                #"si":"Mole & Jameson (1976) sequential insertion heuristic without local search"),
                from classic_heuristics.cheapest_insertion import get_si_algorithm 
                algos.append( ("si",)+get_si_algorithm() )        
            if algo_name in ["mj", "all", "classical"]:
                #"mj":"Mole & Jameson (1976) sequential insertion heuristic with local search"),                
                from classic_heuristics.mole_jameson_insertion import get_mj_algorithm 
                algos.append( ("mj",)+get_mj_algorithm() )
            if algo_name in ["pi", "all"]:
                #"pi":"van Breedam (1995) parallel insertion heuristic"),
                from classic_heuristics.cheapest_insertion import get_pi_algorithm
                algos.append( ("pi",)+get_pi_algorithm() )
            if algo_name in ["mbsa", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [mbsa/DV89-MM] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:
                    #"mbsa":"Desrochers and Verhoog (1989) matching based savings algorithm"
                    from classic_heuristics.matchingvrp import get_mm_algorithm
                    algos.append( ("mbsa",)+get_mm_algorithm())
            if algo_name in ["cmt", "all", "classical"]:
                #"cmt":"Christofides, Mingozzi & Toth (1979) two phase heuristic"
                from classic_heuristics.cmt_2phase import get_cmt2p_algorithm
                algos.append( ("cmt",)+get_cmt2p_algorithm() )
            if algo_name in ["sn", "all"]:
                #"sn":"Sequential Nearest Neighbor construction heuristic"
                from classic_heuristics.nearest_neighbor import get_snn_algorithm
                algos.append( ("sn",)+get_snn_algorithm() )
            if algo_name in ["pn", "all"]:
                #"pn":"Parallel Nearest Neighbor construction heuristic"
                from classic_heuristics.nearest_neighbor import get_pnn_algorithm
                algos.append( ("pn",)+get_pnn_algorithm() )
            if algo_name in ["ty", "all", "classical"]:
                #"ty":"Tyagi (1968) Nearest Neighbor construction heuristic"
                from classic_heuristics.tyagi_nearest_neighbor import get_ty_algorithm
                algos.append( ("ty",)+get_ty_algorithm() ) 
            if algo_name in ["swp", "all"]:
                #"swp":"Plain Sweep algorithm without route improvement heuristics"
               from classic_heuristics.sweep import get_swp_algorithm
               algos.append( ("swp",)+get_swp_algorithm() )
            if algo_name in ["wh", "all", "classical"]:
                #"wh":"Wren and Holliday (1972) Sweep + Local Search algorithm"
               from classic_heuristics.wren_holliday_sweep import get_wh_algorithm
               algos.append( ("wh",)+get_wh_algorithm() ) 
            if algo_name in ["gm", "all", "classical"]:
                #"gm":"Gillett and Miller (1974) Sweep algorithm"
                from classic_heuristics.gillet_miller_sweep import get_gm_algorithm
                algos.append( ("gm",)+get_gm_algorithm() )
            if algo_name in ["rfcs", "all", "classical"]:
                #"rfcs":"Route-first-cluster-second heuristic of Beasley (1983)"
                from classic_heuristics.rfcs import get_rfcs_algorithm
                algos.append( ("rfcs",)+get_rfcs_algorithm()) 
            if algo_name in ["gap", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [gap/FJ81-GAP] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:
                    #"gap":"Fisher & Jaikumar (1981) generalized assignment problem heuristic"
                    from classic_heuristics.gapvrp import get_gap_algorithm
                    algos.append( ("gap",)+get_gap_algorithm() )
            if algo_name in ["ptl", "all", "classical"]:
                if not has_gurobi:
                    print("WARNING: [ptl/FR76-1PTL] heuristic is not available (gurobipy is not installed).", file=sys.stderr)
                else:                    
                    #"ptl":"Foster and Ryan (1976) Petal algorithm"
                    from classic_heuristics.petalvrp import get_ptl_algorithm
                    algos.append( ("ptl",)+get_ptl_algorithm() )  
            if algo_name in ["lr3o", "all", "classical"]:
                #"lr3o":"Stewart and Golden (1982) LR3OPT heuristic"
                from classic_heuristics.lr3opt import get_lr3opt_algorithm
                algos.append( ("lr3o",)+get_lr3opt_algorithm() )  
        else:
            print(algo_name, "is not a valid algorithm name", file=sys.stderr)
    return algos

def _build_algorithm_help():
    algo_shorthands = sorted(list(set(algo_name_aliases.values())-set(["all","classical"])))
    algo_name_descriptors = get_algorithms(algo_shorthands)
    htxt = ""
    for key, algo_name, algo_desc, _ in algo_name_descriptors:
        htxt+=key+" : "+algo_name+" : "+algo_desc+"\n"
    return htxt

def main(overridden_args=None):
    ## 1. parse arguments
    
    parser = ArgumentParser(description="Solve some .vrp problems with the algorithms built into VeRyPy.")
    parser.add_argument('-l', dest='list_algorithms', help="List the available heuristics and quit", action="store_true")
    parser.add_argument('-v', dest='verbosity', help="Set the verbosity level (to completely disable debug output, run this script with 'python -O')", type=int, default=-1)
    parser.add_argument('-a', dest='active_algorithms', help="Algorithm to apply (argument can be set multiple times to enable multiple algorithms, or one can use 'all' or 'classical')", action='append')    
    parser.add_argument('-b', dest='objective', choices=['c', 'cost', 'K', 'vehicles'], help="Primary optimization oBjective (default is cost)", default="cost")
    parser.add_argument('-m', dest='minimal_output', help="Overrides the output options and prints only one line CSV report per solved instance", action="store_true")
    parser.add_argument('-t', dest='print_elapsed_time', help="Print elapsed wall time for each solution attempt", action="store_true")
    parser.add_argument('-c', dest='show_solution_cost', help="Display solution cost instead of solution length", action="store_true")
    parser.add_argument('-D', dest='dist_weight_format', choices=['ROUND', 'EXACT', 'TRUNCATE'], help="Force distance matrix rounding")
    parser.add_argument('-1', dest='use_single_iteration', help="Force the algorithms to use only single iteration.", action="store_true")    
    parser.add_argument('--iinfo', dest='print_instance_info', help="Print the instance info in the collected results", action="store_true")
    parser.add_argument('--routes', dest='print_route_stat', help="Print per route statistics of the final solution", action="store_true")
    parser.add_argument('--vrph', dest='print_vrph_sol', help="Print the final solution in the VRPH format", action="store_true")
    parser.add_argument('--forbid', dest='forbid_algorithms', help="Forbid applying algorithms (argument can set multiple times to forbid multiple algorithms)", action='append')    
    parser.add_argument('--simulate', dest='simulate', help="Do not really invoke algorithms, can be used e.g. to test scripts", action="store_true")
    #TODO: consider adding more LS opts e.g. 2optstart, 3optstart
    parser.add_argument('--post-optimize', dest='local_search_operators', choices=['2opt', '3opt'], help="Do post-optimization with local search operator(s) (can set multiple)", action='append')
    parser.add_argument("problem_file", help="a path of a .vrp problem file, a directory containing .vrp files, or a text file of paths to .vrp files", action='append')
    
    if overridden_args:
        app_args = parser.parse_args(overridden_args)
    elif "-l" in sys.argv:
        print("Select at least one algorithm (with -a) from the list:", file=sys.stderr)
        print(_build_algorithm_help())
        sys.exit(1)
    elif len(sys.argv)==1:
        print("Give at least one .vrp file and use -h to get help.", file=sys.stderr)
        sys.exit(1)
    else:
        app_args = parser.parse_args()
        
    # some further argument validation
    if not app_args.active_algorithms or app_args.list_algorithms:
        print("Select at least one algorithm (with -a) from the list:", file=sys.stderr)
        print(_build_algorithm_help())
        exit()
    if len(app_args.problem_file)==0:
        print("Provide at least one .vrp file to solve", file=sys.stderr)
        exit()
    
    # get .vrp file list
    files_to_solve = shared_cli.get_a_problem_file_list(app_args.problem_file)

    # get algorithms
    algos = get_algorithms(app_args.active_algorithms)
    if app_args.forbid_algorithms:
        forbidden_algos = [algo_name_aliases[algo_name]
                           for algo_name in app_args.forbid_algorithms
                           if (algo_name in app_args.forbid_algorithms)]
        algos = [a for a in algos if (a[0] not in forbidden_algos)]

    # get primary objective
    minimize_K = False
    if app_args.objective=='K' or app_args.objective=='vehicles':
        minimize_K = True
        
    run_single_iteration = False
    if app_args.use_single_iteration:
        run_single_iteration = True
        
    # get post-optimization local search move operators
    ls_ops = []
    ls_algo_names = []
    if app_args.local_search_operators:
        ls_algo_names = app_args.local_search_operators
        for ls_op_name in ls_algo_names:
            if ls_op_name=="2opt":
                ls_ops.append(do_2opt_move)
            if ls_op_name=="3opt":
                ls_ops.append(do_3opt_move)
    
    # verbosity
    if app_args.verbosity >= 0:
        shared_cli.set_logger_level(app_args.verbosity)

        # Magic number :) for very very very detailed debugging.
        if app_args.verbosity >= 10:
            shared_cli.enable_function_tracing()

    # minimal header
    if app_args.minimal_output:
        print("algo;problem;is_feasible;f;K;t")
        
    ## 2. solve 
    results = defaultdict(lambda: defaultdict(float))
    instance_data = dict()

    interrupted = False
    for pfn in files_to_solve:
        bn = path.basename(pfn).replace(".vrp","").replace(".tsp","").replace(".pickle","")
        
        try:
           N, points, dd_points, d, D, C, ewt, K, L, st = pickle.load( open( pfn, "rb" ) )
        except:
           N, points, dd_points, d, D, C, ewt = cvrp_io.read_TSPLIB_CVRP(pfn)
           K, L, st = cvrp_io.read_TSBLIB_additional_constraints(pfn)
        
        # We do not have point coodrinates, but we have D! 
        if points is None:
            if dd_points is not None:
                points = dd_points
            else:
                points, ewt = cvrp_ops.generate_missing_coordinates(D)
        
        if app_args.dist_weight_format == "TRUNCATE":
            D = np.floor(D)
            ewt = "FLOOR_2D"
        if app_args.dist_weight_format == "ROUND":
            D = np.int(D)
            ewt = "EUC_2D"
            
        # Bake service time to D (if needed)
        D_c = cvrp_ops.D2D_c(D, st) if st else D
        
        for algo_abbreviation, algo_name, _, algo_f in algos:
            if not app_args.minimal_output:
                print("Solving %s with %s"%(bn, algo_name))
            start_t = time()
            sol = None
            try:
                if not app_args.simulate:
                    sol = algo_f(points, D_c, d, C, L, st, ewt,
                                 run_single_iteration, minimize_K)
            except (KeyboardInterrupt, Exception) as e:
                if type(e) is KeyboardInterrupt:
                    interrupted = True
                    # if interrupted on initial sol gen, return the best of those
                    if len(e.args)>0 and type(e.args[0]) is list:
                        sol = e.args[0]
                    if not app_args.minimal_output:
                        print("WARNING: Interrupted solving %s with %s"%
                              (bn, algo_abbreviation), file=sys.stderr)
                else:
                    if not app_args.minimal_output:
                        print("ERROR: Failed to solve %s with %s because %s"%
                              (bn, algo_abbreviation, str(e) ), file=sys.stderr)
                    sol = None
            
            if sol:
                sol = cvrp_ops.normalize_solution(sol)
                if app_args.show_solution_cost:
                    sol_q = cvrp_ops.calculate_objective(sol, D_c)
                else:
                    sol_q = cvrp_ops.calculate_objective(sol, D)
                sol_K = sol.count(0)-1
                
                if app_args.local_search_operators:
                    if not app_args.minimal_output:
                        print("Postoptimize with %s ..."%
                              ", ".join(app_args.local_search_operators),end="")
                    sol = do_local_search(ls_ops, sol, D, d, C, L)
                    sol = cvrp_ops.normalize_solution(sol)
                        
                    if app_args.show_solution_cost:
                        ls_sol_q = cvrp_ops.calculate_objective(sol, D_c)
                    else:
                        ls_sol_q = cvrp_ops.calculate_objective(sol, D)
                    if ls_sol_q<sol_q:
                        if not app_args.minimal_output:
                            print(" improved by %.2f%%."%(1-ls_sol_q/sol_q))
                        sol_q = ls_sol_q
                        sol_K = sol.count(0)-1
                    else:
                        if not app_args.minimal_output:
                            print(" did not find improving moves.")
            else:
                sol_q = float('inf')

 
            elapsed_t = time()-start_t
            if app_args.minimal_output:
                print("%s;%s"%(algo_abbreviation, bn),end="")
                timecap_symbol = "*" if interrupted else ""
                if sol:
                    feasible = all( cvrp_ops.check_solution_feasibility( sol,D_c,d,C,L,st) )
                    print(";%s;%.2f;%d;%.2f%s"%
                          (str(feasible), sol_q, sol_K, elapsed_t, timecap_symbol))
                else:
                    print(";False;inf;inf;%.2f%s"%(elapsed_t, timecap_symbol))
                    
            elif sol:
                # Minimal output is not enabled, print like crazy :)
                
                if app_args.print_elapsed_time:
                    print("Algorithm produced a solution in %.2f s\n"%(elapsed_t))
                else:
                    #just a newline
                    print()
                    
                tightness = None
                if C and sol_K:
                    tightness = (sum(d)/(C*sol_K))
                if not bn in instance_data or sol_K<instance_data[bn][1]:
                    #"N K C tightness L st"
                    instance_data[bn] = (N,sol_K,C,"%.3f"%tightness,L,st)
                    
                shared_cli.print_problem_information(
                    points, D_c, d, C, L, st, tightness,
                    verbosity=app_args.verbosity)
                
                solution_print_verbosity = 3 if app_args.print_route_stat else 1
                shared_cli.print_solution_statistics(sol, D, D_c, d,
                                                     C, L, st,
                                                     solution_print_verbosity)
                    
                if app_args.print_vrph_sol:
                    print("SOLUTION IN VRPH FORMAT:")
                    print(" ".join( str(n) for n in cvrp_io.as_VRPH_solution(sol)))
                print("\n")
                
            short_algo_name = algo_name
            results[bn][short_algo_name] = sol_q
            
            if interrupted:
                break # algo loop
        if interrupted:
                break # problem file loop
            
    ## Print collected results 
    sys.stdout.flush()
    sys.stderr.flush()
    if not app_args.minimal_output and (len(results)>1 or len(algos)>1):
        print("\n")
        print_title = True
        ls_label = "+".join(ls_algo_names)
        for problem, algo_results in sorted(results.items()):
            algo_names = [ "%s+%s"%(algo_name,ls_label) if ls_algo_names else (algo_name)\
                           for algo_name in sorted(algo_results.keys())]
            
            if print_title:
                instance_fields = "instance\t"
                if PRINT_INSTANCE_DATA:
                    #"N K C tightness L st"
                    instance_fields+="N\tK*\tC\ttightness\tL\tst\t"
                print(instance_fields+"\t".join(algo_names))
                print_title=False
            print(problem,end="")
            if PRINT_INSTANCE_DATA:            
                print("\t",end="")
                print("\t".join( str(e) for e in instance_data[problem] ),end="")
            for _, result in sorted(algo_results.items()):
                print("\t", result,end="")
            print()
    sys.stdout.flush()
    sys.stderr.flush()

if __name__=="__main__":
    main()
