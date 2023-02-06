# -*- coding: utf-8 -*-
###############################################################################
""" This file implements tests comparing the operation of the actual optimized
local search procedures in inter_route_operations.py/intra_route_operations.py
with their naive implementation counterparts in naive_implementations.py

At the same time it to some extent tests that the do_local_search works 
properly.
NOTE: this does not mean parallel processing!

TODO: write dedicated tests for do_local_search at_lsop_optimal checks (if we
 correclty omit operations that are already at local optima).
"""
###############################################################################


# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

# syslib
import os
import sys
import unittest
import argparse
import logging
from subprocess import check_output
from tempfile import NamedTemporaryFile

# basics
import numpy as np 
from scipy.spatial.distance import pdist, squareform

# project imports
from verypy.cvrp_io import generate_CVRP, read_TSPLIB_CVRP
from verypy.cvrp_io import write_TSPLIB_file, as_OPT_solution
from verypy.cvrp_ops import normalize_solution, validate_solution_feasibility
from random import randint, shuffle
from verypy.local_search.naive_implementations import do_naive_local_search, \
    do_naive_2opt_move, do_naive_2optstar_move, \
    do_naive_1point_move, do_naive_relocate_move, \
    do_naive_exchange_move, do_naive_2point_move
from verypy.local_search import LSOPT, do_local_search
from verypy.local_search.intra_route_operators import do_2opt_move, do_3opt_move,\
    do_relocate_move, do_exchange_move
from verypy.local_search.inter_route_operators import do_2optstar_move,\
    do_1point_move, do_2point_move
                         
from verypy.util import objf, sol2routes, routes2sol
from verypy.config import BENCHMARKS_BASEPATH
from verypy.config import CAPACITY_EPSILON as C_EPS

DEBUG_VRPH_CALL = True
# this is to verify the heuristics against independent implementation
VRPH_SOLVER = r"C:\Users\juherask\Dissertation\Phases\Features\extractors\RKM16\solvers\vrp_load.exe"


def _do_vrph_ls(problem_file_name, solution_file_name, heuristics):
    sover_cmd_w_args = [VRPH_SOLVER,
          '-f', problem_file_name,
          '-s', solution_file_name,
          '-c', #use local search to do intra-route LS with ONE_POINT_MOVE+TWO_POINT_MOVE+TWO_OPT
          ]
    
    
    for h in heuristics:
        sover_cmd_w_args.append('-h')
        sover_cmd_w_args.append(h)          
        #'-h', "ONE_POINT_MOVE",
        #'-h', "TWO_POINT_MOVE",
        #'-h', "TWO_OPT",
        #'-h', "THREE_OPT"

    if DEBUG_VRPH_CALL:
        print(" ".join(sover_cmd_w_args))
    vrph_output = check_output(sover_cmd_w_args)
    
    routes = {}
    next_line_is_route_num = None
    for l in vrph_output.split("\n"):
        if next_line_is_route_num is not None:
            routes[next_line_is_route_num] = [int(n) for n in l.split("-")]
        if "routenum=" in l:
            # e.g.
            # Route 0002(routenum=2)[0-3...7-0, 10398.00, 4, 4]:
            parts = l.replace("(", " (").split()
            next_line_is_route_num = int(parts[1])
        else:
            next_line_is_route_num = None
        
    return routes2sol(routes.values())

class TestVsVRPH(unittest.TestCase):
    def setUp(self):
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(pts, "euclidean") )
        self.d = [1.0]*len(self.D)
        self.d[0] = 0.0
        self.C = float(len(self.D))
        
        self.longMessage = True
        
    def _test_ls_vs_vrhp_w_random_sol(self, vrph_heur, ls_heur,
                        times=10, routes=1):
        for i in range(times):
            initial_sol = list(range(1,len(self.D)))
            shuffle(initial_sol )
            initial_sol = [0]+initial_sol+[0]
            routes_to_create = routes-1
            while routes_to_create>0:
                possible_0_positions = len(initial_sol)-initial_sol.count(0)-2
                insert_0_counter = randint(0, possible_0_positions-1)
                for j in range(0, len(initial_sol)-1):
                    if initial_sol[j]!=0 and initial_sol[j+1]!=0:
                        # valid position to insert
                        
                        if insert_0_counter==0:
                            initial_sol.insert(j+1, 0)
                            break
                        else:
                            insert_0_counter-=1
                routes_to_create-=1
            #initial_sol = [0, 7, 3, 0, 4, 1, 5, 2, 6, 0]
            #print("initial_solution", initial_sol)
            #route = [max(0,int(n)-1) for n in "0-2-5-7-4-6-8-3-0".split("-")]
                            
            with NamedTemporaryFile( delete=False, suffix='.vrp') as tmpfile:
                tsplib_file_path = tmpfile.name
            write_TSPLIB_file(tsplib_file_path, self.D,
                              float_to_int_precision=1000)
            
            with NamedTemporaryFile( delete=False, suffix='.opt') as tmpfile:
                tmpfile.write( as_OPT_solution(objf(initial_sol, self.D), initial_sol) )
                opt_file_path = tmpfile.name

            vrph_sol = _do_vrph_ls(tsplib_file_path, opt_file_path, [vrph_heur])
            
            ls_sol = do_local_search([ls_heur], initial_sol,
                                     self.D, self.d, self.C,
                                     LSOPT.BEST_ACCEPT)
            
            # make sure the routes are in right order
            
            vrph_sol = normalize_solution(vrph_sol)
            ls_sol = normalize_solution(ls_sol)
            
            print("LS on", initial_sol,
                  "vrph_sol = %s (%.2f)"%(str(vrph_sol), objf(vrph_sol,self.D)),
                  "ls_sol = %s (%.2f)"%(str(ls_sol), objf(ls_sol,self.D)))
            self.assertEqual(vrph_sol, ls_sol)
            
            if not DEBUG_VRPH_CALL:
                os.remove(tsplib_file_path)
                os.remove(opt_file_path)
            
    def test_2opt_vs_vrph_from_single_route_random_solution_best_accept(self):
         self._test_ls_vs_vrhp_w_random_sol("TWO_OPT", do_2opt_move)
    
    def test_2opt_vs_vrph_from_multiroute_random_solution_best_accept(self):
        self._test_ls_vs_vrhp_w_random_sol("TWO_OPT", do_2opt_move, routes=2)
        self._test_ls_vs_vrhp_w_random_sol("TWO_OPT", do_2opt_move, routes=3)
    
    # Note, if this fails there is a bug in your VRPH implementation,
    #  the pull request fixing this bug can be found at
    #  https://github.com/geoffleyland/VRPH/pull/1
    def test_3opt_vs_vrph_from_single_route_random_solution_best_accept(self):
         self._test_ls_vs_vrhp_w_random_sol("THREE_OPT", do_3opt_move)
         
    def test_3opt_vs_vrph_from_multiroute_random_solution_best_accept(self):
        self._test_ls_vs_vrhp_w_random_sol("THREE_OPT", do_3opt_move, routes=2)
        self._test_ls_vs_vrhp_w_random_sol("THREE_OPT", do_3opt_move, routes=3)
        
    @unittest.skip("VRPH does not respect the INTER_ROUTE flag (there is a bug)")
    def test_1pm_vs_vrph_from_single_route_random_solution_best_accept(self):
         self._test_ls_vs_vrhp_w_random_sol("ONE_POINT_MOVE", do_relocate_move)
 
    #TestVsVRPH.test_2pm_vs_vrph_from_single_route_random_TSP_solution_best_accept
    #def test_2pm_vs_vrph_from_single_route_random_TSP_solution_best_accept(self):
    #     self._test_ls_vs_vrhp_w_random_sol("TWO_POINT_MOVE", do_exchange_move)

    #def test_2pm_vs_vrph_from_multiroute_random_solution_best_accept(self):
    #     self._test_ls_vs_vrhp_w_random_sol("TWO_POINT_MOVE", do_2point_move, routes=2)
    #     self._test_ls_vs_vrhp_w_random_sol("TWO_POINT_MOVE", do_2point_move, routes=3)
    
    
def _random_cvrp():
    size = randint(5,8)
    C = 50
    K = randint(2, max(3,int(size/3)))
    muC = K*C/size
    sdC = randint(2, max(3,int(C/K)))
    N, points, _, d, D, C, _ = generate_CVRP(size, C, muC, sdC)
    if randint(1,5)!=1:
        d = [int(e) for e in d]
    if randint(1,2)!=1:
        D = D.astype(int)
    return N, points, d, D, C

def _get_random_solution(d,C):
    nodes = list(range(1,len(d)))
    shuffle(nodes)
    sol = [0]
    route_d = 0
    for n in nodes:
        if route_d+d[n]+C_EPS>C:
            sol.append(0)
            route_d = 0
        route_d+=d[n]
        sol.append(n)
    if sol[-1]!=0:
        sol.append(0)
    
    if __debug__:
        print("random sol", sol)
    return sol

def _compare_improved_from_solution(testcase, sol, D,d,C,L,
                                           ls_ops, naive_ops,
                                           operator_strategy=LSOPT.BEST_ACCEPT,
                                           extra_msg=""):
        """ Improves the solution `sol` for the problem `D`,`d`,`C`,`L` using
        the local_search module operators `ls_ops` and naive implementations
        `naive_ops`. 
        """
        
        if __debug__:
            print("THE PROBLEM:")
            print("D,d,C,L =","np.%s"%repr(D),",",d,",",C,",",L)
        
        rls_ops = list(reversed(ls_ops))
        # note: if given multiple operators, the best single move out of among
        #  all moves is chosen at each step due to ITEROPT.BEST_ACCEPT.
        ls_sol_fwdop = do_local_search(ls_ops, sol, D, d, C, L=L,
                                 operator_strategy=operator_strategy,
                                 iteration_strategy=operator_strategy)
        ls_sol_rwdop = do_local_search(rls_ops, sol, D, d, C, L=L,
                                 operator_strategy=operator_strategy,
                                 iteration_strategy=operator_strategy)
        bf_sol = do_naive_local_search(naive_ops, sol, D, d, C, L=L,
                                       operator_strategy=operator_strategy)
        
        ls_sol_fwdop = normalize_solution(ls_sol_fwdop)
        ls_sol_rwdop = normalize_solution(ls_sol_rwdop)
        bf_sol = normalize_solution(bf_sol)
        
        if __debug__:
            print("\nFINAL SOLUTIONS:")
            print("+".join( op.__name__ for op in ls_ops),"alt 1 :",
                  ls_sol_fwdop, "(%.2f)"%objf(ls_sol_fwdop, D))
            print("+".join( op.__name__ for op in ls_ops),"alt 2 :",
                  ls_sol_rwdop, "(%.2f)"%objf(ls_sol_rwdop, D))
            print("+".join( op.__name__ for op in naive_ops),":", bf_sol,
                  "(%.2f)"%objf(bf_sol, D))
        
        testcase.assertTrue(all(validate_solution_feasibility(ls_sol_fwdop, D, d, C, L)))
        testcase.assertTrue(all(validate_solution_feasibility(ls_sol_rwdop, D, d, C, L)))
        testcase.assertTrue(all(validate_solution_feasibility(bf_sol, D, d, C, L)))
        testcase.assertTrue(ls_sol_fwdop==bf_sol or ls_sol_rwdop==bf_sol, extra_msg)   
           
class Test2Opt(unittest.TestCase):
    
    def setUp(self):
        self.N, self.points,self.d, self.D, self.C  = _random_cvrp()
        self.L = None
        self.longMessage = True
        
    def test_regression_cases(self):
        # only test this on the first iteration
        global iteration
        if iteration>0:
            return 
        
        # There was an bug in marking routes locally optimal in do_local_search
        #  that manifested when the search was terminated early i.e. when using
        #  LSOPT.FIRST_ACCEPT 
        # All of the remaining routes were incorrectly marked as optimal for
        #  the ls_op.
        D= np.array([[  0, 142, 322, 168, 176, 285, 160, 198],
                     [142,   0, 185, 298, 238, 298,  19, 153],
                     [322, 185,   0, 459, 352, 439, 172, 209],
                     [168, 298, 459,   0, 150, 416, 317, 275],
                     [176, 238, 352, 150,   0, 461, 255, 146],
                     [285, 298, 439, 416, 461,   0, 297, 441],
                     [160,  19, 172, 317, 255, 297,   0, 163],
                     [198, 153, 209, 275, 146, 441, 163,   0]]) 
        d,C,L = [0, 8, 31, 19, 4, 7, 13, 18], 50, None
        sol = [0, 4, 6, 7, 0, 2, 5, 1, 0, 3, 0]
        _compare_improved_from_solution(
            self, sol, D,d,C,L,
            [do_2opt_move], [do_naive_2opt_move],
            operator_strategy=LSOPT.FIRST_ACCEPT)
        
    def test_improve_random_solution_best_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,self.L,
            [do_2opt_move], [do_naive_2opt_move],
            operator_strategy=LSOPT.BEST_ACCEPT)
        
    def test_improve_random_solution_first_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,self.L,
            [do_2opt_move], [do_naive_2opt_move],
            operator_strategy=LSOPT.FIRST_ACCEPT)
        
class Test2OptStar(unittest.TestCase):
    
    def setUp(self):
        self.N, self.points, self.d, self.D, self.C  = _random_cvrp()
        self.longMessage = True
        
    def test_regression_cases(self):
        # only test this on the first iteration
        global iteration
        if iteration>0:
            return 
        
        ## Case 1
        D = np.array([[  0, 195, 168, 293, 236, 276],
                      [195,   0, 223, 225, 226, 434],
                      [168, 223,   0, 158, 377, 236],
                      [293, 225, 158,   0, 440, 380],
                      [236, 226, 377, 440,   0, 507],
                      [276, 434, 236, 380, 507,   0]])
        d,C = [0, 14, 1, 24, 50, 13], 50
        sol=[0, 5, 3, 0, 4, 0, 1, 2, 0]
        _compare_improved_from_solution(
            self, sol, D, d, C, None,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move])
        
        ## Case 2
        D = np.array([[  0,  65, 133, 299, 214, 169, 211, 299, 143],
                      [ 65,   0, 172, 259, 238, 233, 198, 363, 109],
                      [133, 172,   0, 281, 327, 185, 343, 291, 161],
                      [299, 259, 281,   0, 497, 447, 438, 567, 155],
                      [214, 238, 327, 497,   0, 208, 124, 286, 346],
                      [169, 233, 185, 447, 208,   0, 289, 129, 299],
                      [211, 198, 343, 438, 124, 289,   0, 392, 302],
                      [299, 363, 291, 567, 286, 129, 392,   0, 424],
                      [143, 109, 161, 155, 346, 299, 302, 424,   0]])
        d,C = [0, 26, 18, 21, 18, 24, 16, 16, 27], 50
        sol = [0, 6, 2, 0, 7, 5, 0, 8, 0, 1, 3, 0, 4, 0]
        _compare_improved_from_solution(
            self, sol, D, d, C, None,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move])
    
        
    def test_improve_random_solution_best_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation may coincidentally invert also routes"+
         " which are between the two routes that have their edges removed."+
         " Then subsequent moves may in some cases choose different equally"+
         " good improving moves and the final solutions differ. Both 2-opt*"
         " implementations still work properly.") 
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.BEST_ACCEPT)
        
        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.BEST_ACCEPT)
    
    @unittest.skip("The differences in different search order become more"
                   " pronounced when using FIRST_ACCEPT.")
    def test_improve_random_solution_first_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation may coincidentally invert also routes"+
         " which are between the two routes that have their edges removed."+
         " Then subsequent moves may in some cases choose different equally"+
         " good improving moves and the final solutions differ. Both 2-opt*"
         " implementations still work properly.") 
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.FIRST_ACCEPT)
        
        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_2optstar_move,do_2opt_move], [do_naive_2optstar_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.FIRST_ACCEPT)
    
        
class TestRelocateMove(unittest.TestCase):
    
    def setUp(self):
        self.N, self.points,self.d, self.D, self.C  = _random_cvrp()
        self.longMessage = True
        
    def test_improve_eil101_solution(self):
        # only test this on the first iteration
        global iteration
        if iteration>0:
            return
        
        pp = os.path.join( BENCHMARKS_BASEPATH, "Classic",
                        "ChristofidesEilon1969","10-eil101.vrp")
        N, points, dd_points, d, D, C, _ = read_TSPLIB_CVRP( pp )
        sol = [0, 40, 21, 58, 13, 0, 70, 30, 20, 66, 71, 65, 35, 9, 81, 33, 78, 79, 3, 24, 55, 0, 82, 48, 47, 36, 49, 64, 63, 90, 11, 19, 62, 10, 0, 4, 25, 39, 67, 23, 56, 75, 41, 22, 74, 72, 73, 91, 0, 57, 42, 14, 38, 44, 16, 86, 17, 46, 61, 85, 98, 37, 92, 6, 0, 60, 26, 12, 54, 80, 68, 29, 34, 51, 1, 50, 77, 76, 28, 0, 94, 95, 97, 87, 2, 15, 43, 100, 93, 59, 99, 96, 7, 0, 27, 69, 31, 32, 88, 8, 45, 84, 5, 83, 18, 52, 89, 53, 0]
    
        _compare_improved_from_solution(
            self, sol, D,d,C,None,
            [do_relocate_move, do_1point_move], [do_naive_1point_move])
   
    def test_improve_random_solution_best_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation checks first all insertion positions"+
         " for one customer before moving on to next, while the local search"+
         " Implementation checks all intra route moves first and then the"+
         " inter route moves. If there is equally good move in intra and inter"+
         " route moves, it may be that both implementations do not select the"+
         " same one.")
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_relocate_move, do_1point_move], [do_naive_1point_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.BEST_ACCEPT)
        
        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_relocate_move, do_1point_move], [do_naive_1point_move],
            extra_msg=special_case_fails_message+
            " Problem with L constraint %.2f."%max_L,
            operator_strategy=LSOPT.BEST_ACCEPT)
         
    @unittest.skip("The differences in different search order become more"
                   " pronounced when using FIRST_ACCEPT.")
    def test_improve_random_solution_first_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation checks first all insertion positions"+
         " for one customer before moving on to next, while the local search"+
         " Implementation checks all intra route moves first and then the"+
         " inter route moves. If there is equally good move in intra and inter"+
         " route moves, it may be that both implementations do not select the"+
         " same one.")
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_relocate_move, do_1point_move], [do_naive_1point_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.FIRST_ACCEPT)
        
        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_relocate_move, do_1point_move], [do_naive_1point_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.FIRST_ACCEPT)
        
class TestOnePointMove(unittest.TestCase):
    
     def setUp(self):
        self.N, self.points,self.d, self.D, self.C  = _random_cvrp()
        self.longMessage = True
        
     def test_improve_random_solution_best_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_relocate_move], [do_naive_relocate_move],
            operator_strategy=LSOPT.BEST_ACCEPT)
        
     def test_improve_random_solution_first_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_relocate_move], [do_naive_relocate_move],
            operator_strategy=LSOPT.FIRST_ACCEPT)
                           
class TestTwoPointMove(unittest.TestCase):
    
    def setUp(self):
        self.N, self.points,self.d, self.D, self.C  = _random_cvrp()
        self.longMessage = True
        
    def test_improve_random_solution_best_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation checks first all swap alternatives"+
         " for one customer before moving on to next, while the local search"+
         " Implementation checks all intra route swaps first and then the"+
         " inter route swaps. If there is equally good move in intra and inter"+
         " route moves, it may be that both implementations do not select the"+
         " same one.")
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_exchange_move, do_2point_move], [do_naive_2point_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.BEST_ACCEPT)

        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_exchange_move, do_2point_move], [do_naive_2point_move],
            extra_msg=special_case_fails_message+
            " Problem with L constraint %.2f."%max_L,
            operator_strategy=LSOPT.BEST_ACCEPT)
     
    @unittest.skip("The differences in different search order become more"
               " pronounced when using FIRST_ACCEPT.")
    def test_improve_random_solution_first_accept(self):
        special_case_fails_message = ("NOTE: This may sometimes fail due to"+
         " how the naive implementation checks first all swap alternatives"+
         " for one customer before moving on to next, while the local search"+
         " Implementation checks all intra route swaps first and then the"+
         " inter route swaps. If there is equally good move in intra and inter"+
         " route moves, it may be that both implementations do not select the"+
         " same one.")
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_exchange_move, do_2point_move], [do_naive_2point_move],
            extra_msg=special_case_fails_message,
            operator_strategy=LSOPT.FIRST_ACCEPT)

        # Test with L constraint        
        max_L = max( objf(r,self.D) for r in sol2routes(sol) )
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,max_L,
            [do_exchange_move, do_2point_move], [do_naive_2point_move],
            extra_msg=special_case_fails_message+
            " Problem with L constraint %.2f."%max_L,
            operator_strategy=LSOPT.FIRST_ACCEPT)

class TestExchangeMove(unittest.TestCase):
    def setUp(self):
        self.N, self.points,self.d, self.D, self.C  = _random_cvrp()
        self.longMessage = True
        
    def test_improve_random_solution_best_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_exchange_move], [do_naive_exchange_move],
            operator_strategy=LSOPT.BEST_ACCEPT)  
        
    def test_improve_random_solution_first_accept(self):
        sol = _get_random_solution(self.d,self.C)
        _compare_improved_from_solution(
            self, sol, self.D,self.d,self.C,None,
            [do_exchange_move], [do_naive_exchange_move],
            operator_strategy=LSOPT.FIRST_ACCEPT)

if __name__=="__main__":
    
    if __debug__:
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            stream=sys.stdout,level=logging.DEBUG)
        for lvl in range(1,10):
            logging.addLevelName(lvl, "DEBUG")
        
    # this allows setting the number of repeats (relevant )
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--repeat", dest="repeat",
        help="Repeat tests this many times (improves random solutions)")
    (args, unitargs) = parser.parse_known_args()
    unitargs.insert(0, "placeholder") # unittest ignores first arg
    
    # add more arguments to unitargs here
    repeat = vars(args)["repeat"]
    if repeat is None:
        repeat = 1
    else:
        repeat = int(repeat)
    for iteration in range(repeat):
        print("Test repeat %d/%d"%(iteration+1, repeat))
        wasSuccessful = unittest.main(exit=False, argv=unitargs).result.wasSuccessful()
        if not wasSuccessful:
            break
