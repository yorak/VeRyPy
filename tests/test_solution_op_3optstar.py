# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 19:00:19 2018

@author: juherask
"""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division


import sys
import unittest
from time import time
import numpy as np

from local_search import LSOPT, do_local_search
from local_search.solution_operators import do_3optstar_move
from classic_heuristics.nearest_neighbor import nearest_neighbor_init

from cvrp_io import generate_CVRP
from cvrp_ops import check_solution_feasibility, calculate_objective
from util import sol2routes, routes2sol
from test_intra_route_local_search_operation import Test3Opt

PRINT_ONLY_FINAL_RESULT = True


def _intra_route_3optstar_call(sol, D, strategy=LSOPT.FIRST_ACCEPT):
     return do_3optstar_move(sol, D, [1.0]*len(D), len(D), None, strategy)

class TestIntraRouteMoves3OptStarSolutionOperator(Test3Opt):
    """ This tests all the move operators on a single route. The test reuses
    the test_local_search_operation.Test3Opt unit test.
    """
    def setUp(self):
        super(TestIntraRouteMoves3OptStarSolutionOperator, self).setUp()
        self.move_op = _intra_route_3optstar_call

#TODO: write this
#class TestInterRouteMoves3OptStarSolutionOperator(unittest.TestCase):
# e.g. 6 nodes on a circular formation, depot at the center but tad closer to one

def _strategy_to_str(strategy):
    if LSOPT.FIRST_ACCEPT:
        return "FIRST_ACCEPT"
    elif LSOPT.BEST_ACCEPT:
        return "BEST_ACCEPT"
    else:
        return "N/A"

class TestSmoke3OptStarSolutionOperator(unittest.TestCase):
    def setUp(self):
        self.D = np.array([[  0, 195, 168, 293, 236, 276],
                              [195,   0, 223, 225, 226, 434],
                              [168, 223,   0, 158, 377, 236],
                              [293, 225, 158,   0, 440, 380],
                              [236, 226, 377, 440,   0, 507],
                              [276, 434, 236, 380, 507,   0]])
        self.d = [0, 14, 1, 24, 50, 13]
        self.C = 50
        self.initial_sol=[0, 5, 3, 0, 4, 0, 1, 2, 0]

    def test_smoke(self):
        print("in", self.initial_sol)
        smoke_sol = do_local_search([do_3optstar_move],
                                    self.initial_sol,
                                    self.D, self.d, self.C)
        print("out", smoke_sol)

class TestRandomStressOn3OptStarSolutionOperator(unittest.TestCase):
    # abuse class variable to repeat with different problem sizes
    problem_size = 5
    
    def setUp(self):
        N = TestRandomStressOn3OptStarSolutionOperator.problem_size
        problem = generate_CVRP(N, 50, 10, 5)
        aN,pts,d,D,C = (problem.size, problem.coordinate_points,
                        problem.customer_demands, problem.distance_matrix,
                        problem.capacity_constraint)
        pts = None
        d = [int(dv) for dv in d]
        D = D.astype(int)
        problem = aN,pts,d,D,C
        
        self.naive_sol = routes2sol( [[n] for n in range(1,aN+1)] )
        self.nn_sol = nearest_neighbor_init(D, d, C)
        self.L = max( calculate_objective(r,D) for r in sol2routes(self.nn_sol) )
        self.N = aN
        self.problem = (aN,pts,d,D,C)
        
    def _improve_with_3opt_star(self, problem, solution, strategy):
        if len(problem)==5:
            N,pts,d,D,C = problem
            L = None
        else:
            N,pts,d,D,C,L = problem
        
        sol = list(solution)
        if not PRINT_ONLY_FINAL_RESULT:
            print("\nin", sol)
        total_t = 0
        while True:
            start_t = time()
            out_sol = do_3optstar_move(sol, D, d, C, L, strategy)
            elapsed_t = time()-start_t
            total_t+=elapsed_t
            if not PRINT_ONLY_FINAL_RESULT:
                print("elapsed %.2f s"%elapsed_t)
                print("out (%s)\n"%_strategy_to_str(strategy))
            if out_sol[1] == None:
                print("no more improvements found")
                if PRINT_ONLY_FINAL_RESULT:
                     print("final (%s)"%_strategy_to_str(strategy),
                           sol, calculate_objective(sol, D))
                print("total elapsed %.2f s"%total_t)
                break
            if not PRINT_ONLY_FINAL_RESULT:
                print(out_sol, calculate_objective(out_sol[0], D))
            sol = out_sol[0]
            
            self.assertTrue( all(check_solution_feasibility(sol, D, d, C, L)), "must be feasible")
    
    
    def test_random_problems_with_C_constraints_first_accept_from_naive_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        self._improve_with_3opt_star(self.problem, self.naive_sol, strategy=LSOPT.FIRST_ACCEPT)
        
    def test_random_problems_with_C_constraints_best_accept_from_naive_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        self._improve_with_3opt_star(self.problem, self.naive_sol, strategy=LSOPT.BEST_ACCEPT)
        
    def test_random_problems_with_L_constraints_first_accept_from_naive_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        problem_with_L = tuple( list(self.problem)+[self.L] )
        self._improve_with_3opt_star(problem_with_L, self.naive_sol, strategy=LSOPT.FIRST_ACCEPT)
        
    def test_random_problems_with_L_constraints_best_accept_from_naive_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        problem_with_L = tuple( list(self.problem)+[self.L] )
        self._improve_with_3opt_star(problem_with_L, self.naive_sol, strategy=LSOPT.BEST_ACCEPT)
        
    
    def test_random_problems_with_C_constraints_first_accept_from_nn_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        self._improve_with_3opt_star(self.problem, self.nn_sol, strategy=LSOPT.FIRST_ACCEPT)
        
    def test_random_problems_with_C_constraints_best_accept_from_nn_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        self._improve_with_3opt_star(self.problem, self.nn_sol, strategy=LSOPT.BEST_ACCEPT)
        
    def test_random_problems_with_L_constraints_first_accept_from_nn_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        problem_with_L = tuple( list(self.problem)+[self.L] )
        self._improve_with_3opt_star(problem_with_L, self.nn_sol, strategy=LSOPT.FIRST_ACCEPT)
        
    def test_random_problems_with_L_constraints_best_accept_from_nn_sol(self):
        print("\n\nTEST RANDOM PROBLEM WITH %d CUSTOMERS, NAIVE INITIAL SOLUTION"%self.N)
        problem_with_L = tuple( list(self.problem)+[self.L] )
        self._improve_with_3opt_star(problem_with_L, self.nn_sol, strategy=LSOPT.BEST_ACCEPT)
        
if __name__=="__main__":
    unittest.main()
    if len(sys.argv)<=1 or "TestRandomStressOn3OptStarSolutionOperator" in sys.argv:
        #for N in range(70,101,10):
        for N in range(5,15):
            TestRandomStressOn3OptStarSolutionOperator.problem_size = N
            wasSuccessful = unittest.main(exit=False, argv=["TestRandomStressOn3OptStarSolutionOperator"]).result.wasSuccessful()
            if not wasSuccessful:
                sys.exit(1)
        
