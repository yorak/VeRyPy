# -*- coding: utf-8 -*-

import unittest
from collections import namedtuple
from os import path
from random import shuffle

import numpy as np

from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS
from verypy.classic_heuristics.lr3opt import lr3opt_init, _check_lr3opt_move, _init_with_random
from verypy.cvrp_ops import calculate_objective
    
def _random_init_lr3opt(pts,D,d,C,L,st,times):
    best_sol = None
    best_f = float('inf')
    
    for t in range(times):
        sol = lr3opt_init(D, d, C, L, initialization_algorithm=_init_with_random)
        sol_f = calculate_objective(sol, D)
        if sol_f<best_f:
            best_sol = None
            best_f = None
    return best_sol

def _random_init_lr3opt_once(pts,D,d,C,L,st):
    return lr3opt_init(D, d, C, L, initialization_algorithm=_init_with_random)

LiteratureResult = namedtuple('LiteratureResult', 'obj_f cpu_time')
class TestLR3OPTAlgorithm(unittest.TestCase):
 
    def setUp(self):
        pass
 
    def test_penalty_calculation_fig1_example(self):
        D = np.array([
               [ 0, 19, 39, 51, 66, 59, 42, 22, 30, 40, 54, 68, 73, 62, 41],
               [19,  0, 21, 37, 54, 52, 36, 21, 36, 45, 58, 69, 77, 74, 54],
               [39, 21,  0, 21, 37, 40, 30, 27, 39, 46, 59, 65, 76, 81, 63],
               [51, 37, 21,  0, 17, 20, 19, 31, 34, 37, 48, 49, 61, 75, 60],
               [66, 54, 37, 17,  0, 16, 28, 45, 43, 43, 50, 45, 59, 79, 67],
               [59, 52, 40, 20, 16,  0, 17, 37, 30, 28, 34, 30, 43, 63, 53],
               [42, 36, 30, 19, 28, 17,  0, 19, 15, 18, 30, 34, 45, 55, 41],
               [22, 21, 27, 31, 45, 37, 19,  0, 15, 24, 37, 48, 56, 55, 35],
               [30, 36, 39, 34, 43, 30, 15, 15,  0, 10, 22, 34, 41, 42, 26],
               [40, 45, 46, 37, 43, 28, 18, 24, 10,  0, 13, 25, 32, 37, 24],
               [54, 58, 59, 48, 50, 34, 30, 37, 22, 13,  0, 16, 19, 28, 24],
               [68, 69, 65, 49, 45, 30, 34, 48, 34, 25, 16,  0, 13, 40, 40],
               [73, 77, 76, 61, 59, 43, 45, 56, 41, 32, 19, 13,  0, 32, 39],
               [62, 74, 81, 75, 79, 63, 55, 55, 42, 37, 28, 40, 32,  0, 20],
               [41, 54, 63, 60, 67, 53, 41, 35, 26, 24, 24, 40, 39, 20,  0]])
        C = 140
        #d = [30,15,15,15,15,15,15, #route 3
        #     10,20,20,20,30,30,20] #route 1
        
        #sol = [0,1,2,3,4,5,6,7,0,8,9,10,11,12,13,14,0]
        self.assertAlmostEqual( -10, _check_lr3opt_move(D, C, None, 60, 0,
                     [[0,4],[3,1],[2,5]], 
                     [6,7,8,9,9,10], #end_p
                     [6,7,0,8,8,9], #end_n
                     [105,15,0,10,10,140], #cum_d
                     None,
                     8, None, [2.0, None]))
        
class TestStewartGoldenReplications(ReplicationBase):
    
    def setUp(self):
        self.algorithms = [
            ("lr3opt_det", lambda pts,D,d,C,L,st:\
               lr3opt_init(D, d, C, L)),
            ("lr3opt_ran", _random_init_lr3opt_once)]

        self.problem_names =  [
            "00-CW64_n31_k8c.vrp",
            "05-E051-k5.vrp",
            "06-E076-k15s.vrp",
            "07-E076-k10s.vrp",
            "08-E076-k8s.vrp",
            "09-E076-k7s.vrp",
            "10-E101-k14s.vrp",
            "11-E101-k8.vrp"]
            
        self.targets =  [(1212,521,1058,847,751,692,1117,829), #det
                         (1212,521,1058,847,751,692,1117,829)] #rnd
                         
        self.problem_path = path.join("Classic", "GilletMiller1974")
    
    def test_deterministic_LR3OPT_with_GilletMiller1974_instances(self):   
        avgq, sdq, minq, maxq = self.solve_problems(
            "lr3opt_det", require_K = False,
            round_f_func = np.int,
            cost_compare = False)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 

    #TestStewartGoldenReplications.test_stochastic_LR3OPT_with_GilletMiller1974_instances
    def test_stochastic_LR3OPT_with_GilletMiller1974_instances(self):   
        
        repeats_per_problem = zip(list(range(8)), [10, 10, 7, 8*3, 10, 10*2, 3, 6*2])
        

        
        bestqs = [float('inf')]*8
        for i, repeats in repeats_per_problem: 
            for repeat in range(repeats):
                problem_name = self.problem_names[i]
                print "Repeat %d of %d for %s"%(repeat+1,repeats,problem_name)
                avgq, sdq, minq, maxq = self.solve_problems(
                    "lr3opt_ran", instance_idx = i, require_K = False,
                    round_f_func = np.int,
                    #round_D_func = np.around,
                    cost_compare = False)
                
                if avgq<bestqs[i]:
                    bestqs[i] = avgq
    
        # check the average gap of [10, 10, 7, 8, 10, 10, 3, 6] repeats
        avgq = np.average(bestqs)
        sdq = np.std(bestqs)     
        
        # Usually this assertion succeeds, but because it is stochastic, it 
        #  is possible that by chance some of the results are (much) worse.
        #  Then, it is best to try again on bump the level up to B.
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 


if __name__ == '__main__':
    unittest.main()