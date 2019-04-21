# -*- coding: utf-8 -*-

import unittest
from classic_heuristics.matchingvrp import mbsa_init
import numpy as np
from os import path

from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS
#from tsp_solvers.tsp_solver_ropt import solve_tsp_3opt as solve_tsp
from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp


#
#VRP-15
#[[0,91, 8,12,13,14,15,11,10,9,7,6,5,4,3,1,2,0], # 19->91!
#[0,120,105,106,107,104,103,116,99,101,102,0], #...102,1 -> ...,102,0!
#[0,76,79,53,55,58,56,60,63,66,64,62,61,65,59,57,54,52,0],
#[0,87,92,89,90,18,114,109,115,97,94,93,96,95,0],
#[0,118,108,17,16,91,25,22,24,27,33,30,31,34,36, 29,35,32,28,26,23,20,21,0], 
#[0,100,98,68,73,77,80,78,72,75,74,71,70,69,67,0],
#[0,88,82,111,86,85,112,84,113,83,117,81,119,0],
#[0,110,40,43,45,48,51,50,49,47,46,44,41,42,39,38,37,0]]
#
#VRP-17
#[[0,117,113,83,2,1,3,4,5,6,108,118,18,0],
#[0,107,104,103,69,70,67,98,110,115,116,100,99,0],
#[0,82,81,112,84,85,89,92,91,90,97,94,93,96,0],
#[0,119,120,105,106,101,102,95,87,86,111,88,0],
#[0,53,55,58,54,52,40,42,39,38,37,0],
#[0,114,109,8,12,13,14,15,11,10,9,7,0],
#[0,43,45,48,51,50,49,47,46,44,41,0],
#[0,28,31,27,30,33,34,36,35,32,29,0],
#[0,57,59,65,61,62,64,66,63,60,56,0], # 64, 61, 62, 64, -> 65, 61, 62, 64, !
#[0,73,71,74,72,75,78,80,79,77,76,68,0],
#[0,21,20,26,23,25,24,22,19,16,17,0]]



class TestGaskellReplications(ReplicationBase):
    
    def setUp(self):
        self.algorithms = [
            ("mbsa_vrp", lambda pts,D,d,C,L,st:\
              mbsa_init(D, d, C, L, solve_tsp=solve_tsp))]

        self.problem_path = path.join("Classic", "DesrochersVerhoog1989")
        
    def test_FisheJaikumar_with_Gaskell_instances(self):
        self.problem_names =  [ "01-Gaskell1-E022-k4g.vrp",
                                "02-Gaskell2-E023-k5g.vrp",
                                "03-Gaskell3-E030-k4g.vrp",
                                "04-Gaskell4-E033-k4g.vrp"]
        self.targets =  [(587, 970, 939, 833)]
        
        avgq, sdq, minq, maxq = self.solve_problems("mbsa_vrp",
            round_f_func = np.int, cost_compare = True, require_K=False,
            
            )
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.C_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.C_SD, "There is too much variation between instances") 

    def test_FisheJaikumar_with_CMT_instances(self):
        self.problem_names =  [ "05-Christofides-01.vrp",
                                "06-Christofides-02.vrp",
                                "07-Christofides-03.vrp",
                                "08-Christofides-04.vrp",
                                "09-Christofides-05.vrp",
                                "10-Christofides-06.vrp",
                                "11-Christofides-07.vrp",
                                "12-Christofides-08.vrp",
                                "13-Christofides-09.vrp",
                                "14-Christofides-10.vrp",
                                "15-Christofides-11.vrp",
                                "16-Christofides-12.vrp",
                                "17-Christofides-13.vrp",
                                "18-Christofides-14.vrp"]

        self.targets =  [(586, 885, 889, 1133, 1424,
                          593, 963, 914, 1292, 1559,
                          1058, 828, 1562, 882)]
        
        avgq, sdq, minq, maxq = self.solve_problems("mbsa_vrp",
            round_f_func = np.int, cost_compare = False, require_K = False)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.C_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.C_SD, "There is too much variation between instances") 


if __name__ == '__main__':
    unittest.main()