# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 08:52:30 2017

@author: juherask
"""

import unittest
from scipy.spatial.distance import pdist, squareform
from classic_heuristics.gillet_miller_sweep import gillet_miller_init, _shortest_path_through_nodes
from cvrp_ops import D2D_c, check_solution_feasibility

from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS

from os import path

class TestGilletMillerSweep(unittest.TestCase):
    
    def setUp(self):
        self.longMessage = True
        
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(pts, "euclidean") )
 
    def test_shortest_path_through_nodes_basic_routing(self):
        D = self.D
        path, l = _shortest_path_through_nodes(D, 0, 4, [1,2,3,5,6])
        
        self.assertEqual(path, [0,1,2,3,6,5,4], 
                         "wrong routing or missing some of the points")
        ref_l = D[0,1]+D[1,2]+D[2,3]+D[3,6]+D[6,5]+D[5,4]
        self.assertAlmostEqual(ref_l, l, "returns wrong cost for the path")
    
    def test_shortest_path_through_nodes_degenerate_cases(self):             
        D = self.D
        
        path, l = _shortest_path_through_nodes(D, 0, 1, [4])
        self.assertEqual(path, [0,4,1])
        self.assertEqual(l, D[0,4]+D[4,1])
        
        path, l = _shortest_path_through_nodes(D, 0, 1, [])
        self.assertEqual(path, [0,1])
        self.assertEqual(l, D[0,1])
        
        path, l = _shortest_path_through_nodes(D, 0, 0, [])
        self.assertEqual(path, [0,0])
        self.assertEqual(l, 0.0)
        
    def test_shortest_path_through_nodes_end_index_bigger(self): 
        D = self.D
        
        path, l = _shortest_path_through_nodes(D, 3, 2, [4])
        self.assertEqual(path, [3,4,2])
        self.assertEqual(l, D[3,4]+D[4,2])

        path, l = _shortest_path_through_nodes(D, 3, 2, [6,4,5])
        self.assertEqual(path, [3,4,5,6,2])
        ref_l = D[3,4]+D[4,5]+D[5,6]+D[6,2]
        self.assertEqual(l, ref_l)
    
    #TestGilletMillerSweep.test_regression_small_static_case_with_C_and_L_constraints
    def test_regression_small_static_case_with_C_and_L_constraints(self):
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        D = squareform( pdist(pts, "euclidean") )
        d = [0.0]+[1.0]*(len(D)-1)
        C = 4.0; L = 14.0; st = 2.0
        D_c = D2D_c(D, st)
        
        #import cvrp_io
        #cvrp_io.write_TSPLIB_file("tiny_7_pt_problem.vrp", D, d, C, L)
        
        sol = gillet_miller_init(pts,D_c,d,C,L)
        self.assertTrue(check_solution_feasibility(sol, D_c, d, C, L),"Should produce feasible solution")
        

class TestMillerSolutions(unittest.TestCase):
    """ Because replicating results was not succeedi
    """
    
class TestGilletMillerReplications(ReplicationBase):
    
    def setUp(self):
        reference_results = [
          {
          "name":"01-Gaskell1-E022-k4g",
          # visual solution from Gaskell 1967
          "bks":[0,6,1,2,5,7,9,
                 0,10,8,3,4,11,13,
                 0,12,15,18,20,17,
                 0,14,21,19,16,0],
          "bks_val":585,
          "bks_k":4,
          # best Sweep of Gillet and Miller 1974
          "ref_val":591,
          "ref_k":4
          },
          {
          "name":"02-Gaskell2-E023-k5g",
          # visual solution from Gaskell 1967
          "bks":[0,7,8,4,5,9,13,11,12,
                 0,14,17,15,16,3,2,
                 0,20,22,19,18,
                 0,1,6,10,
                 0,21,0],
          "bks_val":949,
          "bks_k":5, #!
          # best Sweep of Gillet and Miller 1974
          "ref_val":956,
          "ref_k":5
          },
          {
          "name":"03-Gaskell3-E030-k4g",
          # visual solution from Gaskell 1967
          "bks":[0,15,16,13,7,17,9,8,14,21,
                 0,10,11,12,23,18,19,
                 0,22,2,5,4,1,6,3,20,
                 0,26,28,27,29,25,24,0],
          "bks_val":876,
          "bks_k":4,
          # best Sweep of Gillet and Miller 1974
          "ref_sol":[0,18,10,11,12,9,17,7,13,16,15,
                     0,26,28,27,25,24,29,
                     0,20,3,6,1,4,5,2,22,
                     0,23,8,14,21,19,0],
          "ref_val":875,
          "ref_k":4
          },
          {
          "name":"04-Gaskell4-E033-k4g",
          # visual solution from Gaskell 1967
          "bks":[0,29,28,16,27,26,
                 0,15,17,25,24,23,22,20,21,19,18,
                 0,31,30,3,4,2,12,14,
                 0,1,11,5,6,7,32,8,9,10,13,0],
          "bks_val":813,
          "bks_k":4,
          # best Sweep of Gillet and Miller 1974
          "ref_val":810,
          "ref_k":4
          },
          
          
          {
          "name":"05-E051-k5",
          #  best known
          "bks":None,
          "bks_val":521,
          "bks_k":5,
          # best Sweep of Gillet and Miller 1974
          "ref_val":546,
          "ref_k":5
          },
          
          {
          "name":"06-E076-k15s",
          #  best known  (new in Gillet and Miller 1974)
          "bks":None,
          "bks_val":1127,
          "bks_k":15,
          # best Sweep of Gillet and Miller 1974
          "ref_val":1127,
          "ref_k":14
          },
                  
          {
          "name":"07-E076-k10s",
          #  best known
          "bks":None,
          "bks_val":830,
          "bks_k":10,
          # best Sweep of Gillet and Miller 1974
          "ref_val":865,
          "ref_k":10
          },
          
          {
          "name":"08-E076-k8s",
          #  best known (new in Gillet and Miller 1974)
          "bks":None,
          "bks_val":754,
          "bks_k":8,
          # best Sweep of Gillet and Miller 1974
          "ref_val":754,
          "ref_k":8
          },
                  
          {
          "name":"09-E076-k7s",
          #  best known (new in Gillet and Miller 1974)
          "bks":None,
          "bks_val":715,
          "bks_k":7,
          # best Sweep of Gillet and Miller 1974
          "ref_val":715,
          "ref_k":7
          },
                  
          {
          "name":"10-E101-k14s",
          #  best known (new in Gillet and Miller 1974)
          "bks":None,
          "bks_val":1170,
          "bks_k":14,
          # best Sweep of Gillet and Miller 1974
          "ref_val":1170,
          "ref_k":14
          },
                  
          {
          "name":"11-E101-k8",
          #  best known
          "bks":None,
          "bks_val":817,
          "bks_k":8,
          # best Sweep of Gillet and Miller 1974
          "ref_val":862,
          "ref_k":8
          },
                  
          {
          "name":"12-n250-k25",
          #  best known
          "bks":None,
          "bks_val":5794,
          "bks_k":25,
          # best Sweep of Gillet and Miller 1974
          "ref_val":5794,
          "ref_k":25
          },
                          
          ]
        self.algorithms = [
            ("gillet_miller_sweep", lambda pts,D,d,C,L,st:\
              gillet_miller_init(pts,D,d,C,L))]

        self.problem_names =  [ rr["name"]+".vrp" for rr in reference_results ]
        self.targets =  [ tuple( [(rr["ref_k"], rr["ref_val"]) for rr in reference_results] ) ]
        self.problem_path = path.join("Classic","GilletMiller1974")
        
    def test_GilletMiller_with_Gaskell_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("gillet_miller_sweep", round_f_func = lambda x: int(x))
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 

if __name__ == '__main__':
    unittest.main()