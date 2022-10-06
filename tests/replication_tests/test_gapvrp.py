# -*- coding: utf-8 -*-

import unittest
import cvrp_io
from verypy.classic_heuristics.gapvrp import gap_init, _sweep_seed_points
import numpy as np
from collections import namedtuple
from os import path

from math import sqrt

from replicationbase import ReplicationBase

#TODO: use REPRO_QUALITY_LEVELS according to the replication levels of the 
# "Summary of replication results" table.
AVG_QUALITY_REQUREMENT = 0.2 # 0.1%
SD_QUALITY_REQUREMENT = 1.0 # 0.5%
RIGHT_SEED_POS_TOL = 5.0 # 5 % 

def _read_TSPLIB_seed_points(file_name):
    seed_coords = []
    with open(file_name, 'r') as f:
        seed_coord_section = False
        for l in f.readlines():
            if "SECTION" in l:
                seed_coord_section = "SEED_COORD_SECTION" in l
                continue        
            if seed_coord_section:
                parts = l.split()
                seed_coords.append( (float(parts[1]),float(parts[2])) )
    return seed_coords

LiteratureResult = namedtuple('LiteratureResult', 'obj_f cpu_time')
class TestGAPVRPAlgorithm(unittest.TestCase):
 
    def setUp(self):
        pass
 
    def test_seed_point_gen_vs_fig4_example(self):
        pfn = r"fisher_jaikumar_fig4.vrp"
        seed_example_problem_instance = cvrp_io.read_TSPLIB_CVRP(pfn)
        seed_example_seed_points = _read_TSPLIB_seed_points(pfn)
        
        N, points, dd_points, d, D, C, _ = seed_example_problem_instance
        generated_seeds = _sweep_seed_points(points, D, d, C, 3)
        
        #ran = np.amax(points)-np.amin(points)
        #tol = ran*(RIGHT_SEED_POS_TOL/100.0)
        
        avg_rel_d=0.0
        for gx,gy in generated_seeds:
            # find the closest match
            min_d = None
            min_idx = None
            for ri, (rx, ry) in enumerate(seed_example_seed_points):
                
                dd = sqrt( (points[0][0]-rx)**2+(points[0][1]-ry)**2 )
                tol = dd*(RIGHT_SEED_POS_TOL/100.0)
                
                d = sqrt( (gx-rx)**2+(gy-ry)**2 )
                avg_rel_d+=d/dd
                if (min_d is None) or (d<min_d):
                    min_d = d
                    min_idx = ri
            
            #print min_d, tol
            self.assertAlmostEqual(min_d, 0.0, delta=tol,
                msg="The seed point %d is over %.2f %% off" %
                    (min_idx+1, RIGHT_SEED_POS_TOL) )
        avg_rel_d=avg_rel_d/len(generated_seeds)
        print(avg_rel_d)

class TestFisherJaikumarReplications(ReplicationBase):
    
    def setUp(self):
        # note: use peredefined K's
        self.algorithms = [
            ("gapvrp_C_EXACT", lambda pts,D,d,C,L,st,K:\
              gap_init(pts, D, d, C, L, st=st,K=K, minimize_K=True,
              seed_edge_weight_type='EXACT_2D'
              )),
            
            ("gapvrp_CD_EXACT", lambda pts,D,d,C,L,st,K:\
              gap_init(pts, D, d, C, L, st=st,K=K, minimize_K=True,
              seed_edge_weight_type='EXACT_2D'
              )),

            ("gapvrp_C_FLOOR", lambda pts,D,d,C,L,st,K:\
              gap_init(pts, D, d, C, L, st=st,K=K, minimize_K=True,
              seed_edge_weight_type='FLOOR_2D'
              )),
            
            ("gapvrp_CD_FLOOR", lambda pts,D,d,C,L,st,K:\
              gap_init(pts, D, d, C, L, st=st,K=K, minimize_K=True,
              seed_edge_weight_type='FLOOR_2D'
              ))]

        self.problem_names =  [
            "Christofides_01.vrp",
            "Christofides_02.vrp",
            "Christofides_03.vrp",
            "Christofides_04.vrp",
            "Christofides_05.vrp",
            "Christofides_06.vrp",
            "Christofides_07.vrp",
            "Christofides_08.vrp",
            "Christofides_09.vrp",
            "Christofides_10.vrp",
            "Christofides_12.vrp",
            "Christofides_14.vrp"]

        self.targets =  [
                         # C instances
                         (524, 857, 833, 1014, 1420, None, None,
                          None, None, None, 824, None),
                         # CD instances
                         #(None, None, None, None, None, 560, 916,
                         # 885, 1230, 1518, None, 848)]
                         (None, None, None, None, None, None, None,
                          None, None, 1518, None, 848)]
        # duplicate for EXACT and FLOOR targets
	self.targets*=2

        
                         
        self.problem_path = path.join("Classic", "FisherJaikumar1981")
     
    # NOTES: The description of some details of the heuristic are vague in 
    #  (Fisher & Jaikumar 1981). This issue in replication of their 
    #  results has been also noted by Cordeau et al. (2002).
    # 
    # We have made following assumptions:
    #  * even if not explicitly mentioned, the best starting position
    #     for the sweep in the first clustering phase is searched. This
    #     is computationally quite demanding.
    #  * when growing the sector and the cojoined arc in the second
    #     the full (not factorial) weight of the customer nodes 
    #     with factorial weights is considered. This way the automatic
    #     seed point generation is within 2.5% of the seed point
    #     positions of example Figure 4 in (Fisher & Jaikumar 1981).
    #  * Despite euclidiean distances are mentioned, it seems exact 
    #     distances were not used not used. The results are 
    #     given as integers, so possibly some kind of unspecified 
    #     rounding was used. Rounding down (floor) seems to produce
    #     results closest to the ones given in the paper.       
    
    #TestFisherJaikumarReplications.test_FisherJaikumar_with_C_CMT_instances
    def test_FisherJaikumar_with_Floor_C_CMT_instances(self):        
        avgq, sdq, minq, maxq = self.solve_problems(
            "gapvrp_C_FLOOR", require_K = True,
            round_D_func = lambda D: D.astype(int),
            round_f_func = np.int,
            cost_compare = False)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 
    
    #TestFisherJaikumarReplications.test_FisherJaikumar_with_Floor_CD_CMT_instances
    def test_FisherJaikumar_with_Floor_CD_CMT_instances(self):        
        avgq, sdq, minq, maxq = self.solve_problems(
            "gapvrp_CD_FLOOR", require_K = True,
#            instance_idx= [3,4],
            round_D_func = lambda D: D.astype(int),
            round_f_func = np.int,
            cost_compare = False)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 

    #TestFisherJaikumarReplications.test_FisherJaikumar_with_Exact_C_CMT_instances
    def test_FisherJaikumar_with_Exact_C_CMT_instances(self):        
        avgq, sdq, minq, maxq = self.solve_problems(
            "gapvrp_C_EXACT", require_K = True,
            round_f_func = np.int,
            cost_compare = False)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 
    
    #TestFisherJaikumarReplications.test_FisherJaikumar_with_CD_CMT_instances
    def test_FisherJaikumar_with_Exact_CD_CMT_instances(self):        
        avgq, sdq, minq, maxq = self.solve_problems(
            "gapvrp_CD_EXACT", require_K = True,
            round_f_func = np.int,
            cost_compare = False)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 


if __name__ == '__main__':
    unittest.main()
