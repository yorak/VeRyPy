# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from classic_heuristics.wren_holliday_sweep import wren_holliday_init,\
                                _get_route_phi_range
from classic_heuristics.sweep import BEST_ALTERNATIVE
from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS

import numpy as np
from math import pi
from os import path

def _do_ovelap(r1_phi_range, r2_phi_range):
    return (r1_phi_range[0][0] < r2_phi_range[0][1] and
            r2_phi_range[0][0] < r1_phi_range[0][1]) or\
           (r1_phi_range[1][0] < r2_phi_range[1][1] and\
            r2_phi_range[1][0] < r1_phi_range[1][1])

class TestRoutePhiRangeCalculation(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.phis = np.array([
                -pi,     # 0
                -pi*0.9, # 1
                -pi*0.7, # 2
                -pi*0.4, # 3
                -pi*0.1, # 4
                -pi*0.1, # 5
                -pi*0.0, # 6
                pi*0.2,  # 7
                pi*0.3,  # 8
                pi*0.6,  # 9
                pi*0.6,  # 10
                pi*0.9,  # 11
                pi ])    # 12
    
        self.rhos = np.array([10, 3, 6, 9, 2, 3, 7, 7, 5, 4, 9, 5, 9])
        
    def test_self_overlap(self):
        #lower half
        r1_phi_range = _get_route_phi_range(self.phis[[1,2,5]])
        self.assertTrue( _do_ovelap(r1_phi_range, r1_phi_range),
                         "route1 does not overlap with itself" )
        
        #upper half
        r2_phi_range = _get_route_phi_range(self.phis[[6,10]])
        self.assertTrue( _do_ovelap(r2_phi_range, r2_phi_range),
                         "route2 does not overlap with itself" )
        
    def test_on_different_halves(self):
        r1_phi_range = _get_route_phi_range(self.phis[[1,2,5]])
        r2_phi_range = _get_route_phi_range(self.phis[[6,10]])
        
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should NOT overlap" )
        
        
    def test_clear_cut_cases_lower_half(self):
        ## lower circle half, tight 
        r1_phi_range = _get_route_phi_range(self.phis[[1,2,5]])
        r2_phi_range = _get_route_phi_range(self.phis[[3,4]])
        r3_phi_range = _get_route_phi_range(self.phis[[0,1,2]])
        
        self.assertTrue( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertFalse( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should NOT overlap" )
        
    def test_clear_cut_cases_upper_half(self):
        ## upper circle half, tight 
        r1_phi_range = _get_route_phi_range(self.phis[[6,10]])
        r2_phi_range = _get_route_phi_range(self.phis[[8,9,11,12]])
        r3_phi_range = _get_route_phi_range(self.phis[[6,7]])
        
        self.assertTrue( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertFalse( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should NOT overlap" )
    
    def test_concave_cases(self):
        pass
    
    def test_over_origo_case(self):
        pass           
    
    def test_wrapping_case(self):
        # one wrap from upper to lower
        r1_phi_range = _get_route_phi_range(self.phis[[11,2]])
        r2_phi_range = _get_route_phi_range(self.phis[[1,3,4]])
        r3_phi_range = _get_route_phi_range(self.phis[[3,4,5]])
        r4_phi_range = _get_route_phi_range(self.phis[[9,10,8]])
        self.assertTrue( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route3 and route2 should overlap" )
        self.assertFalse( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should NOT overlap" )
        self.assertFalse( _do_ovelap(r1_phi_range, r4_phi_range),
                         "route1 and route4 should NOT overlap" )
        self.assertFalse( _do_ovelap(r2_phi_range, r4_phi_range),
                         "route2 and route4 should NOT overlap" )
        self.assertFalse( _do_ovelap(r3_phi_range, r4_phi_range),
                         "route3 and route4 should NOT overlap" )
        
        # one wrap from lower to upper
        r1_phi_range = _get_route_phi_range(self.phis[[2,11]])
        r2_phi_range = _get_route_phi_range(self.phis[[1,3,4]])
        r3_phi_range = _get_route_phi_range(self.phis[[3,4,5]])
        r4_phi_range = _get_route_phi_range(self.phis[[9,10,8]])
        self.assertTrue( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route3 and route2 should overlap" )
        self.assertFalse( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should NOT overlap" )
        self.assertFalse( _do_ovelap(r1_phi_range, r4_phi_range),
                         "route1 and route4 should NOT overlap" )
        self.assertFalse( _do_ovelap(r2_phi_range, r4_phi_range),
                         "route2 and route4 should NOT overlap" )
        self.assertFalse( _do_ovelap(r3_phi_range, r4_phi_range),
                         "route3 and route4 should NOT overlap" )
        
        #TODO: two wraps by one route from lower to upper
        
        #TODO: two wraps by one route from upper to lower
        
        #TODO: two wrap both routes from upper to lower
        
    def test_same_border_not_overlap(self):
        # two routes share a phi, upper half
        r1_phi_range = _get_route_phi_range(self.phis[[11,9]])
        r2_phi_range = _get_route_phi_range(self.phis[[7,10,8]])
        r3_phi_range = _get_route_phi_range(self.phis[[12,10,9,8]])
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "while sharing a phi value, route1 and route2 should NOT overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should overlap" )
    
        # two routes share a phi, lower half
        r1_phi_range = _get_route_phi_range(self.phis[[1,2,3,4]])
        r2_phi_range = _get_route_phi_range(self.phis[[5,6,7,8]])
        r3_phi_range = _get_route_phi_range(self.phis[[3,7,8,9]])
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "while sharing a phi value, route1 and route2 should NOT overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should overlap" )

        # at origo
        r1_phi_range = _get_route_phi_range(self.phis[[4,5,6,3]])
        r2_phi_range = _get_route_phi_range(self.phis[[10,7,6,8]])
        r3_phi_range = _get_route_phi_range(self.phis[[2,7,6,9]])
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "while sharing a phi value, route1 and route2 should NOT overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should overlap" )

        # over wrap
        r1_phi_range = _get_route_phi_range(self.phis[[1,12,2]])
        r2_phi_range = _get_route_phi_range(self.phis[[11,0,10,8]])
        r3_phi_range = _get_route_phi_range(self.phis[[1,11]])
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "while sharing a phi value, route1 and route2 should NOT overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route1 and route3 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should overlap" )
    
    def test_over_180_case(self):
        r1_phi_range = _get_route_phi_range(self.phis[range(2,11)])
        r2_phi_range = _get_route_phi_range(self.phis[[1,0,12,11]])
        r3_phi_range = _get_route_phi_range(self.phis[range(8,13)+range(0,4)])
        r4_phi_range = _get_route_phi_range(self.phis[[7,6,5,4]])
        
        self.assertFalse( _do_ovelap(r1_phi_range, r2_phi_range),
                         "route1 and route2 should NOT overlap" )
        self.assertFalse( _do_ovelap(r3_phi_range, r4_phi_range),
                         "route3 and route4 should NOT overlap" )
        self.assertFalse( _do_ovelap(r2_phi_range, r4_phi_range),
                         "route1 and route2 should NOT overlap" )
        
        self.assertTrue( _do_ovelap(r1_phi_range, r3_phi_range),
                         "route3 and route4 should overlap" )
        self.assertTrue( _do_ovelap(r1_phi_range, r4_phi_range),
                         "route3 and route4 should overlap" )
        self.assertTrue( _do_ovelap(r2_phi_range, r3_phi_range),
                         "route2 and route3 should overlap" )       
    
    def test_thin_at_wrap_or_origo(self):
        #TODO: implementing this would warrant changes into the comparison
        # criteria. Now a pathological case of a route being on the merge line
        # or at the origo (all points of the route on this line) fails to 
        # be overlapped by a route going over this line, even though it should!
        # However, this is extremely rare case. 
        pass
        
    def test_over_360_case(self):
        # zips back and forth, covering entire area
        r1_phi_range = _get_route_phi_range(self.phis[ [1,6,10,0,2,9] ])
        print(r1_phi_range)
        for i in range(13):
            for j in range(13):
                phi1,phi2 = self.phis[ [i,j] ]
                if phi1==phi2 or (abs(phi1)==pi and abs(phi2)==pi):
                    #ignore this very,very special case
                    continue
                pair_phi_range = _get_route_phi_range([phi1,phi2])
                #print(pair_phi_range)
                self.assertTrue( _do_ovelap(r1_phi_range, pair_phi_range),
                         "route1 is over 2*pi wide and should overlap everyting"+
                         "including i=%d, j=%d"%(i,j))       
    
def fill_missing_pts_as_needed(points, D):
    #print(points, dd_points)
    if points is None:
        # We do not have point coodrinates, but we have D!
        from sklearn import manifold
        mds = manifold.MDS(n_components=2, dissimilarity='precomputed',
                           random_state=42)
        mds_results = mds.fit(D)
        return list( mds_results.embedding_ )
    return points
        
    
        
class TestWrenHollidayReplications(ReplicationBase):
    
    def setUp(self):
        self.algorithms = [
            ("wren_holliday_G", lambda pts, D,d,C,L,st:\
                wren_holliday_init(fill_missing_pts_as_needed(pts,D), D,d,C,L,
                                   full_convergence = True)),
            ("wren_holliday_H", lambda pts, D,d,C,L,st:\
                wren_holliday_init(fill_missing_pts_as_needed(pts,D), D,d,C,L,
                                   full_convergence = False)),
            ("wren_holliday_BEST", lambda pts, D,d,C,L,st:\
                wren_holliday_init(fill_missing_pts_as_needed(pts,D), D,d,C,L,
                                   seed_node=BEST_ALTERNATIVE,
                                   full_convergence = True))]

        self.problem_names =  ["01_Ga67_n37_k5.vrp",
                               #"02-CW64_n31_k8c.vrp",
                               "02-CW64_n32_k9c.vrp",
                               "03-Ga67_n33_k4.vrp",
                               "04-Ga67-n22-k4.vrp",
                               "05-Ga67-n30-k4.vrp",
                               "06-Ga67-n23-k5.vrp",
                               "08-eil51.vrp",
                               "09-eil76.vrp",
                               "10-eil101.vrp"]

        self.targets = [
            #full convergence
            ((5,851),
            (8,1406),#-82*2),
            (4,812),
            (4,593),
            (4,888),
            (5,954),
            #--
            (5,551),
            (10,863),
            (8,851)),

            #terminating according to convergence
            ((5,851),
            (8,1406),#-82*2),
            (4,812),
            (4,593),
            (4,888),
            (5,964),
            #--
            (5,551),
            (10,900),
            (8,851)),
             
            # full convergence, all possible initial solutions
            ((5,842),
            (8,1406),#-82*2),
            (4,812),
            (4,593),
            (4,888),
            (5,954),
            #--
            (5,551),
            (10,863),
            (8,851))]
            
        self.problem_path = path.join("Classic","WrenHolliday1972")
        
    def test_WrenHolliday_full_convergence(self):
        avgq, sdq, minq, maxq = self.solve_problems("wren_holliday_G")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.C_SD, "There is too much variation between instances") 
    
    def test_WrenHolliday_early_termination(self):
        avgq, sdq, minq, maxq = self.solve_problems("wren_holliday_H")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.C_SD, "There is too much variation between instances") 

    #TestWrenHollidayReplications.test_WrenHolliday_all_initial_solutions
    def test_WrenHolliday_all_initial_solutions(self):
        # Target is really given only for the very first
        avgq, sdq, minq, maxq = self.solve_problems("wren_holliday_BEST", instance_idx=0)
        # Compute the others just in case
        all_avgq, all_sdq, all_minq, all_maxq = self.solve_problems("wren_holliday_BEST")
        
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.C_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.C_SD, "There is too much variation between instances") 
        
        
        
if __name__ == '__main__':
    unittest.main()