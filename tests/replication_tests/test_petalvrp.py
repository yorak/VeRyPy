# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from os import path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from verypy.cvrp_io import read_TSPLIB_CVRP, read_TSBLIB_additional_constraints
from verypy.cvrp_ops import check_solution_feasibility, D2D_c
from verypy.util import objf, totald
from replicationbase import ReplicationBase
from verypy.classic_heuristics.petalvrp import petal_init, _remove_multiserved, \
                                        _generate_solution_relaxation_petals
from verypy.config import BENCHMARKS_BASEPATH

MAX_REQUIRE_ITERATIONS = 250

#TODO: use REPRO_QUALITY_LEVELS according to the replication levels of the 
# "Summary of replication results" table.
AVG_QUALITY_REQUREMENT = 0.1 # 0.1%
SD_QUALITY_REQUREMENT = 0.5 # 0.5%

def flatten(routes):
    return [0]+[n for r in routes for n in r[1:] if r!=[0,0]]

class TestMutiservedRemoval(unittest.TestCase):
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
        self.longMessage = True
    
    def test_no_doubles(self):
        target_routes = [[0,1,2,3,4,0],[0,5,6,7,0]]
        modified_routes = [[0,1,2,3,4,0],[0,5,6,7,0]]
        _remove_multiserved(modified_routes, self.D)
        self.assertEqual( flatten(target_routes), flatten(modified_routes),
                          "There is no customers served multiple times, the "+
                          "routes should be left untouched" )
    
    def test_single_doubles(self):
        target_routes = [[0,1,2,3,4,0],[0,5,6,7,0]]
        modified_routes = [[0,1,2,3,4,0],[0,5,1,6,7,0]]
        _remove_multiserved(modified_routes, self.D)
        self.assertEqual( flatten(target_routes), flatten(modified_routes),
                          "The double served should be removed" )
    def test_mutiple_doubles(self):
        target_routes = [[0,1,2,3,4,0],[0,5,6,7,0]]
        modified_routes = [[0,1,2,5,3,4,0],[0,1,5,6,7,0]]
        _remove_multiserved(modified_routes, self.D)
        self.assertEqual( flatten(target_routes), flatten(modified_routes),
                          "The double served should be removed" )
    
    def test_mutiple_doubles_and_triples(self):
        target_routes = [[0,3,4,0],[0,1,2,5,6,7,0]]
        modified_routes = [[0,1,2,3,0],[0,3,4,2,5,0],[0,1,2,5,6,7,0]]
        _remove_multiserved(modified_routes, self.D)
        self.assertEqual( flatten(target_routes), flatten(modified_routes),
                          "All customers served multiple times should be removed" )
    
class TestFosterRyanRelaxationOperation(unittest.TestCase):
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
        self.longMessage = True

    def _make_improving_1pm_move(self, C=None, d=None, L=None):
        D = self.D
        r1 = [0,4,5,6,7,0]
        r2 = [0,1,2,3,0]
        return _generate_solution_relaxation_petals([r1,r2], None, D, C=C, d=d, L=L)
        
    def test_improving_1pm_move(self):
        result= self._make_improving_1pm_move() 
        #print("test_improving_move RESULT :",[rd.route for rd in result],len(result))
        routes = [rd.route for rd in result]
        self.assertTrue( [0,1,2,3,4,0] in routes, "node 4 should be inserted to r2")
        self.assertTrue( [0,5,6,7,0] in routes, "node 4 should be removed from r1")
    
    def test_non_improving_1pm_move(self):
        D = self.D
        r1 = [0,5,6,7,0]
        r2 = [0,1,2,3,4,0]
        result = _generate_solution_relaxation_petals([r1,r2], None, D)
        self.assertTrue( len(result)==0, "r1, r2 configuration should already be at local optima")
        result = _generate_solution_relaxation_petals([r2,r1], None, D)
        self.assertTrue( len(result)==0, "r1, r2 configuration should already be at local optima")
    
    def test_1pm_respects_capacity_constraint(self):
        d=[1]*len(self.D)
        result = self._make_improving_1pm_move(C=3.5, d=d)      
        self.assertTrue( len(result)==0, "no move should be made due to violation of C constraint")
            
    def test_1pm_respects_route_cost_constraint(self):
        r2 = [0,1,2,3,0]
        result = self._make_improving_1pm_move(L=objf(r2,self.D)+1.0)    
        self.assertTrue( len(result)==0, "no move should be made due to violation of C constraint")
    
    def test_1pm_updated_route_cost(self):
        result = self._make_improving_1pm_move()
        for rd in result:
            self.assertAlmostEqual( objf(rd[0],self.D), rd[1], msg="original route cost + savings should match recalculated route cost")
    
    def test_1pm_updated_route_capacity(self):
        # route 2 capacity 3.0 , route 1 capacity 4.0 (in _make_improving_move)
        #  node #4 is moved and with it the demand of 1.4 from r1 to r2
        d=[0.0, 0.7,1.0,1.3, 1.4,1.5,0.5,0.6]
        result = self._make_improving_1pm_move(C=5.0, d=d)
        for rd in result:
            self.assertAlmostEqual( totald(rd[0],d), rd[2], msg="original route cost + savings should match recalculated route cost")

    def _make_improving_2pm_move(self, C=None, d=None, L=None):
        D = self.D
        r1 = [0,4,5,6,0]
        r2 = [0,1,2,3,7,0]
        return _generate_solution_relaxation_petals([r1,r2], None, D, d=d, C=C, L=L)
        
    def test_2pm_with_C_constraint_no_improving(self):
        #     #0   #1   #2   #3   #4   #5   #6   #7
        d = [0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5]
        result = self._make_improving_2pm_move(d=d, C=3.0) 
        print("test_exchage_with_C_constraint_no_improving", [rd.route for rd in result])
        self.assertTrue( len(result)==0, "There should be no improving moves")
    
    def test_2pm_with_C_constraint_improving(self):
        # exhange is initiated only on C violation
        d = [0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5]
        result = self._make_improving_2pm_move(d=d, C=3.5) 
        routes = [rd.route for rd in result]
        self.assertTrue( [0,1,2,3,4,0] in routes, "nodes 4 and 7 should be swapped")
        self.assertTrue( [0,5,6,7,0] in routes, "nodes 4 and 7 should be swapped")
    
    def test_2pm_non_improving_move(self):
        r1 = [0,5,6,7,0]
        r2 = [0,1,2,3,4,0]
        result = _generate_solution_relaxation_petals([r1,r2], None, self.D)
        self.assertTrue( len(result)==0, "r1, r2 configuration should already be at local optima")
        result = _generate_solution_relaxation_petals([r2,r1], None, self.D)
        self.assertTrue( len(result)==0, "r1, r2 configuration should already be at local optima")
        
    def test_2pm_respect_capacity_constraint(self):
        d=[1]*len(self.D)
        result = self._make_improving_2pm_move(C=3.5, d=d)        
        self.assertTrue( len(result)==0, "no move should be made due to violation of C constraint")
        
    def test_2pm_respect_route_cost_constraint(self):
        r2 = [0,1,2,3,0]
        result = self._make_improving_2pm_move(L=objf(r2,self.D)+1.0)    
        self.assertTrue( len(result)==0, "no move should be made due to violation of C constraint")
    
    def test_2pm_updated_route_cost(self):
        result = self._make_improving_exchage_move()
        for rd in result:
            self.assertAlmostEqual( objf(rd[0],self.D), rd[1], msg="original route cost + savings should match recalculated route cost")
    
    def test_2pm_updated_route_capacity(self):
        # route 2 capacity 3.0 , route 1 capacity 4.0 (in _make_improving_move)
        #  node #4 is moved and with it the demand of 1.4 from r1 to r2
        d=[0.0, 0.7,1.0,1.3, 1.4,1.5,0.5,0.6]
        result = self._make_improving_2pm_move(C=5.0, d=d)
        for rd in result:
            self.assertAlmostEqual( totald(rd[0],d), rd[2], msg="original route cost + savings should match recalculated route cost")

     
class TestFosterRyanReplications(ReplicationBase):
    
    def setUp(self):
        self.algorithms = [
             ("petal_lp", lambda pts,D,d,C,L,st,K:\
               petal_init(pts, D, d, C, L,
                          minimize_K = True,
                          relaxe_SCP_solutions=False,
                          required_iterations=1,
                          restricted_route_ratio=None,
                          K_constraint=K)),
             ("relaxations", lambda pts,D,d,C,L,st:\
               petal_init(pts, D, d, C, L,
                          minimize_K = True,
                          relaxe_SCP_solutions=True)),
             ("user_set", lambda pts,D,d,C,L,st:\
               petal_init(pts, D, d, C, L,
                          minimize_K = True,
                          relaxe_SCP_solutions=True,
                          min_iterations=50)),
             ("iterations", lambda pts,D,d,C,L,st, K:\
               #ABUSE K parameter to pass the specific number of iterations
               petal_init(pts, D, d, C, L,
                          minimize_K = True,
                          relaxe_SCP_solutions=True,
                          #restricted_route_ratio=None,
                          required_iterations=K))
             ]

        self.problem_names =  [
            r"01-Ga67-04_n22_k4cd.vrp",
            r"02-Ga67-06_n23_k5cd.vrp",
            r"03-Ga67-05_n30_k4cd.vrp",
            # Omit 04 altogehter, different instance!
            #r"04-CW64_n30a_k8c.vrp",
            #r"04-CW64_n30b_k8c.vrp",
            #r"04-CW64_n31_k9c.vrp",
            r"05-Ga67-03_n33_k4cd.vrp",
            r"06-Ga67-01_n37_k5cd.vrp",
            r"07-E051-k5c.vrp",
            r"08-E076-k15c.vrp",
            r"09-E076-k10c.vrp",
            r"10-E076-k8c.vrp",
            r"11-E076-k7c.vrp",
            r"12-FR076-k10cd.vrp",
            r"13-E101-k14c.vrp",
            r"14-E101-k8c.vrp",
            r"15-FR101-k8cd.vrp"]
        
        self.targets = [
            # petal_lp
           (( 4,607 ),
            ( 5,863 ),
            ( 4,813 ),
            # Omit 04 altogehter, different instance!
            #( 7,1427-82*2),
            #( 8,1427),
            #( 8,1427),
            ( 4,792 ),
            ( 5,841 ),
            ( 5,523 ),
            (14,1084),
            (10,864 ),
            ( 8,768 ),
            ( 7,692 ),
            (10,852 ),
            (14,1174),
            ( 8,825 ),
            ( 8,827 )),
            # relaxations
           (( 4,585 ),  
            ( 5,953 ),  
            ( 4,873 ),  
            # Omit 04 altogehter, different instance!
            #( 7,1377-82*2),  
            #( 8,1377),
            #( 8,1377),
            ( 4,809 ),  
            ( 5,838 ),  
            ( 5,521 ),  
            (14,1081),  
            (10,852 ),  
            ( 8,760 ),  
            ( 7,692 ),  
            (10,865 ),  
            (14,1116),  
            ( 8,825 ),  
            ( 8,826 )),
            # user_set (70 min_iterations)
           (( 4,585 ),  
            ( 5,953 ),  
            ( 4,873 ),  
            # Omit 04 altogehter, different instance!
            #( 7,1377-82*2),  
            #( 8,1377),
            #( 8,1377),
            ( 4,809 ),  
            ( 5,838 ),  
            ( 5,521 ),  
            (14,1081),  
            (10,852 ),  
            ( 8,760 ),  
            ( 7,692 ),  
            (10,865 ),  
            (14,1116),  
            ( 8,825 ),  
            ( 8,826 )),
           # iterations
           #TODO: get from Azure the 06,07,09 if able (1000e test running)
           (( 28,585 ), #01 gap=0.0
            ( 141,953 ),  #02 gap<0.0 (better, no closer result)
            ( 65,873 ),  #03 gap=0.0
            # Omit 04 altogehter, different instance!
            #( 37,1377-82*2), #04 
            #( 8,1377),
            #( 8,1377),
            ( 28,809 ),  #05 gap=0.0 
            ( 1,838 ),   #06 gap>0.0 !
            ( 23,521 ),  #07 gap>0.0 !
            ( 60,1081),  #08 gap<0.0 (better, no closer result)
            ( 40,852 ),  #09 gap>0.0 !
            ( 1,760 ),   #10 gap<0.0 (better, after 1st iter)
            ( 1,692 ),   #11 gap=0.0 
            ( 144,865 ), #12 gap~0.0 (~200 better)
            ( 6,1116),   #13 gap<0.0 (better, 4 worse)
            ( 5,825 ),   #14 gap=0.0
            ( 107,826 ))]#15 gap=0.0
        
        # The solutions from Foster & Ryan 1976 APPENDIX
        self.target_solutions = [
        # 1
        [0,6,1,2,5,7,9,
         0,10,8,3,4,11,13,
         0,17,20,18,15,12,
         0,14,21,19,16,0],
        # 2
        [0,2,3,16,15,17,14,
         0,7,8,5,4,21,
         0,12,1,6,13,9,
         0,11,10,
         0,20,22,19,18,0],
        # 3
        [0,22,2,5,4,1,6,3,20,
         0,15,16,13,7,17,9,12,11,10,18,
         0,26,28,27,25,24,29,
         0,23,8,14,21,19,0],
        # 4, different instance, uncomment to verify
        #[0,24,29,23,
        # 0,17,20,1,2,
        # 0,26,10,19,
        # 0,6,5,13,7,11,16,15,9,14,27,12,
        # 0,25,8,18,
        # 0,21,30,
        # 0,4,28,3,22,
        # 0,31,0],
        # 5
        [0,1,11,5,6,7,9,8,10,32,13,
         0,14,17,25,24,23,20,22,21,18,19,
         0,12,2,4,3,30,31,
         0,15,26,27,28,16,29,0],
        # 6
        [0,9,10,11,18,12,6,5,4,3,
        0,8,13,19,25,31,32,26,20,14,
        0,7,1,2,
        0,16,22,28,34,33,27,21,15,
        0,23,29,35,36,30,24,17,0],
        # 7
        [0,27,48,23,7,43,24,25,14,6,
        0,47,4,17,42,19,40,41,13,18,
        0,11,2,29,21,16,50,34,30,9,38,
        0,46,5,49,10,39,33,45,15,44,37,12,
        0,32,1,22,20,35,36,3,28,31,26,8,0],
        # 8
        [0,17,3,24,49,23,16,51,
        0,63,56,41,43,33,
        0,12,65,66,11,
        0,6,2,28,61,74,
        0,46,8,35,7,67,
        0,75,68,29,45,4,
        0,40,55,25,9,39,
        0,15,57,59,14,53,
        0,36,71,60,70,20,37,5,
        0,44,18,50,32,26,
        0,72,31,10,38,58,
        0,34,52,27,13,54,19,
        0,21,69,47,48,30,
        0,62,22,64,42,1,73,0],
        # 9
        [0,4,45,5,15,57,27,52,
        0,29,37,20,70,60,71,69,36,47,48,
        0,75,62,22,64,42,41,43,1,
        0,7,53,65,38,10,58,72,26,
        0,35,14,59,66,11,
        0,40,9,25,55,31,39,12,
        0,51,16,49,56,23,63,73,33,6,
        0,67,46,8,19,54,13,34,
        0,44,32,50,18,24,3,17,
        0,68,2,28,61,21,74,30,0],
        # 10
        [0,1,43,41,42,64,22,28,62,2,
        0,6,73,33,63,23,56,24,49,16,51,17,
        0,40,9,25,55,18,50,32,44,3,
        0,30,48,47,36,71,60,70,69,61,21,74,68,
        0,12,72,39,31,10,65,38,58,26,
        0,35,14,59,66,11,53,7,
        0,75,4,45,29,5,37,20,15,57,27,
        0,34,52,13,54,19,8,46,67,0],
        # 11
        [0,40,39,9,25,55,31,10,58,72,12,
        0,68,2,62,28,61,21,74,30,75,
        0,63,23,56,41,43,42,64,22,1,73,33,6,
        0,7,53,14,59,66,11,65,38,26,
        0,34,52,27,15,57,13,54,19,35,8,46,67,
        0,45,29,5,37,20,70,60,71,69,36,47,48,4,
        0,17,3,44,32,50,18,24,49,16,51,0],
        # 12
        [0,30,48,47,36,60,71,69,21,74,
        0,6,62,61,28,2,68,75,
        0,29,5,37,70,20,15,45,4,
        0,7,11,66,65,38,26,
        0,9,25,55,18,50,32,17,
        0,51,16,63,23,56,49,24,3,44,
        0,12,72,58,10,31,39,40,
        0,35,19,59,14,53,67,
        0,22,64,42,41,43,1,73,33,
        0,34,52,27,57,13,54,8,46,0],
        # 13
        [0,52,7,48,82,8,45,17,84,18,
        0,61,16,86,38,44,91,98,
        0,13,97,92,100,37,95,94,
        0,87,42,14,43,15,57,2,58,53,
        0,19,11,64,49,36,47,46,83,
        0,26,54,55,25,24,80,68,12,
        0,40,21,72,75,41,22,74,73,
        0,99,93,85,59,96,
        0,9,71,65,35,81,33,50,
        0,88,62,63,90,32,10,31,
        0,4,39,67,23,56,
        0,89,60,5,6,
        0,3,79,78,34,29,77,76,28,
        0,1,51,66,20,30,70,69,27,0],
        #14
        [0,26,4,56,23,67,39,25,55,54,12,
        0,94,95,97,92,98,37,100,14,38,43,42,87,13,
        0,27,31,10,32,90,63,64,49,19,11,62,88,
        0,58,2,57,15,41,22,75,74,72,73,21,40,53,
        0,99,93,61,16,86,44,91,85,59,96,6,
        0,50,33,81,51,9,35,65,71,66,20,30,70,1,69,
        0,28,76,77,3,79,78,34,29,24,80,68,
        0,52,7,82,48,47,36,46,8,45,17,84,5,60,83,18,89,0],
        # 15
        [0,52,7,82,48,47,36,46,8,45,17,84,5,60,83,18,89,
        0,28,76,77,3,79,33,81,78,34,29,24,80,68,
        0,31,10,32,90,63,64,49,19,11,62,88,
        0,99,93,61,16,86,44,91,85,59,96,6,
        0,50,51,9,35,65,71,66,20,30,70,1,69,27,
        0,58,2,57,15,41,22,75,74,72,73,21,40,53,
        0,94,95,97,92,98,37,100,14,38,43,42,87,13,
        0,26,4,56,23,67,39,25,55,54,12,0]]
        
        self.problem_path = path.join("Classic", "FosterRyan1976")
    
    #TestFosterRyanReplications.test_verify_reference_solutions_FosterRyan1976_instances
    def test_verify_reference_solutions_FosterRyan1976_instances(self):
        for problem_idx, problem_name in enumerate(self.problem_names):
            ref_k, ref_f = self.targets[1][problem_idx]
            
            if problem_name==r"04-CW64_n30a_k8c.vrp":
                problem_name=r"04-CW64_n31_k9c.vrp"
                ref_f=1377
                
            pfn = path.join(BENCHMARKS_BASEPATH,self.problem_path, problem_name)
            N, points, dd_points, d, D, C, _ = read_TSPLIB_CVRP(pfn)
            K, L, service_time = read_TSBLIB_additional_constraints(pfn)
            if service_time:
                D_c = D2D_c(D, service_time)
            else:
                D_c = D
            
            ref_sol = self.target_solutions[problem_idx]
            ref_sol_f = int(objf(ref_sol, D_c))
            ref_sol_k = ref_sol.count(0)-1
            
            cover_ok, capa_ok, rlen_ok = check_solution_feasibility(ref_sol, D,d,C,L,True)
            self.assertTrue( cover_ok, "Must be a valid solution")
            self.assertTrue( capa_ok, "Must not violate the C constraint" )
            self.assertTrue( rlen_ok, "Must not violate the L constraint"  )
            
            self.assertEqual(ref_k, ref_sol_k,
             msg=("The appendix solution route count differs from the one given "+
                  "in Table 2 for %s (%d vs %d)"%(problem_name,ref_sol_k,ref_k)))
            self.assertAlmostEqual(ref_f, ref_sol_f,
             msg=("The appendix solution result differs from the one given "+
                  "in Table 2 for %s : %d (ours) vs %d (theirs)"%(problem_name,ref_sol_f,ref_f)))
    
    
    #TestFosterRyanReplications.test_petal_lp_with_FosterRyan1976_instances
    def test_petal_lp_with_FosterRyan1976_instances(self):   
        avgq, sdq, minq, maxq = self.solve_problems(
            "petal_lp", require_K = True,
            round_f_func = np.int,
            cost_compare = True,
            suppress_constraint_check=True)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 

    #TestFosterRyanReplications.test_relaxations_with_FosterRyan1976_instances
    def test_relaxations_with_FosterRyan1976_instances(self):   
        avgq, sdq, minq, maxq = self.solve_problems(
            "relaxations", #require_K = True,
            round_f_func = np.int,
            cost_compare = True)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 

    #TestFosterRyanReplications.test_user_set_with_FosterRyan1976_instances
    def test_user_set_with_FosterRyan1976_instances(self):   
        avgq, sdq, minq, maxq = self.solve_problems(
            "user_set", #require_K = True,
            round_f_func = np.int,
            cost_compare = True)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 

    #TestFosterRyanReplications.test_iterated_with_FosterRyan1976_instances
    @unittest.skip("Not reported in the _Foster Ryan (1976) paper")
    def test_iterated_with_FosterRyan1976_instances(self):   
        avgq, sdq, minq, maxq = self.solve_problems(
            "iterations", require_K = True,
            round_f_func = np.int,
            cost_compare = True)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 



if __name__ == '__main__':
    unittest.main()
