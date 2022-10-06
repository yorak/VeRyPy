# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from os import path

from verypy.classic_heuristics.parallel_savings import parallel_savings_init
from verypy.classic_heuristics.sequential_savings import sequential_savings_init
from verypy.classic_heuristics.gaskell_savings import gaskell_savings_init
from verypy.classic_heuristics.paessens_savings import paessens_savings_init
from verypy.classic_heuristics.suppression_savings import suppression_savings_init

from verypy.local_search import LSOPT, do_local_search
from verypy.local_search.intra_route_operators import do_3opt_move

from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS

class TestSavingsGaskell1967Replications(ReplicationBase):
    
    def setUp(self):
        self.algorithms = [
            ("savings:multiple", lambda pts, D,d,C,L,st:\
                parallel_savings_init(D,d,C,L, minimize_K=True)),
            ("savings:sequential", lambda pts, D,d,C,L,st:\
                sequential_savings_init(D,d,C,L, minimize_K=True,
                                        initialize_routes_with = "farthest")),
            ("lambda:multiple", lambda pts, D,d,C,L,st:\
                gaskell_savings_init(D,d,C,L, minimize_K=True,
                                     savings_method="lambda")),
            ("pi:multiple", lambda pts,D,d,C,L,st:\
                gaskell_savings_init(D,d,C,L, minimize_K=True,
                                     savings_method="pi"))]

        self.problem_names =  ["01_Ga67_n37_k5.vrp",
                          "02-CW64_n32_k8.vrp",
                          "03-Ga67_n33-k4.vrp",
                          "04-Ga67-n22-k4.vrp",
                          "05-Ga67-n30-k4.vrp",
                          "06-Ga67-n23-k5.vrp"]

        self.targets = [
            #savings:multiple
            ((5,923),(8,1427),(5,839),(4,598),(4,963),(5,955)),
            #savings:sequential
            ((5,947),(8,1427),(5,850),(4,648),(5,1017),(5,949)),
            #lambda:multiple
            ((5,913),(8,1434),(5,821),(4,602),(5,979),(5,988)),
            #pi:multiple	    
            ((5,857),(9,1500),(5,850),(4,598),(5,943),(6,1015))]
            
        
        self.problem_path = path.join("Classic","Gaskell1967")
    
    def test_parallel_savings_with_Gaskell1967_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("savings:multiple")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 
    
    def test_sequential_savings_with_Gaskell1967_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("savings:sequential")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq)  
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 
     
    def test_Gaskell_lambda_savings_with_Gaskell1967_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("lambda:multiple")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 
    
    def test_Gaskell_pi_savings_with_Gaskell1967_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("pi:multiple")
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 

class TestSavingsCWEilonEtAl1971Replications(ReplicationBase):
    """ The results were published in Eilon et al 1971 (in the book Distribution management)
    """
    
    def setUp(self):
        self.algorithms = [
            ("savings:multiple", lambda pts, D,d,C,L,st:\
                parallel_savings_init(D,d,C,L)),]

        self.problem_names =  [
            "01-eil7.vrp",
            "02-eil13.vrp",
            "03-Gaskell1-E022-k4g_exact2d.vrp",
            "04-Gaskell2-E023-k5g_exact2d.vrp",
            "05-Gaskell3-E030-k4g_exact2d.vrp",
            "06-CW64_n31_k8c.vrp",
            "07-Gaskell4-E033-k4g_exact2d.vrp",
            "08-eil51_exact2d.vrp",
            "09-eil76_exact2d.vrp",
            "10-eil101_exact2d.vrp"]

        # CW64_n31_k8c is without Birmingham with demand of 140
        #  (it is always served separately with a route of length 82*2)
        self.targets = [
            #savings:multiple
            ((2,119),(4,290),(4,598),(5,955),(5,963),
             (8,1427-82*2),(5,839),(6,585),(10,900),(8,887))
            ]
            
        self.problem_path = path.join("Classic","Beasley1983")
        
    
    def test_parallel_savings_with_Beasley1983_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("savings:multiple", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 
    
    
class TestSavingsPaessensReplications(ReplicationBase):
    """ The results were published in Paessens1988 
    """
    
    def setUp(self):
        self.algorithms = [
            ("clarke_wright_savings", lambda pts, D,d,C,L,st:\
                parallel_savings_init(D,d,C,L)),
            ("paessens_savings_M1", lambda pts, D,d,C,L,st:\
                paessens_savings_init(D,d,C,L,strategy="M1", do_3opt=False)),
            ("paessens_savings_M4", lambda pts, D,d,C,L,st:\
                paessens_savings_init(D,d,C,L,strategy="M4", do_3opt=False)),
            ("clarke_wright_savings_3OPT", lambda pts, D,d,C,L,st:\
                do_local_search([do_3opt_move], parallel_savings_init(D,d,C,L),
                             D, d, C, L, operator_strategy=LSOPT.BEST_ACCEPT)),
            ("paessens_savings_M1_3OPT", lambda pts, D,d,C,L,st:\
                 paessens_savings_init(D,d,C,L,strategy="M1", do_3opt=True)),
            ("paessens_savings_M4_3OPT", lambda pts, D,d,C,L,st:\
                 paessens_savings_init(D,d,C,L,strategy="M4", do_3opt=True))
        ]

        self.problem_names =  [
            "G1.vrp",
            "G2.vrp",
            "G3.vrp",
            "G4.vrp",
            "C1.vrp",
            "C2.vrp",
            "C3.vrp",
            "C4.vrp",
            "C5.vrp",
            "C6.vrp",
            "C7.vrp",
            "C8.vrp",
            "C9.vrp",
            "C10.vrp",
            "C11.vrp",
            "GJ1.vrp" ]
            
        # in Gaskell instances service_time was not included in the table
        # in Christofides et al. service_time was included
        self.targets = [
            #cw
            (599-10*21,956-10*22,964-10*29,841-10*32,
             585,907,889,834,876,1068,1593,1140,1288,1395,1539,5568),
            #M1
            (585-10*21,956-10*22,938-10*29,814-10*32,
             564,866,866,826,870,1065,1584,1102,1222,1370,1486,5380),
            #M4
            (598-10*21,956-10*22,938-10*29,814-10*32,
             571,889,877,826,872,1068,1591,1112,1222,1370,1515,5476),
            #cw+3opt
            (599-10*21,956-10*22,962-10*29,840-10*32,
             579,902,880,824,869,1047,1583,1134,1285,1387,1522,5546),
            #M1+3opt
            (585-10*21,956-10*22,937-10*29,812-10*32,
             557,861,858,822,870,1046,1568,1083,1221,1359,1476,5348),
            #M4+3opt
            (598-10*21,956-10*22,937-10*29,812-10*32,
             570,876,869,823,871,1047,1580,1106,1221,1359,1504,5449),
        ]
                        
        self.problem_path = path.join("Classic","Paessens1988", "literature")
        
    
    def test_CW_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("clarke_wright_savings",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        
        # Is better, give it some slack
        self.assertTrue( abs(avgq)-0.02 < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 

    def test_CW_w_3OPT_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("clarke_wright_savings_3OPT",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        
        # It is in fact, better. But still, B.
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 

    def test_PaessensM1_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("paessens_savings_M1",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        self.assertTrue( abs(avgq)-0.01 < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 

    def test_PaessensM1_w_3OPT_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("paessens_savings_M1_3OPT",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 

    def test_PaessensM4_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("paessens_savings_M4",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 

    def test_PaessensM4_w_3OPT_with_Paessens1988_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems("paessens_savings_M4_3OPT",
                                                    round_f_func = _to_int,
                                                    cost_compare = False)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.A_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.A_SD, "There is too much variation between instances") 

class TestSavingsSuppressionEilonEtAl1971Replications(ReplicationBase):
    """ The results were published in Eilon et al 1971 (in the book Distribution management)
    """
    
    def setUp(self):
        self.algorithms = [
            (r"HP76-PS$\vert$IMS", lambda pts, D,d,C,L,st:\
                suppression_savings_init(D,d,C,L, minimize_K=False)),]

        self.problem_names =  [
            "01-eil7.vrp",
            # some 
            "08-eil51.vrp",
            "09-eil76.vrp",
            "10-eil101.vrp"]

        self.targets = [((2, 114), (5, 573), (10, 886), (8,876))]
            
        self.problem_path = path.join("Classic","ChristofidesEilon1969")
        
    
    def test_suppression_savings_with_ChristofidesEilon1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems(r"HP76-PS$\vert$IMS", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < REPRO_QUALITY_LEVELS.B_AVG, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < REPRO_QUALITY_LEVELS.B_SD, "There is too much variation between instances") 

      
if __name__ == '__main__':
    unittest.main(exit=False)     
            