""" These tests call every one of the classical heuristics included in the
VeRyPy library. It is used as a smoke test to see if every heuristic is working
properly and is able to solve simple problems. If you are able to run this
test, you have installed the VeRyPy and its dependencies correctly.
"""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from scipy.spatial.distance import pdist, squareform

from verypy.cvrp_ops import validate_solution_feasibility, D2D_c
from verypy.util import objf

SMOKE_TEST_VERBOSITY = 0

# TODO: test with a TSP problem

class TestSmokeWithSimple7pProblem(unittest.TestCase):
    def setUp(self):
        self.size = 7
        self.pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(self.pts, "euclidean") )
        self.d = [1.0]*len(self.D)
        self.d[0] = 0.0
        
        self.C = 4.0
        self.L = 14.0
        self.st = 2.0
        
        self.optf_C = objf([0,1,2,3,4,0,5,6,7,0], self.D)
        self.optf_L = objf([0, 3, 4, 0, 5, 6, 7, 0, 1, 2, 0], self.D)+self.st*self.size
        self.worstf_C = objf([0,1,0,2,0,3,0,4,0,5,0,6,0,7,0], self.D)
        self.worstf_L = self.worstf_C+self.st*self.size
        
        self.longMessage = True
        
    def _solve(self, name, desc, algof, C=None, L=None, st=None):
        d = self.d
        if C is None:
            d = None
        if st is not None:
            D_c = D2D_c(self.D, st)
        else:
            D_c = self.D
        
        for minimize_K in [False, True]:
            print(name, "(min_K=%s, C=%s, L=%s)"% (str(minimize_K),"%.1f"%C if C else "None", "%.1f"%L if L else "None"))
            try:
                sol = algof(self.pts, D_c, d, C=C, L=L, st=st, 
                            wtt="EXACT_2D", single=False, minimize_K=minimize_K)
                #print("done")
            except NotImplementedError:
                continue
            
            sol_f = objf(sol, D_c)
           
            print("SOLUTION %s (%.2f)"%(sol,objf(sol,D_c)))
            cover_ok,capa_ok,rlen_ok = validate_solution_feasibility(sol,D_c,d,C,L)
            self.assertTrue( cover_ok, str(sol)+" is not a valid solution")
            self.assertTrue( capa_ok, str(sol)+" violates C constraint" )
            self.assertTrue( rlen_ok, str(sol)+" violates L constraint"  )
            
            
            if L:
                self.assertGreaterEqual(sol_f, self.optf_L,
                    "Cannot be better than the optimal solution")
                self.assertLessEqual(sol_f, self.worstf_L,
                    "Must be better than the worst possible solution")
            else:
                self.assertGreaterEqual(sol_f, self.optf_C,
                    "Cannot be better than the optimal solution")
                self.assertLessEqual(sol_f, self.worstf_C,
                    "Must be better than the worst possible solution")
    
    #TestSmokeWithSimple7pProblem.test_sequential_cheapest_insertion
    def test_sequential_cheapest_insertion(self):
        from verypy.classic_heuristics.cheapest_insertion import get_si_algorithm
        self._solve(*get_si_algorithm(), C=self.C)
        self._solve(*get_si_algorithm(), L=self.L,st=self.st)
        self._solve(*get_si_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.    
    def test_parallel_cheapest_insertion(self):
        from verypy.classic_heuristics.cheapest_insertion import get_pi_algorithm
        self._solve(*get_pi_algorithm(), C=self.C)
        self._solve(*get_pi_algorithm(2), C=self.C)
        self._solve(*get_pi_algorithm(3), C=self.C) 
        
        self._solve(*get_pi_algorithm(), L=self.L,st=self.st)
        self._solve(*get_pi_algorithm(2), L=self.L,st=self.st)
        self._solve(*get_pi_algorithm(3), L=self.L,st=self.st)
        
        self._solve(*get_pi_algorithm(), C=self.C,L=self.L,st=self.st)
        self._solve(*get_pi_algorithm(2), C=self.C,L=self.L,st=self.st)
        self._solve(*get_pi_algorithm(3), C=self.C,L=self.L,st=self.st)
      
    #TestSmokeWithSimple7pProblem.
    def test_cmt_2phase_heuristic(self):
        from verypy.classic_heuristics.cmt_2phase import get_cmt2p_algorithm
        self._solve(*get_cmt2p_algorithm(), C=self.C)
        self._solve(*get_cmt2p_algorithm(), L=self.L,st=self.st)
        self._solve(*get_cmt2p_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.    
    def test_gapvrp_heuristic(self):
        from verypy.classic_heuristics.gapvrp import get_gap_algorithm
        self._solve(*get_gap_algorithm(), C=self.C)
        # "cones" seed generation does requires C, use kmeans instead
        self._solve(*get_gap_algorithm(seed_method="kmeans"), L=self.L,st=self.st)
        self._solve(*get_gap_algorithm(), C=self.C,L=self.L,st=self.st)
     
    #TestSmokeWithSimple7pProblem.
    def test_gaskell_savings(self):
        from verypy.classic_heuristics.gaskell_savings import get_gs_algorithm
        self._solve(*get_gs_algorithm(), C=self.C)
        self._solve(*get_gs_algorithm(), L=self.L,st=self.st)
        self._solve(*get_gs_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.        
    def test_sequential_savings(self):
        from verypy.classic_heuristics.sequential_savings import get_ss_algorithm
        self._solve(*get_ss_algorithm(lambda_multiplier="auto"), C=self.C)
        self._solve(*get_ss_algorithm(lambda_multiplier="auto"), L=self.L,st=self.st)
        self._solve(*get_ss_algorithm(lambda_multiplier="auto"), C=self.C,L=self.L,st=self.st)

    #TestSmokeWithSimple7pProblem.  
    def test_gillet_miller_sweep(self):
        from verypy.classic_heuristics.gillet_miller_sweep import get_gm_algorithm
        self._solve(*get_gm_algorithm(), C=self.C)
        self._solve(*get_gm_algorithm(), L=self.L,st=self.st)
        self._solve(*get_gm_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.  
    def test_lr3opt_sweep(self):
        from verypy.classic_heuristics.lr3opt import get_lr3opt_algorithm
        self._solve(*get_lr3opt_algorithm(), C=self.C)
        self._solve(*get_lr3opt_algorithm(), L=self.L,st=self.st)
        self._solve(*get_lr3opt_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.      
    def test_maximum_matching_heuristic(self):
        from verypy.classic_heuristics.matchingvrp import get_mm_algorithm
        self._solve(*get_mm_algorithm(), C=self.C)
        self._solve(*get_mm_algorithm(), L=self.L,st=self.st)
        self._solve(*get_mm_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.test_mole_jameson_insertion
    def test_mole_jameson_insertion(self):
        from verypy.classic_heuristics.mole_jameson_insertion import get_mj_algorithm
        self._solve(*get_mj_algorithm(), C=self.C)
        self._solve(*get_mj_algorithm(), L=self.L,st=self.st)
        self._solve(*get_mj_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.test_sequential_nearest_neighbor
    def test_sequential_nearest_neighbor(self):
        from verypy.classic_heuristics.nearest_neighbor import get_snn_algorithm
        self._solve(*get_snn_algorithm(), C=self.C)
        self._solve(*get_snn_algorithm(), L=self.L,st=self.st)
        self._solve(*get_snn_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.test_parallel_nearest_neighbor
    def test_parallel_nearest_neighbor(self):
        from verypy.classic_heuristics.nearest_neighbor import get_pnn_algorithm
        self._solve(*get_pnn_algorithm(), C=self.C)
        self._solve(*get_pnn_algorithm(2), C=self.C)
        self._solve(*get_pnn_algorithm(3), C=self.C)

        self._solve(*get_pnn_algorithm(), L=self.L,st=self.st)
        self._solve(*get_pnn_algorithm(2), L=self.L,st=self.st)
        self._solve(*get_pnn_algorithm(3), L=self.L,st=self.st)

        self._solve(*get_pnn_algorithm(), C=self.C,L=self.L,st=self.st)
        self._solve(*get_pnn_algorithm(2), C=self.C,L=self.L,st=self.st)
        self._solve(*get_pnn_algorithm(3), C=self.C,L=self.L,st=self.st)

    #TestSmokeWithSimple7pProblem.  
    def test_paessens_generalized_savings(self):
        from verypy.classic_heuristics.paessens_savings import get_gps_algorithm
        self._solve(*get_gps_algorithm(), C=self.C)
        self._solve(*get_gps_algorithm(), L=self.L,st=self.st)
        self._solve(*get_gps_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.  
    def test_clarke_wright_parallel_savings(self):
        from verypy.classic_heuristics.parallel_savings import get_ps_algorithm
        self._solve(*get_ps_algorithm(), C=self.C) 
        self._solve(*get_ps_algorithm(), L=self.L,st=self.st)
        self._solve(*get_ps_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.test_petal_heuristic
    def test_petal_heuristic(self):
        from verypy.classic_heuristics.petalvrp import get_ptl_algorithm
        self._solve(*get_ptl_algorithm(), C=self.C)  
        self._solve(*get_ptl_algorithm(), L=self.L,st=self.st)
        self._solve(*get_ptl_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.  
    def test_routefirst_clustersecond_heuristic(self):
        from verypy.classic_heuristics.rfcs import get_rfcs_algorithm
        self._solve(*get_rfcs_algorithm(), C=self.C)
        self._solve(*get_rfcs_algorithm(), L=self.L,st=self.st)
        self._solve(*get_rfcs_algorithm(), C=self.C,L=self.L,st=self.st)

    #TestSmokeWithSimple7pProblem.  
    def test_suppression_savings_heuristic(self):
        from verypy.classic_heuristics.suppression_savings import get_ims_algorithm
        self._solve(*get_ims_algorithm(), C=self.C)
        self._solve(*get_ims_algorithm(), L=self.L,st=self.st)
        self._solve(*get_ims_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.  
    def test_basic_sweep_heuristic(self):
        from verypy.classic_heuristics.sweep import get_swp_algorithm
        self._solve(*get_swp_algorithm(), C=self.C)
        self._solve(*get_swp_algorithm(), L=self.L,st=self.st)
        self._solve(*get_swp_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.test_tyagi_nearest_neighbor  
    def test_tyagi_nearest_neighbor(self):
        from verypy.classic_heuristics.tyagi_nearest_neighbor import get_ty_algorithm
        self._solve(*get_ty_algorithm(), C=self.C)
        self._solve(*get_ty_algorithm(), L=self.L,st=self.st)
        self._solve(*get_ty_algorithm(), C=self.C,L=self.L,st=self.st)
    
    #TestSmokeWithSimple7pProblem.  
    def test_wren_holliday_sweep(self):
        from verypy.classic_heuristics.wren_holliday_sweep import get_wh_algorithm
        self._solve(*get_wh_algorithm(), C=self.C) 
        self._solve(*get_wh_algorithm(), L=self.L,st=self.st)
        self._solve(*get_wh_algorithm(), C=self.C,L=self.L,st=self.st)
        
if __name__ == '__main__':
    if SMOKE_TEST_VERBOSITY:
        import logging
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.DEBUG-SMOKE_TEST_VERBOSITY)
        for lvl in range(1,10):
            logging.addLevelName(lvl, "DEBUG")
    unittest.main()