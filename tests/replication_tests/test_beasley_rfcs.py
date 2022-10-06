# -*- coding: utf-8 -*-

from verypy.classic_heuristics.rfcs import route_first_cluster_second_init
from replicationbase import ReplicationBase
from verypy.config import BENCHMARKS_BASEPATH
from verypy.tsp_solvers.tsp_solver_ropt import solve_tsp_ropt 
from verypy.tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp_giant_tour
#from verypy.tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp_giant_tour

from verypy.util import objf
import unittest
from os import path

#TODO: use REPRO_QUALITY_LEVELS according to the replication levels of the 
# "Summary of replication results" table.
# More slack than usual b/c the method is without stochasticity
# for 1 trial these are quadrupled, for 5 trials tripled and for 10 doubled.
# for 25 trials and for deterministic version they are used as is.
AVG_QUALITY_REQUREMENT = 0.5
SD_QUALITY_REQUREMENT = 2.0

def _to_int(x):
    return int(x)
        
def beasley_rfcs_init(D, d, C, L, trials):
    random_2opt_tsp_sol = lambda D, nodes: solve_tsp_ropt(D, nodes,
                   do_shuffle=True, do2opt=True, do3opt=False)
    
    best_sol = None
    best_f = None
    for t in range(trials):
        # 2-opt it
        sol = route_first_cluster_second_init(D, d, C, L, minimize_K=False,
                        tsp_gen_algo=random_2opt_tsp_sol)
        sol_f = objf(sol, D)
        if (best_sol is None) or (sol_f<best_f):
            best_sol = sol
            best_f = sol_f
    return best_sol

def best_solutions_from_the_paper(D, d, C, L):
    N = len(D)-1
    
    if N==6:
        # Problem 1 Routes:
        sol = [0,1,2,3,
        0,4,5,6,0]
    if N==12:
        # Problem 2 Routes:
        sol = [0,1,2,3,4,
        0,5,
        0,6,8,9,
        0,11,12,10,7,0]
    if N==21:
        # Problem 3 Routes:
        sol = [0,9,7,5,2,1,6,
        0,10,8,3,4,11,13,
        0,14,21,19,16,
        0,17,20,18,15,12,0]
    if N==22:
        # Problem 4 Routes:
        sol = [0,7,8,4,5,9,11,12,
        0,10,13,
        0,6,1,2,3,15,16,
        0,14,17,22,20,
        0,18,19,21,0] 
    if N==29:
        # Problem 5 Routes:
        sol = [0,22,2,5,4,1,6,3,20,
        0,21,14,8,9,17,7,13,16,15,
        0,18,23,12,11,10,19,
        0,26,28,27,25,24,29,0]
    if N==30:
        # Problem 6 Routes:
        sol = [0,2,1,20,12,17,
        0,21,30,
        0,18,8,25,
        0,19,10,26,
        0,3,4,6,5,11,16,15,27,23,
        0,22,28,24,
        0,29,13,7,9,14,0]
    if N==32:
        # Problem 7 Routes:
        sol = [0,14,13,10,9,8,32,11,12,2,1,
        0,6,7,5,4,3,30,31,
        0,29,28,16,27,26,
        0,18,19,21,20,22,23,24,25,17,15,0]
    if N==50:
        # Problem 8 Routes:
        sol = [0,5,49,10,45,33,39,30,34,50,9,
        0,12,17,37,15,44,42,19,40,41,4,47,
        0,18,13,25,14,24,43,6,
        0,27,48,23,7,26,8,31,28,22,1,32,46,
        0,11,2,3,36,35,20,29,21,16,38,0]
    if N==75:
        # Problem 9 Routes:
        sol = [0,45,29,27,13,54,52,34,
        0,46,8,19,59,14,35,7,
        0,53,11,66,65,38,
        0,58,10,31,55,25,9,39,72,
        0,50,18,24,49,16,33,
        0,63,23,56,41,64,42,43,1,73,
        0,6,22,62,2,68,75,
        0,51,3,44,32,40,12,17,
        0,26,67,4,
        0,30,48,21,69,61,28,74,
        0,5,47,36,71,60,70,20,37,15,57,0]
        
        # Problem 10 Routes:
    if N==100:
        sol = [0,2,57,42,100,85,91,44,38,14,43,15,41,22,73,21,40,
        0,89,18,60,83,8,45,17,86,16,61,84,5,99,96,6,
        0,94,95,59,93,98,37,92,97,87,13,58,
        0,53,28,26,12,76,50,1,69,27,
        0,52,7,82,48,47,46,36,49,64,11,19,62,88,
        0,31,10,63,90,32,66,65,71,20,30,70,
        0,77,3,79,33,81,51,9,35,34,78,29,24,80,68,
        0,54,4,55,25,39,67,23,56,75,74,72,0]
    
    biggest_n = max(sol)
    for i in range(biggest_n):
        if i not in sol:
            raise Exception("node %d missing from the %d customer solution"%(i,N))
 
    return sol

class TestBeasleyReplications(ReplicationBase):
        
    def setUp(self):        
        self.algorithms = [
            ("best_reference", lambda pts,D,d,C,L,st:\
              best_solutions_from_the_paper(D, d, C, L)),        
            ("deterministic1", lambda pts,D,d,C,L,st:\
              route_first_cluster_second_init(D, d, C, L, minimize_K=False,
                                              tsp_gen_algo=solve_tsp_giant_tour)),
            ("deterministic10", lambda pts,D,d,C,L,st:\
              route_first_cluster_second_init(D, d, C, L, minimize_K=False,
                                              tsp_gen_algo=solve_tsp_giant_tour)),
            ("1_trial", lambda pts,D,d,C,L,st:\
              beasley_rfcs_init(D, d, C, L, 1)),
            ("5_trials", lambda pts,D,d,C,L,st:\
              beasley_rfcs_init(D, d, C, L, 5)),
            ("10_trials", lambda pts,D,d,C,L,st:\
              beasley_rfcs_init(D, d, C, L, 10)),
            ("25_trials", lambda pts,D,d,C,L,st:\
              beasley_rfcs_init(D, d, C, L, 25))]


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
            "10-eil101_exact2d.vrp"
        ]

        # CW64_n31_k8c is without Birmingham with demand of 140
        #  (it is always served separately with a route of length 82*2)
        self.targets =  [
        # These are the best results given in table 1
        ((2,114),(4,290),(4,585),(5,956),(4,875),(7,1444-82*2),(4,822),(5,552),(11,884),(8,873)),
        # the best we can expect from the TSP is the same as 1 trial of
        #  Table 1 of Beasley (1983), but test also against 10 repts
        ((2,114),(4,296),(4,608),(6,1017),(4,879),(8,1524-82*2),(5,848),(5,564),(11,906),(8,902)),
        ((2,114),(4,290),(4,585),(5,968),(4,875),(7,1444-82*2),(5,814),(5,564),(11,895),(8,878)),
        # These are from the Table 1 of Beasley (1983)
        ((2,114),(4,296),(4,608),(6,1017),(4,879),(8,1524-82*2),(5,848),(5,564),(11,906),(8,902)),
        ((2,114),(4,290),(4,585),(6,994),(4,876),(7,1462-82*2),(5,815),(5,564),(11,895),(8,880)),
        ((2,114),(4,290),(4,585),(5,968),(4,875),(7,1444-82*2),(5,814),(5,564),(11,895),(8,878)),
        ((2,114),(4,290),(4,585),(5,956),(4,875),(7,1444-82*2),(4,822),(5,552),(11,884),(8,873))]

                         
        self.problem_path = path.join("Classic","Beasley1983")
    
    def test_reference_solutions_Beasley_RFCS_with_eil_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("best_reference", round_f_func = _to_int)
        self.assertTrue( abs(avgq) == 0.0, "Reference solution quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) == 0.0, "There should be no variation between instances") 
    
    #TestBeasleyReplications.test_Deterministic1_Beasley_RFCS_with_eil_instances
    def test_Deterministic1_Beasley_RFCS_with_eil_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("deterministic1", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*2, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT) )
        
    def test_Deterministic10_Beasley_RFCS_with_eil_instances(self):
        # NOTES: We are happy if the deterministic version is closte to the 
        #  specified quality target with 10 repetitions. The SD does not matter so much.
        avgq, sdq, minq, maxq = self.solve_problems("deterministic10", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*4, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT*2,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT*2) ) 
    
    
    def test_01_trial_Beasley_RFCS_with_eil_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("1_trial", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*4, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT*4,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT*4) ) 
    
    def test_05_trial__Beasley_RFCS_with_eil_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("5_trials", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*3, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT*3,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT*3) ) 
    
    def test_10_trial_Beasley_RFCS_with_eil_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("10_trials", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*2, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT*2,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT*2) ) 
        
    def test_25_trial_Beasley_RFCS_with_eil_instances(self):
        # NOTES: The Beasley 1983 is in fact, stochastic.      
        avgq, sdq, minq, maxq = self.solve_problems("25_trials", round_f_func = _to_int)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq)
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT,
                        "There is too much variation between instances"+\
                        "(%.2f vs. %.2f)"%(sdq, SD_QUALITY_REQUREMENT) ) 

    @unittest.skip("Run only when required")
    def test_stess__Beasley_RFCS_with_7th_eil_instance(self):
        while True:
            problem_name = self.problem_names[6]
            pfn = path.join(BENCHMARKS_BASEPATH,self.problem_path, problem_name)
            sol,sol_f,sol_c = self._solve_instance(
                self.algorithms[3][1],pfn)
            
            print(self.algorithms[3][0], problem_name, sol_c, "\n")
            
if __name__ == '__main__':
    unittest.main()
