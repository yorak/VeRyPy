# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from replicationbase import ReplicationBase
from classic_heuristics.cmt_2phase import cmt_2phase_init
from os import path
import unittest
from util import objf

from random import random

#TODO: use REPRO_QUALITY_LEVELS according to the replication levels of the 
# "Summary of replication results" table.
AVG_QUALITY_REQUREMENT = 0.1
SD_QUALITY_REQUREMENT = 0.5

RANDOM_TRALS = 100
RANDOMIZED_RETRIES = 10

#stochastic_params = {
#    (51,160,None ):(2.60,1.32),
#    (76,140,None ):(2.59,0.87),
#    (151,200,None):(2.99,1.92),
#    (200,200,None):(2.16,1.76),
#    (51,160,200  ):(2.56,1.93),
#    (76,140,160  ):(2.57,1.27),
#    (101,200,230 ):(1.08,0.87),
#    (151,200,200 ):(2.58,1.65),
#    (200,200,200 ):(1.75,0.22),
#    (121,200,None):(2.46,1.81),
#    (101,200,None):(1.54,1.46),
#    (121,200,720 ):(2.91,0.71),
#    (101,200,1040):(2.80,0.89)}


def random_cmt_2phase(D, d, C, L):
    best_sol = None
    best_f = None
    
    best_lambda = None
    best_mu = None
    
    for i in range(RANDOM_TRALS):
        lambda_multiplier=1.0+2.0*random()
        mu_multiplier=0.5+1.5*random()
        
        #pkey = (len(D), C, L)
        #if pkey in stochastic_params:
        #    lambda_multiplier, mu_multiplier = stochastic_params[pkey]
        #    
        #    # little wobble
        #    mu_multiplier+=random()*0.2-0.05
        #    lambda_multiplier+=random()*0.2-0.05
    
        # stochastic with SMAC
        # 8% gap {'pmu': 1.521411377875044, 'plambda': 3.920075418247227}
    
        sol = cmt_2phase_init(D, d, C, L, False,
                           lambda_multiplier, mu_multiplier,
                           phase1_seed_selection_method = "first",
                           phase2_choose_most_associated_route = False,
                           phase2_repeated_association_with_n_routes = 1,
                           number_of_randomized_retries=RANDOMIZED_RETRIES)
        
        sol_f = objf(sol, D)
        
        if (best_sol is None) or (sol_f<best_f):
            
            print("~~~~", "updating the best", "~~~~")
            best_sol = sol
            best_f = sol_f
            best_lambda = lambda_multiplier
            best_mu = mu_multiplier
    
    print("Best parameters for the cmt_2phase after %d trials"%RANDOM_TRALS,
          "lambda=%.2f"%best_lambda, "mu=%.2f"%best_mu)
    return best_sol

class CMT2PhaseReplications(ReplicationBase):
    
    def setUp(self):
        # Configured mu and lambda with SMAC
        
        # deterministic 
        # 1.1 {'pmu': 0.881673806062693, 'plambda': 3.022052284320399}        
        # all but CMT2
        # 0.222565292083 {'pmu': 1.1, 'plambda': 2.1}
        # CMT3,6,13 
        # 1.039674178 {'pmu': 1.3, 'plambda': 1.9}
        # CMT2 
        # 1.16688616136 {'pmu': 0.8, 'plambda': 2.7}        
        # CMT14
        #lambda_multiplier=2.0, mu_multiplier=0.9,
        
        self.algorithms = [
            ("cmt_2phase_stochastic", lambda pts,D,d,C,L,st:\
              random_cmt_2phase(D, d, C, L)),
            ("cmt_2phase_det_ps1", lambda pts,D,d,C,L,st:\
              cmt_2phase_init(D, d, C, L, minimize_K = False,
                                 lambda_multiplier=1.9, mu_multiplier=1.3,
                                 phase2_repeated_association_with_n_routes = 1)),
            ("cmt_2phase_det_ps2", lambda pts,D,d,C,L,st:\
              cmt_2phase_init(D, d, C, L, minimize_K = False,
                                 lambda_multiplier=2.7, mu_multiplier= 0.8,
                                 phase2_repeated_association_with_n_routes = 1)),
            ("cmt_2phase_det_ps3", lambda pts,D,d,C,L,st:\
              cmt_2phase_init(D, d, C, L, minimize_K = False,
                                 lambda_multiplier=2.0, mu_multiplier= 1.0,
                                 phase2_repeated_association_with_n_routes = 1))]

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
            "Christofides_11.vrp",
            "Christofides_12.vrp",
            "Christofides_13.vrp",
            "Christofides_14.vrp",]
        self.targets =  [
               [(5 , 547),
                (11, 883),
                (8 , 851),
                (12, 1093),
                (17, 1418),
                (6 , 565),
                (12, 969),
                (9 , 915),
                (15, 1245),
                (19, 1508),
                (7 , 1066),
                (10, 827),
                (11, 1612),
                (11, 876)],

        
#                # the deterministic are split among ps1, ps2, ps3
               [(5 , 547),
                None,
                None,
                None,
                None,
                (6 , 565),
                (12, 969),
                None,
                None,
                None,
                None,
                None,
                (11, 1612),
                None],
                
                
                [None,
                (11, 883),
                None,
                None,
                None,
                None,
                None,
                (9 , 915),
                None,
                None,
                None,
                None,
                None,
                (11, 876)],
                
                [None,
                None,
                (8 , 851),
                (12, 1093),
                (17, 1418),
                None,
                None,
                None,
                (15, 1245),
                (19, 1508),
                (7 , 1066),
                (10, 827),
                None,
                None]]
                        
        self.problem_path = path.join("Classic","CMT1979")
        
        
    # CMT2PhaseReplications.test_stochastic_CMT2Phase_with_CMT_instances
    def test_stochastic_CMT2Phase_with_CMT_instances(self):
        avgq, sdq, minq, maxq = self.solve_problems("cmt_2phase_stochastic",
                                                    cost_compare=False)
        self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT, "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT, "There is too much variation between instances") 

    # CMT2PhaseReplications.test_deterministic_CMT2Phase_with_CMT_instances
    def test_deterministic_CMT2Phase_with_CMT_instances(self):
        #TODO: currently this does unncecessary computation also on those 
        # instances we know the parameter set (ps1/ps2/ps3) does not work well,
        # change this so that the instance is solved with the right ps.
        avgq, sdq, minq, maxq = self.solve_problems("cmt_2phase_det_ps1",
                                                    cost_compare=False)
        avgq, sdq, minq, maxq = self.solve_problems("cmt_2phase_det_ps2",
                                                    cost_compare=False)
        avgq, sdq, minq, maxq = self.solve_problems("cmt_2phase_det_ps3",
                                                    cost_compare=False)
        
        
        #self.assertTrue( abs(avgq) < AVG_QUALITY_REQUREMENT*2, "Average quality not replicated (%.2f)"%avgq) 
        #self.assertTrue( abs(sdq) < SD_QUALITY_REQUREMENT*5, "There is too much variation between instances") 

if __name__ == '__main__':
    unittest.main()
