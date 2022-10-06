# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from os import path

from verypy.classic_heuristics.mole_jameson_insertion import mole_jameson_insertion_init
from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS

# Mole & Jameson (1976) insertion has five different approaches, each tested
#  separately. However, replicating is hard as the ...
#  1) ... secondary strain sorting criterion is not given and there seems to be 
#          a lot of ties!
#  2) ... the strategy or other details for the local search are not given.

class TestMoleJamesonReplications(ReplicationBase):
    # static storage of results
    results = {}
    
    def setUp(self):  
        self.algorithms = [
            ('proximity_ranking', lambda pts, D,d,C,L,st:\
                   mole_jameson_insertion_init(D,d,C,L,
                       minimize_K = True,
                       strain_criterion='proximity_ranking')
            ),
            ('min_strain', lambda pts, D,d,C,L,st:\
                 mole_jameson_insertion_init(D,d,C,L,
                       minimize_K = True,
                       strain_criterion='min_strain')
            ),
            ('clarke_wright', lambda pts, D,d,C,L,st:\
                 mole_jameson_insertion_init(D,d,C,L,
                       minimize_K = True,
                       strain_criterion='clarke_wright')
            ),
            ('gaskell', lambda pts, D,d,C,L,st:\
                mole_jameson_insertion_init(D,d,C,L,
                        minimize_K = True,
                       strain_criterion='gaskell')
            ),
            ('augumented_min_strain', lambda pts, D,d,C,L,st:\
                mole_jameson_insertion_init(D,d,C,L,
                       minimize_K = True,
                       strain_criterion='augumented_min_strain')
            )]
            
        self.problem_names =  [
            "01-eil7.vrp",
            "02-eil13.vrp",
            "03-Gaskell1-E022-k4g_euc2d.vrp",
            "04-Gaskell2-E023-k5g_euc2d.vrp",
            "05-Gaskell3-E030-k4g_euc2d.vrp",
            "06-CW64_n32_k9c.vrp",
            "07-Gaskell4-E033-k4g_euc2d.vrp",
            "08-eil51.vrp",
            "09-eil76.vrp",
            "10-eil101.vrp"]

        # CW64_n31_k8c is without Birmingham with demand of 140
        #  (it is always served separately with a route of length 82*2)
        self.targets = [
            #with proximity (ranking) criterion
            ((2,119),
            (4,290),
            (4,694),
            (5,994),
            (5,937),
            (8,1422),
            (5,836),
            (5,672),
            (10,1065),
            (8,948)),
            # min.strain
            ((2,121),
            (4,294),
            (4,633),
            (5,994),
            (5,951),
            (8,1465),
            (5,856),
            (5,603),
            (10,1017),
            (8,882)),
            # with savings (ClarkeWright) criterion
            ((2,119),
            (4,290),
            (4,642),
            (5,979),
            (4,923),
            (8,1422),
            (4,834),
            (5,597),
            (10,934),
            (8,955)),
            # with parametrized savings criterion (best of 4)
            ((2,119),
            (4,290),
            (4,598),
            (5,949),
            (4,890),
            (8,1422),
            (4,812),
            (5,590),
            (10,910),
            (8,926)),
             # with parametrized min strain criterion (best of 5)
            ((2,119),
            (4,290),
            (4,594),
            (5,949),
            (4,891),
            (8,1422),
            (4,834),
            (5,575),
            (10,919),
            (8,882)),
            ]
            
        self.problem_path = path.join("Classic","ChristofidesEilon1969")

    
    def test_proximity_criteria_with_CE1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems('proximity_ranking', round_f_func = _to_int)
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        TestMoleJamesonReplications.results['proximity_ranking'] = (avgq, sdq, minq, maxq)

    def test_minstrain_criteria_with_CE1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems('min_strain', round_f_func = _to_int)
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        TestMoleJamesonReplications.results['min_strain'] = (avgq, sdq, minq, maxq)
        
    def test_savings_criteria_with_CE1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems('clarke_wright', round_f_func = _to_int)
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        TestMoleJamesonReplications.results['clarke_wright'] = (avgq, sdq, minq, maxq)

    def test_parametrized_savings_criteria_with_CE1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems('gaskell', round_f_func = _to_int)
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        self.results['gaskell'] = (avgq, sdq, minq, maxq)
        
    def test_parametrized_minstrain_criteria_with_CE1969_instances(self):
        _to_int = lambda x: int(x)
        avgq, sdq, minq, maxq = self.solve_problems('augumented_min_strain', round_f_func = _to_int)
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        TestMoleJamesonReplications.results['augumented_min_strain'] = (avgq, sdq, minq, maxq)
    
    def tearDown(self):
        # There are no subtests in python 2.7. Abuse tearDown and class variable
        #  to check the overall quality at the end.
        target_algo_names,_ = zip(*self.algorithms)
        target_algo_names = set(target_algo_names)
        result_algo_names = set(TestMoleJamesonReplications.results.keys())
        if target_algo_names!=result_algo_names:
            return
            
        avg_avgq=0.0
        avg_sdq=0.0
        for (avgq, sdq, minq, maxq) in TestMoleJamesonReplications.results.values():
            avg_avgq+=avgq
            avg_sdq+=sdq
        avg_avgq = avg_avgq/len(self.results.values())
        avg_sdq = avg_sdq/len(self.results.values())
        
        print("\n In aggregate they are within %.2f %% (SD %.2f %%) of the target."%
              (avg_avgq,avg_sdq))
        self.assertTrue(abs(avgq) < REPRO_QUALITY_LEVELS.D_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) < REPRO_QUALITY_LEVELS.D_SD,
                        "There is too much variation between instances") 
        
if __name__ == '__main__':
    #import logging
    #logging.basicConfig(level=logging.DEBUG-5, filename='debug.log')
    unittest.main()
            