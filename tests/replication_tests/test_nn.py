# -*- coding: utf-8 -*-

from replicationbase import ReplicationBase, REPRO_QUALITY_LEVELS
from classic_heuristics.tyagi_nearest_neighbor import tyagi_init
from os import path
import unittest

class TestNearestNeighbourReplications(ReplicationBase):
    def setUp(self):
        self.algorithms = [
            ("tyagi", lambda pts, D,d,C,L,st:\
                tyagi_init(D,d,C,L))]

        self.problem_names =  ["02-eil13.vrp"]
        self.targets = [((4, 290),)]
        self.problem_path = path.join("Classic","ChristofidesEilon1969")
    def test_tyagi_with_DantzigRamser1959_n13_instance(self):
        avgq, sdq, minq, maxq = self.solve_problems("tyagi")
        self.assertTrue(abs(avgq) <= REPRO_QUALITY_LEVELS.A_AVG,
                        "Average quality not replicated (%.2f)"%avgq) 
        self.assertTrue(abs(sdq) <= REPRO_QUALITY_LEVELS.A_SD,
                        "There is too much variation between instances") 
   
if __name__ == '__main__':
    unittest.main()