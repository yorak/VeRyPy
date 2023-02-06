# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from os import path
from time import time

import numpy as np

from verypy.util import objf, without_empty_routes
from verypy.shared_cli import print_solution_statistics
from verypy.cvrp_ops import validate_solution_feasibility
import verypy.cvrp_io as cvrp_io


from verypy.config import BENCHMARKS_BASEPATH

class REPRO_QUALITY_LEVELS:
    # Allow more leeway as there are many inaccuracies in the implementation details
    A_AVG = 0.2 # 0.2%
    A_SD = 1.0 # 1.0%

    # general level is good, however large variation from instance to instance
    B_AVG = 0.5 # 0.5%
    B_SD = 2.5 # 2.5%
    
    # general level is not so good, and/or very large variation, many tricks needed
    C_AVG = 1.0 # 1.0%
    C_SD = 3.5 # 3.5%
    
    # not reproed and/or many tricks needed
    D_AVG = 2.0 # 2.0%
    D_SD = 5.0 # 5.0%


class ReplicationBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # log every result to a file
        import __main__ as main
        rfn = main.__file__.replace(".pyo", "_results.txt")
        rfn = main.__file__.replace(".pyc", "_results.txt")
        rfn = main.__file__.replace(".py", "_results.txt")
        rfn = path.basename(rfn)
        
        ReplicationBase.result_file = None
        raw_result_path = path.join("replication_results", "raw")
        if path.isdir(raw_result_path):
            rfn = path.join(raw_result_path, rfn)
            nbr = 0
            unique_rfn = rfn
            while path.isfile(unique_rfn):
                unique_rfn = rfn.replace("_results.txt", "_results%d.txt"%nbr)
                nbr+=1
            ReplicationBase.result_file = open(unique_rfn, "w")
        else:
            print("WARNING: the directory \"%s\", which is the target for the replication results, does not exist."%raw_result_path)
    
    @classmethod
    def tearDownClass(cls):
        if (ReplicationBase.result_file):
            ReplicationBase.result_file.close()
        
    def setUp(self):
        """ This is meant to be overridden in derived unittsest class """
        
        self.algorithms = [
            # example algorithm for illustrative purposes, returns a solution 
            # with N feasible routes, each serving a single customer 
            ("algo_name", lambda D,d,C,L:
              [0]+[n for p in zip(range(1,len(D)), [0]*len(D)) for n in p])
            # ...
            ]
    
        self.problem_names = [
            "problem1.vrp",
            "problem2.vrp",
            #..."
            ]

        self.problem_path = path.join("Folder1","Folder2")
        
        # target quality for each algorithm and problem
        self.targets = [
                
            # target solution quality and number of veicles for algo_name on 
            # problem1, problem2, ..
            ((4, 12.3), (6, 45.6), )
            #...
            ]

    def _solve_instance(self, algo,pfn, round_D_func = None, require_K = False,
                        predefined_k = None, suppress_constraint_check = False):
        N, points, dd_points, d, D, C, _ = cvrp_io.read_TSPLIB_CVRP(pfn)
        K, L, service_time = cvrp_io.read_TSBLIB_additional_constraints(pfn)
        
        if round_D_func:
            D = round_D_func(D)
        
        if predefined_k is not None:
            K = predefined_k
        if require_K and (K is None):
            raise IOError("It is required that the VEHICLE field is set in %s"%pfn)
            
        if service_time:
            half_st=service_time/2.0
            if int(half_st)==half_st:
                half_st = int(half_st)
                service_time = int(service_time)
                
            # The service time can be modeled modifying the distance
            #  matrix in a way that any visit to a depot node costs
            #  service_time units.
            D_c = np.copy(D)
            D_c[1:,1:]+=service_time
            D_c[0,:]+=half_st
            D_c[:,0]+=half_st
            np.fill_diagonal(D_c,0.0)
        else:
            D_c = D
            
        if points is None and dd_points is not None:
            points = dd_points
        
        startt = time()
        if require_K:
            sol = algo(points,D_c,d,C,L,service_time,K)
        else:
            sol = algo(points,D_c,d,C,L,service_time)
        endt = time()
        elapsedt = endt - startt
        
        if __debug__:
            print_solution_statistics(sol, D, D_c, d, C, L,service_time)
        
        
        cover_ok, capa_ok, rlen_ok = validate_solution_feasibility(sol, D,d,C,L,True)
        if not suppress_constraint_check:
            self.assertTrue( cover_ok, "Must be a valid solution")
            self.assertTrue( capa_ok, "Must not violate the C constraint" )
            self.assertTrue( rlen_ok, "Must not violate the L constraint"  )
                
        return sol, objf(sol, D), objf(sol, D_c), elapsedt
        
    def solve_problems(self, algo_name, instance_idx='all',
                       round_D_func = None,
                       round_f_func = None,
                       cost_compare = True,
                       require_K = False,
                       suppress_constraint_check=False):
        """
        Solves the problems set up in the setUp using the algo_name. The other
        arguments are:
            
        instance_idx: if 'all' (default), solve all problems. If int, solve the
         problem with that index. If a iterable of ints, solve those problems.
            
        round_D_func: nxn matrix -> nxn matrix, operation that can be used to
         modify the costs of the distance matrix. Usually e.g. np.int_,
         np.around etc.
        
        round_f_func: float -> int/float, the function used to round the result
         (i.e. solution quality).
        
        cost_compare: compare the solution against the target in self.target by
         calculating the quality with the distance or cost matrix (cost matrix
         includes service times).
         
        require_K: gives the required number of vehicles K as a parameter for
         the algorithm.
         
        suppress_constraint_check: can be uesd to disable constraint checking.
        """
        
        algo_idx = (i for i,al in enumerate(self.algorithms)\
                      if al[0]==algo_name).next()
        assert self.algorithms[algo_idx][0]==algo_name
        
        algo_targets = self.targets[algo_idx]
        
        diffs = np.zeros(len(algo_targets))
        
        active_targets = list(enumerate(algo_targets))
        if instance_idx!='all':
            if hasattr(type(instance_idx), '__iter__'):
                for i in instance_idx:
                    active_targets = [active_targets[i]]
            elif type(instance_idx) is int:
                active_targets = [active_targets[instance_idx]]
            else:
                raise ValueError("Invalid problem index value")
            
        for problem_idx, target in active_targets:
            if target is None:
                continue
            elif (type(target) is int) or (len(target)==1):
                target_c = target
                target_k = None
            else:
                target_k, target_c = target
                
            problem_name = self.problem_names[problem_idx]
            pfn = path.join(BENCHMARKS_BASEPATH,self.problem_path, problem_name)
            
            #print(pfn)
            
            sol,sol_f,sol_c,elapsed_t = self._solve_instance(
                self.algorithms[algo_idx][1],pfn,
                round_D_func=round_D_func,
                require_K=require_K,
                predefined_k = target_k,
                suppress_constraint_check=suppress_constraint_check)
            
            sol = without_empty_routes(sol)
            sol_k = sol.count(0)-1
            
            if round_f_func:
                sol_f = round_f_func(sol_f+0.5)
                sol_c = round_f_func(sol_c+0.5)
            
            if cost_compare:
                gap = 100.0*sol_c/target_c-100.0
                print(algo_name, problem_name, sol_c, "VS.", target_c,
                      "(gap %.2f %%)"%gap, "in %.2f s"%elapsed_t)
            else:
                gap = 100.0*sol_f/target_c-100.0
                print(algo_name, problem_name, sol_f, "VS.", target_c,
                      "(gap %.2f %%)"%gap, "in %.2f s"%elapsed_t)
            
            if ReplicationBase.result_file:
                ReplicationBase.result_file.write(";".join( [
                        algo_name,
                        problem_name,
                        str(sol_c) if cost_compare else str(sol_f),
                        str(sol_k),
                        str(target_c),
                        str(target_k),
                        str(gap),
                        str(elapsed_t)]))
                ReplicationBase.result_file.write("\n")
                ReplicationBase.result_file.flush()
                
            diffs[problem_idx] = gap
            
        avg_diff = np.average(diffs)
        sd_diff = np.std(diffs)
        print("On average",algo_name,"is within",
              "%.2f %% (SD %.2f %%) of the target.\n"%(avg_diff, sd_diff), "\n")
              
        return avg_diff, sd_diff, np.min(diffs), np.max(diffs)
    
            
if __name__ == '__main__':
    unittest.main()     
            
