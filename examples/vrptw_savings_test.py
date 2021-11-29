# coding: utf-8

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

import numpy as np

import VeRyPy
import cvrp_io
import cvrp_ops
from util import objf
from classic_heuristics.vrptw_savings import vrptw_savings_init
from classic_heuristics.parallel_savings import parallel_savings_init
from shared_cli import print_solution_statistics

# Show all logging (very verbose)
import logging
logging.basicConfig(level=logging.DEBUG-2)

# Load a sample problem (manually converted Solomon 1987 instance)
N, points, _, d, D, C, _ = cvrp_io.read_TSPLIB_CVRP("examples/C101.25.vrp")
_, _, st, TWs = cvrp_io.read_TSBLIB_additional_constraints("examples/C101.25.vrp")

# model service time with the distance matrix
D_c = cvrp_ops.D2D_c(D, st) if st else D

ctrs = {
    'C': C,
    'TWs': TWs
}

sol = parallel_savings_init(D_c, d, ctrs)
print("The solution using non TW sensitive savings function was:")
print_solution_statistics(sol, D, D_c, d, C, None, TWs, st, verbosity=2 )
print()

#exit()

# Just make sure it works by using savings function with random values
best_random_sol = None
best_random_cost = float('inf')
for i in range(10):
    sol = vrptw_savings_init(D_c, d, ctrs)
    cost = objf(sol, D)
    if cost<best_random_cost:
        best_random_sol = sol
        best_random_cost = cost

print("The solution using TW specific savings function was:")
print_solution_statistics(best_random_sol, D, D_c, d, C, None, TWs, st, verbosity=2 )
print()

# Try with the soon (TM) to be implemented actual TW specific savings function
sol = vrptw_savings_init(D_c, d, ctrs, debug_with_random_savings=False)
print("The solution using TW specific savings function was:")
print_solution_statistics(best_random_sol, D, D_c, d, C, None, TWs, st, verbosity=2 )
print()