#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and demonstrates the simple use case of solving a single
TSPLIB formatted problem instance file with a single heuristic algorithm and
printing the resulting solution route by route."""
###############################################################################

import verypy.cvrp_io as cvrp_io
from verypy.classic_heuristics.parallel_savings import parallel_savings_init
from verypy.util import sol2routes

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
E_n51_k5_path = os.path.join(dir_path, "E-n51-k5.vrp")

problem = cvrp_io.read_TSPLIB_CVRP(E_n51_k5_path)

solution = parallel_savings_init(
    D=problem.distance_matrix, 
    d=problem.customer_demands, 
    C=problem.capacity_constraint)

for route_idx, route in enumerate(sol2routes(solution)):
    print("Route #%d : %s"%(route_idx+1, route))
