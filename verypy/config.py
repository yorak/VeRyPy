# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 15:43:09 2017

@author: juherask
"""

from os import path, name
exe_ext = ".exe" if name == 'nt' else ""

DEBUG_VERBOSITY = 3

COST_EPSILON = 1e-10
CAPACITY_EPSILON = 1e-10

# how many seconds we give to a MIP solver
MAX_MIP_SOLVER_RUNTIME = 60*10 # 10m

MIP_SOLVER_THREADS = 1 # 0 is automatic (parallel computing)

# Set up some paths where to find benchmarks and external solvers
HOME_PATH = path.expanduser("~")

BENCHMARKS_BASEPATH = path.join(HOME_PATH, r"Research/VRPBenchmarks")

LKH_EXE_PATH = path.join(HOME_PATH, r"Research/TSP/LKH-2.0.9/lkh"+exe_ext)
LKH_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 

ACOTSP_EXE_PATH = path.join(HOME_PATH, r"Research/TSP/ACOTSP-1.03-ls/acotsp"+exe_ext)
ACOTSP_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 
