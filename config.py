# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 15:43:09 2017

@author: juherask
"""

import os

DEBUG_VERBOSITY = 3

COST_EPSILON = 1e-10
CAPACITY_EPSILON = 1e-10

# how many seconds we give to a MIP solver
MAX_MIP_SOLVER_RUNTIME = 60*10 # 10m

MIP_SOLVER_THREADS = 1 # 0 is automatic (parallel computing)

# venv does not allow use of ~ for some reason in paths on Ubuntu 20.04.
BENCHMARKS_BASEPATH = os.path.join(os.environ["HOME"], r"Projects/Research/VRPBenchmarks")

LKH_EXE_PATH = os.path.join(os.environ["HOME"], r"Projects/Research/TSP/LKH-2.0.9/LKH")
LKH_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 

ACOTSP_EXE_PATH = os.path.join(os.environ["HOME"], r"Projects/Research/TSP/ACOTSP-master/acotsp")
ACOTSP_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 
