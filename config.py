# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 15:43:09 2017

@author: juherask
"""

DEBUG_VERBOSITY = 3

COST_EPSILON = 1e-10
CAPACITY_EPSILON = 1e-10

# how many seconds we give to a MIP solver
MAX_MIP_SOLVER_RUNTIME = 60*10 # 10m

MIP_SOLVER_THREADS = 1 # 0 is automatic (parallel computing)

BENCHMARKS_BASEPATH = r"C:\Users\juherask\Dissertation\Phases\Benchmarks"

LKH_EXE_PATH = r"C:\Users\juherask\Dissertation\Phases\Selection\solvers\lkh\lkh.exe"
LKH_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 

ACOTSP_EXE_PATH = r"C:\Users\juherask\Dissertation\Phases\Selection\solvers\ACOTSP\acotsp.exe"
ACOTSP_EXACT_DISTANCES_PRECISION_DECIMALS = 1000.0 # of the form 0.123 
