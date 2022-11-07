#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This is an example on how to use all algorithms from VeRyPy to solve a 
problem given in a Python dict from a separate file (not included).
The results are printed to stdout in a TSV (tab separated file) format.
"""
###############################################################################

from __future__ import print_function
from os import path

HOME_PATH = path.expanduser("~")
#TODO: Set your VeRyPy path here
YOUR_VERYPY_PATH = path.join(HOME_PATH, r'Research/VeRyPy')

# TODO: Define your problem in data.py. It is a dictionary with following keys
# See sample_data_dict.py for an example.
"""
Data = {
'distance_matrix' : <2D array given as a list of (int or double) lists.>
'demands' : <List of all demands. At index 0 is the depot with demand of 0.>
'vehicle_capacities' : <List of (equal) capacities for each vehicle.>
'lats' : <Latitudes of the coordinates of the depot (index 0) and customers.>
'longs' : <Longitutdes of the coordinates of the depot (index 0) and customers.>

"""
from sample_data_dict import Data

import sys
sys.path.append(YOUR_VERYPY_PATH)

try:
    from time import clock
except:
    from time import process_time as clock

import numpy as np
from verypy import get_algorithms
from verypy.cvrp_ops import normalize_solution, calculate_objective

def check_symmetric(A, tol=1e-8):
    """ Helper function to check if a matrix is symmetric. """
    return np.all(np.abs(A-A.T) < tol)

### 1. Convert data to VeRyPy format ###
C = Data['vehicle_capacities'][0]
assert \
    all( vC==C for vC in Data['vehicle_capacities']),\
    "VeRyPy requies all vehicles to have the same capacity"

points = list(zip(Data['lats'], Data['longs']))
D = np.array( Data['distance_matrix'] )
d = Data['demands']

assert \
    len(points) == len(d) == D.shape[0] == D.shape[1],\
    "Coordinates, demands, and the distance matrix have to have same number of elements"
    
assert \
    check_symmetric(D),\
    "VeRyPy only supports symmetric distances"

### 2. Solve the problem using all algorithms from VeRyPy ###
algos = get_algorithms(['all'])
print("algorithm", "objf", "K", "sol", "t", sep='\t')
for algo_abbreviation, algo_name, _, algo_f in algos:
    elapsed_t = 0
    try:
        start_t = clock()
        sol = algo_f(points, D, d, C, None, None, "GEO", False, True)
        elapsed_t = clock()-start_t
    except:
        start_t = clock()
        sol = algo_f(points, D, d, C, None, None, "GEO", False, False)
        elapsed_t = clock()-start_t
    
    sol = normalize_solution(sol)
    obj = calculate_objective(sol, D)
    K = sol.count(0)-1
    print(algo_name, obj, K, sol, elapsed_t, sep='\t')