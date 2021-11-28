# coding: utf-8
import VeRyPy
from classic_heuristics.vrptw_savings import vrptw_savings_init

# Show all logging (very verbose)
import logging
logging.basicConfig(level=logging.DEBUG-2)

# Init a super simple 1 customer problem
import numpy as np
D = np.array( [[0,1], [1,0]] )
d = [0, 1]

sol = vrptw_savings_init(D, d, {'C':2, 'TWs':[(0,0), (0,1)]})
print("The solution for our simple test was", sol)

