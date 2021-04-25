# -*- coding: utf-8 -*-

import numpy as np
from cvrp_io import _haversine 

points = [[24.42, 54.44], [24.44, 54.60], [24.42, 54.53]]

n_incl_depot = len(points)
geodesic_D = np.zeros( (n_incl_depot, n_incl_depot) )

for i, p1 in enumerate(points):
    for j, p2 in enumerate(points):
        geodesic_D[i,j] = 0.0 if p1==p2 else _haversine(p1,p2)
    
print(geodesic_D)
# The D can now be given to the VeRyPy *_init construction heuristic functions.