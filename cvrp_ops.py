# -*- coding: utf-8 -*-
################################################################################
""" This file has some convinient operations on normalizing CVRP solutions,
checking their feasibility, and calculating the solution quality.
"""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

from sys import stderr
from itertools import groupby
from config import COST_EPSILON as S_EPS
from config import CAPACITY_EPSILON as C_EPS

import numpy as np

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

################################################################################


def _is_all_integer_array(D):
    return np.all(np.equal(np.mod(D, 1), 0)) # pylint: disable=no-member

def _list_trim(l, e):
    """ works like string trimming but for lists """
    trimmed = list(l)
    while len(trimmed)>0 and trimmed[0]==e:
        del trimmed[0]
    while len(trimmed)>0 and trimmed[-1]==e:
        del trimmed[-1]
    return trimmed
    
def _list_split(l, e):
    """ works like string splitting but for lists """
    return [list(group) for k, group in groupby(l, lambda x:x==e) if not k]  
    
def normalize_solution(sol):
    """ Make sure solution routes are ordered as follows
    1) route start pt < route end pt
    2) prev route start < next route start
    Also removes repeated depots, e.g., [0,0,1,2,0,0,3,0] -> [0,1,2,3,0]
    """
    
    # Already in parts
    if hasattr(sol[0], '__len__'):
        routes = [_list_trim(route, 0) for route in sol]
    else:
        #print(sol
        routes = _list_split(sol, 0)
    nroutes = []
    # 1
    for route in routes:
        if route[-1]<route[0]:
            nroutes.append( list( reversed(route) ) )
        else:
            nroutes.append( list(route) )
    
    # 2
    nroutes.sort()
    
    nsol = [0]
    for nroute in nroutes:
        nsol += nroute + [0]
    
    return nsol

def generate_missing_coordinates(for_D):
    from sklearn import manifold
    mds = manifold.MDS(n_components=2, dissimilarity='precomputed',
                       random_state=42)
    mds_results = mds.fit(for_D)
    points = list( mds_results.embedding_ )
    edge_weight_type = "EUC_2D" if _is_all_integer_array(for_D) else "EXACT_2D"
    return points, edge_weight_type    
    
def check_solution_feasibility(solution, D, d=None,
                               C=None, L=None, print_violations=False):
    """ This checks if the solution is feasible. This feasibility checker also
     supports solutions that do not have 0 as their first node as long as there
     is at least one visit to the depot (0).
    
    Note: This is not performance optimized in any way. Therefore, it is 
     advisable to use this for solution verification purposes only, and not i.e.
     as a building block of an algorithm.
     
    Returns the feasibility status as triple-tuple:
    (covering_feasibility, capacity_feasibility, route_cost_feasibility)
    """
    
    if D is not None:
        N = len(D)
    else:
        N = len(d)
        
    start_node = -1
    tail_c = 0
    tail_l = 0
    covering = [0]*N
    
    covering_feasibility = True
    capacity_feasibility = True
    route_cost_feasibility = True
    
    # some algorighms do not have the 0 as the first node, but start
    #  mid-route point (that is, solution[-1]->solution[0] is an edge)
    prev_node = solution[-1]
    for i,node in enumerate(solution):
        if covering[node]:
            if node!=0:
                if print_violations:
                    print("CONSTRAINT VIOLATION: customer n%d is served twice"%node, file=stderr)
                covering_feasibility = False
        else:
            covering[node]=1
        
        if C:
            tail_c += d[node]
        if L:
            tail_l += D[prev_node, node]
            prev_node = node
        # Depot
        if (node==0):
            start_node = i
            break
    
    c = 0.0
    l = 0.0
    # check second half from first visit to depot to end
    for node in solution[start_node+1:]:
        if covering[node]:
            if node!=0:
                if print_violations:
                    print("CONSTRAINT VIOLATION: customer n%d is served twice"%node, file=stderr)
                covering_feasibility = False
        else:
            covering[node]=1

        if C:
            if node!=0:
                c += d[node]            
            elif c-C_EPS>C:
                if print_violations:
                    print("CONSTRAINT VIOLATION: capacity is exceeded by %.2f"%(c-C), file=stderr)
                capacity_feasibility = False
            
        if L:
            l += D[prev_node, node]
            prev_node = node
            if node==0 and l-S_EPS>L:
                if print_violations:
                    print("CONSTRAINT VIOLATION: maximum route cost is exceeded by %.2f"%(l-L), file=stderr)
                route_cost_feasibility =  False

        if node==0:
            c = 0.0
            l = 0.0
                   
    if sum(covering)!=N:
        if print_violations:
            print("CONSTRAINT VIOLATION: some of the nodes are not served", file=stderr)
        covering_feasibility = False
        
    # check the possible remaining first half from beginning to first depot
    c += tail_c
    if C and c-C_EPS>C:
        if print_violations:
            print("CONSTRAINT VIOLATION: capacity is exceeded by %.2f"%(c-C), file=stderr)
        capacity_feasibility = False
    l += tail_l
    if L and l-S_EPS>L:
        if print_violations:
            print("CONSTRAINT VIOLATION: maximum route cost is exceeded by %.2f"%(l-L), file=stderr)
        route_cost_feasibility = False
    
    return (covering_feasibility, capacity_feasibility, route_cost_feasibility)

def check_route_feasibility(routes, D=None, d=None, C=None, L=None):
    
    N = len(d) if d else len(D) 
    served = [False]*N
    served[0] = True
    
    for route in routes:
        route_demand = 0
        route_cost = 0
        p = 0
        for n in route:
            if n!=0 and served[n]:
                print("ERROR: 'serve only once' constraint violated", file=stderr)
            served[n]=True
            if C:
                route_demand+=d[n]
            if L:
                route_cost+=D[p,n]
                p = n
        if L:
            route_cost+=D[p,0]
            if route_cost-S_EPS > L:
                print("ERROR: maximum route cost violated", file=stderr)
                return False
        if C and route_demand-C_EPS > C:
            print("ERROR: route capacity violated", file=stderr)
            return False
        
    if sum(served)!=len(d):
        print("ERROR: all vertices are not served", file=stderr)
        
    return True
   
def fast_constraint_check(sol,D,d,C,L):
    """ A faster version of the constrain checker. E.g. does not check 
    if all customers are served (it is assumed this constraint is not violated.
    """
    prev_node = None
    c = 0.0
    l = 0.0
    for node in sol:    
        if C and node!=0:
            c += d[node]
            if c-C_EPS>C:
                return False
        if L and prev_node!=None:
            l += D[prev_node,node]
            if l-S_EPS>L:
                return False
        if(node==0):
            c = 0.0
            l = 0.0
        prev_node  = node
    return True

def calculate_objective(sol, D):
    """ The objective is the total cost of all routes in VRP solution. D is 
    the full distance matrix of the points in the problem (also includes the
    depot with index of 0), and solution is a list containing giant tour
    encoded VRP solution, where 0 indicates a visit to the depot OR a list of 
    routes leaving and returning to the depot.
    """
    f = 0.0
    if hasattr(sol[0], '__iter__'):
        # routes are separately
        for route in sol:
            f+= sum( D[route[i-1],route[i]] for i in range(1,len(route)))
            # sometimes a end or start (or both) visit to depot may be missing
            if route[0]!=0:
                f+=D[0,route[0]]
            if route[-1]!=0:
                f+=D[route[-1],0]
    else:
        # giant tour encoding
        f = sum( D[sol[i-1],sol[i]] for i in range(1,len(sol)))
        f+=D[sol[-1], sol[0]]
    return f

def D2D_c(D, st):
    """ The service time can be modeled modifying the distance matrix in a way
    that any visit to a depot node costs service_time units. The D may be int
    or double.  Convert it if modeling requires it.
    
    Args:
        D (np.ndarray) is the NxN distance matrix with integer or real dist.
        st (float or int) is the service time
        
    Return:
        returns a copy of the distance matrix with the service times baked in.
    """
    
    if not st:
        raise ValueError("Service time not set")
    
    tst = D.dtype.type(st)
    if float(tst)!=float(st) or float(tst/2)!=float(st/2.0):
        D_c = D.astype('d')
        halftst = tst/2.0
    else:
        D_c = np.copy(D)    
        halftst = int(tst/2)
    
    D_c[1:,1:]+=tst
    D_c[0,:]+=halftst
    D_c[:,0]+=halftst
    np.fill_diagonal(D_c,0.0)
    
    return D_c
