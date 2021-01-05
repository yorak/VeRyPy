#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
""" This file is a part of the VeRyPy classical vehicle routing problem
heuristic library and provides an implementation of the Fisher&Jaikumar (1981)
heuristic, which generates an approximate solution for a VRP via solving it as
an generalized assignment problem (GAP).

The script is callable and can be used as a standalone solver for TSPLIB 
formatted CVRPs. It has extensive dependencies: MIP solver Gurobi, built-in TSP
solver, and numpy and scipy for reading and preparing the problem instance."""
###############################################################################

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from signal import signal, SIGINT, default_int_handler
from collections import namedtuple
from math import pi, ceil
from logging import log, DEBUG, WARNING
from sys import stderr

import numpy as np
from gurobipy import Model, GRB, LinExpr, GurobiError

try:
    ## For reasonably sized instances you might want to get the optimal TSP solution
    ## with Gurobi. This is also the default used in Rasku et al. (2019) experiments.
    from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi as solve_tsp
    ## For larger instances you might want to use ether of the faster TSP solvers.
    #from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh as solve_tsp
    #from tsp_solvers.tsp_solver_acotsp import solve_tsp_acotsp as solve_tsp
except ImportError:
    print("WARNING: could not use the external TSP solver (probably the executable is not found). "+
          "Relying on internal TSP solver and the results may differ from those that were published.", file=stderr)
    from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt as solve_tsp


from classic_heuristics.sweep import get_sweep_from_cartesian_coordinates, bisect_angle
from cvrp_io import calculate_D
from util import is_better_sol, totald
from config import MAX_MIP_SOLVER_RUNTIME, MIP_SOLVER_THREADS
from config import CAPACITY_EPSILON as C_EPS
from config import COST_EPSILON as S_EPS

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"


# These hard coded parameters define how the relaxation is adjusted if the 
#  GAP solution is not L feasible.
L_MPLR_DEFAULT = 1.0
L_ADAPTIVE_MPLR_INIT = 0.85 
L_ADAPTIVE_MPLR_INC = 0.85
L_ADAPTIVE_MPLR_MAX_TRIES = 3
INCREASE_K_ON_FAILURE_UPTO = 1.1 # = 10% increase to K (or min of 1)

def _decision_variables_to_assignments(m, Y_ik, N, K):
    """ Convert the decision variables in m for keys Y_ik to assignments of 
    customers i==2..N (0 is the depot) to the routes k=1..K . 
    
    TODO: there is probably a neat numpy trick to get node and k indices
      out of the decision variable array. For now just use nested loops,
      but a cleverer way would problably be faster.
     However, "premature optimization is the root of all evil", so profile
      first, and modify only after verifying it to be a real bottleneck."""
      
    assignments = []
    Y_ik_values = m.getAttr('x', Y_ik)
    for k in range(K):
        route_nodes = []
        for i in range(1, N):
            if Y_ik_values[i,k]:
                route_nodes.append(i)
        assignments.append(route_nodes)
    return assignments
                
def _solve_gap(N, D_s, d, C, K, L=None, L_ctr_multipiler=1.0):
    """A helper function that Solves VRP as a Generalized Assignment Problem
    to assign customers to vehicles with a objective function that the delivery
    cost as described in (Fisher & Jaikumar 1981).
    
    D_s is the distance matrix complemented with distances to K seed points.
    That is:
          [D_0]
    D_s = [S  ], where D_0 is the first row of the full distance matrix and
                         S is the distances from seed points to node points
    d is the list of customer demands with d[0]=0 being the depot node
    C is the capacity of the K identical trucks
    
    also, additional (and optional) constraints can be given:
    
    L is the maximum tour cost/duration/length
    L_ctr_multipiler allows iteratively adjusting the max route cost
     approximation constraint in order to avoid producing assignments that are
     ruled infeasible by the feasibility checker. 
    --
    Fisher, M. L. and Jaikumar, R. (1981), A generalized assignment 
    heuristic for vehicle routing. Networks, 11: 109-124.
    """
    
    ## build the cost approximation matrix "insertion_cost"
    
    # it is ~ a insertion cost matrix, where each coefficient is the cost 
    #  of inserting customer i to the route consisting visit to seed k.
    #
    # we assume that distances are symmetric, but if asymmetric
    #  distances are to be used, take min
    # d_{ik} = min(c_{0i}+c_{i{i_k}}+c_{i{i_k}},
    #              c_{0{i_k}}+c_[{i_k}i}+c_{i0})
    #              -(c_{0{i_k}}+c_{{i_k}0})
     
    m = Model("GAPCVRP")
        
    # the order of the keys is important when we interpret the results
    Y_ik_keys = [(i,k) for k in range(K) for i in range(1,N)]
    
    # delivery cost approximation coefficients for the objective function
    insertion_cost = {(i,k): D_s[0,i]+D_s[k,i]-D_s[k,0] \
                              for i,k in Y_ik_keys}
                                  
    # variables and the objective
    Y_ik = m.addVars(Y_ik_keys, obj=insertion_cost, vtype=GRB.BINARY, name='y')    
    
    ## constraints

    # c1, the capacity constraint and optional tour cost constraint cl
    approx_route_cost_constraints = []    
    if C: c1_coeffs = d[1:]
    for k in range(K):
        ck_vars = [Y_ik[i,k] for i in range(1,N)] 
        if C:
            c1_lhs = LinExpr(c1_coeffs,ck_vars)
            #c1_lhs = Y_ik.prod(c1_coeffs, '*', k)
            m.addConstr(c1_lhs <= C, "c1_k%d"%k) 
            
        # ct = optional tour cost constraints
        #  it is a bit hidden, but the additional side constraint can be found
        #  from Fisher & Jaikumar (1981) p121, 2. paragraph.
        # However, for whatever reason, this does not seem to produce the 
        #  same results as reported in their paper as the constraint easily
        #  starts to make the problem infeasible and the exact mechanism to 
        #  recover that is not specified in the paper.
        if L:
            ct_coeffs = [insertion_cost[(i,k)]*L_ctr_multipiler for i in range(1,N)]
            ct_lhs = LinExpr(ct_coeffs,ck_vars)
            #ct_lhs = Y_ik.prod(ct_coeffs, '*', k)
            constr_l = m.addConstr(ct_lhs <= L, "cl_k%d"%k)
            approx_route_cost_constraints.append(constr_l)

    # c2, the assignment constraints 
    for i in range(1,N):
        # c2_1..N every node assigned only to 1 route
        m.addConstr(Y_ik.sum(i, '*') == 1, "c1_i%d"%i) 
        
    ## update the model and solve 
    m._vars = Y_ik  
    m.modelSense = GRB.MINIMIZE
    m.update()
    #m.write("gapvrp_model.lp")
    # disable output
    m.setParam('OutputFlag', 0)    
    m.setParam('Threads', MIP_SOLVER_THREADS)
    # REMOVEME
    m.setParam('MIPFocus', 3)
    m.setParam('TimeLimit', MAX_MIP_SOLVER_RUNTIME)

    m.optimize()

    # restore SIGINT callback handler which is changed by gurobipy
    signal(SIGINT, default_int_handler)
    
    if __debug__:
        log(DEBUG-1,"Gurobi runtime = %.2f"%m.Runtime)
    
    if m.Status == GRB.OPTIMAL:
        return _decision_variables_to_assignments(m, Y_ik, N, K)
    elif m.Status == GRB.INFEASIBLE and L:
        # relax the model and allow violating minimal number of the approximate 
        #  route length constraints
        pens = [1.0]*len(approx_route_cost_constraints)
        m.feasRelax(1, True, None, None, None, approx_route_cost_constraints, pens)
        # TODO: not sure if feasRelax can change Status, test it someday
        if m.Status == GRB.INTERRUPTED:
            raise KeyboardInterrupt() # pass it on
        m.optimize()

        # restore SIGINT callback handler which is changed by gurobipy
        signal(SIGINT, default_int_handler)

        status = m.Status       
        if __debug__:
            log(DEBUG-1, "Relaxed problem Gurobi runtime = %.2f"%m.Runtime)
        if status == GRB.OPTIMAL:
            return _decision_variables_to_assignments(m, Y_ik, N, K)
        elif status == GRB.TIME_LIMIT:
            raise GurobiError(10023, "Gurobi timeout reached when attempting to solve relaxed SCPCVRP")
        elif m.Status == GRB.INTERRUPTED:
            raise KeyboardInterrupt() # pass it on
        return None
    elif m.Status == GRB.TIME_LIMIT:
        raise GurobiError(10023, "Gurobi timeout reached when attempting to solve GAP")
    elif m.Status == GRB.INTERRUPTED:
        raise KeyboardInterrupt() # pass it on
    return None
    
_Cone = namedtuple('_Cone', ['phi1', 'phi2', 'demand', 'nodes'])
def _sweep_seed_points(points, D, d, C, K, trial=0):
    """A seed point generation function that implements the rule used in
    Fisher and Jaikumar (1981) to select the seed customers for the delivery
    cost approximation calculation in their VRP heuristic. It is assumed that
    all customers are located on a plane with euclidean distances D between
    them and that the truck capacity C is  the the same for all vehicles. 
    """
    
    ## Assume planar case and convert to a sweep        
    sweep = get_sweep_from_cartesian_coordinates(points)

    ## Append each of the K customer cones into K groups of concecutive cones
    
    if C:
        alpha = sum(d)/float(K*C)
        group_target = alpha*C # = sum(d)/K
        EPS = C_EPS
    else: #only L set?
        total_sweep_len = sum( D[int(sweep[2][i-1]),int(sweep[2][i])]
                               for i in range(len(sweep[2])) )
        group_target = total_sweep_len/K
        EPS = S_EPS
        
    if __debug__:
        log(DEBUG-2,"Cone group demand/cost target = %.2f"%group_target )

    #for start_cone in range(len(cones)):    
    start_cone_idx = trial
    
    grouped_cones = []
    group_start_ray = None
    group_end_ray = None
    group_cum = 0.0
    group_nodes = []

    prev_node_i = None
    prev_node_rho = None
    prev_node_phi = sweep[0][start_cone_idx-1]
    if start_cone_idx==0:
        prev_node_phi-=2*pi 
    prev_ray = None
    
    # iterate over all (phi,rho,node_idx) staring from start_cone_idx
    #  and doing it twice
    for circle_view in (sweep.T[start_cone_idx:], sweep.T[:start_cone_idx+1]):
        for node_phi,node_rho,i in circle_view:
            i = int(i) # is numpy float
            if (node_phi<prev_node_phi):
                node_phi+=2*pi
            ray = bisect_angle(prev_node_phi,node_phi)
             
            if prev_ray is None:
                group_start_ray = ray
                if __debug__:
                    log(DEBUG-2,"First node %d cone sets group_start_ray=%.2f"%(i,group_start_ray))
            else:
                # calculate if the entire cone (~customer) can be added to the group
                #  or if only a fraction is needed to fill the group.
                if C:
                    cone_fraction = 1.0
                    if d[prev_node_i]!=0:
                        cone_fraction = min(1.0, (group_target-group_cum)/d[prev_node_i])
                    cone_wt = cone_fraction*d[prev_node_i]
                else:
                    cone_fraction = min(1.0, (group_target-group_cum)/(D[prev_node_i,i]))
                    cone_wt = cone_fraction*D[prev_node_i,i]
                
                group_cum+=cone_wt
                group_nodes.append( (prev_node_rho,prev_node_i,
                    d[prev_node_i] if C else D[prev_node_i,i]) ) 
                
                if __debug__:
                    if C:
                        log(DEBUG-3,"Node %d, added %.2f %% of demand (%.2f)" %\
                            (prev_node_i, cone_fraction*100, d[prev_node_i]))
                    else:
                        log(DEBUG-3,"Node %d, added %.2f %% of cost (%.2f)" %\
                            (prev_node_i, cone_fraction*100, 0.5*D[prev_node_i,i]))
                    log(DEBUG-2,"Group %.2f %% full"%\
                        (group_cum/group_target*100.0))
                                    
                if (group_target-group_cum)<EPS:  
                    group_end_ray = bisect_angle(prev_ray, ray, cone_fraction)  
                    # group is full, store it
                    grouped_cones.append( _Cone(group_start_ray,group_end_ray,
                                                group_cum, group_nodes) )
                                               
                    if __debug__:
                        log(DEBUG-2,"Node %d cone sets group_end_ray=%.2f"%\
                            (prev_node_i,group_end_ray))
                        log(DEBUG-2,"Group completed!\n")
                            
                    # next group            
                    group_start_ray = group_end_ray
                    group_nodes = []
                    group_cum = 0
                    
                    if cone_fraction<1.0:
                        if C:
                            rmdr_wt = (1.0-cone_fraction)*d[prev_node_i]
                        else:
                            rmdr_wt = (1.0-cone_fraction)*D[prev_node_i,i]
                        
                        group_cum += rmdr_wt
                        group_nodes.append((prev_node_rho,prev_node_i,
                            d[prev_node_i] if C else D[prev_node_i,i]))

                    if __debug__:
                        if len(grouped_cones)<K:
                            log(DEBUG-2,"Node %d cone sets group_start_ray=%.2f"%\
                                (prev_node_i,group_start_ray))
               
                # the group now spans upto this
                group_end_ray = ray
                
                if __debug__:
                    if len(grouped_cones)<K:
                        log(DEBUG-2,"Node %d cone grows group to ray=%.2f"%\
                            (prev_node_i,group_end_ray))
                    
            prev_ray = ray
            prev_node_i = i
            prev_node_rho = node_rho
            prev_node_phi = node_phi

    ## get seed form the resulting K merged cones
    seed_points = np.zeros((K,2), dtype=np.float64)
    
    depot_x = points[0][0]
    depot_y = points[0][1]
    for k, grouped_cone in enumerate(grouped_cones):
        if __debug__:
            log(DEBUG-3," ===========================================")
            log(DEBUG-3," #%d %s"%(k, str(grouped_cone)))
            log(DEBUG-3," ===========================================\n")
                
        # Find an arc that splits the k-cone in a way that the linear demand 
        #  under the arc is "around" 0.75 (the exact definition is in the
        #  Fisher & Jaikumar (1981) paper. Begin by sorting by distance from
        #  the depot and grow arc as long as weight sum is under the limit.
        seed_rho = 0
        grow_arc_wt = 0
        weight_target = 0.75*group_target # 0.75{\labmda}b
        for cr,ci,cwt in sorted(grouped_cone.nodes):
            if grow_arc_wt+cwt>weight_target:
                # take a fraction of the weight just outside the arc
                seed_rho+=((weight_target-grow_arc_wt)/cwt)*(cr-seed_rho)
                break
            else:
                grow_arc_wt+=cwt
                seed_rho=cr
        
        # Calculate the actual seed point position
        seed_phi = bisect_angle(grouped_cone.phi1,grouped_cone.phi2)
        seed_points[k,0] = depot_x+seed_rho*np.cos(seed_phi)
        seed_points[k,1] = depot_y+seed_rho*np.sin(seed_phi)
    return seed_points.tolist()

def _kmeans_seed_points(points, D, d, C, K, trial=0):
    """A seed point generation function that puts the seed points at customer
    node point cluster centers using k-Means clustering."""
    
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=K, random_state=trial).fit(points[1:])
    return kmeans.cluster_centers_.tolist()
    
def _end_of_thoroughfares_seed_points(points, D, d, C, K, trial=0):
    """A seed point generation function that automates the human assisted 
    idea presented in Fisher and Jaikumar (1981) involving placing the seed 
    points to the end of throughtfares leaving from the depot. A DBSCAN
    clustering is made and the seeds are selected among non-core points. Non-
    core points should be, due to the operating principle of DBSCAN, at the
    ends of long cluster "arms". By selecting the non-core points farthest from
    the depot and previously selected seeds, we should get a set of seed points
    closely following the Fisher and Jaikumar (1981) idea: "customers
    often lie along radial corridors corresponding to major thoroughfares, and
    the most distant ... along these corridors are natural seed customers".
    Fisher and Jaikumar (1981) presented the idea interactive computer systems
    in mind, whereas this implementation is automatic.
    
    TODO: in practice, the results are underwhelming. Instead, one should do
    1d clustering for phis and then choose the farthest point of each
    "Sweep cluster".
    
    parameters:
    - points, D, d, C, K as before
    - trial can be used to get different clusterings from the DBSCAN algorithm.
        the DBSCAN min_size is 2,2,3,3,4,4,... for trial 0,1,2,3,4,5... .
        The inititial eps is determined by getting the median distance of the 
        nn=2., 3., 2., 3., 3., 4., 3,... nearest neightbour of all nodes 
        depending if the trial is 0,1,2,3,4,5,6,7.. following the formula
            
            nn=2+trial%2+int(trial/4))
        
        The seed points are selected among the non-core points S_nc by
        maximizing the squared distances . If it 
        happens that |S_nc|<K, all non-core points are included and the rest
        of the seed points clustered points are 
        enough non-core points are found. 
        
    WARNING: This seed heuristic may return None seeds as the existence of non-
    core points cannot be guranteed.
    """
    
    from sklearn.cluster import DBSCAN
    from util import produce_nn_list
    
    # use a heuristic to get eps that finds all 2. closest nodes and 
    #  uses the median distance of those as the eps
    N = len(d)
    nnD = produce_nn_list(D)
    nn = 2+trial%2+int(trial/4)
    nn2l = [nnS[nn][0] for nnS in nnD]
    nn2l.sort()
    min_size = 3#+int(trial/2)
    eps = nn2l[int(N/2)]

    ## Get non-core DBSCAN points    
    if __debug__:
        log(DEBUG-2,"Doing DBSCAN with eps =", eps, " min_size =",min_size)
    db = DBSCAN(eps=eps, min_samples=min_size).fit(points)   
    outliers_mask = db.labels_ == -1
    clustered_mask = db.labels_ != -1
    core_samples_mask = np.zeros(N, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    # we are interested of the nodes at the fringes of the clusters    
    candidate_mask = clustered_mask^core_samples_mask
    candidate_idxs = np.where(candidate_mask)[0].tolist()
    candidates_type = "cluster non-core"
    
    if __debug__:
        log(DEBUG-3,"DBSCAN labels = %s"%str(list(zip(range(N),db.labels_))))
        log(DEBUG-3,"DBSCAN core = %s"%str(db.core_sample_idxs_))
        log(DEBUG-2,"Select %d seed nodes from non-core nodes %s."%
            (min(len(candidate_idxs),K), str(candidate_idxs)))
    seeds = []
    selected_seeds_mask = np.zeros(N, dtype=bool)
    # make depot like a seed -> maximize distance from it
    selected_seeds_mask[0] = True     
    if len(candidate_idxs)<=K:
        # if all candidates are needed, add them without checking the distances
        for seed_idx in candidate_idxs:
            seeds.append( points[seed_idx] )
            if __debug__:
                log(DEBUG-2,"Selecting n%d (%.2f, %.2f) that is a %s point to be a seed"%
                    (seed_idx,points[seed_idx][0],points[seed_idx][1],candidates_type))
            selected_seeds_mask[seed_idx] = True
        candidate_idxs = []
    
    used_core_points = False
    while len(seeds)<K:
        if not candidate_idxs:
            if not used_core_points:
                # ran out of non-core candidates. Use clustered as candidates
                candidate_mask = core_samples_mask
                candidate_idxs = np.where(core_samples_mask)[0].tolist()
                candidates_type = "cluster core"
                used_core_points = True
                
                if __debug__:
                     log(DEBUG-3,"Ran out of non-core nodes, select %d seed nodes from core nodes %s"%
                         (min(len(candidate_idxs), K-len(seeds)), str(candidate_idxs)))
            else:
                candidate_mask = outliers_mask
                candidate_idxs = np.where(outliers_mask)[0].tolist()
                candidates_type = "outliers"
                
                if __debug__:
                    log(DEBUG-3, "Ran out of core and non-core nodes, select %d seed nodes from outlier nodes %s"%
                        (K-len(seeds), str(candidate_idxs)))
            
        # maximize the distance to other seeds and depot
        if not seeds:
            D_to_seeds = D[selected_seeds_mask,candidate_mask]
        else:
            D_to_seeds = np.sum( np.sqrt((D[selected_seeds_mask,:])[:,candidate_mask]), axis=0)
        seed_idx = candidate_idxs[np.argmax( D_to_seeds )]        
        selected_seeds_mask[seed_idx] = True
        seeds.append( points[seed_idx] )    
        
        if __debug__:
            log(DEBUG-2, "Selecting n%d (%.2f, %.2f) that is a %s point to be a seed"%
                (seed_idx,points[seed_idx][0],points[seed_idx][1], candidates_type))
                     
        # prevent selecting it again 
        candidate_mask[seed_idx] = False
        candidate_idxs.remove(seed_idx)

    return seeds
        
def _large_demand_seed_points(points, D, d, C, K, trial=0):
    """A seed point generation function that automates the human assisted 
    idea presented in Fisher and Jaikumar (1981)
    """
    # make sure we are dealing with np arrays here
    np_d = np.array(d)
    N = len(d)
    
    # we are look mainly the large d nodes where only 1 fits on a route
    can_fit_only_1_mask = np_d > (0.5*C)
    candidate_d_mask = can_fit_only_1_mask.copy()
    candidate_d_idxs = np.where(can_fit_only_1_mask)[0].tolist()
    
    if trial:
        # in addition, add as many OTHER largest d ones as trial is
        not_over_half_idxs = np.where( ~candidate_d_mask )[0].tolist()
        sorted_d = [(d[i], i) for i in not_over_half_idxs]
        sorted_d.sort(reverse=True)
        sorted_d_idxs = list(zip(*sorted_d)[1])
        additional_large_d_idxs = sorted_d_idxs[max(0, trial-N):min(N,trial)]
        candidate_d_idxs+=additional_large_d_idxs
        candidate_d_mask[additional_large_d_idxs] = True
    
    large_d_mask = np.copy(candidate_d_mask)
    
    if __debug__:    
        log(DEBUG-2, "Select %d seed nodes from large demand nodes %s"%
            (min(len(candidate_d_idxs),K), str(candidate_d_idxs)))
    
    seeds = []
    selected_seeds_mask = np.zeros(len(d), dtype=bool)
    # make depot like a seed -> maximize distance from it
    selected_seeds_mask[0] = True
    if len(candidate_d_idxs)<=K:
        # if all candidates are needed, add them without checking the distances
        for seed_idx in candidate_d_idxs:
            seeds.append( points[seed_idx] )    
            selected_seeds_mask[seed_idx] = True        
            if __debug__:
                log(DEBUG-2,"Selecting n%d (%.2f, %.2f) that %s to be a seed"%\
                    (seed_idx,points[seed_idx][0],points[seed_idx][1],
                    "fills over the half of the capacity" if can_fit_only_1_mask[seed_idx]
                    else "is within "+str(trial)+" largest demands"))
        candidate_d_idxs = []
        
    
    select_from_non_large = False
    while len(seeds)<K:
        if not candidate_d_idxs:
            candidate_d_mask = ~large_d_mask
            candidate_d_mask[0]=False
            candidate_d_idxs = np.where(candidate_d_mask)[0].tolist()
            select_from_non_large = True
            
            if __debug__:
                log(DEBUG-2,"Ran out of nodes with large demand, select %d seed nodes from rest of the nodes %s using inter seed distances weighted by the node demand"%
                    (min(len(candidate_d_idxs), K-len(seeds)), str(candidate_d_idxs)))
            
        # maximize the distance to other seeds and depot
        if not seeds:
            D_to_seeds = D[selected_seeds_mask,candidate_d_mask]
        else:
            D_to_seeds = np.sum( np.sqrt((D[selected_seeds_mask,:])[:,candidate_d_mask]), axis=0)
        if select_from_non_large:
            # multiply by demand
            D_to_seeds = np.multiply(D_to_seeds,np_d[candidate_d_mask]/C)
            
        seed_idx = candidate_d_idxs[np.argmax( D_to_seeds )]        
        selected_seeds_mask[seed_idx] = True
        seeds.append( points[seed_idx] )    
        
        if __debug__:
            if can_fit_only_1_mask[seed_idx]:
                candidates_type = "fills over the half of the capacity"
            elif large_d_mask[seed_idx]:
                candidates_type = "is within "+str(trial)+" largest demands"
            else:
                candidates_type = "when weighted by demand has largest distance from other seeds"
            log(DEBUG-2,"Selecting a node n%d (%.2f, %.2f) that %s to be a seed"%\
                  (seed_idx,points[seed_idx][0],points[seed_idx][1], candidates_type))
                  
        # prevent selecting it again 
        candidate_d_mask[seed_idx] = False
        candidate_d_idxs.remove(seed_idx)

    return seeds
    
def gap_init(points, D, d, C, L=None, st=None, K=None, minimize_K=True,
             find_optimal_seeds=True,
             seed_method="cones",
             seed_edge_weight_type='EUC_2D',
             use_adaptive_L_constraint_weights=True,
             increase_K_on_failure=False):
             #REMOVEME, disable!
             #increase_K_on_failure=True):
    """ An implementation of a three phase cluster-first-route-second CVRP
    construction / route initialization algorithm. The first two phases involve
    the clustering. First, a seed point is generated for each route, which is
    then used in approximating customer node service costs in solving
    generalized assignment problem (GAP) relaxation of the VRP. The resulting
    assignments are then routed using a TSP solver. The algorithm has been 
    first proposed in (Fisher and Jaikumar 1981).
    
    The algorithm assumes that the problem is planar and this implementation
    allows seed in two ways: 
    * seed_method="cones", the initialization method of Fisher and Jaikumar
        (1981) which can be described as Sweep with fractional distribution of
        customer demand and placing the seed points approximately to the center
        of demand mass of created sectors.
    * seed_method="kmeans", intialize seed points to k-means cluster centers.
    * seed_method="large_demands", according to Fisher and Jaikumar (1981) 
        "Customers for which d_i > 1/2 C can also be made seed customers". 
        However applying this rule relies on human operator who then decides
        the intuitively best seed points. This implementation selects the
        seed points satisfying the criteria d_i>mC, where m is the fractional
        capacity multipier, that are farthest from the depot and each other.
        The m is made iteratively smaller if there are no at least K seed point
        candidates.
    * seed_method="ends_of_thoroughfares", this option was descibed in 
        (Fisher and Jaikumar 1981) as "Most distant customers at the end of
        thoroughfares leaving from the depot are natural seed customers". They
        relied on human operator. To automate this selection we make a 
        DBSCAN clustering with eps = median 2. nearest neighbor of all nodes
        and min_samples of 3. 
        
    
    The other parameters are:
    * points is a list of x,y coordinates of the depot [0] and the customers.
    * D is a numpy ndarray (or equvalent) of the full 2D distance matrix.
       including the service times (st/2.0 for leaving and entering nodes).
    * d is a list of demands. d[0] should be 0.0 as it is the depot.
    * C is the capacity constraint limit for the identical vehicles.
    * L is the optional constraint for the maximum route length/duration/cost.
    * st is the service time. However, also the D should be modified with 
       service times to allow straight computation of the TSP solutions (see
       above)
    * K is the optional parameter specifying the required number of vehicles.
       The algorithm is only allowed to find solutions with this many vehicles.
    * minimize_K, if set to True (default), makes the minimum number of routes
       the primary and the solution cost the secondary objective. If set False
       the algorithm optimizes for mimimum solution / route cost by increasing
       K as long as it seems beneficial. WARNING: the algorithm suits this use
       case (cost at the objective) poorly and setting this option to False may
       significantly increase the required CPU time.
       
    * find_optimal_seeds if set to True, tries all possible Sweep start
       positions / k-Means with N different seeds. If False, only one sweep 
       from the node closest to the depot is done / k-Means clustering is done
       only once with one random seed value.
    * seed_edge_weight_type specifies how to round off the distances from the
       customer nodes (points) to the seed points. Supports all TSPLIB edge
       weight types.
       
    Note1: The GAP is optimized using Gurobi solver. If L constraint is set,
     the side constraints may make the GAP instance tricky to solve and it 
     is advisable to set a sensible timeout with config.MAX_MIP_SOLVER_RUNTIME
    * use_adaptive_L_constraint_weights if set True, and the L constraint is 
     set, the algorithm adaptively adjusts the route cost approximation of the
     relevant side constraints so that a solution which is not L infeasible or
     GAP infeasible is found. The exact handling of L consraint is vague in
     (Fisher and Jaikumar 1981) and this was our best guess on how the
     feasible region of the problem can be found. Note that if GAP solver is
     terminated due to a timeout, the adaptive multipier is increased and 
     GAP solution is attempted again. However, if increase_K_on_failure is set,
     (see below) it takes priority over this.
    * increase_K_on_failure (default False) is another countermeasure against
     long running GAP solving attempts for problem instances without L 
     constraint (if there is L constraint, and use_adaptive_L_constraint_-
     weights is enabled, this is ignored) or instances where K estimation 
     does not work and it takes excessively long time to check all initial 
     seed configurations before increasing K. If Gurobi timeout is encountered
     or the solution is GAP infeasible, and this option is enabled, the K is
     temporately increased, new seeds points generated for current sweep start
     location and another GAP solution attempt is made. K is allowed to
     increased temporarely up to 10% of the mimimum K allowed (or 1, whichever
     is larger).
    
    Note2: logger controls the debug level but running the script with
     Python -O option disables all debug output.
    
    Fisher, M. L. and Jaikumar, R. (1981), A generalized assignment heuristic
    for vehicle routing. Networks, 11: 109-124. doi:10.1002/net.3230110205
    """   #TODO: other alternatives
    # customers with maximum demand or most distant customer from origin
    
    if seed_method=="cones":
        seed_f = _sweep_seed_points
    if seed_method=="kmeans":
        seed_f = _kmeans_seed_points
    if seed_method=="large_demands":
        if not C: raise ValueError("""The "large_demands" seed initialization method requires demands and C constraint to be known.""")
        seed_f = _large_demand_seed_points
    if seed_method=="ends_of_thoroughfares":
        seed_f = _end_of_thoroughfares_seed_points
    
    int_dists = issubclass(D.dtype.type, np.integer)
    if seed_edge_weight_type=="EXPLICIT":
        seed_edge_weight_type = "EUC_2D" if int_dists else "EXACT_2D"
    
    if not points:
        raise ValueError("The algorithm requires 2D coordinates for the points")
    N = len(D)    
    if K:
        startK = K
        maxK = K
    else:
        # start from the smallest K possible
        if C:
            startK = int(ceil(sum(d)/C))
        elif L:
            # find a lower bound by checking how many visits from the TSP
            #  tour need to add to have any chance of making this L feasible.
            _,tsp_f = solve_tsp(D, list(range(1,N)))
            shortest_depot_edges = list(D[0,1:])
            shortest_depot_edges.sort()
            startK = int(ceil(tsp_f/L))
            while True:
                if tsp_f+sum(shortest_depot_edges[:startK*2])<=startK*L:
                    break
                startK+=1
        else:
            raise ValueError("If C and L have not been set, K is required")   
        maxK = N-1
    
    # We only need first row of the distance matrix to calculcate insertion 
    #  costs for GAP objective function
    D_0 = np.copy( D[0,:] )
    
    best_sol = None
    best_f = None
    best_K = None
    seed_trial = 0
    incK = 0
    maxKinc = max(startK+1, int(startK*INCREASE_K_ON_FAILURE_UPTO))
    
    L_ctr_multipiler = L_MPLR_DEFAULT
    if L and use_adaptive_L_constraint_weights:
        # Adaptive L constraint multipier 
        L_ctr_multipiler = L_ADAPTIVE_MPLR_INIT
        L_ctr_multipiler_tries = 0
    
    try:
        for currentK in range(startK, maxK+1):
            found_improving_solution_for_this_K = False
            seed_trial=0
            while True:
                if __debug__:
                    log(DEBUG, "ITERATION:K=%d, trial=%d, L_ctr_mul=%.6f\n"%
                        (currentK+incK,seed_trial,L_ctr_multipiler))
                    log(DEBUG-1, "Getting %d seed points...\n"%(currentK+incK))
                
                # Get seed points
                seed_points = seed_f(points, D, d, C, currentK+incK, seed_trial)
                if __debug__:
                    log(DEBUG-1, "...got seed points %s\n"%str(seed_points))
                
                # Extend the distance matrix with seed distances
                
                S = calculate_D(seed_points, points, seed_edge_weight_type)
                if st:
                    # include the "leaving half" of the service_time in the 
                    #  distances (the other half is already added to the D
                    #  prior to gapvrp_init)
                    halftst = int(st/2) if int_dists else st/2.0
                    S[:,1:] += halftst
                D_s = np.vstack( (D_0, S) )
        
                GAP_infeasible = False        
                L_infeasible = False
                solution = [0]
                sol_f = 0
                solved = False
                sol_K = 0
                take_next_seed = False
                try:
                    # Distribute the nodes to vehicles using the approxmate 
                    # service costs in D_s and by solving it as GAP
                    #
                    #TODO: the model has the same dimensions for all iterations
                    # with the same K and only the weights differ. Consider
                    # replacing the coefficient matrix  e.g. via C interface
                    #https://stackoverflow.com/questions/33461329
                    assignments = _solve_gap(N, D_s, d, C, currentK+incK, L,
                                             L_ctr_multipiler)
                    if not assignments:
                        if __debug__:
                            log(DEBUG, "INFEASIBILITY: GAP infeasible solution")
                            corrective_action = "try with another seed = %d"%seed_trial
                        GAP_infeasible = True                
                    else:
                        if __debug__:
                            log(DEBUG-1, "Assignments = %s"%str(assignments))
                        
                        # Due to floating point inaccuracies in L constrained
                        #  cases the feasrelax may be used, which, in turn, can
                        #  in some corner cases return solutions that are not
                        #  really feasible. Make sure it is not the case
                        if L: served = set([0])
                                                
                        for route_nodes in assignments:
                            if not route_nodes:
                                continue
                            route,route_l = solve_tsp(D, [0]+route_nodes)
                        
                            # Check for feasibility violations due to feasrelax
                            if L:
                                served |= set(route_nodes)
                                if C and d and totald(route,d)-C_EPS>C:
                                    if __debug__: 
                                        log(DEBUG, "INFEASIBILITY: feasRelax "+
                                            "caused GAP infeasible solution "+
                                            " (capacity constraint violation)")
                                    GAP_infeasible = True
                                    break # the route loop
                                        
                            solution += route[1:]
                            sol_f += route_l
                            sol_K += 1
                            
                            if __debug__:
                                log(DEBUG-2, "DEBUG: Got TSP solution %s (%.2f)"%
                                    (str(route),route_l))
                                
                            if L and route_l-S_EPS>L:
                                if __debug__:
                                    log(DEBUG, "INFEASIBILITY: L infeasible solution")
                                L_infeasible = True
                                break # break route for loop
                                
                        # Check for feasibility violations due to feasrelax.
                        #  Have all customers been served?
                        if not GAP_infeasible and not L_infeasible and\
                           L and len(served)<len(D):
                            if __debug__: 
                                log(DEBUG, "INFEASIBILITY: feasRelax caused GAP "+
                                           "infeasible solution (all customers "+
                                           "are not served)")
                            GAP_infeasible = True                        
                        
                    if not GAP_infeasible and not L_infeasible:
                        if __debug__:
                            log(DEBUG, "Yielded feasible solution = %s (%.2f)"%(str(solution), sol_f))
                        solved = True
                        
                except GurobiError as grbe:
                    if __debug__: log(WARNING, str(grbe))

                    if L and use_adaptive_L_constraint_weights and \
                         L_ctr_multipiler_tries<L_ADAPTIVE_MPLR_MAX_TRIES:
                       L_ctr_multipiler+=L_ADAPTIVE_MPLR_INC
                       L_ctr_multipiler_tries+=1
                       if __debug__: corrective_action = "Gurobi timeout, try with another L_ctr_multipiler = %.2f"%L_ctr_multipiler
                    elif increase_K_on_failure and currentK+incK+1<=maxKinc:
                        if L and use_adaptive_L_constraint_weights and\
                           L_ctr_multipiler_tries>=L_ADAPTIVE_MPLR_MAX_TRIES:
                            # try with all multiplier values for larger K
                            L_ctr_multipiler = L_ADAPTIVE_MPLR_INIT
                            L_ctr_multipiler_tries = 0
                        incK+=1
                        if __debug__: corrective_action = "Gurobi timeout, temporarely increase K by %d"%incK
                    elif find_optimal_seeds:
                       take_next_seed = True
                    else:
                        grbe.message+=", consider increasing the MAX_MIP_SOLVER_RUNTIME in config.py"
                        raise grbe
                else:
                    if L and use_adaptive_L_constraint_weights:
                        ## Adaptive GAP/L constraint multiplier reset 
                        # reset multiplier in case it the L feasibility was not violated
                        #  or it has reached the max_value. 
                        if solved or L_ctr_multipiler_tries>=L_ADAPTIVE_MPLR_MAX_TRIES:
                            L_ctr_multipiler = L_ADAPTIVE_MPLR_INIT
                            L_ctr_multipiler_tries = 0
                            take_next_seed = True
                            if not solved and increase_K_on_failure and currentK+incK+1<=maxKinc:
                                incK+=1
                                take_next_seed = False
                                if __debug__: corrective_action = "temporarely increase K by %d"%incK
                            else:
                                if __debug__: corrective_action = "try with another seed = %d"%seed_trial
                        ## Adaptive GAP/L constraint multiplier update
                        else:
                            L_ctr_multipiler+=L_ADAPTIVE_MPLR_INC
                            L_ctr_multipiler_tries+=1
                            if __debug__: corrective_action = "try with another L_ctr_multipiler = %.2f"%L_ctr_multipiler
                    else:
                        if not solved and increase_K_on_failure and currentK+incK+1<=maxKinc:
                            incK+=1
                            if __debug__: corrective_action = "temporarely increase K by %d"%incK
                        else:
                            take_next_seed = True


                # Store the best so far
                if solved:
                    if is_better_sol(best_f, best_K, sol_f, sol_K, minimize_K):
                        best_sol = solution
                        best_f = sol_f
                        best_K = sol_K
                        found_improving_solution_for_this_K = True
                else:
                    # No feasible solution was found for this trial (max route cost 
                    #  or capacity constraint was violated).            
                    if __debug__:
                        if GAP_infeasible or L_infeasible:
                            log(DEBUG, "Constraint is violated, "+corrective_action)
                        else:
                            log(DEBUG, "Continuing search, "+corrective_action)
    
                if take_next_seed:
                    incK = 0
                    seed_trial+=1
                    if not find_optimal_seeds: 
                        break # seed loop, possibly try next K
                if seed_trial==N:
                    incK = 0
                    break # seed loop, possibly try next K
            
            if minimize_K:
                # do not try different K if we found a solution
                if best_sol:
                    break # K loop
            else: # not minimize_K
                # We already have an feasible solution for K<K_current, and could 
                # not find a better solution than that on K_current. Therefore, it 
                # is improbable we will find one even if we increase K and we
                # should stop here.
                if best_sol and not found_improving_solution_for_this_K:
                    break
    except KeyboardInterrupt: #or SIGINT
        #  pass on the current best_sol
        raise KeyboardInterrupt(best_sol)
        
    return best_sol                

# ---------------------------------------------------------------------
# Wrapper for the command line user interface (CLI)
def get_gap_algorithm(seed_method="cones"):
    algo_name = "FJ81-GAP"
    algo_desc = "Fisher & Jaikumar (1981) generalized assignment problem heuristic"
    def call_init(points, D, d, C, L, st, wtt, single, minimize_K):
        return gap_init(points, D, d, C, L=L, st=st,
                        K=None, minimize_K=minimize_K,
                        seed_edge_weight_type=wtt,
                        find_optimal_seeds=(not single),
                        seed_method=seed_method)
    return (algo_name, algo_desc, call_init)
    
if __name__=="__main__":
    from shared_cli import cli
    cli(*get_gap_algorithm())
