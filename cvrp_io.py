# -*- coding: utf-8 -*-
################################################################################
""" This file implements the necessary functionality for reading TSPLIB CVRP
problem instance files, additional constraints from the said files, and
generating new random instances.
"""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
from builtins import range

import os
import re
import random

from collections import namedtuple
from math import pi, radians, cos, sin, asin, sqrt, acos, modf
from itertools import groupby   
from sys import stderr

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2018, Jussi Rasku"
__credits__ = ["Jussi Rasku"]
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = "Jussi Rasku"
__email__ = "jussi.rasku@jyu.fi"
__status__ = "Development"

################################################################################

k_re = re.compile(r"-k([0-9]+)[\.-]")

def _haversine(pt1, pt2):
    """from http://stackoverflow.com/questions/4913349/
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    The distance should be within ~0.3% of the correct value.
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [pt1[0], pt1[1], pt2[0], pt2[1]])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

def _degrees_and_minutes_to_radians(x):
    """ Adapted from Reneilt 1991 TSPLIB article / TSPFAQ """
    PI = 3.141592
    mins, degs = modf(x)
    return (degs+100/60.0*mins)*PI/180.0

def _geo(pt1, pt2):
    """ Adapted from Reneilt 1991 TSPLIB article / TSPFAQ
    this togehter with the _degrees_and_minutes_to_radians conversion produces
    the same results than the optimal solution on the original GEO TSP files."""
    RRR = 6378.388

    latitude_i_rads, longitude_i_rads = pt1
    latitude_j_rads, longitude_j_rads = pt2
    
    q1 = cos(longitude_i_rads - longitude_j_rads)
    q2 = cos(latitude_i_rads - latitude_j_rads)
    q3 = cos(latitude_i_rads + latitude_j_rads)
    return int( RRR*acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0 )

def _att(pt1,pt2):
    dx = pt1[0]-pt2[0]
    dy = pt1[1]-pt2[1]
    r = sqrt(dx**2+dy**2)/10.0
    t = int(r)
    return t+1 if t<r else t

def calculate_D(pts, opts=None, tsplib_distances_type='EUC_2D'):
    pdtype = 'euclidean'
    postprocess = lambda M: M
    
    if tsplib_distances_type=='MAX_2D':
        pdtype = 'chebyshev'
    elif tsplib_distances_type=='MAN_2D':
        pdtype = 'cityblock'
    elif tsplib_distances_type=='CEIL_2D':
        postprocess = lambda D: np.ceil(D).astype(int)
    elif tsplib_distances_type=='FLOOR_2D':
        postprocess = lambda D: np.floor(D).astype(int)
    elif tsplib_distances_type=='EUC_2D':
        postprocess = lambda D: np.round(D).astype(int)
    elif tsplib_distances_type=='ATT':
        pdtype = lambda v,w : _att(v, w)
    elif tsplib_distances_type=='GEO':
        pdtype = lambda v,w : _geo(v, w)
    elif tsplib_distances_type=='EXACT_2D':
        pass
    else:
        raise ValueError("Unknown distance method")
    if opts is None:
        return postprocess(squareform(pdist(pts, pdtype)))
    else:
        return postprocess(cdist(pts, opts, pdtype))

def read_OPT_CVRP(file_name):    
    solution = [0]
    opt_f = None
    opt_k = None
    
    re_k = k_re.findall(file_name)
    if re_k:
        opt_k = int(re_k[0])
    
    
    file_ext = os.path.splitext(file_name)[1]

    count_k = 0
    with open(file_name, "r") as f:
        for l in f.readlines():
            if file_ext == ".opt":
                if "route" in l.lower():
                    if not opt_k:
                        count_k+1
                    _, routestring = l.split(":")
                    p_idxs = [int(s) for s in routestring.split()]
                    first_node = True
                    for p_idx in p_idxs:
                        if first_node and solution[-1]!=0:
                            solution.append(0)
                        solution.append(p_idx)
                        first_node = False
                if "cost" in l.lower():
                    _, coststring = l.split()
                    
                    # tries to convert to int and if it fails to float
                    opt_f = None
                    try:
                        opt_f = int(coststring)
                    except ValueError:
                        opt_f = float(coststring)
            else:
                raise NotImplementedError("This solution file is not supported (yet)")
    if len(solution)>1:
        solution.append(0)
        
    if count_k or not opt_k:
        opt_k = count_k
    elif opt_k!=count_k:
        print("WARNING: the vehicle count in file name and solution differ", file=stderr)

    return solution, opt_f, opt_k
         
       
ProblemDefinition = namedtuple('ProblemDefinition',
    ['size', 'coordinate_points', 'display_coordinate_points',
     'customer_demands', 'distance_matrix', 'capacity_constraint', 'edge_weight_type'])
def read_TSPLIB_CVRP(file_name):
    """ Returns a namedtuple (N, points, dd_points, demands, D, C, ewt) where
    * N is the size of the problem,
    * points has the coordinates of the depot (index 0) and customers,
        note: points can be None if the file does not have NODE_COORD_SECTION
    * dd_points has the DISPLAY coordinates,
        note: is usually None as files containing DISPLAY_DATA_SECTION are rare
    * demands is a list of demands with the depot demand (index 0) set to 0
    * D is the distance matrix as a numpy 2D ndarray,
    * C is the vehicle capacity constraint, can be None if it is not set
    * ewt is the EDGE_WEIGHT_TYPE
    
    The reader supports following TSPLIB (Reinelt, 1991) fields:
        NAME
        TYPE
        DIMENSION
        CAPACITY
        EDGE_WEIGHT_FORMAT (FUNCTION/FULL_MATRIX/
                            LOWER_ROW/LOWER_DIAG_ROW/
                            UPPER_ROW/UPPER_DIAG_ROW/
                            LOWER_COL)
        EDGE_WEIGHT_TYPE (MAX_2D/MAN_2D/EXACT_2D/CEIL_2D/EUC_2D/EXPLICIT/GEO/ATT)
        NODE_COORD_TYPE
        
    and sections:
        EDGE_WEIGHT_SECTION
        NODE_COORD_SECTION
        DEMAND_SECTION
        DEPOT_SECTION
        DISPLAY_DATA_SECTION
        
    However, these are ignored (but see read_TSBLIB_additional_constraints):
        SVC_TIME_SECTION
        DISTANCE
        SERVICE_TIME

    Reinelt, G. (1991). Tsplib a traveling salesman problem library. ORSA 
        journal on computing, 3(4):376-384
    """
    with open(file_name, "r") as f:
        # pylint: disable=unsubscriptable-object
        
        section = None
        section_pos = 0
        ij_section_pos = {'i':0,'j':0}
        N=0
        C=None
        
        points = None  
        dd_points = None
        demands = None
        D = None
        D_needs_update = False     
        edge_weight_type = None
        edge_weight_format = None
        
        depot_ids = []
        
        while 1:
            line = f.readline().strip()
            if not line:
                continue
            
            # Parse fields
            if ':' in line:
                field, value = line.split(":",1)
                field = field.strip()
                
                if 'TYPE' == field:
                    if not 'CVRP' in value and not 'TSP' in value:
                        raise IOError("Only CVRP TSPLIB files are supported")
                elif 'DIMENSION' in field:
                    N = int(value)-1 # depot excluded
                elif 'CAPACITY' in field:
                    C = int(value)
                elif 'EDGE_WEIGHT_TYPE' in field:
                    edge_weight_type = value.strip()
                    if edge_weight_type not in ["MAX_2D", "MAN_2D", "EXACT_2D",
                                                "CEIL_2D", "FLOOR_2D", "EUC_2D",
                                                "EXPLICIT", "GEO", "ATT"]:
                        raise IOError("Only matrix and euclidian distance notation is supported")
                elif 'EDGE_WEIGHT_FORMAT' in field:
                    edge_weight_format  = value.strip()
            
            # Section handling
            else:
                if 'EOF' in line:
                    break
            
                if 'EDGE_WEIGHT_SECTION' in line:
                    section = 'EDGE_WEIGHT_SECTION'
                    D = np.zeros((N+1,N+1))
                    ij_section_pos = {'i':0,'j':0}
                    if (edge_weight_format=="LOWER_ROW"):
                        ij_section_pos['j']=1
                    elif (edge_weight_format=="UPPER_ROW" or
                        edge_weight_format=="LOWER_COL"):
                        ij_section_pos['i']=1
                elif 'DEMAND_SECTION' in line:
                    demands = [None]*(N+1)
                    section = 'DEMAND_SECTION'
                    section_pos = 0
                elif 'DEPOT_SECTION' in line:
                    section = 'DEPOT_SECTION'
                    section_pos = 0
                elif 'NODE_COORD_SECTION' in line:
                    section = 'NODE_COORD_SECTION'
                    points = [ [None, None] for i in range(N+1) ]
                    if edge_weight_type!='EXPLICIT':
                        # sometimes coordinates are incorrectly not given in 
                        #  DISPLAY_DATA_SECTION even if a matrix is defined.
                        D_needs_update = True
                    section_pos = 0
                elif 'DISPLAY_DATA_SECTION' in line:
                    if points is None:
                        section = 'DISPLAY_DATA_SECTION'
                        dd_points = [ [None, None] for i in range(N+1) ]
                        D_needs_update = False
                        section_pos = 0
                    else:
                        section = ''
                elif 'SVC_TIME_SECTION' in line:
                    section = 'SVC_TIME_SECTION'
                else:
                    if section == 'EDGE_WEIGHT_SECTION':
                        distances = line.split() 
                        #print distances, section_pos, edge_weight_format
                        for d in distances: 
                            
                            D[ij_section_pos['i']][ij_section_pos['j']] = float(d)
                            D[ij_section_pos['j']][ij_section_pos['i']] = float(d)   
                            
                            if (edge_weight_format=="LOWER_ROW"):
                                # incrementer
                                ij_section_pos['i']+=1
                                if ij_section_pos['i']==ij_section_pos['j']:
                                    ij_section_pos['i'] = 0
                                    ij_section_pos['j'] += 1                            
                            elif (edge_weight_format=="UPPER_ROW" or 
                                  edge_weight_format=="LOWER_COL"):
                                # incrementer
                                ij_section_pos['i']+=1
                                if ij_section_pos['i']==len(D):
                                    ij_section_pos['j'] += 1
                                    ij_section_pos['i'] = ij_section_pos['j']+1                                    
                            elif (edge_weight_format=="FULL_MATRIX"):
                                # incrementer
                                ij_section_pos['i']+=1
                                if ij_section_pos['i']==len(D):
                                    ij_section_pos['j'] += 1
                                    ij_section_pos['i'] = 0
                            elif (edge_weight_format=="LOWER_DIAG_ROW"):
                                # incrementer
                                ij_section_pos['i']+=1
                                if ij_section_pos['i']==ij_section_pos['j']+1:
                                    ij_section_pos['i'] = 0
                                    ij_section_pos['j'] += 1
                            elif (edge_weight_format=="UPPER_DIAG_ROW"):
                                # incrementer
                                ij_section_pos['i']+=1
                                if ij_section_pos['i']==len(D):
                                    ij_section_pos['j'] += 1
                                    ij_section_pos['i'] = ij_section_pos['j']                            
                            
                            
                    elif section == 'NODE_COORD_SECTION':
                        coords = line.split()
                        x = float( coords [1] )
                        y = float( coords [2] )
                        
                        # According to TSPLIB format spec. the GEO coordinates
                        #  are of  format degrees.minutes. Convert to radians
                        #  BUT FIX THE ISSUE WITH THE NEGATIVE MINUTES THE
                        #  ORIGINAL SPEC HAS!
                        if edge_weight_type=='GEO':
                            x = _degrees_and_minutes_to_radians(x)
                            y = _degrees_and_minutes_to_radians(y)
                            #print("lat, lon (in rads) : %.2f, %.2f"%(x,y))
                        
                        points[section_pos][0] = x
                        points[section_pos][1] = y
                        section_pos+=1
                    elif  section == 'DISPLAY_DATA_SECTION':
                        coords = line.split()
                        x = float( coords [1] )
                        y = float( coords [2] )
                        dd_points[section_pos][0] = x
                        dd_points[section_pos][1] = y
                        section_pos+=1
                    elif  section == 'DEMAND_SECTION':
                        demand = line.split()
                        c = float( demand[1] )
                        demands[section_pos] = c # pylint: disable=unsupported-assignment-operation
                        section_pos+=1
                    elif  section == 'DEPOT_SECTION': 
                        value = int(line)
                        if value>0:
                            depot_ids.append(value)
                            if len(depot_ids)>1:
                                raise IOError("multi depot problems not supported")
                    
        f.close()
        
        if edge_weight_type=='EXPLICIT' and not (
                ((edge_weight_format in ['FULL_MATRIX', 'LOWER_ROW', 'LOWER_DIAG_ROW']) and \
                    ij_section_pos['i']==0 and ij_section_pos['j']==len(D)) or\
                  (edge_weight_format in ['UPPER_ROW','LOWER_COL'] and \
                    ij_section_pos['i']==len(D) and ij_section_pos['j']==len(D)-1) or\
                  (edge_weight_format == 'UPPER_DIAG_ROW' and \
                    ij_section_pos['i']==len(D) and ij_section_pos['j']==len(D))
                ):
            #print edge_weight_format, ij_section_pos 
            raise IOError("Explicit distance matrix did not have enough values")
    
    if D_needs_update:
        D = calculate_D(points, None, edge_weight_type )
    
    if edge_weight_type == "EXPLICIT":
        # check if the matrix had integer dinstances (as they often have)
        D_int = D.astype(int)
        if np.all((D - D_int) == 0):
            D = D_int
    
    # depot is not node 0!       
    if depot_ids and depot_ids[0]>1:
        # make sure depot is the 0
        idx_0 = depot_ids[0]-1
        row_col_permutation = [idx_0]+list(range(0,idx_0))+list(range(idx_0+1,len(D)))
        for i in range(N):
            D[:,i] = D[row_col_permutation,i]
        for i in range(N):
            D[i,:] = D[i,row_col_permutation]
        if demands is not None and len(demands)>0:
            demands = [demands[idx_0]]+demands[:idx_0]+demands[idx_0+1:]
        if points is not None and len(points)>0:
            points = [points[idx_0]]+points[:idx_0]+points[idx_0+1:]
        if dd_points is not None and len(dd_points)>0:
            dd_points = [dd_points[idx_0]]+dd_points[:idx_0]+dd_points[idx_0+1:]    
        
    if edge_weight_type=="GEO":
        dd_points = points
        points = None
    
    return ProblemDefinition(N, points, dd_points, demands, D, C, edge_weight_type)
 
namedtuple('AdditionalConstraints',
    'vehicle_count_constraint maximum_route_cost_constraint service_time_at_customer')
def read_TSBLIB_additional_constraints(custom_tsplib_file):
    """ An unofficial/custom and optional way of storing route cost/length/
    duration constraint in a TSBLIB file as an additional DISTANCE, VEHICLES
    and SERVICE_TIME fields (e.g. in CMT instances).
    
    Also SVC_TIME_SECTION is supported but only if the service time is set to
    the same value for all customers.
    """    
    K = None
    L = None
    ST = None
    reading_service_time_section = False
    with open(custom_tsplib_file) as fh:
        for l in fh.readlines():
            if reading_service_time_section:
                nid, nst = l.split()
                if "." in nst:
                    nst = float(nst)
                else:
                    nst = int(nst)
                if ST is not None and nst!=ST:
                    raise IOError("Only single (same) service time for all customers is supported")
                elif nid!=1:
                    ST = nst
            if "DISTANCE" in l:
                if "." in l:
                    L = float( l.split()[-1] )
                else:
                    L = int( l.split()[-1] )
            if "SERVICE_TIME" in l:
                if "." in l:
                    ST = float( l.split()[-1] )
                else:
                    ST = int( l.split()[-1] )
            if "VEHICLES" in l:
                K = int( l.split()[-1] )
            if "SVC_TIME_SECTION" in l:
                reading_service_time_section = True
                    
    return K, L, ST
 
def generate_CVRP(N, C, muC, sdC, regular=False, R=200.0):
    """ Generate new random CVRP with N customer points and capacity of C.
    Demand of customers is randomly generated with mean of muC and standard
    deviation sdC.
    returns (N, points,demands, D, C)
    """
    points = []
    demands = []
    points.append((0.0,0.0)) # Depot at 0,0
    demands.append(0)
    sumc = 0.0
    alpha = pi/4.0
    for _ in range(N):
        if regular:
            alpha+=(2*pi/N)
            r = R
        else:
            # Random angle
            alpha = random.random()*2*pi
            r = R*random.gauss(1.0, 0.33)
        pt_x = r*cos(alpha)
        pt_y = r*sin(alpha)
        c = min(C, max(1.0, random.gauss(muC, sdC)))
        sumc+=c
        points.append((pt_x, pt_y))
        demands.append(c)
    #points[0][2] = -sumc
        
    D = calculate_D(points)
           
    return ProblemDefinition(N,points,None,demands,D,C,None)

def as_VRPH_solution(sol):
    """ Return a string containing the solution in the format used by VRPH
    (Groër et al 2010) """
    
    vrph_sol = []
    vrph_sol.append(max(sol)+1)
    visit_depot = False
    for node in sol:
        if node==0:
            visit_depot = True
        elif visit_depot:
            vrph_sol.append(-node)
            visit_depot = False
        else:
            vrph_sol.append(node)
    vrph_sol.append(0)
    return vrph_sol

def write_OPT_file(opt_file_path, D, sol):
    routes = [[0]+list(r)+[0] for x, r in groupby(sol, lambda z: z == 0) if not x]    
    with open(opt_file_path, 'w') as opt_file:
        for ri, route in enumerate(routes):
            opt_file.write("Route #%d: "%(ri+1))
            opt_file.write("\t".join( str(n) for n in route if n!=0))
            opt_file.write("\n")
        
        cost = sum(( D[sol[i-1],sol[i]] for i in range(1,len(sol))))
        if cost == int(cost):
            opt_file.write("Cost : %d\n"%int(cost))
        else:
            opt_file.write("Cost : %.2f\n"%cost)
            
    return opt_file_path

def write_TSPLIB_file(tsplib_file_path, D,
                      d=None, C=None, L=None, selected_idxs=None,
                      float_to_int_precision=None):
    
    if not selected_idxs:
        selected_idxs=list(range(len(D)))
    write_cvrp = False
    if tsplib_file_path[-4:].lower()==".vrp":
        write_cvrp = True
        
    with open(tsplib_file_path, 'w') as problem_file:
        problem_file.write("NAME: temporary\n")
        if write_cvrp:
            problem_file.write("TYPE: CVRP\n")
            if C:
                problem_file.write("CAPACITY: %d\n"%C)
            else:
                problem_file.write("CAPACITY: %d\n"%len(D))
            if L:
                problem_file.write("DISTANCE: %d\n"%L)
        else:    
            problem_file.write("TYPE: TSP\n")
        problem_file.write("COMMENT: temporary CVRP or TSP problem\n")
        problem_file.write("DIMENSION: %d\n" % len(selected_idxs))
        problem_file.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
        problem_file.write("EDGE_WEIGHT_FORMAT: UPPER_ROW\n")
        problem_file.write("EDGE_WEIGHT_SECTION\n")
        for ii, i in enumerate(selected_idxs):
            for j in selected_idxs[ii+1:]:
                if float_to_int_precision is not None:
                    problem_file.write(str(int(D[i,j]*float_to_int_precision)))
                else:
                    problem_file.write(str(D[i,j]))
                problem_file.write(" ")
            if ii!=len(selected_idxs)-1:
                problem_file.write("\n")
        if write_cvrp:
            problem_file.write("DEMAND_SECTION\n1 0\n")
            if d:
                for i in range(2,len(d)+1):
                    problem_file.write("%d %d\n"%(i, int(d[i-1])))
            else:
                for i in range(2,len(D)+1):
                    problem_file.write("%d 1\n"%i)
            problem_file.write("DEPOT_SECTION\n")
            problem_file.write("1\n")
            problem_file.write("-1\n")
        problem_file.write("EOF")
