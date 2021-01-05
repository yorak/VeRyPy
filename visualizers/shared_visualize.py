# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import subprocess
from export_dot import BB_HALF_SIZE, output_dot
from os import makedirs, path
from shutil import copyfile
from sys import stderr, argv, exit, version_info

GRAPHVIZ_NEATO_EXE = r"C:\Program Files (x86)\Graphviz\bin\neato.exe"
GIFSICLE_EXE = r"C:\MyTemp\LocalApplications\gifsicle-1.88-win64\gifsicle.exe"
PYTHON_EXE = r"C:\Anaconda2\python.exe"

class VISUALIZE:
    # flags that can be given to visualize_procedure. e.g. EMPTY|SOLUTION
    EMPTY = 1
    SEARCH = 2
    TSP = 4
    SOLUTION = 8
    ALL = 15

START_STATE_ANIM_FRAMES = 3
INTERMEDIATE_STATE_ANIM_FRAMES = 3
END_STATE_ANIM_FRAMES = 4
MAX_STEPS = 100

def visualize_procedure(algo_name, selector=VISUALIZE.ALL, make_anim=True, 
              process_debug_line_callback = None):
                  
    file_path = argv[-1]
    file_ext = file_path.lower()[-4:]
    if not (path.isfile(file_path) and (file_ext==".vrp"  or file_ext==".txt")):
        print(r"usage: visualize_X <path/to/problem_file.vrp|/path/to/output_file.txt>", file=stderr)
        print(r" outputs to EXISTING folder ./output/X/, where", file=stderr)
        print(r" X is the name of the algorithm (e.g. sweep).", file=stderr)
        exit()
    
    # "-v", "3" set the verbosity to print EVERYTHING
    # "-1" print only one try
    
    if file_ext==".vrp":
        algo_path = path.join("..","classic_heuristics", "%s.py"%algo_name)
        subpargs = [PYTHON_EXE, algo_path, "-1", "-v", "3", argv[-1]]
        proc = subprocess.Popen(subpargs , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # In Python3 bytes go in and come out. Special handling is required.
        if version_info[0] >= 3:
            out = p.communicate()[0].decode('ascii')
        else:
            out = p.communicate()[0]

        #print(out)
        #print(" ".join(subpargs))
    else:
        with open(file_path, 'r' ) as precalculated_output_file:
            out = precalculated_output_file.read()
    
        
    output_folder = "output/%s"%algo_name
    if not path.exists(output_folder):
        makedirs(output_folder) 
            
    points = None
    sol = None
    final_K = 1
    
    additional_points = []
    
    # read the problem specifiers (or load from file)
    for line in out.split("\n"):
        if "POINTS:" in line:
            points = eval(line.split(":")[1])
            max(zip(points))
        if "SOLUTION:" in line:
            #print(line)
            sol = eval(line.split(":")[1])
        
        # sometimes some special points may be outlside 
        if "points = " in line:
            point_coordinates_in_a_string = line.split("=")[-1]
            additional_points.extend( list( eval(point_coordinates_in_a_string ) ))
            
            
    if not points:
        print(r"ERROR: the vrp file does not have coordinates!", file=stderr)
    
    dot_files = []
   
    nparams = get_normalization_params(points+additional_points, (-BB_HALF_SIZE,BB_HALF_SIZE), keep_aspect_ratio=True)
    npts = normalize_to_rect(points, nparams)
    
    if selector&VISUALIZE.EMPTY:
        # write the initial state
        print("write step 0 (initial state)")
        with open(output_folder+"/i00_0000_gapvrp_initial.dot", "w") as initial_fh:
            output_dot(initial_fh, npts, rays=[])
            dot_files.append(initial_fh.name)
    
    step = 1
    trial = 0        
    K = 1
    active_nodes = []
    rays = []
    active_ray_idxs = []
    points_of_interest = []
    
    candidate_routes = []
    infeasible_routes = []
    complete_routes = []
    
    labels = []
    
    if process_debug_line_callback:  
        for line in out.split("\n"):
            print(line)
            newK = None
            if "DEBUG:" in line and process_debug_line_callback:
                changed, newK = process_debug_line_callback(
                     #inputs
                     line, nparams, K,
                     #outpupts
                     rays, active_nodes, active_ray_idxs, points_of_interest,
                     candidate_routes, infeasible_routes, complete_routes,
                     labels)
                
                if changed and selector&VISUALIZE.SEARCH:
                    if K is not None:
                        state_steps = 1
                        if newK:
                            state_steps = INTERMEDIATE_STATE_ANIM_FRAMES
                        for i in range(state_steps):
                            print("write step %d"%step)
                            dfn = output_folder+"/k%02d"%K+\
                                  ("" if trial==0 else ("_t%04d"%trial))+\
                                  "_%06d_gapvrp_step.dot"%step
                            
                            with open(dfn, "w") as step_fh:
                                output_dot(step_fh, npts, active_point_idxs=active_nodes,
                                       rays=rays, active_ray_idxs=active_ray_idxs,
                                       points_of_interest=points_of_interest,
                                       gray_routes=candidate_routes,
                                       red_routes=infeasible_routes,
                                       black_routes=complete_routes,
                                       label=", ".join(labels))
                                dot_files.append(step_fh.name)
                            step+=1
                if newK is not None:
                    if K==newK:
                        trial+=1
                    K = newK
                    final_K = newK
                    step = 1
                    active_nodes = []
                    rays = []
                    active_ray_idxs = []
                    points_of_interest = []
                    candidate_routes = []
                    infeasible_routes = []
                    complete_routes = []
            
            if step>MAX_STEPS:
                print("WARNING: Maximum number of steps (%d) for the visualization reached."%MAX_STEPS)
                break
                        
                        
    
    ## generate graphs for TSP sols
    if sol:
        from itertools import groupby
        routes = [[0]+list(group)+[0] for k, group in groupby(sol, lambda x:x==0) if not k]
        if selector&VISUALIZE.SOLUTION:
            print("write final solution")
            with open(output_folder+"/k%02d_%04d_gapvrp_final.dot"%(final_K,step+1), "w") as step_fh:
                output_dot(step_fh, npts, active_point_idxs=None,
                       rays=rays, active_ray_idxs=active_ray_idxs,
                       points_of_interest=points_of_interest,
                       black_routes=routes)
                dot_files.append(step_fh.name)
                step_fh.close()
    print("\n")
    
    # duplicate the start and end frames to view those a bit longer
    start_file = dot_files[0]
    end_file = dot_files[-1]
    
    if make_anim:
        if selector&VISUALIZE.EMPTY:
            for i in range(1,START_STATE_ANIM_FRAMES+1):
                start_copy = start_file.replace("0000", "%04d"%i)
                print("copy", start_file, start_copy)
                if start_file!=start_copy:
                    copyfile(start_file, start_copy)
                    dot_files.insert(0, start_copy)
    
        for i in range(1,END_STATE_ANIM_FRAMES+1):
            print(end_file)
            last_frame_number = int(end_file.split("_")[-3])
            end_copy = end_file.replace(
                "%04d"%last_frame_number,
                "%04d"%(last_frame_number+i)).replace("step","sol")
            print("copy", end_file, end_copy)
            if end_file!=end_copy:
                copyfile(end_file, end_copy)
                dot_files.append(end_copy)
    print("\n")
    
    ## convert dot to png
    for step, dfn in enumerate(dot_files):
        print("convert plot %d to gif "%step)
        ofn = dfn.replace(".dot", ".gif")
        subprocess.call([GRAPHVIZ_NEATO_EXE, "-Tgif", "-o",  path.abspath(ofn), path.abspath(dfn)])
    print("\n")        
    
    if make_anim:
        ## to gif anim
        print
        if not path.exists(output_folder+"/anims"):
            makedirs(output_folder+"/anims") 
        with open(output_folder+"/anims/%s_anim.gif"%algo_name, "w") as animf:
            print("convert gifs to anim")
            subprocess.call([GIFSICLE_EXE, "--delay=33", "--loop",  output_folder+"/*.gif"], stdout=animf)
            
    print("\nall done, ktxhbye")
    
def get_normalization_params(pts, to_range, keep_aspect_ratio=False):
    xl, yl  = zip(*pts)
    minx = min(xl)
    maxx = max(xl)
    rangex = maxx-minx
    miny = min(yl)
    maxy = max(yl)
    rangey = maxy-miny
    minr = min(to_range)
    maxr = max(to_range)
    ranger = maxr-minr
    
    if keep_aspect_ratio:
        rangex = max(rangex, rangey)
        rangey = rangex
    
    return (minx,miny,rangex,rangey,minr,ranger)
    
def normalize_to_rect(pts, normalization_params):
    """ pts elements are (x,y,d), where d is demand
    to_range is a 2-tuple for the desired output range min and max"""
    minx,miny,rangex,rangey,minr,ranger = normalization_params
    xl, yl  = zip(*pts)
    
    new_xl = [(x-minx)/float(rangex)*(ranger)+minr for x in xl]
    new_yl = [(y-miny)/float(rangey)*(ranger)+minr for y in yl]
    
    return zip(new_xl, new_yl)    
    