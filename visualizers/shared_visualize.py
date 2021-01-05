# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import subprocess
import argparse
from export_dot import BB_HALF_SIZE, output_dot
from os import makedirs, path
from shutil import copyfile
from sys import stderr, argv, exit, version_info

# Exoect these to be found from the PATH, edit if not
GRAPHVIZ_NEATO_EXE = r"neato"
GIFSICLE_EXE = r"gifsicle"
PYTHON_EXE = r"python"

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

def visualize_cli(algo_name):

    parser = argparse.ArgumentParser(description="Produces animated gif to a folder ./output")
    parser.add_argument('file_path', help="/path/to/problem_file.vrp OR /path/to/algo_raw_output_file.txt")
    parser.add_argument('-k', dest='keep_files', help="Keep intermediate dot and gif files.", action='store_true')
    args = parser.parse_args()

    file_base, file_ext = path.splitext(args.file_path)
    file_ext = file_ext.lower()
    file_base = path.basename(file_base)

    if not (path.isfile(args.file_path) or (file_ext!=".vrp" and file_ext!=".txt")):
        parser.print_help(stderr)
        exit()

    # "-v", "3" set the verbosity to print EVERYTHING
    # "-1" print only one try
    
    if file_ext==".vrp":
        pathname = path.dirname(argv[0])        
        algo_path = path.abspath(path.join(pathname, "..", "classic_heuristics", "%s.py"%algo_name))
        subpargs = [PYTHON_EXE, algo_path, "-1", "-v", "3", args.file_path]
        #print("CALL", " ".join(subpargs))
        p = subprocess.Popen(subpargs , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # In Python3 bytes go in and come out. Special handling is required.
        if version_info[0] >= 3:
            out = p.communicate()[0].decode('ascii')
        else:
            out = p.communicate()[0]
        
        #print(out)
        #print(" ".join(subpargs))
    else:
        with open(args.file_path, 'r' ) as precalculated_output_file:
            out = precalculated_output_file.read()
    
    return (out, file_base, args.keep_files)

def visualize_procedure(algo_output, algo_name, data_name, selector=VISUALIZE.ALL,
              make_anim=True, keep_files=False,
              process_debug_line_callback = None):

    output_folder = path.join("output",algo_name)
    if not path.exists(output_folder):
        makedirs(output_folder) 
            
    points = None
    sol = None
    final_K = 1
    
    additional_points = []
    
    # read the problem specifiers (or load from file)
    for line in algo_output.split("\n"):
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
    gif_files = []
    
    nparams = get_normalization_params(points+additional_points, (-BB_HALF_SIZE,BB_HALF_SIZE), keep_aspect_ratio=True)
    scaled_points = list(normalize_to_rect(points, nparams))
    
    if selector&VISUALIZE.EMPTY:
        # write the initial state
        print("write step 0 (initial state)")
        with open(path.join(output_folder, "i00_000000_gapvrp_initial.dot"), "w") as initial_fh:
            output_dot(initial_fh, scaled_points, rays=[], label="START")
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
        for line in algo_output.split("\n"):
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
                            dfn = ("k%02d"%K+
                                   ("" if trial==0 else ("_t%04d"%trial))+
                                   "_%06d_gapvrp_step.dot"%step)
                            dfnp = path.join(output_folder, dfn)
                                            
                            
                            with open(dfnp, "w") as step_fh:
                                output_dot(step_fh, scaled_points, active_point_idxs=active_nodes,
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
            dotfp = path.join(output_folder,"k%02d_%06d_gapvrp_final.dot"%(final_K,step+1))
            with open(dotfp, "w") as step_fh:
                output_dot(step_fh, scaled_points, active_point_idxs=None,
                       rays=rays, active_ray_idxs=active_ray_idxs,
                       points_of_interest=points_of_interest,
                       black_routes=routes, label="STOP")
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
        gif_files.append(ofn)
    print("\n")        
    
    if make_anim:
        ## to gif anim
        print()
        anims_dir = "output"
        if not path.exists(anims_dir):
            makedirs(anims_dir) 
        
        gif_input_files_path = path.join(output_folder, "*.gif") # all outputted gif files. 
        anim_output_file_path = path.join(anims_dir, "%s_%s_anim.gif"%(algo_name, data_name))
        print("Converting gifs to anim ...", end="")
        call_gifsicle = [GIFSICLE_EXE, "--delay=33", "--loop", 
                         path.abspath(gif_input_files_path),
                         "-o", path.abspath(anim_output_file_path)]
        #print(" ".join(call_gifsicle))
        # Call with shell=Trye to use shell wildcard expands.
        if subprocess.call(" ".join(call_gifsicle), shell=True)==0:
            print("DONE!")
            print("Final gif animation written to file", anim_output_file_path)
        else:
            print("FAIL!")
        
        if not keep_files:
            removed_cnt = 0
            for tempf in dot_files+gif_files:
                if path.isfile(tempf):
                    try:
                        os.remove(tempf)
                    except:
                        pass
                    finally:
                        removed_cnt+=1
            print("Cleaned up %d of %d temp gif/dot files."%
                  (removed_cnt,len(dot_files)+len(gif_files)))

            
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
    