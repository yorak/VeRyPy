# -*- coding: utf-8 -*-
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE
from visualize_sweep import _process_debug_line as _process_sweep_debug_line
import re

MAKE_ANIM = True

node_re = re.compile("n([0-9]+)")
KII_node_re = re.compile("KII=n([0-9]+)")
JJX_node_re = re.compile("JJX=n([0-9]+)")
JII_node_re = re.compile("JII=n([0-9]+)")
route_re = re.compile("(\[[0-9, ]+\])")
              
def _process_gms_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):
                           
    changed, newK = _process_sweep_debug_line(line, normalization_parameters,
                       currentK, rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels)

    if "Trying to replace" in line:
        # active node is the one that will be removed
        changed = True
        active_nodes[:] = []
        
    KII_match = KII_node_re.search(line)    
    KII = int(KII_match.group(1)) if KII_match else None

    JJX_match = JJX_node_re.search(line)    
    JJX = int(JJX_match.group(1)) if JJX_match else None

    JII_match = JII_node_re.search(line)    
    JII = int(JII_match.group(1)) if JII_match else None
    
    if KII is not None:
        active_nodes.append(KII)
    if JJX is not None:
        active_nodes.append(JJX)
    if JII is not None:
        active_nodes.append(JII)
        
    # not improving move
    if "which forms a route" in line or\
       "would have formed a route" in line:
        route_match = route_re.search(line)
        
        if candidate_routes:
            del candidate_routes[-1]
        
        if route_match:
            route = eval(route_match.group(1))
            
        is_infeasible = "would have formed a route" in line
        
        if is_infeasible:
            infeasible_routes[:] = [route]
        else:
            candidate_routes[:] = [route]
        changed = True
        
        print("Active nodes on reject/accept", active_nodes)

    if len(candidate_routes)>1:
        candidate_routes[:] = [candidate_routes[-1]]
    
    return changed, newK

if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("gillet_miller_sweep")
    visualize_procedure(algo_output, "gillet_miller_sweep", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_gms_debug_line)