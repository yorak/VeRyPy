# -*- coding: utf-8 -*-

import re
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE

MAKE_ANIM = True

node_re = re.compile("node\s([0-9]+),?\s")
seeds_re = re.compile("getting\s([0-9]+)\sseed")
tsp_re = re.compile("Got TSP solution \((-?[0-9]+\.[0-9]+)\) = (\[.*?\])")
select_seed_re = re.compile("Selecting a node n([0-9]+) (\([0-9]+\.[0-9]+, [0-9]+\.[0-9]+\)) that is .*? point to be a seed")
seed_candidates_re = re.compile("[Ss]elect [0-9]+ seed nodes from .+? nodes = (\[.*?\])")

#TODO: 
# - does not show infeasibility violations (when L constraint is enabled)
# - shows TSP solutions at the end, even though they are calculated as soon as
#    a cone gets ready. For problems with L constraint the visualization does
#    not illustrate the reason for the incoplete sweep (the TSP sol violated 
#    the L constraint)

def _process_gap_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):
    newK = None
    changed = False
    
    se = select_seed_re.search(line)
    sc = seed_candidates_re.search(line)
    
    mo = seeds_re.search(line)
    if mo:
        newK = int(mo.group(1))
        changed = True
    elif sc:
        new_seed_candidates = eval(sc.group(1))
        active_nodes[:] = new_seed_candidates
        changed = True
    elif se:
        new_seed_node = int(se.group(1))
        active_nodes.remove(new_seed_node)
        node_coord = normalize_to_rect( [eval(se.group(2))], normalization_parameters)[0]       
        points_of_interest.append( node_coord )
        changed = True
    elif "group_start_ray" in line:
        start_ray = float(line.split("=")[1])
        active_ray_idxs.append(len(rays))
        rays.append(start_ray)
        changed = True
    elif "cone grows group to ray" in line:
        inc_ray = float(line.split("=")[1])
        rays.append(inc_ray)
        changed = True
    elif  "...got seed points"  in line:
        new_seed_points = eval(line.split("points")[1])
        new_seed_points = normalize_to_rect(new_seed_points, normalization_parameters)
        # empty the list and add in all new seed points
        points_of_interest[:] = []
        points_of_interest.extend(new_seed_points)
        changed = True

    mo = node_re.search(line)
    if mo:
        node = int(mo.group(1))
        if not active_nodes or active_nodes[0]!=node:
            # empty the list, but keep the list object and add the node as 
            #  the single element
            active_nodes[:] = []
            active_nodes.append( node)
            changed = True
            
    tsp_match = tsp_re.search(line)
    if tsp_match:
        route = eval(tsp_match.group(2))
        complete_routes.append(route)
        changed = True
    
    #print "VGV:", line,
    #print "VGV:", changed, newK
    return changed, newK
    
if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("gapvrp")
    visualize_procedure(algo_output, "gapvrp", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_gap_debug_line)