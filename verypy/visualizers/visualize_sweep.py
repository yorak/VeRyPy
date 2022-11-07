# -*- coding: utf-8 -*-

import re
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE
from collections import OrderedDict
from math import pi 

MAKE_ANIM = True

included_re = re.compile("Included n([0-9]+) to route set ([0-9]+)")
rays_re = re.compile("Considering n([0-9]+) between rays (-?[0-9]+\.[0-9]+), (-?[0-9]+\.[0-9]+)")
added_re = re.compile("Added n([0-9]+) to the route set")
tsp_re = re.compile("Got TSP solution (\[.*?\]) \((-?[0-9]+\.[0-9]+)\)")
route_set_re = re.compile("Route set .*?(\[.*?\]) full")

L_violation_re = re.compile("L constraint was violated, removed n([0-9]+) from the route set")

# global state
ray_cache = OrderedDict()

def _process_sweep_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):                           
                         
    global ray_cache 
        
    if "Do a sweep from" in line:
        # start a new sweep
        ray_cache = OrderedDict()
        return True, 0
    
    rays_match = rays_re.search(line)
    if rays_match:
        node = float( rays_match.group(1) )
        l = float( rays_match.group(2) )
        r = float( rays_match.group(3) )
        ray_cache[node] = r
            
        if len(rays)==0:
            if (l>r):
                l-=2*pi
            rays.append(l)
            active_ray_idxs.append(0)
            # first ray, show
            return True, None
        else:
            # do not update, because this node might not belong to the current group
            return False, None
    
    added_match = added_re.search(line)
    if added_match:
        node = int( added_match.group(1) )
        
        active_nodes[:] = []
        active_nodes.append(node)
    
        if active_ray_idxs[-1]==len(rays)-1:
            rays.append(ray_cache[node])
        else:
            rays[-1] = ray_cache[node]
        
        return True, None
    
    rs_match = route_set_re.search(line)
    if rs_match:
        active_nodes[:] = eval( rs_match.group(1) )
        active_ray_idxs.append(len(rays)-1)
        return True, None

    tsp_match = tsp_re.search(line)
    if tsp_match:
        #rays[-1] = ray_cache[ray_cache.keys()[-1]]
        active_nodes[:] = []
        infeasible_routes[:] = []
        route = eval(tsp_match.group(1))
        candidate_routes.append(route)
        return True, None

    L_violation_match = L_violation_re.search(line)
    if L_violation_match and len(candidate_routes)>0:
        node_to_remove = int(L_violation_match.group(1))
        active_nodes[:] = [node_to_remove]
        violating_route = candidate_routes[-1]
        del candidate_routes[-1]
        infeasible_routes.append(violating_route)
        
        # sweeps back
        ray_to_prev = False
        for node, ray in reversed(ray_cache.items()):
            if ray_to_prev:
                rays[-1] = ray
                break
            if node == node_to_remove:
                ray_to_prev = True
                
        return True, None
        
    if "Route completed" in line: 
        complete_route = eval( line[line.find('['):] )
        complete_routes.append( complete_route )
        candidate_routes[:] = []
        infeasible_routes[:] = []
        active_ray_idxs.append( len(rays)-1 )
        return True, None
    
    #print "VSW:", line,
    #print "VSW:", changed, newK
    return False, None
    
if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("sweep")
    visualize_procedure(algo_output, "sweep", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_sweep_debug_line)