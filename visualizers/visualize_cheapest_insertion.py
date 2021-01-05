# -*- coding: utf-8 -*-

import re
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE

MAKE_ANIM = True

node_re = re.compile("\sn([0-9]+)\s")
route_init_re = re.compile("Initialized a new route \#([0-9]+) (\[.*?\])")
route_select_re = re.compile("Inserting to route \#([0-9]+) (\[.*?\])")
#insert_attempt_re = "Try insertion ([0-9]+)-([0-9]+)-([0-9]+) with l_delta ([0-9]+\.[0-9]+)"
insert_re = re.compile("Chose to insert n([0-9]+) resulting in route (\[.*?\]) \(([0-9]+\.[0-9]+)\)")
route_complete_re =  re.compile("Route \#([0-9]+) finished as (\[.*?\]) \(([0-9]+\.[0-9]+)\)")
#constraint_re = re.compile("Insertion would break ([CL]) constraint.")

def _set_candidate_route_based_on_match(matcho, candidate_routes):
    route_id = int(matcho.group(1))
    candidate_route = eval(matcho.group(2))
    while len(candidate_routes)<route_id+1:
        candidate_routes.append([])
    candidate_routes[route_id] = candidate_route
    return route_id

route_id = None  
def _process_ci_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):
    global route_id
    
    newK = None
    changed = False
    
    # a route is selected for insertion
    rso = route_select_re.search(line)
    if rso:
        route_id = _set_candidate_route_based_on_match(rso, candidate_routes)

    # a new route is created
    rio = route_init_re.search(line)
    if rio:
        route_id = _set_candidate_route_based_on_match(rio, candidate_routes)
        changed = True
        labels.append("N/A")
        
    imo = insert_re.search(line)    
    if imo:
        inserted_node = int(imo.group(1))
        active_nodes[:] = [inserted_node]
        candidate_routes[route_id] = eval(imo.group(2))
        labels[-1] = imo.group(3)
        changed = True
    
    #cmo = constraint_re.search(line)    
    #elif cmo:
    #    # the candidate route was infeasible
    #    infeasible_routes[:] = [candidate_routes[-1]]
    #    candidate_routes[:] = []
    #    changed = True
    
    rco = route_complete_re.search(line) 
    if rco:
        complete_route = eval(rco.group(2))
        infeasible_routes[:] = []
        candidate_routes[:] = []
        active_nodes[:] = []
        complete_routes.append( complete_route )
        labels[-1] = rco.group(3)
        changed = True

    return changed, newK
    
if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("cheapest_insertion")
    visualize_procedure(algo_output, "cheapest_insertion", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_ci_debug_line)