# -*- coding: utf-8 -*-

import re
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE

MAKE_ANIM = True

node_re = re.compile("\sn([0-9]+)\s")
#tsp_re = re.compile("Got TSP solution \((-?[0-9]+\.[0-9]+)\) = (\[.*?\])")
init_re = re.compile("Initialize route #([0-9]+) with n([0-9]+), resulting to (\[.*?\])")
complete_re = re.compile("storing route #([0-9]+) as (\[.*?\])")
insert_re = re.compile("Inserted n([0-9]+) to create a route (\[.*?\]).")
parallel_re = re.compile("Parallel route building with seeds (\[.*?\])")

def _process_nn_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):
    newK = None
    changed = False
    
    nmo = init_re.search(line)
    imo = None #insert_re.search(line)
    cmo = None #constraint_re.search(line)
    rmo = None #complete_re.search(line)
    pmo = None #parallel_re.search(line)
    raise NotImplementedError # TODO: the re's must match to the nn output!

    if nmo:
        candidate_route = eval(nmo.group(3))
        candidate_routes.append( candidate_route )
        changed = True
    elif imo:
        #inserted_node = int(imo.group(1))
        candidate_route = eval(imo.group(2))
        candidate_routes[:] = [candidate_route+[0]]
        changed = True
    elif cmo:
        # the candidate route was infeasible
        infeasible_routes[:] = [candidate_routes[-1]]
        candidate_routes[:] = []
        changed = True
    elif rmo:
        complete_route = eval(rmo.group(1))
        infeasible_routes[:] = []
        candidate_routes[:] = []
        active_nodes[:] = []
        complete_routes.append( complete_route )
        changed = True        
    elif pmo:
        seeds = eval(pmo.group(1))
        newK = len(seeds)
        changed = True
    #elif tmo:
    #    prev_tsp_sol = eval(tmo.group(2))
    elif "Initialize route" in line:
        mo = node_re.search(line)
        seed_node = int(mo.group(1))
        active_nodes[:] = [seed_node]
        candidate_routes.append( [0, seed_node, 0] )
        changed = True
    elif "Check feasibility of inserting" in line:
        mo = node_re.search(line)
        candidate_routes[-1].insert(-1, int(mo.group(1)))
        changed = False
    elif "Sequential route bulding" in line:
        newK=0
        
    #print "VGV:", line,
    #print "VGV:", changed, newK
    return changed, newK
    
if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("nearest_neighbor")
    visualize_procedure(algo_output, "nearest_neighbor", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_nn_debug_line)