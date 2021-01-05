# -*- coding: utf-8 -*-

import re
from shared_visualize import visualize_cli, visualize_procedure, VISUALIZE
from util import sol2routes

MAKE_ANIM = True

initial_re = re.compile("Start from initial solution (\[.*?\]) \(([0-9]+\.[0-9]+)\), and with l1=([0-9]+\.[0-9]+), l2=([0-9]+\.[0-9]+)")
#move_start_re = re.compile("Finding a LR3OPT move for (\[.*?\]) \(([0-9]+\.[0-9]+)\)")
move_done_re = re.compile("Found improving LR3OPT move leading to (\[.*?\]) \(([0-9]+\.[0-9]+)\)")
remain_infeasible_re = re.compile("However, routes (\[.*?\]) remain infeasible.")
feasible_found_re = re.compile("Reached feasible solution (\[.*?\]) \(([0-9]+\.[0-9]+)\)")
inc_lambda_re =  re.compile("No improving moves left, increasing lambda to l1=([0-9]+\.[0-9]+), l2=([0-9]+\.[0-9]+)")

def _process_lr3opt_debug_line(line, normalization_parameters, currentK,
                       #output here
                       rays, active_nodes, active_ray_idxs,
                       points_of_interest,
                       candidate_routes, infeasible_routes, complete_routes,
                       labels):    
    newK = None
    changed = False
    
    ism = initial_re.search(line)
    #msm = move_start_re.search(line)
    mdm = move_done_re.search(line)
    rim = remain_infeasible_re.search(line)
    ffm = feasible_found_re.search(line)
    ilm = inc_lambda_re.search(line)
    
    if ism:
        sol = eval(ism.group(1))
        infeasible_routes[:] = sol2routes(sol)
        labels[:] = ["l1_C=%s, l1_L=%s"%(ism.group(3), ism.group(4))]
        changed = True
   
    if mdm:
        sol = eval(mdm.group(1))
        candidate_routes[:] = sol2routes(sol)
        changed = False
   
    if rim:
        infeasible_routes[:] = eval(rim.group(1))
        print("REMAIN INFEASIBLE", infeasible_routes[:])
        candidate_routes[:] = [r for r in candidate_routes if r not in infeasible_routes] 
        changed = True
   
    if ffm:
        infeasible_routes[:] = []
        sol = eval(ffm.group(1))
        complete_routes[:] = sol2routes(sol)
        changed = True
    
    if ilm:
        labels[:] = ["l1_C=%s, l1_L=%s"%(ilm.group(1), ilm.group(2))]
    
    return changed, newK
    
if __name__=="__main__":
    algo_output, problem_name, keep_files = visualize_cli("lr3opt")
    visualize_procedure(algo_output, "lr3opt", problem_name, selector=VISUALIZE.ALL,
              make_anim=MAKE_ANIM, keep_files=keep_files,
              process_debug_line_callback = _process_lr3opt_debug_line)