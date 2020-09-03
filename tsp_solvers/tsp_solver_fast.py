# -*- coding: utf-8 -*-
###############################################################################
""" Tries to be as fast as possible in solving TSPs as possible using Gurobi
and LKH as experimentally it was found out to be time-optimal."""
###############################################################################

from tsp_solvers.tsp_solver_ropt import solve_tsp_ropt
from tsp_solvers.tsp_solver_lkh import solve_tsp_lkh
from tsp_solvers.tsp_solver_gurobi import solve_tsp_gurobi

def solve_tsp_fast(D, selected_idxs):
    if len(selected_idxs)<5:
        return solve_tsp_ropt(D, selected_idxs)
    elif len(selected_idxs)<20:
        return solve_tsp_gurobi(D, selected_idxs)
    else:
        return solve_tsp_lkh(D, selected_idxs, num_runs=1)
    
if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("fast", solve_tsp_fast)