# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division


from util import objf
from local_search.intra_route_operators import do_2opt_move, do_3opt_move
from random import shuffle
    
def solve_tsp_2opt(D, selected_idxs):
    return solve_tsp_ropt(D, selected_idxs,
               do_shuffle=False, do2opt=True, do3opt=False)

def solve_tsp_3opt(D, selected_idxs):
    return solve_tsp_ropt(D, selected_idxs,
               do_shuffle=False, do2opt=False, do3opt=True)
    
    
def solve_tsp_ropt(D, selected_idxs,
                   do_shuffle=False, do2opt=True, do3opt=True):    
    # r-Opt (r \in {2,3} )
    endp = selected_idxs[0]
    
    if do_shuffle:
        shuffled_idxs = list(selected_idxs[1:])
        shuffle(shuffled_idxs)
        new_route = [endp]+shuffled_idxs+[endp]
    elif selected_idxs[-1]!=endp:
        new_route = selected_idxs+[endp]
    else:
        new_route = selected_idxs
        
    new_route_cost = objf(new_route, D)
    
    # make first 2-optimal
    if do2opt:
        improved = True
        while improved:
            improved = False
            improved_route, delta = do_2opt_move(new_route, D, 1)
            if improved_route is not None:
                new_route = improved_route 
                new_route_cost+=delta
                improved = True
    
    # then 3-optimal (do not waste time on "easy" 2-opt
    #  operations if the route has already been made 2-optimal
    if do3opt:
        improved = True
        while improved:
            improved = False
            improved_route, delta = do_3opt_move(new_route, D, 1)
            if improved_route is not None:
                new_route = improved_route 
                new_route_cost+=delta
                improved = True
            
    return new_route, new_route_cost

if __name__=="__main__":
    from shared_cli import tsp_cli
    tsp_cli("ropt", solve_tsp_ropt)
    