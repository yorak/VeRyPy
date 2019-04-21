# -*- coding: utf-8 -*-

from scipy.spatial.distance import pdist, squareform
from local_search.intra_route_operators import do_3opt_move

# this is just to get pretty pictures out of 3-opt moves
SHOW_DEBUG_GRAPHVIZ_VISUALIZATIONS = True

import sys
import itertools
from shutil import copyfile
sys.path.append(r'C:\Users\juherask\Dissertation\Phases\Features\extractors\RKM16\tools')
from vrpopt2dot import write_dot_to_stream, run_graphviz_and_show_image

def _set_weight(D,n1,n2,wt):
    """ A helper shorthand to set the symmetric distance matrix weights. This
    makes the D manipulation code esier to read and less error prone to write. """
    D[n1,n2]=wt
    D[n2,n1]=wt
    
def _visualize_result(fname, sol, pts, D):
    routes = [[0]+list(r)+[0] for x, r in itertools.groupby(sol, lambda z: z == 0) if not x]
    dotfile = open(fname, "w")
    
    write_dot_to_stream(dotfile, pts, routes=routes, D=D)
    dotfile.close()
    if SHOW_DEBUG_GRAPHVIZ_VISUALIZATIONS:
        tmpfile = run_graphviz_and_show_image(dotfile.name)
        copyfile(tmpfile, fname.replace(".dot","")+".png")

def visualize_3opt_move_s1_s2_r3(D, pts, optimum):
    """ Test the 1st 3-opt alternative (actually corresponding a 2-opt
    move), where the last  segment (segment3) is reversed."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,8,9,100)
    _set_weight(D,16,17,100)
    _set_weight(D,8,16,1)
    _set_weight(D,9,17,1)
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)

    _visualize_result("output/3opt_move1_s1_s2_r3.dot", sol, pts, D)
        
def visualize_3opt_move_s1_r2_s3(D, pts, optimum):
    """ Test the 2nd 3-opt alternative (actually corresponding a 2-opt
    move), where the middle segment (segment2) is reversed."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,8,9,100)
    _set_weight(D,2,8,1)
    _set_weight(D,3,9,1)
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move2_s1_r2_s3.dot", sol, pts, D)
            
def visualize_3opt_move_s1_r3_r2(D, pts, optimum):
    """ Test the 3rd 3-opt alternative (actually corresponding a 2-opt
    move), where the first segment (segment1) is reversed. However, in
    symmetric case this is equal to  reversing the order and direction of
    the segments 2 and 3 and this is the expected move operation here."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,16,17,100)
    _set_weight(D,2,16,1)
    _set_weight(D,3,17,1)   
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move3_s1_r3_r2.dot", sol, pts, D)
            
            
def visualize_3opt_move_s1_r2_r3(D, pts, optimum):
    """ Test the 4th 3-opt alternative, where the middle AND last segments
    (segment2 and segment3) are both reversed but their order kept."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,8,9,100)
    _set_weight(D,16,17,100)
    _set_weight(D,2,8,1)
    _set_weight(D,3,16,1)   
    _set_weight(D,9,17,1)   
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move4_s1_r2_r3.dot", sol, pts, D)
        
        
def visualize_3opt_move_s1_s3_s2(D, pts, optimum):
    """ Test the 5th 3-opt alternative, where the middle AND last segments
    (segment2 and segment3) are swapped, but their traverse direction kept."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,8,9,100)
    _set_weight(D,16,17,100)
    _set_weight(D,2,9,1)
    _set_weight(D,16,3,1)   
    _set_weight(D,8,17,1)   
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move5_s1_s3_s2.dot", sol, pts, D)
        
def visualize_3opt_move_s1_s3_r2(D, pts, optimum):
    """Test the 6th 3-opt alternative, where the middle AND last segments
    (segment2 and segment3) are swapped, and the segment2 is reversed."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,8,9,100)
    _set_weight(D,16,17,100)
    _set_weight(D,2,9,1)
    _set_weight(D,16,8,1)   
    _set_weight(D,3,17,1)   
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move6_s1_s3_r2.dot", sol, pts, D)
        
def visualize_3opt_move_s1_r3_s2(D, pts, optimum):
    """Test the 7th and final 3-opt alternative, where the middle AND last segments
    (segment2 and segment3) are swapped, and the segment3 is reversed."""
    
    # modify the distance matrix to force a recombination move
    D = D.copy()
    _set_weight(D,2,3,100)
    _set_weight(D,8,9,100)
    _set_weight(D,16,17,100)
    _set_weight(D,2,16,1)
    _set_weight(D,9,3,1)   
    _set_weight(D,8,17,1)   
    
    initial_sol = optimum
    sol = do_3opt_move(initial_sol,D)
    
    _visualize_result("output/3opt_move7_s1_r3_s2.dot", sol, pts, D)
    
def main():
    """ specify a problem and plot all 3-opt moves """
    # a problem that is a circle with a radius of ~7, the depot on the rim
    pts = [(-4,4), #0, the depot
           (-3,5), #1
           (-1,6), #2
           
           (1,6), #3               
           (3,5), #4
           (4,4), #5
           (5,3), #6
           (6,1), #7               
           (6,-1), #8
           
           (5,-3), #9
           (4,-4), #10
           (3,-5), #11
           (1,-6), #12               
           (-1,-6), #13
           (-3,-5), #14
           (-4,-4), #15
           (-5,-3), #16
           
           (-6,-1), #17               
           (-6,1), #18
           (-5,3), #19
          ]
    
    D = squareform( pdist(pts, "euclidean") )  
    optimum = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0]

    visualize_3opt_move_s1_s2_r3(D, pts, optimum)
    visualize_3opt_move_s1_r2_s3(D, pts, optimum)
    visualize_3opt_move_s1_r3_r2(D, pts, optimum)
    visualize_3opt_move_s1_r2_r3(D, pts, optimum)
    visualize_3opt_move_s1_s3_s2(D, pts, optimum)
    visualize_3opt_move_s1_s3_r2(D, pts, optimum)
    visualize_3opt_move_s1_r3_s2(D, pts, optimum)
    
if __name__=="__main__":
    main()