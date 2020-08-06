# -*- coding: utf-8 -*-


# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from scipy.spatial.distance import pdist, squareform
from local_search import LSOPT
from local_search.intra_route_operators import do_2opt_move, do_3opt_move,\
                            do_relocate_move, do_exchange_move
from cvrp_io import generate_CVRP

# Helpers, so simple that they are sure to work right
from util import routes2sol, sol2routes
from util import objf

def set_weight(D,n1,n2,wt):
    """ A helper shorthand to set the symmetric distance matrix weights. This
    makes the D manipulation code esier to read and less error prone to write. """
    D[n1,n2]=wt
    D[n2,n1]=wt
   
def _normalise_route_order(sol):
    routes = sol2routes(sol)
    # make sure route nodes are directed so that smaller end point is first
    for route in routes:
        if route[1]>route[-2]:
            route.reverse()
    # make sure routes are from smallest first node first
    routes.sort()
    return routes2sol(routes)
    
class Test2Opt(unittest.TestCase):
 
    def setUp(self):
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(pts, "euclidean") )
        self.longMessage = True
 
    def test_empty_route(self):
        self.assertEqual(do_2opt_move([],self.D), (None, None))
        self.assertEqual(do_2opt_move([0,0],self.D), (None, None))
 
    def test_2opt_none_crossed(self):
        route = [0,1,2,3,0]
        sol, delta_f = do_2opt_move(route,self.D)
        self.assertEqual(sol, None, "Route was already 2-optimal, improvements are not possible")
       
    def test_2opt_one_crossed(self):
        route = [0,1,6,2,7,0]
        
        initial_f =objf(route,self.D)
        sol, delta_f = do_2opt_move(route,self.D, LSOPT.BEST_ACCEPT)
        do_2opt_f = objf(sol,self.D)
        self.assertEqual(sol ,[0,1,2,6,7,0],
         "chose invalid move, initial %f, optimized %f" % (initial_f,do_2opt_f))
        self.assertEqual(initial_f+delta_f,do_2opt_f, "The delta based and recalculated objective functions differ")
         
    
    def test_2opt_two_crossed(self):
        route = [0,1,6,3,4,5,2,7,0]
        required_moves = 2
        for i in range(required_moves ):
            #print "op %d:"%i+1
            #print "initial f =",f(route,self.D)
            opt_route, opt_f = do_2opt_move(route,self.D, LSOPT.BEST_ACCEPT)
            #print "2-opt'd f =",f(opt_route,self.D)
            if opt_route is None:
                return route
            else:
                route = opt_route
        self.assertEqual( i+1,required_moves )
        self.assertEqual( opt_route,[0,1,2,3,4,5,6,7,0] )

class Test3Opt(unittest.TestCase):
    def setUp(self):
        # a problem that is a circle with a radius of ~7, the depot on the rim
        self.pts = [(-4,4), #0, the depot
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

        self.D = squareform( pdist(self.pts, "euclidean") )  
        
        self.optimum = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0]
        self.segment1_head = [0, 1, 2]
        self.segment1_tail = [17, 18, 19, 0]
        self.segment1 = [17, 18, 19, 0, 1, 2]
        self.segment2 = [3, 4, 5, 6, 7, 8]
        self.segment3 = [9, 10, 11, 12, 13, 14, 15, 16]
        
        self.longMessage = True
        
        self.move_op = do_3opt_move

    def test_empty_route(self):
        self.assertEqual(self.move_op([],self.D), (None,None))
        self.assertEqual(self.move_op([0,0],self.D), (None,None))
        
    def test_3opt_already_optimal(self):
        initial_sol = list(range(20))
        #print "initial f =",f(initial_sol,self.D)
        sol, sol_f = self.move_op(initial_sol,self.D)
        #print "3-opt'd f =",f(initial_sol,self.D) if sol is None else objf(sol,self.D)
        self.assertEqual(sol, None )
    
    def test_3opt_move_s1_s2_r3(self):
        """ Test the 1st 3-opt alternative (actually corresponding a 2-opt
        move), where the last  segment (segment3) is reversed."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,8,9,100)
        set_weight(D,16,17,100)
        set_weight(D,8,16,1)
        set_weight(D,9,17,1)
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +
            self.segment2 +
            self.segment3[::-1] +
            self.segment1_tail) 
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
        

    def test_3opt_move_s1_r2_s3(self):
        """ Test the 2nd 3-opt alternative (actually corresponding a 2-opt
        move), where the middle segment (segment2) is reversed."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,8,9,100)
        set_weight(D,2,8,1)
        set_weight(D,3,9,1)
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol, D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +
            self.segment2[::-1] +
            self.segment3 +
            self.segment1_tail) 
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
            
    def test_3opt_move_s1_r3_r2(self):
        """ Test the 3rd 3-opt alternative (actually corresponding a 2-opt
        move), where the first segment (segment1) is reversed. However, in
        symmetric case this is equal to  reversing the order and direction of
        the segments 2 and 3 and this is the expected move operation here."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,16,17,100)
        set_weight(D,2,16,1)
        set_weight(D,3,17,1)   
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +
            self.segment3[::-1] +
            self.segment2[::-1]  +
            self.segment1_tail)
        self.assertEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
            
            
    def test_3opt_move_s1_r2_r3(self):
        """ Test the 4th 3-opt alternative, where the middle AND last segments
        (segment2 and segment3) are both reversed but their order kept."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,8,9,100)
        set_weight(D,16,17,100)
        set_weight(D,2,8,1)
        set_weight(D,3,16,1)   
        set_weight(D,9,17,1)   
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +
            self.segment2[::-1]  +
            self.segment3[::-1] +            
            self.segment1_tail)
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
           
            
    def test_3opt_move_s1_s3_s2(self):
        """ Test the 5th 3-opt alternative, where the middle AND last segments
        (segment2 and segment3) are swapped, but their traverse direction kept."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,8,9,100)
        set_weight(D,16,17,100)
        set_weight(D,2,9,1)
        set_weight(D,16,3,1)   
        set_weight(D,8,17,1)   
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +            
            self.segment3 +            
            self.segment2 +
            self.segment1_tail)
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
        
            
    def test_3opt_move_s1_s3_r2(self):
        """Test the 6th 3-opt alternative, where the middle AND last segments
        (segment2 and segment3) are swapped, and the segment2 is reversed."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,8,9,100)
        set_weight(D,16,17,100)
        set_weight(D,2,9,1)
        set_weight(D,16,8,1)   
        set_weight(D,3,17,1)   
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +            
            self.segment3 +            
            self.segment2[::-1] +
            self.segment1_tail)
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
        
            
    def test_3opt_move_s1_r3_s2(self):
        """Test the 7th and final 3-opt alternative, where the middle AND last segments
        (segment2 and segment3) are swapped, and the segment3 is reversed."""
        
        # modify the distance matrix to force a recombination move
        D = self.D.copy()
        set_weight(D,2,3,100)
        set_weight(D,8,9,100)
        set_weight(D,16,17,100)
        set_weight(D,2,16,1)
        set_weight(D,9,3,1)   
        set_weight(D,8,17,1)   
        
        initial_sol = self.optimum
        initial_sol_f = objf(initial_sol,D)
        sol, delta_f = self.move_op(initial_sol,D,
                                  strategy=LSOPT.BEST_ACCEPT)
        do_3opt_f = objf(sol,D)
        self.assertEqual(sol,
            self.segment1_head +            
            self.segment3[::-1] +            
            self.segment2 +
            self.segment1_tail) 
        self.assertAlmostEqual(initial_sol_f+delta_f, do_3opt_f, msg=\
          "The delta based and recalculated objective function values differ")
        
            
    def smoke_test_3opt(self):
        N, points, _, _, D, _, _ = generate_CVRP(50, 40, 10, 5)
        current_sol = list(range(N))
        while True:
            sol = self.move_op(current_sol ,D)
            if sol:
                current_sol = sol
            else:
                break
            
        
class TestRelocate(unittest.TestCase):
 
    def setUp(self):
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(pts, "euclidean") )
        self.longMessage = True
    
    def test_empty_route(self):
        self.assertEqual(do_relocate_move([],self.D), (None, None))
        self.assertEqual(do_relocate_move([0,0],self.D), (None, None))
 
    def test_no_improvements(self):
        route = [0,1,2,3,0]
        sol, delta_f = do_relocate_move(route,self.D)
        self.assertEqual(sol, None, "Route was already optimal, improvements are not possible")
       
    def test_one_move(self):
        route = [0,1,3,4,2,0]
        
        initial_f = objf(route,self.D)
        sol, delta_f = do_relocate_move(route,self.D, LSOPT.BEST_ACCEPT)
        do_1pm_f = objf(sol,self.D)
        self.assertEqual(sol ,[0,1,2,3,4,0],
         "chose invalid move, initial %f, optimized %f" % (initial_f,do_1pm_f))
        self.assertEqual(initial_f+delta_f,do_1pm_f, "The delta based and recalculated objective functions differ")
         
    def test_many_moves(self):
        route = [0,1,6,3,4,5,2,7,0]
        required_moves = 2
        for i in range(required_moves ):
            opt_route, opt_f = do_relocate_move(route,self.D, LSOPT.BEST_ACCEPT)
            if opt_route is None:
                return route
            else:
                route = opt_route
        self.assertEqual( i+1,required_moves )
        self.assertEqual( opt_route,[0,1,2,3,4,5,6,7,0] )
        
        
class TestExchange(unittest.TestCase):
 
    def setUp(self):
        pts = [(0,0), #0
               (1,1), #1
               (1,2), #2
               (1,3), #3
               (0,4), #4
               (-2,3),#5
               (-2,2),#6
               (-2,1)]#7
        self.D = squareform( pdist(pts, "euclidean") )
        self.longMessage = True
    
    def test_empty_route(self):
        self.assertEqual(do_exchange_move([],self.D), (None, None))
        self.assertEqual(do_exchange_move([0,0],self.D), (None, None))
 
    def test_no_improvements(self):
        route = [0,1,2,3,0]
        sol, delta_f = do_relocate_move(route,self.D)
        self.assertEqual(sol, None, "Route was already optimal, improvements are not possible")
       
    def test_one_move(self):
        route = [0,4,2,3,1,0]
        
        initial_f = objf(route,self.D)
        sol, delta_f = do_exchange_move(route,self.D, LSOPT.BEST_ACCEPT)
        do_2pm_f = objf(sol,self.D)
        
        sol = _normalise_route_order(sol)
        self.assertEqual(sol ,[0,1,2,3,4,0],
         "chose invalid move, initial %f, optimized %f" % (initial_f,do_2pm_f))
        self.assertAlmostEqual(initial_f+delta_f,do_2pm_f, msg="The delta based and recalculated objective functions differ")
         
    def test_many_moves(self):
        route = [0,1,6,3,4,5,2,7,0]
        required_moves = 2
        for i in range(required_moves ):
            opt_route, opt_f = do_exchange_move(route,self.D, LSOPT.BEST_ACCEPT)
            if opt_route is None:
                return route
            else:
                route = opt_route
        self.assertEqual( i+1,required_moves )
        self.assertEqual( opt_route,[0,1,2,3,4,5,6,7,0] )
        
    def test_a_tricky_move(self):
        route = [0, 7, 6, 4, 1, 2, 3, 5, 0]
        
        initial_f = objf(route,self.D)
        sol, delta_f = do_exchange_move(route,self.D, LSOPT.BEST_ACCEPT)
        do_2pm_f = objf(sol,self.D)
        
        sol = _normalise_route_order(sol)
        self.assertEqual(sol, [0, 4, 3, 2, 1, 5, 6, 7, 0],
         "chose invalid move, initial %f, optimized %f" % (initial_f,do_2pm_f))
        self.assertAlmostEqual(initial_f+delta_f,do_2pm_f, msg="The delta based and recalculated objective functions differ")
        
    
if __name__ == '__main__':
    unittest.main()
