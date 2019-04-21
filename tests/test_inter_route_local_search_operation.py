# -*- coding: utf-8 -*-


# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

import unittest
from scipy.spatial.distance import pdist, squareform
from local_search import LSOPT
from local_search.inter_route_operators import do_2optstar_move,\
                            do_insert_move,do_redistribute_move,\
                            do_2point_move, do_1point_move
                            
from routedata import RouteData

# Helpers, so simple that they are sure to work right
from util import routes2sol, sol2routes
from util import objf as route_l
from util import totald as route_d

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
    
class Test2OptStar(unittest.TestCase):
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
        self.d = [1.0]*len(self.D)
        self.d[0] = 0.0
        self.longMessage = True
    
    def _make_improving_move(self, r1, r2, C=None, d=None, L=None):
        D = self.D
        r1rd = RouteData(r1, route_l(r1,D), route_d(r1, d), None)
        r2rd = RouteData(r2, route_l(r2,D), route_d(r2, d), None)
        return do_2optstar_move(r1rd, r2rd, D, d=d, C=C, L=L,
                               strategy=LSOPT.BEST_ACCEPT)
        
    def test_improving_move(self):
        new_r1d, new_r2d, f_delta = self._make_improving_move(
                r1=[0,1,2,7,0], r2=[0,4,3,6,5,0], C=4.0, d = self.d) 
        new_r1d.normalize()
        new_r2d.normalize()
        
        #print(new_r1d[0], new_r2d[0], route_l(new_r1d[0], self.D)+route_l(new_r2d[0], self.D))
        #print([0,1,2,3,4,0], [0,5,6,7,0], route_l([0,1,2,3,4,0], self.D)+route_l([0,5,6,7,0], self.D))
        self.assertEqual( new_r1d[0], [0,1,2,3,4,0], "nodes 4 and 7 should be swapped")
        self.assertEqual( new_r2d[0], [0,5,6,7,0], "nodes 4 and 7 should be swapped")


class TestInsertion(unittest.TestCase):
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
        self.d = [1.0]*len(self.D)
        self.d[0] = 0.0
        self.longMessage = True
        
        
    def test_insert_from_empty(self):
        rd1 = RouteData()
        r2 = [0,1,3,4,0]
        rd2 = RouteData(r2, route_l(r2, self.D), route_d(r2, self.d), None)
        
        _, result_rd, result_delta = do_insert_move(rd1, rd2, self.D)
        self.assertEqual( result_rd.route, [0,1,3,4,0], "All should be inserted")
        
    def test_insert_to_empty(self):
        r1 = [0,1,3,4,0]
        rd1 = RouteData(r1, route_l(r1, self.D), route_d(r1, self.d), None)
        rd2 = RouteData()
        
        _, result_rd, result_delta = do_insert_move(rd1, rd2, self.D)
        self.assertEqual( result_rd.route, [0,1,3,4,0], "All should be inserted")
        
    def test_insert_from_empty_to_empty(self):
        _, result_rd, result_delta = do_insert_move(RouteData(), RouteData(), self.D)
        self.assertEqual( result_rd.route, [0,0], "Should remain empty")
        
    def _insert_one(self, strategy, result_target, msg_if_fail, C=None, L=None, to_insert=2):
        r1 = [0,1,3,4,0]
        rd1 = RouteData(r1, route_l(r1, self.D), route_d(r1, self.d), None)
        _, result_rd, result_delta = do_insert_move(to_insert, rd1, self.D,
                                    C=C, d=self.d, L=L, 
                                    strategy=strategy)
        result_route = None if (result_rd is None) else result_rd.route
        self.assertEqual( result_route, result_target, msg_if_fail)
        if result_delta is not None:
            self.assertAlmostEqual(rd1.cost+result_delta,result_rd.cost, msg="The delta based and stored objective functions differ")
            self.assertAlmostEqual(rd1.cost+result_delta,route_l(result_rd.route,self.D), msg="The delta based and recalculated objective functions differ")
            #print("one base",route_l(rd1.route,self.D))
            #print("one result",route_l(result_rd.route,self.D))

    
    def test_insert_one(self):
        self._insert_one(LSOPT.FIRST_ACCEPT,
                         [0,1,3,4,2,0],
                         "FIRST_ACCEPT should just append")
        self._insert_one(LSOPT.BEST_ACCEPT,
                         [0,1,2,3,4,0],
                         "BEST_ACCEPT should insert at the best position")

    def _insert_many(self, strategy, result_target, msg_if_fail, C=None, L=None):
        to_insert = [3,2]
        r1 = [0]+to_insert+[0]
        rd1 = RouteData(r1, route_l(r1, self.D), route_d(r1, self.d), None)
        r2 = [0,1,4,0]
        rd2 = RouteData(r2, route_l(r2, self.D), route_d(r2, self.d), None)
        
        _, list_input_result_rd, list_input_result_delta = do_insert_move(
                                        to_insert, rd2, self.D,
                                        C=C, d=self.d, L=L, 
                                        strategy=strategy)
        _, rd_input_result_rd, rd_input_result_delta = do_insert_move(
                                        rd1, rd2, self.D,
                                        C=C, d=self.d, L=L, 
                                        strategy=strategy)
        list_input_result_route = None if (list_input_result_rd is None) else  list_input_result_rd.route
        rd_input_result_route = None if (rd_input_result_rd is None) else  rd_input_result_rd.route
        self.assertEqual( list_input_result_route, result_target, msg_if_fail)
        self.assertEqual( list_input_result_route, rd_input_result_route,
                         "Inserting from a list or RouteData should lead to same result")
        if list_input_result_delta is not None:
            self.assertAlmostEqual(rd2.cost+list_input_result_delta,list_input_result_rd.cost, msg="The delta based and stored objective functions differ")
            self.assertAlmostEqual(rd2.cost+list_input_result_delta,route_l(list_input_result_rd.route,self.D), msg="The delta based and recalculated objective functions differ")
        if rd_input_result_delta is not None:
            self.assertAlmostEqual(rd2.cost+rd_input_result_delta,rd_input_result_rd.cost, msg="The delta based and stored objective functions differ")
            self.assertAlmostEqual(rd2.cost+rd_input_result_delta,route_l(rd_input_result_rd.route,self.D), msg="The delta based and recalculated objective functions differ")
            #print("multi base",route_l(rd2.route,self.D))
            #print("multi result",route_l(rd_input_result_rd.route,self.D))

    def test_insert_many(self):
        self._insert_many(LSOPT.FIRST_ACCEPT,
                          [0,1,4,3,2,0],
                          "FIRST_ACCEPT should just append")
        self._insert_many(LSOPT.BEST_ACCEPT,
                         [0,1,2,3,4,0],
                         "BEST_ACCEPT should insert at the best position")
    
    def test_insert_success_with_C_costraint(self):
        self._insert_many(LSOPT.FIRST_ACCEPT,
                          [0,1,4,3,2,0],
                          "both of the inserted should fit with FIRST_ACCEPT",
                          C=5.5)
        self._insert_many(LSOPT.BEST_ACCEPT,
                          [0,1,2,3,4,0],
                          "both of the inserted should fit with BEST_ACCEPT",
                          C=5.5)
        self._insert_one(LSOPT.FIRST_ACCEPT,
                         [0,1,3,4,2,0],
                         "the inserted should fit with FIRST_ACCEPT",
                         C=5.5)
        self._insert_one(LSOPT.BEST_ACCEPT,
                         [0,1,2,3,4,0],
                         "the inserted should fit with BEST_ACCEPT",
                         C=5.5)
    
    def test_insert_fail_due_to_C_costraint(self):
        self._insert_many(LSOPT.FIRST_ACCEPT,
                          None,
                          "both of the inserted should not fit with FIRST_ACCEPT",
                          C=3.5)
        self._insert_many(LSOPT.BEST_ACCEPT,
                          None,
                          "both of the inserted should not fit with BEST_ACCEPT",
                          C=3.5)
        self._insert_one(LSOPT.FIRST_ACCEPT,
                         None,
                         "the inserted should not fit with FIRST_ACCEPT",
                         C=3.5)
        self._insert_one(LSOPT.BEST_ACCEPT,
                         None,
                         "the inserted should not fit with BEST_ACCEPT",
                         C=3.5)
    
    def test_insert_success_with_L_costraint(self):
        self._insert_many(LSOPT.FIRST_ACCEPT,
                          [0, 2, 1, 3, 4, 0],
                          "both of the inserted should fit with FIRST_ACCEPT",
                          L=20)
        self._insert_many(LSOPT.BEST_ACCEPT,
                          [0, 1, 2, 3, 4, 0],
                          "both of the inserted should fit with BEST_ACCEPT",
                          L=20)
        self._insert_one(LSOPT.FIRST_ACCEPT,
                         [0, 5, 1, 3, 4, 0],
                         "the inserted should fit with FIRST_ACCEPT",
                         L=20,
                         to_insert=5)
        self._insert_one(LSOPT.BEST_ACCEPT,
                         [0, 1, 3, 4, 5, 0],
                         "the inserted should fit with BEST_ACCEPT",
                         L=20,
                         to_insert=5)
        
    def test_insert_fail_due_to_L_costraint(self):
        self._insert_many(LSOPT.FIRST_ACCEPT,
                          None,
                          "both of the inserted should not fit with FIRST_ACCEPT",
                          L=8.8)
        self._insert_many(LSOPT.BEST_ACCEPT,
                          None,
                          "both of the inserted should not fit with BEST_ACCEPT",
                          L=8.8)
        self._insert_one(LSOPT.FIRST_ACCEPT,
                         None,
                         "the inserted should not fit with FIRST_ACCEPT",
                         L=9.0,
                         to_insert=5)
        self._insert_one(LSOPT.BEST_ACCEPT,
                         None,
                         "the inserted should not fit with BEST_ACCEPT",
                         L=9.0,
                         to_insert=5)

class TestRedistribute(unittest.TestCase):
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
        self.d = [1.0]*len(self.D)
        self.d[0] = 0.0
        self.longMessage = True

    def test_redistribute_to_one_route(self):
        r1= [0, 1, 4, 6, 0]
        rd1_redistribute = RouteData(r1, route_l(r1, self.D),
                                     route_d(r1, self.d), None)
        r2= [0, 2, 3, 5, 7, 0]
        rd2_recieving = RouteData(r2, route_l(r2, self.D),
                                  route_d(r2, self.d), None)
    
        result = do_redistribute_move(rd1_redistribute, rd2_recieving,
                            self.D, strategy=LSOPT.FIRST_ACCEPT)
        self.assertEqual( len(result), 3, "The redistribute operator should return the redistributed and the new combined routes and the delta")
        self.assertEqual( result[1].route, [0,2,3,5,7,1,4,6,0], "It should be possible to insert all and they sould be appended to the route")
        
        result = do_redistribute_move(rd1_redistribute, rd2_recieving,
                            self.D, strategy=LSOPT.BEST_ACCEPT)
        self.assertEqual( len(result), 3, "The redistribute operator should return the redistributed and the new combined routes and the delta")
        self.assertEqual( result[1].route, [0,1,2,3,4,5,6,7,0], "It should be possible to insert all")
      

    def test_redistribute_to_two_routes(self):
        r1= [0, 1, 4, 6, 0]
        rd1_redistribute = RouteData(r1, route_l(r1, self.D),
                                     route_d(r1, self.d), None)
        r2= [0, 2, 3, 0]
        rd2_recieving = RouteData(r2, route_l(r2, self.D),
                                  route_d(r2, self.d), None)
    
        r3= [0, 7, 5, 0]
        rd3_recieving = RouteData(r3, route_l(r3, self.D),
                                  route_d(r3, self.d), None)
        
        # depending on how the recombinations are made, the results differ
        FI = LSOPT.FIRST_ACCEPT
        BE = LSOPT.BEST_ACCEPT
        tests = [
        (0, "FIRST", FI, ([0, 2, 3, 1, 4, 0],[0, 7, 5, 6, 0])),
        (0, "BEST",  BE, ([0, 1, 2, 3, 4, 0],[0, 7, 6, 5, 0])),
        (1, "FIRST", FI, ([0, 2, 3, 4, 1, 0],[0, 7, 5, 6, 0])),
        (1, "BEST",  BE, ([0, 1, 2, 3, 4, 0],[0, 7, 6, 5, 0])),
        (2, "FIRST", FI, ([0, 2, 3, 1, 4, 0],[0, 7, 5, 6, 0])),
        (2, "BEST",  BE, ([0, 1, 2, 3, 4, 0],[0, 7, 6, 5, 0])),
        (3, "FIRST", FI, ([0, 2, 3, 4, 1, 0],[0, 7, 5, 6, 0])),
        (3, "BEST",  BE, ([0, 1, 2, 3, 4, 0],[0, 7, 6, 5, 0]))]
        for recombination_level, strategy_name, strategy, target_result in tests:
            result = do_redistribute_move(rd1_redistribute, [rd2_recieving, rd3_recieving],
                                self.D, C=4.0, d=self.d, strategy=strategy,
                                recombination_level=recombination_level)
            #print("QUALITY", strategy_name, "/", recombination_level, "delta", result[-1])
            self.assertEqual( len(result), 4, "The redistribute operator should return the redistributed and the new combined routes and the delta")
            self.assertEqual( result[1].route, target_result[0], "It should be possible to insert all and they sould be appended to the route on recombination level %d with strategy %s"%(recombination_level,strategy_name))
            self.assertEqual( result[2].route, target_result[1], "It should be possible to insert all and they sould be appended to the route on recombination level %d with strategy %s"%(recombination_level,strategy_name))
    
    def test_redistribute_does_not_fit(self):
        demands = [0.0, 1, 1, 1, 1, 1, 1, 1]
        r1= [0, 1, 2, 0]
        rd1_redistribute = RouteData(r1, route_l(r1, self.D),
                                     route_d(r1, demands), None)
        r2= [0, 3, 4, 5, 6, 7, 0]
        rd2_recieving = RouteData(r2, route_l(r2, self.D),
                                  route_d(r2, demands), None)
        
        result = do_redistribute_move(rd1_redistribute,
                    [rd2_recieving],
                    self.D, C=6.0, d=self.d, strategy=LSOPT.BEST_ACCEPT,
                    recombination_level=1)
        
        # Fails -> None,None...None is returned (no partial fits)
        self.assertTrue( all(rp==None for rp in result), "The redistribute operator should return the redistributed and the new combined routes and the delta")
        
        
        
    def test_redistribute_find_best_fit(self):
        demands = [0.0, 1, 1, 2, 3, 1, 3, 1]
        r1= [0, 1, 4, 6, 0]
        rd1_redistribute = RouteData(r1, route_l(r1, self.D),
                                     route_d(r1, demands), None)
        r2= [0, 2, 3, 0]
        rd2_recieving = RouteData(r2, route_l(r2, self.D),
                                  route_d(r2, demands), None)
    
        r3= [0, 5, 0]
        rd3_recieving = RouteData(r3, route_l(r3, self.D),
                                  route_d(r3, demands), None)
        
        r4= [0, 7, 0]
        rd4_recieving = RouteData(r4, route_l(r4, self.D),
                                  route_d(r4, demands), None)
        
        result = do_redistribute_move(rd1_redistribute,
                    [rd2_recieving, rd3_recieving, rd4_recieving],
                    self.D, C=4.0, d=demands, strategy=LSOPT.BEST_ACCEPT,
                    recombination_level=1)
        
        self.assertEqual( len(result), 5, "The redistribute operator should return the redistributed and the new combined routes and the delta")
        self.assertEqual( result[0].route, [0,0], "It should be possible to insert all customers" )
        self.assertEqual( _normalise_route_order(result[1].route), [0, 1, 2, 3, 0], "n1 should be redistributed to first route")
        self.assertEqual( _normalise_route_order(result[2].route), [0, 4, 5, 0], "n4 should be redistributed to first route")
        self.assertEqual( _normalise_route_order(result[3].route), [0, 6, 7, 0], "n6 should be redistributed to first route")
            
        
    #def test_redistribute_to_many_routes(self):
    #    pass
    #TODO: write more tests
    
    
    
class TestOnePointMove(unittest.TestCase):
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

    def _make_improving_move(self, C=None, d=None, L=None):
        D = self.D
        r1 = [0,4,5,6,7,0]
        r2 = [0,1,2,3,0]
        r1d = RouteData(r1, route_l(r1,D), 4, None)
        r2d = RouteData(r2, route_l(r2,D), 3, None)
        return do_1point_move(r1d, r2d, D, d=d, C=C, L=L,
                                strategy=LSOPT.BEST_ACCEPT)
 
    def test_improving_move(self):
        new_r1d, new_r2d, delta = self._make_improving_move() 
        self.assertEqual( new_r2d[0], [0,1,2,3,4,0], "node 4 should be inserted to r2")
        self.assertEqual( new_r1d[0], [0,5,6,7,0], "node 4 should be removed from r1")
    
    def test_until_no_improvements_move(self):
        D = self.D
        r1 = [0,1,2,3,5,6,7,0]
        r2 = [0,4,0]
        r1d = (r1, route_l(r1,D), 6, None)
        r2d = (r2, route_l(r2,D), 1, None)
        
        loop_count = 0
        while True:        
            result = do_1point_move(r1d, r2d, D)
            if result[2] is None:
                break
            # unpack result
            r1d, r2d, delta = result
            loop_count+=1
            
            self.assertTrue( loop_count < 8, "Too many operations, potential "+
             "infinite improvement loop")
        
        # unconstrained should run until other route is empty
        self.assertEqual( r2d[0], [0,1,2,3,4,0], "should reach a local optima")
    
    def test_non_improving_move(self):
        D = self.D
        r1 = [0,5,6,7,0]
        r2 = [0,1,2,3,4,0]
        r1d = RouteData(r1, route_l(r1,D), 3, None)
        r2d = RouteData(r2, route_l(r2,D), 4, None)
        result = do_1point_move(r1d, r2d, D)
        self.assertEqual( result, (None, None, None), "r1, r2 configuration should already be at local optima")
        result = do_1point_move(r2d, r1d, D)
        self.assertEqual( result, (None, None, None), "r1, r2 configuration should already be at local optima")
    
    def test_respect_capacity_constraint(self):
        d=[1]*len(self.D)
        result = self._make_improving_move(C=3.5, d=d)        
        self.assertEqual( result, (None,None,None), "no move should be made due to violation of C constraint")
            
    def test_respect_route_cost_constraint(self):
        r2 = [0,1,2,3,0]
        result = self._make_improving_move(L=route_l(r2,self.D)+1.0)    
        self.assertEqual( result, (None,None,None), "no move should be made due to violation of C constraint")
    
    def test_updated_route_cost(self):
        new_r1d, new_r2d, delta = self._make_improving_move()
        self.assertEqual( route_l(new_r1d[0],self.D), new_r1d[1], "original route cost + savings should match recalculated route cost")
        self.assertEqual( route_l(new_r2d[0],self.D), new_r2d[1], "original route cost + savings should match recalculated route cost")
    
    def test_updated_route_capacity(self):
        # route 2 capacity 3.0 , route 1 capacity 4.0 (in _make_improving_move)
        #  node #4 is moved and with it the demand of 1.4 from r1 to r2
        d=[0.0, 0.7,1.0,1.3, 1.4,1.5,0.5,0.6]
        new_r1d, new_r2d, delta = self._make_improving_move(C=5.0, d=d)
        self.assertEqual( route_d(new_r1d[0],d), new_r1d[2], "original route cost + modification should match recalculated route demand")
        self.assertEqual( route_d(new_r2d[0],d), new_r2d[2], "original route cost + modification should match recalculated route demand")

class TestTwoPointMove(unittest.TestCase):
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
    
    def _make_improving_move(self, C=None, d=None, L=None):
        D = self.D
        r1 = [0,5,6,4,0]
        r2 = [0,1,2,3,7,0]
        r1rd = RouteData(r1, route_l(r1,D), route_d(r1, d), None)
        r2rd = RouteData(r2, route_l(r2,D), route_d(r2, d), None)
        return do_2point_move(r1rd, r2rd, D, d=d, C=C, L=L,
                               strategy=LSOPT.BEST_ACCEPT)
        
    def test_improving_move(self):
        new_r1d, new_r2d, f_delta = self._make_improving_move() 
        self.assertEqual( new_r2d[0], [0,1,2,3,4,0], "nodes 4 and 7 should be swapped")
        self.assertEqual( new_r1d[0], [0,5,6,7,0], "nodes 4 and 7 should be swapped")

    def test_with_C_constraint_no_improving(self):
        C = 3.0
        #     #0   #1   #2   #3   #4   #5   #6   #7
        d = [0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5]
        new_r1d, new_r2d, f_delta = self._make_improving_move(d=d, C=C)   
        self.assertEqual( f_delta, None, "There should be no improving moves")
    
    def test_with_C_constraint_improving(self):
        C = 3.0
        #     #0   #1   #2   #3   #4   #5   #6   #7
        d = [0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 0.5]
        new_r1d, new_r2d, f_delta = self._make_improving_move(d=d, C=C) 
        self.assertEqual( new_r2d[0], [0,1,2,6,7,0], "nodes 3 and 7 should be swapped")
        self.assertEqual( new_r1d[0], [0,5,3,4,0], "nodes 3 and 7 should be swapped")

    #def test_with_L_constraint(self):
    #    raise NotImplementedError()
        
    def test_until_no_improvements_move(self):
        D = self.D
        r1 = [0,4,6,1,0]
        r2 = [0,7,2,3,5,0]
        r1d = RouteData(r1, route_l(r1,D), -1, r1[:-1])
        r2d = RouteData(r2, route_l(r2,D), -1, r2[:-1])
        while True:        
            result = do_2point_move(r1d, r2d, D)
            if result[0] is None:
                break
            # unpack result
            r1d, r2d, delta  = result
        
        # unconstrained should run until other route is empty
        self.assertEqual( r2d[0], [0,1,2,3,4,0], "should reach a local optima")
    
    def test_non_improving_move(self):
        D = self.D
        r1 = [0,5,6,7,0]
        r2 = [0,1,2,3,4,0]
        r1d = RouteData(r1, route_l(r1,D), 3, r1[:-1])
        r2d = RouteData(r2, route_l(r2,D), 4, r2[:-1])
        result = do_2point_move(r1d, r2d, D)
        self.assertEqual( result, (None,None,None), "r1, r2 configuration should already be at local optima")
        result = do_2point_move(r2d, r1d, D)
        self.assertEqual( result, (None,None,None), "r1, r2 configuration should already be at local optima")
        
    def test_respect_capacity_constraint(self):
        d=[1]*len(self.D)
        result = self._make_improving_move(C=3.5, d=d)        
        self.assertEqual( result, (None,None,None), "no move should be made due to violation of C constraint")
        
    def test_respect_route_cost_constraint(self):
        r2 = [0,1,2,3,0]
        result = self._make_improving_move(L=route_l(r2,self.D)+1.0)    
        self.assertEqual( result, (None,None,None), "no move should be made due to violation of C constraint")
    
    def test_updated_route_cost(self):
        new_r1d, new_r2d, delta = self._make_improving_move()
        self.assertAlmostEqual( route_l(new_r1d[0],self.D), new_r1d[1], msg="original route cost + savings should match recalculated route cost")
        self.assertAlmostEqual( route_l(new_r2d[0],self.D), new_r2d[1], msg="original route cost + savings should match recalculated route cost")
    
    def test_updated_route_capacity(self):
        # route 2 capacity 3.0 , route 1 capacity 4.0 (in _make_improving_move)
        #  node #4 is moved and with it the demand of 1.4 from r1 to r2
        d=[0.0, 0.7,1.0,1.3, 1.4,1.5,0.5,0.6]
        new_r1d, new_r2d, delta = self._make_improving_move(C=5.0, d=d)
        self.assertAlmostEqual( route_d(new_r1d[0],d), new_r1d[2], msg="original route cost + modification should match recalculated route demand")
        self.assertAlmostEqual( route_d(new_r2d[0],d), new_r2d[2], msg="original route cost + modification should match recalculated route demand")

# The test_local_search_parallel does this
#class TestLSAPI(unittest.TestCase):
#    def test_2opt_from_poor_solution(self):
#        raise NotImplementedError()
        
    
if __name__ == '__main__':
    unittest.main()