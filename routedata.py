# -*- coding: utf-8 -*-
from util import sol2routes, objf

class RouteData:
    def __init__(self, route=None, cost=0.0, demand=0.0, node_set=None,
                 generate_node_set=False):
        self.route = [0,0] if (route is None) else route
        self.cost = cost
        self.demand = demand
        self.node_set = node_set
        self.aux_data_updated = False
        self.fwd_l = None
        self.rwd_l = None
        self.fwd_d = None
        self.rwd_d = None
        
        if generate_node_set:
            assert route is not None, "No route to generate it from"
            self.node_set = RouteData._route_to_nodeset(route)

    def __str__(self):
        return "%s (d=%.2f, f=%.2f)"%(self.route, self.demand, self.cost)
        
    def __iter__(self):
        """ allows unpacking the route_data e.g. """
        for e in (self.route,self.cost,self.demand,self.node_set):
            yield e
    
    def __getitem__(self, i):
        if (i==0):
            return self.route
        elif (i==1):
            return self.cost
        elif (i==2):
            return self.demand
        elif (i==3):
            return self.node_set
        else:
            raise IndexError("0=route, 1=cost, 2=demand, 3=node_set (can be None) and that's it!")
    
    def is_empty(self):
        if not self.route or self.route == [0,0]:
            return True
        else:
            return False
        
    def update_node_set(self):
        self.node_set = RouteData._route_to_nodeset(self.route)
        
    def update_auxiliary_data(self, D, d, direction=None):
        """ This data makes checking C and L constraints linear time-
            direction can be -1 or 1. None means both """
        if not direction:
            self.update_auxiliary_data(D, d, direction=1)
            self.update_auxiliary_data(D, d, direction=-1)
            return
            
        route_d = [0.0]*len(self.route)
        route_l = [0.0]*len(self.route)
        
        cum_d = 0.0
        cum_l = 0.0
        prev_n = 0
        start_i = 1 if direction==1 else -2
        i = start_i
        
        #TODO: there might be a faster (numpy?) way of doing this
        # e.g. using cumsum
        for n in self.route[start_i::direction]:
            if d:
                cum_d += d[n]
                route_d[i] = cum_d
                
            cum_l += D[prev_n,n]
            route_l[i] = cum_l
            prev_n = n
            i+=direction
            
        if direction>0:
            self.fwd_l = route_l
            self.fwd_d = route_d
        elif direction<0:
            self.rwd_l = route_l
            self.rwd_d = route_d
            
    def normalize(self):
        """ The smaller route start/end node comes first """
        if self.route[1]>self.route[-2]:
            self.route.reverse()
    
    @staticmethod
    def _route_to_nodeset(route):
        if not route or route==[0,0]:
            return set()
        start_from = 0
        while route[start_from]==0 and start_from<len(route)-1:
            start_from += 1
        if route[-1]==0:
            end_to = -2
            while route[end_to]==0 and abs(end_to)<len(route):
                end_to -= 1 
            r_nodes = set(route[start_from:end_to+1])
        else:
            r_nodes = set(route[start_from:])
            
        return r_nodes
        
    @staticmethod
    def from_routes(routes, D, d):
        # consruct route data array
        route_datas = []
        for r in routes:
            r_cost = objf(r, D)
            r_demand = sum( d[node] for node in r ) if d else 0
            r_nodes = RouteData._route_to_nodeset(r)
            route_datas.append( RouteData(r, r_cost, r_demand, r_nodes) )    
        return route_datas
        
    @staticmethod
    def from_solution(solution, D, d):
        routes = sol2routes(solution)
        return RouteData.from_routes(routes, D, d)
        
    @staticmethod
    def to_solution(route_datas):
        return [0]+[n for rd in route_datas for n in rd.route[1:]]
        
        
        
        