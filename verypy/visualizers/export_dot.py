# -*- coding: utf-8 -*-

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division

from math import sin, cos, log

# used to compute rays and their bounding box (bb)
BB_PADDING = 0.2
BB_HALF_SIZE = 5.0

#~defines, for readability
X = 0 
Y = 1

def _segment_length_squared(line1):
    xdiff = line1[0][0] - line1[1][0]
    ydiff = line1[0][1] - line1[1][1]
    return xdiff**2+ydiff**2
    
    
def _segment_intersect_border(line1, lim_x=None, lim_y=None):
    """ Instead of using shapely just adapt the simplest code snipplet kindly
    posted by Stackoverflow users Paul Draper and zidik.
    https://stackoverflow.com/a/20677983
    """
    
    x1, y1 = line1[0]
    x2, y2 = line1[1]
    dx = x2-x1
    dy = y2-y1
    
    if lim_x is not None and dx!=0 and min(x1,x2)<lim_x and max(x1,x2)>lim_x:
        return lim_x, y1+dy/dx*(lim_x-x1)
    elif lim_y is not None and dy!=0 and min(y1,y2)<lim_y and max(y1,y2)>lim_y:
        return x1+dx/dy*(lim_y-y1), lim_y
    else:
        return None

def print_solution_edges(routes, output_handle, color=None):
    for route in routes:
        prev_v = None
        for v in route:
            if prev_v is not None:
                if color:
                    print('    n%04d--n%04d [color="%s"];'%(prev_v, v, color),
                          file=output_handle)  
                else:
                    print('    n%04d--n%04d;'%(prev_v, v),
                          file=output_handle)  
            prev_v = v
        
def output_dot(output_handle, npts, active_point_idxs=None,
               rays=None, active_ray_idxs=None, points_of_interest = None,
               gray_routes=None, red_routes=None, black_routes=None,
               label=None):
    
    
    xs, ys = zip(*npts)
    min_x = min(xs)
    max_x = max(xs)
    min_y = min(ys)
    max_y = max(ys)
    
    # N=1-9->1, N=10-99->2, N=10-99->2, ...
    label_length = int((log(abs(len(npts)),10)+1e-15))+1
    label_format = "%0"+str(label_length)+"d"

    print("""Graph G {""", file=output_handle)  
    if label:
        print('label="%s"'%label, file=output_handle)
    print("node [shape=circle width=.4 fixedsize=true];", file=output_handle)
    print("n0000 [fillcolor=black style=filled fontcolor=white", 
          file=output_handle, end=" ") 
    print('label="'+label_format%0+'" pos="%f,%f!"];'%(npts[0][X], npts[0][Y]),
          file=output_handle)
    for i in range(1,len(npts)):
        if active_point_idxs and i in active_point_idxs:
            print('    n%04d [shape=doublecircle label="'%i+
                  label_format%i+'" pos="%f,%f!"];'%(npts[i][X], npts[i][Y]),
                  file=output_handle)    
        else:
            print('    n%04d [label="'%i+label_format%i+
                  '" pos="%f,%f!"];'%(npts[i][X], npts[i][Y]),
                  file=output_handle)      
    print("", file=output_handle)
        
    if rays is not None:
        print("""ur [shape=none label="" pos="%f,%f!" ];"""%
              (max_x+BB_PADDING, max_y+BB_PADDING), file=output_handle)
        print("""lr [shape=none label="" pos="%f,%f!" ];"""%
              (max_x+BB_PADDING, min_y-BB_PADDING), file=output_handle)
        print("""ul [shape=none label="" pos="%f,%f!" ];"""%
              (min_x-BB_PADDING, max_y+BB_PADDING), file=output_handle)
        print("""ll [shape=none label="" pos="%f,%f!" ];"""%
              (min_x-BB_PADDING, min_y-BB_PADDING), file=output_handle)
                
        b1 = (max_x+BB_PADDING,None)
        b2 = (min_x-BB_PADDING,None)
        b3 = (None, max_y+BB_PADDING)
        b4 = (None, min_y-BB_PADDING)
        
        for ri, ray_alpha in enumerate(rays):
            # find a ray that gets clipped by the bounding box
            from_pt = (npts[0][X], npts[0][Y])
            to_pt = (npts[0][X]+3*(BB_HALF_SIZE)*cos(ray_alpha),
                     npts[0][Y]+3*(BB_HALF_SIZE)*sin(ray_alpha))
            ray_len_sqrd = _segment_length_squared( (from_pt, to_pt) )
            long_ray = [from_pt,to_pt]
            
            intersection_n = 0
            for x_lim, y_lim in [b1,b2,b3,b4]:
                isect = _segment_intersect_border(long_ray, x_lim, y_lim)
                
                if isect:
                    intersection_n+=1
                    #print("""is%d [shape=point label="" pos="%f,%f!" ];"""%(intersection_n, isect[0], isect[1]), file=output_handle)
                    
                    isect_ray_len_sqrd = _segment_length_squared( (from_pt, isect) )
                    if isect_ray_len_sqrd<ray_len_sqrd:
                        to_pt = isect
                        ray_len_sqrd = isect_ray_len_sqrd
            
            #print("ray_endpoint", ri, ray_alpha, "x,y", to_pt[X], to_pt[Y])
            
            # hidden endpoint of the ray
            print("""    r%04d [shape=none label="" pos="%f,%f!" ];"""%
                  (ri, to_pt[X], to_pt[Y]), file=output_handle)
            # draw ray
            if ri in active_ray_idxs:
                print("    n%04d--r%04d [style=dashed];"%(0,ri), file=output_handle)
            else:
                print("    n%04d--r%04d [style=dotted];"%(0,ri), file=output_handle)
           
    if points_of_interest:
        for pi, poi in enumerate(points_of_interest):
            print("""    poi%04d [fontsize=20 shape=none label=<<b>×</b>> pos="%f,%f!" ];"""%
                  (pi, poi[X], poi[Y]), file=output_handle)
      
    if gray_routes:
        print_solution_edges(gray_routes, output_handle, "gray")
    if red_routes:
        print_solution_edges(red_routes, output_handle, "red")  
    if black_routes:
        print_solution_edges(black_routes, output_handle)
    
    print("}", file=output_handle)