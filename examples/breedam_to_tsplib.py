#!/usr/bin/env python

"""Converts van Breedam VRP problem files to the TSPLIB95 format.

Execute as script, provide the globbed inputs and an output folder, 
and see it go.

Can also be studied as an example of generating .vrp files.
"""

# Try to support Python2
from __future__ import print_function
from __future__ import division

from os import path
import glob
import string

""" Converts a single problem file from format to another. """
def breedam_to_tsplib95(in_filepath, out_filepath):
    data = {'N':0, # number of customers
            'Q':0, # capacity constraint
            'xys':[], # x,y coordinates
            'qs':[]}  # demands (qs[0]=0 b/c it is the depot)

    # The capacity constraint varies between the van Breedam sets,
    #  and is set based on the problem name.
    bn = path.basename(in_filepath)
    if 'G1' in bn or "1PP" in bn or "2PP" in bn:
        data['Q']=100
    elif 'G2' in bn or "3PP" in bn or "4PP" in bn:
        data['Q']=50
    elif 'G3' in bn or "5PP" in bn or "6PP" in bn:
        data['Q']=200
        
    # read
    bf = open(in_filepath, 'r')  
    N = 0
    for l in bf.readlines():
        # strip non-printable non-ascii characters
        l = ''.join(filter(lambda c: c in string.printable, l))
        parts = [int(e) for e in l.split()]
        
        if len(parts)!=10:
            continue
        
        # van Breedam problem format. Each line has some or all of following info:
        # <id>, <x_coord>, <y_coord>, <tw1_a>, <tw1_b>, <tw2_a>, <tw2_b>, <demand>, <service>, <delivery_or_pickup>
        idx, x, y, _, _, _, _, q, _, _  = parts
        # As one can see, we stay in the realm of CVRP
        
        data['xys'].append( (idx+1,x,y) )
        data['qs'].append( (idx+1,q) )
        N+=1
    bf.close()
    
    print("read", in_filepath)
    
    # write
    tf = open(out_filepath, 'w')  
    tf.write("NAME : %s\n"%bn)
    tf.write("COMMENT : Automatically converted from a van Breedam DAT file.\n")
    tf.write("TYPE : CVRP\n")
    tf.write("DIMENSION : %d\n" % N)
    tf.write("EDGE_WEIGHT_TYPE : EUC_2D\n")
    tf.write("CAPACITY : %d\n" % data['Q'])
    tf.write("NODE_COORD_SECTION\n")
    for coodinate_data in data['xys']:
        tf.write("%d %d %d\n" % coodinate_data)
    tf.write("DEMAND_SECTION\n")
    for demand_data in data['qs']:
        tf.write("%d %d\n" % demand_data)
    tf.write("DEPOT_SECTION\n")
    tf.write(" 1\n -1\n")
    tf.write("EOF\n")
    tf.close()
    
    print("wrote", out_filepath, "\n")


def main(in_files_glob, output_folder):
    for in_fn in glob.glob(in_files_glob):
        in_bn = path.basename(in_fn)
        out_fn = path.join(output_folder, path.splitext(in_bn)[0]+".vrp")
        breedam_to_tsplib95(in_fn, out_fn)


if __name__=="__main__":
    def folder_type(astring):
        if not path.isdir(astring):
            raise ValueError
        return astring
        
    import argparse
    parser = argparse.ArgumentParser(description="Convert van Breedam .DAT VRP problem files to TSPLIB95 format")
    parser.add_argument('inglob', help='the input file or glob (pathname pattern expansion)')
    parser.add_argument('outdir', help='the output directory', type=folder_type)
    args = parser.parse_args()

    main(args.inglob, args.outdir)
