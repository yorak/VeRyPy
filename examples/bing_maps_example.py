#!/usr/bin/env python

"""This script calculates distance matrix for given addresses using Bing maps API.

More specifically:
1) it reads adresses from a file (tab separated fields, one per line format):
   
   addressLine \t postalCode \t locality  \t countryRegion [\t goodsDemand [\t comment]] 
   
   where addressLine is, e.g., the street and number, 
   locality is, e.g., the city/town name, 
   countryRegion is ISO country code, e.g. "US", and
   optional goodsDemand is the demand for CVRP
   if goodsDemand is given, also a commend may be added to it.

2) geocodes them (gets coordinates) using Bing maps API,

3) and calculates a distance matrix using Bing maps API.

4) Optionally dups the Numpy distance matrix to file.

5) Interpretes the read and computed data as a capacitated vehicle routing problem (CVRP)
   and solves it using Clarke and Wright savings heuristic.

To use it you need Bing maps API key. Get it here: https://www.bingmapsportal.com
and place it into a text file named MyBingMapApi.key, where the script reads it."""

# Written in Python 2.7, but try to maintain Python 3+ compatibility
from __future__ import print_function
from __future__ import division
try:
    from urllib.parse import urlparse, urlencode, quote
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode, quote
    from urllib2 import urlopen, Request, HTTPError

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2020"
__license__ = "MIT"
__version__ = "1.0.0"
__status__ = "Example"

import json
import sys
import pickle
import argparse
import os

from pprint import pprint
import numpy as np


# Due to the limitation of the free api, one has to batch the distance queries.
BING_DM_API_CUSTOMER_LIMIT = 50 # 50 x 50 = 2500
# Use a key loaded from the file "MyBingMapApi.key" containing only *your* key.
BING_MAP_API_KEY = ""
try:
    BING_MAP_API_KEY = open('MyBingMapApi.key', mode='r').read().strip()
except:
    sys.stderr.write("ERROR: No file called \"MyBingMapApi.key\" found in \"%s\"\n"%os.getcwd() )
    sys.stderr.write("It is a text file that should contain your personal Microsoft Bing Maps API key.\n")
    sys.stderr.write("You can get one from here: https://www.bingmapsportal.com/\n")
    sys.stderr.write("Exiting.\n")
    sys.exit(1)

def geocode_data(address_data, bing_maps_key):
    req_data = dict(address_data)
    req_data['locality'] = quote(address_data['locality'].strip(), safe='')
    req_data['addressLine'] = quote(address_data['addressLine'].strip(), safe='')
    req_data['bing_maps_key'] = bing_maps_key

    geogode_url =  ("http://dev.virtualearth.net/REST/v1/Locations"+
                   "?countryRegion="+str(req_data['countryRegion'])+
                   "&locality="+str(req_data['locality'])+
                   "&postalCode="+str(req_data['postalCode'])+
                   "&addressLine="+str(req_data['addressLine'])+
                   "&maxResults=1"+
                   "&key="+str(req_data['bing_maps_key']))
    #print(geogodeUrl)
    request = Request(geogode_url)
    response = urlopen(request)
    result = json.loads(response.read().decode(encoding='utf-8'))
    lat, lon = result['resourceSets'][0]['resources'][0]['point']['coordinates']

    ret_data = dict(address_data)
    ret_data['latitude'] = lat
    ret_data['longitude'] = lon
    return ret_data

def request_distance_matrix(address_list, travel_mode, bing_maps_key):
    coordinateList = []
    for a in address_list:
        coordinateList.append( str(a['latitude'])+","+str(a['longitude']) )
    coordinates = ";".join(coordinateList)

    distanceUrl = ("https://dev.virtualearth.net/REST/v1/Routes/DistanceMatrix"+
        "?origins="+coordinates+"&destinations="+coordinates+
        "&travelMode="+str(travel_mode)+
        "&distanceUnit=km"+
        "&key="+str(bing_maps_key))
    #print(distanceUrl)
    request = Request(distanceUrl)
    response = urlopen(request)
    result = json.loads(response.read().decode(encoding='utf-8'))
    #pprint(result)

    N = len(address_list)
    D = np.zeros( (N,N) )
    for cell in result['resourceSets'][0]['resources'][0]['results']:
        D[cell['originIndex'], cell['destinationIndex']] =  cell['travelDistance'] # could also be cell['travelDuration'] 
    return D

def main(args):
    ## Read addresses from the file, geocode them, and compute the distance matrix using Bing maps API ##
    locations = []
    if args.address_file:
        for l in args.address_file.readlines():
            l = l.strip()
            if len(l)==0:
                continue
            parts = l.split('\t')
            if len(parts)<4 or len(parts)>6:
                sys.stderr.write("ERROR: Invalid address data line \""+l+"\"\n")
                sys.stderr.write("Data lines should be in format: addressLine \\t postalCode \\t locality \\t countryRegion [\\t goodsDemand [\\t comment]]\n")
                sys.stderr.write("Exiting.\n")
                sys.exit(1)
            location = {
                "addressLine":parts[0],
                "postalCode":parts[1],
                "locality":parts[2],
                "countryRegion":parts[3],
                "demand":1 if len(locations)>0 else 0
            }
            if len(parts)==5:
                location["demand"]=float(parts[4])
            if len(parts)==6:
                location["comment"]=parts[5]

            if args.verbosity:
                sys.stderr.write("INFO: Geocoding \""+location["addressLine"]+"\"\n")
            locations.append( geocode_data(location, BING_MAP_API_KEY) )

        if args.verbosity:
            sys.stderr.write("INFO: Filling %d x %d distance matrix.\n"%(len(locations), len(locations)))

        if (len(locations)>BING_DM_API_CUSTOMER_LIMIT):
            sys.stderr.write("ERROR: Bing Distance Matrix API is limited to a query with %d\n"%(BING_DM_API_CUSTOMER_LIMIT*BING_DM_API_CUSTOMER_LIMIT))
            sys.stderr.write("distances. Hence, at most %d coordinates can be queried in one go.\n")
            sys.stderr.write("Exiting.\n")
            sys.exit(1)
            
        D = request_distance_matrix(locations, args.travel_mode, BING_MAP_API_KEY)
        if args.output_file:
            pickle.dump(D, args.output_file)
        print("D = ", end="");pprint(D)
    elif args.distance_matrix_file:
        D = pickle.load(args.distance_matrix_file)
    else:
        assert(False) # argparse should take care that we never arrive here
        pass

    ## Solve the related CVRP using VeRyPy ##
    if args.C:
        #Got an error? Remember to set your PYTHON_PATH to point to VeRyPy folder.
        from verypy.util import sol2routes
        from verypy.classic_heuristics.parallel_savings import parallel_savings_init

        if locations:
            d = [loc["demand"] for loc in locations]
        else:
            d = [0.]+[1.0]*(len(D)-1)
            if args.verbosity:
                sys.stderr.write("WARNING: Setting all goods demands to be 1.0 (no data).\n")
        
        if args.verbosity:
            sys.stderr.write("INFO: Solving a %d customer CVRP with Savings heuristic.\n"%(len(locations)-1))
        
        # Solve and print the resulting routes
        solution = parallel_savings_init(D=D, d=d, C=args.C)
        print("\nCorresponding CVRP solution is")
        for route_idx, route in enumerate(sol2routes(solution)):
            print("Route #%d : %s"%(route_idx+1, route))
    else:
        sys.stderr.write("WARNING: No capacity given, no CVRP to solve (give one with -C <float> command line argument).\n")
        

if __name__ == "__main__":
    ## CLI specification ##
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', type=argparse.FileType('r'), dest='address_file', help="the tsv (tab separated values) file to load addresses from")
    group.add_argument('-t', type=str, dest='travel_mode', default='Driving', help="travel mode, i.e., 'Driving' (default), 'Walking' or 'Transit' as per Bing Maps API" )
    group.add_argument('-D', type=argparse.FileType('rb'), dest="distance_matrix_file", help="pickled Numpy distance matrix file")
    parser.add_argument('-C', type=float, help="capacity of the vehicles (i.e. C in CVRP)")
    parser.add_argument('-o', type=argparse.FileType('wb'), dest='output_file', help="output file for pickled Numpy distance matrix")
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbosity', help="verbose output")
    args = parser.parse_args()
    main(args)