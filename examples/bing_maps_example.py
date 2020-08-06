#!/usr/bin/env python3

"""This script calculates distance matrix for given addresses using Bing maps API.

More specifically:
1) it reads adresses from a file (tab separated fields, one per line format):
   
   addressLine \t postalCode \t locality  \t countryRegion [\t goodsDemand]
   
   where addressLine is, e.g., the street and number, 
   locality is, e.g., the city/town name, 
   countryRegion is ISO country code, e.g. "US", and
   optional goodsDemand is the demand for CVRP

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
    
import json
import sys
import pickle
import argparse

__author__ = "Jussi Rasku"
__copyright__ = "Copyright 2020"
__license__ = "MIT"
__version__ = "1.0.0"
__status__ = "Example"

# TODO: include as a CLI parameter
TRAVEL_MODE = 'driving' #'walking'

from pprint import pprint
import numpy as np

def completeDataWithGeocoding(addressData, BingMapsKey):
    reqData = dict(addressData)
    reqData['locality'] = quote(addressData['locality'].strip(), safe='')
    reqData['addressLine'] = quote(addressData['addressLine'].strip(), safe='')
    reqData['BingMapsKey'] = BingMapsKey

    geogodeUrl =  ("http://dev.virtualearth.net/REST/v1/Locations"+
                   "?countryRegion="+str(reqData['countryRegion'])+
                   "&locality="+str(reqData['locality'])+
                   "&postalCode="+str(reqData['postalCode'])+
                   "&addressLine="+str(reqData['addressLine'])+
                   "&maxResults=1"+
                   "&key="+str(reqData['BingMapsKey']))
    #print(geogodeUrl)
    request = Request(geogodeUrl)
    response = urlopen(request)
    result = json.loads(response.read().decode(encoding='utf-8'))
    lat, lon = result['resourceSets'][0]['resources'][0]['point']['coordinates']

    retData = dict(addressData)
    retData['latitude'] = lat
    retData['longitude'] = lon
    return retData

def getDistanceMatrix(addressDataList, travelMode, BingMapsAPIKey):
    coordinateList = []
    for ad in addressDataList:
        coordinateList.append( str(ad['latitude'])+","+str(ad['longitude']) )
    coordinates = ";".join(coordinateList)

    distanceUrl = ("https://dev.virtualearth.net/REST/v1/Routes/DistanceMatrix"+
        "?origins="+coordinates+"&destinations="+coordinates+
        "&travelMode="+str(travelMode)+
        "&distanceUnit=km"+
        "&key="+str(BingMapsAPIKey))
    #print(distanceUrl)
    request = Request(distanceUrl)
    response = urlopen(request)
    result = json.loads(response.read().decode(encoding='utf-8'))
    #pprint(result)

    N = len(addressDataList)
    D = np.zeros( (N,N) )
    for cell in result['resourceSets'][0]['resources'][0]['results']:
        D[cell['originIndex'], cell['destinationIndex']] =  cell['travelDistance'] # could also be cell['travelDuration'] 
    return D

# CLI
parser = argparse.ArgumentParser(description=__doc__)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', type=argparse.FileType('r'), metavar="address_file", help='address file')
group.add_argument('-D', type=argparse.FileType('rb'), metavar="distance_matrix_file", help='pickled Numpy distance matrix file')
parser.add_argument('-C', type=float, metavar="capacity", help='solve as CVRP with vehicles with this capacity')
parser.add_argument('-o', type=argparse.FileType('wb'), metavar="output", help='output file for pickled Numpy distance matrix')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
args = parser.parse_args()

locations = []
if args.f:
    mapApiKey = ""
    try:
        mapApiKey = open('MyBingMapApi.key', mode='r').read().strip()
    except:
        sys.stderr.write("ERROR: No file called \"MyBingMapApi.key\"\n")
        sys.stderr.write("It is a text file that should contain your personal Microsoft Bing Maps API key.\n")
        sys.stderr.write("You can get one from here: https://www.bingmapsportal.com/\n")
        sys.stderr.write("Exiting.\n")
        sys.exit(1)

    for l in args.f.readlines():
        l = l.strip()
        if len(l)==0:
            continue
        parts = l.split('\t')
        if len(parts)<4 or len(parts)>5:
            sys.stderr.write("ERROR: Invalid address data line \""+l+"\"\n")
            sys.stderr.write("Data lines should be in format: addressLine \\t postalCode \\t locality \\t countryRegion [\\t goodsDemand]\n")
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

        if args.verbose:
            sys.stderr.write("INFO: Geocoding \""+location["addressLine"]+"\"\n")
        locations.append( completeDataWithGeocoding(location, mapApiKey) )

    if args.verbose:
        sys.stderr.write("INFO: Filling %d x %d distance matrix."%(len(locations), len(locations)))

    D = getDistanceMatrix(locations, TRAVEL_MODE, mapApiKey)
    if args.o:
        pickle.dump(D, args.o)
    pprint(D)
elif args.D:
    D = pickle.load(args.D)
    
else:
    assert(False) # argparse should take care that we never arrive here
    pass

if args.C:
    from classic_heuristics.parallel_savings import parallel_savings_init
    from util import sol2routes

    if locations:
        d = [loc["demand"] for loc in locations]
    else:
        d = [0.]+[1.0]*(len(D)-1)
        if args.verbose:
            sys.stderr.write("Warning: Setting all goods demands to be 1.0 (no data).\n")
    
    if args.verbose:
        sys.stderr.write("INFO: Solving a %d customer CVRP with Savings heuristic.\n"%(len(locations)-1))

    solution = parallel_savings_init(D=D, d=d, C=args.C)
    print("\nCVRP solution\n")
    for route_idx, route in enumerate(sol2routes(solution)):
        print("Route #%d : %s"%(route_idx+1, route))
