# CORE
import sys; print("python",sys.version)
#2.7.12 (default, Dec  4 2017, 14:50:18)
#[GCC 5.4.0 20160609]
import numpy; print("numpy", numpy.__version__)
#1.11.0
import scipy; print("scipy", scipy.__version__)
#0.17.0
import sklearn; print("sklearn", sklearn.__version__)
#'0.17'

# PyPI
import natsort; print ("natsort", natsort.__version__)
#5.3.2
import orderedset; print ("orderedset", orderedset.__version__)
#2.0.1
#import llist; print("llist",llist.__version__)
#Traceback (most recent call last):
#  File "<stdin>", line 1, in <module>
#AttributeError: 'module' object has no attribute '__version__'
print("llist version from https://pypi.org/project/llist/", "0.6")
#0.6

import gurobipy; print("gurobipy", gurobipy.gurobi.version())
#(7, 5, 2)
print("LKH version must be checked manually (the executable has not --version flag)", "Version 2.0.7 - November 2012")
#Version 2.0.7 - November 2012

