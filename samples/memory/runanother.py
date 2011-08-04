# Simple script to assess memory usage of eGFRD

# Modules
# ===============================
import sys
#Also, set relative python path
sys.path.append('../../')
import os
sys.path.append('/storage1/wehrens/lib/python/')
sys.path.append('/storage1/wehrens/lib/python2.6/site-packages')
import guppy
import gc
from dozer import Dozer

#from memorymonitor import MemoryMonitor
import memoryusage
from egfrd import *

import model
import gfrdbase
import _gfrd
#import logger
#import dumper
#from bd import *

import cluster


# Constants
# ===============================

username = 'wehrens'
outFilename = 'anotherreport.txt'

# runs
runs = 100

# Clustersize
N = 50

# Initializing
# ===============================

# GC
gc.set_debug(gc.DEBUG_LEAK)
sys.stderr = open('bugreport.txt', 'w')
print str(gc.get_debug())
# Dozer
wsgi_app = Dozer(wsgi_app)

# Functions
# ================================

def cleaningcrew(themess):
    # print str(themess)

#    for mess in themess:
#        del(mess)
#    for key, item in cluster.__dict__.items():
#        tokill = getattr(cluster, key)
#        print str(tokill)
#        del(tokill)

    # print "AFTER CLEANING:" + str(cluster.__dict__)

#    methodToCall = getattr(foo, 'bar')
#    result = methodToCall()
    pass

# print "GC enabled?-: " + str(gc.isenabled())

# print "Initializing cyclops."

#import Cyclops
#z = Cyclops.CycleFinder()

# print "Starting."

# Go
outFile = open(outFilename, 'w')
outFile.close()
outFile = open(outFilename, 'a')

for M in range(runs):
    print "run " + str(M),
    #print "Run " + str(M),
    cluster.single_run(N, False)
    outFile.write(str(memoryusage.memory()) + "\n")
    outFile.flush    
    #print "(" + str(memoryusage.memory()) + ").\n"
    #print "CREW:"
#    dirt = []
#    cleaningcrew(dir()) 
    #print str(gc.garbage)
    print "Annoying items: " + str(gc.collect())
    #print "END CREW"
#    del cluster.w
#    del cluster.s # If this one is thrown in, the runs suddenly stop?!
#    del cluster.m
    print "(" + str(memoryusage.memory()) + ")."    

outFile.close()

#print "Memory usage after exiting loop & closing file: " + str(memoryusage.memory()) + "."

#print "Looking for cycles"

#print "z.find_cycles() returned: " + str(z.find_cycles())

#print "Finished."

# sys.stderr.close()






















