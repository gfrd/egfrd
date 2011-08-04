# Simple script to assess memory usage of eGFRD

# Modules
# ===============================
import sys
#Also, set relative python path
sys.path.append('../../')
import os

#from memorymonitor import MemoryMonitor
import memoryusage
from egfrd import *

import model
import gfrdbase
import _gfrd
#import logger
#import dumper
#from bd import *


# Constants
# ===============================
username = 'wehrens'
outFilename = 'report.txt'

# Number of runs and number of timesteps
runs = 1000
N = 1000

# Number of particles
Np = 100 

# General constants
sigma = 1e-9        # Diameter particle; more realistic would be e-9
D = 1e-12           # Diffusion constant
world_size = 1e-3   # Lengths of simulation box

# Reaction constants
kf = 1e-18
kb = 1


# Memory Monitor
# ===============================
# memory_mon = MemoryMonitor(username)


def singlerun(N):
    """ Creates simple simulation set up and runs N timesteps """

    # Basic set up simulator
    # ===============================
    # Model
    m = model.ParticleModel(world_size)
    # Species
    A = model.Species('A', D, sigma/2)                                   
    m.add_species_type(A) 
    B = model.Species('B', D, sigma/2)                                   
    m.add_species_type(B) 
    C = model.Species('C', D, sigma/2)                                   
    m.add_species_type(C) 
    X = model.Species('X', D, sigma/2)                                   
    m.add_species_type(X) 
    # Reaction rules
    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)
    r2 = model.create_unbinding_reaction_rule(A, B, C, kb)
    m.network_rules.add_reaction_rule(r2)

    # World
    w = gfrdbase.create_world(m, 3)
    # Simulator
    s = EGFRDSimulator(w, myrandom.rng)

    # Throw in particles
    throw_in_particles(w, A, Np)
    throw_in_particles(w, B, Np)
    throw_in_particles(w, C, Np)
    throw_in_particles(w, X, Np)

    # Running
    # ===============================
    for t in range(N):
        s.step()

    s.stop(s.t)

print "Starting."

# Go
outFile = open(outFilename, 'w')
outFile.close()
outFile = open(outFilename, 'a')

for M in range(runs):
    print "Run " + str(M),
    singlerun(N)
    outFile.write(str(memoryusage.memory()) + "\n")
    outFile.flush    
    print "(" + str(memoryusage.memory()) + ")."

outFile.close()

print "Finished."



















