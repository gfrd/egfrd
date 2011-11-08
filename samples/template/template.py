#!/usr/bin/python

# Template file to copy and paste basic code from.
#
#
# Make sure that the egfrd system is added to your PYTHONPATH
# This means, in bash for example:
# $ export PYTHONPATH=$HOME/egfrd
#

# Modules
# ===============================
import sys
import os

from egfrd import *

import model
import gfrdbase
import _gfrd
#import logger
#import dumper
#from bd import *


# Constants
# ===============================
sigma = 1e-4        # Diameter particle; more realistic would be e-9
D = 1e-12           # Diffusion constant
N = 1000            # Number of steps simulation will last
world_size = 1e-3   # Lengths of simulation box


# Basic set up simulator
# ===============================
# Model
m = model.ParticleModel(world_size)
# Species
X = model.Species('X', D, sigma/2)                                   
m.add_species_type(X) 
# Reaction rules
### Binding: A + B -> C
r1 = model.create_binding_reaction_rule(A, B, C, k)
m.network_rules.add_reaction_rule(r1)
### Unbinding: C -> A + B
r2 = model.create_unbinding_reaction_rule(C, A, B, k)
m.network_rules.add_reaction_rule(r2)
# Unimolecular: A - > B
r3 = model.create_unimolecular_reaction_rule(A, B, k)
m.network_rules.add_reaction_rule(r3)
### Decay: A - > 0
r4 = model.create_decay_reaction_rule(A, k)
m.network_rules.add_reaction_rule(r4)
### Annihilation: A + B - > 0
r5 = model.create_annihilation_reaction_rule(A, B, k)
m.network_rules.add_reaction_rule(r5)

# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

# Throw in particles
throw_in_particles(w, X, numberToThrowIn)
place_particle(w, X, [0,0,0])

# Running 1
# ===============================
while (s.get_next_time() < end_time):
    s.step()

s.stop(end_time)
    

# Running 2
# ===============================
for t in range(N):
    s.step()

s.stop(s.t)


# Additional code
# ===============================

# VTK Logger
# ===============================
vtk_output_directory = 'VTK_out'
if (os.path.exists(vtk_output_directory)):
    print '** WARNING: VTK output directory already exists, possible '+ \
            'present old VTK files might lead to errors.'
    
from visualization import vtklogger
l = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 

# log a step
l.log()
# stop logger, writing away essential data
l.stop()

