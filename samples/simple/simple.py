#!/usr/bin/env python

#
# Make sure that the egfrd system is added to your PYTHONPATH
# This means, in bash for example:
# $ export PYTHONPATH=$HOME/egfrd
#

# Modules
import sys

from egfrd import *
import model
import gfrdbase

import logger


# Constants
size = 1e-6

# Model
m = model.ParticleModel(size)
# Species
P = model.Species('P', 1e-12, 3e-9)                                   
m.add_species_type(P) 

# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

# Throw in particles
throw_in_particles(w, P, 60)

# Logger
l = logger.Logger('simple')
interrupter = logger.FixedIntervalInterrupter(s, 3.33e-4, l.log)

# Simulation
l.start(s)
while s.t < .1:
    interrupter.step()
