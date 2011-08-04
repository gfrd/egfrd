#!/usr/bin/env python

#
# Make sure that the egfrd system is added to your PYTHONPATH
# This means, in bash for example:
# $ export PYTHONPATH=$HOME/egfrd
#
#
# Run with bash command:
# $ python tf.py

import sys
import os

from egfrd import *

import model
import gfrdbase
import _gfrd
import logger
import dumper
#from bd import *

# Constants
k_fR = 6e9 * 1000 / N_A
k_bR = 0.1  # 1 - 0.01
k_f_rp = 38
k_b_rp = 0.5
k_OC = 1 # 0.3 - 3
t_clear = 1  # should not be poisson
t_elon = 50 # 50-100
k_dm = 0.019
k_ribo = 5 * k_dm
k_dp = 2.4e-4
t_trans = 30

# Model
m = model.ParticleModel(1e-6)
# Species
O = model.Species('O', 0, 1e-8)
R = model.Species('R', 1e-12, 1e-8)
P = model.Species('R', 1e-12, 1e-8)
OR = model.Species('OR', 0, 1e-8)
ORp = model.Species('ORp', 0, 1e-8)
ORpa = model.Species('ORpa', 0, 1e-8)
T = model.Species('T', 1e-12, 1e-8)
M = model.Species('M', 1e-12, 1e-8)
Mribo = model.Species('Mribo', 1e-12, 1e-8)

m.add_species_type(O) 
m.add_species_type(R) 
m.add_species_type(P) 
m.add_species_type(OR) 
m.add_species_type(ORp) 
m.add_species_type(ORpa) 
m.add_species_type(T) 
m.add_species_type(M) 
m.add_species_type(Mribo) 

# Reaction rules
r1 = model.create_binding_reaction_rule(O, R, OR, k_fR)
m.network_rules.add_reaction_rule(r1)
r2 = model.create_unbinding_reaction_rule(OR, O, R, k_bR)
m.network_rules.add_reaction_rule(r2)
r3 = model.create_unimolecular_reaction_rule(O, ORp, k_f_rp)
m.network_rules.add_reaction_rule(r3)
r4 = model.create_unimolecular_reaction_rule(ORp, O, k_b_rp)
m.network_rules.add_reaction_rule(r4)
r5 = model.create_unimolecular_reaction_rule(ORp, ORpa, k_OC)
m.network_rules.add_reaction_rule(r5)
r6 = model.create_unbinding_reaction_rule(ORpa, T, O, 1/t_clear)
m.network_rules.add_reaction_rule(r6)
r7 = model.create_decay_reaction_rule(M, k_dm)
m.network_rules.add_reaction_rule(r7)
r8 = model.create_unbinding_reaction_rule(M, M, Mribo, k_ribo)
m.network_rules.add_reaction_rule(r8)
r9 = model.create_unimolecular_reaction_rule(Mribo, P, 1/t_trans)
m.network_rules.add_reaction_rule(r9)
r10 = model.create_decay_reaction_rule(P, k_dp)
m.network_rules.add_reaction_rule(r10)

# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

#####

#EMPTY = model.Species('EMPTY', 2e-12, 5e-8)

#  1 2 O + R <-> OR
#  3 4 O     <-> ORp
#  5   ORp    -> ORpa
#  6   ORpa   -> T + O
#  7   M      -> EMPTY
#  8   M      -> M + Mribo
#  9   Mribo  -> P
# 10   P      -> EMPTY


place_particle(w, O, [0,0,0])

#s.throw_in_particles(R, 50, box1)

l = logger.Logger('pushpull')
interrupter = logger.FixedIntervalInterrupter(s, 1e-3, l.log)

l.start(s)
while s.t < 1000:
    interrupter.step()
    dumper.dump_particles(s)
