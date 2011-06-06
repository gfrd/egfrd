#!/usr/bin/env python

# Call with:
# PYTHONPATH=../../ python pushpull.py [Keq] [koff_ratio] [N_S_total] [N_K] [N_P] [V] [mode] [T]

import sys
import os

from egfrd import *
import gfrdbase
import _gfrd
import logger
import dumper
import model

from fractionS import *


# Args:
# Keq
# koff_ratio
# N_K
# N_P
# V (liter)
# mode:  'normal' 'immobile' 'localized' 'single' 'clustered'
# T

Keq_str = sys.argv[1]

koff_ratio_str = sys.argv[2]
N_S_total = int(sys.argv[3])
N_K = int(sys.argv[4])
N_P = int(sys.argv[5])
V_str = sys.argv[6]
mode = sys.argv[7]
T_str = sys.argv[8]


Keq = float(Keq_str)
koff_ratio = float(koff_ratio_str)
V = float(V_str)
T = float(T_str)

radius = 2.5e-9
sigma = radius * 2
D1 = 1.0e-12



if mode == 'normal':
    D2 = D1
elif mode == 'immobile' or mode == 'localized' or mode == 'single':
    D2 = 0
else:
    raise 'invalid mode'


L = (V * 1e-3) ** (1.0 / 3.0)


N = N_S_total * 1.1
matrix_size = min(max(3, int((3 * N) ** (1.0/3.0))), 60)
print 'matrix_size=', matrix_size


# Model
m = model.ParticleModel(L)
# Species
S = model.Species('S', D1, radius)
P = model.Species('P', D2, radius)
K = model.Species('K', D2, radius)
KS = model.Species('KS', D2, radius)
Sp = model.Species('Sp', D1, radius)
PSp = model.Species('PSp', D2, radius)

m.add_species_type(S) 
m.add_species_type(P) 
m.add_species_type(K) 
m.add_species_type(KS) 
m.add_species_type(Sp) 
m.add_species_type(PSp) 
# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

#s.set_dt_factor(1e-5)
print V, L

print C2N(498e-9, V)

#TODO Correct sequence?
#(id, corner, unit_x, unit_y, length_x, length_y)
plain1 = _gfrd.create_planar_surface('plain1', [0,0,0], [1,0,0], [0,1,0], L, L)
plain2 = _gfrd.create_planar_surface('plain2', [L/2,0,0], [1,0,0], [0,1,0], L, L)


#fracS = fraction_S(N_K, N_P, Keq)
fracS = 1


S_conc = N_S_total / V * 1e3   # in #/m^3

N_S = N_S_total * fracS
N_Sp = N_S_total - N_S

Dtot = D1 + D2

#Dtot_ref = 1e-12

#ka = k_a(kon, k_D(Dtot, sigma))
#ka = 9e9 / N_A / 1e3 # 1/M s -> m^3/s

kD = k_D(Dtot, sigma)
#ka = k_a(kon, kD)
#kon = per_M_to_m3(0.03e9)

ka = 7e-19
kon = k_on(ka, kD)


Keq_S = Keq * S_conc

kcatkoff = Keq_S * kon
koff = kcatkoff * koff_ratio
kcat = kcatkoff - koff

if mode == 'single':
    kcat1 = kcat * float(N_K) / float(N_P)
    koff1 = kcatkoff - kcat1
    kcat2 = kcat
    koff2 = koff
else:
    kcat1 = kcat2 = kcat
    koff1 = koff2 = koff


kd1 = k_d(koff, kon, kD)
kd2 = k_d(koff2, kon, kD)

print 'ka', ka, 'kD', kD, 'kd1', kd1, 'kd2', kd2
print 'kon m^3/s', kon, '1/M s', kon * N_A * 1e3
print 'koff1 1/s ', koff1
print 'kcat1 1/s ', kcat1
print 'koff2 1/s ', koff2
print 'kcat2 1/s ', kcat2

assert koff2 >= 0



print 'S mol conc', S_conc / 1e3 / N_A

print (koff1 + kcat1)/kon/S_conc


#sys.exit(0)

if mode == 'normal' or mode == 'immobile':
    s.throw_in_particles(w, K, N_K)
    s.throw_in_particles(w, P, N_P)
elif mode == 'localized':
    
    P = model.Species('P', D2, radius,'plain1') 
    K = model.Species('K', D2, radius,'plain2') 
    m.add_species_type(P) 
    m.add_species_type(K)

    s.throw_in_particles(K, N_K, plain1)
    s.throw_in_particles(P, N_P, plain2)
elif mode == 'single':
    x = L/2
    yz = L/2
    tl = L/4
    place_particle(w, K, [tl, tl, tl])
    place_particle(w, K, [tl, tl, yz+tl])
    place_particle(w, K, [tl, yz+tl, tl])
    place_particle(w, K, [tl, yz+tl, yz+tl])
    place_particle(w, P, [x+tl, tl, tl])
    place_particle(w, P, [x+tl, tl, yz+tl])
    place_particle(w, P, [x+tl, yz+tl, tl])
    place_particle(w, P, [x+tl, yz+tl, yz+tl])
else:
    assert False

throw_in_particles(w, Sp, N_Sp)
throw_in_particles(w, S, N_S)

# Stir before actually start the sim.

stir_time = 1e-7
while 1:
    s.step()
    next_time = s.get_next_time()
    if next_time > stir_time:
        s.stop(stir_time)
        break

s.reset()

#  1 2 S + K  <-> KS
#  3   KS      -> K + Sp
#  4 5 Sp + P <-> PSp
#  6   PSp     -> P + S




# Create reaction rules (can this be done after rest is already performed?) TODO
r1 = model.create_binding_reaction_rule(S, K, KS, ka)
m.network_rules.add_reaction_rule(r1)
r2 = model.create_unbinding_reaction_rule(KS, S, K, kd1)
m.network_rules.add_reaction_rule(r2)
r3 = model.create_unbinding_reaction_rule(KS, K, Sp, kcat1)
m.network_rules.add_reaction_rule(r3)
r4 = model.create_binding_reaction_rule(Sp, P, PSp, ka)
m.network_rules.add_reaction_rule(r4)
r5 = model.create_unbinding_reaction_rule(PSp, Sp, P, kd2)
m.network_rules.add_reaction_rule(r5)
r6 = model.create_unbinding_reaction_rule(PSp, P, S, kcat2)
m.network_rules.add_reaction_rule(r6)


model = 'pushpull'

# 'pushpull-Keq-koff_ratio-N_K-N_P-V-mode.dat'
l = logger.Logger(logname = model + '_' + '_'.join(sys.argv[1:8]) + '_' # +\
#               os.environ['SGE_TASK_ID']
           ,comment = '@ model=\'%s\'; Keq=%s; koff_ratio=%s\n' %
           (model, Keq_str, koff_ratio_str) +
           '#@ V=%s; N_K=%s; N_P=%s; mode=\'%s\'; T=%s\n' % 
           (V_str, N_K, N_P, mode, T_str) +
           '#@ kon=%g; koff1=%g; koff2=%g; N_S_total=%s\n' %
           (kon, koff1, koff2, N_S_total) +
           '#@ kcat1=%g; kcat2=%g\n' %
           (kcat1, kcat2) +
           '#@ ka=%g; kd1=%g; kd2=%g\n' %
           (ka, kd1, kd2))
interrupter = logger.FixedIntervalInterrupter(s, 1e-7, l.log)

l.start(s)
while s.t < T:
    interrupter.step()

    if s.last_reaction:
        #log.info(dumper.dump_particles(s))
        l.log(s, s.t)
