#!/usr/bin/env python

#
# Make sure that the egfrd system is added to your PYTHONPATH
# This means, in bash for example:
# $ export PYTHONPATH=$HOME/egfrd
#

'''
LOGLEVEL=ERROR python -O run.py rev.3.out 1.25e-2 5000000 &
LOGLEVEL=ERROR python -O run.py rev.2.out 1.25e-3 4000000 &
LOGLEVEL=ERROR python -O run.py rev.1.out 1.25e-4 2000000 &
LOGLEVEL=ERROR python -O run.py rev.0.out 1.25e-5 2000000 &
LOGLEVEL=ERROR python -O run.py rev.-1.out 1.25e-6 2000000 &
LOGLEVEL=ERROR python -O run.py rev.-2.out 1.25e-7 2000000 &
LOGLEVEL=ERROR python -O run.py rev.-3.out 1.25e-8 1000000 &
'''

import sys
from egfrd import *
from bd import *
import model
import gfrdbase
import _gfrd

def run(outfilename, T, N):
    print outfilename

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun(T)
        outfile.write('%.18g\n' % d)
        outfile.flush()
        #print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun(T):

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    D_tot = D * 2

    tau = sigma * sigma / D_tot

    kf = 100 * sigma * D_tot
    koff = 0.1 / tau

    m = model.ParticleModel(1e-3)

    A = model.Species('A', D, sigma/2)
    B = model.Species('B', D, sigma/2)
    C = model.Species('C', D, sigma/2)

    m.add_species_type(A)
    m.add_species_type(B)
    m.add_species_type(C)

    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    r2 = model.create_unbinding_reaction_rule(C, A, B, koff)
    m.network_rules.add_reaction_rule(r2)

    w = gfrdbase.create_world(m, 3)
    nrw = _gfrd.NetworkRulesWrapper(m.network_rules)
    s = EGFRDSimulator(w, myrandom.rng, nrw)

    place_particle(w, A, [0,0,0])
    place_particle(w, B, [(float(A['radius']) + float(B['radius']))+1e-13,0,0])

    end_time = T
    s.step()

    while 1:
        next_time = s.get_next_time()
        if next_time > end_time:
            s.stop(end_time)
            break
        s.step()

    
    if len(s.world.get_particle_ids(C.id)) != 0:
        return 0, s.t

    distance = w.distance(s.get_position(A.id), s.get_position(B.id))

    return distance, s.t
    
if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))
