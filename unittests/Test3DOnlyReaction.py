# Simulation of the reactions:
# A + B -(ka)-> C
# C     -(kd)-> A + B


# ### Importing necessary files ###
from egfrd import *
import model
import gfrdbase
import _gfrd
import os
import sys


# Extra functions
def check_existance( species_id ):
    pids = w.get_particle_ids(species_id)
    if(len(pids) > 0): return True
    else: return False


def give_pid( species_id ):
    pids = w.get_particle_ids(species_id)
    for pid in pids: return pid


def constants():
    '''Define the simulation constants.'''
    global L, sigma, D, ka, kd, N, rnd_seed, tend, M_PI, P_c

    M_PI = numpy.pi  # Pi
    L     = 1e-6     # Lengths of simulation
    V     = L*L*L    # Volume of simulation
    sigma = 4.4e-9   # Diameter particle
    D     = 3e-12    # Diffusion constant (3D)
    ka    = .001 * 4 * M_PI * (2*sigma)**2 # Association constant
    kd    = 10.      # Complex dissociation constant
    N     = 100      # Starting number of particles
    rnd_seed = 0     # seed for randomizer
    tend  = 5        # Simulated time
    c_bar = (N - 1)/V
    P_c = ka*c_bar/(ka*c_bar + kd) #=t_c/t_tot

def setup_model():
    '''Define the model (m) with three particle species (A,B,C)
    and two reaction rules (r1,r2).'''
    global m, A, B, C, r1, r2

    m = model.ParticleModel(L)
    A = model.Species('A', 0.,   sigma )
    B = model.Species('B', D ,   sigma )
    C = model.Species('C', 0., 2*sigma )

    m.add_species_type(A)
    m.add_species_type(B)
    m.add_species_type(C)

    r1 = model.create_binding_reaction_rule(A, B, C, ka)
    m.network_rules.add_reaction_rule(r1)

    r2 = model.create_unbinding_reaction_rule(C, A, B, kd)
    m.network_rules.add_reaction_rule(r2)


def setup_simulator():
    '''SetUp the simulator (s), given the world (w) and the model (m).'''
    global w, s

    w = gfrdbase.create_world(m, 3)

    myrandom.seed( rnd_seed )
    s = EGFRDSimulator(w, myrandom)

    place_particle(w, A, [0.5*L,0.5*L,0.5*L])
    throw_in_particles(w, B, N - 1)



# ### Set everything up: ###
constants()
setup_model()
setup_simulator()


# ### Let's roll: ###
t_current  = 0.0
t_bind     = 0.0
t_decay    = 0.0
oldt_decay = 0.0
Nbinds     = 0
step_cntr  = 0
lastC_pid  = give_pid( A )
CatLastReaction = False;

print "P_c: {0}".format(P_c)

while(t_current < tend):
    s.step()

    #print "time: {0}".format(s.t)

    if len(w.get_particle_ids(C)) is not 0:
        print "step#: %s mln, t = %s, #rebinds = %s, #A = %s, #B = %s" %( step_cntr/1e6, s.t, Nbinds, len(w.get_particle_ids(A)), len(w.get_particle_ids(B)) )

    if step_cntr%1000 == 0:
        print "step#: %s mln, t = %s, #rebinds = %s, #A = %s, #B = %s" %( step_cntr/1e6, s.t, Nbinds, len(w.get_particle_ids(A)), len(w.get_particle_ids(B)) )


    if check_existance(C):
        print "got C particle! time: {0}".format(s.t)

    if s.last_reaction:
        #Check is last_reaction was the formation of a NEW C particle.
        if check_existance( C ) and give_pid( C ) != lastC_pid:
            lastC_pid = give_pid( C )
            t_bind = s.t
            Nbinds += 1

            print "time: {0}".format(s.t)

        #Check if last_reaction was the dissociation of a C particle.
        if CatLastReaction and not check_existance( C ):
            oldt_decay = t_decay
            t_decay = s.t
            outfile_t.write('%.15g %.15g %.15g\n' %(t_bind, t_decay, t_bind - oldt_decay ) )
            print "Bind @: %s, Decay @: %s, Rebind time: %s " %(t_bind, t_decay, t_bind - oldt_decay)
            outfile_t.flush()

        CatLastReaction = check_existance( C )
    step_cntr += 1
    t_current = s.get_next_time()

# Synchronize all domains to tend.
s.stop(tend)


# ### Quitting: ###
print '* Bye.'
