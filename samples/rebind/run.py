#!/usr/bin/env python

# See README for a description

import sys

from egfrd import *
from bd import *
import gfrdbase
import model

def run(outfilename, N_B, N_X, N, T_list):
    """ 
    Function that loops over single runs of binding/rebinding
    simulations. (Plus some overhead code.)

    Arguments:
        - outfilename:  This string can be found back in the
                        output filename.
        - T_list:       List of times at which distance between 
                        particles is evaluated and logged.
        - N_B:          Number of instances ("molecules") of 
                        particle B.
        - N_X:          Number of instances ("molecules") of 
                        particle X.

    """ 
    global DX_factor, V, L, matrix_size, D, radius, DX, kf, tau   

    print "Printing output to: " + outfilename

    ################# Set some constants 
    radius = 2.5e-9
    sigma = radius * 2
    D = 1e-12
    D_tot = D * 2

    tau = sigma**2 / D_tot

    DX_factor = 1

    V = 1e-18 # m^3
    L = V ** (1.0/3.0) 

    kf = 0.092e-18

    DX = D * DX_factor

    matrix_size = min(max(3, int((9 * (N_X+N_B)) ** (1.0/3.0))), 60)
    print 'matrix_size=', matrix_size

    # Clear the distance output files
    for j in T_list:
        outfile_r = open(outfilename + '_r_'+str(j)+'.dat', 'a')
        outfile_r.close()

    # Open the time output file
    outfile_t = open(outfilename + '_t.dat', 'w')

    ################# Loop through the multiple (N) simulations runs
    for i in range(N):

        # Perform a single simulation run
        r_list, t_list = singlerun(T_list, N_B, N_X)

        # Write away time data
        for t in t_list:
            outfile_t.write('%g\n' % t)
    
        # Write away distance data, write data of different 
        # measurement points to different files. (Ending up
        # with a list of distances corresponding to measure-
        # ments at times from T_list.)
        for j in range(len(r_list)):
            outfile_r = open(outfilename + '_r_-' + str(T_list[j]) +'.dat', 'a')
            outfile_r.write(str(r_list[j]) + '\n')
            outfile_r.close()

        # Make sure script writes data NOW
        outfile_t.flush()

    # Close time output file
    outfile_t.close()

def singlerun(T_list, N_B, N_X):
    """ 
    Function that performs one simulation run of binding/
    rebinding. 

    Arguments:
        - T_list:       List of times at which distance between 
                        particles is evaluated and logged.
        - N_B:          Number of instances ("molecules") of 
                        particle B.
        - N_X:          Number of instances ("molecules") of 
                        particle X.
    """ 

    ################# Define the model:
    
    # Create model class
    m = model.ParticleModel(L)

    # Particle species classes
    A = model.Species('A', D, radius)
    m.add_species_type(A)
    B = model.Species('B', D, radius)
    m.add_species_type(B)
    C = model.Species('C', D, radius)
    m.add_species_type(C)
    X = model.Species('X', DX, radius)
    m.add_species_type(X)
    
    # Reaction rules
    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    r2 = model.create_unbinding_reaction_rule(C, A, B, 1e3)
    m.network_rules.add_reaction_rule(r2)

    ################# Set up simulator

    w = gfrdbase.create_world(m, matrix_size)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = EGFRDSimulator(w, myrandom.rng, nrw)

    ################# Place particles 

    # Throw in some X particles and stirr
    if N_X != 0:
        gfrdbase.throw_in_particles(w, X, N_X)

        end_time = tau * 1
        while 1:
            s.step()
            next_time = s.get_next_time()
            if next_time > end_time:
                s.stop(end_time)
                break

    # Reset simulation, so stirring has no effect on output data
    s.reset()

    # Define positions for particles to place.
    A_pos = [0,0,0]
    B_pos = [(float(A['radius']) + float(B['radius']))+1e-23,0,0]

    # Clear an area at position A_pos and place particle A there
    # 
    # How this should be done:
    #    - find particles at position A_pos within species A radius
    #    - delete those
    #    - add a number of X particles _randomly to the box_, the 
    #      number being equal to # removed particles.
    #    - if a particle is accidently placed within the "cleared"
    #      area, repeat the process (in very crowded box, this 
    #      leads to infinite loop)
    #
    # Currently, however, it doesn't get the particles within a certain radius, 
    # but the particles which domain lies within a certain radius. 
    #
    # Therefore, you might be removing particles that after bursting are not 
    # within the radius anymore, but whose (now bursted) domains were.
    # 
    # TODO
    # This code should thus be adapted to solve this problem.
    # 
    # Perhaps it would be more easy to just place all the particles at
    # the beginning, but give A and B initially a diffusion constant of
    # 0, then stirr, and then set the diffusion constant to their
    # desired values.."""

    while 1:

        # Burst domains within a certain volume and collect their ids
        dd = s.clear_volume(A_pos, float(A['radius'])) 

        if not dd:
            break
        for d in dd:
            for p in d.particle_id_pair:
                s.remove_particle(p)
        s.throw_in_particles(X, len(pp), box1)

    place_particle(w, A, A_pos)

    # Idem for a particle B:
    while 1:
    
        # Burst domains within a certain volume and collect their ids
        dd = s.clear_volume(B_pos, float(B['radius'])) 

        if not dd:
            break
        for d in dd:
            for p in d.particle_id_pair:
                s.remove_particle(p)
        s.throw_in_particles(X, len(pp), box1)

    place_particle(w, B, B_pos) 

    # Place rest of B at random positions
    if N_B > 1:
        gfrdbase.throw_in_particles(w, B, N_B-1)

    # Create some vars and perform first simulation step
    r_list = []
    t_list = []
    t_last = 0

    s.step()

    next_stop = T_list[0] 

    i_T = 0

    #  ### Start simulating
    while 1:

        # What happens in the following if-statement: 
        # Create a list of the durations particles were in bound state.
        #
        # If there was a reaction:
        #     - Binding: Indicated by the fact that there are no C-particles.
        #       In this case: record the current time (s.t) to t_last.
        #     - Unbinding: Indicated because this is the only other case when
        #       a reaction has taken place. 
        #       In this case: record the time binding lasted, i.e. dt between
        #       current time (s.t) and last reaction time (t_last).
        if s.last_reaction: # aka if there was a reaction

            if len(s.world.get_particle_ids(C.id)) == 0:  #A,B
                print '    - set t_last', str(s.t)
                # set t_last to "current time" in simulator
                t_last = s.t  
            else:    
                print '    - reaction: ', str(s.t - t_last)
                t_list.append(s.t - t_last)

        # If it's time to log, then log distance between particles. 
        next_time = s.get_next_time()
        if next_time > next_stop:
            print '* Measuring distances at ', str(next_stop),\
                            '(measurement #', str(i_T), ')'
            s.stop(next_stop)
            print '* stopped'
            # If there is a particle C, the distance is "0".
            if len(s.world.get_particle_ids(C.id)) != 0:  #A,B
                r_list.append(0)
            # If particles are manifested as A and B, log distance
            # inbetween.
            else: # C = 0, i.e. there are particles A and B

                # Calculate distances:
                A_pos = s.get_position(A.id)
                distance_list = []
                # If there are >1 B particles, we want to log all   
                # distances seperately, hence the loop:
                for particle_id in w.get_particle_ids(B.id):
                    B_pos = s.get_position(particle_id)
                    distance = w.distance(A_pos, B_pos)
                    distance_list.append(distance)

                # Write list to distances table
                r_list.append(distance_list)

            i_T += 1
            next_stop = T_list[i_T]
        
        # If there are no measuring moments on the list any more, and
        # a duration of being bounded has been recorded, then stop
        if (next_stop == INF) and (len(t_list) != 0):
            print 'break', s.t
            break

        s.step()

    return r_list, t_list
    

if __name__ == '__main__':

    import os
    
    # convert last user arguments to array T_list
    T_list = []
    for measurement_time in sys.argv[4:]:    
        T_list.append(float(measurement_time))

    run('data/rebind', float(sys.argv[1]), 
        int(sys.argv[2]), int(sys.argv[3]), T_list)

    print "Done."











