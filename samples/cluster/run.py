
# Rebinding to a cluster

# Run by:
# $ python run.py [N] [runs] [outFilename] [Logmode, default=False]
#
# Arguments:
#   - N: Number of particles in cluster
#   - runs: Number of simulation runs
#   - outFilename: Name of output file
#   - Logmode: false by default, if True, only 1 VTK-logged run is
#     performed.    
#
# E.g.:
# $ python run.py 7 1000 data.out
# Or:
# $ python run.py 7 1 data.out True

# Modules
# ===============================
import sys
#Also, set relative egfrd directory path
sys.path.append('../../')
import os
import shutil
import datetime

import math
from egfrd import *
import model
import gfrdbase
import _gfrd
from visualization import vtklogger

if __name__ == "__main__":
    # Constants
    # ===============================
    # Number of particles in the cluster
    N = int(sys.argv[1])
    # Number of runs
    runs = int(sys.argv[2])
    # Output file
    outFilename = sys.argv[3]
    # LOGGING mode
    try: 
        LOGGING = bool(sys.argv[4])
        if (LOGGING==True): 
            print "* Performing only 1 logging run."
            runs = 1
    except: LOGGING = False

    # Particle constants
    sigma = 1e-5        # Diameter particle; big:1e-5
    D = 1e-8            # Diffusion constant; big:1e-8
    world_size = 1e-3   # Lengths of simulation box; normal: 1e-3

    k1 = 1e-10
    k2 = 1e2

    # Spacing inbetween cluster particles AND cluster/B-particle
    spacing = sigma/1e5

    # Create "unique" seed
    # ===============================
    currenttime = (long(datetime.datetime.now().year*3600*24*365+
    datetime.datetime.now().month*30.5*24*3600+
    datetime.datetime.now().day*24*3600+datetime.datetime.now().hour*3600+
    datetime.datetime.now().minute*60+datetime.datetime.now().second))
    myrandom.seed(currenttime)

    print str('Seed: '+str(currenttime))


# Functions
# ===============================
def cartesian(g1, g2, ng1, ng2):
    """ Converts lattice coordinates (ng1, ng2) to cartesian coordinates
    
    Arguments:
        - Lattice vectors g1=[x1,y1] and g2=[x2,y2]
        - Lattice coordinates (ng1, ng2)

    """
    cartesian_x = (ng1 * (g1[0]) + ng2 * (g2[0]))
    cartesian_y = (ng2 * (g2[1]) + ng1 * (g1[1]))
    return cartesian_x, cartesian_y

def distance(x1, y1, x2, y2):
    """ Calculates the distance between points (x1, y1) and (x2, y2) """
    return math.sqrt(math.pow((x1-x2),2)+math.pow((y1-y2),2))

def generate_possible_positions(N):
    """ Generates a list with possible coordinates and distances to center

    This function does the following:
        - It loops over square area of an hexagonal lattice (of size 
          calculated to facilitate N particles) 
        - Lattice coordinates and distance to the center of the area 
          are stored if they are within circular boundaries. (The 
          latter check is only done to be absolutely sure a circular 
          shape is formed.)

    This information can be used in a later function. If the list of
    generated coordinates is sorted by distance to the center, N 
    particles can be placed within a circular geometry by looping over
    the coordinates.

    Arguments:
        - N: number of particles to place.

    """
    global w, A, spacing, sigma #, g1, g2, cartesian_x, cartesian_y, ng1, ng2

    # Lattice vectors
    d1 = [1,0]
    d2 = [math.cos((math.pi)/3.0),math.sin((math.pi)/3.0)]
    # Scaled lattice vectors
    lengthLatticeVector = (sigma+spacing)
    g1 = [lengthLatticeVector*d1[0], lengthLatticeVector*d1[1]]
    g2 = [lengthLatticeVector*d2[0], lengthLatticeVector*d2[1]]

    # Factors to calculate square area of lattice needed
    """
    factorNonSquareGrid = math.sqrt(1/(d1[0]*d2[1]-d1[1]*d2[0]))
    factorSquareToCircle = math.pi/2 #math.sqrt(math.pow((math.pi*0.5),2)) 
    # Calculate diameter of sphere
    L = (spacing+sigma)*(math.ceil(math.sqrt(N))+2)*factorNonSquareGrid*factorSquareToCircle
        # "(math.ceil(math.sqrt(N))+2)"; ceil and +2: some margin
    diamondSide = (L/2)*math.tan(pi/6)+(L/2)/math.tan(pi/6)
    diamondDiagonal = sin(pi/3)*L #2*sin(pi/3)*L/2
    """
    # d1x*d2y-d1y*d2x
    diamondToSquare = math.sqrt((4*d2[0]*d2[1])/(0.25*math.pi))
    vectorsOnDiameter = math.ceil(math.sqrt(N)+2)
    numberOfVectors = diamondToSquare * vectorsOnDiameter

    diameterCircle = math.ceil(math.sqrt(N)+2)*lengthLatticeVector
    radiusCircle = diameterCircle/2

    # Cartesian coordinates
    cartesian_x = 0
    cartesian_y = 0
    cartesian_z = 0

    # The list with possible coordinates and distances to center
    clusterParticleCoordinates = []

    for ng1 in range(int(-math.ceil(numberOfVectors/2)), 
                            int(math.ceil(numberOfVectors/2))):
        for ng2 in range(int(-math.ceil(numberOfVectors/2)), 
                                int(math.ceil(numberOfVectors/2))):
            cartesian_x, cartesian_y = cartesian(g1, g2, ng1, ng2)
            distanceToCenter = distance(cartesian_x, cartesian_y, 0, 0)
            if (distanceToCenter < radiusCircle):
                clusterParticleCoordinates.append([
                            distanceToCenter,
                            cartesian_x, 
                            cartesian_y, 
                            cartesian_z
                            ])
    

    return clusterParticleCoordinates

def make_cluster(N, coord_x, coord_y, coord_z): 
    """ Places N particles on hexagonal lattice in spherical symmetry

    Arguments:
        - N: number of particles to place
        - coord_x, coord_y, coord_z: Coordinates where to place cluster
    """
    # Generate list of positions and there distance to the center
    clusterParticleCoordinates = generate_possible_positions(N)
    # Sort this list
    clusterParticleCoordinates.sort(key=lambda x: x[0])

    # Place center particle
    place_particle(w, C,
         [clusterParticleCoordinates[0][1]+coord_x, 
          clusterParticleCoordinates[0][2]+coord_y,
          clusterParticleCoordinates[0][3]+coord_z])
    # Place N particles starting with the particle closest to the
    # center, then placing the second closest particle, etc..
    for i in range(1,N):
        try:
            place_particle(w, A,
                [clusterParticleCoordinates[i][1]+coord_x, 
                 clusterParticleCoordinates[i][2]+coord_y,
                 clusterParticleCoordinates[i][3]+coord_z])
        except: 
            print "ERROR: couldn't place particle."

def single_run(N, LOGGING):
    """ Single run of simulation """
    global w, A, C, k1, k2
    # Basic set up simulator
    # ===============================
    # Model
    m = model.ParticleModel(world_size)
    # Species
    A = model.Species('A', 0, sigma/2)                                   
    m.add_species_type(A) 
    B = model.Species('B', D, sigma/2)                                   
    m.add_species_type(B) 
    C = model.Species('C', 0, sigma/2)                                   
    m.add_species_type(C) 
    # Reaction rules
    r1 = model.create_binding_reaction_rule(A, B, C, k1)
    m.network_rules.add_reaction_rule(r1)
    r2 = model.create_unbinding_reaction_rule(C, A, B, k2)
    m.network_rules.add_reaction_rule(r2)

    # World
    w = gfrdbase.create_world(m, 3)
    # Simulator   
    s = EGFRDSimulator(w, myrandom.rng)

    # Put in cluster
    make_cluster(N, world_size/2, world_size/2, world_size/2)

    # Put in reactants
    # place_particle(w, B, [world_size/2, world_size/2, world_size/2+sigma+spacing])    

    # Enable VTK Logger
    # ===============================
    if (LOGGING == True):
        vtk_output_directory = 'VTK_out'
        if (os.path.exists(vtk_output_directory)):
            print '** Warning: VTK directory already exists.'
        l = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 
    
    # Running 
    # ===============================
    numberDetected = 0

    if (LOGGING == True):
        while 1:
            l.log() # log
            s.step() # and make eGFRD step
            if s.last_reaction:
                numberDetected = numberDetected+1
                if (numberDetected == 2):
                    # print "2nd Reaction detected at: " + str(s.t) + "(" + str(s.last_reaction) + ")"
                    reaction_time = s.t - previous_time
                    break  
                else: previous_time = s.t
        l.stop() 
    else:
        while 1:
            s.step() # make eGFRD step
            if s.last_reaction:
                numberDetected = numberDetected+1
                if (numberDetected == 2):
                    # print "2nd Reaction detected at: " + str(s.t) + "(" + str(s.last_reaction) + ")"
                    reaction_time = s.t - previous_time
                    break  
                else: previous_time = s.t

    s.stop(s.t)
    return (reaction_time)
    
# Main part
# ===============================

if __name__ == "__main__":
    # Output file
    outFile = open(outFilename, 'w')
    outFile.close()
    outFile = open(outFilename, 'a')

    for M in range(runs):
        outFile.write(str(single_run(N, LOGGING)) + '\n')
        outFile.flush    

    outFile.close()

    timetaken = (long(datetime.datetime.now().year*3600*24*365+
    datetime.datetime.now().month*30.5*24*3600+
    datetime.datetime.now().day*24*3600+datetime.datetime.now().hour*3600+
    datetime.datetime.now().minute*60+datetime.datetime.now().second)-currenttime)
    print "Done in "+ str(timetaken) +" seconds."


# #### OLD CODE

#def in_circle(cartesian_x, cartesian_y, r_x, r_y, R):
#    if (math.sqrt(math.pow((cartesian_x-r_x),2)+math.pow((cartesian_y-r_y),2)) < R):
#        return True
#    else:
#        return False


    # while end time hasn't come and #timesteps still below limit
    # while ((s.get_next_time() < endTime) and (n < EARLY_STOP)): # TODO


    # Terminate the simulation at appropriate time

    #"""
    #if (n >= EARLY_STOP):
    #    print "Early stop (timestep limit)."
    #    s.stop(s.t)
    #    reaction_time = INF
    #else:
    #    print "Early stop (time limit)."
    #    s.stop(s.get_next_time())
    #    reaction_time = INF
    #"""

    # Make sure logger writes away information.
    # l.stop() #TODO REMOVED LOGGING












