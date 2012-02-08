#!/usr/env python

import math
import sys

import numpy
import scipy


import _gfrd
from utils import *

import os
import logging
import logging.handlers

import myrandom

import model

__all__ = [
    'log',
    'setup_logging',
    'p_free',
    'throw_in_particles',
    'place_particle',
    'NoSpace',
    'create_world',
    'ParticleSimulatorBase',
    'get_surfaces',
    'get_closest_surface',
    'create_network_rules_wrapper',
    'place_particle'
    ]

World = _gfrd.World

log = None

def setup_logging():
    global log 
    log = logging.getLogger('ecell')

    if 'LOGFILE' in os.environ:
        if 'LOGSIZE' in os.environ and int(os.environ['LOGSIZE']) != 0:
            handler = logging.handlers.\
                RotatingFileHandler(os.environ['LOGFILE'], mode='w',
                                    maxBytes=int(os.environ['LOGSIZE']))
        else:
            handler = logging.FileHandler(os.environ['LOGFILE'], 'w', )
            
        if 'LOGLEVEL' in os.environ:
            levelvalue = getattr(logging, os.environ['LOGLEVEL'])
            handler.setLevel(levelvalue)
            log.setLevel(levelvalue)
        else:
            handler.setLevel(logging.INFO)
            log.setLevel(logging.INFO)
    else:
        handler = _gfrd.CppLoggerHandler(_gfrd.Logger.get_logger("ecell"))
        if 'LOGLEVEL' in os.environ:
            levelvalue = getattr(logging, os.environ['LOGLEVEL'])
            handler.logger.manager.level = _gfrd.CppLoggerHandler.translateLevelValue(levelvalue)
            log.setLevel(levelvalue)

    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)

    log.addHandler(handler)

setup_logging()


def p_free(r, t, D):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp(- rsq / Dt4) / math.sqrt(Pi4Dt * Pi4Dt * Pi4Dt)

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    
def get_surfaces(world, pos, ignore=[]):
    
    surface_distances = []

    for surface in world.structures:
        if isinstance(surface, _gfrd.Surface) and surface.id not in ignore:
            pos_transposed = \
                world.cyclic_transpose(pos, surface.shape.position)
            distance = world.distance(surface.shape, pos_transposed)
            surface_distances.append((surface, distance))

    return surface_distances

def get_closest_surface(world, pos, ignore=[]):
    """Return
      - closest surface
      - distance to closest surface
    
    We can not use matrix_space, it would miss a surface if the 
    origin of the surface would not be in the same or neighboring 
    cells as pos."""

    surface_distances = get_surfaces(world, pos, ignore)

    if surface_distances:
        return min(surface_distances, key=lambda surface_dist: surface_dist[1])
    else:
        return None, numpy.inf

def create_world(m, matrix_size=10):
    """Create a world object.
    
    The world object keeps track of the positions of the particles
    and the protective domains during an eGFRD simulation.

    Arguments:
        - m
            a ParticleModel previously created with model.ParticleModel.
        - matrix_size
            the number of cells in the MatrixSpace along the x, y and z 
            axis. Leave it to the default number if you don't know what 
            to put here.

    The simulation cube "world" is divided into (matrix_size x matrix_size 
    x matrix_size) cells. Together these cells form a MatrixSpace. The 
    MatrixSpace keeps track in which cell every particle and protective 
    domain is at a certain point in time. To find the neigherest 
    neighbours of particle, only objects in the same cell and the 26 
    (3x3x3 - 1) neighbouring cells (the simulation cube has periodic
    boundary conditions) have to be taken into account.

    The matrix_size limits the size of the protective domains. If you 
    have fewer particles, you want a smaller matrix_size, such that the 
    protective domains and thus the eGFRD timesteps can be larger. If 
    you have more particles, you want a larger matrix_size, such that 
    finding the neigherest neighbours is faster.
    
    Example. In samples/dimer/dimer.py a matrix_size of
    (N * 6) ** (1. / 3.) is used, where N is the average number of 
    particles in the world.

    """
    m.set_all_repulsive()
    world_region = m.get_structure("world")
    if not isinstance(world_region, _gfrd.CuboidalRegion):
        raise TypeError("the world should be a CuboidalRegion")

    if not numpy.all(world_region.shape.half_extent ==
                     world_region.shape.half_extent[0]):
        raise NotImplementedError("non-cuboidal world is not supported")

    world_size = world_region.shape.half_extent[0] * 2

    world = _gfrd.World(world_size, matrix_size)

    for st in m.species_types:
        try:
            structure = st["structure"]
        except _gfrd.NotFound:
            structure = "world"
        world.add_species(
            _gfrd.SpeciesInfo(st.id, 
                              float(st["D"]), 
                              float(st["radius"]), 
                              structure,
                              float(st["v"])))

    for r in m.structures.itervalues():
        world.add_structure(r)

    world.model = m
    return world

def create_network_rules_wrapper(model):
    return _gfrd.NetworkRulesWrapper(model.network_rules)

def throw_in_particles(world, sid, n):
    """Add n particles of a certain Species to the specified world.

    Arguments:
        - sid
            a Species previously created with the function 
            model.Species.
        - n
            the number of particles to add.

    Make sure to first add the Species to the model with the method
    model.ParticleModel.add_species_type.

    """
    species = world.get_species(sid)
    structure = world.get_structure(species.structure_id)

    if __debug__:
        name = world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        log.info('\n\tthrowing in %s particles of type %s to %s' %
                 (n, name, structure.id))

    # This is a bit messy, but it works.
    i = 0
    while i < int(n):
        position = structure.random_position(myrandom.rng)
        position = apply_boundary(position, world.world_size)

        # Check overlap.
        if not world.check_overlap((position, species.radius)):
            create = True
            # Check if not too close to a neighbouring structures for 
            # particles added to the world, or added to a self-defined 
            # box.
            if isinstance(structure, _gfrd.CuboidalRegion):
                surface, distance = get_closest_surface(world, position, [])
                if(surface and
                   distance < surface.minimal_distance(species.radius)):
                    if __debug__:
                        log.info('\t%d-th particle rejected. Too close to '
                                 'surface. I will keep trying.' % i)
                    create = False
            if create:
                # All checks passed. Create particle.
                p = world.new_particle(sid, position)
                i += 1
                if __debug__:
                    log.info('(%s,\n %s' % (p[0], p[1]))
        elif __debug__:
            log.info('\t%d-th particle rejected. I will keep trying.' % i)

def place_particle(world, sid, position):
    """Place a particle of a certain Species at a specific position in 
    the specified world.

    Arguments:
        - sid
            a Species previously created with the function 
            model.Species.
        - position
            a position vector [x, y, z]. Units: [meters, meters, meters].

    Make sure to first add the Species to the model with the method 
    model.ParticleModel.add_species_type.

    """
    species = world.get_species(sid)
    structure = world.get_structure(species.structure_id)
    radius = species.radius

    if world.check_overlap((position, radius)):
        raise NoSpace, 'overlap check failed'

    # Check if not too close to a neighbouring structures for particles 
    # added to the world, or added to a self-defined box.
    if isinstance(structure, _gfrd.CuboidalRegion):
        surface, distance = get_closest_surface(world, position, [])
        if(surface and
           distance < surface.minimal_distance(species.radius)):
            raise RuntimeError('Placing particle failed: %s %s. '
                               'Too close to surface: %s.' %
                               (sid, position, distance))

    if __debug__:
        species = world.get_species(sid)
        structure = world.get_structure(species.structure_id)
        name = world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        log.info('\n\tplacing particle of type %s to %s at position %s' %
                 (name, structure.id, position))

    particle = world.new_particle(sid, position)
    return particle

class ParticleSimulatorBase(object):
    def __init__(self, world, rng, network_rules):
        self.world = world
        self.rng = rng
        self.network_rules = network_rules

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0

        self.dissociation_retry_moves = 1
        
        self.dt_limit = 1e-3
        self.dt_max = self.dt_limit

        # counters
        self.rejected_moves = 0
        self.reaction_events = 0

        self.max_matrix_size = 0

        self._distance_sq = None
        self._disatnce_sq_array = None

        self.last_reaction = None

    def initialize(self):
        pass

    def get_species(self):
        """
        Return an iterator over the Species in the simulator. 
        To be exact, it returns an iterator that returns 
        all SpeciesInfo instances defined in the simulator.
        
        Arguments: 
            - sim an EGFRDSimulator. 

        More on output:
            - SpeciesInfo is defined in the C++ code, and has
              the following attributes that might be of interest
              to the user:
                   * SpeciesInfo.id
                   * SpeciesInfo.radius
                   * SpeciesInfo.structure_id
                   * SpeciesInfo.D
                   * SpeciesInfo.v
                
        Example:

            for species in s.get_species():
                print str(species.radius)

        Note:
            Dumper.py also holds a copy of this function.        

        """ #TODO: Added by wehrens@amolf.nl; please revise
        return self.world.species

    def get_first_pid(self, sid):
        """
        Returns the first particle ID of a certain species.

        Arguments:
            - sid: a(n) (eGFRD) simulator.
        """ #TODO: Added by wehrens@amolf.nl; please revise
        return iter(self.world.get_particle_ids(sid)).next()

    def get_position(self, object):
        """
        Function that returns particle position. Can take multiple
        sorts of arguments as input. 

        Arguments:
            - object: can be:
                * pid_particle_pair tuple 
                  (of which pid_particle_pair[0] is _gfrd.ParticleID).
                * An instance of _gfrd.ParticleID.

                * An instance of _gfrd.Particle.
                * An instance of _gfrd.SpeciesID.
                    In the two latter cases, the function returns the 
                    first ID of the position of the first particle of 
                    the given species.
        """ #TODO: Added by wehrens@amolf.nl; please revise
        if type(object) is tuple and type(object[0]) is _gfrd.ParticleID:
            pid = object[0]
        elif type(object) is _gfrd.ParticleID:
            pid = object
        elif type(object) is _gfrd.Particle:
            pid = self.get_first_pid(object.sid)
        elif type(object) is _gfrd.SpeciesID:
            pid = self.get_first_pid(object)

        return self.world.get_particle(pid)[1].position


    def clear(self):
        self.dt_max = self.dt_limit
        self.dt = self.dt_limit

    def check_particle_matrix(self):
        total = sum(len(self.world.get_particle_ids(s.id))
                    for s in self.world.species)

        if total != self.world.num_particles:
            raise RuntimeError('total number of particles %d != '
                               'self.world.num_particles %d' %
                               (total, self.world.num_particles))

    def check_particles(self):
        for pid_particle_pair in self.world:
            pid = pid_particle_pair[0]
            pos = pid_particle_pair[1].position
            if (pos >= self.world.world_size).any() or (pos < 0.0).any():
                raise RuntimeError('%s at position %s out of the world '
                                   '(world size=%g).' %
                                   (pid, pos, self.world.world_size))


    def check(self):
        self.check_particle_matrix()
        self.check_particles()

    def print_report(self):
        pass

class NoSpace(Exception):
    pass

