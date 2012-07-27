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
    'create_box',
    'ParticleSimulatorBase',
    'get_all_surfaces',
    'get_neighbor_surfaces',
    'create_network_rules_wrapper',
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

def get_neighbor_surfaces(world, pos, current_struct_id, ignores=[]):
    
    if ignores:
        ignore = ignores[0]
    else:
        ignore = world.get_def_structure_id()   # smt random

    surface_distances = [(structure, distance) for ((id, structure), distance) in world.get_close_structures(pos, current_struct_id, ignore)]

    return surface_distances

def get_all_surfaces(world, pos, ignores=[]):

    surface_distances = []

    for surface in world.structures:
        if isinstance(surface, _gfrd.Surface) and surface.id not in ignores:
            pos_transposed = \
                world.cyclic_transpose(pos, surface.shape.position)
            distance = world.distance(surface.shape, pos_transposed)
            surface_distances.append((surface, distance))

    return sorted(surface_distances, key=lambda surface_dist: surface_dist[1])

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

    # make the world
    world = _gfrd.World(m.world_size, matrix_size)

    # copy all structure_types to the world (including the default one) and set the default id.
    for st in m.structure_types:
        world.add_structure_type(st)
    world_structure_type_id = m.get_def_structure_type_id()
    world.set_def_structure_type_id(world_structure_type_id)

    # The model (a ParticleModel) now only holds SpeciesTypes
    # The SpeciesTYpes hold the information that is required for spatial simulations
    # in auxiliairy string fields. This information is here converted to the right
    # SpeciesInfo object and added to the world.
    for st in m.species_types:
        if st["structure_type"] == "world":
            # The default structure_type
            structure_type_id = world_structure_type_id
        else:
            # Find the corresponding structure_type between all the structure_types known
            # FIXME: ugly hack
            for stt in m.structure_types:
                if st["structure_type"] == stt['name']:
                    structure_type_id = stt.id
                    break
            else:
                raise RuntimeError("StructureType used in Species not found in Model when initializing World.")

        world.add_species(
            _gfrd.SpeciesInfo(st.id, 
                              structure_type_id,
                              float(st["D"]), 
                              float(st["radius"]), 
                              float(st["v"])))

    # make the default structure: A cuboidal region of the default structure_type
    # This needs to be done after making the structure_types or add_structure can't
    # find the proper structure_type.
    x = numpy.repeat(m.world_size / 2, 3)
    region = _gfrd.CuboidalRegion("world", world_structure_type_id, world.get_def_structure_id(), _gfrd.Box(x, x))
    world.set_def_structure(region)

    world.model = m
    return world

def create_box(world, structure_type, center, size):
    # Creates a box of PlanarSurface sections and adds it to the world.
    # size is a vector in x, y, z

    # assert that the center and size is ok
    center = numpy.array(center)
    size   = numpy.array(size)
    assert all(0 < center) and all(center < world.world_size)
    assert all(0 < size)   and all(size < world.world_size)
    assert all(0 < center - size/2.0) and all(center + size/2.0 < world.world_size)

    sid = structure_type.id
    name = 'box'
    def_struct_id = world.get_def_structure_id()
    
    front  = model.create_planar_surface(sid, name+'_front', [center[0] - size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2], [0, 0, 1], [1, 0, 0], size[2], size[0], def_struct_id)
    back   = model.create_planar_surface(sid, name+'_back',  [center[0] - size[0]/2, center[1] + size[1]/2, center[2] - size[2]/2], [1, 0, 0], [0, 0, 1], size[0], size[2], def_struct_id)
    right  = model.create_planar_surface(sid, name+'_right', [center[0] + size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2], [0, 0, 1], [0, 1, 0], size[2], size[1], def_struct_id)
    left   = model.create_planar_surface(sid, name+'_left',  [center[0] - size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2], [0, 1, 0], [0, 0, 1], size[1], size[2], def_struct_id)
    top    = model.create_planar_surface(sid, name+'_top',   [center[0] - size[0]/2, center[1] - size[1]/2, center[2] + size[2]/2], [0, 1, 0], [1, 0, 0], size[1], size[0], def_struct_id)
    bottom = model.create_planar_surface(sid, name+'_bottom',[center[0] - size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2], [1, 0, 0], [0, 1, 0], size[0], size[1], def_struct_id)
    world.add_structure(front)
    world.add_structure(back)
    world.add_structure(right)
    world.add_structure(left)
    world.add_structure(top)
    world.add_structure(bottom)
    world.connect_structures( front, 3, top, 2)
    world.connect_structures( front, 2, bottom, 1)
    world.connect_structures( front, 1, left, 2)
    world.connect_structures( front, 0, right, 1)
    world.connect_structures(  back, 3, right, 0)
    world.connect_structures(  back, 2, left, 3)
    world.connect_structures(  back, 1, bottom, 0)
    world.connect_structures(  back, 0, top, 3)
    world.connect_structures(   top, 1, left, 0)
    world.connect_structures(   top, 0, right, 3)
    world.connect_structures(bottom, 3, right, 2)
    world.connect_structures(bottom, 2, left, 1)


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
    structure_type = world.get_structure_type(species.structure_type_id)
    structure_list = list(world.get_structure_ids(structure_type))

    if __debug__:
        name = world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        log.info('\n\tthrowing in %s particles of type %s to %s' %
                 (n, name, structure_type.id))

    i = 0
    while i < int(n):
        # This is messy, but it works.
        myrandom.shuffle(structure_list)
        structure = world.get_structure(structure_list[0])
        position = structure.random_position(myrandom.rng)
        position, structure_id = world.apply_boundary((position, structure.id))

        # Check overlap. TODO put in 'if' statement for improved efficiency?
        particle_overlaps = world.check_overlap((position, species.radius))
        surface_overlaps  = world.check_surface_overlap((position, species.radius*MINIMAL_SEPARATION_FACTOR),
                                                        position, structure_id, species.radius)

        if (not particle_overlaps) and (not surface_overlaps):
            # All checks passed. Create particle.
            p = world.new_particle(sid, structure_id, position)
            i += 1
            if __debug__:
                log.info('particle accepted: (%s,\n %s)' % (p[0], p[1]))
        elif __debug__:
            if particle_overlaps:
                log.info('\t%d-th particle rejected. Too close to particle. I will keep trying.' % i)
            if surface_overlaps:
                log.info('\t%d-th particle rejected. Too close to surface. I will keep trying.' % i)

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
    radius = species.radius

    # FIXME This is a mess!!

    # check that particle doesn't overlap with other particles.
    if world.check_overlap((position, radius)):
        raise NoSpace, 'Placing particle failed: overlaps with other particle.'

    surfaces = get_all_surfaces(world, position, [])
    if surfaces:
        surface, distance = surfaces[0]
    else:
        surface, distance = None, numpy.inf

    # Check if not too close to a neighbouring structures for particles 
    # added to the world, or added to a self-defined box.
    if species.structure_type_id == world.get_def_structure_type_id():
        structure_id = world.get_def_structure_id()
        if world.check_surface_overlap((position, radius*MINIMAL_SEPARATION_FACTOR),
                                       position, structure_id, radius):#(surface and distance < surface.minimal_distance(species.radius)):
            raise RuntimeError('Placing particle failed: %s %s. '
                               'Too close to surface: %s.' %
                               (sid, position, distance))
    else:
        # If the particle lives on a surface then the position should be in the closest surface.
        # The closest surface should also be of the structure_type associated with the species.
        structure_ids = world.get_structure_ids(world.get_structure_type(species.structure_type_id))
        if not (surface and 
                distance < TOLERANCE*radius and 
                surface.id in structure_ids):
            raise RuntimeError('Placing particle failed: %s %s. Position should be in structure of structure_type \"%s\".' %
                                (sid, position, world.get_structure_type(species.structure_type_id)['name']))
        else:
            structure_id = surface.id

    if __debug__:
        name = world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        log.info('\n\tplacing particle of type %s to structure \"%s\" at position %s' %
                 (name, world.get_structure(structure_id).name, position))

    position, structure_id = world.apply_boundary((position, structure_id))
    particle = world.new_particle(sid, structure_id, position)

    if __debug__:
        log.info('(%s,\n %s' % (particle[0], particle[1]))

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
                   * SpeciesInfo.structure_type_id
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

