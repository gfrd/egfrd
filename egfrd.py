#!/usr/env python

from weakref import ref
import math

import numpy

from _gfrd import (
    Event,
    EventScheduler,
    Particle,
    SphericalShell,
    CylindricalShell,
    DomainIDGenerator,
    ShellIDGenerator,
    DomainID,
    ParticleContainer,
    CuboidalRegion,
    CylindricalSurface,
    PlanarSurface,
    Surface,
    _random_vector,
    Cylinder,
    Sphere,
    NetworkRulesWrapper,
    )

from gfrdbase import *
from single import *
from pair import *
from multi import *
from utils import *
from constants import *
from shellcontainer import ShellContainer

import logging
import os

from bd import DEFAULT_DT_FACTOR

log = logging.getLogger('ecell')

if __debug__:
    PRECISION = 3
    FORMAT_DOUBLE = '%.' + str(PRECISION) + 'g'

def create_default_single(domain_id, pid_particle_pair, shell_id, rt, structure):
    if isinstance(structure, CuboidalRegion):
        return SphericalSingle(domain_id, pid_particle_pair,
                               shell_id, rt, structure)
    elif isinstance(structure, CylindricalSurface):
        return CylindricalSurfaceSingle(domain_id, pid_particle_pair, 
                                        shell_id, rt, structure)
    elif isinstance(structure, PlanarSurface):
        return PlanarSurfaceSingle(domain_id, pid_particle_pair, 
                                   shell_id, rt, structure)


def create_default_pair(domain_id, com, single1, single2, shell_id, 
                        r0, shell_size, rt, structure):
    if isinstance(structure, CuboidalRegion):
        return SphericalPair(domain_id, com, single1, single2,
                             shell_id, r0, shell_size, rt, structure)
    elif isinstance(structure, CylindricalSurface):
        return CylindricalSurfacePair(domain_id, com, single1, single2,
                                      shell_id, r0, shell_size, rt, structure)
    elif isinstance(structure, PlanarSurface):
        return PlanarSurfacePair(domain_id, com, single1, single2,
                                 shell_id, r0, shell_size, rt, structure)


class DomainEvent(Event):
    __slot__ = ['data']
    def __init__(self, time, domain):
        Event.__init__(self, time)
        self.data = domain.domain_id		# Store the domain_id key refering to the domain
						# in domains{} in the scheduler

class Delegate(object):
    def __init__(self, obj, method, arg):
        self.ref = ref(obj)
        self.method = method
        self.arg = arg

    def __call__(self):
        return self.method(self.ref(), self.arg)


class EGFRDSimulator(ParticleSimulatorBase):
    """The eGFRDSimulator implements the asynchronous egfrd scheme of performing
    the diffusing and reaction of n particles. The eGFRDsimulator acts on a 'world'
    object containing particles and structures, and can be attached and detached from
    this 'world'.
    """
    def __init__(self, world, rng=myrandom.rng, network_rules=None):
        """Create a new EGFRDSimulator.

        Arguments:
            - world
                a world object created with the function 
                gfrdbase.create_world.
            - rng
                a random number generator. By default myrandom.rng is 
                used, which uses Mersenne Twister from the GSL library.
                You can set the seed of it with the function 
                myrandom.seed.
            - network_rules
                you don't need to use this, for backward compatibility only.

        """
        if network_rules == None:
            network_rules = NetworkRulesWrapper(world.model.network_rules)
        ParticleSimulatorBase.__init__(self, world, rng, network_rules)

        self.domain_id_generator = DomainIDGenerator(0)
        self.shell_id_generator = ShellIDGenerator(0)

	# some constants
        self.MULTI_SHELL_FACTOR = 1.05		# This is the threshold for when the algorithm switches from
						# NonInteractionSingles to a Multi and also defines the Multi
						# shell size.
        self.SINGLE_SHELL_FACTOR = 2.0		# This is the threshold for when the algorithm switches from
						# NonInteractionSingles to a Pair or Interaction. It also defines
						# the radius in which the NonInteractionSingle will burst.
	self.MAX_NUM_DT0_STEPS = 10000

	# used datastructrures
        self.scheduler = EventScheduler()	# contains the events. Note that every domains has exactly one event

        self.domains = {}			# a dictionary containing references to the domains that are defined
						# in the simulation system. The id of the domain (domain_id) is the key.
						# The domains can be a single, pair or multi of any type.

	# other stuff
        self.is_dirty = True			# The simulator is dirty if the state if the simulator is not
						# consistent with the content of the world that it represents
						# (or if we don't know for sure)

        self.reset()				# The Simulator is only initialized at the first step, allowing
						# modifications to the world to be made before simulation starts

    def get_next_time(self):
        """ 
        Returns the time it will be when the next egfrd timestep
        is completed.        
        """ #~MW
        if self.scheduler.size == 0:
            return self.t			# self.t is the current time of the simulator
	else:
            return self.scheduler.top[1].time

    ###
    # The next methods control the general functioning of the simulator
    ###
    def reset(self):
        """
        This function resets the "records" of the simulator. This means
        the simulator time is reset, the step counter is reset, events
        are reset, etc.
        Can be for example usefull when users want to first do an equilibration
	run before starting the "real experiment".
        """ #~ MW
        self.t = 0.0
        self.dt = 0.0
        self.step_counter = 0
        self.single_steps = {EventType.SINGLE_ESCAPE:0,
                             EventType.SINGLE_REACTION:0}
        self.interaction_steps = {EventType.IV_INTERACTION:0,
                                  EventType.IV_ESCAPE:0}
        self.pair_steps = {EventType.SINGLE_REACTION:0,
                           EventType.IV_REACTION:0,
                           EventType.IV_ESCAPE:0,
                           EventType.COM_ESCAPE:0}
        self.multi_steps = {EventType.MULTI_ESCAPE:0,
                            EventType.MULTI_UNIMOLECULAR_REACTION:0,
                            EventType.MULTI_BIMOLECULAR_REACTION:0, 3:0}
        self.zero_steps = 0
        self.rejected_moves = 0
        self.reaction_events = 0
        self.last_event = None
        self.last_reaction = None

        self.is_dirty = True		# simulator needs to be re-initialized

    def initialize(self):
	"""Initialize the eGFRD simulator

	This method (re-)initializes the simulator with the current state of the 'world'
	that it represents.
	Call this method if you change the 'world' using methods outside of the Simulator
	thereby invalidating the state of the eGFRD simulator (the simulator is dirty), or
	in an other case that the state of the simulator is inconsistent with the state of
	the world.
	"""
	# 1. (re-)initialize previously set datastructures
        ParticleSimulatorBase.initialize(self)
        self.scheduler.clear()
        self.domains = {}

	# create/clear other datastructures
	self.geometrycontainer = ShellContainer(self.world)

	# 2. Couple all the particles in 'world' to a new 'single' in the eGFRD simulator
        # Fix order of adding particles (always, or at least in debug mode).
        singles = []
        pid_particle_pairs = list(self.world)
        pid_particle_pairs.sort()

        for pid_particle_pair in pid_particle_pairs:
            single = self.create_single(pid_particle_pair)
            if __debug__:
                log.debug('%s as single %s' %
                          (pid_particle_pair[0], single.domain_id))
            singles.append(single)
        assert len(singles) == self.world.num_particles
        for single in singles:
            self.add_domain_event(single)

	# 3. The simulator is now consistent with 'world'
        self.is_dirty = False

    def stop(self, t):
        """Bring the simulation to a full stop, which synchronizes all
	particles at time t. The state is similar to the state after
	initialization.

        With eGFRD, particle positions are normally updated 
        asynchronously. This method bursts all protective domains so that 
        the position of each particle is known.

        Arguments:
            - t
                the time at which to synchronize the particles. Usually 
                you will want to use the current time of the simulator: 
                EGFRDSimulator.t.

        This method is called stop because it is usually called at the 
        end of a simulation. It is possible to call this method at an 
        earlier time. For example the Logger module does this, because 
        it needs to know the positions of the particles at each log 
        step.

        """
        if __debug__:
            log.info('stop at %s' % (FORMAT_DOUBLE % t))

        if self.t == t:				# We actually already stopped?
            return				# FIXME: is this accurate? Probably use feq from utils.py

        if t >= self.scheduler.top[1].time:	# Can't schedule a stop later than the next event time
            raise RuntimeError('Stop time >= next event time.')

        if t < self.t:
            raise RuntimeError('Stop time < current time.')

        self.t = t

        non_single_list = []

        # first burst all Singles, and put Pairs and Multis in a list.
        for id, event in self.scheduler:
            obj = self.domains[event.data]
            if isinstance(obj, Pair) or isinstance(obj, Multi):
                non_single_list.append(obj)
            elif isinstance(obj, Single):
                if __debug__:
                    log.debug('burst %s, last_time= %s' % 
                              (obj, FORMAT_DOUBLE % obj.last_time))
                self.burst_single(obj)
            else:
                assert False, 'object from scheduler was no Single, Pair or Multi'


        # then burst all Pairs and Multis.
        if __debug__:
            log.debug('burst %s' % non_single_list)
        self.burst_objs(non_single_list)

        self.dt = 0.0

    def step(self):
        """Execute one eGFRD step.

        """
        self.last_reaction = None

        if self.is_dirty:
            self.initialize()
            
        if __debug__:
            if int("0" + os.environ.get("ECELL_CHECK", ""), 10):
                self.check()
        
        self.step_counter += 1

        if __debug__:
            if self.scheduler.size == 0:
                raise RuntimeError('No Events in scheduler.')
	
	# 1. Get the next event from the scheduler
	#
        id, event = self.scheduler.pop()
	domain = self.domains[event.data]
        self.t, self.last_event = event.time, event

        if __debug__:
            domain_counts = self.count_domains()
            log.info('\n\n%d: t=%s dt=%s\t' %
                     (self.step_counter, FORMAT_DOUBLE % self.t,
                      FORMAT_DOUBLE % self.dt) + 
                     'Singles: %d, Pairs: %d, Multis: %d\n' % domain_counts + 
                     'event=#%d reactions=%d rejectedmoves=%d' %
                     (id, self.reaction_events, self.rejected_moves))
       
	# 2. Use the correct method to process (fire) the shell that produced the event
	#
        # Dispatch is a dictionary (hash) of what function to use to fire
        # different classes of shells (see bottom egfrd.py) 
        # event.data holds the object (Single, Pair, Multi) that is associated with the next event.
        # e.g. "if class is Single, then fire_single" ~ MW
        for klass, f in self.dispatch:
            if isinstance(domain, klass):
                f(self, domain)		# fire the correct method for the class (e.g. fire_singel(self, Single))

        if __debug__:
            if self.scheduler.size == 0:
                raise RuntimeError('Zero events left.')

	# 3. Adjust the simulation time
	#
        next_time = self.scheduler.top[1].time
        self.dt = next_time - self.t

        # assert if not too many successive dt=0 steps occur.
        if __debug__:
            if self.dt == 0:
                log.warning('dt=zero step, working in s.t >> dt~0 Python limit.')
                self.zero_steps += 1
                # TODO Changed from 10 to 10000, because only a problem 
                # when reaching certain magnitude.
                if self.zero_steps >= max(self.scheduler.size * 3, self.MAX_NUM_DT0_STEPS): 
                    raise RuntimeError('too many dt=zero steps. '
                                       'Simulator halted?'
                                    'dt= %.300g-%.300g' % (self.scheduler.top[1].time, self.t))
            else:
                self.zero_steps = 0



    def create_single(self, pid_particle_pair):
	# 1. generate identifiers for the domain and shell. The event_id is
	# generated by the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        rrs = self.network_rules.query_reaction_rule(pid_particle_pair[1].sid)
        # Get structure (region or surface) where the particle lives.
        species = self.world.get_species(pid_particle_pair[1].sid)
        structure = self.world.get_structure(species.structure_id)

        # 2. Create and register the single domain.
	# The type of the single that will be created 
        # depends on the structure (region or surface) this particle is 
        # in/on. Either SphericalSingle, PlanarSurfaceSingle, or 
        # CylindricalSurfaceSingle.
        single = create_default_single(domain_id, pid_particle_pair, 
                                       shell_id, rrs, structure)
        assert isinstance(single, NonInteractionSingle)
        single.initialize(self.t)
        self.domains[domain_id] = single

	# 3. update the proper shell container
        self.geometrycontainer.move_shell(single.shell_id_shell_pair)

        if __debug__:
            # Used in __str__.
            single.world = self.world
        return single

    def create_interaction(self, pid_particle_pair, surface, shell_center,
			   shell_radius, shell_half_length, shell_unit_z):

##        assert single.dt == 0.0 and single.get_mobility_radius() == 0.0
#	assert that the particle is not already associate with another domain

	# 1. generate identifiers for the domain and shell. event_id is generated by
	# the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        species_id = pid_particle_pair[1].sid
	species = self.world.get_species(species_id)

        reaction_rules = self.network_rules.query_reaction_rule(species_id)
        #interaction_type = self.query_interaction_rule(species_id, surface)
        # TODO.
        interaction_rules = None
        structure = self.world.get_structure(species.structure_id)

	particle_pos = pid_particle_pair[1].position

    	if isinstance(surface, CylindricalSurface):
            particle_pos = self.world.cyclic_transpose(particle_pos, surface.shape.position)
            projected_point, r0 = surface.projected_point(particle_pos)
            shell_unit_r = normalize(particle_pos - projected_point)

            z0 = numpy.dot (shell_unit_z, (projected_point - shell_center))

            interaction = CylindricalSurfaceInteraction(domain_id, pid_particle_pair,
                                                        reaction_rules, structure,
                                                        shell_id, shell_center, shell_radius,
                                                        shell_half_length, shell_unit_z, z0,
                                                        shell_unit_r, r0, interaction_rules, surface)
	elif isinstance(surface, PlanarSurface):
            particle_pos = self.world.cyclic_transpose(particle_pos, shell_center)
            z0 = numpy.dot (shell_unit_z, (particle_pos - shell_center))

            interaction = PlanarSurfaceInteraction(domain_id, pid_particle_pair,
                                                   reaction_rules, structure,
                                                   shell_id, shell_center, shell_radius,
                                                   shell_half_length, shell_unit_z,
                                                   z0, interaction_rules, surface)


        assert isinstance(interaction, InteractionSingle)
        interaction.initialize(self.t)
        self.domains[domain_id] = interaction

	# 3. update the shell containers 
        self.geometrycontainer.move_shell(interaction.shell_id_shell_pair)

        if __debug__:
            # Used in __str__.
            interaction.world = self.world

        return interaction

    def create_pair(self, single1, single2, com, r0, shell_size):
        assert single1.dt == 0.0
        assert single2.dt == 0.0

	pid_particle_pair1 = single1.pid_particle_pair
	pid_particle_pair2 = single2.pid_particle_pair

	# 1. generate needed identifiers
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # Select 1 reaction type out of all possible reaction types between the two particles.
        rrs = self.network_rules.query_reaction_rule(
                pid_particle_pair1[1].sid,
                pid_particle_pair2[1].sid)
        k_array = numpy.add.accumulate([rr.k for rr in rrs])
        k_max = k_array[-1]
        rnd = myrandom.uniform()
        i = numpy.searchsorted(k_array, rnd * k_max)
        rr = rrs[i]
        # The probability for this reaction to happen is proportional to 
        # the sum of the rates of all the possible reaction types. 
        rr.ktot = k_max

        # Get structure (region or surface) where particle1 lives (assuming particle2
	# also lives here -> TODO).
        species1 = self.world.get_species(pid_particle_pair1[1].sid)
        structure = self.world.get_structure(species1.structure_id)


        # 2. Create pair. The type of the pair that will be created depends
        # on the structure (region or surface) the particles are in/on.  
        # Either SphericalPair, PlanarSurfacePair, or 
        # CylindricalSurfacePair.
        pair = create_default_pair(domain_id, com, single1, single2, shell_id, 
                                   r0, shell_size, rr, structure)

	assert isinstance(pair, Pair)
        pair.initialize(self.t)
        self.domains[domain_id] = pair

	# 3. update the shell containers with the new shell
        self.geometrycontainer.move_shell(pair.shell_id_shell_pair)

        if __debug__:
            # Used in __str__.
            pair.world = self.world
        return pair

    def create_multi(self):
	# 1. generate necessary id's
        domain_id = self.domain_id_generator()

        if __debug__:
            try:
                # Option to make multis run faster for nicer visualization.
                dt_factor = DEFAULT_DT_FACTOR * self.bd_dt_factor
            except AttributeError:
                dt_factor = DEFAULT_DT_FACTOR 
        else:
            dt_factor = DEFAULT_DT_FACTOR

	# 2. Create and register domain object
        multi = Multi(domain_id, self, dt_factor)
        self.domains[domain_id] = multi

	# 3. no shells are yet made, since these are added later
	# -> a multi can have 0 shells
        return multi


    def move_single(self, single, position, radius=None):
        self.move_single_shell(single, position, radius)
        self.move_single_particle(single, position)

    def move_single_shell(self, single, position, radius=None):
        if radius == None:
            # By default, don't change radius.
            radius = single.shell.shape.radius

        # Reuse shell_id and domain_id.
        shell_id = single.shell_id
        domain_id = single.domain_id

        # Replace shell.
        shell = single.create_new_shell(position, radius, domain_id)
        shell_id_shell_pair = (shell_id, shell) 

        single.shell_id_shell_pair = shell_id_shell_pair
        self.geometrycontainer.move_shell(shell_id_shell_pair)

    def move_single_particle(self, single, position):
        new_pid_particle_pair = (single.pid_particle_pair[0],
                          Particle(position,
                                   single.pid_particle_pair[1].radius,
                                   single.pid_particle_pair[1].D,
                                   single.pid_particle_pair[1].sid))
        single.pid_particle_pair = new_pid_particle_pair

        self.world.update_particle(new_pid_particle_pair)

    def remove_domain(self, obj):
	# Removes all the ties to a domain (single, pair, multi) from the system.
	# Note that the particles that it represented still exist
	# in 'world' and that obj also still persits (?!)
        if __debug__:
            log.info("remove: %s" % obj)
	# TODO assert that the domain is not still in the scheduler
        del self.domains[obj.domain_id]
        for shell_id_shell_pair in obj.shell_list:
	    self.geometrycontainer.remove_shell(shell_id_shell_pair)

### TODO These methods can be made methods to the scheduler class
    def add_domain_event(self, domain):
    # This method makes an event for domain 'domain' in the scheduler.
    # The event will have the domain_id pointing to the appropriate domain in 'domains{}'.
    # The domain will have the event_id pointing to the appropriate event in the scheduler.
	event_time = self.t + domain.dt
        event_id = self.scheduler.add(
            DomainEvent(event_time, domain))
        if __debug__:
            log.info('add_event: %s, event=#%d, t=%s' %
                     (domain.domain_id, event_id,
                      FORMAT_DOUBLE % (event_time)))
        domain.event_id = event_id			# FIXME side effect programming -> unclear!!

    def remove_event(self, event):
        if __debug__:
            log.info('remove_event: event=#%d' % event.event_id)
        del self.scheduler[event.event_id]

    def update_domain_event(self, t, domain):
        if __debug__:
            log.info('update_event: %s, event=#%d, t=%s' %
                     (domain.domain_id, domain.event_id, FORMAT_DOUBLE % t))
        self.scheduler.update((domain.event_id, DomainEvent(t, domain)))

    def burst_domain(self, domain):
    # Reduces and domain (Single, Pair or Multi) to Singles with the zero shell, and dt=0
        if __debug__:
            log.info('burst_domain: %s' % domain)

        if isinstance(domain, Single):
            # TODO. Compare with gfrd.
            domain = self.burst_single(domain)
            bursted = [domain, ]
        elif isinstance(domain, Pair):  # Pair
            single1, single2 = self.burst_pair(domain)
            # Don't schedule events in burst/propagate_pair, because 
            # scheduling is different after a single reaction in 
            # fire_pair.
            self.add_domain_event(single1)
            self.add_domain_event(single2)
            self.remove_event(domain)
            bursted = [single1, single2]
        else:  # Multi
#            bursted = self.burst_multi(domain)
            bursted = self.break_up_multi(domain)	# Multi's can't really be 'bursted' since the
							# positions of the particles are known. They
							# are broken up to singles with a dt=0 shell instead.
            self.remove_event(domain)

        if __debug__:
            # After a burst, InteractionSingles should be gone.
            assert all(not isinstance(b, InteractionSingle) for b in bursted)
            log.info('bursted = %s' % ',\n\t  '.join(str(i) for i in bursted))

        return bursted

    def burst_objs(self, domains):
    # bursts all the domains that are in the list 'domains'
        bursted = []
        for domain in domains:
            d = self.burst_domain(domain)
            bursted.extend(d)

        return bursted

    def clear_volume(self, pos, radius, ignore=[]):
        """ Burst domains within a certain volume and give their ids.

        This function actually has a confusing name, as it only bursts 
        domains within a certain radius, and gives their ids. It doesn't
        remove the particles or something like that.

        (Bursting means it propagates the particles within the domain until
        the current time, and then creates a new, minimum-sized domain.)

        Arguments:
            - pos: position of area to be "bursted"
            - radius: radius of area to be "bursted"
            - ignore: domains that should be ignored, none by default.
        """ # ~ MW
	# TODO clear_volume now always assumes a spherical volume to burst.
	# It is inefficient to burst on the 'other' side of the membrane or to
	# burst in a circle, when you want to make a new shell in the membrane.
	# Then we may want to burst in a cylindrical volume, not a sphere
        neighbor_ids = self.geometrycontainer.get_neighbors_within_radius_no_sort(pos, radius,
                                                                               ignore)
	neighbors = [self.domains[domain_id] for domain_id in neighbor_ids]
        return self.burst_objs(neighbors)

    def burst_non_multis(self, domains):
        bursted = []

        for domain in domains:
            if not isinstance(domain, Multi):
		# domain is Single or Pair of some subclass
                single_list = self.burst_domain(domain)
                bursted.extend(single_list)
            else:
		# domain is a Multi
                bursted.append(domain)

        return bursted

    def fire_single_reaction(self, single, interaction_flag=False):
	# This takes care of the identity change when a single particle decays into
	# one or a number of other particles
	# TODO This can also be used when a particle falls off a surface.
        if interaction_flag == True:
            raise NotImplementedError
        reactant_species_radius = single.pid_particle_pair[1].radius
        oldpos = single.pid_particle_pair[1].position
        current_structure = single.structure
        
        rt = single.draw_reaction_rule()

        if len(rt.products) == 0:
            
            self.world.remove_particle(single.pid_particle_pair[0])

            self.last_reaction = (rt, (single.pid_particle_pair[1], None), [])

            
        elif len(rt.products) == 1:
            
            product_species = self.world.get_species(rt.products[0])

            if reactant_species_radius < product_species.radius:
                self.clear_volume(oldpos, product_species.radius)

            if self.world.check_overlap((oldpos, product_species.radius),
                                        single.pid_particle_pair[0]):
                if __debug__:
                    log.info('no space for product particle.')
                raise NoSpace()

            self.world.remove_particle(single.pid_particle_pair[0])
            newparticle = self.world.new_particle(product_species.id, oldpos)
            newsingle = self.create_single(newparticle)
            if __debug__:
                log.info('product = %s' % newsingle)
            self.add_domain_event(newsingle)

            self.last_reaction = (rt, (single.pid_particle_pair[1], None),
                                  [newparticle])


            
        elif len(rt.products) == 2:
            product_species1 = self.world.get_species(rt.products[0])
            product_species2 = self.world.get_species(rt.products[1])
            
            D1 = product_species1.D
            D2 = product_species2.D
            D12 = D1 + D2
            
            particle_radius1 = product_species1.radius
            particle_radius2 = product_species2.radius
            particle_radius12 = particle_radius1 + particle_radius2

            # clean up space.
            rad = max(particle_radius12 * (D1 / D12) + particle_radius1,
                      particle_radius12 * (D2 / D12) + particle_radius2)

            self.clear_volume(oldpos, rad)

            for _ in range(self.dissociation_retry_moves):
                vector = _random_vector(current_structure, particle_radius12 *
                                        MINIMAL_SEPARATION_FACTOR, self.rng)
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?

                while 1:
                    newpos1 = oldpos + vector * (D1 / D12)
                    newpos2 = oldpos - vector * (D2 / D12)
                    newpos1 = self.world.apply_boundary(newpos1)
                    newpos2 = self.world.apply_boundary(newpos2)

                    if(self.world.distance(newpos1, newpos2) >= 
                       particle_radius12):
                        break

                    vector *= 1.0 + 1e-7


                # accept the new positions if there is enough space.
                if(not self.world.check_overlap((newpos1, particle_radius1),
                                                 single.pid_particle_pair[0])
                   and
                   not self.world.check_overlap((newpos2, particle_radius2),
                                                 single.pid_particle_pair[0])):
                    break
            else:
                if __debug__:
                    log.info('no space for product particles.')
                raise NoSpace()

            self.world.remove_particle(single.pid_particle_pair[0])

            particle1 = self.world.new_particle(product_species1.id, newpos1)
            particle2 = self.world.new_particle(product_species2.id, newpos2)
            newsingle1 = self.create_single(particle1)
            newsingle2 = self.create_single(particle2)

            if __debug__:
                log.info('product1 = %s\nproduct2 = %s' % 
                     (newsingle1, newsingle2))

            self.add_domain_event(newsingle1)
            self.add_domain_event(newsingle2)

            self.last_reaction = (rt, (single.pid_particle_pair[1], None),
                                  [particle1, particle2])


        else:
            raise RuntimeError('num products >= 3 not supported.')

        self.reaction_events += 1

    def propagate_single(self, single):
        # The difference between a burst and a propagate is that a burst 
        # always takes place before the actual scheduled event for the 
        # single, while propagate_single can be called for an escape event.

        # Another subtle difference is that burst_single always 
        # reschedules (update_event) the single, while just calling 
        # propagate does not.  So whoever calls propagate_single 
        # directly should reschedule the single afterwards.
        if __debug__:
            log.debug("single.dt=%s, single.last_time=%s, self.t=%s" %
                      (FORMAT_DOUBLE % single.dt,
                       FORMAT_DOUBLE % single.last_time,
                       FORMAT_DOUBLE % self.t))

        newpos = single.draw_new_position(single.dt, single.event_type) 
        newpos = self.world.apply_boundary(newpos)

        if __debug__:
            log.debug("propagate %s: %s => %s" %
                      (single, single.pid_particle_pair[1].position, newpos))

            if self.world.check_overlap((newpos,
                                        single.pid_particle_pair[1].radius),
                                        single.pid_particle_pair[0]):
                raise RuntimeError('propagate_single: check_overlap failed.')

        if(single.event_type == EventType.SINGLE_REACTION or
           single.event_type == EventType.IV_INTERACTION):
            # Single reaction or interaction, and not a burst. Only 
            # update particle, single is removed anyway.
            self.move_single_particle(single, newpos)
            return single
        else:
            if isinstance(single, InteractionSingle):
                # If for an interaction single a single reaction or iv 
                # interaction occurs, we create a new single and get 
                # rid of the old interactionSingle in 
                # fire_single_reaction.
                # For escapes and bursts of interaction singles we do 
                # it here.
                self.move_single_particle(single, newpos)
                newsingle = self.create_single(single.pid_particle_pair)
                self.remove_domain(single)
                if __debug__:
                    log.debug('    *New %s.\n'
                              '        radius = %.3g. dt = %.3g.' %
                              (newsingle, newsingle.shell.shape.radius,
                               newsingle.dt))
                return newsingle
            else:
                single.initialize(self.t)
                self.move_single(single, newpos,
                                 single.pid_particle_pair[1].radius)
                return single

    def fire_single(self, single):
	### 1. Process event produced by the single
	### 1.1 Special cases (shortcuts)
        # In case nothing is scheduled to happen: do nothing; 
        # results also in disappearance from scheduler. (?)
        if single.dt == numpy.inf:
            return 
	# Else continue

        # check timeline (??)
        assert (abs(single.dt + single.last_time - self.t) <= 1e-18 * self.t)

        # If the single had a decay reaction.
        if single.event_type == EventType.SINGLE_REACTION:
            if __debug__:
                log.info('%s' % single.event_type)
                log.info('reactant = %s' % single)

            self.single_steps[single.event_type] += 1


            single = self.propagate_single(single)

            try:
                self.remove_domain(single)
                self.fire_single_reaction(single)
            except NoSpace:
                self.reject_single_reaction(single)

            return

	# Else the event was either an ESCAPE (NonInteractionSingle or InteractionSingle),
	# or IV_EVENT (InteractionSingle only), NO BURST event comes up here!


        if __debug__:
            log.info('%s' % single.event_type)
            log.info('single = %s' % single)

        # Handle immobile case (what event has taken place to get here?!).
        if single.getD() == 0:
	    # The domain has not produced an decay reaction event and is immobile
            # ->no propagation, no other changes just calculate next reaction time
	    # and reschedule the single. So this also takes care of the domain making part
            single.dt, single.event_type = single.determine_next_event() 
            single.last_time = self.t
            self.add_domain_event(single)
            return
        

	### 1.2 General case
	# There is a non-trivial event time for the domain
        if single.dt != 0.0:
            # Propagate this particle to the exit point on the shell.
            single = self.propagate_single(single)


        # An interaction event is similar to a single reaction.
        if single.event_type == EventType.IV_INTERACTION:
            try:
                self.remove_domain(single)
                self.fire_single_reaction(single, interaction_flag=True)
            except NoSpace:
                self.reject_single_reaction(single)

            return

	self.make_new_domain(single)
	return


    def make_new_domain(self, single):
	### Make a new domain out of a NonInteractionSingle that was
	### just fired

	# can only make new domain from NonInteractionSingle
	assert isinstance(single, NonInteractionSingle)
	
	# 2.1
	# We prefer to make NonInteractionSingles for efficiency.
	# But if objects (Shells and surfaces) get close enough (closer than
	# reaction_threshold) we want to try Interaction and/or Pairs

	single_pos = single.pid_particle_pair[1].position
	single_radius = single.pid_particle_pair[1].radius

	# Check if there are shells with the burst radius (reaction_threshold)
	# of the particle (intruders). Note that we approximate the reaction_volume
	# with a sphere (should be cylinder for 2D or 1D particle)
        reaction_threshold = single.pid_particle_pair[1].radius * \
                             self.SINGLE_SHELL_FACTOR

        intruder_ids, closest_id, closest_distance = \
            self.geometrycontainer.get_intruders(single_pos, reaction_threshold,
                                                 ignore=[single.domain_id, ])
	intruders = [self.domains[domain_id] for domain_id in intruder_ids]

        if __debug__:
            log.debug("intruders: %s" %
                      (', '.join(str(i) for i in intruders)))


        # 2.2 Burst the shells with the following conditions
	# -shell is in the burst radius (intruder)
	# -shell is burstable (not Multi)
	# -shell is not newly made
	# -shell is not already a zero shell (just bursted)
	burst = []
        if intruders:
#            burst = self.burst_non_multis(intruders)
	    for domain in intruders:
        	if not isinstance(domain, Multi) and \
                   not self.t == domain.last_time: # and \
#		   not domain.dt == 0.0:
                    # domain is Single or Pair of some subclass
            	    single_list = self.burst_domain(domain)
                    burst.extend(single_list)
                elif isinstance(domain, Multi):
                    # domain is a Multi
                    burst.append(domain)

	    # burst now contains all the Domains resulting from the burst.
	    # These are NonInteractionSingles and Multis


	# 2.3 get the closest object (a Domain or surface) which can be a
	# partner for a new domain. This can be:
	# zero-shell NonInteractionSingle -> Pair or Multi
	# Multi -> Multi
	# Surface -> Interaction

	# Potential partners are for now just the bursted domains
	partners = burst

	dists = []
	if partners:
            # sort partners by distance
	    dists = self.obj_distance_array(single_pos, partners)
            if len(dists) >= 2:
        	n = dists.argsort()
        	partners = numpy.take(partners, n)	# sort the potential partners using the index
        	dists = dists.take(n)			# sort the distances using the index

	    closest_domain = partners[0]
	    domain_distance = dists[0]
	    rest_domains = partners[1:]
	    rest_dists = dists[1:]
	else:
	    closest_domain = None
	    domain_distance = numpy.inf
	    rest_domains = []
	    rest_dists = []
	
	# Potential partners are also surfaces
        closest_surface, surface_distance = \
	get_closest_surface(self.world, single_pos, ignore=[single.structure.id])



	# 2.4 calculate the thresholds (horizons) for trying to form the various domains.
	# Note that we do not differentiate between directions. This means that we
	# look around in a sphere, the horizons are spherical
	surface_horizon = single_radius * self.SINGLE_SHELL_FACTOR
	multi_horizon = single_radius * self.MULTI_SHELL_FACTOR

	# 2.4 check that the object is a potential partner for an
	# -Interaction 	(a surface)
	# -Pair		(a zero dt NonInteractionSingle)
	# -Multi	(a zero dt NonInteractionSingle or Multi)
	if isinstance(closest_domain, NonInteractionSingle) and \
	   closest_domain.is_reset():
		pair_horizon = (single_radius + \
                       closest_domain.pid_particle_pair[1].radius) * self.SINGLE_SHELL_FACTOR
		if domain_distance < pair_horizon:
	            # try making a Pair (can still be Mixed Pair or Normal Pair)
		    domain = self.try_pair (single, closest_domain, rest_domains)
		else:
		    domain = None

	elif isinstance(single, SphericalSingle) and \
	     surface_distance < domain_distance and \
	     surface_distance < surface_horizon:
		# try making an Interaction
		domain = self.try_interaction (single, closest_surface)
	else:
		domain = None


	if not domain:
	    # No Pair or Interaction could be formed
	    # Now choose between NonInteractionSingle and Multi

	    if domain_distance > multi_horizon and \
	       surface_distance > multi_horizon:
	        # just make a normal NonInteractionSingle
	        self.update_single(single)
	        self.add_domain_event(single)
	    else:
		# The closest object was too close to make a NonInteractionSingle
		domain = self.form_multi(single, partners, dists)

	return domain


    def reject_single_reaction(self, single):
        if __debug__:
            log.info('single reaction; placing product failed.')
        self.domains[single.domain_id] = single
        self.geometrycontainer.move_shell(single.shell_id_shell_pair)
        self.rejected_moves += 1
        single.initialize(self.t)
        self.add_domain_event(single)

    def calculate_simplepair_shell_size(self, single1, single2, burst):
        # 1. Determine min shell size for the SimplePair.

	assert single1.structure == single2.structure

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

        sigma = radius1 + radius2

        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D12 = D1 + D2

#        assert (pos1 - single1.shell.shape.position).sum() == 0	# TODO Not sure why this is here
        r0 = self.world.distance(pos1, pos2)
        distance_from_sigma = r0 - sigma
        assert distance_from_sigma >= 0, \
            'distance_from_sigma (pair gap) between %s and %s = %s < 0' % \
            (single1, single2, FORMAT_DOUBLE % distance_from_sigma)

        shell_size1 = r0 * D1 / D12 + radius1
        shell_size2 = r0 * D2 / D12 + radius2
        shell_size_margin1 = radius1 * 2
        shell_size_margin2 = radius2 * 2
        shell_size_with_margin1 = shell_size1 + shell_size_margin1
        shell_size_with_margin2 = shell_size2 + shell_size_margin2
        if shell_size_with_margin1  >= shell_size_with_margin2:
            min_shell_size = shell_size1
            shell_size_margin = shell_size_margin1
        else:
            min_shell_size = shell_size2
            shell_size_margin = shell_size_margin2

        # 2. Check if min shell size for the Pair not larger than max shell size or 
        # sim cell size.
        min_shell_size_with_margin = min_shell_size + shell_size_margin
        max_shell_size = min(self.geometrycontainer.get_max_shell_size(),
                             distance_from_sigma * 100 +
                             sigma + shell_size_margin)

        if min_shell_size_with_margin >= max_shell_size:
            if __debug__:
                log.debug('%s not formed: min_shell_size %s >='
                          'max_shell_size %s' %
                          ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                             single2.pid_particle_pair[0]),
                           FORMAT_DOUBLE % min_shell_size_with_margin,
                           FORMAT_DOUBLE % max_shell_size))
            return None, None, None


        # 3. Check if bursted Singles can still make a minimal shell.
	# TODO This should be changed to a more general algorithm
	# that checks if the closest nearest relevant shell is a
	# NonInteractionSingle or smt else.
	# Now it only processes the shells that were actually bursted before.
	#
	# In case it's a NonInteractionSingle observe a distance of reaction_radius
	# In case it's smt else no extra distance has to be held

        com = self.world.calculate_pair_CoM(pos1, pos2, D1, D2)
        com = self.world.apply_boundary(com)
        closest, closest_shell_distance = None, numpy.inf
        for b in burst:
            if isinstance(b, Single):
                bpos = b.shell.shape.position
                d = self.world.distance(com, bpos) - \
                    b.pid_particle_pair[1].radius * self.SINGLE_SHELL_FACTOR
                if d < closest_shell_distance:
                    closest, closest_shell_distance = b, d

        if closest_shell_distance <= min_shell_size_with_margin:
            if __debug__:
                log.debug('%s not formed: squeezed by burst neighbor %s' %
                          ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                             single2.pid_particle_pair[0]),
                           closest))
            return None, None, None

        assert closest_shell_distance > 0

        # 4. Determine shell size and check if closest object not too 
        # close (squeezing).
        c, d = self.geometrycontainer.get_closest_obj(com, self.domains, ignore=[single1.domain_id,
                                                      single2.domain_id],
                                                      ignores=[single1.structure.id])

        if d < closest_shell_distance:
            closest, closest_shell_distance = c, d

        if __debug__:
            log.debug('Pair closest neighbor: %s %s, '
                      'min_shell_with_margin %s' %
                      (closest, FORMAT_DOUBLE % closest_shell_distance,
                       FORMAT_DOUBLE % min_shell_size_with_margin))

        assert closest_shell_distance > 0

        if isinstance(closest, Single):

            D_closest = closest.pid_particle_pair[1].D
            D_tot = D_closest + D12
            closest_particle_distance = self.world.distance(
                    com, closest.pid_particle_pair[1].position)

            closest_min_radius = closest.pid_particle_pair[1].radius
            closest_min_shell = closest_min_radius * self.SINGLE_SHELL_FACTOR

            # options for shell size:
            # a. ideal shell size
            # b. closest shell is from a bursted single
            # c. closest shell is closer than ideal shell size 
            shell_size = min((D12 / D_tot) *
                            (closest_particle_distance - min_shell_size 
                             - closest_min_radius) + min_shell_size,
                            closest_particle_distance - closest_min_shell,
                            closest_shell_distance)

            shell_size /= SAFETY
            assert shell_size < closest_shell_distance

        else:
            assert isinstance(closest, (Pair, Multi, Surface, None.__class__))

            shell_size = closest_shell_distance / SAFETY

        if shell_size <= min_shell_size_with_margin:
            if __debug__:
                log.debug('%s not formed: squeezed by %s' %
                          ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                             single2.pid_particle_pair[0]),
                           closest))
            return None, None, None


        # 5. Check if singles would not be better.
	# TODO clear this up -> strange rule
        d1 = self.world.distance(com, pos1)
        d2 = self.world.distance(com, pos2)

        if shell_size < max(d1 + single1.pid_particle_pair[1].radius *
                            self.SINGLE_SHELL_FACTOR, \
                            d2 + single2.pid_particle_pair[1].radius * \
                            self.SINGLE_SHELL_FACTOR) * 1.3:
            if __debug__:
                log.debug('%s not formed: singles are better' %
                          'Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                            single2.pid_particle_pair[0]))
            return None, None, None

        shell_size = min(shell_size, max_shell_size)

        assert closest_shell_distance == numpy.inf or \
               shell_size < closest_shell_distance
        assert shell_size >= min_shell_size_with_margin
        assert shell_size <= max_shell_size

	if __debug__:
	    log.info('SimplePair shell made. shell_size=%s, '
                     'closest_shell_distance=%s,\nclosest = %s' %
                     (FORMAT_DOUBLE % shell_size, FORMAT_DOUBLE % closest_shell_distance, closest))

	return com, r0, shell_size


    def calculate_single_shell_size(self, single, closest, 
                                 distance, shell_distance):
	# only if the closest shell is also a 'simple' Single
        assert isinstance(closest, NonInteractionSingle)

        min_radius1 = single.pid_particle_pair[1].radius
        D1 = single.getD()

        if D1 == 0:
            return min_radius1

        D2 = closest.getD()
        min_radius2 = closest.pid_particle_pair[1].radius
        min_radius12 = min_radius1 + min_radius2
        sqrtD1 = math.sqrt(D1)
            
        shell_size = min(sqrtD1 / (sqrtD1 + math.sqrt(D2))
                        * (distance - min_radius12) + min_radius1,
                        shell_distance / SAFETY)
        if shell_size < min_radius1:
            shell_size = min_radius1

        return shell_size

    def update_single(self, single): 
        assert isinstance(single, NonInteractionSingle)	# This only works for 'simple' Singles

	singlepos = single.pid_particle_pair[1].position
	# TODO we need the closest RELEVANT domain, i.e. the domain in the
	# direction in which we can scale the domain
	closest, distance_to_shell = \
            self.geometrycontainer.get_closest_obj(singlepos, self.domains, ignore=[single.domain_id],
               	                                   ignores=[single.structure.id])

        if isinstance(closest, NonInteractionSingle):
            closestpos = closest.shell.shape.position
            distance_to_closest = self.world.distance(singlepos, closestpos)
            new_shell_size = self.calculate_single_shell_size(single, closest, 
                                                      distance_to_closest,
                                                      distance_to_shell)
        else:  # Pair or Multi or Surface
            new_shell_size = distance_to_shell / SAFETY
            new_shell_size = max(new_shell_size,
                                 single.pid_particle_pair[1].radius)

        new_shell_size = min(new_shell_size, self.geometrycontainer.get_max_shell_size())

        # Resize shell, don't change position.
        # Note: this should be done before determine_next_event.
        self.move_single_shell(single, singlepos, new_shell_size)        

        single.dt, single.event_type = single.determine_next_event()
        single.last_time = self.t

    def fire_pair(self, pair):
        assert self.check_obj(pair)

        single1 = pair.single1
        single2 = pair.single2
        particle1 = single1.pid_particle_pair
        particle2 = single2.pid_particle_pair
        pos1 = particle1[1].position
        pos2 = particle2[1].position
        
        if pair.event_type == EventType.IV_EVENT:
            # Draw actual pair event for iv at very last minute.
            r0 = self.world.distance(pos1, pos2)
            pair.event_type = pair.draw_iv_event_type(r0)

        self.pair_steps[pair.event_type] += 1

        if __debug__:
            log.info('FIRE PAIR: %s' % pair.event_type)
            log.info('single1 = %s' % pair.single1)
            log.info('single2 = %s' % pair.single2)


        old_com = pair.com

        # Four cases:
        #  1. Single reaction
        #  2. Pair reaction
        #  3a. IV escape
        #  3b. com escape

        #
        # 1. Single reaction
        #
        if pair.event_type == EventType.SINGLE_REACTION:
            reactingsingle = pair.reactingsingle

            if reactingsingle == single1:
                theothersingle = single2
            else:
                theothersingle = single1

            self.burst_pair(pair)	# TODO this should be a propagate -> cleaner

            self.add_domain_event(theothersingle)

            if __debug__:
                log.info('reactant = %s' % reactingsingle)
            try:
                self.remove_domain(reactingsingle)
                self.fire_single_reaction(reactingsingle)
            except NoSpace:
                self.reject_single_reaction(reactingsingle)

            return
        
        #
        # 2. Pair reaction
        #
        if pair.event_type == EventType.IV_REACTION:
            self.world.remove_particle(single1.pid_particle_pair[0])
            self.world.remove_particle(single2.pid_particle_pair[0])

            if len(pair.rt.products) == 0:
		# no product particles, nothing new to be made
                product = []

	    elif len(pair.rt.products) == 1:
                species3 = self.world.get_species(pair.rt.products[0])
                # calculate new position
                event_type = pair.event_type
                new_com = pair.draw_new_com(pair.dt, event_type)

                if __debug__:
                    shell_size = pair.get_shell_size()
                    assert self.world.distance(old_com, new_com) < \
                           shell_size - species3.radius

                new_com = self.world.apply_boundary(new_com)

		# make product particle
                particle = self.world.new_particle(species3.id, new_com)

		# create a new domain for the new particle
                product = self.create_single(particle)
		# schedule domain
                self.add_domain_event(product)

		# update counters
                self.reaction_events += 1
                self.last_reaction = (pair.rt, (particle1, particle2),
                                      [particle])
            else:
                raise NotImplementedError('num products >= 2 not supported.')

            if __debug__:
                log.info('product = %s' % product)
            self.remove_domain(pair)

            return

        #
        # 3a. Escaping through a_r.
        # 3b. Escaping through a_R.
        #
        elif(pair.event_type == EventType.IV_ESCAPE or
             pair.event_type == EventType.COM_ESCAPE):
            dt = pair.dt
            event_type = pair.event_type
            single1, single2 = self.propagate_pair(pair, dt, event_type)
            self.add_domain_event(single1)
            self.add_domain_event(single2)
        else:
            raise SystemError('Bug: invalid event_type.')

        return

    def fire_multi(self, multi):
        self.multi_steps[3] += 1  # multi_steps[3]: total multi steps
        multi.step()

        if __debug__:
            log.info('FIRE MULTI: %s' % multi.last_event)

        if(multi.last_event == EventType.MULTI_UNIMOLECULAR_REACTION or
           multi.last_event == EventType.MULTI_BIMOLECULAR_REACTION):
            self.reaction_events += 1
            self.last_reaction = multi.last_reaction

        if multi.last_event is not None:		# if an event took place
            self.break_up_multi(multi)
            self.multi_steps[multi.last_event] += 1
        else:
            self.add_domain_event(multi)

    def break_up_multi(self, multi):
    #Dissolves a multi in singles with a zero shell (dt=0)
        singles = []
        for pid_particle_pair in multi.particles:
            single = self.create_single(pid_particle_pair)
            self.add_domain_event(single)
            singles.append(single)

        self.remove_domain(multi)
        return singles

#    def burst_multi(self, multi):
#        #multi.sim.sync()
#        assert isinstance(multi, Multi)
#        singles = self.break_up_multi(multi)
#
#        return singles

    def burst_single(self, single):
        # Sets next event time of single domain in such a way it will end 
        # up at current time if fired, and then outputs newly created 
        # single that is result of firing old single

        # Check correct timeline ~ MW
        assert self.t >= single.last_time
        assert self.t <= single.last_time + single.dt

        # record important single data ~ MW
#        oldpos = single.shell.shape.position
#        old_shell_size = single.get_shell_size()

        particle_radius = single.pid_particle_pair[1].radius

        # Override dt, burst happens before single's scheduled event.
        single.dt = self.t - single.last_time
        # Override event_type. Always call gf.drawR on BURST.
        single.event_type = EventType.BURST
        newsingle = self.propagate_single(single)

        newpos = newsingle.pid_particle_pair[1].position
        # Check if stays within domain ~MW
#        assert self.world.distance(newpos, oldpos) <= \
#               old_shell_size - particle_radius
        # Displacement check is in NonInteractionSingle.draw_new_position.

        if isinstance(single, InteractionSingle):
            # Removing the event has to be done for *bursting* 
            # *Interaction*Singles, not for propagating 
            # InteractionSingles nor for bursting 
            # NonInteractionSingles.
            self.remove_event(single)
            self.add_domain_event(newsingle)
        else:
            assert single == newsingle
            self.update_domain_event(self.t, single)

        assert newsingle.shell.shape.radius == particle_radius

        # Returned single is different from original single in the case 
        # of an InteractionSingle only.
        return newsingle

    def burst_pair(self, pair):
	# TODO this should also schedule the resulting Singles and remove the
	# Pair event -> cleaner
        if __debug__:
            log.debug('burst_pair: %s', pair)

        assert self.t >= pair.last_time
        assert self.t <= pair.last_time + pair.dt

        dt = self.t - pair.last_time 
        # Override event_type. Always call sgf.drawR and pgf.drawR on BURST.
        event_type = EventType.BURST
        single1, single2 = self.propagate_pair(pair, dt, event_type)

        return single1, single2

    def propagate_pair(self, pair, dt, event_type):
        single1 = pair.single1
        single2 = pair.single2

        particle1 = single1.pid_particle_pair
        particle2 = single2.pid_particle_pair

        pos1 = particle1[1].position
        pos2 = particle2[1].position

        if dt > 0.0:
            D1 = particle1[1].D
            D2 = particle2[1].D

            pos2t = self.world.cyclic_transpose(pos2, pos1)
            old_inter_particle = pos2t - pos1
            r0 = self.world.distance(pos1, pos2)
            assert feq(r0, length(old_inter_particle))

            old_com = pair.com

            newpos1, newpos2 = pair.draw_new_positions(dt, r0, 
                                                     old_inter_particle, 
                                                     event_type)

            newpos1 = self.world.apply_boundary(newpos1)
            newpos2 = self.world.apply_boundary(newpos2)
            assert not self.world.check_overlap((newpos1, particle1[1].radius),
                                                particle1[0], particle2[0])
            assert not self.world.check_overlap((newpos2, particle2[1].radius),
                                                particle1[0], particle2[0])
            assert self.check_pair_pos(pair, newpos1, newpos2, old_com,
                                       pair.get_shell_size())
        else:
            newpos1 = particle1[1].position
            newpos2 = particle2[1].position

        if __debug__:
            log.debug("fire_pair: #1 { %s: %s => %s }" %
                      (single1, pos1, newpos1))
            log.debug("fire_pair: #2 { %s: %s => %s }" %
                      (single2, str(pos2), str(newpos2)))

        single1.initialize(self.t)
        single2.initialize(self.t)
        
        self.remove_domain(pair)
        assert single1.domain_id not in self.domains
        assert single2.domain_id not in self.domains
        self.domains[single1.domain_id] = single1
        self.domains[single2.domain_id] = single2
        self.move_single(single1, newpos1, particle1[1].radius)
        self.move_single(single2, newpos2, particle2[1].radius)

	# TODO This is not very clear or clean now. Also doesn't support MixedPair
        if __debug__:
            container = self.geometrycontainer.get_container(single1.shell)
            assert container[single1.shell_id].shape.radius == \
                   single1.shell.shape.radius
            assert container[single2.shell_id].shape.radius == \
                   single2.shell.shape.radius

            if type(single1.shell) is CylindricalShell:
                assert container[single1.shell_id].shape.half_length == \
                       single1.shell.shape.half_length
                assert container[single2.shell_id].shape.half_length == \
                       single2.shell.shape.half_length

        assert single1.shell.shape.radius == particle1[1].radius
        assert single2.shell.shape.radius == particle2[1].radius

        assert self.check_obj(single1)
        assert self.check_obj(single2)

        return single1, single2

    def try_interaction(self, single, surface):
        # Try to form an interaction between the 'single' particle and the 'surface'.

        pid_particle_pair = single.pid_particle_pair
	particle = pid_particle_pair[1]

	# calculate the Projected_point and distance to the surface
        # Cyclic transpose needed when calling surface.projected_point!
        pos_transposed = self.world.cyclic_transpose(particle.position, 
                                                     surface.shape.position) 
        projected_point, projection_distance = \
                surface.projected_point(pos_transposed)

	# projection_distance is relative to the unit_z of the surface -> can be negative
        particle_distance = abs(projection_distance)

        # For interaction with a planar surface, decide orientation. 
	# Note that the orientation_vector is normalized
        orientation_vector = cmp(projection_distance, 0) * surface.shape.unit_z 

        if __debug__:
           log.debug('trying to form Interaction(%s, %s)' % (particle, surface))


	### 1. See if there is enough space for the shell


	# Make sure the maximal cylinder fits in the maximal sphere. Matrix space
	# doesn't allow to check for shells outside the maximal sphere.
	max_cylinder_radius      = self.geometrycontainer.get_max_shell_size()/math.sqrt(2)
	max_cylinder_half_length = max_cylinder_radius

        # Initialize dr, dz_left, dz_right to maximum allowed values.
        # And decide minimal dr, dz_left, dz_right.
        if isinstance(surface, PlanarSurface):
            dr = max_cylinder_radius
            # Leave enough for the particle itself to the left.
            dz_left = particle.radius
            dz_right = max_cylinder_half_length * 2 - dz_left

            min_dr = particle.radius * self.SINGLE_SHELL_FACTOR
            min_dz_left = dz_left
            min_dz_right = particle_distance + particle.radius * self.SINGLE_SHELL_FACTOR

        elif isinstance(surface, CylindricalSurface):
            dr = max_cylinder_radius
            dz_left = max_cylinder_half_length
            dz_right = max_cylinder_half_length

            min_dr = particle_distance + particle.radius * self.SINGLE_SHELL_FACTOR
            min_dz_left = particle.radius * self.SINGLE_SHELL_FACTOR
            min_dz_right = particle.radius * self.SINGLE_SHELL_FACTOR

        # Miedema's algorithm.
        dr, dz_left, dz_right = \
            self.calculate_max_cylinder(single, surface,
                                        projected_point,
                                        particle_distance,
                                        orientation_vector,
                                        dr, dz_left, dz_right)

        dr /= SAFETY
        dz_right /= SAFETY
#	 This will break the conditions below for membrane interaction
#        dz_left /= SAFETY

        # Decide if interaction domain is possible.
        if dr < min_dr or dz_left < min_dz_left or dz_right < min_dz_right:
            if __debug__:
                log.debug('        *Interaction not possible:\n'
                          '            %s +\n'
                          '            %s.\n'
                          '            dr = %.3g. min_dr = %.3g.\n'
                          '            dz_left = %.3g. min_dz_left = %.3g.\n'
                          '            dz_right = %.3g. min_dz_right = %.3g.' %
                          (single, surface, dr, min_dr, dz_left, min_dz_left, 
                           dz_right, min_dz_right))
            return None


	### 2. The shell can be made. Now do what's necessary to make it
        # Compute origin, radius and half_length of cylinder.
        origin = projected_point + ((dz_right - dz_left)/2.0) * orientation_vector
        half_length = (dz_left + dz_right) / 2.0
        radius = dr

        interaction = self.create_interaction(pid_particle_pair, surface,
					      origin, radius, half_length,
					      orientation_vector)

        interaction.dt, interaction.event_type, = \
            interaction.determine_next_event()
        assert interaction.dt >= 0

        self.last_time = self.t

        self.remove_domain(single)
        # the event associated with the single will be removed by the scheduler.

        #assert self.check_obj(interaction)
        self.add_domain_event(interaction)

        if __debug__:
            log.debug('        *create_interaction\n'
                      '            dr = %s. dz_left = %s. dz_right = %s.\n' %
                      (FORMAT_DOUBLE % dr, FORMAT_DOUBLE % dz_left,
                       FORMAT_DOUBLE % dz_right))

        return interaction


    def calculate_max_cylinder(self, single, surface, projected_point, 
                                   particle_distance, orientation_vector,
                                   dr, dz_left, dz_right):
        # Find optimal cylinder around particle and surface, such that 
        # it is not interfering with other shells.
	#
	# To determine the maximal cylindrical shell a starting cylinder is defined
	# * The projected_point is a reference point along the axis of the cylinder,
	#   the position, z, can be chosen freely
	# * orientation_vector is a vector along the z-axis providing orientation
	#   It usually points in the general direction of the particle
	#
        # * dr is the radius of the cylinder.
        # * dz_right is the distance between the projected_point and the face of the
	#   cylinder on the side of the orientation_vector 
        # * dz_left is the distance from the projected_point to the face of the cylinder
	#   on the side oposite of the orientation_vector
	#
	# * particle_distance is the distance from the center of the particle to the
	#   projected_point

	# the search point is the center of the sphere that surrounds the
	# maximal cylinder
        search_point = projected_point + ((dz_right - dz_left)/2.0) * orientation_vector
        all_neighbor_ids = \
            self.geometrycontainer.get_neighbors_within_radius_no_sort(search_point, 
                                                     self.geometrycontainer.get_max_shell_size(),
                                                     ignore=[single.domain_id])
	all_neighbors = [self.domains[domain_id] for domain_id in all_neighbor_ids]

        for domain in all_neighbors:
            if isinstance(domain, Multi):
                for _, shell in domain.shell_list:
                    shell_position = shell.shape.position
                    shell_size = shell.shape.radius
                    dr, dz_left, dz_right = \
                        self.miedema_algorithm(shell_position, shell_size,
                                               projected_point, 
                                               orientation_vector, dr, 
                                               dz_left, dz_right, 
                                               surface, particle_distance)
            else:
                # Get shell size, which normally is the radius of the 
                # shell, but in case of a cylinder around the same 
                # CylindricalSurface it should be the half_length of 
                # the cylinder.
                shell_position = domain.shell.shape.position
                shell_size = domain.get_shell_size()

                # Make bursted singles look bigger,
                # because the size of their shell is only 
                # particle.radius (not yet multiplied by 
                # SINGLE_SHELL_FACTOR)
                # (and no we can not do that immediately after they are 
                # bursted, singles might start overlapping).
                if domain.dt == 0.0 and domain.getD() > 0:
                    # This is one of the bursted singles.
                    # Or a particle that just escaped it's multi.
                    shell_size *= self.SINGLE_SHELL_FACTOR

                dr, dz_left, dz_right = \
                    self.miedema_algorithm(shell_position, shell_size, 
                                           projected_point, 
                                           orientation_vector, dr, 
                                           dz_left, dz_right, surface,
                                           particle_distance)


        return dr, dz_left, dz_right

    def miedema_algorithm(self, shell_position, shell_size, projected_point,
                          orientation_vector, dr, dz_left, dz_right,
                          surface, particle_distance):
        # Find optimal cylinder around particle and surface, such that 
        # it does not interfere with the shell at position 
        # shell_position and with size (radius or half_length) 
        # shell_size.

	shell_position = self.world.cyclic_transpose(shell_position, projected_point)
        shell_vector = shell_position - projected_point

        # Calculate zi and ri for this shell.
        zi = numpy.dot(shell_vector, orientation_vector)
        z_vector = zi * numpy.array(orientation_vector)
        r_vector = numpy.array(shell_vector) - numpy.array(z_vector)
        ri = numpy.linalg.norm(r_vector)

        # Calculate dr_i for this shell.
        dr_i = ri - shell_size

	if isinstance(surface, PlanarSurface):
	    dz_right -= particle_distance
        elif isinstance(surface, CylindricalSurface):
            # Run Miedema's algorithm in the r direction 
            # relative to the particle's position, not to 
            # projected point (r=0).
            dr_i -= particle_distance
	    dr -= particle_distance
	else:
	    raise SystemError('Error in Miedema\'s algorithm: surface is not a Surface')

        # Calculate dz_left_i or dz_right_i (both will usually 
        # be positive values).
        if zi < 0:
            # Calculate dz_left_i for this shell.
            dz_left_i = - zi - shell_size

            # Miedema's algorithm left side.
            if dz_left_i < dz_left and dr_i < dr:
                if dz_left_i > dr_i:
                    dz_left = dz_left_i
                else:
                    dr = dr_i
        else:
            # Calculate dz_right_i for this shell.
            dz_right_i = zi - shell_size

            if isinstance(surface, PlanarSurface):
                # On the particle side (right side), run 
                # Miedema's algorithm in the z direction 
                # relative to the particle's position, not the 
                # projected point (z=0).
                dz_right_i -= particle_distance

            # Miedema's algorithm right side.
            if dz_right_i < dz_right and dr_i < dr:
                if dz_right_i > dr_i:
                    dz_right = dz_right_i
                else:
                    dr = dr_i

	if isinstance(surface, PlanarSurface):
	    dz_right += particle_distance
	elif isinstance(surface, CylindricalSurface):
	    dr += particle_distance

        return dr, dz_left, dz_right


    def try_pair(self, single1, single2, burst):
        if __debug__:
           log.debug('trying to form Pair(%s, %s)' %
                     (single1.pid_particle_pair, single2.pid_particle_pair))

        assert single1.is_reset()
        assert single2.is_reset()

        # Try forming a Pair only if singles are on same structure.
        if single1.structure != single2.structure:
	    # TODO implement shell making method for MixedPair

            if __debug__:
                log.debug('Pair(%s, %s) not formed: not on same structure.' %
                          (single1.pid_particle_pair[0],
                           single2.pid_particle_pair[0]))

            center, r0, shell_size = None, None, None
#	    shell = self.make_mixedpair_shell (single1, single2)
	else:
	    # particles are on the same structure
	    center, r0, shell_size = self.calculate_simplepair_shell_size (single1, single2, burst)


	if shell_size:
            # A shell could be made and makes sense. Create a Pair
            pair = self.create_pair(single1, single2, center, r0, shell_size)

            pair.dt, pair.event_type, pair.reactingsingle = \
            pair.determine_next_event(r0)

            assert pair.dt >= 0

            self.last_time = self.t

            self.remove_domain(single1)
            self.remove_domain(single2)

            # single1 will be removed by the scheduler.
            self.remove_event(single2)

            assert self.check_obj(pair)
            self.add_domain_event(pair)

            if __debug__:
                log.info('%s,\ndt=%s, r0=%s, shell_size=%s, ' %
                         (pair, FORMAT_DOUBLE % pair.dt, FORMAT_DOUBLE % r0, 
                          FORMAT_DOUBLE % shell_size)) 

            return pair
	else:
	    return None
    

    def form_multi(self, single, neighbors, dists):
	# form a Multi with the 'single'
	# The 'neighbors' are neighboring NonInteractionSingles and Multi which
	# can be added to the Multi (this can also be empty)
	# 'dists' are the distances of the 'neighbors'

	# Filter out relevant neighbors if present
	# only consider neighboring bursted domains that are within the Multi horizon
	min_shell = single.pid_particle_pair[1].radius * self.MULTI_SHELL_FACTOR
	dists = numpy.array(dists)	# FIXME Not sure why this is necessary, dists should already be array
        neighbors = [neighbors[i] for i in (dists <= min_shell).nonzero()[0]]
	if neighbors:
            closest = neighbors[0]
	else:
	    closest = None


	# 1. Make new Multi if Necessary
        if isinstance(closest, Multi):
        # if the closest to this Single is a Multi, reuse the Multi.
	    multi = closest
	    neighbors = neighbors[1:]	# the rest of the neighbors are added
        else: 
        # the closest is not a Multi. Can be NonInteractionSingle, surface or
	# nothing. Create new Multi
            multi = self.create_multi()


	# 2. Add the single and neighbors to the Multi
        self.add_to_multi(single, multi)
        self.remove_domain(single)
        for neighbor in neighbors:
            self.add_to_multi_recursive(neighbor, multi)


	# 3. Initialize and (re-)schedule
        multi.initialize(self.t)
        if isinstance(closest, Multi):
            self.update_domain_event(self.t + multi.dt, multi)
	else:
            self.add_domain_event(multi)

        return multi


    def add_to_multi_recursive(self, domain, multi):
        if isinstance(domain, NonInteractionSingle):
            if multi.has_particle(domain.pid_particle_pair[0]):
                # Already in the Multi.
                return
            assert domain.is_reset()
            objpos = domain.shell.shape.position
            
            self.add_to_multi(domain, multi)
            self.remove_domain(domain)
            self.remove_event(domain)

            radius = domain.pid_particle_pair[1].radius * \
                self.MULTI_SHELL_FACTOR
            neighbor_ids = self.geometrycontainer.get_neighbors_within_radius_no_sort(
                    objpos, radius, ignore=[domain.domain_id])
	    neighbors = [self.domains[domain_id] for domain_id in neighbor_ids]

            burst = self.burst_non_multis(neighbors)
            neighbor_dists = self.obj_distance_array(objpos, burst)
            neighbors = [burst[i] for i
                                  in (neighbor_dists <= radius).nonzero()[0]]

            for domain in neighbors:
                self.add_to_multi_recursive(domain, multi)

        elif isinstance(domain, Multi):
            for pp in multi.particles:
                if domain.has_particle(pp[0]):
                    if __debug__:
                        log.debug('%s already added. skipping.' % domain)
                    break
            else:
                self.merge_multis(domain, multi)
        else:
            assert False, 'do not reach here.'  # Pairs are burst

    def add_to_multi(self, single, multi):
        if __debug__:
            log.info('add to multi:\n  %s\n  %s' % (single, multi))

	shell_id = self.shell_id_generator()
	shell = multi.create_new_shell(single.pid_particle_pair[1].position,
                single.pid_particle_pair[1].radius * self.MULTI_SHELL_FACTOR)
	shell_id_shell_pair = (shell_id, shell)
        self.geometrycontainer.move_shell(shell_id_shell_pair)
        multi.add_shell(shell_id_shell_pair)
        multi.add_particle(single.pid_particle_pair)

    def merge_multis(self, multi1, multi2):
        # merge multi1 into multi2. multi1 will be removed.
        if __debug__:
            log.info('merging %s to %s' % (multi1.domain_id, multi2.domain_id))
            log.info('  %s' % multi1)
            log.info('  %s' % multi2)

            try:
                particle_of_multi1 = iter(multi1.particle_container).next()
                assert particle_of_multi1[0] not in \
                        multi2.particle_container
            except:
                pass

        for sid_shell_pair in multi1.shell_list:
            sid_shell_pair[1].did = multi2.domain_id
            self.geometrycontainer.move_shell(sid_shell_pair)
            multi2.add_shell(sid_shell_pair)

        for pid_particle_pair in multi1.particles:
            multi2.add_particle(pid_particle_pair)

        del self.domains[multi1.domain_id]
        self.remove_event(multi1)

    def domain_distance(self, pos, domain):
	# calculates the shortest distance from 'pos' to A shell (a Multi
	# can have multiple) of 'domain'.
	# Note: it returns the distance between pos and the CENTER of the shell
        return min(self.world.distance(shell.shape, pos)
                   for i, (_, shell) in enumerate(domain.shell_list))

    def obj_distance_array(self, pos, domains):
        dists = numpy.array([self.domain_distance(pos, domain) for domain in domains])
        return dists
            

    #
    # statistics reporter
    #

    def print_report(self, out=None):
        """Print various statistics about the simulation.
        
        Arguments:
            - None

        """
        report = '''
t = %g
steps = %d 
\tSingle:\t%d\t(escape: %d, reaction: %d)
\tInteraction: %d\t(escape: %d, interaction: %d)
\tPair:\t%d\t(escape r: %d, R: %d, reaction pair: %d, single: %d)
\tMulti:\t%d\t(escape: %d, reaction pair: %d, single: %d)
total reactions = %d
rejected moves = %d
''' \
            % (self.t, self.step_counter,
               numpy.array(self.single_steps.values()).sum(),
               self.single_steps[EventType.SINGLE_ESCAPE],
               self.single_steps[EventType.SINGLE_REACTION],
               numpy.array(self.interaction_steps.values()).sum(),
               self.interaction_steps[EventType.IV_ESCAPE],
               self.interaction_steps[EventType.IV_INTERACTION],
               numpy.array(self.pair_steps.values()).sum(),
               self.pair_steps[EventType.IV_ESCAPE],
               self.pair_steps[EventType.COM_ESCAPE],
               self.pair_steps[EventType.IV_REACTION],
               self.pair_steps[EventType.SINGLE_REACTION],
               self.multi_steps[3], # total multi steps
               self.multi_steps[EventType.MULTI_ESCAPE],
               self.multi_steps[EventType.MULTI_BIMOLECULAR_REACTION],
               self.multi_steps[EventType.MULTI_UNIMOLECULAR_REACTION],
               self.reaction_events,
               self.rejected_moves
               )

        print >> out, report

    #
    # consistency checkers
    #

    def check_obj(self, obj):
        obj.check()

        if isinstance(obj, Multi):
            # Ignore all surfaces, multi shells can overlap with 
            # surfaces.
            ignores = [s.id for s in self.world.structures]
	elif isinstance(obj, InteractionSingle):
	    # Ignore surface of the particle and interaction surface
	    ignores = [obj.structure.id, obj.surface.id]
        else:
            ignores = [obj.structure.id]

        for shell_id, shell in obj.shell_list:
            closest, distance = self.geometrycontainer.get_closest_obj(shell.shape.position,
								       self.domains,
                                                                       ignore=[obj.domain_id],
                                                                       ignores=ignores)
	#TODO
            if(type(shell.shape) is Cylinder and
               closest and type(closest.shell.shape) is Sphere):
                # Note: this case is special.
                # Note: only checking if cylinder doesn't overlap with 
                # closest sphere, like we do here, is not really 
                # sufficient (but checking all spheres is too much 
                # work).
                shell_size = shell.shape.half_length
                # Reverse overlap calculation: from closest sphere to 
                # cylinder is easier than the other way around, because 
                # the distance calculation from a point to a cylinder 
                # is already implemented.
                sphere = closest.shell.shape
                diff = self.world.distance(shell.shape, sphere.position) - \
                       sphere.radius
            else:
                if(type(obj) is CylindricalSurfaceSingle or
                   type(obj) is CylindricalSurfacePair or
                   type(obj) is CylindricalSurfaceInteraction):
                    # On CylindricalSurface, use half_lenghts.
                    # (assume all nearby other cylinders are on the 
                    # same surface)
                    shell_size = shell.shape.half_length
                else:
                    # Normally compare radii.
                    shell_size = shell.shape.radius
                diff = distance - shell_size

            assert shell_size <= self.geometrycontainer.get_user_max_shell_size(), \
                '%s shell size larger than user-set max shell size' % \
                str(shell_id)

            assert shell_size <= self.geometrycontainer.get_max_shell_size(), \
                '%s shell size larger than simulator cell size / 2' % \
                str(shell_id)

            assert diff >= 0.0, \
                '%s overlaps with %s. (shell: %s, dist: %s, diff: %s.' % \
                (str(obj), str(closest), FORMAT_DOUBLE % shell_size,
                 FORMAT_DOUBLE % distance,
                 FORMAT_DOUBLE % diff)

        return True

    def check_obj_for_all(self):
        for id, event in self.scheduler:
	    domain = self.domains[event.data]
            self.check_obj(domain)

    def check_event_stoichiometry(self):
    # checks if the number of particles in the world is equal to the number
    # of particles represented in the EGFRD simulator
	world_population = self.world.num_particles
        domain_population = 0
        for id, event in self.scheduler:
	    domain = self.domains[event.data]
            domain_population += domain.multiplicity

        if world_population != domain_population:
            raise RuntimeError('population %d != domain_population %d' %
                               (world_population, domain_population))

    def check_shell_matrix(self):
	did_map, shell_map = self.geometrycontainer.get_dids_shells()

        shell_population = 0
        for id, event in self.scheduler:
	    domain = self.domains[event.data]
            shell_population += domain.num_shells
            shell_ids = did_map[domain.domain_id]
            if len(shell_ids) != domain.num_shells:
                diff = set(sid for (sid, _)
                               in domain.shell_list).difference(shell_ids)
                for sid in diff:
                    print shell_map.get(sid, None)

                raise RuntimeError('number of shells are inconsistent '
                                   '(%d != %d; %s) - %s' %
                                   (len(shell_ids), domain.num_shells, 
                                    domain.domain_id, diff))

        matrix_population = self.geometrycontainer.get_total_num_shells()
        if shell_population != matrix_population:
            raise RuntimeError('num shells (%d) != matrix population (%d)' %
                               (shell_population, matrix_population))

    def check_domains(self):
	# make set of event_id that are stored in all the domains
        event_ids = set(domain.event_id
                        for domain in self.domains.itervalues())
	# check that all the event_id in the scheduler are also stored in a domain
        for id, event in self.scheduler:
            if id not in event_ids:
                raise RuntimeError('Event %s in EventScheduler has no domain in self.domains' %
                                   event)
	    else:
                event_ids.remove(id)

        # self.domains always include a None  --> this can change in future
        if event_ids:
            raise RuntimeError('following domains in self.domains are not in '
                               'Event Scheduler: %s' % str(tuple(event_ids)))

    def check_pair_pos(self, pair, pos1, pos2, com, radius):
        particle1 = pair.single1.pid_particle_pair[1]
        particle2 = pair.single2.pid_particle_pair[1]

        old_com = com
        
        # debug: check if the new positions are valid:
        new_distance = self.world.distance(pos1, pos2)
        particle_radius12 = particle1.radius + particle2.radius

        # check 1: particles don't overlap.
        if new_distance <= particle_radius12:
            if __debug__:
                log.info('rejected move: radii %s, particle distance %s',
                         (FORMAT_DOUBLE % particle1.radius + particle2.radius,
                          FORMAT_DOUBLE % new_distance))
            if __debug__:
                log.debug('DEBUG: pair.dt %s, pos1 %s, pos2 %s' %
                          (FORMAT_DOUBLE % pair.dt, FORMAT_DOUBLE % pos1,
                           FORMAT_DOUBLE % pos2))
            raise RuntimeError('New particles overlap')

        # check 2: particles within mobility radius.
        d1 = self.world.distance(old_com, pos1) + particle1.radius
        d2 = self.world.distance(old_com, pos2) + particle2.radius
        if d1 > radius or d2 > radius:
            raise RuntimeError('New particle(s) out of protective sphere. ' 
                               'radius = %s, d1 = %s, d2 = %s ' %
                               (FORMAT_DOUBLE % radius, FORMAT_DOUBLE % d1,
                                FORMAT_DOUBLE % d2))

        return True




    def check(self):
        ParticleSimulatorBase.check(self)

        assert self.scheduler.check()

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.check_shell_matrix()
        self.check_domains()
        self.check_event_stoichiometry()
        
        self.check_obj_for_all()

    #
    # methods for debugging.
    #

    def dump_scheduler(self):
        """Dump scheduler information.

        """
        for id, event in self.scheduler:
            print id, event

    def dump(self):
        """Dump scheduler and event information.

        """
        for id, event in self.scheduler:
            print id, event, self.domains[event.data]

    def count_domains(self):
        # Returns a tuple (# Singles, # Pairs, # Multis).

        num_singles = 0
        num_pairs = 0
        num_multis = 0
        for d in self.domains.itervalues():
            if isinstance(d, Single):
                num_singles += 1
            elif isinstance(d, Pair):
                num_pairs += 1
            elif isinstance(d, Multi):
                num_multis += 1
            else:
                raise RuntimeError('DO NOT GET HERE')

        return (num_singles, num_pairs, num_multis)

    dispatch = [
        (Single, fire_single),
        (Pair, fire_pair),
        (Multi, fire_multi)
        ]


