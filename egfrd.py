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
    Plane,
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

log = logging.getLogger('ecell')

if __debug__:
    PRECISION = 7
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


def create_default_pair(domain_id, single1, single2, shell_id, com,
                        shell_size, r0, rrs, structure):
    if isinstance(structure, CuboidalRegion):
        return SphericalPair(domain_id, com, single1, single2,
                             shell_id, r0, shell_size, rrs, structure)
    elif isinstance(structure, CylindricalSurface):
        return CylindricalSurfacePair(domain_id, com, single1, single2,
                                      shell_id, r0, shell_size, rrs, structure)
    elif isinstance(structure, PlanarSurface):
        return PlanarSurfacePair(domain_id, com, single1, single2,
                                 shell_id, r0, shell_size, rrs, structure)


class DomainEvent(Event):
    __slot__ = ['data']
    def __init__(self, time, domain):
        Event.__init__(self, time)
        self.data = domain.domain_id            # Store the domain_id key refering to the domain
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
        self.MULTI_SHELL_FACTOR = math.sqrt(2)  # This stems from the fact that there vacant space in the cylinder
                                                # NonInteractionSingles to a Multi and also defines the Multi
                                                # shell size.
        self.SINGLE_SHELL_FACTOR = 2.0          # This is the threshold for when the algorithm switches from
                                                # NonInteractionSingles to a Pair or Interaction. It also defines
                                                # the radius in which the NonInteractionSingle will burst.
                                                # NOTE THAT ABOVE CONSTANTS ARE ALSO IN domain.py
        self.MAX_NUM_DT0_STEPS = 1000

        self.MAX_TIME_STEP = 10

        self.DEFAULT_DT_FACTOR = 1e-5           # Diffusion time prefactor in oldBD algortithm to determine time step.
        
        self.DEFAULT_STEP_SIZE_FACTOR = 0.05    # The maximum step size in the newBD algorithm is determined as DSSF * sigma_min.
                                                # Make sure that DEFAULT_STEP_SIZE_FACTOR < MULTI_SHELL_FACTOR, or else the 
                                                # reaction volume sticks out of the multi.

        # used datastructrures
        self.scheduler = EventScheduler()       # contains the events. Note that every domains has exactly one event

        self.domains = {}                       # a dictionary containing references to the domains that are defined
                                                # in the simulation system. The id of the domain (domain_id) is the key.
                                                # The domains can be a single, pair or multi of any type.

        # other stuff
        self.is_dirty = True                    # The simulator is dirty if the state if the simulator is not
                                                # consistent with the content of the world that it represents
                                                # (or if we don't know for sure)

        self.reset()                            # The Simulator is only initialized at the first step, allowing
                                                # modifications to the world to be made before simulation starts

    def get_next_time(self):
        """ 
        Returns the time it will be when the next egfrd timestep
        is completed.        
        """ #~MW
        if self.scheduler.size == 0:
            return self.t                       # self.t is the current time of the simulator
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
                             EventType.SINGLE_REACTION:0,
                             EventType.BURST:0}
        self.interaction_steps = {EventType.IV_INTERACTION:0,
                                  EventType.IV_ESCAPE:0,
                                  EventType.BURST:0}
        self.pair_steps = {EventType.SINGLE_REACTION:0,
                           EventType.IV_REACTION:0,
                           EventType.IV_ESCAPE:0,
                           EventType.COM_ESCAPE:0,
                           EventType.BURST:0}
        self.multi_steps = {EventType.MULTI_ESCAPE:0,
                            EventType.MULTI_UNIMOLECULAR_REACTION:0,
                            EventType.MULTI_BIMOLECULAR_REACTION:0, 3:0}
        self.zero_steps = 0
        self.rejected_moves = 0
        self.reaction_events = 0
        self.last_event = None
        self.last_reaction = None

        self.is_dirty = True            # simulator needs to be re-initialized

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

        if self.t == t:                         # We actually already stopped?
            return                              # FIXME: is this accurate? Probably use feq from utils.py

        if t >= self.scheduler.top[1].time:     # Can't schedule a stop later than the next event time
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
        if event.time == numpy.inf:
            self.t, self.last_event = self.MAX_TIME_STEP, event
        else:
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
        # e.g. "if class is Single, then process_single_event" ~ MW
        for klass, f in self.dispatch:
            if isinstance(domain, klass):
                f(self, domain)         # fire the correct method for the class (e.g. fire_singel(self, Single))

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
                log.debug('dt=zero step, working in s.t >> dt~0 Python limit.')
                self.zero_steps += 1
                # TODO Changed from 10 to 10000, because only a problem 
                # when reaching certain magnitude.
                if self.zero_steps >= max(self.scheduler.size * 3, self.MAX_NUM_DT0_STEPS): 
                    raise RuntimeError('too many dt=zero steps. '
                                       'Simulator halted?'
                                    'dt= %.10g-%.10g' % (self.scheduler.top[1].time, self.t))
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
        single.initialize(self.t)               # set the initial event and event time
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
#       assert that the particle is not already associate with another domain

        # 1. generate identifiers for the domain and shell. event_id is generated by
        # the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        species_id = pid_particle_pair[1].sid
        species = self.world.get_species(species_id)
        # the structure on which the particle lives
        structure = self.world.get_structure(species.structure_id)
        # the surface_id if the interaction surface
        surface_id = surface.sid

        # get unimolecular reaction rules
        reaction_rules = self.network_rules.query_reaction_rule(species_id)
        # get reaction rules for interaction
        interaction_rules = self.network_rules.query_reaction_rule(species_id, surface_id)


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

    def create_pair(self, single1, single2, shell_center, shell_radius, r0, shell_half_length=0,
                    shell_orientation_vector=None):

        assert single1.dt == 0.0
        assert single2.dt == 0.0

        pid_particle_pair1 = single1.pid_particle_pair
        pid_particle_pair2 = single2.pid_particle_pair

        # 1. generate needed identifiers
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # Select 1 reaction type out of all possible reaction types between the two particles.
        rrs = self.network_rules.query_reaction_rule(pid_particle_pair1[1].sid,
                                                     pid_particle_pair2[1].sid)

        # Get structure (region or surface) where particles live.
        species1 = self.world.get_species(pid_particle_pair1[1].sid)
        species2 = self.world.get_species(pid_particle_pair2[1].sid)
        structure1 = self.world.get_structure(species1.structure_id)
        structure2 = self.world.get_structure(species2.structure_id)

        # 2. Create pair. The type of the pair that will be created depends
        # on the structure (region or surface) the particles are in/on.  
        if structure1 == structure2:
            # Either SphericalPair, PlanarSurfacePair, or 
            # CylindricalSurfacePair.
            pair = create_default_pair(domain_id, single1, single2, shell_id, shell_center,
                                       shell_radius, r0, rrs, structure1)
        else:
            # MixedPair (3D/2D)
            assert shell_orientation_vector != None
            assert shell_half_length != 0
            pair = MixedPair(domain_id, single1, single2, shell_id, shell_center, shell_radius,
                             shell_half_length, shell_orientation_vector, r0, rrs)

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

        # 2. Create and register domain object
        multi = Multi(domain_id, self, self.DEFAULT_STEP_SIZE_FACTOR)
        self.domains[domain_id] = multi

        # 3. no shells are yet made, since these are added later
        # -> a multi can have 0 shells
        return multi


#    def move_single(self, single, position, radius=None):
#        single.pid_particle_pair = self.move_particle(single.pid_particle_pair, position)
#        self.update_single_shell(single, position, radius)

    def update_single_shell(self, single, position, radius=None):
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

    def move_particle(self, pid_particle_pair, position):
        # moves a particle in world based on an existing particle

        new_pid_particle_pair = (pid_particle_pair[0],
                                 Particle(position,
                                          pid_particle_pair[1].radius,
                                          pid_particle_pair[1].D,
                                          pid_particle_pair[1].sid))

        self.world.update_particle(new_pid_particle_pair)

        return new_pid_particle_pair

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
        domain.event_id = event_id                      # FIXME side effect programming -> unclear!!

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
    # return a list of Singles that was the result of the bursting
        if __debug__:
            log.info('burst_domain: %s' % domain)

        if isinstance(domain, Single):
            # TODO. Compare with gfrd.
            bursted = self.burst_single(domain)
        elif isinstance(domain, Pair):  # Pair
            bursted = self.burst_pair(domain)
        else:  # Multi
#            bursted = self.burst_multi(domain)
            bursted = self.break_up_multi(domain)       # Multi's can't really be 'bursted' since the
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

    def burst_volume(self, pos, radius, ignore=[]):
        """ Burst domains within a certain volume and give their ids.

        (Bursting means it propagates the particles within the domain until
        the current time, and then creates a new, minimum-sized domain.)

        Arguments:
            - pos: position of area to be "bursted"
            - radius: radius of area to be "bursted"
            - ignore: domains that should be ignored, none by default.
        """ # ~ MW
        # TODO burst_volume now always assumes a spherical volume to burst.
        # It is inefficient to burst on the 'other' side of the membrane or to
        # burst in a circle, when you want to make a new shell in the membrane.
        # Then we may want to burst in a cylindrical volume, not a sphere
        neighbor_ids = self.geometrycontainer.get_neighbors_within_radius_no_sort(pos, radius,
                                                                               ignore)
        neighbors = [self.domains[domain_id] for domain_id in neighbor_ids]
        return self.burst_objs(neighbors)

    def burst_non_multis(self, domain_distances, burstradius):
        # Bursts the domains in 'domains' if within burstradius
        # Returns the updated list of domains with their distances

        domain_distances_updated = []

        for domain, distance in domain_distances:
            if distance <= burstradius and \
               not isinstance(domain, Multi) and \
               not (isinstance(domain, NonInteractionSingle) and domain.is_reset()) and \
               self.t != domain.last_time:
                # Burst domain only if (AND):
                # -within radius
                # -domain is not a multi
                # -domain is not zero dt NonInteractionSingle
                # -domain has time passed
                single_list = self.burst_domain(domain)
                single_distances = [(single, 0) for single in single_list]  # TODO? calculate real distance, now just 0
                domain_distances_updated.extend(single_distances)
            else:
                # Don't burst domain if (OR):
                # -domain is farther than burst radius
                # -domain is a multi
                # -domain is already a burst NonInteractionSingle (zero dt NonInteractionSingle)
                # -domain is domain in which no time has passed
                domain_distances_updated.append((domain, distance))

        return domain_distances_updated

    def fire_single_reaction(self, single, reactant_pos):
        # This takes care of the identity change when a single particle decays into
        # one or a number of other particles
        # It performs the reactions:
        # A(any structure) -> 0
        # A(any structure) -> B(same structure or 3D)
        # A(any structure) -> B(same structure) + C(same structure or 3D)

        # 0. get reactant info
        reactant           = single.pid_particle_pair
        reactant_radius    = reactant[1].radius
        reactant_structure = single.structure
        rr = single.reactionrule

        # 1. remove the particle

        if len(rr.products) == 0:
            
            # 1. No products, no info
            # 2. No space required
            # 3. No space required
            # 4. process the changes (remove particle, make new ones)
            self.world.remove_particle(reactant[0])
            products = []

            # 5. No new single to be made
            # 6. Log the change

            
        elif len(rr.products) == 1:
            
            # 1. get product info
            product_species = self.world.get_species(rr.products[0])
            product_radius  = product_species.radius
            product_structure = self.world.get_structure(product_species.structure_id) 

            # 1.5 produce a number of new possible positions for the product particle
            # TODO make these generators for efficiency
            if product_structure != reactant_structure:
                assert isinstance(product_structure, CuboidalRegion)

                if isinstance(reactant_structure, PlanarSurface):
                    a = myrandom.choice(-1, 1)
                    directions = [-a,a]
                    vector_length = (product_radius + 0.0) * MINIMAL_SEPARATION_FACTOR  # the thickness of the membrane is 0.0
                    product_pos_list = [reactant_pos + vector_length * reactant_structure.shape.unit_z * direction \
                                        for direction in directions]

                elif isinstance(reactant_structure, CylindricalSurface):
                    vector_length = (product_radius + reactant_structure.shape.radius) * MINIMAL_SEPARATION_FACTOR
                    product_pos_list = []
                    for _ in range(self.dissociation_retry_moves):
                        unit_vector3D = random_unit_vector()
                        unit_vector2D = normalize(unit_vector3D - numpy.dot(unit_vector3D, reactant_structure.shape.unit_z))
                        vector = reactant_pos + vector_length * unit_vector2D
                        product_pos_list.append(vector)
                else:
                    # cannot decay from 3D to other structure
                    raise RuntimeError('fire_single_reaction: Can not decay from 3D to other structure')
            else:
                product_pos_list = [reactant_pos]               # no change of position is required if structure doesn't change


            for product_pos in product_pos_list:
                product_pos = self.world.apply_boundary(product_pos)

                # 2. make space for the products.
#                if product_pos != reactant_pos or product_radius > reactant_radius:
                self.burst_volume(product_pos, product_radius)

                # 3. check that there is space for the products 
                if (not self.world.check_overlap((product_pos, product_radius), reactant[0])):

                    # 4. process the changes (remove particle, make new ones)
                    self.world.remove_particle(reactant[0])
                    newparticle = self.world.new_particle(product_species.id, product_pos)
                    products = [newparticle]

                    # 5. update counters
                    self.reaction_events += 1
                    self.last_reaction = (rr, (reactant[1], None), products)

                    # 6. Log the change
                    if __debug__:
                        log.info('product (%s) = %s' % (len(products), products))

                    # exit the loop, we have found a new position
                    break

            else:
                    # 4. Process (the lack of) change
                    moved_reactant = self.move_particle(reactant, reactant_pos)
                    products = [moved_reactant]     # no change

                    # 5. update counters
                    self.rejected_moves += 1

                    # 6. Log the event
                    if __debug__:
                        log.info('single reaction; placing product failed.')

            
        elif len(rr.products) == 2:
            # 1. get product info
            product1_species = self.world.get_species(rr.products[0])
            product2_species = self.world.get_species(rr.products[1])
            product1_structure = self.world.get_structure(product1_species.structure_id) 
            product2_structure = self.world.get_structure(product2_species.structure_id) 
            
            D1 = product1_species.D
            D2 = product2_species.D
            
            product1_radius = product1_species.radius
            product2_radius = product2_species.radius
            particle_radius12 = product1_radius + product2_radius

            # calculate the displace of the particles according to the ratio D1:D2
            # this way, species with small D only have small displacement from the reaction point.
            if D1 == D2:
            # special case (this should also catch D1==D2==0.0)
                product1_displacement = particle_radius12 * 0.5
                product2_displacement = particle_radius12 * 0.5
            else:
                D12 = D1 + D2
                product1_displacement = particle_radius12 * (D1 / D12)
                product2_displacement = particle_radius12 * (D2 / D12)

            product1_displacement *= MINIMAL_SEPARATION_FACTOR
            product2_displacement *= MINIMAL_SEPARATION_FACTOR


            if product1_structure != product2_structure:
                # make sure that the reactant was on a 2D or 1D surface
                assert not isinstance(reactant_structure, CuboidalRegion)
                # make sure that only either one of the target structures is the 3D ('^' is exclusive-OR)
                assert (isinstance(product1_structure, CuboidalRegion) ^ isinstance(product2_structure, CuboidalRegion))

                # figure out which product stays in the surface and which one goes to 3D
                if isinstance(product2_structure, CuboidalRegion):
                    # product2 goes to 3D and is now particleB (product1 is particleA)
                    productA_radius = product1_radius
                    productB_radius = product2_radius
                    DA = D1
                    DB = D2
                    default = True      # we like to think of this as the default
                else:
                    # product1 goes to 3D and is now particleB (product2 is particleA)
                    productA_radius = product2_radius
                    productB_radius = product1_radius
                    DA = D2
                    DB = D1
                    default = False
                # Note that A is a particle in the surface and B is in the 3D

                if isinstance(reactant_structure, PlanarSurface):

                    # draw a number of new positions for the two product particles
                    # TODO make this into a generator

                    product_pos_list = []
                    for _ in range(self.dissociation_retry_moves):
                        # draw the random angle for the 3D particle relative to the particle left in the membrane
                        # not all angles are available because particle cannot intersect with the membrane

                        # do the backtransform with a random iv with length such that the particles are at contact
                        # Note we make the iv slightly longer because otherwise the anisotropic transform will produce illegal
                        # positions
                        iv = random_vector(particle_radius12 * MixedPair.calc_z_scaling_factor(DA, DB))
                        iv *= MINIMAL_SEPARATION_FACTOR

                        unit_z = reactant_structure.shape.unit_z * myrandom.choice(-1, 1)
                        newposA, newposB = MixedPair.do_back_transform(reactant_pos, iv, DA, DB,
                                                                       productA_radius, productB_radius,
                                                                       reactant_structure, unit_z)

                        newposA = self.world.apply_boundary(newposA)
                        newposB = self.world.apply_boundary(newposB)

                        assert (self.world.distance(newposA, newposB) >= particle_radius12)
                        assert (self.world.distance(reactant_structure.shape, newposB) >= productB_radius)

                        if default:
                            newpos1, newpos2 = newposA, newposB
                        else:
                            newpos1, newpos2 = newposB, newposA

                        product_pos_list.append((newpos1, newpos2))

                elif isinstance(reactant_structure, CylindricalSurface):

                    product_pos_list = []
                    for _ in range(self.dissociation_retry_moves):

                        iv = random_vector(particle_radius12 * MixedPair3D1D.calc_r_scaling_factor(DA, DB))
                        iv *= MINIMAL_SEPARATION_FACTOR

                        newposA, newposB = MixedPair3D1D.do_back_transform(reactant_pos, iv, DA, DB,
                                                                           productA_radius, productB_radius, reactant_structure)

                        newposA = self.world.apply_boundary(newposA)
                        newposB = self.world.apply_boundary(newposB)

                        assert (self.world.distance(newposA, newposB) >= particle_radius12)
                        assert (self.world.distance(reactant_structure.shape, newposB) >= productB_radius)

                        if default:
                            newpos1, newpos2 = newposA, newposB
                        else:
                            newpos1, newpos2 = newposB, newposA

                        product_pos_list.append((newpos1, newpos2))

                else:
                    # cannot decay from 3D to other structure
                    raise RuntimeError('fire_single_reaction: Can not decay from 3D to other structure')

            else:
                # The two particles stay on the same structure.

                # generate new positions in the structure
                # TODO Make this into a generator
                product_pos_list = []
                for _ in range(self.dissociation_retry_moves):
                    # FIXME for particles on the cylinder there are only two possibilities
                    # calculate a random vector in the structure with unit length
                    orientation_vector = _random_vector(reactant_structure, 1.0, self.rng)
            
                    newpos1 = reactant_pos + product1_displacement * orientation_vector
                    newpos2 = reactant_pos - product2_displacement * orientation_vector
                    newpos1 = self.world.apply_boundary(newpos1)
                    newpos2 = self.world.apply_boundary(newpos2)

                    assert (self.world.distance(newpos1, newpos2) >= particle_radius12)
                    product_pos_list.append((newpos1, newpos2))

            # 2. make space for the products. 
            # TODO
            #    calculate the sphere around the two product particles
            rad = max(product1_displacement + product1_radius,
                      product2_displacement + product2_radius)
            self.burst_volume(reactant_pos, rad)

            # 3. check that there is space for the products (try different positions if possible)
            # accept the new positions if there is enough space.
            for newpos1, newpos2 in product_pos_list:

                if (not self.world.check_overlap((newpos1, product1_radius), reactant[0])
                   and
                   not self.world.check_overlap((newpos2, product2_radius), reactant[0])):
                    
                    # 4. process the changes (remove particle, make new ones)
                    self.world.remove_particle(reactant[0])
                    product1_pair = self.world.new_particle(product1_species.id, newpos1)
                    product2_pair = self.world.new_particle(product2_species.id, newpos2)
                    products = [product1_pair, product2_pair]

                    # 5. update counters
                    self.reaction_events += 1
                    self.last_reaction = (rr, (reactant[1], None), products)

                    # 6. Log the change
                    if __debug__:
                        log.info('products (%s) = %s' % (len(products), products))

                    # exit the loop, we have found new positions
                    break
            else:
                if __debug__:
                    log.info('single reaction; placing products failed.')
                moved_reactant = self.move_particle(reactant, reactant_pos)
                products = [moved_reactant]
                self.rejected_moves += 1

        else:
            raise RuntimeError('num products >= 3 not supported.')

        return products

    def fire_interaction(self, single, reactant_pos):
        # This takes care of the identity change when a particle associates with a surface
        # It performs the reactions:
        # A(3D) + surface -> 0
        # A(3D) + surface -> B(surface)

        # 0. get reactant info
        reactant        = single.pid_particle_pair
        reactant_radius    = reactant[1].radius
        rr = single.interactionrule

        # 1. remove the particle

        if len(rr.products) == 0:
            
            # 1. No products, no info
            # 2. No space required
            # 3. No space required
            # 4. process the changes (remove particle, make new ones)
            self.world.remove_particle(reactant[0])
            products = []

            # 5. No new single to be made
            # 6. Log the change

            
        elif len(rr.products) == 1:
            
            # 1. get product info
            product_species = self.world.get_species(rr.products[0])
            product_radius  = product_species.radius
            product_surface = single.surface

            # 1.5 get new position of particle
            transposed_pos = self.world.cyclic_transpose(reactant_pos, product_surface.shape.position)
            product_pos, _ = product_surface.projected_point(transposed_pos)
            product_pos = self.world.apply_boundary(product_pos)        # not sure this is still necessary

            # 2. burst the volume that will contain the products.
            #    Note that an interaction domain is already sized such that it includes the
            #    reactant particle moving into the surface.
            if product_radius > reactant_radius:
                self.burst_volume(product_pos, product_radius)

            # 3. check that there is space for the products 
            if self.world.check_overlap((product_pos, product_radius), reactant[0]):
                if __debug__:
                    log.info('interaction; placing product failed.')
                moved_reactant = self.move_particle(reactant, reactant_pos)
                products = [moved_reactant]     # no change 
                self.rejected_moves += 1

            else:
                # 4. process the changes (remove particle, make new ones)
                self.world.remove_particle(reactant[0])
                newparticle = self.world.new_particle(product_species.id, product_pos)
                products = [newparticle]

                # 5. update counters
                self.reaction_events += 1
                self.last_reaction = (rr, (reactant[1], None), products)

                # 6. Log the change
                if __debug__:
                     log.info('product (%s) = %s' % (len(products), products))

        else:
            raise RuntimeError('fire_interaction: num products > 1 not supported.')

        return products


    def fire_pair_reaction(self, pair, reactant1_pos, reactant2_pos):
        # This takes care of the identity change when two particles react with each other
        # It performs the reactions:
        # A(any structure) + B(same structure) -> 0
        # A(any structure) + B(same structure or 3D) -> C(same structure)

        # 0. get reactant info
        pid_particle_pair1 = pair.pid_particle_pair1
        pid_particle_pair2 = pair.pid_particle_pair2
        reactant1_species_id = pid_particle_pair1[1].sid
        reactant2_species_id = pid_particle_pair2[1].sid

        rr = pair.draw_reaction_rule(pair.interparticle_rrs)
        # 1. remove the particles

        if len(rr.products) == 0:

            # 1. no product particles, no info
            # 2. No space required
            # 3. No space required
            # 4. process the change (remove particles, make new ones)
            self.world.remove_particle(pid_particle_pair1[0])
            self.world.remove_particle(pid_particle_pair2[0])
            products = []

            # 5. No new single to be made
            # 6. Log the change
            self.last_reaction = (rr, (pid_particle_pair1[1], pid_particle_pair2[1]), products)

        elif len(rr.products) == 1:

            # 1. get product info
            product_species = self.world.get_species(rr.products[0])
            product_radius  = product_species.radius

            # 1.5 get new position for product particle
            product_pos = pair.draw_new_com (pair.dt, pair.event_type)
            product_pos = self.world.apply_boundary(product_pos)

            # 2.  make space for products
            # 2.1 if the product particle sticks out of the shell
#            if self.world.distance(old_com, product_pos) > (pair.get_shell_size() - product_radius):
            self.burst_volume(product_pos, product_radius)

            # 3. check that there is space for the products ignoring the reactants
            if self.world.check_overlap((product_pos, product_radius), pid_particle_pair1[0],
                                                                       pid_particle_pair2[0]):
                if __debug__:
                    log.info('fire_pair_reaction: no space for product particle.')
                moved_reactant1 = self.move_particle(pid_particle_pair1, reactant1_pos)
                moved_reactant2 = self.move_particle(pid_particle_pair2, reactant2_pos)
                products = [moved_reactant1, moved_reactant2]     # no change 
                self.rejected_moves += 1

            else:
                # 4. process the changes (remove particle, make new ones)
                self.world.remove_particle(pid_particle_pair1[0])
                self.world.remove_particle(pid_particle_pair2[0])
                newparticle = self.world.new_particle(product_species.id, product_pos)
                products = [newparticle]

                # 5. update counters
                self.reaction_events += 1
                self.last_reaction = (rr, (pid_particle_pair1, pid_particle_pair2),
                                      products)

                # 6. Log the change
                if __debug__:
                    log.info('fire_pair_reaction: product (%s) = %s' % (len(products), products))

        else:
            raise NotImplementedError('num products >= 2 not supported.')

        return products


    def fire_move(self, single, reactant_pos, ignore=None):
        # No reactions/Interactions have taken place -> no identity change of the particle
        # Just only move the particles

        reactant = single.pid_particle_pair

        if ignore:
            if self.world.check_overlap((reactant_pos, reactant[1].radius),
                                        reactant[0], ignore[0]):
                raise RuntimeError('fire_move: particle overlap failed.')
        else:
            if self.world.check_overlap((reactant_pos, reactant[1].radius),
                                        reactant[0]):
                raise RuntimeError('fire_move: particle overlap failed.')

        moved_reactant = self.move_particle(reactant, reactant_pos)
        return [moved_reactant]



    def make_new_domain(self, single):
        ### Make a new domain out of a NonInteractionSingle that was
        ### just fired

        # can only make new domain from NonInteractionSingle
        assert isinstance(single, NonInteractionSingle)
        
        # Special case: When D=0, nothing needs to change and only new event
        # time is drawn
        # TODO find nicer construction than just this if statement
        if single.getD() == 0:
            single.dt, single.event_type = single.determine_next_event() 
            single.last_time = self.t
            self.add_domain_event(single)
            return single


        # 0. Get generic info
        single_pos = single.pid_particle_pair[1].position
        single_radius = single.pid_particle_pair[1].radius

        # 1.0 Get neighboring domains and surfaces
        neighbors = self.geometrycontainer.get_neighbor_domains(single_pos, self.domains, ignore=[single.domain_id, ])
        # Get also surfaces but only if the particle is in 3D
        surface_distances = []
        if isinstance(single, SphericalSingle):
            surface_distances = self.geometrycontainer.get_neighbor_surfaces(single_pos, ignores=[single.structure.id])


        # Check if there are shells with the burst radius (reaction_threshold)
        # of the particle (intruders). Note that we approximate the reaction_volume
        # with a sphere (should be cylinder for 2D or 1D particle)
        reaction_threshold = single_radius * self.SINGLE_SHELL_FACTOR


        # TODO redundant
        intruders = self.geometrycontainer.filter_distance(neighbors, reaction_threshold)
        if __debug__:
            log.debug("intruders: %s" %
                      (', '.join(str(i) for i in intruders)))


        # 2.2 Burst the shells with the following conditions
        # -shell is in the reaction_threshold (intruder)
        # -shell is burstable (not Multi)
        # -shell is not newly made
        # -shell is not already a zero shell (just bursted)
        neighbor_distances = self.burst_non_multis(neighbors, reaction_threshold)


        # We prefer to make NonInteractionSingles for efficiency.
        # But if objects (Shells and surfaces) get close enough (closer than
        # reaction_threshold) we want to try Interaction and/or Pairs

        # 2.3 From the updated list of neighbors, we calculate the thresholds (horizons) for trying to form
        #     the various domains. The domains can be:
        # -zero-shell NonInteractionSingle -> Pair or Multi
        # -Surface -> Interaction
        # -Multi -> Multi
        # Note that we do not differentiate between directions. This means that we
        # look around in a sphere, the horizons are spherical
        pair_interaction_partners = []
        for domain, _ in neighbor_distances:
            if (isinstance (domain, NonInteractionSingle) and domain.is_reset()):
                # distance from the center of the particles/domains
                pair_distance = self.world.distance(single_pos, domain.shell.shape.position)
                pair_horizon  = (single_radius + domain.pid_particle_pair[1].radius) * self.SINGLE_SHELL_FACTOR
                pair_interaction_partners.append((domain, pair_distance - pair_horizon))
        
        for surface, surface_distance in surface_distances:
            surface_horizon = single_radius * self.SINGLE_SHELL_FACTOR
            pair_interaction_partners.append((surface, surface_distance - surface_horizon))

        pair_interaction_partners = sorted(pair_interaction_partners, key=lambda domain_overlap: domain_overlap[1])

        log.debug('pair_interaction_partners: %s' % str(pair_interaction_partners))

        # For each particles/particle or particle/surface pair, check if a pair or interaction can be
        # made. Note that we check the closest one first and only check if the object are within the horizon
        domain = None
        for obj, hor_overlap in pair_interaction_partners:

            if hor_overlap > 0.0 or domain:
                # there are no more potential partners (everything is out of range)
                # or a domain was formed successfully previously
                break

            if isinstance(obj, NonInteractionSingle):
                # try making a Pair (can be Mixed Pair or Normal Pair)
                domain = self.try_pair (single, obj)

            elif isinstance(obj, CylindricalSurface) or isinstance(obj, PlanarSurface):
                # try making an Interaction
                domain = self.try_interaction (single, obj)


        if not domain:
            # No Pair or Interaction domain was formed
            # Now choose between NonInteractionSingle and Multi

            # make a new list of potential partners
            # Is now zero-dt NonInteractionSingles and Multi's
            # TODO: combine with first selection such that we have to do this loop only once
            multi_partners = []
            for domain, dist_to_shell in neighbor_distances:
                if (isinstance (domain, NonInteractionSingle) and domain.is_reset()):
                    multi_horizon = (single_radius + domain.pid_particle_pair[1].radius) * self.MULTI_SHELL_FACTOR
                    distance = self.world.distance(single_pos, domain.shell.shape.position)
                    multi_partners.append((domain, distance - multi_horizon))

                elif isinstance(domain, Multi):
                    # The dist_to_shell = dist_to_particle - multi_horizon_of_target_particle
                    # So only the horizon and distance of the current single needs to be taken into account
                    # Note: this is built on the assumption that the shell of a Multi has the size of the horizon.
                    multi_horizon = (single_radius * self.MULTI_SHELL_FACTOR)
                    multi_partners.append((domain, dist_to_shell - multi_horizon))

            # Also add surfaces
            for surface, distance in surface_distances:
                surface_horizon = single_radius * self.MULTI_SHELL_FACTOR
                multi_partners.append((surface, distance - surface_horizon))


            # get the closest object (if there)
            if multi_partners:
                multi_partners = sorted(multi_partners, key=lambda domain_overlap: domain_overlap[1])
#                log.debug('multi_partners: %s' % str(multi_partners))
                closest_overlap = multi_partners[0][1]
            else:
                # In case there is really nothing
                closest_overlap = numpy.inf

            log.debug('Single or Multi: closest_overlap: %s' % (FORMAT_DOUBLE % closest_overlap))

            # If the closest partner is within the multi horizon we do Multi, otherwise Single
            if closest_overlap > 0.0: 
                # just make a normal NonInteractionSingle
                self.update_single(single)
                self.add_domain_event(single)
            else:
                # An object was closer than the Multi horizon
                # Form a multi with everything that is in the multi_horizon
                domain = self.form_multi(single, multi_partners)

        return domain


    def calculate_simplepair_shell_size(self, single1, single2):
        assert single1.structure == single2.structure

        # 0. Get some necessary information
        com, iv = SimplePair.do_transform(single1, single2, self.world)
        r0 = length (iv)

        # 1. Get the minimal possible shell size (including margin?)
        min_shell_size = SimplePair.get_min_shell_size(single1, single2, self.geometrycontainer)

        # 2. Get the maximum possible shell size
        max_shell_size = SimplePair.get_max_shell_size(com, single1, single2,
                                                       self.geometrycontainer, self.domains)

        # 3. Calculate the maximum based on some other criterium (convergence?)
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius
        sigma = radius1 + radius2
        distance_from_sigma = r0 - sigma        # the distance between the surfaces of the particles
#        convergence_max = distance_from_sigma * 100 + sigma + shell_size_margin                # FIXME
#       max_shell_size = min(max_shell_size, convergence_max)


        # 4. Check if min shell size for the Pair not larger than max shell size or 
        # sim cell size.
        if min_shell_size >= max_shell_size:
            if __debug__:
                log.debug('%s not formed: min_shell_size %s >='
                          'max_shell_size %s' %
                          ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                             single2.pid_particle_pair[0]),
                           FORMAT_DOUBLE % min_shell_size,
                           FORMAT_DOUBLE % max_shell_size))
            return None, None, None

        # 5. Calculate the 'ideal' shell size, it matches the expected first-passage time of the closest particle
        # with the expected time of the particles in the pair.
#        D1 = single1.pid_particle_pair[1].D
#        D2 = single2.pid_particle_pair[1].D
#        D12 = D1 + D2
#
#        D_closest = closest.pid_particle_pair[1].D
#        D_tot = D_closest + D12
#        ideal_shell_size = min((D12 / D_tot) *
#                            (closest_particle_distance - min_shell_size - closest_min_radius) +
#                            min_shell_size,
        ideal_shell_size = numpy.inf
        shell_size = min(max_shell_size, ideal_shell_size)


        # 6. Check if singles would not be better.
        # TODO clear this up -> strange rule
        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
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

        assert shell_size >= min_shell_size
        assert shell_size <= max_shell_size

        if __debug__:
            log.info('SimplePair shell can be made. shell_size=%s, '
                     'closest_shell_distance=%s,\nclosest = %s' %
                     (FORMAT_DOUBLE % shell_size, FORMAT_DOUBLE % 0, None))

        return com, r0, shell_size



    def calculate_mixedpair_shell_dimensions (self, single2D, single3D):


        # 1. Get the minimal possible shell size (including margin?)
        r_min, z_left_min, z_right_min = MixedPair.get_min_shell_dimensions(single2D, single3D, self.geometrycontainer)

        # 2. Get the maximum possible shell size
        r, z_left, z_right = MixedPair.get_max_shell_dimensions(single2D, single3D, self.geometrycontainer, self.domains)

        # Decide if MixedPair domain is possible.
        if r < r_min or z_left < z_left_min or z_right < z_right_min:
            if __debug__:
                log.debug('        *MixedPair not possible:\n'
                          '            %s +\n'
                          '            %s.\n'
                          '            r = %.3g. r_min = %.3g.\n'
                          '            z_left = %.3g. z_left_min = %.3g.\n'
                          '            z_right = %.3g. z_right_min = %.3g.' %
                          (single2D, single3D, r, r_min, z_left, z_left_min, 
                           z_right, z_right_min))
            return None, None, None, None, None
        else:

            # Compute origin, radius and half_length of cylinder.
            com, iv = MixedPair.do_transform(single2D, single3D, self.geometrycontainer.world)
            reference_point = com
            orientation_vector = normalize(single2D.structure.shape.unit_z * numpy.dot (iv, single2D.structure.shape.unit_z))
            r0 = length(iv)

            center = reference_point + ((z_right - z_left)/2.0) * orientation_vector
            center = self.world.apply_boundary(center)
            shell_half_length = (z_left + z_right) / 2.0
            shell_radius = r

            return center, shell_radius, shell_half_length, r0, orientation_vector



    def update_single(self, single): 
        assert isinstance(single, NonInteractionSingle) # This only works for 'simple' Singles

        min_shell_size = single.get_min_shell_size()
        max_shell_size = single.get_max_shell_size(self.geometrycontainer, self.domains)

        # Make sure that the new shell size is not too small or big
        new_shell_size = max(max_shell_size, min_shell_size)

        # Resize shell, don't change position.
        # Note: this should be done before determine_next_event.
        singlepos = single.pid_particle_pair[1].position
        self.update_single_shell(single, singlepos, new_shell_size)        

        single.dt, single.event_type = single.determine_next_event()
        single.last_time = self.t


    def process_single_event(self, single):

        ### log Single event
        if __debug__:
            log.info('FIRE SINGLE: %s' % single.event_type)
            log.info('single = %s' % single)



        if single.is_reset():
        ### If no event event really happened and just need to make new domain
            domains = [self.make_new_domain(single)]
            # domain is already scheduled in make_new_domain

        ### 1.1 Special cases (shortcuts)
        # In case nothing is scheduled to happen: do nothing and just reschedule
        elif single.dt == numpy.inf:
            self.add_domain_event(single)
            domains = [single]

        else:
        ### 1. Process 'normal' event produced by the single

            # check that the event time of the single (last_time + dt) is equal to the
            # simulator time
            assert (abs(single.last_time + single.dt - self.t) <= TIME_TOLERANCE * self.t), \
                'Timeline incorrect. single.last_time = %s, single.dt = %s, self.t = %s' % \
                (FORMAT_DOUBLE % single.last_time, FORMAT_DOUBLE % single.dt, FORMAT_DOUBLE % self.t)


            pid_particle_pair = single.pid_particle_pair
            # in case of an interaction domain: determine real event before doing anything
            if single.event_type == EventType.IV_EVENT:
                single.event_type = single.draw_iv_event_type()


            # get the (new) position
            if single.getD() != 0 and single.dt > 0.0:
                # If the particle had the possibility to diffuse
                newpos = single.draw_new_position(single.dt, single.event_type)
                newpos = self.world.apply_boundary(newpos)
            else:
                # no change in position has taken place
                newpos = pid_particle_pair[1].position
            # TODO? Check if the new positions are within domain

            # newpos now hold the new position of the particle (not yet committed to the world)
            # if we would here move the particles and make new shells, then it would be similar to a propagate

            self.remove_domain(single)

            # If the single had a decay reaction or interaction.
            if single.event_type == EventType.SINGLE_REACTION or \
               single.event_type == EventType.IV_INTERACTION:
                if __debug__:
                    log.info('%s' % single.event_type)
                    log.info('reactant = %s' % single)

                if single.event_type == EventType.SINGLE_REACTION:
                    self.single_steps[single.event_type] += 1       # TODO counters should also be updated for escape events
                    particles = self.fire_single_reaction(single, newpos)

                else:
                    self.interaction_steps[single.event_type] += 1  # TODO similarly here
                    particles = self.fire_interaction(single, newpos)

            else:
                particles = self.fire_move(single, newpos)

            # 5. Make a (new) domain (reuse) domain for each particle(s)
            #    (Re)schedule the (new) domain
            domains = []
            for pid_particle_pair in particles:
#               single.pid_particle_pair = pid_particle_pair    # reuse single
                single = self.create_single(pid_particle_pair)  # TODO re-use NonInteractionSingle domain if possible
                self.add_domain_event(single)                   # TODO re-use event if possible
                domains.append(single) 

            # 6. Log change?
                
# NOTE Code snippets to re-use the domain/event
#
#        if isinstance(single, InteractionSingle):
#            # When bursting an InteractionSingle, the domain changes from Interaction
#           # NonInteraction domain. This needs to be reflected in the event 
#            self.remove_event(single)
#            self.add_domain_event(newsingle)
#        else:
#            assert single == newsingle
#            self.update_domain_event(self.t, single)
#
#        particle_radius = single.pid_particle_pair[1].radius
#        assert newsingle.shell.shape.radius == particle_radius

        return domains


    def process_pair_event(self, pair):
        assert self.check_domain(pair)

        if __debug__:
            log.info('FIRE PAIR: %s' % pair.event_type)
            log.info('single1 = %s' % pair.single1)
            log.info('single2 = %s' % pair.single2)

        # check that the event time of the single (last_time + dt) is equal to the
        # simulator time
        assert (abs(pair.last_time + pair.dt - self.t) <= TIME_TOLERANCE * self.t), \
            'Timeline incorrect. pair.last_time = %s, pair.dt = %s, self.t = %s' % \
            (FORMAT_DOUBLE % pair.last_time, FORMAT_DOUBLE % pair.dt, FORMAT_DOUBLE % self.t)



        single1 = pair.single1
        single2 = pair.single2
        assert single1.domain_id not in self.domains
        assert single2.domain_id not in self.domains

        pid_particle_pair1 = pair.pid_particle_pair1
        pid_particle_pair2 = pair.pid_particle_pair2
        pos1 = pid_particle_pair1[1].position
        pos2 = pid_particle_pair2[1].position

        # TODO store old_iv, old_com, r0 in pair object -> useless to recalculate these things all the time
        old_com, old_iv = pair.do_transform(single1, single2, self.world)
        r0 = length(old_iv)

        if pair.event_type == EventType.IV_EVENT:
            # Draw actual pair event for iv at very last minute.
            pair.event_type = pair.draw_iv_event_type(r0)

        self.pair_steps[pair.event_type] += 1


        # Get new position of particles
        if pair.dt > 0.0:
            newpos1, newpos2 = pair.draw_new_positions(pair.dt, r0, old_iv, pair.event_type)
            newpos1 = self.world.apply_boundary(newpos1)
            newpos2 = self.world.apply_boundary(newpos2)

            # check that the particles do not overlap with any other particles in the world
            assert not self.world.check_overlap((newpos1, pid_particle_pair1[1].radius),
                                                pid_particle_pair1[0], pid_particle_pair2[0])
            assert not self.world.check_overlap((newpos2, pid_particle_pair2[1].radius),
                                                pid_particle_pair1[0], pid_particle_pair2[0])

        else:
            newpos1 = pid_particle_pair1[1].position
            newpos2 = pid_particle_pair2[1].position
        # TODO? Check if the new positions are within domain
        # TODO? some more consistency checking of the positions
##            assert self.check_pair_pos(pair, newpos1, newpos2, old_com,
##                                       pair.get_shell_size())


        # newpos1/2 now hold the new positions of the particles (not yet committed to the world)
        # if we would here move the particles and make new shells, then it would be similar to a propagate

        self.remove_domain(pair)

        # If identity changing processes have taken place
        # Four cases:
        #  1. Single reaction
        if pair.event_type == EventType.SINGLE_REACTION:
            reactingsingle = pair.reactingsingle

            if reactingsingle == single1:
                theothersingle = single2
                # first move the non-reacting particle
                particles = self.fire_move(single2, newpos2, pid_particle_pair1)
                # don't ignore the (moved) non-reacting particle
                particles.extend(self.fire_single_reaction(single1, newpos1))
            else:
                theothersingle = single1
                particles = self.fire_move(single1, newpos1, pid_particle_pair2)
                particles.extend(self.fire_single_reaction(single2, newpos2))

            if __debug__:
                log.info('reactant = %s' % reactingsingle)

            
            domains = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                single = self.create_single(pid_particle_pair)  # TODO reuse the non-reacting single domain
                self.add_domain_event(single)
                domains.append(single)

        #
        # 2. Pair reaction
        #
        elif pair.event_type == EventType.IV_REACTION:

            particles = self.fire_pair_reaction (pair, newpos1, newpos2)

            domains = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                single = self.create_single(pid_particle_pair)
                self.add_domain_event(single)
                domains.append(single)

        # Just moving the particles
        #  3a. IV escape
        #  3b. com escape
        elif(pair.event_type == EventType.IV_ESCAPE or
             pair.event_type == EventType.COM_ESCAPE or
             pair.event_type == EventType.BURST):

            particles = self.fire_move (single1, newpos1, pid_particle_pair2)
            particles.extend(self.fire_move (single2, newpos2, pid_particle_pair1))

            # make new NonInteractionSingeles and reschedule domains
            domains = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                single = self.create_single(pid_particle_pair)  # TODO reuse domains that were cached in the pair
                self.add_domain_event(single)
                domains.append(single)

        else:
            raise SystemError('process_pair_event: invalid event_type.')

# NOTE code snippit to reuse domains that were cached in pair domain
#        # re-use single domains for particles
#        # normally we would take the particles + new positions from the old pair
#        # and use them to make new singles. Now we re-use the singles stored in the
#        # pair
#        single1 = pair.single1
#        single2 = pair.single2
#        assert single1.domain_id not in self.domains
#        assert single2.domain_id not in self.domains
#        single1.pid_particle_pair = pid_particle_pair1  # this is probably redundant
#        single2.pid_particle_pair = pid_particle_pair2
#
#        # 'make' the singles
#        single1.initialize(self.t)
#        single2.initialize(self.t)
#        
#        self.update_single_shell(single1, newpos1, pid_particle_pair1[1].radius)
#        self.update_single_shell(single2, newpos2, pid_particle_pair2[1].radius)
#
#
#        self.domains[single1.domain_id] = single1
#        self.domains[single2.domain_id] = single2
#
#        # Check the dimensions of the shells of the singles with the shell in the container
#        if __debug__:
#            container1 = self.geometrycontainer.get_container(single1.shell)
#            assert container1[single1.shell_id].shape.radius == \
#                   single1.shell.shape.radius
#            if type(single1.shell) is CylindricalShell:
#                assert container1[single1.shell_id].shape.half_length == \
#                       single1.shell.shape.half_length
#
#            container2 = self.geometrycontainer.get_container(single2.shell)
#            assert container2[single2.shell_id].shape.radius == \
#                   single2.shell.shape.radius
#            if type(single2.shell) is CylindricalShell:
#                assert container2[single2.shell_id].shape.half_length == \
#                       single2.shell.shape.half_length
#
#        assert single1.shell.shape.radius == pid_particle_pair1[1].radius
#        assert single2.shell.shape.radius == pid_particle_pair2[1].radius
#        # even more checking
#        assert self.check_domain(single1)
#        assert self.check_domain(single2)
#        # Now finally we are convinced that the singles were actually made correctly


        # Log the event
        if __debug__:
            log.debug("process_pair_event: #1 { %s: %s => %s }" %
                      (single1, str(pos1), str(newpos1)))
            log.debug("process_pair_event: #2 { %s: %s => %s }" %
                      (single2, str(pos2), str(newpos2)))

        return domains

    def process_multi_event(self, multi):
        self.multi_steps[3] += 1  # multi_steps[3]: total multi steps
        multi.step()

        if __debug__:
            log.info('FIRE MULTI: %s' % multi.last_event)

        if(multi.last_event == EventType.MULTI_UNIMOLECULAR_REACTION or
           multi.last_event == EventType.MULTI_BIMOLECULAR_REACTION):
            self.reaction_events += 1
            self.last_reaction = multi.last_reaction

        if multi.last_event is not None:                # if an event took place
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
        # Sets next event time of 'single' to current time and event to burst event
        # returns a new single that is the result of firing the old single
        # new single also has proper zero shell and is rescheduled

        if __debug__:
            log.debug('burst single: %s', single)

        # Check correct timeline ~ MW
        assert single.last_time <= self.t
        assert self.t <= single.last_time + single.dt

        # Override dt, burst happens before single's scheduled event.
        # Override event_type. 
        single.dt = self.t - single.last_time
        single.event_type = EventType.BURST
        self.remove_event(single)

        newsingles = self.process_single_event(single)

        return newsingles

    def burst_pair(self, pair):
        # Sets next event time of 'pair' to current time and event to burst event
        # returns two new singles that are the result of firing the pair
        # new singles also have proper zero shell and are rescheduled

        if __debug__:
            log.debug('burst_pair: %s', pair)

        # make sure that the current time is between the last time and the event time
        assert pair.last_time <= self.t
        assert self.t <= pair.last_time + pair.dt

        # Override event time and event_type. 
        pair.dt = self.t - pair.last_time 
        pair.event_type = EventType.BURST
        self.remove_event(pair)         # remove the event -> was still in the scheduler

        newsingles = self.process_pair_event(pair)

        return newsingles


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
        # Don't tweak dz_left, this will break the conditions below for membrane interaction

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
        origin = self.world.apply_boundary(origin)
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

        assert self.check_domain(interaction)
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


    def try_pair(self, single1, single2):
        if __debug__:
            log.debug('trying to form Pair(%s, %s)' %
                (single1.pid_particle_pair, single2.pid_particle_pair))

        assert single1.is_reset()
        assert single2.is_reset()


        # Try forming a Pair only if singles are on same structure.
        if single1.structure == single2.structure:
            # particles are on the same structure

            center, r0, shell_size = self.calculate_simplepair_shell_size (single1, single2)
            if shell_size:
                # A shell could be made and makes sense. Create a Pair
                pair = self.create_pair(single1, single2, center, shell_size, r0)
                if __debug__:
                    log.debug('Created: %s, shell_size = %.3g, r0 = %.3g' %
                              (pair, shell_size, r0))

            else:
                pair = None

        elif (isinstance(single1.structure, PlanarSurface) and isinstance(single2.structure, CuboidalRegion)) ^ \
             (isinstance(single2.structure, PlanarSurface) and isinstance(single1.structure, CuboidalRegion)):
            # one particle in 2D, the other in 3D

            if __debug__:
                log.debug('Mixed pair')

            # make sure that we know which single is on what structure
            if isinstance (single1.structure, PlanarSurface):
                single2D = single1
                single3D = single2
            else:
                single2D = single2
                single3D = single1

            center, shell_radius, shell_half_length, r0, shell_orientation_vector = \
            self.calculate_mixedpair_shell_dimensions (single2D, single3D)

            if shell_radius:
                pair = self.create_pair(single2D, single3D, center, shell_radius, r0,
                                        shell_half_length, shell_orientation_vector)
                if __debug__:
                    log.debug('Created: %s, shell_radius = %.3g, shell_half_length = %.3g, r0 = %.3g' %
                              (pair, shell_radius, shell_half_length, r0))

            else:
                pair = None

        else:
            # a 1D/3D pair was supposed to be formed -> unsupported
            if __debug__:
                log.debug('Pair(%s, %s) not formed: combination of structures not supported.' %
                          (single1.pid_particle_pair[0],
                           single2.pid_particle_pair[0]))
            pair = None


        # if a pair has been formed
        if pair:

            pair.dt, pair.event_type, pair.reactingsingle = \
            pair.determine_next_event(r0)

            assert pair.dt >= 0

            self.last_time = self.t

            self.remove_domain(single1)
            self.remove_domain(single2)

            # single1 will be removed by the scheduler.
            self.remove_event(single2)

            assert self.check_domain(pair)
            self.add_domain_event(pair)

            if __debug__:
                log.info('dt=%s' % (FORMAT_DOUBLE % pair.dt)) 

            return pair
        else:
            return None
    

    def form_multi(self, single, multi_partners):
        # form a Multi with the 'single'

        # The 'multi_partners' are neighboring NonInteractionSingles, Multis and surfaces which
        # can be added to the Multi (this can also be empty)
        for partner, overlap in multi_partners:
            assert ((isinstance(partner, NonInteractionSingle) and partner.is_reset()) or \
                    isinstance(partner, Multi) or \
                    isinstance(partner, PlanarSurface) or \
                    isinstance(partner, CylindricalSurface)), \
                    'multi_partner %s was not of proper type' % (partner)


        # Only consider objects that are within the Multi horizon
        # assume that the multi_partners are already sorted by overlap
        neighbors = [obj for (obj, overlap) in multi_partners if overlap < 0]

        # TODO there must be a nicer way to do this
        if neighbors:
            closest = neighbors[0]
        else:
            closest = None


        # 1. Make new Multi if Necessary (optimization)
        if isinstance(closest, Multi):
        # if the closest to this Single is a Multi, reuse the Multi.
            multi = closest
            neighbors = neighbors[1:]   # don't add the closest multi with add_to_multi_recursive below
        else: 
        # the closest is not a Multi. Can be NonInteractionSingle, surface or
        # nothing. Create new Multi
            multi = self.create_multi()

        # 2. Add the single to the Multi
        self.add_to_multi(single, multi)
        self.remove_domain(single)


        # Add all partner domains (multi and dt=0 NonInteractionSingles) in the horizon to the multi
        for partner in neighbors:
            if (isinstance(partner, NonInteractionSingle) and partner.is_reset()) or \
               isinstance(partner, Multi):
                self.add_to_multi_recursive(partner, multi)
            else:
                # neighbor is a surface
                log.debug('add_to_multi: object(%s) is surface, not added to multi.' % partner)


        # 3. Initialize the multi and (re-)schedule
        multi.initialize(self.t)
        if isinstance(closest, Multi):
            # Multi existed before
            self.update_domain_event(self.t + multi.dt, multi)
        else:
            # In case of new multi
            self.add_domain_event(multi)

        return multi


    def add_to_multi_recursive(self, domain, multi):
    # adds 'domain' (which can be zero-dt NonInteractionSingle or Multi) to 'multi'


        if __debug__:
            log.debug('add_to_multi_recursive.\n domain: %s, multi: %s' % (domain, multi))

        if isinstance(domain, NonInteractionSingle):
            assert domain.is_reset()        # domain must be zero-dt NonInteractionSingle

            # check that the particles were not already added to the multi previously
            if multi.has_particle(domain.pid_particle_pair[0]):
                # Already in the Multi.
                if __debug__:
                    log.debug('%s already in multi. skipping.' % domain)
                return

            
            # 1. Get the neighbors of the domain
            dompos = domain.shell.shape.position
            neighbor_distances = self.geometrycontainer.get_neighbor_domains(dompos, self.domains,
                                                                             ignore=[domain.domain_id, multi.domain_id])

            # 2. Burst the domains that interfere with making the multi shell
            burstradius = domain.pid_particle_pair[1].radius * self.MULTI_SHELL_FACTOR
            neighbor_distances = self.burst_non_multis(neighbor_distances, burstradius)
            # This bursts only domains in which time has passed, it is assumed that other domains
            # have been made 'socially', meaning that they leave enough space for this particle to make a multi
            # shell without overlapping. 

            # neighbor_distances has same format as before
            # It contains the updated list of domains in the neighborhood with distances

            # 3. add domain to the multi
            self.add_to_multi(domain, multi)
            self.remove_domain(domain)
            self.remove_event(domain)

            # 4. Add recursively to multi, neighbors that:
            #    - have overlapping 'multi_horizons'
            #    - are Multi or
            #    - are just bursted (initialized) NonInteractionSingles
            for neighbor, dist_to_shell in neighbor_distances:
                if (isinstance (neighbor, NonInteractionSingle) and neighbor.is_reset()):
                    multi_horizon = (domain.pid_particle_pair[1].radius + neighbor.pid_particle_pair[1].radius) * \
                                    self.MULTI_SHELL_FACTOR
                    # distance from the center of the particles/domains
                    distance = self.world.distance(dompos, neighbor.shell.shape.position)
                    overlap = distance - multi_horizon

                elif isinstance(neighbor, Multi):
                    # The dist_to_shell = dist_to_particle - multi_horizon_of_target_particle
                    # So only the horizon and distance of the current single needs to be taken into account
                    # Note: this is built on the assumption that the shell of a Multi has the size of the horizon.
                    multi_horizon = (domain.pid_particle_pair[1].radius * self.MULTI_SHELL_FACTOR)
                    overlap = dist_to_shell - multi_horizon
                else:
                    # neighbor is not addible to multi
                    overlap = numpy.inf

                # 4.2 if the domains are within the multi horizon
                if overlap < 0:
                    self.add_to_multi_recursive(neighbor, multi)


        elif isinstance(domain, Multi):

            # check that the particles were not already added to the multi previously
            for pp in multi.particles:
                if domain.has_particle(pp[0]):
                    if __debug__:
                        log.debug('%s already added. skipping.' % domain)
                    break
            else:
                # 1.   No bursting is needed, since the shells in the multi already exist and no space has be made
                # 2.   Add the domain (the partner multi) to the existing multi
                self.merge_multis(domain, multi)

                # 3/4. No neighbors are searched and recursively added since everything that was in the multi
                #      horizon of the domain was already in the multi.

        else:
            # only NonInteractionSingles and Multi should be selected
            assert False, 'Trying to add non Single or Multi to Multi.'


    def add_to_multi(self, single, multi):
    # This method is similar to the 'create_' methods for pair and interaction, except it's now for a multi
    # adding a single to the multi instead of creating a pair or interaction out of single(s)

        if __debug__:
            log.info('add to multi:\n  %s\n  %s' % (single, multi))

        shell_id = self.shell_id_generator()
        shell = multi.create_new_shell(single.pid_particle_pair[1].position,
                single.pid_particle_pair[1].radius * self.MULTI_SHELL_FACTOR)
        shell_id_shell_pair = (shell_id, shell)
        self.geometrycontainer.move_shell(shell_id_shell_pair)
        multi.add_shell(shell_id_shell_pair)
        multi.add_particle(single.pid_particle_pair)

        # TODO maybe also cache the singles in the Multi for reuse?

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
        # Note: it returns the distance between pos and the surface of the shell
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
\tSingle:\t%d\t(escape: %d, reaction: %d, bursted: %d)
\tInteraction: %d\t(escape: %d, interaction: %d, bursted: %d)
\tPair:\t%d\t(escape r: %d, R: %d, reaction pair: %d, single: %d, bursted: %d)
\tMulti:\t%d\t(escape: %d, reaction pair: %d, single: %d)
total reactions = %d
rejected moves = %d
''' \
            % (self.t, self.step_counter,
               numpy.array(self.single_steps.values()).sum(),
               self.single_steps[EventType.SINGLE_ESCAPE],
               self.single_steps[EventType.SINGLE_REACTION],
               self.single_steps[EventType.BURST],
               numpy.array(self.interaction_steps.values()).sum(),
               self.interaction_steps[EventType.IV_ESCAPE],
               self.interaction_steps[EventType.IV_INTERACTION],
               self.interaction_steps[EventType.BURST],
               numpy.array(self.pair_steps.values()).sum(),
               self.pair_steps[EventType.IV_ESCAPE],
               self.pair_steps[EventType.COM_ESCAPE],
               self.pair_steps[EventType.IV_REACTION],
               self.pair_steps[EventType.SINGLE_REACTION],
               self.pair_steps[EventType.BURST],
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

    def check_domain(self, domain):
        domain.check()

        if isinstance(domain, Multi):
            # Ignore all surfaces, multi shells can overlap with 
            # surfaces.
            ignores = [s.id for s in self.world.structures]
        elif isinstance(domain, InteractionSingle) or isinstance(domain, MixedPair):
            # Ignore surface of the particle and interaction surface
            ignores = [domain.structure.id, domain.surface.id]
        else:
            ignores = [domain.structure.id]



        for shell_id, shell in domain.shell_list:

            # TODO should be replace by list of surface
            surface, distance = get_closest_surface(self.world, shell.shape.position, ignores)
            if surface:
                surfaces = [(surface, distance), ]
                for surface, distance in surfaces:
                    assert self.check_surface_overlap(shell, surface), \
                        '%s (%s) overlaps with %s.' % \
                        (str(domain), str(shell), str(surface))


            neighbors = self.geometrycontainer.get_neighbor_domains(shell.shape.position,
                                                                    self.domains, ignore=[domain.domain_id])
            # TODO maybe don't check all the shells, this takes a lot of time

            # testing overlap criteria
            for neighbor, _ in neighbors:
                for _, neighbor_shell in neighbor.shell_list:
                    overlap = self.check_shell_overlap(shell, neighbor_shell)
                    assert overlap >= 0.0, \
                        '%s (%s) overlaps with %s (%s) by %s.' % \
                        (domain, str(shell), str(neighbor), str(neighbor_shell), FORMAT_DOUBLE % overlap)


            # checking wether the shell don't exceed the maximum size
            if (type(shell.shape) is Cylinder):
                # the cylinder should at least fit in the maximal sphere
                shell_size = math.sqrt(shell.shape.radius**2 + shell.shape.half_length**2)
            elif (type(shell.shape) is Sphere):
                shell_size = shell.shape.radius
            else:
                raise RuntimeError('check_domain error: Shell shape was not Cylinder of Sphere')

            assert shell_size <= self.geometrycontainer.get_user_max_shell_size(), \
                '%s shell size larger than user-set max shell size, shell_size = %s, max = %s.' % \
                (str(shell_id), FORMAT_DOUBLE % shell_size, FORMAT_DOUBLE % self.geometrycontainer.get_user_max_shell_size())

            assert shell_size <= self.geometrycontainer.get_max_shell_size(), \
                '%s shell size larger than simulator cell size / 2, shell_size = %s, max = %s.' % \
                (str(shell_id), FORMAT_DOUBLE % shell_size, FORMAT_DOUBLE % self.geometrycontainer.get_max_shell_size())

        return True

    def check_shell_overlap(self, shell1, shell2):
    # Returns True if the shells DO NOT overlap
    # Returns False if the shells DO overlap
    # Return the amount of overlap (positive number means no overlap, negative means overlap)

        # overlap criterium when both shells are spherical
        if (type(shell1.shape) is Sphere) and (type(shell2.shape) is Sphere):
            distance = self.world.distance(shell1.shape.position, shell2.shape.position)
            return distance - (shell1.shape.radius + shell2.shape.radius)

        # overlap criterium when both shells are cylindrical 
        elif (type(shell1.shape) is Cylinder) and (type(shell2.shape) is Cylinder):
            # assuming that the cylinders are always parallel
            assert feq(abs(numpy.dot (shell1.shape.unit_z, shell2.shape.unit_z)), 1.0), \
                'shells are not oriented parallel, dot(unit_z1, unit_z2) = %s' % \
                (abs(numpy.dot (shell1.shape.unit_z, shell2.shape.unit_z)))

            shell1_pos = shell1.shape.position
            shell2_pos = shell2.shape.position
            shell2_post = self.world.cyclic_transpose(shell1_pos, shell2_pos)
            inter_pos = shell1_pos - shell2_pos
            inter_pos_z = shell1.shape.unit_z * numpy.dot(inter_pos, shell1.shape.unit_z)
            inter_pos_r = inter_pos - inter_pos_z
            overlap_r = length(inter_pos_r) - (shell1.shape.radius + shell2.shape.radius)
            overlap_z = length(inter_pos_z) - (shell1.shape.half_length + shell2.shape.half_length)
            if (overlap_r < 0.0) and (overlap_z < 0.0):
                return -1           # TODO: find better number for overlap measure
            else:
                return 1


        # overlap criterium when one shell is spherical and the other is cylindrical
        elif (type(shell1.shape) is Sphere) and (type(shell2.shape) is Cylinder):
            distance = self.world.distance(shell2.shape, shell1.shape.position)
            return distance - shell1.shape.radius

        # the other way around
        elif (type(shell1.shape) is Cylinder) and (type(shell2.shape) is Sphere):
            distance = self.world.distance(shell1.shape, shell2.shape.position)
            return distance - shell2.shape.radius

        # something was wrong (wrong type of shell provided)
        else:
            raise RuntimeError('check_shell_overlap: wrong shell type(s) provided.')


    def check_surface_overlap(self, shell, surface):
        if (type(shell.shape) is Sphere):
            distance = self.world.distance(surface.shape, shell.shape.position)
            return shell.shape.radius < distance

        elif (type(shell.shape) is Cylinder):
            distance = self.world.distance(surface.shape, shell.shape.position)

            # Only take shells into account that have unit_z perpendicular to unit_z of surface
            if (type(surface.shape) is Plane):
                return shell.shape.radius < distance
            elif (type(surface.shape) is Cylinder):
                return shell.shape.half_length < distance
            else:
                raise RuntimeError('check_surface_overlap: Surface was not Plane of Cylinder.')

        # something was wrong (wrong type of shell provided)
        else:
            raise RuntimeError('check_surface_overlap: wrong shell type provided.')


    def check_domain_for_all(self):
        for id, event in self.scheduler:
            domain = self.domains[event.data]
            self.check_domain(domain)

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
    # checks that the events in the scheduler are consistent with the domains

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
        particle1 = pair.pid_particle_pair1[1]
        particle2 = pair.pid_particle_pair2[1]

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
        
        self.check_domain_for_all()

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
        (Single, process_single_event),
        (Pair, process_pair_event),
        (Multi, process_multi_event)
        ]


