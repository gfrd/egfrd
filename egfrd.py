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
    SphericalSurface,
    CylindricalSurface,
    DiskSurface,
    PlanarSurface,
    Surface,
    _random_vector,
    Cylinder,
    Disk,
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
from shells import (
    testShellError,
    testPair,
    testInteractionSingle,
    hasSphericalShell,
    hasCylindricalShell,
    SphericalSingletestShell,
    SphericalPairtestShell,
    PlanarSurfaceSingletestShell,
    PlanarSurfacePairtestShell,
    PlanarSurfaceTransitionSingletestShell,
    PlanarSurfaceTransitionPairtestShell,
    CylindricalSurfaceSingletestShell,
    CylindricalSurfacePairtestShell,
    DiskSurfaceSingletestShell,
    PlanarSurfaceInteractiontestShell,
    CylindricalSurfaceInteractiontestShell,
    CylindricalSurfaceCapInteractiontestShell,
    CylindricalSurfaceSinktestShell,
    MixedPair2D3DtestShell,
    MixedPair1DCaptestShell,
    )

import logging
import os

log = logging.getLogger('ecell')

### Singles
def create_default_single(domain_id, shell_id, pid_particle_pair, structure, reaction_rules, geometrycontainer, domains):
    if isinstance(structure, CuboidalRegion):
        # first make the test shell
        testSingle = SphericalSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
        return SphericalSingle          (domain_id, shell_id, testSingle, reaction_rules)
    elif isinstance(structure, PlanarSurface):
        # first make the test shell
        testSingle = PlanarSurfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
        return PlanarSurfaceSingle      (domain_id, shell_id, testSingle, reaction_rules)
    elif isinstance(structure, CylindricalSurface):
        # first make the test shell
        testSingle = CylindricalSurfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
        return CylindricalSurfaceSingle (domain_id, shell_id, testSingle, reaction_rules)
    elif isinstance(structure, DiskSurface):
        # first make the test shell
        testSingle = DiskSurfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
        return DiskSurfaceSingle (domain_id, shell_id, testSingle, reaction_rules)

### Interactions
def try_default_testinteraction(single, target_structure, geometrycontainer, domains):
    if isinstance(single.structure, CuboidalRegion):
        if isinstance(target_structure, PlanarSurface):
            return PlanarSurfaceInteractiontestShell(single, target_structure, geometrycontainer, domains)
        elif isinstance(target_structure, CylindricalSurface):
            return CylindricalSurfaceInteractiontestShell(single, target_structure, geometrycontainer, domains)
        else:
            raise testShellError('(Interaction). Combination of (3D particle, target_structure) is not supported')
    elif isinstance(single.structure, PlanarSurface):
        raise testShellError('(Interaction). Combination of (2D particle, target_structure) is not supported')
    elif isinstance(single.structure, CylindricalSurface):
        if isinstance(target_structure, DiskSurface):
            return CylindricalSurfaceCapInteractiontestShell (single, target_structure, geometrycontainer, domains)
        elif isinstance(target_structure, CylindricalSurface):
            return CylindricalSurfaceSinktestShell (single, target_structure, geometrycontainer, domains)
        else:
            raise testShellError('(Interaction). Combination of (1D particle, target_structure) is not supported')
    else:
        raise testShellError('(Interaction). structure of particle was of invalid type')

def create_default_interaction(domain_id, shell_id, testShell, reaction_rules, interaction_rules):
    if isinstance(testShell, CylindricalSurfaceInteractiontestShell):
        return CylindricalSurfaceInteraction    (domain_id, shell_id, testShell, reaction_rules, interaction_rules)
    elif isinstance(testShell, PlanarSurfaceInteractiontestShell):
        return PlanarSurfaceInteraction         (domain_id, shell_id, testShell, reaction_rules, interaction_rules)
    elif isinstance(testShell, CylindricalSurfaceCapInteractiontestShell):
        return CylindricalSurfaceCapInteraction (domain_id, shell_id, testShell, reaction_rules, interaction_rules)
    elif isinstance(testShell, CylindricalSurfaceSinktestShell):
        return CylindricalSurfaceSink           (domain_id, shell_id, testShell, reaction_rules, interaction_rules)

### Transitions
def try_default_testtransition(single, target_structure, geometrycontainer, domains):
    if isinstance(single.structure, CuboidalRegion):        
            raise testShellError('(Transition). Combination of (3D particle, target_structure) is not supported')
    elif isinstance(single.structure, PlanarSurface):
        if isinstance(target_structure, PlanarSurface):
            return PlanarSurfaceTransitionSingletestShell(single, target_structure, geometrycontainer, domains)
        else:
            raise testShellError('(Transition). Combination of (2D particle, target_structure other than plane) is not supported')
    elif isinstance(single.structure, CylindricalSurface):        
            raise testShellError('(Transition). Combination of (1D particle, target_structure) is not supported')
    else:
        raise testShellError('(Transition). structure of particle was of invalid type')

def create_default_transition(domain_id, shell_id, testShell, reaction_rules):
    if isinstance(testShell, PlanarSurfaceTransitionSingletestShell):
        return PlanarSurfaceTransitionSingle       (domain_id, shell_id, testShell, reaction_rules)

### Pairs
def try_default_testpair(single1, single2, geometrycontainer, domains):
    if single1.structure == single2.structure:
        if isinstance(single1.structure, CuboidalRegion):
            return SphericalPairtestShell           (single1, single2, geometrycontainer, domains)
        elif isinstance(single1.structure, PlanarSurface):
            return PlanarSurfacePairtestShell       (single1, single2, geometrycontainer, domains)
        elif isinstance(single1.structure, CylindricalSurface):
            return CylindricalSurfacePairtestShell  (single1, single2, geometrycontainer, domains)
    elif (isinstance(single1.structure, PlanarSurface) and isinstance(single2.structure, PlanarSurface)):
        return PlanarSurfaceTransitionPairtestShell (single1, single2, geometrycontainer, domains) 
    #elif (isinstance(single1.structure, PlanarSurface) and isinstance(single2.structure, CuboidalRegion)):
        #return MixedPair2D3DtestShell               (single1, single2, geometrycontainer, domains) 
    #elif (isinstance(single2.structure, PlanarSurface) and isinstance(single1.structure, CuboidalRegion)):
        #return MixedPair2D3DtestShell(single2, single1, geometrycontainer, domains)
    elif (isinstance(single1.structure, CylindricalSurface) and isinstance(single2.structure, DiskSurface)):
        return MixedPair1DCaptestShell(single1, single2, geometrycontainer, domains)
    elif (isinstance(single2.structure, CylindricalSurface) and isinstance(single1.structure, DiskSurface)):
        return MixedPair1DCaptestShell(single2, single1, geometrycontainer, domains)
    else:
        # another mixed pair was supposed to be formed -> unsupported
        raise testShellError('(MixedPair). combination of structures not supported')

def create_default_pair(domain_id, shell_id, testShell, reaction_rules):
    # Either SphericalPair, PlanarSurfacePair, or CylindricalSurfacePair.
    if   isinstance(testShell, SphericalPairtestShell):
        return SphericalPair               (domain_id, shell_id, testShell, reaction_rules)
    elif isinstance(testShell, PlanarSurfacePairtestShell):
        return PlanarSurfacePair           (domain_id, shell_id, testShell, reaction_rules)
    elif isinstance(testShell, PlanarSurfaceTransitionPairtestShell):
        return PlanarSurfaceTransitionPair (domain_id, shell_id, testShell, reaction_rules)
    elif isinstance(testShell, CylindricalSurfacePairtestShell):
        return CylindricalSurfacePair      (domain_id, shell_id, testShell, reaction_rules)
    # or MixedPair (3D/2D or 1D/Cap)
    elif isinstance(testShell, MixedPair2D3DtestShell):
        return MixedPair2D3D               (domain_id, shell_id, testShell, reaction_rules)
    elif isinstance(testShell, MixedPair1DCaptestShell):
        return MixedPair1DCap              (domain_id, shell_id, testShell, reaction_rules)

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
        self.MAX_NUM_DT0_STEPS = 1000

        self.MAX_TIME_STEP = 10

        self.DEFAULT_DT_FACTOR = 1e-5           # Diffusion time prefactor in oldBD algortithm to determine time step.
        
        self.DEFAULT_STEP_SIZE_FACTOR = 0.05    # The maximum step size in the newBD algorithm is determined as DSSF * sigma_min.
                                                # Make sure that DEFAULT_STEP_SIZE_FACTOR < MULTI_SHELL_FACTOR, or else the 
                                                # reaction volume sticks out of the multi. 

        self.BD_ONLY_FLAG = False               # Will force the algorithm into Multi-creation, i.e. always to use BD
                                                # Take care: This is for testing only! Keep this 'False' for normal sims!

        self.BD_DT_HARDCORE_MIN = -1e-9         # This is to define a hardcore lower bound for the timestep that will be
                                                # dynamically determined by the new BD scheme. It will prevent the algorithm
                                                # to calculate ridiculously small timesteps, but will break detail balance.
                                                # Take care: This is for testing only! Keep this at a negative value for normal sims!

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

    #####################################################
    #### METHODS FOR GENERAL SIMULATOR FUNCTIONALITY ####
    #####################################################
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
                            EventType.MULTI_BIMOLECULAR_REACTION:0,
                            EventType.MULTI_DIFFUSION:0,
                            EventType.BURST:0,
                            3:0}
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

        #non_single_list = []

        # Burst all the domains that we know of.
        self.burst_all_domains()

        ##first burst all Singles, and put Pairs and Multis in a list.
        #for id, event in self.scheduler:
        #ignore = []
        #for domain_id, domain in self.domains.items():
            #if domain_id not in ignore:
                  ##domain = self.domains[event.data]
                #if isinstance(domain, Pair) or isinstance(domain, Multi):
                    #non_single_list.append(domain)
                #elif isinstance(domain, Single):
                    #if __debug__:
                        #log.debug('burst %s, last_time= %s' % 
                                  #(domain, FORMAT_DOUBLE % domain.last_time))
                    #ignore.append(domain.domain_id)
                    #ignore = self.burst_single(domain, ignore)
                #else:
                    #assert False, 'domain from domains{} was no Single, Pair or Multi'


        ## then burst all Pairs and Multis.
        #if __debug__:
            #log.debug('burst %s' % non_single_list)
        #non_single_ids = [non_single.domain_id for non_single in non_single_list]
        #self.burst_domains(non_single_ids, ignore)

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
                f(self, domain, [domain.domain_id])         # fire the correct method for the class (e.g. process_single_event(self, Single))

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
                self.zero_steps += 1
                # TODO What is best solution here? Usually no prob, -> just let 
                # user know?
                if self.zero_steps >= max(self.scheduler.size * 3, self.MAX_NUM_DT0_STEPS): 
                    raise RuntimeError('too many dt=zero steps. '
                                       'Simulator halted?'
                                    'dt= %.10g-%.10g' % (self.scheduler.top[1].time, self.t))
                    log.warning('dt=zero step, working in s.t >> dt~0 Python limit.')
            else:
                self.zero_steps = 0


    #####################################
    #### METHODS FOR DOMAIN CREATION ####
    #####################################
    def create_single(self, pid_particle_pair):
    # Create a new single domain from a particle.
    # The interaction can be any NonInteractionSingle (SphericalSingle, PlanarSurface or CylindricalSurface
    # NonInteractionSingle).

        # 1. generate identifiers for the domain and shell. The event_id is
        # generated by the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # get unimolecular reaction rules
        reaction_rules = self.network_rules.query_reaction_rule(pid_particle_pair[1].sid)
        # Get structure (region or surface) where the particle lives.
        species = self.world.get_species(pid_particle_pair[1].sid)
        structure = self.world.get_structure(pid_particle_pair[1].structure_id)

        # 2. Create and register the single domain.
        # The type of the single that will be created 
        # depends on the structure (region or surface) this particle is 
        # in/on. Either SphericalSingle, PlanarSurfaceSingle, or 
        # CylindricalSurfaceSingle.
        single = create_default_single(domain_id, shell_id, pid_particle_pair, 
                                       structure, reaction_rules,
                                       self.geometrycontainer, self.domains)

        assert isinstance(single, NonInteractionSingle)
        single.initialize(self.t)               # set the initial event and event time
        self.domains[domain_id] = single

        # 3. update the proper shell container
        self.geometrycontainer.move_shell(single.shell_id_shell_pair)

        #if __debug__:
            ## Used in __str__.
        single.world = self.world

        return single


    def create_interaction(self, testShell):
    # Create a new interaction domain from a testShell.
    # The interaction can be any interaction (PlanarSurface or CylindricalSurface). The creation of the testShell was
    # succesful, so here we can just make the domain.

        assert isinstance(testShell, testInteractionSingle)
        # TODO? Assert that the particle is not already associated with another domain

        # 1. generate identifiers for the domain and shell. event_id is generated by
        # the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # get unimolecular reaction rules
        species_id = testShell.pid_particle_pair[1].sid
        reaction_rules = self.network_rules.query_reaction_rule(species_id)
        # get reaction rules for interaction
        structure_type_id = testShell.target_structure.sid
        interaction_rules = self.network_rules.query_reaction_rule(species_id, structure_type_id)

        # 2. Create and register the interaction domain.
        # The type of the interaction that will be created 
        # depends on the surface (planar or cylindrical) the particle is 
        # trying to associate with. Either PlanarSurfaceInteraction or 
        # CylindricalSurfaceInteraction.
        interaction = create_default_interaction(domain_id, shell_id, testShell,
                                                 reaction_rules, interaction_rules)

        assert isinstance(interaction, InteractionSingle)
        interaction.initialize(self.t)
        self.domains[domain_id] = interaction

        # 3. update the shell containers 
        self.geometrycontainer.move_shell(interaction.shell_id_shell_pair)

        #if __debug__:
            ## Used in __str__.
        interaction.world = self.world

        return interaction


    def create_transition(self, testShell):
    # Create a new transition domain from a testShell.
    # Currently only plane-plane transitions are supported.
    # We assume here that the testShell was successfully created.

        assert isinstance(testShell, PlanarSurfaceTransitionSingletestShell) # TODO should be generalized

        # 1. generate identifiers for the domain and shell. event_id is generated by
        # the scheduler
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # get unimolecular reaction rules
        species_id = testShell.pid_particle_pair[1].sid    
        reaction_rules = self.network_rules.query_reaction_rule(species_id)

        # 2. Create and register the transition domain.
        transition = create_default_transition(domain_id, shell_id, testShell, reaction_rules)

        assert isinstance(transition, Single)
        transition.initialize(self.t)
        self.domains[domain_id] = transition

        # 3. update the shell containers 
        self.geometrycontainer.move_shell(transition.shell_id_shell_pair)
        # TODO What precisely happens here?

        transition.world = self.world

        return transition


    def create_pair(self, testShell):
    # Create a new pair domain from a testShell.
    # The pair can be any pair (Simple or Mixed). The creation of the testShell was
    # succesful, so here we can just make the domain.

        assert isinstance(testShell, testPair)

        # 1. generate needed identifiers
        domain_id = self.domain_id_generator()
        shell_id = self.shell_id_generator()

        # Select 1 reaction type out of all possible reaction types between the two particles.
        reaction_rules = self.network_rules.query_reaction_rule(testShell.pid_particle_pair1[1].sid,
                                                                testShell.pid_particle_pair2[1].sid)


        # 2. Create pair. The type of the pair that will be created depends
        # on the structure (region or surface) the particles are in/on.  
        pair = create_default_pair(domain_id, shell_id, testShell, reaction_rules)

        assert isinstance(pair, Pair)
        pair.initialize(self.t)
        self.domains[domain_id] = pair

        # 3. update the shell containers with the new shell
        self.geometrycontainer.move_shell(pair.shell_id_shell_pair)

        #if __debug__:
            ## Used in __str__.
        pair.world = self.world

        return pair


    def create_multi(self):
        # 1. generate necessary id's
        domain_id = self.domain_id_generator()

        # 2. Create and register domain object
        multi = Multi(domain_id, self, self.DEFAULT_STEP_SIZE_FACTOR, self.BD_DT_HARDCORE_MIN)
        self.domains[domain_id] = multi

        # 3. no shells are yet made, since these are added later
        # -> a multi can have 0 shells
        return multi


    ##############################################
    #### METHODS FOR PARTICLE / DOMAIN UPDATE ####
    ##############################################
   #def move_single(self, single, position, radius=None):
       #single.pid_particle_pair = self.move_particle(single.pid_particle_pair, position)
       #self.update_single_shell(single, position, radius)

   #def update_single_shell(self, single, shell):#position, radius=None):
       ## Reuse shell_id.
       #shell_id = single.shell_id

       ## Replace shell.
       #shell_id_shell_pair = (shell_id, shell) 

       #single.shell_id_shell_pair = shell_id_shell_pair
       #self.geometrycontainer.move_shell(shell_id_shell_pair)


    def move_particle(self, pid_particle_pair, position, structure_id):
        # Moves a particle in World and performs the required change of structure_id
        # based on an existing particle.

        new_pid_particle_pair = (pid_particle_pair[0],
                                 Particle(position,
                                          pid_particle_pair[1].radius,
                                          pid_particle_pair[1].D,
                                          pid_particle_pair[1].v,
                                          pid_particle_pair[1].sid,
                                          structure_id))

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


    #######################################
    #### METHODS FOR EVENT BOOKKEEPING ####
    #######################################
    # TODO This method can be made a method to the scheduler class
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


    # TODO This method can be made a method to the scheduler class
    def remove_event(self, event):
        if __debug__:
            log.info('remove_event: event=#%d' % event.event_id)
        del self.scheduler[event.event_id]


    # TODO This method can be made a method to the scheduler class
    def update_domain_event(self, t, domain):
        if __debug__:
            log.info('update_event: %s, event=#%d, t=%s' %
                     (domain.domain_id, domain.event_id, FORMAT_DOUBLE % t))
        self.scheduler.update((domain.event_id, DomainEvent(t, domain)))


    #####################################
    #### METHODS FOR DOMAIN BURSTING ####
    #####################################
    def burst_domain(self, domain, ignore):
    # Reduces 'domain' (Single, Pair or Multi) to 'zero_singles', singles
    # with the zero shell, and dt=0.
    # returns:
    # - list of zero_singles that was the result of the bursting
    # - updated ignore list

        domain_id = domain.domain_id

        if __debug__:
            log.info('burst_domain: %s' % domain)

        if isinstance(domain, Single):  # Single
            # TODO. Compare with gfrd.
            zero_singles, ignore = self.burst_single(domain, ignore)
        elif isinstance(domain, Pair):  # Pair
            zero_singles, ignore = self.burst_pair(domain, ignore)
        else:                           # Multi
            assert isinstance(domain, Multi)
            zero_singles, ignore = self.burst_multi(domain, ignore)

        if __debug__:
            # After a burst, the domain should be gone and should be on the ignore list.
            # Also the zero_singles should only be NonInteractionSingles
            assert domain_id not in self.domains
            assert domain_id in ignore
            assert all(isinstance(b, NonInteractionSingle) for b in zero_singles)
            log.info('zero_singles = %s' % ',\n\t  '.join(str(i) for i in zero_singles))

        return zero_singles, ignore


    def burst_single(self, single, ignore):
        # Bursts the 'single' domain and updates the 'ignore' list.
        # Returns:
        # - zero_singles that are the result of bursting the single (can be multiple
        #   because the burst is recursive). The zero_single have zero shell and are scheduled.
        # - the updated ignore list

        if __debug__:
            log.debug('burst single: %s', single)

        # Check correct timeline ~ MW
        assert single.last_time <= self.t
        assert self.t <= single.last_time + single.dt

        # Override dt, the burst happens before the single's scheduled event.
        # Override event_type. 
        single.dt = self.t - single.last_time
        single.event_type = EventType.BURST
        # with a burst there is always an associated event in the scheduler.
        # to simulate the natural occurence of the event we have to remove it from the scheduler
        self.remove_event(single)

        return self.process_single_event(single, ignore)


    def burst_pair(self, pair, ignore):
        # Bursts the 'pair' domain and updates the 'ignore' list.
        # Returns:
        # - zero_singles that are the result of bursting the pair (can be more than two
        #   since the burst is recursive). The zero_single have zero shell and are scheduled.
        # - the updated ignore list

        if __debug__:
            log.debug('burst_pair: %s', pair)

        # make sure that the current time is between the last time and the event time
        if __debug__:
            assert pair.last_time <= self.t
            assert self.t <= pair.last_time + pair.dt

        # Override event time and event_type. 
        pair.dt = self.t - pair.last_time 
        pair.event_type = EventType.BURST
        # with a burst there is always an associated event in the scheduler.
        # to simulate the natural occurence of the event we have to remove it from the scheduler
        self.remove_event(pair)

        return self.process_pair_event(pair, ignore)


    def burst_multi(self, multi, ignore):
        # Bursts the 'multi' domain and updates the 'ignore' list.
        # Returns:
        # - zero_singles that are the result of bursting the multi (note that the burst is
        #   recursive). The zero_single have zero shell and are scheduled.
        # - the updated ignore list

        if __debug__:
            log.debug('burst multi: %s', multi)

        #multi.sim.sync()
        # make sure that the current time is between the last time and the event time
        if __debug__:
            assert multi.last_time <= self.t
            assert self.t <= multi.last_time + multi.dt

        # with a burst there is always an associated event in the scheduler.
        # to simulate the natural occurence of the event we have to remove it from the scheduler
        self.remove_event(multi)        # The old event was still in the scheduler

        return self.break_up_multi(multi, ignore)


    def burst_volume(self, pos, radius, ignore=[]):
        """ Burst domains within a certain volume and give their ids.

        (Bursting means it propagates the particles within the domain until
        the current time, and then creates a new, minimum-sized domain.)

        Arguments:
            - pos: position of area to be "bursted"
            - radius: radius of spherical area to be "bursted"
            - ignore: ids of domains that should be ignored, none by default.
        """
        # TODO burst_volume now always assumes a spherical volume to burst.
        # It is inefficient to burst on the 'other' side of the membrane or to
        # burst in a circle, when you want to make a new shell in the membrane.
        # Then we may want to burst in a cylindrical volume, not a sphere.
        # Same for the 1D particles on rods.

        neighbor_ids = self.geometrycontainer.get_neighbors_within_radius_no_sort(pos, radius, ignore)

        return self.burst_domains(neighbor_ids, ignore)


    def burst_all_domains(self):
        # Bursts all the domains in the simulator

        all_domains = self.domains.items()
        all_domain_ids =[domain_id for domain_id, _ in all_domains]
        zero_singles, ignore = self.burst_domains(all_domain_ids, [])

        # Check if the bursting procedure was correct.
        if __debug__:
            assert len(zero_singles) == self.world.num_particles
            assert len(ignore) == (len(all_domains) + self.world.num_particles)


    def burst_domains(self, domain_ids, ignore=[]):
        # bursts all the domains that are in the list 'domain_ids'
        # ignores all the domain ids on 'ignore'
        # returns:
        # - list of zero_singles that was the result of the bursting
        # - updated ignore list

        zero_singles = []
        for domain_id in domain_ids:
            if domain_id not in ignore:
                domain = self.domains[domain_id]

                # add the domain_id to the list of domains that is already bursted (ignore list),
                # and burst the domain.
                ignore.append(domain_id)
                zero_singles2, ignore = self.burst_domain(domain, ignore)
                # add the resulting zero_singles from the burst to the total list of zero_singles.
                zero_singles.extend(zero_singles2)

        return zero_singles, ignore


    def burst_non_multis(self, pos, radius, ignore):
        # Recursively bursts the domains within 'radius' centered around 'pos'
        # Similarly to burst_volume but now doesn't burst all domains.
        # Returns:
        # - the updated list domains that are already bursted -> ignore list
        # - zero_singles that were the result of the burst

        # get the neighbors in the burstradius that are not already bursted.
        neighbor_ids = self.geometrycontainer.get_neighbors_within_radius_no_sort(pos, radius, ignore)

        zero_singles = []
        for domain_id in neighbor_ids:
            if domain_id not in ignore:                 # the list 'ignore' may have been updated in a
                domain = self.domains[domain_id]        # previous iteration of the loop

                if isinstance(domain, NonInteractionSingle) and domain.is_reset():
                    # If the domain was already bursted, but was not on the ignore list yet, put it there.
                    # NOTE: we assume that the domain has already bursted around it when it was made!
                    ignore.append(domain_id)

                elif (not isinstance(domain, Multi)) and (self.t != domain.last_time):
                    # Burst domain only if (AND):
                    # -within burst radius
                    # -not already bursted before (not on 'ignore' list)
                    # -domain is not zero dt NonInteractionSingle
                    # -domain is not a multi
                    # -domain has time passed
 
                    # add the domain_id to the list of domains that is already bursted (ignore list),
                    # and burst the domain.
                    ignore.append(domain_id)
                    zero_singles2, ignore = self.burst_domain(domain, ignore)
                    # add the resulting zero_singles from the burst to the total list of zero_singles.
                    zero_singles.extend(zero_singles2)

                #else:
                    # Don't burst domain if (OR):
                    # -domain is farther than burst radius
                    # -domain is already a burst NonInteractionSingle (is on the 'ignore' list)
                    # -domain is a multi
                    # -domain is domain in which no time has passed

                    # NOTE that the domain id is NOT added to the ignore list since it may be bursted at a later
                    # time

        return zero_singles, ignore


    ####################################
    #### METHODS FOR EVENT HANDLING ####
    ####################################
    def fire_single_reaction(self, single, reactant_pos, reactant_structure_id, ignore):
        # This takes care of the identity change when a single particle decays into
        # one or a number of other particles
        # It performs the reactions:
        # A(any structure) -> 0
        # A(any structure) -> B(same structure or 3D)
        # A(any structure) -> B(same structure) + C(same structure or 3D)
        if __debug__:
            assert isinstance(single, Single)

        # Note that the reactant_structure_id is the id of the structure on which the particle was located at the time of the reaction.
        if __debug__:
            assert isinstance(single, Single)

        # 0. get reactant info
        reactant                    = single.pid_particle_pair
        reactant_radius             = reactant[1].radius
        species                     = self.world.get_species(reactant[1].sid)
        reactant_structure_type_id  = species.structure_type_id
        reactant_structure          = self.world.get_structure(reactant_structure_id)
        rr = single.reactionrule

        # The zero_singles are the NonInteractionSingles that are the total result of the recursive
        # bursting process.
        zero_singles = []

        if len(rr.products) == 0:
            
            # 1. No products, no info
            # 1.5 No new position/structure_id required
            # 2. No space required
            # 3. No space required
            # 4. process the changes (remove particle, make new ones)
            self.world.remove_particle(reactant[0])
            products = []
            # zero_singles/ignore is unchanged.

            # 5. Update counters
            self.reaction_events += 1
            self.last_reaction = (rr, (reactant, None), products)

            # 6. Log the change
            if __debug__:
                 log.info('single reaction: product = None.')
            
        elif len(rr.products) == 1:
            
            # 1. get product info
            product_species           = self.world.get_species(rr.products[0])
            product_radius            = product_species.radius
            product_structure_type_id = product_species.structure_type_id

            # 1.5 get new position and structure_id of particle
            # If the particle falls off a surface (non CuboidalRegion)
            if product_structure_type_id != reactant_structure_type_id:
                assert (product_structure_type_id == self.world.get_def_structure_type_id())

                # produce a number of new possible positions for the product particle
                # TODO make these generators for efficiency
                product_pos_list = []
                if isinstance(reactant_structure, PlanarSurface):
                    if reactant_structure.shape.is_one_sided:
                        # unbinding always in direction of unit_z
                        directions = [1]
                    else:
                        # randomize unbinding direction
                        a = myrandom.choice(-1, 1)
                        directions = [-a,a]
                    # place the center of mass of the particle 'at contact' with the membrane
                    vector_length = (product_radius + 0.0) * (MINIMAL_SEPARATION_FACTOR - 1.0)  # the thickness of the membrane is 0.0
                    product_pos_list = [reactant_pos + vector_length * reactant_structure.shape.unit_z * direction \
                                        for direction in directions]

                elif isinstance(reactant_structure, CylindricalSurface):
                    vector_length = (product_radius + reactant_structure.shape.radius) * MINIMAL_SEPARATION_FACTOR
                    for _ in range(self.dissociation_retry_moves):
                        unit_vector3D = random_unit_vector()
                        unit_vector2D = normalize(unit_vector3D - 
                                        (reactant_structure.shape.unit_z * numpy.dot(unit_vector3D, reactant_structure.shape.unit_z)))
                        vector = reactant_pos + vector_length * unit_vector2D
                        product_pos_list.append(vector)

                elif isinstance(reactant_structure, DiskSurface):
                    # unbinding in direction of disk unit vector
                    #vector_length = 1.0*product_radius * MINIMAL_SEPARATION_FACTOR
                    #vector        = reactant_pos + vector_length * reactant_structure.shape.unit_z
                    #product_pos_list.append(vector)
                    # unbinding perpendicularly to disk unit vector (i.e. like on cylinder)
                    vector_length = (product_radius + reactant_structure.shape.radius) * MINIMAL_SEPARATION_FACTOR
                    for _ in range(self.dissociation_retry_moves):
                        unit_vector3D = random_unit_vector()
                        unit_vector2D = normalize(unit_vector3D - 
                                        (reactant_structure.shape.unit_z * numpy.dot(unit_vector3D, reactant_structure.shape.unit_z)))
                        vector = reactant_pos + vector_length * unit_vector2D
                        product_pos_list.append(vector)

                else:
                    # cannot decay from 3D to other structure
                    raise RuntimeError('fire_single_reaction: Can not decay from 3D to other structure')

                product_structure_id = self.world.get_def_structure_id()    # product structure is always the default structure (bulk)

            # If decay happens to same structure_type
            else:
                product_pos_list     = [reactant_pos]           # no change of position is required if structure_type doesn't change
                product_structure_id = reactant_structure_id    # The product structure is the structure where the reaction takes place.


            # 2. make space for the products (kinda brute force, but now we only have to burst once).
            if not ((product_pos_list[0] == reactant_pos).all()) or product_radius > reactant_radius:
                product_pos = self.world.apply_boundary(product_pos_list[0])
                radius = self.world.distance(product_pos, reactant_pos) + product_radius
                zero_singles, ignore = self.burst_volume(product_pos, radius, ignore)

            # 3. check that there is space for the products (try different positions if possible)
            # accept the new positions if there is enough space.
            for product_pos in product_pos_list:
                product_pos = self.world.apply_boundary(product_pos)

                # check that there is space for the products 
                # Note that we do not check overlap with surfaces (we assume that its no problem)
                if (not self.world.check_overlap((product_pos, product_radius), reactant[0])):

                    # 4. process the changes (remove particle, make new ones)
                    self.world.remove_particle(reactant[0])
                    newparticle = self.world.new_particle(product_species.id, product_structure_id, product_pos)
                    products = [newparticle]

                    # 5. update counters
                    self.reaction_events += 1
                    self.last_reaction = (rr, (reactant, None), products)

                    # 6. Log the change
                    if __debug__:
                        log.info('single reaction: product (%s) = %s' % (len(products), products))

                    # exit the loop, we have found a new position
                    break

            else:
                    # 4. Process (the lack of) change
                    moved_reactant = self.move_particle(reactant, reactant_pos, reactant_structure_id)
                    products = [moved_reactant]

                    # 5. update counters
                    self.rejected_moves += 1

                    # 6. Log the event
                    if __debug__:
                        log.info('single reaction: placing product failed.')

            
        elif len(rr.products) == 2:
            # 1. get product info
            product1_species            = self.world.get_species(rr.products[0])
            product2_species            = self.world.get_species(rr.products[1])
            product1_structure_type_id  = product1_species.structure_type_id 
            product2_structure_type_id  = product2_species.structure_type_id 
            product1_radius             = product1_species.radius
            product2_radius             = product2_species.radius
            D1 = product1_species.D
            D2 = product2_species.D
            particle_radius12           = product1_radius + product2_radius


            # 1.5 Get new positions and structure_ids of particles
            #     If one of the particle is not on the surface of the reactant.
            if product1_structure_type_id != product2_structure_type_id:
                # make sure that the reactant was on a 2D or 1D surface
                assert reactant_structure_type_id != self.world.get_def_structure_type_id()
                # make sure that only either one of the target structures is the 3D ('^' is exclusive-OR)
                assert ((product1_structure_type_id == self.world.get_def_structure_type_id()) ^ \
                        (product2_structure_type_id == self.world.get_def_structure_type_id()))

                # figure out which product stays in the surface and which one goes to 3D
                # Note that A is a particle in the surface and B is in the 3D
                if (product2_structure_type_id == self.world.get_def_structure_type_id()):
                    # product2 goes to 3D and is now particleB (product1 is particleA and is on the surface)
                    product1_structure_id = reactant_structure_id               # TODO after the displacement the structure can change!
                    product2_structure_id = self.world.get_def_structure_id()
                    productA_radius = product1_radius
                    productB_radius = product2_radius
                    DA = D1
                    DB = D2
                    default = True      # we like to think of this as the default
                else:
                    # product1 goes to 3D and is now particleB (product2 is particleA)
                    product2_structure_id = reactant_structure_id
                    product1_structure_id = self.world.get_def_structure_id()
                    productA_radius = product2_radius
                    productB_radius = product1_radius
                    DA = D2
                    DB = D1
                    default = False

                if isinstance(reactant_structure, PlanarSurface):
                    # draw a number of new positions for the two product particles
                    # TODO make this a generator

                    product_pos_list = []
                    for _ in range(self.dissociation_retry_moves):
                        # draw the random angle for the 3D particle relative to the particle left in the membrane

                        # do the backtransform with a random iv with length such that the particles are at contact
                        # Note we make the iv slightly longer because otherwise the anisotropic transform will produce illegal
                        # positions
                        iv = random_vector(particle_radius12 * MixedPair2D3D.calc_z_scaling_factor(DA, DB))
                        iv *= MINIMAL_SEPARATION_FACTOR

                        # determine the side of the membrane the dissociation takes place
                        if(reactant_structure.shape.is_one_sided):
                            unit_z = reactant_structure.shape.unit_z
                        else:
                            unit_z = reactant_structure.shape.unit_z * myrandom.choice(-1, 1)

                        # calculate the new positions and structure IDs
                        newposA, newposB, sidA, sidB = MixedPair2D3D.do_back_transform(reactant_pos, iv, DA, DB,
                                                                                       productA_radius, productB_radius,
                                                                                       reactant_structure, reactant_structure,
                                                                                       unit_z, self.world)
                        # the second reactant_structure parameter passed is ignored here

                        if default:
                            newpos1, newpos2 = newposA, newposB
                        else:
                            newpos1, newpos2 = newposB, newposA
                        product_pos_list.append((newpos1, newpos2))

                elif isinstance(reactant_structure, CylindricalSurface):

                    product_pos_list = []
                    for _ in range(self.dissociation_retry_moves):
                        iv = random_vector(particle_radius12 * MixedPair1D3D.calc_r_scaling_factor(DA, DB))
                        iv *= MINIMAL_SEPARATION_FACTOR

                        unit_z = reactant_structure.shape.unit_z

                        newposA, newposB, sidA, sidB = MixedPair1D3D.do_back_transform(reactant_pos, iv, DA, DB,
                                                                                       productA_radius, productB_radius,
                                                                                       reactant_structure, reactant_structure,
                                                                                       unit_z, self.world)
                        # the second reactant_structure parameter passed is ignored here
                        # unit_z is passed here for completeness but not used in the method

                        newposB = self.world.apply_boundary(newposB)
                        if __debug__:
                            assert (self.world.distance(reactant_structure.shape, newposB) >= productB_radius)

                        if default:
                            newpos1, newpos2 = newposA, newposB
                        else:
                            newpos1, newpos2 = newposB, newposA
                        product_pos_list.append((newpos1, newpos2))

                else:
                    # cannot decay from 3D to other structure
                    raise RuntimeError('fire_single_reaction: Can not decay from 3D to other structure')

            # 1.5 Get new positions and structure_ids of particles
            #     If the two particles stay on the same structure type as the reactant
            #     (but not necessarily on the same structure!)
            else:
                # Case product1_structure_type_id == product2_structure_type_id
                # Generate new positions in the structure
                # TODO Make this into a generator
                product_pos_list = []
                for _ in range(self.dissociation_retry_moves):
                    # FIXME for particles on the cylinder there are only two possibilities
                    # calculate a random vector in the structure with unit length
                    iv = _random_vector(reactant_structure, particle_radius12, self.rng)
                    iv *= MINIMAL_SEPARATION_FACTOR

                    unit_z = reactant_structure.shape.unit_z    # not used but passed for completeness

                    # If reactant_structure is a finite PlanarSurface the new position potentially can
                    # lie outside of the plane. The two product particles then will end up on different
                    # planes of the same structure_type. Therefore this case is treated separately.
                    # We use a specific do_back_transform method here which automatically takes care
                    # of the potential structure change and potentially resulting overlaps.
                    if isinstance(reactant_structure, PlanarSurface):

                        # Recalculate the back transform taking into account the deflection
                        newpos1, newpos2, sid1, sid2 = PlanarSurfaceTransitionPair.do_back_transform(reactant_pos, iv,
                                                                                                      D1, D2,
                                                                                                      product1_radius, product2_radius,
                                                                                                      reactant_structure, reactant_structure,
                                                                                                      unit_z, self.world)
                    else:
                        
                        # Calculate the standard back transform
                        newpos1, newpos2, sid1, sid2 = SimplePair.do_back_transform(reactant_pos, iv, D1, D2,
                                                                                    product1_radius, product2_radius,
                                                                                    reactant_structure, reactant_structure,
                                                                                    unit_z, self.world)
                        # for the SimplePair do_back_transform() requires that the two structures passed are the same
                        # because it will check for this!

                # Pass on the newly generated information
                product_pos_list.append((newpos1, newpos2))
                product1_structure_id = sid1
                product2_structure_id = sid2


            # 2. make space for the products. 
            #    calculate the sphere around the two product particles
            product1_pos = self.world.apply_boundary(product_pos_list[0][0])
            product2_pos = self.world.apply_boundary(product_pos_list[0][1])
            radius = max(self.world.distance(product1_pos, reactant_pos) + product1_radius,
                         self.world.distance(product2_pos, reactant_pos) + product2_radius)
            zero_singles, ignore = self.burst_volume(reactant_pos, radius, ignore)

            # 3. check that there is space for the products (try different positions if possible)
            # accept the new positions if there is enough space.
            for newpos1, newpos2 in product_pos_list:
                newpos1 = self.world.apply_boundary(newpos1)
                newpos2 = self.world.apply_boundary(newpos2)

                if __debug__:
                    assert (self.world.distance(newpos1, newpos2) >= particle_radius12)

                # check that there is space for the products 
                # Note that we do not check overlap with surfaces -> TODO
                if (not self.world.check_overlap((newpos1, product1_radius), reactant[0])
                   and
                   not self.world.check_overlap((newpos2, product2_radius), reactant[0])):
                    
                    # 4. process the changes (remove particle, make new ones)
                    self.world.remove_particle(reactant[0])
                    product1_pair = self.world.new_particle(product1_species.id, product1_structure_id, newpos1)
                    product2_pair = self.world.new_particle(product2_species.id, product2_structure_id, newpos2)
                    products = [product1_pair, product2_pair]

                    # 5. update counters
                    self.reaction_events += 1
                    self.last_reaction = (rr, (reactant, None), products)

                    # 6. Log the change
                    if __debug__:
                        log.info('single reaction: products (%s) = %s' % (len(products), products))

                    # exit the loop, we have found new positions
                    break
            else:
                # Log the lack of change
                if __debug__:
                    log.info('single reaction: placing products failed.')

                # 4. process the changes (move the particle and update structure)
                moved_reactant = self.move_particle(reactant, reactant_pos, reactant_structure_id)
                products = [moved_reactant]

                # 5. update counters.
                self.rejected_moves += 1

        else:
            raise RuntimeError('fire_single_reaction: num products >= 3 not supported.')

        return products, zero_singles, ignore


    def fire_interaction(self, single, reactant_pos, reactant_structure_id, ignore):
        # This takes care of the identity change when a particle associates with a surface.
        # It performs the reactions:
        # A(structure) + surface -> 0
        # A(structure) + surface -> B(surface)

        # Note that the reactant_structure_id is the id of the structure on which the particle was located at the time of the reaction.

        if __debug__:
            assert isinstance(single, InteractionSingle)

        # 0. get reactant info
        reactant              = single.pid_particle_pair
        reactant_radius       = reactant[1].radius
        rr = single.interactionrule

        # The zero_singles are the NonInteractionSingles that are the total result of the recursive
        # bursting process.
        zero_singles = []

        if len(rr.products) == 0:
            
            # 1. No products, no info
            # 1.5 No new position/structure_id required
            # 2. No space required
            # 3. No space required
            # 4. process the changes (remove particle, make new ones)
            self.world.remove_particle(reactant[0])
            products = []
            # zero_singles/ignore is unchanged.

            # 5. Update counters
            self.reaction_events += 1
            self.last_reaction = (rr, (reactant, None), products)

            # 6. Log the change
            if __debug__:
                 log.info('interaction: product = None.')
            
        elif len(rr.products) == 1:
            
            # 1. get product info
            product_species      = self.world.get_species(rr.products[0])
            product_radius       = product_species.radius
            product_structure    = single.target_structure

            # TODO make sure that the product species lives on the same type of the interaction surface
           #product_structure_type = self.world.get_structure(product_species.structure_id)
           #assert (single.target_structure.sid == self.world.get_structure(product_species.structure_id)), \
                  #'Product particle should live on the surface of interaction after the reaction.'

            # 1.5 get new position and structure_id of particle
            transposed_pos       = self.world.cyclic_transpose(reactant_pos, product_structure.shape.position)
            product_pos, _       = product_structure.project_point(transposed_pos)
            product_pos          = self.world.apply_boundary(product_pos)        # not sure this is still necessary
            product_structure_id = product_structure.id

            # 2. burst the volume that will contain the products.
            #    Note that an interaction domain is already sized such that it includes the
            #    reactant particle moving into the surface.
            #    Note that the burst is recursive!!
            if product_radius > reactant_radius:
                zero_singles, ignore = self.burst_volume(product_pos, product_radius, ignore)

            # 3. check that there is space for the products 
            # Note that we do not check for interfering surfaces (we assume this is no problem)
            if self.world.check_overlap((product_pos, product_radius), reactant[0]):

                # Log the (absence of) change
                if __debug__:
                    log.info('interaction: no space for product particle.')

                # 4. process the changes (only move particle, stay on same structure)
                moved_reactant = self.move_particle(reactant, reactant_pos, reactant_structure_id)
                products = [moved_reactant]

                # 5. update counters
                self.rejected_moves += 1

            else:
                # 4. process the changes (remove particle, make new ones)
                self.world.remove_particle(reactant[0])
                newparticle = self.world.new_particle(product_species.id, product_structure_id, product_pos)
                products = [newparticle]

                # 5. update counters
                self.reaction_events += 1
                self.last_reaction = (rr, (reactant, None), products)

                # 6. Log the change
                if __debug__:
                     log.info('interaction: product (%s) = %s' % (len(products), products))

        else:
            raise RuntimeError('fire_interaction: num products > 1 not supported.')

        return products, zero_singles, ignore


    def fire_pair_reaction(self, pair, reactant1_pos, reactant2_pos, reactant1_structure_id, reactant2_structure_id, ignore):
        # This takes care of the identity change when two particles react with each other
        # It performs the reactions:
        # A(any structure) + B(same structure) -> 0
        # A(any structure) + B(same structure or 3D) -> C(same structure)
        if __debug__:
            assert isinstance(pair, Pair)

        # Note that the reactant_structure_ids are the ids of the structures on which the particles were located at the time of the reaction.

        # 0. get reactant info
        pid_particle_pair1     = pair.pid_particle_pair1
        pid_particle_pair2     = pair.pid_particle_pair2
        rr = pair.draw_reaction_rule(pair.interparticle_rrs)

        # The zero_singles are the NonInteractionSingles that are the Total result of the recursive
        # bursting process.
        zero_singles = []

        if len(rr.products) == 0:

            # 1. no product particles, no info
            # 1.5 No new position/structure_id required
            # 2. No space required
            # 3. No space required
            # 4. process the change (remove particles, make new ones)
            self.world.remove_particle(pid_particle_pair1[0])
            self.world.remove_particle(pid_particle_pair2[0])
            products = []
            # zero_singles/ignore is unchanged.

            # 5. update counters
            self.reaction_events += 1
            self.last_reaction = (rr, (pid_particle_pair1, pid_particle_pair2), products)

            # 6. Log the change
            if __debug__:
                 log.info('pair reaction: product = None.')

        elif len(rr.products) == 1:

            # 1. get product info
            product_species = self.world.get_species(rr.products[0])
            product_radius  = product_species.radius


            # 1.5 get new position and structure_id for product particle
            product_pos = pair.draw_new_com (pair.dt, pair.event_type)
            product_pos = self.world.apply_boundary(product_pos)

            # select the structure_id for the product particle to live on.
            # TODO doesn't work for all needed cases
            product_structure_ids = self.world.get_structure_ids(self.world.get_structure_type(product_species.structure_type_id))
            if reactant1_structure_id in product_structure_ids:
                product_structure_id = reactant1_structure_id
            elif reactant2_structure_id in product_structure_ids:
                product_structure_id = reactant2_structure_id

            # 2.  make space for products
            # 2.1 if the product particle sticks out of the shell

            #if self.world.distance(old_com, product_pos) > (pair.get_shell_size() - product_radius): # TODO What is that?
            zero_singles, ignore = self.burst_volume(product_pos, product_radius, ignore)

            # 3. check that there is space for the products ignoring the reactants
            # Note that we do not check for interfering surfaces (we assume this is no problem)
            if self.world.check_overlap((product_pos, product_radius), pid_particle_pair1[0],
                                                                       pid_particle_pair2[0]):
                # Log the (lack of) change
                if __debug__:
                    log.info('pair reaction: no space for product particle.')

                # 4. process the changes (move particles, change structure)
                moved_reactant1 = self.move_particle(pid_particle_pair1, reactant1_pos, reactant1_structure_id)
                moved_reactant2 = self.move_particle(pid_particle_pair2, reactant2_pos, reactant2_structure_id)
                products = [moved_reactant1, moved_reactant2]

                # 5. update counters
                self.rejected_moves += 1

            else:
                # 4. process the changes (remove particle, make new ones)
                self.world.remove_particle(pid_particle_pair1[0])
                self.world.remove_particle(pid_particle_pair2[0])
                newparticle = self.world.new_particle(product_species.id, product_structure_id, product_pos)
                products = [newparticle]

                # 5. update counters
                self.reaction_events += 1
                self.last_reaction = (rr, (pid_particle_pair1, pid_particle_pair2),
                                      products)

                # 6. Log the change
                if __debug__:
                    log.info('pair reaction: product (%s) = %s' % (len(products), products))

        else:
            raise NotImplementedError('fire_pair_reaction: num products >= 2 not supported.')

        return products, zero_singles, ignore


    def fire_move(self, single, reactant_pos, reactant_structure_id, ignore_p=None):
        # No reactions/Interactions have taken place -> no identity change of the particle
        # Just only move the particles and process putative structure change

        # Note that the reactant_structure_id is the id of the structure on which the particle was located at the time of the move.

        if __debug__:
            assert isinstance(single, Single)

        # 0. get reactant info
        reactant = single.pid_particle_pair
        # 1. No product->no info
        # 1.5 No additional positional/structure change is needed
        # 2. Position is in protective shell

        # 3. check that there is space for the reactants
        if ignore_p:
            if self.world.check_overlap((reactant_pos, reactant[1].radius),
                                        reactant[0], ignore_p[0]):
                raise RuntimeError('fire_move: particle overlap failed.')
        else:
            if self.world.check_overlap((reactant_pos, reactant[1].radius),
                                        reactant[0]):
                raise RuntimeError('fire_move: particle overlap failed.')

        # 4. process the changes (move particles, change structure)
        moved_reactant = self.move_particle(reactant, reactant_pos, reactant_structure_id)
        # 5. No counting
        # 6. No Logging
        return [moved_reactant]


    #########################################
    #### METHODS FOR NEW DOMAIN CREATION ####
    #########################################
    def make_new_domain(self, single):
        ### Make a new domain out of a NonInteractionSingle that was
        ### just fired

        # can only make new domain from NonInteractionSingle
        assert isinstance(single, NonInteractionSingle)
        assert single.is_reset()
        
        # Special case: When D=0, nothing needs to change and only new event
        # time is drawn
        # TODO find nicer construction than just this if statement
       #if single.getD() == 0:
           #single.dt, single.event_type = single.determine_next_event() 
           #single.last_time = self.t
           #self.add_domain_event(single)
           #return single


        # 0. Get generic info
        single_pos = single.pid_particle_pair[1].position
        single_radius = single.pid_particle_pair[1].radius
        reaction_threshold = single_radius * SINGLE_SHELL_FACTOR

        # 1.0 Get neighboring domains and surfaces
        neighbor_distances = self.geometrycontainer.get_neighbor_domains(single_pos, self.domains, ignore=[single.domain_id, ])
        # Get also surfaces but only if the particle is in 3D
        surface_distances = get_neighbor_surfaces(self.world, single_pos, single.structure.id, ignores=[])

        # 2.3 We prefer to make NonInteractionSingles for efficiency.
        #     But if objects (Shells and surfaces) get close enough (closer than
        #     reaction_threshold) we want to try Interaction and/or Pairs. 

        #     From the list of neighbors, we calculate the thresholds (horizons) for trying to form
        #     the various domains. The domains can be:
        #     -zero-shell NonInteractionSingle -> Pair or Multi
        #     -Surface -> Interaction
        #     -Multi -> Multi
        #
        # Note that we do not differentiate between directions. This means that we
        # look around in a sphere, the horizons are spherical.
        pair_interaction_partners = []
        for domain, _ in neighbor_distances:
            if (isinstance (domain, NonInteractionSingle) and domain.is_reset()):
                # distance from the center of the particles/domains
                pair_distance = self.world.distance(single_pos, domain.shell.shape.position)
                pair_horizon  = (single_radius + domain.pid_particle_pair[1].radius) * SINGLE_SHELL_FACTOR
                pair_interaction_partners.append((domain, pair_distance - pair_horizon))
        
        for surface, surface_distance in surface_distances:            
            if isinstance(surface, PlanarSurface):# or isinstance(surface, DiskSurface):
                # with a planar surface it is the center of mass that 'looks around'
                surface_horizon = single_radius * (SINGLE_SHELL_FACTOR - 1.0)            
            else:
                # with a cylindrical surface it is the surface of the particle
                surface_horizon = single_radius * SINGLE_SHELL_FACTOR

            pair_interaction_partners.append((surface, surface_distance - surface_horizon))

        pair_interaction_partners = sorted(pair_interaction_partners, key=lambda domain_overlap: domain_overlap[1])

        # For each particles/particle or particle/surface pair, check if a pair or interaction can be
        # made. Note that we check the closest one first and only check if the object are within the horizon
        domain = None
        for obj, hor_overlap in pair_interaction_partners:

            if hor_overlap > 0.0 or domain or self.BD_ONLY_FLAG :
                # there are no more potential partners (everything is out of range)
                # or a domain was formed successfully previously
                # or the user wants to force the system into BD mode = Multi creation
                break

            if isinstance(obj, NonInteractionSingle):
                # try making a Pair (can be Mixed Pair or Normal Pair)
                domain = self.try_pair (single, obj)

            elif isinstance(obj, PlanarSurface) and isinstance(single.structure, PlanarSurface):
                domain = self.try_transition(single, obj)

            elif ( isinstance(obj, CylindricalSurface) or isinstance(obj, PlanarSurface) or isinstance(obj, DiskSurface) ):
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
                    multi_horizon = (single_radius + domain.pid_particle_pair[1].radius) * MULTI_SHELL_FACTOR
                    distance = self.world.distance(single_pos, domain.shell.shape.position)
                    multi_partners.append((domain, distance - multi_horizon))

                elif isinstance(domain, Multi):
                    # The dist_to_shell = dist_to_particle - multi_horizon_of_target_particle
                    # So only the horizon and distance of the current single needs to be taken into account
                    # Note: this is built on the assumption that the shell of a Multi has the size of the horizon.
                    multi_horizon = (single_radius * MULTI_SHELL_FACTOR)
                    multi_partners.append((domain, dist_to_shell - multi_horizon))

            # Also add surfaces
            for surface, distance in surface_distances:
                if isinstance(surface, PlanarSurface):
                    # with a planar surface it is the center of mass that 'looks around'
                    surface_horizon = single_radius * (MULTI_SHELL_FACTOR - 1.0)
                else:
                    # with a cylindrical surface it is the surface of the particle
                    surface_horizon = single_radius * MULTI_SHELL_FACTOR

                multi_partners.append((surface, distance - surface_horizon))


            # get the closest object (if there)
            if multi_partners:
                multi_partners = sorted(multi_partners, key=lambda domain_overlap: domain_overlap[1])
                #log.debug('multi_partners: %s' % str(multi_partners))
                closest_overlap = multi_partners[0][1]
            else:
                # In case there is really nothing
                closest_overlap = numpy.inf

            if __debug__:
                log.debug('Single or Multi: closest_overlap: %s' % (FORMAT_DOUBLE % closest_overlap))

            # If the closest partner is within the multi horizon we do Multi, otherwise Single
            if closest_overlap > 0.0 and not self.BD_ONLY_FLAG : 
                # just make a normal NonInteractionSingle
                self.update_single(single)
            else:
                # An object was closer than the Multi horizon
                # Form a multi with everything that is in the multi_horizon
                domain = self.form_multi(single, multi_partners)

        return domain


    #### Redundant but maybe useful for parts
    def calculate_simplepair_shell_size(self, single1, single2):
        # 3. Calculate the maximum based on some other criterium (convergence?)
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius
        sigma = radius1 + radius2
        distance_from_sigma = r0 - sigma        # the distance between the surfaces of the particles
       #convergence_max = distance_from_sigma * 100 + sigma + shell_size_margin                # FIXME
      #max_shell_size = min(max_shell_size, convergence_max)

        #5. Calculate the 'ideal' shell size, it matches the expected first-passage time of the closest particle
        #with the expected time of the particles in the pair.
       #D1 = single1.pid_particle_pair[1].D
       #D2 = single2.pid_particle_pair[1].D
       #D12 = D1 + D2

       #D_closest = closest.pid_particle_pair[1].D
       #D_tot = D_closest + D12
       #ideal_shell_size = min((D12 / D_tot) *
                           #(closest_particle_distance - min_shell_size - closest_min_radius) +
                           #min_shell_size,
        ideal_shell_size = numpy.inf
        shell_size = min(max_shell_size, ideal_shell_size)


        # 6. Check if singles would not be better.
        # TODO clear this up -> strange rule
        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        d1 = self.world.distance(com, pos1)
        d2 = self.world.distance(com, pos2)

        if shell_size < max(d1 + single1.pid_particle_pair[1].radius *
                            SINGLE_SHELL_FACTOR, \
                            d2 + single2.pid_particle_pair[1].radius * \
                            SINGLE_SHELL_FACTOR) * 1.3:
            if __debug__:
                log.debug('%s not formed: singles are better' %
                          'Pair(%s, %s)' % (single1.pid_particle_pair[0], 
                                            single2.pid_particle_pair[0]))
            return None, None, None


    #################################################
    #### METHODS FOR EVENT CONSEQUENCES HANDLING ####
    #################################################
    def update_single(self, single): 
    # updates a NonInteractionSingle given that the single already holds its new position

        # The method only works for NonInteractionSingles
        assert isinstance(single, NonInteractionSingle)

        singlepos = single.pid_particle_pair[1].position    # TODO get the position as an argument


        # create a new updated shell
        new_shell = single.create_updated_shell(singlepos)
        assert new_shell, 'single.create_updated_shell() returned None.'

        # Replace shell in domain and geometrycontainer.
        # Note: this should be done before determine_next_event.
        # Reuse shell_id.
        shell_id_shell_pair = (single.shell_id, new_shell) 
        single.shell_id_shell_pair = shell_id_shell_pair
        self.geometrycontainer.move_shell(shell_id_shell_pair)
        if __debug__:
            log.info('        *Updated single: single: %s, new_shell = %s' % \
                     (single, str(new_shell)))

        single.dt, single.event_type = single.determine_next_event()
        assert single.dt >= 0
        if __debug__:
            log.info('dt=%s' % (FORMAT_DOUBLE % single.dt)) 

        single.last_time = self.t

        # add to scheduler
        self.add_domain_event(single)
        # check everything is ok
        if __debug__:
            assert self.check_domain(single)

        return single


    def process_single_event(self, single, ignore):
    # This method handles the things that need to be done when the current event was
    # produced by a single. The single can be a NonInteractionSingle or an InteractionSingle.
    # Note that this method is also called when a single is bursted, in that case the event
    # is a BURST event.
    #
    # Note that 'ignore' should already contain the id of 'single'.
    #
    # This method returns:
    #
    # - the 'domains' that were created as a result of the processing. This can be a newly made domains
    #   (since make_new_domain is also called from here) or zero_singles
    #   Note that these could be many since the event can induce recursive bursting. 
    # - the updated ignore list. This contains the id of the 'single' domain, IDs of domains
    #   bursted in the process and the IDs of zero-dt singles.

        if __debug__:
           ## TODO assert that there is no event associated with this domain in the scheduler
           #assert self.check_domain(single)
            assert (single.domain_id in ignore), \
                   'Domain_id should already be on ignore list before processing event.'

        ### SPECIAL CASE 1:
        # If just need to make new domain.
        if single.is_reset():
            if __debug__:
                log.info('FIRE SINGLE: make_new_domain()')
                log.info('single = %s' % single)
            domains = [self.make_new_domain(single)]
            # domain is already scheduled in make_new_domain

        ### SPECIAL CASE 2:
        # In case nothing is scheduled to happen: do nothing and just reschedule
        elif single.dt == numpy.inf:
            if __debug__:
                log.info('FIRE SINGLE: Nothing happens-> single.dt=inf')
                log.info('single = %s' % single)
            self.add_domain_event(single)
            domains = [single]

        ### NORMAL:
        # Process 'normal' event produced by the single
        else:
            if __debug__:
                log.info('FIRE SINGLE: %s' % single.event_type)
                log.info('single = %s' % single)

            ### 1. check that everything is ok
            if __debug__:
                if single.event_type != EventType.BURST:
                    # The burst of the domain may be caused by an overlapping domain
                    assert self.check_domain(single)

                # check that the event time of the single (last_time + dt) is equal to the
                # simulator time
                assert (abs(single.last_time + single.dt - self.t) <= TIME_TOLERANCE * self.t), \
                       'Timeline incorrect. single.last_time = %s, single.dt = %s, self.t = %s' % \
                        (FORMAT_DOUBLE % single.last_time, FORMAT_DOUBLE % single.dt, FORMAT_DOUBLE % self.t)


            ### 2. get some necessary information
            pid_particle_pair = single.pid_particle_pair
            # in case of an interaction domain: determine real event before doing anything
            if single.event_type == EventType.IV_EVENT:
                single.event_type = single.draw_iv_event_type()


            # get the (new) position and structure on which the particle is located.
            if single.getD() != 0 and single.dt > 0.0:
                # If the particle had the possibility to diffuse
                newpos, struct_id = single.draw_new_position(single.dt, single.event_type)
                newpos, struct_id = self.world.apply_boundary((newpos,struct_id))
            else:
                # no change in position has taken place
                newpos    = pid_particle_pair[1].position
                struct_id = pid_particle_pair[1].structure_id
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
                    self.single_steps[single.event_type] += 1
                    particles, zero_singles_b, ignore = self.fire_single_reaction(single, newpos, struct_id, ignore)
                    # The 'zero_singles_b' list now contains the resulting singles of a putative burst
                    # that occured in fire_single_reaction

                else:
                    self.interaction_steps[single.event_type] += 1
                    particles, zero_singles_b, ignore = self.fire_interaction(single, newpos, struct_id, ignore)

            else:
                particles = self.fire_move(single, newpos, struct_id)
                zero_singles_b = []   # no bursting takes place, ignore list remains unchanged
                # Update statistics
                if single.event_type == EventType.SINGLE_ESCAPE or\
                   single.event_type == EventType.BURST:
                    self.single_steps[single.event_type] += 1                
                elif __debug__:
                    log.warning('Omitting to count a single event with unforeseen event type (%s).' % single.event_type)
                    

            domains = zero_singles_b

            ### 5. Make a new domain (or reuse the old one) for each particle(s)
            #      (Re)schedule the (new) domain
            zero_singles = []
            for pid_particle_pair in particles:
#               single.pid_particle_pair = pid_particle_pair         # TODO reuse single
                zero_single = self.create_single(pid_particle_pair)  # TODO re-use NonInteractionSingle domain if possible
                self.add_domain_event(zero_single)                   # TODO re-use event if possible
                zero_singles.append(zero_single)       # Add the newly made zero-dt singles to the list
                ignore.append(zero_single.domain_id)   # Ignore these newly made singles (they should not be bursted)
            domains.extend(zero_singles)

            ### 6. Recursively burst around the newly made zero-dt NonInteractionSingles that surround the particles.
            #      For each zero_single the burst radius equals the particle radius * SINGLE_SHELL_FACTOR
            #      Add the resulting zero_singles to the already existing list.
            for zero_single in zero_singles:
                more_zero_singles, ignore = self.burst_non_multis(zero_single.pid_particle_pair[1].position,
                                                                  zero_single.pid_particle_pair[1].radius*SINGLE_SHELL_FACTOR,
                                                                  ignore)
                domains.extend(more_zero_singles)

            if __debug__:
                # check that at least all the zero_singles are on the ignore list
                assert all(zero_single.domain_id in ignore for zero_single in domains)

            ### 7. Log change?
                
       #NOTE Code snippets to re-use the domain/event
       #if isinstance(single, InteractionSingle):
           ## When bursting an InteractionSingle, the domain changes from Interaction
          ## NonInteraction domain. This needs to be reflected in the event 
           #self.remove_event(single)
           #self.add_domain_event(newsingle)
       #else:
           #assert single == newsingle
           #self.update_domain_event(self.t, single)

       #particle_radius = single.pid_particle_pair[1].radius
       #assert newsingle.shell.shape.radius == particle_radius

        return domains, ignore


    def process_pair_event(self, pair, ignore):
    # This method handles the things that need to be done when the current event was
    # produced by a pair. The pair can be any type of pair (Simple or Mixed).
    # Note that this method is also called when a pair is bursted, in that case the event
    # is a BURST event.
    # Note that ignore should already contain the id of 'pair'.
    # Returns:
    # - the 'domains' that were the results of the processing. These are zero_singles that are
    #   the result of the process. Note that these could be many since the event can induce
    #   recursive bursting.
    # - the updated ignore list. This contains the id of the 'pair' domain, ids of domains
    #   bursted in the process and the ids of zero-dt singles.

        if __debug__:
            log.info('FIRE PAIR: %s' % pair.event_type)
            log.info('single1 = %s' % pair.single1)
            log.info('single2 = %s' % pair.single2)

        ### 1. check that everything is ok
        if __debug__:
            assert (pair.domain_id in ignore), \
                   'Domain_id should already be on ignore list before processing event.'
            assert self.check_domain(pair)
            assert pair.single1.domain_id not in self.domains
            assert pair.single2.domain_id not in self.domains
            # TODO assert that there is no event associated with this domain in the scheduler

            # check that the event time of the pair (last_time + dt) is equal to the
            # simulator time
            assert (abs(pair.last_time + pair.dt - self.t) <= TIME_TOLERANCE * self.t), \
                    'Timeline incorrect. pair.last_time = %s, pair.dt = %s, self.t = %s' % \
                    (FORMAT_DOUBLE % pair.last_time, FORMAT_DOUBLE % pair.dt, FORMAT_DOUBLE % self.t)

        ### 2. get some necessary information
        single1 = pair.single1
        single2 = pair.single2
        pid_particle_pair1 = pair.pid_particle_pair1
        pid_particle_pair2 = pair.pid_particle_pair2
        oldpos1 = pid_particle_pair1[1].position
        oldpos2 = pid_particle_pair2[1].position


        if pair.event_type == EventType.IV_EVENT:
            # Draw actual pair event for iv at very last minute.
            pair.event_type = pair.draw_iv_event_type(pair.r0)


        ### 3. Process the event produced by the pair

        # Count the event in the event statistics
        self.pair_steps[pair.event_type] += 1

        ### 3.1 Get new position and current structures of particles
        if pair.dt > 0.0:
            newpos1, newpos2, struct1_id, struct2_id = pair.draw_new_positions(pair.dt, pair.r0, pair.iv, pair.event_type)
            newpos1, struct1_id = self.world.apply_boundary((newpos1, struct1_id))
            newpos2, struct2_id = self.world.apply_boundary((newpos2, struct2_id))

            # check that the particles do not overlap with any other particles in the world
            assert not self.world.check_overlap((newpos1, pid_particle_pair1[1].radius),
                                                pid_particle_pair1[0], pid_particle_pair2[0])
            assert not self.world.check_overlap((newpos2, pid_particle_pair2[1].radius),
                                                pid_particle_pair1[0], pid_particle_pair2[0])

        else:
            newpos1 = pid_particle_pair1[1].position
            newpos2 = pid_particle_pair2[1].position
            struct1_id = pid_particle_pair1[1].structure_id
            struct2_id = pid_particle_pair2[1].structure_id
        # TODO? Check if the new positions are within domain
        # TODO? some more consistency checking of the positions
        #     assert self.check_pair_pos(pair, newpos1, newpos2, old_com,
        #                                pair.get_shell_size())


        # newpos1/2 now hold the new positions of the particles (not yet committed to the world)
        # if we would here move the particles and make new shells, then it would be similar to a propagate

        self.remove_domain(pair)

        ### 3.2 If identity changing processes have taken place
        #       Four cases:
        #       1. Single reaction
        if pair.event_type == EventType.SINGLE_REACTION:
            reactingsingle = pair.reactingsingle

            if reactingsingle == single1:
                theothersingle = single2
                # first move the non-reacting particle
                particles = self.fire_move(single2, newpos2, struct2_id, pid_particle_pair1)
                # don't ignore the (moved) non-reacting particle
                more_particles, zero_singles_b, ignore = self.fire_single_reaction(single1, newpos1, struct1_id, ignore)
                particles.extend(more_particles)
            else:
                theothersingle = single1
                particles = self.fire_move(single1, newpos1, struct1_id, pid_particle_pair2)
                more_particles, zero_singles_b, ignore = self.fire_single_reaction(single2, newpos2, struct2_id, ignore)
                particles.extend(more_particles)

            if __debug__:
                log.info('reactant = %s' % reactingsingle)

            # Make new NonInteractionSingle domains for every particle after the reaction.
            zero_singles = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                zero_single = self.create_single(pid_particle_pair)  # TODO reuse the non-reacting single domain
                self.add_domain_event(zero_single)
                zero_singles.append(zero_single)

        #
        #       2. Pair reaction
        #
        elif pair.event_type == EventType.IV_REACTION:

            particles, zero_singles_b, ignore = self.fire_pair_reaction (pair, newpos1, newpos2, struct1_id, struct2_id, ignore)

            # Make new NonInteractionSingle domains for every particle after the reaction.
            zero_singles = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                zero_single = self.create_single(pid_particle_pair)
                self.add_domain_event(zero_single)
                zero_singles.append(zero_single)

        # Just moving the particles
        #       3a. IV escape
        #       3b. CoM escape
        elif(pair.event_type == EventType.IV_ESCAPE or
             pair.event_type == EventType.COM_ESCAPE or
             pair.event_type == EventType.BURST):

            particles      = self.fire_move (single1, newpos1, struct1_id, pid_particle_pair2)
            particles.extend(self.fire_move (single2, newpos2, struct2_id, pid_particle_pair1))
            zero_singles_b = []
            # ignore is unchanged

            # Make new NonInteractionSingle domains for every particle after the reaction.
            zero_singles = []
            for pid_particle_pair in particles:
                # 5. make a new single and schedule
                zero_single = self.create_single(pid_particle_pair)  # TODO reuse domains that were cached in the pair
                self.add_domain_event(zero_single)
                zero_singles.append(zero_single)

        else:
            raise SystemError('process_pair_event: invalid event_type.')

       ##NOTE code snippit to reuse domains that were cached in pair domain
       ## re-use single domains for particles
       ## normally we would take the particles + new positions from the old pair
       ## and use them to make new singles. Now we re-use the singles stored in the
       ## pair
       #single1 = pair.single1
       #single2 = pair.single2
       #assert single1.domain_id not in self.domains
       #assert single2.domain_id not in self.domains
       #single1.pid_particle_pair = pid_particle_pair1  # this is probably redundant
       #single2.pid_particle_pair = pid_particle_pair2

       ## 'make' the singles
       #single1.initialize(self.t)
       #single2.initialize(self.t)
       
       #self.update_single_shell(single1, newpos1, pid_particle_pair1[1].radius)
       #self.update_single_shell(single2, newpos2, pid_particle_pair2[1].radius)


       #self.domains[single1.domain_id] = single1
       #self.domains[single2.domain_id] = single2

       ## Check the dimensions of the shells of the singles with the shell in the container
       #if __debug__:
           #container1 = self.geometrycontainer.get_container(single1.shell)
           #assert container1[single1.shell_id].shape.radius == \
                  #single1.shell.shape.radius
           #if type(single1.shell) is CylindricalShell:
               #assert container1[single1.shell_id].shape.half_length == \
                      #single1.shell.shape.half_length

           #container2 = self.geometrycontainer.get_container(single2.shell)
           #assert container2[single2.shell_id].shape.radius == \
                  #single2.shell.shape.radius
           #if type(single2.shell) is CylindricalShell:
               #assert container2[single2.shell_id].shape.half_length == \
                      #single2.shell.shape.half_length

       #assert single1.shell.shape.radius == pid_particle_pair1[1].radius
       #assert single2.shell.shape.radius == pid_particle_pair2[1].radius
       ## even more checking
       #assert self.check_domain(single1)
       #assert self.check_domain(single2)
       ## Now finally we are convinced that the singles were actually made correctly


        # Put the zero singles resulting from putative bursting in fire_... in the final list.
        zero_singles_fin = zero_singles_b
        # Ignore the zero singles that were made around the particles (cannot be bursted)
        for zero_single in zero_singles:
            ignore.append(zero_single.domain_id)
        zero_singles_fin.extend(zero_singles)       # Add the zero-dt singles around the particles

        ### 6. Recursively burst around the newly made zero-dt NonInteractionSingles that surround the particles.
        for zero_single in zero_singles:
            more_zero_singles, ignore = self.burst_non_multis(zero_single.pid_particle_pair[1].position,
                                                              zero_single.pid_particle_pair[1].radius*SINGLE_SHELL_FACTOR,
                                                              ignore)
            # Also add the zero-dt singles from this bursting to the final list
            zero_singles_fin.extend(more_zero_singles)

        # Some usefull checks in debugging mode
        # Check that at least all the zero_singles are on the ignore list
        if __debug__:
            assert all(zero_single.domain_id in ignore for zero_single in zero_singles_fin)

        ### 7. Log the event
        if __debug__:
            log.debug("process_pair_event: #1 { %s: %s => %s }" %
                      (single1, str(oldpos1), str(newpos1)))
            log.debug("process_pair_event: #2 { %s: %s => %s }" %
                      (single2, str(oldpos2), str(newpos2)))

        return zero_singles_fin, ignore


    def process_multi_event(self, multi, ignore):
    # This method handles the things that need to be done when the current event was
    # produced by a multi.
    # Returns:
    # - Nothing if the particles were just propagated.
    # Or:
    # - The zero_singles that are the result when the multi was broken up due to an event.
    #   Note that these could be many since the breakup can induce recursive bursting.
    # - The updated ignore list in case of an event. This contains the id of the 'multi'
    #   domain, ids of domains bursted in the process and the ids of zero-dt singles of the
    #   particles.

        if __debug__:
            log.info('FIRE MULTI: %s' % multi.last_event)

        if __debug__:
            assert (multi.domain_id in ignore), \
                   'Domain_id should already be on ignore list before processing event.'

            # check that the event time of the multi (last_time + dt) is equal to the
            # simulator time
            assert (abs(multi.last_time + multi.dt - self.t) <= TIME_TOLERANCE * self.t), \
                    'Timeline incorrect. multi.last_time = %s, multi.dt = %s, self.t = %s' % \
                    (FORMAT_DOUBLE % multi.last_time, FORMAT_DOUBLE % multi.dt, FORMAT_DOUBLE % self.t)

        self.multi_steps[multi.last_event] += 1
        self.multi_steps[3] += 1  # multi_steps[3]: total multi steps
        multi.step()

        if(multi.last_event == EventType.MULTI_UNIMOLECULAR_REACTION or
           multi.last_event == EventType.MULTI_BIMOLECULAR_REACTION):
            self.reaction_events += 1
            self.last_reaction = multi.last_reaction

        if multi.last_event is not EventType.MULTI_DIFFUSION:                # if an event took place
            zero_singles, ignore = self.break_up_multi(multi, ignore)
            #self.multi_steps[multi.last_event] += 1
        else:
            multi.last_time = self.t
            self.add_domain_event(multi)
            zero_singles = []
            # ignore is unchanged

        return zero_singles, ignore


    def break_up_multi(self, multi, ignore):
        # Dissolves 'multi' into zero_singles, single with a zero shell (dt=0)
        # - 'ignore' contains domain ids that should be ignored
        # returns:
        # - updated ignore list
        # - zero_singles that were the product of the breakup

        self.remove_domain(multi)

        zero_singles = []
        for pid_particle_pair in multi.particles:
            zero_single = self.create_single(pid_particle_pair)
            self.add_domain_event(zero_single)
            zero_singles.append(zero_single)
            ignore.append(zero_single.domain_id)

        zero_singles_fin = zero_singles     # Put the zero-dt singles around the particles into the final list

        ### Recursively burst around the newly made zero-dt NonInteractionSingles that surround the particles.
        for zero_single in zero_singles:
            zero_singles2, ignore = self.burst_non_multis(zero_single.pid_particle_pair[1].position,
                                                          zero_single.pid_particle_pair[1].radius*SINGLE_SHELL_FACTOR,
                                                          ignore)
            zero_singles_fin.extend(zero_singles2)

        # check that at least all the zero_singles are on the ignore list
        if __debug__:
            assert all(zero_single.domain_id in ignore for zero_single in zero_singles_fin)

        return zero_singles_fin, ignore


    ###################################################
    #### METHODS FOR NEW SHELL CONSTRUCTION TRIALS ####
    ###################################################
    def try_interaction(self, single, target_structure):
    # Try to form an interaction between the 'single' particle and the 'target_structure'. The interaction can be a:
    # -CylindricalSurfaceInteraction
    # -PlanarSurfaceInteraction

        if __debug__:
           log.debug('trying to form Interaction(%s, %s)' % (single.pid_particle_pair[1], target_structure))


        ### 1. Attempt to make a testShell. If it fails it will throw an exception.
        try:
            testShell = try_default_testinteraction(single, target_structure, self.geometrycontainer, self.domains)
        except testShellError as e:
            testShell = None
            if __debug__:
                log.debug('%s not formed %s' % \
                          ('Interaction(%s, %s)' % (single.pid_particle_pair[0], target_structure),
                          str(e) ))


        ### 2. The testShell was made succesfully. Now make the complete domain
        if testShell:
            # make the domain
            interaction = self.create_interaction(testShell)
            if __debug__:
                log.info('        *Created: %s' % (interaction))

            # get the next event time
            interaction.dt, interaction.event_type, = interaction.determine_next_event()
            if __debug__:
                assert interaction.dt >= 0.0
                log.info('dt=%s' % (FORMAT_DOUBLE % interaction.dt)) 

            self.last_time = self.t

            self.remove_domain(single)
            # the event associated with the single will be removed by the scheduler.

            # add to scheduler
            self.add_domain_event(interaction)
            # check everything is ok
            if __debug__:
                assert self.check_domain(interaction)

            return interaction
        else:
            return None


    def try_transition(self, single, target_structure):
    # Try to form a transition between the 'single' particle on a structure and another structure.
    # This is called after an interaction trial has failed.
    # Currently only transitions between two planar surfaces are supported.

        if __debug__:
           log.debug('trying to form Transition(%s, %s)' % (single.pid_particle_pair[1], target_structure))


        ### 1. Attempt to make a testShell. If it fails it will throw an exception.
        try:
            testShell = try_default_testtransition(single, target_structure, self.geometrycontainer, self.domains)
        except testShellError as e:
            testShell = None
            if __debug__:
                log.debug('%s not formed %s' % \
                          ('Transition(%s, %s)' % (single.pid_particle_pair[0], target_structure),
                          str(e) ))


        ### 2. The testShell was made succesfully. Now make the complete domain
        if testShell:
            # make the domain
            transition = self.create_transition(testShell)
            if __debug__:
                log.info('        *Created: %s' % (transition))

            # get the next event time
            transition.dt, transition.event_type, = transition.determine_next_event()
            if __debug__:
                assert transition.dt >= 0.0
                log.info('dt=%s' % (FORMAT_DOUBLE % transition.dt)) 

            self.last_time = self.t

            self.remove_domain(single)
            # the event associated with the single will be removed by the scheduler.

            # add to scheduler
            self.add_domain_event(transition)
            # check everything is ok
            if __debug__:
                #assert self.check_domain(transition) # HACK
                pass

            return transition
        else:
            return None


    def try_pair(self, single1, single2):
    # Try to make a pair domain out of the two singles. A pair domain can be a:
    # -SimplePair (both particles live on the same structure)
    # -MixedPair (the particles live on different structures -> MixedPair2D3D)
    # Note that single1 is always the single that initiated the creation of the pair
    #  and is therefore no longer in the scheduler

        if __debug__:
            log.debug('trying to form Pair(%s, %s)' %
                (single1.pid_particle_pair, single2.pid_particle_pair))

        ### 1. Attempt to make a testShell. If it fails it will throw an exception.
        try:
            testShell = try_default_testpair(single1, single2, self.geometrycontainer, self.domains)
        except testShellError as e:
            testShell = None
            if __debug__:
                log.debug('%s not formed %s' % \
                          ('Pair(%s, %s)' % (single1.pid_particle_pair[0], single2.pid_particle_pair[0]),
                          str(e) ))


        ### 2. If the testShell could be formed, make the complete domain.
        if testShell:
            pair = self.create_pair(testShell)
            if __debug__:
                log.info('        *Created: %s' % (pair))

            pair.dt, pair.event_type, pair.reactingsingle = pair.determine_next_event(pair.r0)
            if __debug__:
                assert pair.dt >= 0
                log.info('dt=%s' % (FORMAT_DOUBLE % pair.dt)) 

            self.last_time = self.t

            self.remove_domain(single1)
            self.remove_domain(single2)

            # single1 will be removed by the scheduler, since it initiated the creation of the
            # pair.
            self.remove_event(single2)

            # add to scheduler
            self.add_domain_event(pair)
            # check everything is ok
            if __debug__:
                assert self.check_domain(pair)

            return pair
        else:
            return None
    

    ##############################################
    #### METHODS FOR MULTI SHELL CONSTRUCTION ####
    ##############################################
    def form_multi(self, single, multi_partners):
        # form a Multi with the 'single'

        # The 'multi_partners' are neighboring NonInteractionSingles, Multis and surfaces which
        # can be added to the Multi (this can also be empty)
        for partner, overlap in multi_partners:
            assert ((isinstance(partner, NonInteractionSingle) and partner.is_reset()) or \
                    isinstance(partner, Multi) or \
                    isinstance(partner, PlanarSurface) or \
                    isinstance(partner, DiskSurface) or \
                    isinstance(partner, CylindricalSurface)) or \
                    isinstance(partner, SphericalSurface), \
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

        if __debug__:
            assert self.check_domain(multi)

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

            # neighbor_distances has same format as before
            # It contains the updated list of domains in the neighborhood with distances
            neighbor_distances = self.geometrycontainer.get_neighbor_domains(dompos, self.domains,
                                                                             ignore=[domain.domain_id, multi.domain_id])

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
                                    MULTI_SHELL_FACTOR
                    # distance from the center of the particles/domains
                    distance = self.world.distance(dompos, neighbor.shell.shape.position)
                    overlap = distance - multi_horizon

                elif isinstance(neighbor, Multi):
                    # The dist_to_shell = dist_to_particle - multi_horizon_of_target_particle
                    # So only the horizon and distance of the current single needs to be taken into account
                    # Note: this is built on the assumption that the shell of a Multi has the size of the horizon.
                    multi_horizon = (domain.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR)
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
                                       single.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR,
                                       multi.domain_id)
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


    ##############################
    #### SOME HELPER METHODS  ####
    ##############################
    def domain_distance(self, pos, domain):
        # calculates the shortest distance from 'pos' to A shell (a Multi
        # can have multiple) of 'domain'.
        # Note: it returns the distance between pos and the surface of the shell
        return min(self.world.distance(shell.shape, pos)
                   for i, (_, shell) in enumerate(domain.shell_list))


    def obj_distance_array(self, pos, domains):
        dists = numpy.array([self.domain_distance(pos, domain) for domain in domains])
        return dists
            

    ###########################
    #### STATISTICS REPORT ####
    ###########################
    def print_report(self, out=None):
        """Print various statistics about the simulation.
        
        Arguments:
            - None

        """
        total_steps = self.step_counter
        single_steps = numpy.array(self.single_steps.values()).sum()
        interaction_steps = numpy.array(self.interaction_steps.values()).sum()
        pair_steps = numpy.array(self.pair_steps.values()).sum()
        multi_steps = self.multi_steps[3] # total multi steps

        report = '''
t = %g
steps = %d 
\tSingle:\t%d\t(%.2f %%)\t(escape: %d, reaction: %d, bursted: %d)
\tInteraction: %d\t(%.2f %%)\t(escape: %d, interaction: %d, bursted: %d)
\tPair:\t%d\t(%.2f %%)\t(r-escape: %d, R-escape: %d, reaction pair: %d, single: %d, bursted: %d)
\tMulti:\t%d\t(%.2f %%)\t(diffusion: %d, escape: %d, reaction pair: %d, single: %d, bursted: %d)
total reactions = %d
rejected moves = %d
''' \
            % (self.t, total_steps,
               single_steps,
               (100.0*single_steps) / total_steps,
               self.single_steps[EventType.SINGLE_ESCAPE],
               self.single_steps[EventType.SINGLE_REACTION],
               self.single_steps[EventType.BURST],
               interaction_steps,
               (100.0*interaction_steps) / total_steps,
               self.interaction_steps[EventType.IV_ESCAPE],
               self.interaction_steps[EventType.IV_INTERACTION],
               self.interaction_steps[EventType.BURST],
               pair_steps,
               (100.0*pair_steps) / total_steps,
               self.pair_steps[EventType.IV_ESCAPE],
               self.pair_steps[EventType.COM_ESCAPE],
               self.pair_steps[EventType.IV_REACTION],
               self.pair_steps[EventType.SINGLE_REACTION],
               self.pair_steps[EventType.BURST],
               multi_steps,
               (100.0*multi_steps) / total_steps,
               self.multi_steps[EventType.MULTI_DIFFUSION],
               self.multi_steps[EventType.MULTI_ESCAPE],
               self.multi_steps[EventType.MULTI_BIMOLECULAR_REACTION],
               self.multi_steps[EventType.MULTI_UNIMOLECULAR_REACTION],
               self.multi_steps[EventType.BURST],
               self.reaction_events,
               self.rejected_moves
               )

        print >> out, report


    ##########################################
    #### METHODS FOR CONSISTENCY CHECKING ####
    ##########################################
    def check_domain(self, domain):
        domain.check()

        # Construct ignore list for surfaces and the structures that the domain is associated with
        # __
        # This is necessary because structures which are involved in the 
        # mathematical solution need to be excluded from the overlap 
        # check.
        if isinstance(domain, Multi):
            # Ignore all surfaces, multi shells can overlap with 
            # surfaces.
            ignores = [s.id for s in self.world.structures]
            associated = []
        elif isinstance(domain, SphericalSingle) or isinstance(domain, SphericalPair) or isinstance(domain, PlanarSurfaceSingle):
            # 3D NonInteractionSingles can overlap with planar surfaces but not with rods
            ignores = [s.id for s in self.world.structures if isinstance(s, PlanarSurface)]
            associated = []
        elif isinstance(domain, DiskSurfaceSingle):
            # Particles bound to DiskSurfaces ignore all rods for now. TODO Only ignore neighboring structure
            ignores = [s.id for s in self.world.structures if isinstance(s, CylindricalSurface)]
            associated = [domain.structure.id]
        elif isinstance(domain, CylindricalSurfaceInteraction):
            # Ignore surface of the particle and interaction surface and all DiskSurfaces for now.
            ignores = [s.id for s in self.world.structures if isinstance(s, DiskSurface)]
            associated = [domain.origin_structure.id, domain.target_structure.id]
        elif isinstance(domain, InteractionSingle) or isinstance(domain, MixedPair2D3D) or isinstance(domain, MixedPair1DCap):
            # Ignore surface of the particle and interaction surface
            ignores = []
            associated = [domain.origin_structure.id, domain.target_structure.id]
        elif isinstance(domain, PlanarSurfaceTransitionSingle) or isinstance(domain, PlanarSurfaceTransitionPair):
            # Ignore surface of the particle and interaction surface
            ignores = []
            associated = [domain.structure1.id, domain.structure2.id]
        else:
            # Ignore the structure that the particles are associated with
            ignores = []
            associated = [domain.structure.id]

        ###### check shell consistency
        # check that shells of the domain don't overlap with surfaces that are not associated with the domain.
        # check that shells of the domain DO overlap with the structures that it is associated with.
        # check that shells of the domain do not overlap with shells of other domains
        # check that shells are not bigger than the allowed maximum

        for shell_id, shell in domain.shell_list:

            ### check that shells do not overlap with non associated surfaces
            surfaces = get_all_surfaces(self.world, shell.shape.position, ignores + associated)
            for surface, _ in surfaces:
                assert self.check_shape_overlap(shell.shape, surface.shape) >= 0.0, \
                    '%s (%s) overlaps with %s.' % \
                    (str(domain), str(shell), str(surface))

            ### check that shells DO overlap with associated surfaces.
            surfaces = get_all_surfaces(self.world, shell.shape.position, ignores)
            for surface, _ in surfaces:
                if surface.id in associated:
                    assert self.check_shape_overlap(shell.shape, surface.shape) < 0.0, \
                    '%s (%s) should overlap with associated surface %s.' % \
                    (str(domain), str(shell), str(surface))

            ### Check shell overlap with other shells
            neighbors = self.geometrycontainer.get_neighbor_domains(shell.shape.position,
                                                                    self.domains, ignore=[domain.domain_id])
            # TODO maybe don't check all the shells, this takes a lot of time
            # testing overlap criteria
            for neighbor, _ in neighbors:
                # note that the shell of a MixedPair or Multi that has have just been bursted can stick into each other.
                if not (((isinstance(domain, hasCylindricalShell)   and isinstance(domain, NonInteractionSingle)   and domain.is_reset()) and \
                         (isinstance(neighbor, hasSphericalShell)   and isinstance(neighbor, NonInteractionSingle) and domain.is_reset()) ) or \
                        ((isinstance(domain, hasSphericalShell)     and isinstance(domain, NonInteractionSingle)   and domain.is_reset()) and \
                         (isinstance(neighbor, hasCylindricalShell) and isinstance(neighbor, NonInteractionSingle) and domain.is_reset()) )):

                    for _, neighbor_shell in neighbor.shell_list:
                        overlap = self.check_shape_overlap(shell.shape, neighbor_shell.shape)
                        assert overlap >= 0.0, \
                            '%s (%s) overlaps with %s (%s) by %s.' % \
                            (domain, str(shell), str(neighbor), str(neighbor_shell), FORMAT_DOUBLE % overlap)


            ### checking if the shell don't exceed the maximum size
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


    def check_shape_overlap(self, shape1, shape2):
    # Return the distance between two shapes (positive number means no overlap, negative means overlap)

        # overlap criterium when one of the shapes is a sphere
        if (type(shape1) is Sphere):
            distance = self.world.distance(shape2, shape1.position)
            return distance - shape1.radius
        elif (type(shape2) is Sphere):
            distance = self.world.distance(shape1, shape2.position)
            return distance - shape2.radius

        # overlap criterium when one of the shapes is a cylinder and the other is a plane
        elif ((type(shape1) is Cylinder) and (type(shape2) is Plane)) or \
             ((type(shape2) is Cylinder) and (type(shape1) is Plane)):

            # putatively swap to make sure that shape1 is Cylinder and shape2 is Plane
            if (type(shape2) is Cylinder) and (type(shape1) is Plane):
                temp = shape2
                shape2 = shape1
                shape1 = temp

            if math.sqrt(shape1.half_length**2 + shape1.radius**2) < self.world.distance(shape2, shape1.position):
                # if the plane is outside the sphere surrounding the cylinder.
                return 1
            else: 
                shape1_pos = shape1.position
                shape2_pos = shape2.position
                shape1_post = self.world.cyclic_transpose(shape1_pos, shape2_pos)
                inter_pos = shape1_post - shape2_pos
                relative_orientation = numpy.dot(shape1.unit_z, shape2.unit_z)
                # if the unit_z of the plane is perpendicular to the axis of the cylinder.
                if feq(relative_orientation, 0):
                    if abs(numpy.dot(inter_pos, shape2.unit_x)) - shape1.half_length > shape2.half_extent[0] or \
                       abs(numpy.dot(inter_pos, shape2.unit_y)) - shape1.half_length > shape2.half_extent[1] or \
                       abs(numpy.dot(inter_pos, shape2.unit_z)) > shape1.radius:
                        return 1
                    else:
                        # NOTE this is still incorrect. By the above if statement the cylinder is made into a cuboid
                        return -1

                elif feq(abs(relative_orientation), 1):
                    if abs(numpy.dot(inter_pos, shape2.unit_x)) - shape1.radius > shape2.half_extent[0] or \
                       abs(numpy.dot(inter_pos, shape2.unit_y)) - shape1.radius > shape2.half_extent[1] or \
                       abs(numpy.dot(inter_pos, shape2.unit_z)) > shape1.half_length:
                        return 1
                    else:
                        # NOTE this is still incorrect. By the above if statement the cylinder is made into a cuboid
                        return -1
                else:
                    raise RuntimeError('check_shape_overlap: planes and cylinders should be either parallel or perpendicular.')

        # FIXME FIXME separate Disk overlap criterium in separate check (probably clearer/faster than combining)
        # overlap criterium when both shells are cylindrical or disks (= cylinders with half_length = 0.0)
        elif (type(shape1) is Cylinder) and (type(shape2) is Cylinder) or \
             (type(shape1) is Cylinder) and (type(shape2) is Disk) or \
             (type(shape1) is Disk)     and (type(shape2) is Cylinder) or \
             (type(shape1) is Disk)     and (type(shape2) is Disk) :

            # Disk has no property half_length, so check type before any calculation
            shape1_hl = 0.0
            shape2_hl = 0.0
            if type(shape1) is Cylinder:
                shape1_hl = shape1.half_length
            if type(shape2) is Cylinder:                
                shape2_hl = shape2.half_length

            if math.sqrt(shape1_hl**2 + shape1.radius**2) < self.world.distance(shape2, shape1.position) or \
               math.sqrt(shape2_hl**2 + shape2.radius**2) < self.world.distance(shape1, shape2.position):
                # if cylinder2 is outside the sphere surrounding the cylinder.
                return 1
            else: 
                shape1_pos = shape1.position
                shape2_pos = shape2.position
                shape1_post = self.world.cyclic_transpose(shape1_pos, shape2_pos)
                inter_pos = shape1_post - shape2_pos

                relative_orientation = numpy.dot(shape1.unit_z, shape2.unit_z)
                # if the unit_z of cylinder1 (disk1) is perpendicular to the axis of cylinder2 (disk2).
                if feq(relative_orientation, 0.0):
                    return 1                                      # TODO Why does this always return 1 ?


                elif feq(abs(relative_orientation), 1.0):
                    inter_pos_z = shape2.unit_z * numpy.dot(inter_pos, shape2.unit_z)
                    inter_pos_r = inter_pos - inter_pos_z
                    overlap_r = length(inter_pos_r) - (shape1.radius + shape2.radius)
                    overlap_z = length(inter_pos_z) - (shape1_hl + shape2_hl)
                    if (overlap_r <= 0.0) and (overlap_z <= 0.0):
                        if (overlap_r == 0.0) or (overlap_z == 0.0):
                            # Special treatment of "direct overlap", i.e. e.g. two overlapping disks
                            return -1
                        else:
                            return -(overlap_r * overlap_z)        # TODO: find better number for overlap measure
                    else:
                        return 1
                else:
                    raise RuntimeError('check_shape_overlap: cylinders are not oriented parallel or perpendicular.')

        # If both shapes are planes.
        elif (type(shape1) is Plane) and (type(shape2) is Plane):
            # TODO
            return 1

        # something was wrong (wrong type of shape provided)
        else:
            raise RuntimeError('check_shape_overlap: wrong shape type(s) provided.')


    def check_domain_for_all(self):
    # For every domain in domains, checks the consistency
        for domain in self.domains.itervalues():
            self.check_domain(domain)


    def check_event_stoichiometry(self):
        ### check scheduler
        # -check that every event on scheduler is associated with domain or is a non domain related event
        # -check that every event with domain is also in domains{}

        SMT = 6                 # just some temporary definition to account for non-domain related events
        self.other_objects = {} # ditto
        for id, event in self.scheduler:
            domain_id = event.data
            if domain_id != SMT:
                # in case the domain_id is not special, it needs to be in self.domains
                assert domain_id in self.domains
            else:
                # if it is special in needs to be in self.other_objects
                assert domain_id in self.other_objects


        ### performs two checks:
        # -check that every domain in domains{} has one and only one event in the scheduler
        # -check that dt, last_time of the domain is consistent with the event time on the scheduler
        already_in_scheduler = set()        # make a set of domain_ids that we found in the scheduler
        for domain_id, domain in self.domains.iteritems():
            # find corresponding event in the scheduler
            for id, event in self.scheduler:
                if event.data == domain_id:
                    # We found the event in the scheduler that is associated with the domain
                    assert (abs(domain.last_time + domain.dt - event.time) <= TIME_TOLERANCE * event.time)
                    break
            else:
                raise RuntimeError('Event for domain_id %s could not be found in scheduler' % str(domain_id) )

            # make sure that we had not already found it earlier
            assert domain_id not in already_in_scheduler
            already_in_scheduler.add(domain_id)


    def check_particle_consistency(self):
        ### Checks the particles consistency between the world and the domains of
        #   the simulator

        # First, check num particles in domains is equal to num particles in the world
        world_population = self.world.num_particles
        domain_population = 0
        for domain in self.domains.itervalues():
            domain_population += domain.multiplicity

        if world_population != domain_population:
            raise RuntimeError('population %d != domain_population %d' %
                               (world_population, domain_population))

       ## Next, check that every particle in the domains is once and only once in the world.
       #already_in_domain = set()
       #world_particles = list(self.world)
       #for domain in self.domains.itervalues():
           #for domain_particle in domain.particles:
               #for world_particle in world_particles:
                   #if domain_particle == world_particle:
                       #break
               #else:
                   #raise RuntimeError('particle %s not in world' % str(domain_particle))

               #assert domain_particle not in already_in_domain, 'particle %s twice in domains' % str(domain_particle)
               #already_in_domain.add(domain_particle)

       ## This now means that all the particles in domains are represented in the world and that
       ## all the particles in the domains are unique (no particle is represented in two domains)

       #assert already_in_domain == len(world_particles)
       ## Since the number of particles in the domains is equal to the number of particles in the
       ## world this now also means that all the particles in the world are represented in the
       ## domains.

    def check_shell_matrix(self):
        ### checks the consistency between the shells in the geometry container and the domains
        # -check that shells of a domain in domains{} are once and only once in the shell matrix
        # -check that all shells in the matrix are once and only once in a domain

        did_map, shell_map = self.geometrycontainer.get_dids_shells()

        shell_population = 0
        for domain in self.domains.itervalues():
            shell_population += domain.num_shells

            # get the shells associated with the domain according to the geometrycontainer
            shell_ids = did_map[domain.domain_id]
            for shell_id, shell in domain.shell_list:
                # check that all shells in the domain are in the matrix and have the proper domain_id
                assert shell_id in shell_ids

            # check that the number of shells in the shell_list is equal to the number of shells the domain
            # says it has is equal to the number of shell the geometrycontainer has registered to the domain.
            assert domain.num_shells == len(list(domain.shell_list))
            assert domain.num_shells == len(shell_ids)
            if len(shell_ids) != domain.num_shells:
                diff = set(sid for (sid, _)
                               in domain.shell_list).difference(shell_ids)
                for sid in diff:
                    print shell_map.get(sid, None)

                raise RuntimeError('number of shells are inconsistent '
                                   '(%d != %d; %s) - %s' %
                                   (len(shell_ids), domain.num_shells, 
                                    domain.domain_id, diff))

        # check that total number of shells in domains==total number of shells in shell container
        matrix_population = self.geometrycontainer.get_total_num_shells()
        if shell_population != matrix_population:
            raise RuntimeError('num shells (%d) != matrix population (%d)' %
                               (shell_population, matrix_population))
        # Since the number of shells in the domains is equal to the number of shells in the
        # geometrycontainer this now also means that all the shells in the geometrycontainer are
        # represented in the domains.

        # This now also means that the shells of all the cached singles in Pairs are not in the containers.


#####################################################
#### TODO: METHODS PROTOTYPES FOR FURTHER CHECKS ####
#####################################################
### check that all the cached singles in Pairs are not in domains{}.

### check world
# check that particles do not overlap with cylinderical surfaces unless the domain is associated with the surface

### check domains
# check that testShell is of proper type
# check that the shell(s) of domain have the proper geometry (sphere/cylinder)
# check that shell is within minimum and maximum
# check that shell obeys scaling properties
# check that the particles are within the shells of the domain
# check that the position of the shell(s) and particle(s) are within the structures that the particles live on
# check that dt != 0 unless NonInteractionSingle.is_reset()
# check associated structures are of proper type (plane etc, but also of proper StructureTypeID)
# check that Green's functions are defined.

### check NonInteractionSingle
# check if is_reset() then shellsize is size particle, dt==0, event_type==ESCAPE
# check shell.unit_z == structure.unit_z
# check drift < inf
# check inner space < shell

### check Pair
# check number of particles==2
# check D_R, D_r > 0
# check r0 is between sigma and a (within bounds of greens function)
# check sigma

### check SimplePair
# check CoM, iv lies in the structure of the particles
# check shell.unit_z == structure.unit_z

### check MixedPair
# check CoM is in surface, IV is NOT

### check Interaction
# check shell.unit_z = +- surface.unit_z

### check consistency of the model
# check that the product species of an interaction lives on the interaction surface.
# check that the product species of a bimolecular reaction lives on either structure of reactant species
# check that the product species of a unimolecular reaction lives on structure of the reactant species or the bulk.
# structures cannot overlap unless substructure

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

    ##########################
    #### CHECKLIST METHOD ####
    ##########################
    def check(self):
        ParticleSimulatorBase.check(self)

        assert self.scheduler.check()

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.check_shell_matrix()
        self.check_domains()
        self.check_event_stoichiometry()
        self.check_particle_consistency()
        
        self.check_domain_for_all()


    ###############################
    #### METHODS FOR DEBUGGING ####
    ###############################
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


