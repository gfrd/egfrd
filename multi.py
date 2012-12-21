from weakref import ref

from gfrdbase import *
from utils import *
import itertools

import _gfrd
from _gfrd import (SphericalShell, Sphere)
from constants import EventType
from domain import Domain
from shells import (hasSphericalShell, Others)

import os

import logging
log = logging.getLogger('ecell')


class Multi(Domain, hasSphericalShell, Others):
    def __init__(self, domain_id, main, step_size_factor, dt_hardcore_min = -1):
        Domain.__init__(self, domain_id)
        # note: hasSphericalShell.__init__ is not called and that the Multi has no self.testShell
        #       Only the methods of the hasSphericalShell class are used

        self.main = ref(main)
        self.last_event = EventType.MULTI_DIFFUSION     #None
#        self.event_type = EventType.MULTI_DIFFUSION

        self.sphere_container = _gfrd.SphericalShellContainer(main.world.world_size, 3)
        self.particle_container = _gfrd.MultiParticleContainer(main.world)
        self.escaped = False
        self.step_size_factor = step_size_factor
        self.dt_hardcore_min = dt_hardcore_min
        self.last_reaction = None
        self.reaction_length = 0

    def initialize(self, t):
        self.last_time = t
        self.start_time = t
        main = self.main()
        self.set_dt_and_reaction_length()

    def calculate_bd_dt(species_list):
        D_list = []
        radius_list = []
        D_max = 0.
        radius_min = numpy.inf
        for species in species_list:
            if D_max < species.D:
                D_max = species.D
            if radius_min > species.radius:
                radius_min = species.radius
        return (radius_min * 2) ** 2 / (D_max * 2)

    def get_multiplicity(self):
        return self.particle_container.num_particles
    multiplicity = property(get_multiplicity)

    def within_shell(self, pp):
        # test if pid_particle_pair 'pp' is within at least one shell.
        shells = self.sphere_container.get_neighbors_within_radius(pp[1].position, -(pp[1].radius + self.reaction_length))
        return len(shells) >= 1

    def add_shell(self, shell_id_shell_pair):
        if __debug__:
            log.info("add shell to multi:\t (%s,\t%s)" %
                     (shell_id_shell_pair[0], shell_id_shell_pair[1]))
        self.sphere_container.update(shell_id_shell_pair)

    def add_particle(self, pid_particle_pair):
        if __debug__:
            log.info("add particle to multi:\t (%s,\t\t%s)" % 
                     (pid_particle_pair[0], pid_particle_pair[1]))
        self.particle_container.update_particle(pid_particle_pair)
        
    def set_dt_and_reaction_length(self):
        main = self.main()
        self.dt, self.reaction_length = \
            self.particle_container.determine_dt_and_reaction_length(main.network_rules, self.step_size_factor, self.dt_hardcore_min)

    def step(self):
        self.escaped = False
        tx = self.particle_container.create_transaction()
        main = self.main()

        class check_reaction(object):
            def __init__(self):
                self.reactions = []

            def __call__(self, ri):
                self.reactions.append(ri)

        class clear_volume(object):
            def __init__(self, outer):
                self.outer_ = outer

            def __call__(self, shape, ignore0, ignore1=None):

                within_radius = self.outer_.sphere_container.get_neighbors_within_radius(
                               shape.position, -(shape.radius + self.outer_.reaction_length) )

                if list(within_radius) == []:
                   # Convert to list first because of tricky C++/Python binding here
                   # Type casting bool(within_radius) caused segfaults in some simulations
                    main = self.outer_.main()
                    if self.outer_.last_event == EventType.MULTI_DIFFUSION: #None:
                        self.outer_.last_event = EventType.MULTI_ESCAPE
                    main.burst_volume(shape.position, shape.radius, 
                                      ignore=[self.outer_.domain_id, ])
                    if ignore1 is None:
                        return not main.world.check_overlap(
                            (shape.position, shape.radius), ignore0)
                    else:
                        return not main.world.check_overlap(
                            (shape.position, shape.radius), ignore0, ignore1)

                return True

        cr = check_reaction()
        vc = clear_volume(self)

        ppg = _gfrd.newBDPropagator(tx, main.network_rules,
                     myrandom.rng, self.dt, main.dissociation_retry_moves, self.reaction_length,
                     cr, vc, [pid for pid, _ in self.particle_container])

        self.last_event = EventType.MULTI_DIFFUSION     #None
        while ppg():
            if cr.reactions:
                self.last_reaction = cr.reactions[-1]
                if len(self.last_reaction.reactants) == 1:
                    self.last_event = EventType.MULTI_UNIMOLECULAR_REACTION
                else:
                    self.last_event = EventType.MULTI_BIMOLECULAR_REACTION
                break
                
        # for pid_particle_pair in itertools.chain(
        #         tx.modified_particles, tx.added_particles):
        #     overlapped = main.world.check_overlap(pid_particle_pair[1].shape, pid_particle_pair[0])
        #     if overlapped:
        #         if __debug__:
        #             log.info("collision occurred between particles of a multi and the outside: %s - %s.  moves will be rolled back." % (pid_particle_pair, list(overlapped)))
        #         tx.rollback()
        #         return

        #     if not self.within_shell(pid_particle_pair):
        #         if self.last_event == EventType.MULTI_DIFFUSION   #None:
        #             self.last_event = EventType.MULTI_ESCAPE
        #         main.clear_volume(
        #             pid_particle_pair[1].position,
        #             pid_particle_pair[1].radius, ignore=[self.domain_id, ])

    def check(self):
        # check Multi num particles >= 1
        assert self.multiplicity >= 1, 'Multi needs to have at least one particle'
        # check number of shells==num of particles in Multi
        assert self.multiplicity == self.num_shells, \
            'number of particles is not equal to number of shells'
        # Event is DIFFUSION, SINGLE_REACTION, BIMOL_REACTION, ESCAPE
        assert self.last_event == EventType.MULTI_DIFFUSION or \
               self.last_event == EventType.MULTI_ESCAPE or \
               self.last_event == EventType.MULTI_UNIMOLECULAR_REACTION or \
               self.last_event == EventType.MULTI_BIMOLECULAR_REACTION, \
                    'event_type of Multi was not of proper type'

        # check that the shells of the multi are contiguous
        # FIXME: This can be done more efficiently because we already have an
        #        exhaustive list of the shells to check
        def check_connected(shell_id_shell_pair, already_checked):

            shell_id, shell = shell_id_shell_pair
            if shell_id not in already_checked:
                already_checked.add(shell_id)
                # TODO get only the neighbors that are not in already_checked
                neighbor_distances = self.sphere_container.get_neighbors_within_radius(shell.shape.position, shell.shape.radius)
                for shell_id_shell_pair, _ in neighbor_distances:
                    already_checked = check_connected(shell_id_shell_pair, already_checked)

            return already_checked

        connected_shell_ids = check_connected(self.shell_list.next(), set([]))
        assert len(connected_shell_ids) == self.num_shells, \
            'shells of %s are not contiguous.' % str(self)


        # check that all particles are within one of the shells.
        for pid_particle_pair in self.particle_container:
            assert self.within_shell(pid_particle_pair),\
                'not all particles within the shell.'

        # check that all the shells in the private multi shell container
        # are also in the shell container of the Simulator
        main = self.main()
        for shell_id, shell in self.shell_list:
            container = main.geometrycontainer.get_container(shell)
            if not container.contains(shell_id):
                raise RuntimeError,\
                    'self.sim.main.sphere_container does not contain %s'\
                    % str(shell_id)
        # check that all the shells in the container of the Simulator that belong to this multi
        # domain are also in the private shell container.
        for shell_id, shell in main.geometrycontainer.containers[0]:
            if shell.did == self.domain_id:
                if not self.sphere_container.contains(shell_id):
                    raise RuntimeError,\
                        'self.sphere_container does not contain %s'\
                        % str(shell_id)

    def __repr__(self):
        return 'Multi[domain_id=%s, event_id=%s,\n    %s,\n    %s]' % (
            self.domain_id, self.event_id,
            ',\n    '.join('(%s,\n     %s)' % (p[0], repr(p[1])) for p in self.particle_container),
            ',\n    '.join('(%s,\n     %s)' % (s[0], repr(s[1])) for s in self.sphere_container)
            )

    def has_particle(self, pid):
        return pid in [pid_particle_pair[0] for pid_particle_pair in self.particles]
        # try:
        #     self.particle_container.get_particle(pid)
        #     return True
        # except:
        #     return False

    def get_particles(self):
        return iter(self.particle_container)
    particles = property(get_particles)

    def num_shells(self):
        return len(self.sphere_container)
    num_shells = property(num_shells)    

    def shell_list(self):
    # returns a list of (shell_id, shell) pairs of the Multi
        return iter(self.sphere_container)
    shell_list = property(shell_list)
