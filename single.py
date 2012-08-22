from _gfrd import (
    SphericalShell,
    CylindricalShell,
    Sphere,
    Cylinder,
    Disk,
    CuboidalRegion,
    PlanarSurface,
    CylindricalSurface,
    DiskSurface,
    )
from _greens_functions import *
from greens_function_wrapper import *
from constants import EventType
from utils import *
from transitiontools import *

from domain import (
    Domain,
    ProtectiveDomain)

from shells import *

__all__ = [
    'CylindricalSurfaceSingle',
    'DiskSurfaceSingle',
    'PlanarSurfaceSingle',
    'SphericalSingle',
    'Single',
    'NonInteractionSingle',
    'InteractionSingle',
    'CylindricalSurfaceInteraction',
    'CylindricalSurfaceCapInteraction',
    'CylindricalSurfaceSink',
    'PlanarSurfaceInteraction',
    'TransitionSingle',
    'PlanarSurfaceTransitionSingle',
    ]


class Single(ProtectiveDomain):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see 
    coordinate.py. For each coordinate the Green's function is 
    specified.

    """
    def __init__(self, domain_id, shell_id, reactionrules):
        ProtectiveDomain.__init__(self, domain_id, shell_id)

        assert self.testShell                           # make sure that the hasShell class was initialized earlier
        assert isinstance(self.testShell, testSingle)

        self.multiplicity = 1                           # all singles have only one particle

        self.pid_particle_pair = self.testShell.pid_particle_pair      # the particle in the single
        self.structure = self.testShell.structure                      # the structure in/on which the particle lives

        self.reactionrules = reactionrules              # the reaction rules for mono molecular reactions
        self.rrule = None                               # the reaction rule of the actual reaction (property reactionrule)
        self.k_tot = self.calc_ktot(reactionrules)

    def initialize(self, t):
        self.dt = 0.0
        self.last_time = t
        self.event_type = None

    def getD(self):
        return self.pid_particle_pair[1].D
    D = property(getD)

    def get_reaction_rule(self):
        if self.rrule == None:
            self.rrule = self.draw_reaction_rule(self.reactionrules)
        return self.rrule
    reactionrule = property(get_reaction_rule)

    def get_inner_a(self):
        # inner_a is the outer radius of the 'inner shell', which demarkates
        # the space where the CoM of the particle can go. The space can be
        # one, two or three dimensional and can be linear, round, spherical or
        # cylindric in form.
        return self.shell.shape.radius - self.pid_particle_pair[1].radius

    def draw_new_position(self, dt, event_type):
        '''Draw a new position of the particle in the coordinate system given that an
        event of event_type has occured. No periodic boundary conditions are applied yet
        to the coordinates, as the Domains do not understand this.
        The method also returns the structure_id of the structure on which the particle
        is located at the time of the event.

        Note also that it calculates the positional information PRIOR to any identity or
        structure change. In case of a structure change the structure_id therefor is a
        reference to the original structure.

        When for example a particles dissociates from the membrane (a structure change),
        it is placed next to the membrane. 'draw_new_position' will only determine the
        location in the membrane where it dissociates, not the position where the particle
        should be placed.
        '''
        pass

    def is_reset(self):
        return self.dt == 0.0 and self.event_type == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple(self):
        # MOVE to ProtectiveDomain class as this is actually a property of a single particle
        # and not of a Single domain. It should be a method to a particle (or species)
        # but that's too much hassle. Now make it a function (!) of Domain class for
        # easy access. This mechanism is used everywhere a particle can decay.
        """Return a (reaction time, event type)-tuple for a unimolecular reaction.

        """
        if self.k_tot <= 0:
            dt = numpy.inf
        elif self.k_tot == numpy.inf:
            dt = 0.0
        else:
            rnd = myrandom.uniform()
            if rnd == 0:
                dt = numpy.inf
            else:
                dt = (1.0 / self.k_tot) * (- math.log(rnd)) # log(1/x) == - log(x)
        return dt, EventType.SINGLE_REACTION

    def draw_escape_time_tuple(self):
        """Return an (escape time, event type)-tuple.

        In case of an interaction single, handles only the escape in the
        non-interacting coordinate.
        
        """
        if self.getD() == 0:
            dt = numpy.inf
        else:
            dt = draw_time_wrapper(self.greens_function())

        event_type = EventType.SINGLE_ESCAPE
        return dt, event_type

    def get_particles(self):
        return [self.pid_particle_pair]
    particles = property(get_particles)

    def check(self):
        # performs an internal consistency check
        assert ProtectiveDomain.check(self)
        return True

    def __str__(self):
        pid = self.pid_particle_pair[0]
        sid = self.pid_particle_pair[1].sid
        name = self.world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        return 'Single[%s: %s, ST%s, %s]' % (self.domain_id, pid, name,
                                             self.pid_particle_pair[1].position)


class NonInteractionSingle(Single, NonInteractionSingles):
    """1 Particle inside a shell, no other particles around. 

    There are 3 types of NonInteractionSingles:
        * SphericalSingle: spherical shell, 3D movement.
        * PlanarSurfaceSingle: cylindrical shell, 2D movement.
        * CylindricalSurfaceSingle: cylindrical shell, 1D movement.

    NonInteractionSingle hold two independent processes:
        * Diffusion of the particle in n (1-3) dimensions
        * Decay process of the particle

    """
    def __init__(self, domain_id, shell_id, reactionrules):
        Single.__init__(self, domain_id, shell_id, reactionrules)
        # Note: for the NonInteractionSingles superclass nothing is to be initialized.

    def initialize(self, t):
        Single.initialize(self, t)
        self.event_type = EventType.SINGLE_ESCAPE


    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple())

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos
        else:
            # calculate r
            if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  

                r = self.get_inner_a()
            else:
                gf = self.greens_function()
                r = draw_r_wrapper(gf, dt, self.get_inner_a())
                # Note that in case of 1D diffusion r has a direction. It is
                # the lateral displacement and can be positive or negative. 
                # In other cases r is always positive and denotes
                # radial displacement from the center.

            # calculate other coordinates
            displacement = self.create_position_vector(r)

            # This should be checked in the unit test of random_vector.
            # if __debug__:
            #     scale = self.pid_particle_pair[1].radius
            #     if feq(length(displacement), abs(r), typical=scale) == False:
            #         raise AssertionError('displacement != abs(r): %g != %g.' % 
            #                              (length(displacement), abs(r)))

            # Add displacement to shell.shape.position, not to particle.position.  
            # This distinction is important only in the case of an 
            # asymmetric 1D domain (r0 != 0, or drift), since draw_r always 
            # returns a distance relative to the centre of the shell (r=0), 
            # not relative to r0.

            # note that we need to make sure that the shell.shape.position and displacement vector
            # are in the structure to prevent the particle leaving the structure due to numerical errors
            newpos = self.shell.shape.position + displacement
    
        # The structure on which the particle ended is always the same as on which it began.
        structure_id = self.structure.id

        return newpos, structure_id

### TODO move this in some way to the shell making stuff in shells.py
    def calculate_shell_size_to_single(self, closest, distance_to_shell, geometrycontainer):
        # only if the closest shell is also a 'simple' Single
        assert isinstance(closest, NonInteractionSingle)

        min_radius1 = self.get_min_shell_size()

        D1 = self.getD()
        if D1 == 0:
            return min_radius1

        min_radius2 = closest.pid_particle_pair[1].radius
        D2 = closest.getD()
        min_radius12 = min_radius1 + min_radius2
        sqrtD1 = math.sqrt(D1)


        singlepos = self.pid_particle_pair[1].position
        closestpos = closest.shell.shape.position
        distance = geometrycontainer.distance(singlepos, closestpos)

        shell_size = min(sqrtD1 / (sqrtD1 + math.sqrt(D2))
                        * (distance - min_radius12) + min_radius1,
                        distance_to_shell / SAFETY)
        if shell_size < min_radius1:
            shell_size = min_radius1

        return shell_size


class SphericalSingle(NonInteractionSingle, hasSphericalShell):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector:
          theta, phi.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):
        assert isinstance(testShell, SphericalSingletestShell)

        assert isinstance(testShell, SphericalSingletestShell)
        hasSphericalShell.__init__(self, testShell, domain_id)           # here also the shell is created
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def greens_function(self):
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_inner_a())

    # specific shell methods for the SphericalSingle
    # these return potentially corrected dimensions
    # For explanation see NonInteractionSingles and Others in shells.py
    def shell_list_for_single(self):
        # make sure that there is at least space for the multi shell
        min_radius = self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR
        if  self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            fake_shell = self.create_new_shell(position, min_radius, self.domain_id)
            # TODO maybe just returning the dimension instead of actually making the shell is cheaper.
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        # keeps all the other domains at a distance of the burst radius
        min_radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            fake_shell = self.create_new_shell(position, min_radius, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        try:
            radius = self.testShell.determine_possible_shell(self.structure.id, [self.domain_id], [])
            return self.create_new_shell(position, radius, self.domain_id)

        except ShellmakingError as e:
            raise Exception('SphericalSingle, create_updated_shell failed: %s' %
                            (str(e)))

    def create_position_vector(self, r):
        return random_vector(r)

    def check(self):
        assert ProtectiveDomain.check(self)
        # check shell stuff
        assert isinstance(self.testShell, SphericalSingletestShell)
        assert isinstance(self.shell.shape, Sphere)
        assert self.shell.shape.radius <= self.testShell.get_max_radius()
        if not self.is_reset():
            assert self.shell.shape.radius*SAFETY >= self.testShell.get_min_radius()
        if self.is_reset():
            assert self.shell.shape.radius == self.pid_particle_pair[1].radius

        assert self.is_reset() ^ (self.dt != 0.0)
        assert isinstance(self.structure, CuboidalRegion)
        assert self.greens_function

        return True

    def __str__(self):
        return 'Spherical' + Single.__str__(self)


class PlanarSurfaceSingle(NonInteractionSingle, hasCylindricalShell):
    """1 Particle inside a (cylindrical) shell on a PlanarSurface. (Hockey 
    pucks).

        * Particle coordinates on surface: x, y.
        * Domain: radial r. (determines x and y together with theta).
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        assert isinstance(testShell, PlanarSurfaceSingletestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def greens_function(self):
        return GreensFunction2DAbsSym(self.D,
                                          self.get_inner_a())

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    # For explanation see NonInteractionSingles and Others in shells.py
    def shell_list_for_single(self):
        # The shell should always be larger that the bare minimum for a test shell
        # Note that in case of a cylindrical shell this must fit inside of the spherical multi shell because
        # this demarks the boundary outside of which domains can be.
        min_radius = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.radius < min_radius:
            # This should only happen if the domain is just initialized.
            position = self.shell.shape.position
            radius = self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR
            # There should always be space to make the spherical multi shell.
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            # When the shell has proper dimensions -> just return the real shell.
            return self.shell_list

    def shell_list_for_other(self):
        # set the threshold below which the shell is adjusted.
        min_radius = self.pid_particle_pair[1].radius * math.sqrt(SINGLE_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR     # The burst radius
            # keep everything at a burstradius distance
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            # When the shell is larger than the burst volume -> return the real shell.
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        try:

            dr, dz_right, dz_left = self.testShell.determine_possible_shell(self.structure.id, [self.domain_id], [])
            center, radius, half_length = self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                                             self.testShell.get_orientation_vector(),
                                                                             dr, dz_right, dz_left) 
            
            return self.create_new_shell(center, radius, half_length, self.domain_id)
        except ShellmakingError as e:
            raise Exception('PlanarSurfaceSingle, create_updated_shell failed: %s' % str(e) )

    def create_position_vector(self, r):
        # project the vector onto the surface unit vectors to make sure that the coordinates are in the surface
        x, y = random_vector2D(r)
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def check(self):
        assert ProtectiveDomain.check(self)
        # check shell stuff
        assert isinstance(self.testShell, PlanarSurfaceSingletestShell)
        assert isinstance(self.shell.shape, Cylinder)

        max_radius, _, _ = self.testShell.get_max_dr_dzright_dzleft()
        assert self.shell.shape.radius <= max_radius
        if not self.is_reset():
            min_radius, _, _ = self.testShell.get_min_dr_dzright_dzleft()
            assert self.shell.shape.radius*SAFETY >= min_radius
        assert self.shell.shape.half_length == self.testShell.z_right(self.shell.shape.radius) 
        assert self.shell.shape.half_length == self.testShell.z_left (self.shell.shape.radius)
        if self.is_reset():
            assert self.shell.shape.radius == self.pid_particle_pair[1].radius

        assert self.is_reset() ^ (self.dt != 0.0)
        assert isinstance(self.structure, PlanarSurface)
        assert self.greens_function
        assert (self.shell.shape.unit_z == self.structure.shape.unit_z).all()

        assert feq(numpy.dot(self.pid_particle_pair[1].position - self.structure.shape.position,\
                    self.structure.shape.unit_z), 0.0, typical=self.shell.shape.half_length)

        return True

    def __str__(self):
        return 'PlanarSurface' + Single.__str__(self)


class CylindricalSurfaceSingle(NonInteractionSingle, hasCylindricalShell):
    """1 Particle inside a (cylindrical) shell on a CylindricalSurface. 
    (Rods).

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        assert isinstance(testShell, CylindricalSurfaceSingletestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def getv(self):
        return self.pid_particle_pair[1].v
    v = property(getv)

    def get_inner_a(self):
        # Here the free coordinate is not the radius of the domain as in the
        # PlanarSurfaceSingle and SphericalSingle, although it is treated as a
        # radius in NonInteractionSingle.draw_new_position for simplicity
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        inner_half_length = self.get_inner_a()
        return GreensFunction1DAbsAbs(self.D, self.v, 0.0, -inner_half_length, inner_half_length)

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    # For explanation see NonInteractionSingles and Others in shells.py
    def shell_list_for_single(self):
        # The shell should always be larger that the bare minimum for a test shell
        # Note that in case of a cylindrical shell this fits inside of the spherical multi shell
        min_half_length = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.half_length < min_half_length:
            # This should only happen if the domain is just initialized.
            position = self.shell.shape.position
            # There should always be space to make the spherical multi shell.
            radius = self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_half_length = self.pid_particle_pair[1].radius * math.sqrt(SINGLE_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.half_length < min_half_length:
            position = self.shell.shape.position
            radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
            # Keep all the other shells outside the burst radius (which is spherical!)
            # Note that this makes the cylindrical shell have a slightly smaller radius than the burst radius.
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            # When the shell is larger than the burst volume -> return the real shell.
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        # maybe update the reference_point of the shell before updating the shell
        try:
            dr, dz_right, dz_left = self.testShell.determine_possible_shell(self.structure.id, [self.domain_id], [])
            dz_right = min(dz_right, dz_left)       # make sure the domain is symmetric around the particle
            dz_left  = dz_right                     # This is not necessary but it's assumed later
            center, radius, half_length = self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                                             self.testShell.get_orientation_vector(),
                                                                             dr, dz_right, dz_left)
            return self.create_new_shell(center, radius, half_length, self.domain_id)
        except ShellmakingError as e:
            raise Exception('CylindricalSurfaceSingle, create_updated_shell failed: %s' % str(e) )

    def create_position_vector(self, z):
        # 'z' can be interpreted in two different ways here, it may a coordinate in the z direction or it may
        # be a displacement from the origin. In the first case it will already contain the drift information.
        if self.v == 0.0:
            # if there is no drift then we regard 'z' as a displacement from the origin. Although
            # it actually represents a coordinate when drawR was called. It is a displacement when z=a.
            # We randomize the direction.
            z = myrandom.choice(-1, 1) * z 

        elif self.v != 0.0 and feq(z, self.get_inner_a()):
            # When there is drift and z=a, the 'z' actually represent the displacement from the origin and a
            # boundary must be chosen.

            # The Escape can be either to the left or to the right.
            # The side of escape should be decided on by the flux through both boundaries at the escape time
            gf = self.greens_function()
            event_kind = draw_event_type_wrapper(gf, self.dt)
            if event_kind == PairEventKind.IV_REACTION:     # IV_REACTION -> ESCAPE through 'left' boundary
                z = -z
            elif event_kind == PairEventKind.IV_ESCAPE:     # IV_ESCAPE -> ESCAPE through 'right' boundary
                #z = z
                pass 
            else:
                raise NotImplemented()
        else:
            # When there was drift and the particle was not at the boundary.
            # -> In this case the 'z' actually signifies a coordinate and nothing has to be done.
            pass

        return z * self.structure.shape.unit_z

    def check(self):
        assert ProtectiveDomain.check(self)
        # check shell stuff
        assert isinstance(self.testShell, CylindricalSurfaceSingletestShell)
        assert isinstance(self.shell.shape, Cylinder)
#        assert self.shell.shape.radius <= self.testShell.get_max_radius()
#        if not self.is_reset():
#            assert self.shell.shape.radius >= self.testShell.get_min_radius()
        if self.is_reset():
            assert self.shell.shape.half_length == self.pid_particle_pair[1].radius

        assert self.is_reset() ^ (self.dt != 0.0)
        assert isinstance(self.structure, CylindricalSurface)
        assert self.greens_function
        assert (self.shell.shape.unit_z == self.structure.shape.unit_z).all()
        assert self.v < numpy.inf

        return True

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


class DiskSurfaceSingle(NonInteractionSingle, hasCylindricalShell):
    """1 Particle inside a (cylindrical) shell on a DiskSurface. 

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        assert isinstance(testShell, DiskSurfaceSingletestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def get_inner_a(self):
        return self.pid_particle_pair[1].radius

    def greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        # The Green's function should not be used for anything on the DiskSurface
        # but is defined for completeness here. # TODO Do we really have to do that?
        inner_half_length = self.get_inner_a()
        return GreensFunction1DAbsAbs(0.0, 0.0, 0.0, -inner_half_length, inner_half_length)

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    # For explanation see NonInteractionSingles and Others in shells.py
    def shell_list_for_single(self):
        # The shell should always be larger that the bare minimum for a test shell
        # Note that in case of a cylindrical shell this fits inside of the spherical multi shell
        min_half_length = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.half_length < min_half_length:
            # This should only happen if the domain is just initialized.
            position = self.shell.shape.position
            # There should always be space to make the spherical multi shell.
            radius = self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_half_length = self.pid_particle_pair[1].radius * math.sqrt(SINGLE_SHELL_FACTOR**2 - 1.0)
        if self.shell.shape.half_length < min_half_length:
            position = self.shell.shape.position
            radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
            # Keep all the other shells outside the burst radius (which is spherical!)
            # Note that this makes the cylindrical shell have a slightly smaller radius than the burst radius.
            fake_shell = SphericalShell(self.domain_id, Sphere(position, radius))
            return [(self.shell_id, fake_shell), ]
        else:
            # When the shell is larger than the burst volume -> return the real shell.
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        # maybe update the reference_point of the shell before updating the shell
        try:
            dr, dz_right, dz_left = self.testShell.determine_possible_shell(self.structure.id, [self.domain_id], [])
            dz_right = min(dz_right, dz_left)       # make sure the domain is symmetric around the particle
            dz_left  = dz_right                     # This is not necessary but it's assumed later
            center, radius, half_length = self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                                             self.testShell.get_orientation_vector(),
                                                                             dr, dz_right, dz_left)
            return self.create_new_shell(center, radius, half_length, self.domain_id)
        except ShellmakingError as e:
            raise Exception('DiskSurfaceSingle, create_updated_shell failed: %s' % str(e) )

    def create_position_vector(self, z):
        # Particles on DiskStructures are immobile; always keep them in the center of the structure
        return self.structure.shape.position

    def check(self):
        assert ProtectiveDomain.check(self)
        # check shell stuff
        assert isinstance(self.testShell, DiskSurfaceSingletestShell)
        assert isinstance(self.shell.shape, Cylinder)

        if self.is_reset():
            assert self.shell.shape.half_length == self.pid_particle_pair[1].radius

        assert self.is_reset() ^ (self.dt != 0.0)
        assert isinstance(self.structure, DiskSurface)
        #assert self.greens_function
        assert (self.shell.shape.unit_z == self.structure.shape.unit_z).all()

        return True

    def __str__(self):
        return 'DiskSurface' + Single.__str__(self)


class InteractionSingle(Single, hasCylindricalShell, Others):
    """Interactions singles are used when a particle is close to a 
    surface.

    There are 2 types of InteractionSingles:
        * PlanarSurfaceInteraction.
        * CylindricalSurfaceInteraction.

    The InteractionSingle holds three independent processes:
        * Diffusion in the direction that allows the interaction (radial or z)
        * Diffusion in the direction that is normal the the direction of interaction
        * Decay of the particle into a different or no particle
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):

        hasCylindricalShell.__init__(self, testShell, domain_id)
        Single.__init__(self, domain_id, shell_id, reactionrules)
        # Note: for the Others superclass nothing is to be initialized.

        self.interactionrules = interactionrules        # the surface interaction 'reactions'
        self.intrule = None
        self.interaction_ktot = self.calc_ktot(interactionrules)

        self.origin_structure = self.testShell.origin_structure
        self.target_structure = self.testShell.target_structure 
             # The surface with which the particle is trying to interact                         

    def get_interaction_rule(self):
        if self.intrule == None:
            self.intrule = self.draw_reaction_rule(self.interactionrules)
        return self.intrule
    interactionrule = property(get_interaction_rule)

        # There are defined here but are overloaded later
        # These ProtectiveDomains have a Cylindrical shell that is scaled in r and z
        # Below methods go together with this property
    def get_inner_dz_left(self):
        pass

    def get_inner_dz_right(self):
        pass

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple(),
                   self.draw_iv_event_time_tuple())     # above two events also occur in NonInteractingSingle
                                                        # IV event is added in InteractionSingle

    def draw_iv_event_time_tuple(self):
        """Return an (interaction/escape time, event type)-tuple.
        Note that an interaction is modelled similarly to an escape.
 
        Note: we are not calling single.draw_iv_event_type() just yet, but 
        postpone it to the very last minute (when this event is executed 
        in fire_single). So IV_EVENT can still be an iv escape or an iv 
        interaction.

        """
        dt = draw_time_wrapper(self.iv_greens_function())

        return dt, EventType.IV_EVENT

    def draw_iv_event_type(self):
        gf = self.iv_greens_function()
        event_kind = draw_event_type_wrapper(gf, self.dt)
        if event_kind == PairEventKind.IV_REACTION:     # Note REACTION -> INTERACTION
            return EventType.IV_INTERACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def __str__(self):
        pass

class PlanarSurfaceInteraction(InteractionSingle):
    """1 Particle close to a PlanarSurface, inside a cylindrical shell 
    placed against the surface with its flat face.

        * Particle coordinates inside shell: r, theta, z.
        * Coordinates: radial r, cartesian z.
        * Initial position: r = 0, z = z0 (relative to center shell).
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):
        # shell_orientation_vector and z0 are supplied by arguments because Domains can't do
        # spatial calculations (they don't know about the world)
        # NOTE shell_orientation_vector should be normalized and pointing from the surface to the particle!

        assert isinstance(testShell, PlanarSurfaceInteractiontestShell)
        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        # z0 is the initial position of the particle relative to the center
        # of the shell in the direction normal to the surface (so z0 is usually negative).
        self.z0 = -self.get_inner_dz_left() + self.testShell.particle_surface_distance

    def get_inner_dz_right(self):
        # This is the distance from the center of the shell
        # to the flat boundary away from the surface
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_left(self):
        # This is the distance from the center of the shell to the surface
        # Note that it is implied that the shell sticks out at the other end
        # of the planar surface by the particle radius, and that the surface
        # has no thickness.
        dist_to_surface = self.shell.shape.half_length - self.pid_particle_pair[1].radius
        return dist_to_surface #- self.pid_particle_pair[1].radius

    def greens_function(self):
        # This is the Greens function that does not produce REACTION event
        #return GreensFunction3DAbsSYm(self.getD(), get_inner_a())
        return GreensFunction2DAbsSym(self.D, self.get_inner_a())

    def iv_greens_function(self):
        # The green's function that also modelles the association of the particle
        # with the planar surface.
        return GreensFunction1DRadAbs(self.D, self.interaction_ktot, self.z0,
                                      -self.get_inner_dz_left(), self.get_inner_dz_right())

    def draw_new_position(self, dt, event_type):
        """Draws a new position for the particle for a given event type at a given time
        Note that it also 'moves the particle into the surface' when an interaction
        takes place.
        """
        oldpos = self.pid_particle_pair[1].position
        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos             # just return oldpos, particle will be removed anyway
        else:
            # calculate r vector
            a = self.get_inner_a()
            if event_type == EventType.SINGLE_ESCAPE:
                r = a
            else:
                gf = self.greens_function()
                r = draw_r_wrapper(gf, dt, a)

            x, y = random_vector2D(r)
            # express the r vector in the unit vectors of the surface to make sure the particle is
            # parallel to the surface (due to numerical problem)
            vector_r = x * self.target_structure.shape.unit_x + y * self.target_structure.shape.unit_y

            # calculate z vector
            z_surface = -self.get_inner_dz_left()       # left side of the inner domain
            z_not_surface = self.get_inner_dz_right()   # right side of the inner domain (away from the surface)

            # If the event was not yet fully specified
            if event_type == EventType.IV_EVENT:        
                self.event_type = self.draw_iv_event_type()
                event_type = self.event_type


            if event_type == EventType.IV_INTERACTION:
                z = z_surface 
            elif event_type == EventType.IV_ESCAPE:
                z = z_not_surface
            else:
                gf_iv = self.iv_greens_function()
                z = draw_r_wrapper(gf_iv, dt, z_not_surface, z_surface) # last two args= a, sigma

            # express the vector_z in the unit vectors of the surface to prevent the particle from
            # leaving the surface due to numerical problem
            vector_z = z * self.shell.shape.unit_z

            # The new position is relative to the center of the shell
            # note that we need to make sure that the shell.shape.position and vector_r and vector_z
            # are correct (in the structure) to prevent the particle leaving the structure due to numerical errors
            newpos = self.shell.shape.position + vector_r + vector_z

        # The structure on which the particle ended is always the same as on which it began.
        structure_id = self.origin_structure.id

        return newpos, structure_id

    def __str__(self):
        return ('PlanarSurfaceInteraction' + Single.__str__(self) + \
               'radius = %s, half_length = %s' %
                (self.shell.shape.radius, self.shell.shape.half_length))


class CylindricalSurfaceInteraction(InteractionSingle):
    """1 Particle close to a CylindricalSurface, inside a cylindrical shell 
    that surrounds the surface.

        * Particle coordinates inside shell: r, theta, z.
        * Domains: composite r-theta, cartesian z.
        * Initial position: r = r, theta = 0, z = z.
        * Selected randomly when drawing displacement vector: none.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):

        assert isinstance(testShell, CylindricalSurfaceInteractiontestShell)
        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        # unit_r is the vector from the cylindrical surface to the particle
        # r0 is the distance from the center of the cylinder to the center of the particle
        # z0 is known to be zero (the particle being in the center of the shell in the z direction)
        self.unit_r = self.testShell.reference_vector    
        self.r0     = self.testShell.particle_surface_distance + self.target_structure.shape.radius
        self.z0     = 0.0

    def get_inner_dz_left(self):
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_right(self):
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_sigma(self):
        return self.target_structure.shape.radius + self.pid_particle_pair[1].radius

    def greens_function(self):
        # The greens function not used for the interaction but for the other coordinate
        # Free diffusion in z direction, drift is zero. Note that also z0 = 0.0
        return GreensFunction1DAbsAbs(self.D, 0.0, self.z0,
                                      -self.get_inner_dz_left(),
                                      self.get_inner_dz_right())

    def iv_greens_function(self):
        # Green's function used for the interaction
        # Interaction possible in r direction.
        return GreensFunction2DRadAbs(self.D, self.interaction_ktot, self.r0, self.get_inner_sigma(),
                                      self.get_inner_a())

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos
        else:
            # 1) Draw z.
            if event_type == EventType.SINGLE_ESCAPE:
                # Moving this checks to the Green's functions is not a good 
                # idea, because then you'd draw an unused random number.  
                # The same yields for the draw_new_com and draw_new_iv.  
                # Assume for now that z0 = 0 -> both boundaries are equally likely
                z = myrandom.choice(-1, 1) * self.get_inner_dz_left()
            else:
                gf = self.greens_function()
                z = draw_r_wrapper(gf, dt, self.get_inner_dz_left())

            # Direction matters. We assume that the direction of the shell is equal to the direction of
            # the surface.
            z_vector = z * self.shell.shape.unit_z


            # 2) Draw r and theta.
            # If the event was not yet fully specified
            if event_type == EventType.IV_EVENT:        
                self.event_type = self.draw_iv_event_type()
                event_type = self.event_type

            sigma = self.get_inner_sigma()
            a_r = self.get_inner_a()
            gf_iv = self.iv_greens_function()

            if event_type == EventType.IV_ESCAPE:
                r = a_r
            elif event_type == EventType.IV_INTERACTION:
                r = sigma
            else:
                r = draw_r_wrapper(gf_iv, dt, a_r, sigma)
            theta = draw_theta_wrapper(gf_iv, r, dt)

            r_vector = r * rotate_vector(self.unit_r, self.shell.shape.unit_z,
                                     theta)

            # Add displacement to shell.shape.position, not to particle.position.  
            newpos = self.shell.shape.position + z_vector + r_vector

        # The structure on which the particle ended is always the same as on which it began.
        structure_id = self.origin_structure.id

        return newpos, structure_id

    def __str__(self):
        return ('CylindricalSurfaceInteraction' + Single.__str__(self) + \
               'radius = %s, half_length = %s' %
                (self.shell.shape.radius, self.shell.shape.half_length))


class CylindricalSurfaceSink(InteractionSingle):
    """1 Particle inside a (Cylindrical) shell on a CylindricalSurface.
        Inside the domain is a sink.

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):

        assert isinstance(testShell, CylindricalSurfaceSinktestShell)
        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        # zsink is the position of the sink relative to the center of the shell.
        # Note that we assume that the particle is at the center of the shell, and the unit_z of the shell
        # points from the sink to the particle.
        self.zsink = -self.testShell.particle_surface_distance

    def get_inner_a(self):
        # This is the distance that the particle can travel from the center to one of the sides.
        # NOTE: it is implied that the particle starts in the center of the shell.
        # NOTE: although abstractly defined in InteractionSingle inner_dz_right/left are not
        #       specified here.
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def determine_next_event(self):
        """Return an (event_time, event_type)-tuple.

        """
        return min(self.draw_reaction_time_tuple(),
                   self.draw_iv_event_time_tuple())

    def iv_greens_function(self):
        # The Green's function used for both ESCAPE and IV (sink) events.
        # Particle allways starts in the middle for now.
        inner_half_length = self.get_inner_a()

        return GreensFunction1DAbsSinkAbs(self.D, self.interaction_ktot, 0.0, self.zsink,
                                          -inner_half_length, inner_half_length)

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0 or \
                (event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0):
            newpos = oldpos
        else:
            gf = self.iv_greens_function()

            # If an IV event took place, and the event was not yet fully specified
            if event_type == EventType.IV_EVENT:
                self.event_type = self.draw_iv_event_type()
                event_type = self.event_type

            if event_type == EventType.IV_ESCAPE:
                # The particle left the domain through either boundary

                # Determine whether the particle escaped through the left (at s) or right (at a) boundary.
                flux_a = abs( gf.flux_leavea( dt ) )
                flux_s = abs( gf.flux_leaves( dt ) )
                rnd = myrandom.uniform() * ( flux_a + flux_s )

                inner_half_length = self.get_inner_a()

                if( rnd < flux_a ):
                    # The 'right' boundary (boundary to which the shell.shape.unit_z point)
                    z = inner_half_length
                else:
                    # The 'left' boundary
                    z = -inner_half_length

            elif event_type == EventType.IV_INTERACTION:
                # The particle left the domain through the sink
                z = self.zsink
            else:
                # Some other event took place and the particle didn't escape
                z = draw_r_wrapper(gf, dt, self.get_inner_a(), -self.get_inner_a())

            # Add displacement to shell.shape.position, not to particle.position.
            # The direction of the shell is leading (not direction of the structure)
            z_vector = self.shell.shape.unit_z * z
            newpos = self.shell.shape.position + z_vector

        # The structure on which the particle ended is always the same as on which it began.
        structure_id = self.origin_structure.id

        return newpos, structure_id

    def __str__(self):
        return 'CylindricalSurfaceSink ' + Single.__str__(self)


class CylindricalSurfaceCapInteraction(InteractionSingle):
    """1 Particle inside a (Cylindrical) shell on a CylindricalSurface
       limited by a cap. The cap is a reactive surface to the particle
       and defines the exit point from the cylinder.

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none. # TODO What does this mean?
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):

        assert isinstance(testShell, CylindricalSurfaceCapInteractiontestShell)
        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        self.zcap = self.testShell.get_referencepoint()

    def getv(self):
        return self.pid_particle_pair[1].v
    v = property(getv)

    def get_inner_dz_left(self):
        # This is the distance that the particle can travel from its initial position
        # towards the reactive cap.
        # For the cap it is the center of particle that "reacts" with it.
        # The testShell is taking this into account at construction to ensure that
        # there is enough space around the cap in case of reaction.
        # Note that the reference point of the shell is the cap position.
        return self.testShell.particle_surface_distance

    def get_inner_dz_right(self):
        # This is the distance that the particle can travel from its initial position
        # towards the absorbing cylinder boundary opposite of the cap.
        # Note that the reference point of the shell is the cap position,
        # therefore we have to subtract particle_surface_distance here.
        return self.testShell.dz_right - self.testShell.particle_surface_distance \
                                       - self.pid_particle_pair[1].radius
        # TODO is shell.dz_right always positive?

    def determine_next_event(self):
        """Return an (event_time, event_type)-tuple.

        """
        return min(self.draw_reaction_time_tuple(),
                   self.draw_iv_event_time_tuple())

    def iv_greens_function(self):

        # zcap = position of the cap relative to the initial particle position.
        # zabs = position of the absorbing cylinder boundary opposite of the cap
        dz_cap = -self.get_inner_dz_left()
        dz_abs = self.get_inner_dz_right()

        return GreensFunction1DRadAbs(self.D, self.v, self.interaction_ktot, 0.0, dz_cap, dz_abs)

    def draw_new_position(self, dt, event_type):

        oldpos = self.pid_particle_pair[1].position

        if self.D == 0 or \
                (event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0):
            newpos = oldpos
        else:
            gf = self.iv_greens_function()

            if event_type == EventType.IV_INTERACTION:
                # The particle left the domain through the reactive boundary, i.e. bound to the cap.
                dz = -self.get_inner_dz_left()
                structure_id = self.target_structure.id

            elif event_type == EventType.IV_ESCAPE:
                # The particle left the domain through the absorbing boundary,
                # i.e. opposite of the cap.
                dz = self.get_inner_dz_right()
                structure_id = self.origin_structure.id

            else:
                # Some other event took place and the particle didn't escape;
                # as a consequence it stays on the origin structure.
                # NOTE: draw_r_wrapper(): attention, order of arguments inverted!
                # Last two arguments are a, sigma (in this order!)
                dz = draw_r_wrapper(gf, dt, self.get_inner_dz_right(), -self.get_inner_dz_left())
                structure_id = self.origin_structure.id

            # Add displacement to particle.position.
            # The direction of the shell is leading (not direction of the structure)
            # The convention is that shell.shape.unit_z points from the cap towards
            # the particle's initial position.
            dz_vector = self.shell.shape.unit_z * dz
            newpos = oldpos + dz_vector

        return newpos, structure_id

    def __str__(self):
        return 'CylindricalSurfaceCapInteraction ' + Single.__str__(self)


class TransitionSingle(Single, TransitionTools, Others):
    """TransitionSingles are used when a particle changes from
       one structure to another. The difference to an InteractionSingle
       is that the transition process is instantaneous, i.e. the 
       transition propensity is infinitely high.    

       TransitionSingles feature:
        * a function process_new_position_vector which checks whether a newly
          drawn position is still on the original structure of the particle.
          If this is the case it will not be changed.
          If it is on the new structure then the new position will be transformed
          accordingly and the structure_id of the particle will be updated.
        * several properties that are needed for the above transform.
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        Single.__init__(self, domain_id, shell_id, reactionrules)
        TransitionTools.__init__(self, testShell, self.pid_particle_pair[1].position)
        # Note: for the Others superclass nothing is to be initialized.

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple())     # above two events also occur in NonInteractingSingle

    def __str__(self):
        pass


class PlanarSurfaceTransitionSingle(TransitionSingle, hasSphericalShell):
    """1 Particle inside a (spherical) shell at the edge of two planar surfaces

        * Particle coordinates on surface: x, y.
        * Domain: radial r. (determines x and y together with theta).
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        assert isinstance(testShell, PlanarSurfaceTransitionSingletestShell)
        hasSphericalShell.__init__(self, testShell, domain_id)
        TransitionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules)


    # The same Greens function is used as for the normal PlanarSurfaceSingle;
    # If the position is off the surface of origin, it will be transformed towards the target surface
    def greens_function(self):
        return GreensFunction2DAbsSym(self.D, self.get_inner_a())

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0:
            newpos, new_structure_id = oldpos, self.origin_structure.id
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos, new_structure_id = oldpos, self.origin_structure.id
        else:
            # Calculate r
            if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  

                r = self.get_inner_a()
            else:
                gf = self.greens_function()
                r = draw_r_wrapper(gf, dt, self.get_inner_a())
                # Note that in case of 1D diffusion r has a direction. It is
                # the lateral displacement and can be positive or negative. 
                # In other cases r is always positive and denotes
                # radial displacement from the center.

            # Calculate the displacement from oldpos in the surface of origin
            displacement = self.create_position_vector(r)

            # Now check whether the new position is in the surface of origin;
            # if not, transform it to the target surface (function inherited from TransitionTools)
            newpos, new_structure_id = self.process_new_position_vector(oldpos, displacement)

        return newpos, new_structure_id

    def create_position_vector(self, r):
        # project the vector onto the surface unit vectors to make sure
        # that the coordinates are in the surface
        x, y = random_vector2D(r)
        return x * self.origin_structure.shape.unit_x + y * self.origin_structure.shape.unit_y

    def __str__(self):
        return 'PlanarSurfaceTransition' + Single.__str__(self)


class DummySingle(object):
    def __init__(self):
        self.multiplicity = 1
        self.half_length = 0.0
        self.shellList = [Sphere(NOWHERE, 0.0)]

    def getMinSize(self):
        return 0.0

    def getD(self):
        return 0.0

    def getPos(self):
        return NOWHERE
    pos = property(getPos)

    def __str__(self):
        return 'DummySingle()'

