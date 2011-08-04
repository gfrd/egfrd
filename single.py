from _gfrd import *
from _greens_functions import *
from greens_function_wrapper import *
from constants import EventType
import utils

__all__ = [
    'CylindricalSurfaceSingle',
    'PlanarSurfaceSingle',
    'SphericalSingle',
    'Single',
    'InteractionSingle',
    'CylindricalSurfaceInteraction',
    'PlanarSurfaceInteraction',
    ]

class Single(object):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see 
    coordinate.py. For each coordinate the Green's function is 
    specified.

    """
    def __init__(self, domain_id, pid_particle_pair, reactiontypes, surface):
        self.multiplicity = 1
        self.num_shells = 1

        self.pid_particle_pair = pid_particle_pair
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.last_time = 0.0
        self.dt = 0.0
        self.event_type = None

        self.surface = surface

        self.event_id = None

        self.domain_id = domain_id

        self.updatek_tot()

    def getD(self):
        return self.pid_particle_pair[1].D
    D = property(getD)

    def getv(self):
        return 0 # TODO.
        return self.pid_particle_pair[1].v
    v = property(getv)

    def get_shell_id(self):
        return self.shell_list[0][0]
    shell_id = property(get_shell_id)

    def get_shell(self):
        return self.shell_list[0][1]
    shell = property(get_shell)

    def get_shell_id_shell_pair(self):
        return self.shell_list[0]
    def set_shell_id_shell_pair(self, value):
        self.shell_list[0] = value
    shell_id_shell_pair = property(get_shell_id_shell_pair, 
                                   set_shell_id_shell_pair)

    def get_mobility_radius(self):
        return self.get_shell_size() - self.pid_particle_pair[1].radius

    def get_shell_size(self):
        return self.shell_list[0][1].shape.radius

    def initialize(self, t):
        '''
        Reset the Single.

        Radius (shell size) is shrunken to the actual radius of the 
        particle.  self.dt is reset to 0.0.  Do not forget to reschedule 
        this Single after calling this method.
        '''
        self.dt = 0.0
        self.last_time = t
        self.event_type = EventType.SINGLE_ESCAPE

    def is_reset(self):
        return self.dt == 0.0 and self.event_type == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple(self):
        """Return a (reaction time, event type)-tuple.

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

        In case of an interaction single, does not handle escapes in 
        the interaction coordinate.
        
        """
        if self.getD() == 0:
            dt = numpy.inf
        else:
            dt = draw_time_wrapper(self.greens_function())

        event_type = EventType.SINGLE_ESCAPE
        return dt, event_type

    def updatek_tot(self):
        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k

    def draw_reaction_rule(self):
        k_array = [rt.k for rt in self.reactiontypes]
        k_array = numpy.add.accumulate(k_array)
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted(k_array, rnd * k_max)

        return self.reactiontypes[i]

    def check(self):
        pass

    def __str__(self):
        pid = self.pid_particle_pair[0]
        sid = self.pid_particle_pair[1].sid
        name = self.world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        return 'Single[%s: %s, ST%s, %s]' % (self.domain_id, pid, name,
                                             self.pid_particle_pair[1].position)


class NonInteractionSingle(Single):
    """1 Particle inside a shell, no other particles around. 

    There are 3 types of NonInteractionSingles:
        * SphericalSingle: spherical shell, 3D movement.
        * PlanarSurfaceSingle: cylindrical shell, 2D movement.
        * CylindricalSurfaceSingle: cylindrical shell, 1D movement.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        Single.__init__(self, domain_id, pid_particle_pair, reactiontypes,
                        surface)

        # Create shell.
        shell = self.create_new_shell(pid_particle_pair[1].position,
                                      pid_particle_pair[1].radius, domain_id)

        self.shell_list = [(shell_id, shell), ]

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple())

    def draw_new_position(self, dt, event_type):
        if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  
            r = self.get_mobility_radius()
        else:
            gf = self.greens_function()
            r = draw_r_wrapper(gf, dt, self.get_mobility_radius())

        displacement = self.create_position_vector(r)

        # This should be checked in the unit test of random_vector.
        # if __debug__:
        #     scale = self.pid_particle_pair[1].radius
        #     if feq(length(displacement), abs(r), typical=scale) == False:
        #         raise AssertionError('displacement != abs(r): %g != %g.' % 
        #                              (length(displacement), abs(r)))

        # Add displacement to shape.position, not to particle.position.  
        # This distinction is important only in the case of an 
        # asymmetric 1D domain (r0 != 0, or drift), since draw_r always 
        # returns a position relative to the centre of the shell (r=0), 
        # not relative to r0.
        return self.shell.shape.position + displacement


class SphericalSingle(NonInteractionSingle):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector:
          theta, phi.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_mobility_radius())

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    def create_position_vector(self, r):
        return random_vector(r)

    def __str__(self):
        return 'Spherical' + Single.__str__(self)


class PlanarSurfaceSingle(NonInteractionSingle):
    """1 Particle inside a (cylindrical) shell on a PlanarSurface. (Hockey 
    pucks).

        * Particle coordinates on surface: x, y.
        * Domain: radial r. (determines x and y together with theta).
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        return GreensFunction2DAbsSym(self.getD(),
                                          self.get_mobility_radius())

    def create_new_shell(self, position, radius, domain_id):
        # The half_length (thickness) of a hockey puck is not more than 
        # it has to be (namely the radius of the particle), so if the 
        # particle undergoes an unbinding reaction we still have to 
        # clear the target volume and the move may be rejected (NoSpace 
        # error).
        orientation = normalize(
            utils.crossproduct(self.surface.shape.unit_x,
                               self.surface.shape.unit_y))
        half_length = self.pid_particle_pair[1].radius
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def __str__(self):
        return 'PlanarSurface' + Single.__str__(self)


class CylindricalSurfaceSingle(NonInteractionSingle):
    """1 Particle inside a (cylindrical) shell on a CylindricalSurface. 
    (Rods).

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.getD(), self.getv(), 0.0, -self.get_mobility_radius(), self.get_mobility_radius())

    def create_new_shell(self, position, half_length, domain_id):
        # The radius of a rod is not more than it has to be (namely the 
        # radius of the particle), so if the particle undergoes an 
        # unbinding reaction we still have to clear the target volume 
        # and the move may be rejected (NoSpace error).
        radius = self.pid_particle_pair[1].radius
        orientation = self.surface.shape.unit_z
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, z):
        if utils.feq(z, self.get_mobility_radius()):
            # Escape, can be either to the left or to the right.
            z = myrandom.choice(-1, 1) * z 
        return z * self.shell_list[0][1].shape.unit_z

    def get_shell_size(self):
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell_list[0][1].shape.half_length

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


class InteractionSingle(Single):
    """Interactions singles are used when a particle is close to a 
    surface.

    There are 2 types of InteractionSingles:
        * PlanarSurfaceInteraction.
        * CylindricalSurfaceInteraction.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, shell,  
                 reactiontypes, interaction_type, surface):

        Single.__init__(self, domain_id, pid_particle_pair, reactiontypes,
                        surface)
        self.interaction_type = interaction_type
        self.shell_list = [(shell_id, shell), ]

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_iv_interaction_or_escape_time_tuple(),
                   self.draw_reaction_time_tuple())

    def draw_iv_interaction_or_escape_time_tuple(self):
        """Return an (interaction/escape time, event type)-tuple.
 
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
        if event_kind == PairEventKind.IV_REACTION:
            return EventType.IV_REACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()


class PlanarSurfaceInteraction(InteractionSingle):
    """1 Particle close to a PlanarSurface, inside a cylindrical shell 
    placed on top of the surface.

        * Particle coordinates inside shell: r, theta, z.
        * Coordinates: radial r, cartesian z.
        * Initial position: r = 0, z = z.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes,
                 surface, interaction_type, origin, orientation, half_length, 
                 particle_offset, projected_point, size_of_domain):
        InteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                   shell_id, reactiontypes, surface, 
                                   interaction_type, origin, orientation,
                                   half_length, particle_offset,
                                   projected_point, size_of_domain)

        # Compute from dr, dz_left and dz_right the length of the new 
        # domain (only that part of the cylinder where the particle can 
        # diffuse to) and the radius and half_length of the new 
        # cylinder.

        # a
        # sigma
        sigma = self.surface.shape.Lz + self.pid_particle_pair[1].radius
        # r0
        if isinstance(surface, PlanarSurface):
            # Heads up: length of domain (the part of the cylinder 
            # where the particle can get) is a bit different from the 
            # length of the cylinder because of dz_left.
            # Assume surface.Lz == 0
            length_of_domain = particle_distance + dz_right
            radius = dr
            half_length = (particle_distance + dz_left + dz_right) / 2

        # TODO. This is not correct anymore.
        if isinstance(surface, PlanarSurface):
            # Compute new particle offset relative to origin of the 
            # domain in the z-direction (z=0).
            # The surface boundary is at z=-length_of_domain / 2, the 
            # other end in the middles*at* the surface boundary. (z=L) 
            # is at the end of the cylinder.
            # (min_radius correction in the constructor).
            particle_offset = [0, particle_distance - length_of_domain.Lz]

    def greens_function_r(self):
        # Todo. 2D gf Abs Sym.
        #gf = FirstPassageGreensFunction2D(self.getD())
        return FirstPassageGreensFunction(self.getD(),
                                          self.get_mobility_radius())

    def greens_function_z(self):
        # Todo. 1D gf Rad Abs should be sigma to a.
        #gf = FirstPassageGreensFunction1DRad(self.D_tot, self.rt.k)

        return FirstPassagePairGreensFunction(self.getD(),
                                              self.interaction_type.k, self.r0,
                                              self.sigma, self.sigma)
        gf.seta(a)
        return gf

        # TODO.

        # Cartesian domain of half_length L. Correction with getMinRadius().
        #zDomain = CartesianDomain(particle_offset[1] - self.getMinRadius(), 
        #                           size_of_domain - 2 * self.getMinRadius(), 
        #                           gfz)

    def createNewShell(self, position, radius, domain_id):
        orientation = self.orientation
        half_length = self.half_length
        return CylindricalShell(position, radius, orientation, half_length,
                                domain_id)

    def draw_new_positions(self, dt, eventType):
        gf_r = self.greens_function_r()
        a = self.get_mobility_radius()
        r = draw_displacement_wrapper(gf_r, dt, eventType, a)
        x, y = randomVector2D(r)

        # Todo. Cartesian coordinate will return absolute position.
        gf_z = self.greens_function_z()
        # Heads up. size_of_domain is different from the half_length of 
        # the cylinder, because the domain ends where the surface 
        # starts, while the cylinder is continuing into the surface.
        a = 1 #TODO 
        r0 = self.particle_offset
        sigma = self.sigma
        z = draw_displacement_wrapper(gf_z, dt, eventType, a, r0, sigma)

        # Add displacement vector to self.particle.pos, not self.pos.
        return self.particle.pos + x * self.interactionSurface.unitX + \
               y * self.interactionSurface.unitY + z * self.shellList[0].unitZ

    def __str__(self):
        return 'PlanarSurfaceInteraction' + Single.__str__(self)


class CylindricalSurfaceInteraction(InteractionSingle):
    """1 Particle close to a CylindricalSurface, inside a cylindrical shell 
    that surrounds the surface.

        * Particle coordinates inside shell: r, theta, z.
        * Domains: composite r-theta, cartesian z.
        * Initial position: r = r, theta = 0, z = z.
        * Selected randomly when drawing displacement vector: none.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes,
                 interaction_type, surface, projected_point, particle_distance,
                 orientation_vector, dr, dz_left, dz_right):

        self.particle_distance = particle_distance
        self.dz_left = dz_left

        # Only needed for this type of Interaction.
        self.unit_r = normalize(pid_particle_pair[1].position - projected_point)

        shell = self.create_shell(domain_id, projected_point,
                                  particle_distance, orientation_vector,
                                  dr, dz_left, dz_right)

        InteractionSingle.__init__(self, domain_id, pid_particle_pair,
                                   shell_id, shell, reactiontypes,
                                   interaction_type, surface)

    def create_shell(self, domain_id, projected_point, particle_distance,
                     orientation_vector, dr, dz_left, dz_right):

        # Compute origin, radius and half_length of cylinder.
        # This stuff is complicated.
        half_length = (dz_left + dz_right) / 2
        shiftZ = half_length - dz_left
        origin = projected_point + shiftZ * orientation_vector
        radius = dr + particle_distance

        # Return shell.
        return CylindricalShell(domain_id, Cylinder(origin, radius, 
                                                    orientation_vector, 
                                                    half_length))

    def greens_function(self):
        z0 = self.dz_left - self.shell.shape.half_length

        # Free diffusion in z direction.
        return GreensFunction1DAbsAbs(self.getD(), self.getv(), z0,
                                      -self.get_mobility_radius(),
                                      self.get_mobility_radius())

    def iv_greens_function(self):
        # Interaction possible in r direction.
        # Todo.
        #k = self.interaction_type.k
        k = 0
        r0 = self.particle_distance
        sigma = self.surface.shape.radius + self.pid_particle_pair[1].radius
        a_r = self.shell.shape.radius - self.pid_particle_pair[1].radius
        # Todo 2D.
        return GreensFunction3DRadAbs(self.getD(), k, r0, sigma, a_r)

    def draw_new_position(self, dt, event_type):
        # 1) Draw z.
        if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  
            z = self.get_mobility_radius()
        else:
            gf = self.greens_function()
            z = draw_r_wrapper(gf, dt, self.get_mobility_radius())

        # Orientation matters, so use shell.shape.unit_z instead of 
        # surface.shape.unit_z.
        z_vector = z * self.shell.shape.unit_z

        # 2) Draw r and theta.
        sigma = self.surface.shape.radius + self.pid_particle_pair[1].radius
        a_r = self.shell.shape.radius - self.pid_particle_pair[1].radius
        iv_gf = self.iv_greens_function()
        if event_type == EventType.IV_ESCAPE:
            r = a_r
        elif event_type == EventType.IV_INTERACTION:
            r = sigma
        else:
            r = draw_r_wrapper(iv_gf, dt, a_r, sigma)
        theta = draw_theta_wrapper(iv_gf, r, dt)

        r_vector = r * rotate_vector(self.unit_r, self.shell.shape.unit_z,
                                     theta)

        # Add displacement to shape.position, not to particle.position.  
        return self.shell.shape.position + z_vector + r_vector

    def get_shell_size(self):
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell_list[0][1].shape.half_length

    def __str__(self):
        return 'CylindricalSurfaceInteraction' + Single.__str__(self)

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

