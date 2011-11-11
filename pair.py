from _gfrd import *
from constants import EventType
from _greens_functions import *
from greens_function_wrapper import *

from domain import (
    Domain,
    ProtectiveDomain)
from single import (
    NonInteractionSingle,
    InteractionSingle)
from multi import (
    Multi)

__all__ = [
    'CylindricalSurfacePair',
    'PlanarSurfacePair',
    'SphericalPair',
    'Pair',
    'SimplePair',
    'MixedPair',
    ]

if __debug__:
    PRECISION = 3
    FORMAT_DOUBLE = '%.' + str(PRECISION) + 'g'


class Pair(ProtectiveDomain):
    """There are 3 types of pairs:
        * SphericalPair
        * PlanarSurfacePair
        * CylindricalSurfacePair

    """
    # CUTOFF_FACTOR is a threshold to choose between the real and 
    # approximate Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, domain_id, single1, single2, shell_id, r0, 
                 rrs):
        ProtectiveDomain.__init__(self, domain_id)

        self.multiplicity = 2

        self.single1 = single1
        self.single2 = single2 
        self.pid_particle_pair1 = single1.pid_particle_pair
        self.pid_particle_pair2 = single2.pid_particle_pair

        self.interparticle_rrs = rrs
        self.interparticle_ktot = self.calc_ktot(rrs)

        self.shell_id = shell_id

    def __del__(self):
        if __debug__:
            log.debug('del %s' % str(self))

    def get_D_tot(self):
        return self.pid_particle_pair1[1].D + \
               self.pid_particle_pair2[1].D
    D_tot = property(get_D_tot)
    D_r   = property(get_D_tot)

    def get_D_R(self):
        return (self.pid_particle_pair1[1].D *
                self.pid_particle_pair2[1].D) / self.D_tot
    D_R = property(get_D_R)

    def initialize(self, t):
        # re-initialize the time parameters of the Domain to reuse it.
        # THIS IS PROBABLY USELESS
        self.last_time = t
        self.dt = 0
        self.event_type = None

    # draws a first event time and the (not fully specified) type of event
    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(r0))  # uses the full solution to draw IV first event time
        if dt_com < dt_iv:
            return dt_com, EventType.COM_ESCAPE, None
        else:
            # Note: we are not calling pair.draw_iv_event_type yet, but 
            # postpone it to the very last minute (when this event is 
            # executed in fire_pair). So IV_EVENT can still be an iv 
            # escape or an iv reaction.
            return dt_iv, EventType.IV_EVENT, None

    def draw_single_reaction_time_tuple(self):
        """Return a (reaction time, event type, reactingsingle)-tuple.

        """
        dt_reaction1, event_type1 = self.single1.draw_reaction_time_tuple()
        dt_reaction2, event_type2 = self.single2.draw_reaction_time_tuple()
        if dt_reaction1 < dt_reaction2:
            return dt_reaction1, event_type1, self.single1
        else:
            return dt_reaction2, event_type2, self.single2

    def determine_next_event(self, r0):
        """Return a (event time, event type, reactingsingle)-tuple.

        """
        return min(self.draw_com_escape_or_iv_event_time_tuple(r0), 
                   self.draw_single_reaction_time_tuple()) 

    def draw_iv_event_type(self, r0):
        gf = self.iv_greens_function(r0)
        event_kind = draw_event_type_wrapper(gf, self.dt)
        if event_kind == PairEventKind.IV_REACTION:
            return EventType.IV_REACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def draw_new_positions(self, dt, r0, old_iv, event_type):
        """Calculate new positions of the pair particles using a new 
        center-of-mass, a new inter-particle vector, and an old 
        inter-particle vector.

        """
#        # If the two particles have reacted this returns the
#        # new position twice
#        # The positions are relative to the center of the shell
        new_com = self.draw_new_com(dt, event_type)

#        if event_type == EventType.IV_REACTION:
#            newpos1 = new_com
#            newpos2 = new_com
#        else:
        new_iv = self.draw_new_iv(dt, r0, old_iv, event_type)
        newpos1, newpos2 = self.do_back_transform(new_com, new_iv)

        return newpos1, newpos2

    def draw_new_com(self, dt, event_type):
        # draws a new coordinate for the CoM in world coordinates
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        return self.com + self.create_com_vector(r)

        # Add displacement to old CoM. This assumes (correctly) that 
        # r0=0 for the CoM. Compare this to 1D singles, where r0 is not  
        # necesseraly 0.
#        return displacement

    def draw_new_iv(self, dt, r0, old_iv, event_type):
        gf = self.choose_pair_greens_function(r0, dt)
        if event_type == EventType.IV_ESCAPE:
            r = self.a_r
        elif event_type == EventType.IV_REACTION:
            r = self.sigma      # maybe this should be zero (product particle is at CoM)
        else:
            r = draw_r_wrapper(gf, dt, self.a_r, self.sigma)

        return self.create_interparticle_vector(gf, r, dt, r0, old_iv)

    def check(self):
        pass

    def __str__(self):
        sid1 = self.pid_particle_pair1[1].sid
        sid2 = self.pid_particle_pair2[1].sid
        name1 = self.world.model.get_species_type_by_id(sid1)["name"]
        name2 = self.world.model.get_species_type_by_id(sid2)["name"]
        return 'Pair[%s: %s, %s, (%s, %s)]' % (
            self.domain_id,
            self.pid_particle_pair1[0],
            self.pid_particle_pair2[0],
            name1, name2)

class SimplePair(Pair):

    @classmethod
    def get_min_shell_size(cls, single1, single2, geometrycontainer):
        # 1. Determine min shell size for a putative SimplePair.
        # TODO Note that the calculation is dimension independent. Also
        # 'size' is not the right term here

        assert single1.structure == single2.structure

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

#        assert (pos1 - single1.shell.shape.position).sum() == 0        # TODO Not sure why this is here
        sigma = radius1 + radius2
        r0 = geometrycontainer.world.distance(pos1, pos2)
#        distance_from_sigma = r0 - sigma        # the distance between the surfaces of the particles
        assert r0 >= sigma, \
            'distance_from_sigma (pair gap) between %s and %s = %s < 0' % \
            (single1, single2, FORMAT_DOUBLE % (r0 - sigma))

        # calculate the minimum shell size for the interparticle domain
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D12 = D1 + D2
        dist_from_com1 = r0 * D1 / D12                  # particle distance from CoM
        dist_from_com2 = r0 * D2 / D12
        iv_shell_size1 = dist_from_com1 + radius1       # the shell should surround the particles
        iv_shell_size2 = dist_from_com2 + radius2

        # also fix the minimum shellsize for the CoM domain
        com_shell_size = max(radius1, radius2)

        # calculate total radii including the margin for the burst volume for the particles
        shell_size = max(iv_shell_size1 + com_shell_size + radius1 * (1 - Domain.SINGLE_SHELL_FACTOR),
                         iv_shell_size2 + com_shell_size + radius2 * (1 - Domain.SINGLE_SHELL_FACTOR))

        return shell_size


    @classmethod
    def get_max_shell_size(cls, shell_center, single1, single2, geometrycontainer, domains):
        # returns the maximum shell size based on geometric constraints
        # It makes sure that surrounding NonInteractionSingles have at least the
        # 'bursting volume' of space.

#        pos1 = single1.pid_particle_pair[1].position
#        pos2 = single2.pid_particle_pair[1].position
#        D1 = single1.pid_particle_pair[1].D
#        D2 = single2.pid_particle_pair[1].D
#        D12 = D1 + D2

        # TODO it makes no sence to calculate the CoM in the world object
#        com  = geometrycontainer.world.calculate_pair_CoM(pos1, pos2, D1, D2)
#        shell_center = geometrycontainer.world.apply_boundary(com)

        # TODO get_closest_obj searches the spherical surrounding although this doesn't make sence for
        # a 2D or 1D pair
        closest, closest_distance = geometrycontainer.get_closest_obj(shell_center, domains, ignore=[single1.domain_id,
                                                                      single2.domain_id],
                                                                      ignores=[single1.structure.id])
        if __debug__:
            log.debug('Pair closest neighbor: %s %s' % \
                      (closest, FORMAT_DOUBLE % closest_distance))

        assert closest_distance > 0



        # 3. Check if bursted Singles can still make a minimal shell.
        # TODO This should be changed to a more general algorithm
        # that checks if the closest nearest relevant shell is a
        # NonInteractionSingle or smt else.
        # Now it only processes the shells that were actually bursted before.

#        closest, closest_shell_distance = None, numpy.inf
#        for b in burst:
#            if isinstance(b, NonInteractionSingle):
#                bpos = b.shell.shape.position
#                d = self.world.distance(shell_center, bpos) - \
#                    b.pid_particle_pair[1].radius * self.SINGLE_SHELL_FACTOR
#                if d < closest_shell_distance:
#                    closest, closest_shell_distance = b, d
#


        if isinstance(closest, NonInteractionSingle):
            closest_particle_distance = geometrycontainer.world.distance(
                    shell_center, closest.pid_particle_pair[1].position)

            closest_min_radius = closest.pid_particle_pair[1].radius
            closest_min_shell = closest_min_radius * Domain.SINGLE_SHELL_FACTOR

            # options for shell size:
            # a. upto the closest shell assuming that this shell is bigger than the minimum
            # b. The distance to a bursted single including its minimal shell
            shell_size = min(closest_distance, closest_particle_distance - closest_min_shell)
            assert shell_size <= closest_distance

        # TODO this is not correct. If the closest domain is not a NonInteractionSingle
        # we still have to make sure that the (non closest) NonInteractionSingles have
        # the required space
        else:
            assert isinstance(closest, (InteractionSingle, Pair, Multi, Surface, None.__class__))
            shell_size = closest_distance 


        return min(geometrycontainer.get_max_shell_size(),
                   shell_size)


    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
                      r0, shell_size, rrs, structure):
        Pair.__init__(self, domain_id, single1, single2, shell_id,
                      r0, rrs)

        self.structure = structure              # structure on which both particles live

        self.a_R, self.a_r = self.determine_radii(r0, shell_size)
        # set the radii of the inner domains as a function of the outer protective domain

        self.shell = self.create_new_shell(shell_center, shell_size, domain_id)


    def getSigma(self):
        return self.pid_particle_pair1[1].radius + \
               self.pid_particle_pair2[1].radius
    sigma = property(getSigma)

    def get_com(self):
        return self.shell.shape.position
    com = property(get_com)

    @classmethod
    def do_transform(cls, single1, single2, world):

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D

        com = world.calculate_pair_CoM(pos1, pos2, D1, D2)
        com = world.apply_boundary(com)

        pos2t = world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        return com, iv

    def do_back_transform(self, com, iv):
        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        pos1 = com - iv * (D1 / self.D_tot)
        pos2 = com + iv * (D2 / self.D_tot)

        return pos1, pos2

    def get_shell_size(self):
        return self.shell.shape.radius

    def determine_radii(self, r0, shell_size):
        """Determine a_r and a_R from the size of the protective domain.

        Todo. Make dimension (1D/2D/3D) specific someday. Optimization only.

        """
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        LD_MAX = 20 # temporary value
        a_r_max = LD_MAX * (r0 - self.sigma) + self.sigma
        # a_r_max is the maximum size of a_r for a given maximum ratio of l/delta

        # Make sure that D1 != 0 to avoid division by zero in the followings.
        if D1 == 0:
            D1, D2 = D2, D1     # shouldn't we also here swap the particle radii if we're swapping diffusion contants?

        shell_size /= SAFETY    #??

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shell_size +\
                math.sqrt(D2 / D1) * (radius1 - shell_size) - radius2 >= 0:
            Da = D1
            Db = D2
            radiusa = radius1
            radiusb = radius2
        else:
            Da = D2
            Db = D1
            radiusa = radius2
            radiusb = radius1


        #aR
        a_R = (D_geom * (Db * (shell_size - radiusa) + \
                         Da * (shell_size - r0 - radiusa))) /\
              (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        a_r = (D_geom * r0 + D_tot * (shell_size - radiusa)) / (Da + D_geom)


        # Now if the planned domainsize for r is too large for proper convergence of the Green's
        # functions, make it the maximum allowed size
        if (a_r > a_r_max):
            a_r = a_r_max
            a_R = shell_size - radiusa - a_r * Da / D_tot
            if __debug__:
                log.info('domainsize changed for convergence.')

        assert a_R + a_r * Da / D_tot + radiusa >= \
               a_R + a_r * Db / D_tot + radiusb

        assert abs(a_R + a_r * Da / D_tot + radiusa - shell_size) \
            < 1e-12 * shell_size                        # here the shell_size is the relevant scale


        if __debug__:
          log.debug('a %g, r %g, R %g r0 %g' % 
                 (shell_size, a_r, a_R, r0))
        if __debug__:
            tr = ((a_r - r0)**2) / (6 * self.D_r)       # the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (6 * self.D_R)          # the expected escape time of the CoM
            log.debug('tr %g, tR %g' % (tr, tR))


        assert a_r > 0
        assert a_r > r0, '%g %g' % (a_r, r0)
        assert a_R > 0 or (a_R == 0 and (D1 == 0 or D2 == 0))

        return a_R, a_r


class SphericalPair(SimplePair):
    """2 Particles inside a (spherical) shell not on any surface.

    """
    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
                 r0, shell_size, rrs, structure):
        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
                      r0, shell_size, rrs, structure)

    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Green's function for interparticle vector inside absorbing 
        # sphere.  This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                              self.sigma, self.a_r)

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    # selects between the full solution or an approximation where one of
    # the boundaries is ignored
    def choose_pair_greens_function(self, r0, t):
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(6.0 * self.D_r * t)

        # if sigma reachable
        if distance_from_sigma < threshold_distance:
        
            # if shell reachable
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction3DRadAbs
                if __debug__:
                    log.debug('GF: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction3DRadInf
                if __debug__:
                    log.debug('GF: only sigma')
                return GreensFunction3DRadInf(self.D_r, self.interparticle_ktot, r0,
                                               self.sigma)

        # sigma unreachable
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF: only a')
                return GreensFunction3DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF: free')
                return GreensFunction3D(self.D_r, r0)

    def create_com_vector(self, r):
        return random_vector(r)

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        # FIXME code doesn't handle r=0 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

#        new_inter_particle_s = numpy.array([r, theta, phi])
#        new_iv = spherical_to_cartesian(new_inter_particle_s)

        # calculate the old_iv theta (theta is the angle with the zenith (z-axis) of the system)
#        old_iv_theta = vector_angle_against_z_axis(old_iv)
        # the old_iv phi is taken to be zero. This means that the theta is in the plane of
        # the z-axis and the old_iv

            # We now rotate the new interparticle vector in the plane of
            # the old_iv and the z axis. This menas along the
            # rotation axis that is perpendicular to both the
            # z-axis and the original interparticle vector. The rotation
            # angle is the angle between the old_iv and z axis.
        
            # the rotation axis is a normalized cross product of
            # the z-axis and the original vector.
            # rotation_axis = crossproduct([0,0,1], inter_particle)

#            rotation_axis = crossproduct_against_z_axis(old_iv)
#            rotation_axis = normalize(rotation_axis)
#            new_iv = rotate_vector(new_iv, rotation_axis, old_iv_theta)

        if theta == 0.0:
            # no rotation is necessary -> new_iv = new_iv
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # alternative calculation
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi and adjust length
            phi = myrandom.uniform() * 2 * Pi
            old_iv = normalize(old_iv)
            new_iv = rotate_vector(new_iv, old_iv, phi)
        else:
            # angle == numpi.pi -> just mirror the old_iv 
            new_iv = -old_iv

        new_iv = (r/r0) * new_iv

        return new_iv

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


class PlanarSurfacePair(SimplePair):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. 
    (Hockey pucks).

    """
    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
                 r0, shell_size, rrs, structure):
        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
                      r0, shell_size, rrs, structure)

    def com_greens_function(self):
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
	    # TODO still doesn't work with 2D Green's functions
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                      self.sigma, self.a_r)

    def create_new_shell(self, position, radius, domain_id):
        # The half_length (thickness/2) of a hockey puck is not more 
        # than it has to be (namely the radius of the particle), so if 
        # the particle undergoes an unbinding reaction we still have to 
        # clear the target volume and the move may be rejected (NoSpace 
        # error).
        orientation = crossproduct(self.structure.shape.unit_x,
                                   self.structure.shape.unit_y)
        half_length = max(self.pid_particle_pair1[1].radius,
                          self.pid_particle_pair2[1].radius)
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def choose_pair_greens_function(self, r0, t):
        # selects between the full solution or an approximation where one of
        # the boundaries is ignored
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(4.0 * self.D_r * t)

        if distance_from_sigma < threshold_distance:
        
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction3DRadAbs
                if __debug__:
                    log.debug('GF2D: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction3DRadInf
                if __debug__:
                    log.debug('GF2D: only sigma')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DRadInf(self.D_r, self.interparticle_ktot, r0, self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF2D: only a')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF2D: free')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2D(self.D_r, r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        # note that r0 = length (old_iv)

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

#        unit_x = self.structure.shape.unit_x
#        unit_y = self.structure.shape.unit_y
#        angle = vector_angle(unit_x, old_iv)
        # FIXME THIS IS WRONG: angle returns the angle between the vectors regardsless of order
        #       now we do not see the different if the old_iv is mirrored in the unit_x
        # Todo. Test if nothing changes when theta == 0.
#        new_angle = angle + theta
#
#        new_iv = r * math.cos(new_angle) * unit_x + \
#                 r * math.sin(new_angle) * unit_y

        new_iv = (r/r0) * rotate_vector(old_iv, self.structure.shape.unit_z, theta)
        # note that unit_z can point two ways rotating the vector clockwise or counterclockwise
        # Since theta is symmetric this doesn't matter.

        return new_iv

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


class CylindricalSurfacePair(SimplePair):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface.  
    (Rods).

    """
    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
                 r0, shell_size, rrs, structure):
        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
                      r0, shell_size, rrs, structure)

    def get_v_tot(self):
        return self.pid_particle_pair2[1].v - \
               self.pid_particle_pair1[1].v
    v_tot = property(get_v_tot)
    
    v_r = property(get_v_tot)

    def get_v_R(self):
        return (self.pid_particle_pair1[1].v * 
                self.pid_particle_pair2[1].D +
                self.pid_particle_pair2[1].v *
                self.pid_particle_pair1[1].D) / self.D_tot
    v_R = property(get_v_R)

    def com_greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.D_R, self.v_R, 0.0, -self.a_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction1DRadAbs(self.D_r, self.v_r, self.interparticle_ktot, r0, self.sigma, self.a_r)

    def create_new_shell(self, position, half_length, domain_id):
        # The radius of a rod is not more than it has to be (namely the 
        # radius of the biggest particle), so if the particle undergoes 
        # an unbinding reaction we still have to clear the target volume 
        # and the move may be rejected (NoSpace error).
        radius = max(self.pid_particle_pair1[1].radius,
                     self.pid_particle_pair2[1].radius)
        orientation = self.structure.shape.unit_z
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    def create_com_vector(self, r):
        return r * self.structure.shape.unit_z

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        # note that r0 = length (old_iv)
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        # Note: using self.structure.shape.unit_z here might accidently 
        # interchange the particles.
        # That's why we would use self.shell.shape.unit_z right?!
        return (r/r0) * old_iv

    def get_shell_size(self):
        # Heads up.
        return self.shell.shape.half_length

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)


class MixedPair(Pair):

    @classmethod
    def get_min_shell_dimensions(cls, single1, single2, geometrycontainer):

        assert isinstance (single1.structure, PlanarSurface)
        assert isinstance (single2.structure, CuboidalRegion)

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D12 = D1 + D2

        # calculate the interparticle vector and project it onto the coordinate system of the structure
        # of particle1
        pos2_t = geometrycontainer.world.cyclic_transpose(pos2, pos1)
        iv = pos2_t - pos1
        iv_x =     numpy.dot(iv, single1.structure.shape.unit_x)
        iv_y =     numpy.dot(iv, single1.structure.shape.unit_y)
        iv_z = abs(numpy.dot(iv, single1.structure.shape.unit_z))

        z_scaling = cls.calc_z_scaling_factor (D1, D2)
        r0 = math.sqrt( iv_x**2 + iv_y**2 + (iv_z*z_scaling)**2)

        # calculate the minimal height z_right1 of the shell including burst radius
        # with the accompanying radius r1
        z_right1 = iv_z + radius2 * Domain.SINGLE_SHELL_FACTOR
        r1       = cls.r_right(single1, single2, r0, z_right1)

        len_iv_projected = math.sqrt(iv_x*iv_x + iv_y*iv_y)
        dist_r_from_com1 = len_iv_projected * D1 / D12
        dist_r_from_com2 = len_iv_projected * D2 / D12
        iv_shell_radius1 = dist_r_from_com1 + radius1
        iv_shell_radius2 = dist_r_from_com2 + radius2

        # fix the minimum shell size for the CoM domain
        com_shell_radius = max(radius1, radius2)

        # calculate the minimal dimensions of the protective domain including space for the
        # burst volumes of the particles
        r2 = max(iv_shell_radius1 + com_shell_radius + radius1 * (1 - Domain.SINGLE_SHELL_FACTOR),
                 iv_shell_radius2 + com_shell_radius + radius2 * (1 - Domain.SINGLE_SHELL_FACTOR))
        z_right2 = cls.z_right(single1, single2, r0, r2)

        # of both alternatives pick the largest one
        z_right, r = max((z_right1, r1),(z_right2, r2))

        # z_left is always the same
        z_left = radius1

        return r, z_left, z_right

    @classmethod
    def get_max_shell_dimensions(cls, single1, single2, geometrycontainer, domains):
        # Properties of the MixedPair system

        assert isinstance (single1.structure, PlanarSurface)
        assert isinstance (single2.structure, CuboidalRegion)

        com, iv = cls.do_transform(single1, single2, geometrycontainer.world)
        reference_point = com
        orientation_vector = normalize(single1.structure.shape.unit_z * numpy.dot (iv, single1.structure.shape.unit_z))
        r0 = length(iv)

        # Make sure the maximal cylinder fits in the maximal sphere. Matrix space
        # doesn't allow to check for shells outside the maximal sphere.
        max_cylinder_radius      = geometrycontainer.get_max_shell_size()/math.sqrt(2)
        max_cylinder_half_length = max_cylinder_radius

        # set initial values of shell dimensions
        z_left = single1.pid_particle_pair[1].radius
        z_right = max_cylinder_half_length * 2 - z_left
        r = max_cylinder_radius
        
        # get all the domains in the neighborhood
        search_point = reference_point + ((z_right - z_left)/2.0) * orientation_vector
        all_neighbor_ids = \
            geometrycontainer.get_neighbors_within_radius_no_sort(search_point,
                                                     geometrycontainer.get_max_shell_size(),
                                                     ignore=[single1.domain_id, single2.domain_id])
        all_neighbors = [domains[domain_id] for domain_id in all_neighbor_ids]

        world = geometrycontainer.world

        # adjust the maximal shell dimensions for every shell in the neighborhood
        for domain in all_neighbors:
            if isinstance(domain, Multi):
                for _, shell in domain.shell_list:
                    shell_position = shell.shape.position
                    shell_size = shell.shape.radius
                    r, z_left, z_right = cls.laurens_algorithm(single1, single2, r0, shell_position, shell_size,
                                                                     reference_point, orientation_vector, r, z_left,
                                                                     z_right, world)
            else:
                shell_position = domain.shell.shape.position
                shell_size = domain.get_shell_size()

                if domain.dt == 0.0 and domain.getD() > 0:
                    # This is one of the bursted singles.
                    # Or a particle that just escaped it's multi.
                    shell_size *= Domain.SINGLE_SHELL_FACTOR

                r, z_left, z_right = cls.laurens_algorithm(single1, single2, r0, shell_position, shell_size,
                                                                 reference_point, orientation_vector, r, z_left,
                                                                 z_right, world)



        return r, z_left, z_right

    @classmethod
    def laurens_algorithm(cls, single1, single2, r0, shell_position, shell_size, reference_point,
                          orientation_vector, r, z_left, z_right, world):

        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D


        # determine on what side the midpoint of the shell is relative to the reference_point
        shell_position_t = world.cyclic_transpose (shell_position, reference_point)
        direction = numpy.dot(shell_position_t - reference_point, orientation_vector)

        # if the shell is on the side of orientation_vector -> use z_right
        # also use this one if the shell is in plane with the reference_point
        if direction >= 0:
            phi = cls.calc_right_scalingangle(D1, D2)
            scalecenter_h0 = cls.z_right(single1, single2, r0, 0.0)   # r = 0.0
            r1_function =  cls.r_right
            z1_function = cls.z_right
            z2_function = cls.z_left
            z1 = z_right
            z2 = z_left
            direction = 1           # improve direction specification such that following calculations still work
                                    # Otherwise problems when direction = 0
        # else -> use z_left
        else:
            phi = cls.calc_left_scalingangle(D1, D2)
            scalecenter_h0 = cls.z_left(single1, single2, r0, 0.0)    # r = 0.0
            r1_function =  cls.r_left
            z1_function = cls.z_left
            z2_function = cls.z_right
            z1 = z_left
            z2 = z_right
            direction = -1

        # calculate the center from which linear scaling will take place
        # TODO this should be generalized to accomodate scaling in just z
        # TODO Also the scale center is still strangly defined.
        if phi == 0.0:
            scalecenter_r0 = r  # This is an approximation, since the h0 is infinitely far away of phi == 0,
        else:                   # we introduce an r0
            scalecenter_r0 = 0
        # TODO check this
        scale_center = reference_point + scalecenter_h0 * direction * orientation_vector
        scale_center = world.apply_boundary(scale_center)

        # calculate the vector from the scale center to the center of the shell
        shell_position_t = world.cyclic_transpose (shell_position, scale_center)
        shell_scale_center = shell_position_t - scale_center
        shell_scalecenter_z = numpy.dot( shell_scale_center, direction * orientation_vector )
        shell_scalecenter_r = length(shell_scale_center - (shell_scalecenter_z*direction*orientation_vector))
        # calculate the angle theta of the vector from the scale center to the shell with the vector
        # to the scale center (which is +- the orientation_vector)
        theta = vector_angle (shell_scale_center, direction * orientation_vector)

        psi = theta - phi

        ### check what situation arrises
        # The shell can hit the cylinder on the flat side (situation 1),
        #                                on the edge (situation 2),
        #                                or on the round side (situation 3).
        # I think this also works for phi == 0 and phi == Pi/2
        if psi <= -phi:
        # The midpoint of the shell lies on the axis of the cylinder
            situation = 1

        elif -phi < psi and psi < 0:
        # The (spherical) shell can touch the cylinder on its flat side or its edge
            if phi == Pi/2.0:
                r_tan_phi = 0
            else:
                r_tan_phi = shell_scalecenter_r/math.tan(phi)   # phi == 0 should not get here

            a_thres = shell_scalecenter_z - r_tan_phi
            if shell_size < a_thres:
                situation = 1
            else:
                situation = 2

        elif 0 <= psi and psi < (Pi/2.0 - phi):
        # The (spherical) shell can touch the cylinder on its edge or its radial side
            tan_phi = math.tan(phi)                             # phi == Pi/2 should not get here
            a_thres = shell_scalecenter_r - shell_scalecenter_z * tan_phi
            if shell_size > a_thres:
                situation = 2
            else:
                situation = 3

        elif (Pi/2.0 - phi) <= psi:
        # The shell is always only on the radial side of the cylinder
            situation = 3

        else:
        # Don't know how we would get here, but it shouldn't happen
            raise RuntimeError('Error: psi was not in valid range. psi = %s, phi = %s, theta = %s' %
                               (FORMAT_DOUBLE % psi, FORMAT_DOUBLE % phi, FORMAT_DOUBLE % theta))

        log.debug('situation = %s' % (situation))
        print "h0", scalecenter_h0
        print "scale_center" , scale_center
        print "theta = " ,theta
        print "phi= ", phi
        print "psi= ", psi
        print "shell_scale_center= ", shell_scale_center
        print "shell_scalecenter_r= ", shell_scalecenter_r
        print "shell_scalecenter_z= ", shell_scalecenter_z
#        print "a_thres=", a_thres
        print "shell_size=", shell_size


        ### Get the right values for z and r for the given situation
        if situation == 1:      # shell hits on the flat side
            z1_new = min(z1, scalecenter_h0 + shell_scalecenter_z - shell_size)
            r_new  = min(r,  r1_function(single1, single2, r0, z1_new))
            z2_new = min(z2, z2_function(single1, single2, r0, r_new))

        elif situation == 2:    # shell hits on the edge
            a_sq = shell_size*shell_size
            shell_scalecenter_len = length(shell_scale_center)
            ss_sq = shell_scalecenter_len*shell_scalecenter_len
            sin_phi = math.sin(phi)
            sin_psi = math.sin(psi)
            cos_psi = math.cos(psi)
            scalecenter_shell_dist = (shell_scalecenter_len * cos_psi - math.sqrt(a_sq - ss_sq*sin_psi*sin_psi) )
            # FIXME UGLY FIX BELOW
            scalecenter_shell_dist /= 1.1

            if phi <= Pi/4:
                print "scale z first."
                z1_new = min(z1, scalecenter_h0 + cos_phi * scalecenter_shell_dist)
                r_new  = min(r, r1_function(single1, single2, r0, z1_new))
                z2_new = min(z2, z2_function(single1, single2, r0, r_new))
            else:
                print "scale r first."
                r_new = min(r, scalecenter_r0 + sin_phi * scalecenter_shell_dist)
                z1_new = min(z1, z1_function(single1, single2, r0, r_new))
                z2_new = min(z2, z2_function(single1, single2, r0, r_new))

        elif situation == 3:    # shell hits on the round side
            r_new = min(r, scalecenter_r0 + shell_scalecenter_r - shell_size)
            z1_new = min(z1, z1_function(single1, single2, r0, r_new))
            z2_new = min(z2, z2_function(single1, single2, r0, r_new))
        else:
            raise RuntimeError('Bad situation for MixedPair shell making')

        ## DEBUGGING
        h_l = (z1_new + z2_new)/2.0
        print "h_l = ", h_l
        g = z1_new - h_l - scalecenter_h0
        print "g = ", g
        x = length(shell_scale_center) 
        print "distance = ", math.sqrt(g*g + x*x - 2.0*g*x*math.cos(theta))


        # switch the z values in case it's necessary. r doesn't have to be switched.
        r = r_new
        if direction >= 0.0:
            z_right = z1_new
            z_left  = z2_new
        else:
            z_right = z2_new
            z_left  = z1_new

        ## DEBUGGING
        print "r = " , r
        print "z_left = ", z_left
        print "z_right = ", z_right

        return r, z_left, z_right

    ### class methods that are for the right side of the cylindrical shell
    @classmethod
    def calc_right_scalingangle(cls, D1, D2):
        D_r = D1 + D2
        D_R = (D1*D2)/(D1 + D2)

        z_scaling_factor = cls.calc_z_scaling_factor(D1, D2)

        angle = math.atan( ( z_scaling_factor * ( math.sqrt( (2*D_R)/(3*D_r) ) + \
                                                      max( D1/D_r, D2/D_r )
                                                    )))
        return angle

    @classmethod
    def z_right(cls, single1, single2, r0, r):
        # if the radius is known and we want to determine the height z_right
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D_r = D1 + D2
        D_R = (D1*D2)/(D1 + D2)

        # calculate a_r such that the expected first-passage for the CoM and IV are equal
        sqrt_DRDr = math.sqrt((2*D_R)/(3*D_r))
        a_r1 = (r - radius1 + r0*sqrt_DRDr ) / (sqrt_DRDr + (D1/D_r) )
        a_r2 = (r - radius2 + r0*sqrt_DRDr ) / (sqrt_DRDr + (D2/D_r) )
        # take the smallest that, if entered in the function for r below, would lead to this z_right
        a_r = min (a_r1, a_r2)

        z_scaling_factor = cls.calc_z_scaling_factor(D1, D2)
        z_right = (a_r/z_scaling_factor) + (radius2)

        return z_right

    @classmethod
    def r_right(cls, single1, single2, r0, z_right):
        # if the z_right is known and we want to know the radius r
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D_r = D1 + D2
        D_R = (D1*D2)/(D1 + D2)

        # we first calculate the a_r, since it is the only radius that depends on z_right only.
        z_scaling_factor = cls.calc_z_scaling_factor(D1, D2)
        a_r = (z_right - (radius2)) * z_scaling_factor

        # We equalize the estimated first passage time for the CoM (2D) and IV (3D) for a given a_r
        # Note that for the IV we only take the distance to the outer boundary into account.
        a_R = (a_r - r0)*math.sqrt((2*D_R)/(3*D_r))

        # We calculate the maximum space needed for particles A and B based on maximum IV displacement
        iv_max = max( (D1/D_r * a_r + radius1),
                      (D2/D_r * a_r + radius2))

        r = a_R + iv_max
        return r

    ### class methods for the left side of the cylindrical shell
    @classmethod
    def calc_left_scalingangle(cls, D1, D2):
        return Pi/2

    @classmethod
    def z_left(cls, single1, single2, r0, r):
        return single1.pid_particle_pair[1].radius

    @classmethod
    def r_left(cls, single1, single2, r0, z_left):
        # FIXME this should return the maximum of the shell container
        return numpy.inf


    ### class methods that are general
    @classmethod
    def calc_z_scaling_factor(cls, D2d, D3d):
        # calculates the scaling factor to make the anisotropic diffusion problem into a isotropic one
        return math.sqrt( (D2d + D3d) / D3d)

    ################

    def __init__(self, domain_id, single1, single2, shell_id, shell_center, shell_radius,
                 shell_half_length, shell_orientation_vector, r0, rrs):

        Pair.__init__(self, domain_id, single1, single2, shell_id, r0, rrs)

        self.surface = single1.structure    # the surface on which particle1 lives
        self.structure = single1.structure    # the surface on which particle1 lives
        self.r0 = r0

        self.a_R, self.a_r = self.determine_radii(shell_half_length, shell_radius)

        self.shell = self.create_new_shell (shell_center, shell_radius, shell_half_length, shell_orientation_vector,
                                            domain_id)

    def determine_radii(self, half_length, radius):
        # determines the dimensions of the domains used for the Green's functions from the dimensions
        # of the cylindrical shell.
        # Note that the function assumes that the cylinder is dimensioned properly
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        # use class methods to check dimensions of the cylinder
        # note that we assume that only the z_right is scalable
        r_check = self.r_right (self.single1, self.single2, self.r0, half_length*2 - radius1)
        hl2_check = self.z_right (self.single1, self.single2, self.r0, radius) + \
                    self.z_left  (self.single1, self.single2, self.r0, radius)
        if r_check > radius:
            # radius should have been larger, adjust z_right instead
            if __debug__:
                log.debug('MixedPair: half_length was too high for radius:'
                          'radius = %.3g, r_check = %.3g, 2half_length = %.3g, 2hl_check = %.3g' %
                          (radius, r_check, half_length*2, hl2_check))

            half_length = hl2_check/2

        elif hl2_check > (half_length*2):
            # half_length should have been larger, adjust radius instead
            if __debug__:
                log.debug('MixedPair: radius was too high for half_length:'
                          'radius = %.3g, r_check = %.3g, 2half_length = %.3g, 2hl_check = %.3g' %
                          (radius, r_check, half_length*2, hl2_check))
            radius = r_check
        else:
            # dimensions were ok
            pass

        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D
        D12 = D1 + D2

        a_r = ((half_length * 2) - radius1 - (radius2)) * self.z_scaling_factor

        # calculate the maximal displacement of a particle given an a_r. Also include its radius
        space_for_iv = max( (D1/D12) * a_r + radius1,
                            (D2/D12) * a_r + radius2)
        a_R = radius - space_for_iv


        if __debug__:
            tr = ((a_r - self.r0)**2) / (6 * self.D_r)  # the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (4 * self.D_R)          # the expected escape time of the CoM
            log.debug('a_r= %s, tr= %s, a_R= %s, tR= %s' % \
                      (FORMAT_DOUBLE % a_r, FORMAT_DOUBLE % tr,
                       FORMAT_DOUBLE % a_R, FORMAT_DOUBLE % tR))

        assert (self.sigma < a_r) and (a_r < half_length*2)
        assert (0 < a_R) and (a_R < radius)

        return a_R, a_r

    def get_com(self):
        # Note that it is implied that the center of the shell is chosen such that its projection on the
        # membrane is the CoM of the particles.
        dist_to_surface = self.shell.shape.half_length - self.pid_particle_pair1[1].radius
        return self.shell.shape.position - dist_to_surface * self.shell.shape.unit_z
    com = property(get_com)

    def get_scaling_factor(self):
        D_2 = self.pid_particle_pair2[1].D      # particle 2 is in 3D and is the only contributor to diffusion
                                                # normal to the plane
        return math.sqrt(self.D_r/D_2)
    z_scaling_factor = property(get_scaling_factor)

    @classmethod
    def do_transform(cls, single1, single2, world):

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        D_1 = single1.pid_particle_pair[1].D
        D_2 = single2.pid_particle_pair[1].D

        surface = single1.structure
        # we assume that the particle of single1 lives on the membrane
        assert isinstance(surface, PlanarSurface)


        # the CoM is calculated in a similar way to a normal 3D pair
        com = (D_2 * pos1 + D_1 * pos2) / (D_1 + D_2)
        # and then projected onto the plane
        com = world.cyclic_transpose(com, surface.shape.position)
        com, _ = surface.projected_point (com)
        com = world.apply_boundary(com)

        # calculate the interparticle vector
        pos2t = world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        # calculate the scaling factor due to anisotropic diffusion of iv
        z_scaling_factor = math.sqrt((D_1 + D_2)/D_2)

        # move and recale the iv in the axis normal to the plane
        iv_z = surface.shape.unit_z * numpy.dot (iv, surface.shape.unit_z)
        iv_z_length = length(iv_z)

        new_iv_z_length = (iv_z_length) * z_scaling_factor
        new_iv_z = (new_iv_z_length / iv_z_length) * iv_z

        iv = iv - iv_z + new_iv_z

        return com, iv

    def do_back_transform(self, com, iv):
        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        weight1 = D1 / self.D_tot
        weight2 = D2 / self.D_tot

        # get the coordinates of the iv relative to the system of the surface (or actually the shell)
        iv_x = self.surface.shape.unit_x * numpy.dot(iv, self.surface.shape.unit_x)
        iv_y = self.surface.shape.unit_y * numpy.dot(iv, self.surface.shape.unit_y)

        # reflect the coordinates in the unit_z direction back to the side of the membrane
        # where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        # plane of the membrane
        iv_z_length = abs(numpy.dot(iv, self.shell.shape.unit_z))
        # do the reverse scaling
        iv_z_length = iv_z_length / self.z_scaling_factor

        # if the particle is overlapping with the membrane, make sure it doesn't
        if iv_z_length < self.pid_particle_pair2[1].radius:
            iv_z_length = self.pid_particle_pair2[1].radius * MINIMAL_SEPARATION_FACTOR

        iv_z = self.shell.shape.unit_z * iv_z_length

        pos1 = com - weight1 * (iv_x + iv_y)
        pos2 = com + weight2 * (iv_x + iv_y) + \
               iv_z #/self.z_scaling_factor + \
#               self.shell.shape.unit_z * self.pid_particle_pair2[1].radius

        return pos1, pos2

    def get_sigma(self):
        # rescale sigma to correct for the rescaling of the coordinate system
        # This is the sigma that is used for the evaluation of the Green's function and is in this case slightly
        # different than the sums of the radii of the particleso
        xi = self.z_scaling_factor
        xi_inv = 1.0/xi
        alpha = math.acos(xi_inv)
        sigma = self.pid_particle_pair1[1].radius + self.pid_particle_pair2[1].radius
        rho = abs(sigma * math.sqrt(0.5 + (alpha * xi/(2.0*math.sin(alpha)))))
        return rho
    sigma = property (get_sigma)

    def choose_pair_greens_function(self, r0, t):
        return self.iv_greens_function(r0)

    def com_greens_function(self):
        # diffusion of the CoM of a mixed pair is only two dimensional
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # diffusion of the interparticle vector is three dimensional
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                      self.sigma, self.a_r)

    def create_new_shell(self, position, radius, half_length, orientation_vector, domain_id):
        # It is assumed that the parameter are correct -> the 'bottom' of the cylinder
        # sticks through the membrane enough to accomodate the particle in the membrane
        # assert that orientation_vector is on same line as surface.shape.unit_z
        return CylindricalShell(domain_id, Cylinder(position, radius, orientation_vector,
                                                    half_length))

    def get_shell_size(self):
        return self.shell.shape.radius

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv):
        # Note: old_iv is actually the r0_vector
        # FIXME This is actually the same method as for SphericalPair

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        theta = draw_theta_wrapper(gf, r, dt)   # draw_theta return a value between -Pi and Pi

        if theta == 0.0:
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi
            phi = myrandom.uniform() * 2 * Pi
            old_iv = normalize(old_iv)
            new_iv = rotate_vector(new_iv, old_iv, phi)
        else:
            # theta == Pi -> just mirror the old_iv
            new_iv = -old_iv

        #adjust length of the vector to new r
        new_iv = (r/r0) * new_iv

        return new_iv

    def __str__(self):
        return 'Mixed' + Pair.__str__(self)

