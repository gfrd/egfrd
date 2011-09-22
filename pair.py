from _gfrd import *
from constants import EventType
from _greens_functions import *
from greens_function_wrapper import *
from domain import *

__all__ = [
    'CylindricalSurfacePair',
    'PlanarSurfacePair',
    'SphericalPair',
    'Pair',
    'SimplePair',
    'MixedPair',
    ]

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
                 rt):
	ProtectiveDomain.__init__(self, domain_id)

        self.multiplicity = 2

        self.single1 = single1
        self.single2 = single2 
	self.pid_particle_pair1 = single1.pid_particle_pair
	self.pid_particle_pair2 = single2.pid_particle_pair

        self.rt = rt

	self.shell_id = shell_id

    def __del__(self):
        if __debug__:
            log.debug('del %s' % str(self))

    def get_D_tot(self):
        return self.pid_particle_pair1[1].D + \
               self.pid_particle_pair2[1].D
    D_tot = property(get_D_tot)

    def get_D_R(self):
        return (self.pid_particle_pair1[1].D *
                self.pid_particle_pair2[1].D) / self.D_tot
    D_R = property(get_D_R)

    def initialize(self, t):
	# re-initialize the time parameters of the Domain to reuse it.
        self.last_time = t
        self.dt = 0
        self.event_type = None

    # draws a first event time and the (not fully specified) type of event
    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(r0))	# uses the full solution to draw IV first event time
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
        new_com = self.draw_new_com(dt, event_type)
        new_iv = self.draw_new_iv(dt, r0, old_iv, event_type)

        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        newpos1 = new_com - new_iv * (D1 / self.D_tot)
        newpos2 = new_com + new_iv * (D2 / self.D_tot)
        return newpos1, newpos2

    def draw_new_com(self, dt, event_type):
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        displacement = self.create_com_vector(r)

        # Add displacement to old CoM. This assumes (correctly) that 
        # r0=0 for the CoM. Compare this to 1D singles, where r0 is not  
        # necesseraly 0.
        return self.com + displacement

    def draw_new_iv(self, dt, r0, old_iv, event_type):
        gf = self.choose_pair_greens_function(r0, dt)
        if event_type == EventType.IV_ESCAPE:
            r = self.a_r
        elif event_type == EventType.IV_REACTION:
            r = self.sigma
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
    def __init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, structure):
	Pair.__init__(self, domain_id, single1, single2, shell_id,
                      r0, rt)

        self.structure = structure		# structure on which both particles live

        self.a_R, self.a_r = self.determine_radii(r0, shell_size)
	# set the radii of the inner domains as a function of the outer protective domain

        self.shell = self.create_new_shell(com, shell_size, domain_id)


    def getSigma(self):
        return self.pid_particle_pair1[1].radius + \
               self.pid_particle_pair2[1].radius
    sigma = property(getSigma)

    def get_com(self):
        return self.shell.shape.position
    com = property(get_com)

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
            D1, D2 = D2, D1	# shouldn't we also here swap the particle radii if we're swapping diffusion contants?

        shell_size /= SAFETY	#??

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

        assert a_R + a_r * Da / D_tot + radius1 >= \
               a_R + a_r * Db / D_tot + radius2

        assert abs(a_R + a_r * Da / D_tot + radiusa - shell_size) \
            < 1e-12 * shell_size			# here the shell_size is the relevant scale


        if __debug__:
          log.debug('a %g, r %g, R %g r0 %g' % 
                 (shell_size, a_r, a_R, r0))
        if __debug__:
            tr = ((a_r - r0)**2) / (6 * self.D_tot)	# the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (6 * self.D_R)		# the expected escape time of the CoM
            log.debug('tr %g, tR %g' % (tr, tR))


        assert a_r > 0
        assert a_r > r0, '%g %g' % (a_r, r0)
        assert a_R > 0 or (a_R == 0 and (D1 == 0 or D2 == 0))

        return a_R, a_r


class SphericalPair(SimplePair):
    """2 Particles inside a (spherical) shell not on any surface.

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, structure):
        SimplePair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, structure)

    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Green's function for interparticle vector inside absorbing 
        # sphere.  This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_tot, self.rt.ktot, r0,
                                              self.sigma, self.a_r)

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(self.domain_id, Sphere(position, radius))

    # selects between the full solution or an approximation where one of
    # the boundaries is ignored
    def choose_pair_greens_function(self, r0, t):
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(6.0 * self.D_tot * t)

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
                return GreensFunction3DRadInf(self.D_tot, self.rt.ktot, r0,
                                               self.sigma)

        # sigma unreachable
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF: only a')
                return GreensFunction3DAbs(self.D_tot,
                                                                 r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF: free')
                return GreensFunction3D(self.D_tot, r0)

    def create_com_vector(self, r):
        return random_vector(r)

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        new_inter_particle_s = numpy.array([r, theta, 
                                         myrandom.uniform() * 2 * Pi])
        new_iv = spherical_to_cartesian(new_inter_particle_s)

        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotation_axis = crossproduct([0,0,1], inter_particle)
        angle = vector_angle_against_z_axis(old_iv)
        if angle % numpy.pi != 0.0:
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            rotated = rotate_vector(new_iv, rotation_axis, angle)
        elif angle == 0.0:
            rotated = new_iv
        else:
            rotated = numpy.array([new_iv[0], new_iv[1], - new_iv[2]])
        return rotated

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


class PlanarSurfacePair(SimplePair):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. 
    (Hockey pucks).

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, structure):
        SimplePair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, structure)

    def com_greens_function(self):
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
	# TODO still doesn't work with 2D Green's functions
        return GreensFunction3DRadAbs(self.D_tot, self.rt.ktot, r0,
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
            math.sqrt(4.0 * self.D_tot * t)

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
                # return GreensFunction2DRadInf(self.D_tot, self.rt.ktot, r0, self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF2D: only a')
                return self.iv_greens_function(r0)
	        # Todo
                # return GreensFunction2DAbs(self.D_tot, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF2D: free')
                return self.iv_greens_function(r0)
        	# Todo
                # return GreensFunction2D(self.D_tot, r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        #FIXME: need better handling of angles near zero and pi?
        unit_x = self.structure.shape.unit_x
        unit_y = self.structure.shape.unit_y
        angle = vector_angle(unit_x, old_iv)
        # Todo. Test if nothing changes when theta == 0.
        new_angle = angle + theta

        new_iv = r * math.cos(new_angle) * unit_x + \
                 r * math.sin(new_angle) * unit_y

        return new_iv

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


class CylindricalSurfacePair(SimplePair):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface.  
    (Rods).

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, structure):
        SimplePair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, structure)

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
	# TODO Fix ugly hack to avoid k=0 below
        #return GreensFunction1DRadAbs(self.D_tot, self.v_r, self.rt.ktot, r0, self.sigma, self.a_r)
        return GreensFunction1DRadAbs(self.D_tot, self.v_r, self.rt.ktot + 1e-10, r0, self.sigma, self.a_r)

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
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        # Note: using self.structure.shape.unit_z here might accidently 
        # interchange the particles.
        return r * normalize(old_iv)

    def get_shell_size(self):
        # Heads up.
        return self.shell.shape.half_length

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)


class MixedPair(Pair):

    def __init__(self, domain_id, single1, single2, shell_id, shell_center, shell_radius,
		 shell_half_length, shell_orientation_vector,
		 unit_r, r0, pair_reaction_rules):
	Pair.__init__(self, domain_id, single1, single2, shell_id, r0, 
                 pair_reaction_rules)


	self.surface = surface	# the surface on which particle1 lives

    def determine_radii(self):
	#TODO INCORRECT
	a_r = self.shell.shape.half_length*2.0 - self.pid_particle_pair1[1].radius
	a_R = self.shell.shape.radius - a_r*0.5
	return 

    def get_com(self):
	# TODO INCORRECT
	return self.shell.shape.position
    com = property(get_com)

    def get_sigma(self):
	# rescale sigma to correct for the rescaling of the coordinate system
	# TODO check if this is actually correct
	D_1 = self.pid_particle_pair1[1].D
	D_2 = self.pid_particle_pair2[1].D
	xi_inv = math.sqrt(D_2/(D_1 + D_2))
	alpha = math.acos(xi_inv)
	sigma = self.pid_particle_pair1[1].radius + self.pid_particle_pair2[1].radius
	rho = abs(sigma*math.sqrt(0.5 + (alpha/(2.0*xi_inv*math.sin(alpha)))))
	return rho
    sigma = property (get_sigma)

    def com_greensfunction(self):
	# diffusion of the CoM of a mixed pair is only two dimensional
	return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greensfunction(self):
	# diffusion of the interparticle vector is three dimensional
	return GreensFunction3DRadAbs(self.D_tot, self.rt.ktot, r0,
				      self.sigma, self.a_r)

    def create_new_shell(self, position, radius, orientation_vector, half_length):
	# It is assumed that the parameter are correct -> the 'bottom' of the cylinder
	# sticks through the membrane enough to accomodate the particle in the membrane
	# assert that orientation_vector is on same line as surface.shape.unit_z
	return CylindricalShell(self.domain_id, Cylinder(position, radius, orientation_vector,
							 half_length))

    def create_com_vector(self, r):
	x, y = random_vector2D(r)
	return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def create_iv_vector(self, gf, r, dt, r0, old_iv):
	# TODO This is not yet correct -> stub
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        new_inter_particle_s = numpy.array([r, theta,
                                         myrandom.uniform() * 2 * Pi])
        new_iv = spherical_to_cartesian(new_inter_particle_s)

        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.

        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotation_axis = crossproduct([0,0,1], inter_particle)
        angle = vector_angle_against_z_axis(old_iv)
        if angle % numpy.pi != 0.0:
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            rotated = rotate_vector(new_iv, rotation_axis, angle)
        elif angle == 0.0:
            rotated = new_iv
        else:
            rotated = numpy.array([new_iv[0], new_iv[1], - new_iv[2]])
        
	# TODO
	# rescale z component by 1/xi (z component becomes slightly smaller, since D = D_A + D_B
	# was used for calculation, whereas in reality D is only D_B
	xi_inv = math.sqrt((D_2)/D_1 + D_2)
	z /= xi
	# position is reflected back in the membrane if necessary
	return rotated 
