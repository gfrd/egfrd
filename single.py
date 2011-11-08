from _gfrd import *
from _greens_functions import *
from greens_function_wrapper import *
from constants import EventType
import utils

from domain import (
    Domain,
    ProtectiveDomain)

__all__ = [
    'CylindricalSurfaceSingle',
    'PlanarSurfaceSingle',
    'SphericalSingle',
    'Single',
    'NonInteractionSingle',
    'InteractionSingle',
    'CylindricalSurfaceInteraction',
    'PlanarSurfaceInteraction',
    ]


class Single(ProtectiveDomain):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see 
    coordinate.py. For each coordinate the Green's function is 
    specified.

    """
    def __init__(self, domain_id, pid_particle_pair, reactionrules, structure):
	ProtectiveDomain.__init__(self, domain_id)

        self.multiplicity = 1				# all singles have only one particle

        self.pid_particle_pair = pid_particle_pair	# the particle in the single
        self.structure = structure			# the structure in/on which the particle lives

        self.reactionrules = reactionrules		# the reaction rules for mono molecular reactions
	self.rrule = None				# the reaction rule of the actual reaction (property reactionrule)
        self.k_tot = self.calc_ktot(reactionrules)
        #self.updatek_tot()				# RENAME to calc_ktot and MOVE to ProtectiveDomain

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

	Note also that it calculates the positional information PRIOR to any identity or
	structure change.

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

        event_type = EventType.SINGLE_ESCAPE	# TODO Event is not an ESCAPE when in 1D there is drift
        return dt, event_type

#    def updatek_tot(self):
#	# MOVE to ProtectiveDomain class. See also reasoning above for
#	# draw_reaction_time_tuple
#        self.k_tot = 0
#
#        if not self.reactionrules:
#            return
#
#        for rr in self.reactionrules:
#            self.k_tot += rr.k
#
#    def draw_reaction_rule(self):
#	# MOVE to ProtectiveDomain class
#        k_array = [rr.k for rr in self.reactionrules]
#        k_array = numpy.add.accumulate(k_array)
#        k_max = k_array[-1]
#
#        rnd = myrandom.uniform()
#        i = numpy.searchsorted(k_array, rnd * k_max)
#
#        return self.reactionrules[i]

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

    NonInteractionSingle hold two independent processes:
	* Diffusion of the particle in n (1-3) dimensions
	* Decay process of the particle

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactionrules, 
                 structure):
        Single.__init__(self, domain_id, pid_particle_pair, reactionrules,
                        structure)

        # Create shell.
	self.shell_id = shell_id
        self.shell = self.create_new_shell(pid_particle_pair[1].position,
                                      pid_particle_pair[1].radius, domain_id)

    def get_shell_size(self):
	# Note: this method only means something for Protective Domains that can
	# only be sized in one direction
	# This is predominantly used when making Interaction domains
        return self.shell.shape.radius


    def initialize(self, t):
	# self.shell.shape.radius = self.pid_particle_pair[1].radius
        self.dt = 0.0
        self.last_time = t
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

	    # TODO Either include check for non-zero drift here, or overload
	    # this method in CylindricalSurfaceSingle to check exit interface
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
	    newpos = self.shell.shape.position + displacement

        return newpos

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

    def get_min_shell_size(self):
	return self.pid_particle_pair[1].radius

    def get_max_shell_size(self, geometrycontainer, domains):
	# This calculates the maximum shell size that the single can have in the current
	# situation.
	# TODO This is now general for all NonInteractionSingles and should be different
	# for the different dimensions since they have different domains (spheres vs cylinders)
	# and scale the domains in different directions (r or z).

	# 1. Calculate the maximal dimensions of the shell
	# TODO We now just check in a sphere around the particle and size accordingly. For
	# NonInteractionSingles on the dna and membrane this is quite inefficient-> want to
	# size in the appropriate coordinate (r or z).
	singlepos = self.pid_particle_pair[1].position
        closest, distance_to_shell = \
            geometrycontainer.get_closest_obj(singlepos, domains, ignore=[self.domain_id],
                                              ignores=[self.structure.id])

	# 2. In case the nearest neighbor was a NonInteractionSingle we need to tweak the size
	# TODO shell size has no real meaning for cylinderical domains
        if isinstance(closest, NonInteractionSingle):
            new_shell_size = self.calculate_shell_size_to_single(closest, distance_to_shell,
						                 geometrycontainer)
	# TODO This is incorrect. The closest can be a Pair but still NonInteractionSingles may need
	# more space
        else:  # Pair or Multi or Surface
            new_shell_size = distance_to_shell / SAFETY

	return min(new_shell_size, geometrycontainer.get_max_shell_size())


class SphericalSingle(NonInteractionSingle):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector:
          theta, phi.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactionrules, 
                 structure):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactionrules, structure)

    def greens_function(self):
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_inner_a())

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(self.domain_id, Sphere(position, radius))

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
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactionrules, 
                 structure):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactionrules, structure)

    def greens_function(self):
        return GreensFunction2DAbsSym(self.D,
                                          self.get_inner_a())

    def create_new_shell(self, position, radius, domain_id):
        # The half_length (thickness) of a hockey puck is not more than 
        # it has to be (namely the radius of the particle), so if the 
        # particle undergoes an unbinding reaction we still have to 
        # clear the target volume and the move may be rejected (NoSpace 
        # error).
        orientation = normalize(
            utils.crossproduct(self.structure.shape.unit_x,
                               self.structure.shape.unit_y))
        half_length = self.pid_particle_pair[1].radius
        return CylindricalShell(self.domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

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
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactionrules, 
                 structure):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactionrules, structure)

    def getv(self):
        return 0 # TODO, include drift.
        return self.pid_particle_pair[1].v
    v = property(getv)

    def get_shell_size(self):
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell.shape.half_length

    def get_inner_a(self):
	# Here the free coordinate is not the radius of the domain as in the
	# PlanarSurfaceSingle and SphericalSingle, although it is treated as a
	# radius in NonInteractionSingle.draw_new_position for simplicity
	return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
	inner_half_length = self.get_inner_a()
        return GreensFunction1DAbsAbs(self.D, self.v, 0.0, -inner_half_length, inner_half_length)

    def create_new_shell(self, position, half_length, domain_id):
        # The radius of a rod is not more than it has to be (namely the 
        # radius of the particle), so if the particle undergoes an 
        # unbinding reaction we still have to clear the target volume 
        # and the move may be rejected (NoSpace error).
        radius = self.pid_particle_pair[1].radius
        orientation = self.structure.shape.unit_z
        return CylindricalShell(self.domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, z):
        if utils.feq(z, self.get_inner_a()):
            # Escape, can be either to the left or to the right.
	    # The side of escape should be decided on by the flux through
	    # both boundaries at the escape time
	    # TODO, include drift
            z = myrandom.choice(-1, 1) * z 
        return z * self.shell.shape.unit_z

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


class InteractionSingle(Single):
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
    def __init__(self, domain_id, pid_particle_pair, reactionrules, structure, shell_id,   
                 interactionrules, surface):

        Single.__init__(self, domain_id, pid_particle_pair, reactionrules,
                        structure)

	self.shell_id = shell_id
        self.interactionrules = interactionrules	# the surface interaction 'reactions'
	self.intrule = None
	self.interaction_ktot = self.calc_ktot(interactionrules)
	self.surface = surface				# The surface with which the particle is
							# trying to interact

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

    def initialize(self, t):	
	# initialize the domain object with the appropriate time to allow reuse of
	# the Domain at different times
	self.dt = 0
	self.last_time = t
	self.event_type = None

	# This is defined here, but only used in the subclasses
    def create_new_shell(self, position, radius, orientation_vector, half_length):
        return CylindricalShell(self.domain_id, Cylinder(position, radius,
                                                    orientation_vector, half_length))

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple(),
                   self.draw_iv_event_time_tuple())	# above two events also occur in NonInteractingSingle
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
        if event_kind == PairEventKind.IV_REACTION:	# Note REACTION -> INTERACTION
            return EventType.IV_INTERACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()


class PlanarSurfaceInteraction(InteractionSingle):
    """1 Particle close to a PlanarSurface, inside a cylindrical shell 
    placed against the surface with its flat face.

        * Particle coordinates inside shell: r, theta, z.
        * Coordinates: radial r, cartesian z.
        * Initial position: r = 0, z = z0 (relative to center shell).
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, pid_particle_pair, reactionrules, structure,
                 shell_id, shell_center, shell_radius, shell_half_length, shell_orientation_vector, z0,
		 interactionrules, surface):
	# shell_orientation_vector and z0 are supplied by arguments because Domains can't do
	# spatial calculations (they don't know about the world)
	# NOTE shell_orientation_vector should be normalized and pointing from the surface to the particle!

	InteractionSingle.__init__(self, domain_id, pid_particle_pair, reactionrules,
				   structure, shell_id, interactionrules, surface)

	# z0 is the initial position of the particle relative to the center
	# of the shell in the direction normal to the surface.
	self.z0 = z0
	self.shell = self.create_new_shell(shell_center, shell_radius, 
					   shell_orientation_vector, shell_half_length)

    def get_inner_dz_right(self):
	# This is the distance the particle can travel from its initial position
	# to the flat boundary away from the surface
	return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_left(self):
	# This is the distance the particle can traval towards the surface
	# Note that it is implied that the shell sticks out at the other end
	# of the planar surface by the particle radius, and that the surface
	# has no thickness.
	dist_to_surface = self.shell.shape.half_length - self.pid_particle_pair[1].radius
	return dist_to_surface - self.pid_particle_pair[1].radius

    def greens_function(self):
        # This is the Greens function that does not produce REACTION event
        #return GreensFunction3DAbsSYm(self.getD(), get_inner_a())
        return GreensFunction2DAbsSym(self.D, self.get_inner_a())

    def iv_greens_function(self):
	# The green's function that also modelles the association of the particle
	# with the planar surface.
	# TODO choose interactionrule
        return GreensFunction1DRadAbs(self.D, self.interaction_ktot, self.z0,
				      -self.get_inner_dz_left(), self.get_inner_dz_right())
#        return GreensFunction1DRadAbs(self.D, 0, self.z0,
#				      -self.get_inner_dz_left(), self.get_inner_dz_right())

    def draw_new_position(self, dt, event_type):
	"""Draws a new position for the particle for a given event type at a given time
	Note that it also 'moves the particle into the surface' when an interaction
	takes place.
	"""
	oldpos = self.pid_particle_pair[1].position
	if self.D == 0:
	    newpos = oldpos
	elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
	    newpos = oldpos		# just return oldpos, particle will be removed anyway
	else:
	    # calculate r vector
            a = self.get_inner_a()
	    if event_type == EventType.SINGLE_ESCAPE:
		r = a
	    else:
		gf = self.greens_function()
        	r = draw_r_wrapper(gf, dt, a)

	    x, y = random_vector2D(r)
	    vector_r = x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

	    # calculate z vector
            # Todo. Cartesian coordinate will return absolute position.
            # Heads up. size_of_domain is different from the half_length of 
            # the cylinder, because the domain ends where the surface 
            # starts, while the cylinder is continuing into the surface.
            z_surface = -self.get_inner_dz_left()	# left side of the inner domain
	    z_not_surface = self.get_inner_dz_right()	# right side of the inner domain (away from the surface)

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
        	z = draw_r_wrapper(gf_iv, dt, z_not_surface, z_surface)	# last two args= a, sigma

	    vector_z = z * self.shell.shape.unit_z

	    # The new position is relative to the center of the shell
	    newpos = self.shell.shape.position + vector_r + vector_z

        return newpos


    def get_shell_size(self):
	# REMOVE this method, it doesn't mean anything here.
	# THis method is only used for making an Interaction
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell.shape.radius

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
    def __init__(self, domain_id, pid_particle_pair, reactionrules, structure,
		 shell_id, shell_center, shell_radius, shell_half_length, shell_orientation_vector, z0,
		 unit_r, r0, interactionrules, surface):

        InteractionSingle.__init__(self, domain_id, pid_particle_pair, reactionrules,
				   structure, shell_id, interactionrules, surface)

	# z0 is implied to be zero (the particle being in the center of the shell in the z direction)
	self.z0 = z0
	self.unit_r = unit_r	# This is the vector from the cylinder to the particle
	self.r0 = r0		# r0 is the distance from the center of the cylinder to the
				# center of the particle
        self.shell = self.create_new_shell(shell_center, shell_radius,
                                           shell_orientation_vector, shell_half_length)


    def get_inner_dz_left(self):
	return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_right(self):
	return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_sigma(self):
	return self.surface.shape.radius + self.pid_particle_pair[1].radius

    def greens_function(self):
	# The greens function not used for the interaction but for the other coordinate
        # Free diffusion in z direction, drift is zero.
        return GreensFunction1DAbsAbs(self.D, 0.0, self.z0,
                                      -self.get_inner_dz_left(),
                                      self.get_inner_dz_right())

    def iv_greens_function(self):
	# Green's function used for the interaction
        # Interaction possible in r direction.
        # TODO 2D.
        return GreensFunction3DRadAbs(self.D, self.interaction_ktot, self.r0, self.get_inner_sigma(),
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

            # Direction matters, so use shell.shape.unit_z instead of 
            # structure.shape.unit_z.
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

	return newpos

    def get_shell_size(self):
	# REMOVE this method, it doesn't mean anything here.
	# THis method is only used for making an Interaction
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell.shape.half_length

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

