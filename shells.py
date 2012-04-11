#!/usr/bin/python

import math

from _gfrd import (
    Sphere,
    SphericalShell,
    Cylinder,
    CylindricalShell,
    CylindricalSurface,
    Plane,
    PlanarSurface,
    CuboidalRegion,
    )

from utils import *

__all__ = [
    'ShellmakingError',
    'testShellError',
    'Others',
    'NonInteractionSingles',
    'testSingle',
    'testNonInteractionSingle',
    'testInteractionSingle',
    'testTransitionSingle',
    'testPair',
    'testSimplePair',
    'testMixedPair2D3D',
    'hasSphericalShell',
    'hasCylindricalShell',   
    'SphericalSingletestShell',
    'SphericalPairtestShell',
    'PlanarSurfaceSingletestShell',
    'PlanarSurfacePairtestShell',
    'PlanarSurfaceTransitionSingletestShell',
    'CylindricalSurfaceSingletestShell',
    'CylindricalSurfacePairtestShell',
    'PlanarSurfaceInteractiontestShell',
    'CylindricalSurfaceInteractiontestShell',
    'CylindricalSurfaceSinktestShell',
    'MixedPair2D3DtestShell',
    ]

# These are exceptions that are raise when trying to make a new domain
class ShellmakingError(Exception):
    # This ShellmakingError is raised when a Spherical or Cylindrical shell could
    # not be made, usually due to the fact that the possible shell dimensions are
    # smaller than the minimal dimensions.
    pass

class protoDomainError(Exception):
    # The protoDomainError is raised when something in the initialization of the variables
    # that is NOT related to the shell goes wrong. This will lead later to a testShellError
    pass

class testShellError(Exception):
    # The testShellError is raised when the testShell could not be made, usually
    # due to a ShellmakingError
    pass

####################################

# To make sure that Pairs and Interaction domain keep some distance from NonInteractionSingle domains,
# the domain that is being made asks the neighboring domains what their size is. This allows the
# neighboring domain to change its appeared size depending on who asks.
# To this end, the two classes below define methods that allow an existing domain to return its shell
# dimensions dependent on the one who asks. The asker and the askee are now both a NonInteractionSingle
# or an Other.
# The main difference between the NonInteractionSingles and Others is that the NonInteractionSingle
# keep the Other domains at a distance such that the probability that the Other domain get bursted is
# minimized.
class Others(object):
# Subclasses of Others are: all Pair, Interaction and Multi domains.

    def __init__(self):
        pass    # do nothing

    # The next two methods are only relevant for HASSHELL subclasses that derive from this class.
    def shell_list_for_single(self):
        # The shell dimensions for the Other classes are unmodified, regardless of the asker.
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list


    # The next method is only relevant for TESTSHELL subclasses that derive from this class.
    def get_neighbor_shell_list(self, neighbor):
        # This gets the shells of an existing domain given that the testShell is of the type Other
        return neighbor.shell_list_for_other()

class NonInteractionSingles(object):
# Subclasses of NonInteractionSingles are: all NonInteractionSingle domains (Spherical, etc)

    def __init__(self):
        pass    # do nothing

    # The next two methods are only relevant for HASSHELL subclasses that derive from this class.
    # NOTE below two methods have quite a different purpose!
    def shell_list_for_single(self):
        # In case the asker is also a NonInteractionSingle, the shell dimensions of the current domain
        # is modified such that the domain is at least the size of the Multi shell.

        # NOTE: the Multi shell is spherical whereas the real shell of the NonInteractionSingle can
        # be cylinderical. This method should therefore always return at the Multi shell.
        # In practice it will return the multi shell when the real shell fits inside the multi shell (typically when
        # the domain is just initialized), and will return the real shell when the domain exceeds the multi shell.
        # This means that there may not be space for the multi shell when the domain has non-null dimensions.
        # This should be ok, because no multi will be directly made when the domain has non-null dimensions but
        # will only be made after first bursting nearby domains.

        pass    # needs to be specified in the subclass


    def shell_list_for_other(self):
        # In case the asker is an Other, the shell dimensions of the current NonInteractionSingle domain
        # is modified such that the domain is at least its minimum size (SINGLE_SHELL_FACTOR).
        # This is to reduce unnecessary bursting of the (expensive) Other domains.

        # NOTE: Since the burst radius is spherical, a spherical shell of the size of the burst radius is returned
        # when the real shell is small enough to fit inside the burst radius. When the real shell no longer fits inside
        # the burst radius the real shell is returned.
        # This means that more bursting will take place than necessary because cylindrical domain do not keep the domains
        # that are out of the 'diffusion dimension' at the proper distance -> TODO

        pass    # needs to be specified in the subclass


    # The next method is only relevant for TESTSHELL subclasses that derive from this class.
    def get_neighbor_shell_list(self, neighbor):
        # This gets the shells of an existing domain given that the testShell is of the type Other
        return neighbor.shell_list_for_single()

##################################
# Below classes define the minimal information needed to try to make a domain of a certain type.
# They implement the non geometric aspect of the testShells, and are inherited by the respective testShell
# classes.
class testSingle(object):

    def __init__(self, pid_particle_pair, structure):
        self.pid_particle_pair = pid_particle_pair
        self.structure = structure

class testNonInteractionSingle(testSingle, NonInteractionSingles):

    def __init__(self, pid_particle_pair, structure):
        testSingle.__init__(self, pid_particle_pair, structure)
        # Note: for the NonInteractionSingles superclass nothing is to be initialized.

class testInteractionSingle(testSingle, Others):

    def __init__(self, single, structure, surface):

        # assert that the single is reset
        assert single.is_reset()

        testSingle.__init__(self, single.pid_particle_pair, structure)
        # Note: for the Others superclass nothing is to be initialized.

        self.single  = single
        self.surface = surface

        ### initialize some constants that are common for the testInteractionSingle domains.

        # calculate the Projected_point and distance to the surface
        # Cyclic transpose needed when calling surface.projected_point!
        pos_transposed = self.world.cyclic_transpose(self.pid_particle_pair[1].position,
                                                     self.surface.shape.position)
        self.reference_point, projection_length = self.surface.projected_point(pos_transposed)

        # projection_distance can be negative in case of plane
        # note that the distance is to the center of the cylinder in case of cylinders
        self.particle_surface_distance = self.world.distance(self.surface.shape, self.pid_particle_pair[1].position)

        # The reference_vector is the normalized vector from the reference_point to the particle
        self.reference_vector = normalize(pos_transposed - self.reference_point)

class testTransitionSingle(testSingle, Others):

    def __init__(self, single, target_structure):

        testSingle.__init__(self, single.pid_particle_pair, single.structure)

        assert (self.structure.sid == target_structure.sid)     # Both structures are of same structure_type
        self.single = single

        self.origin_structure = self.structure
        self.target_structure = target_structure

        # Note that we assume the structures are exactly connecting
        self.distance_to_target_structure = self.world.distance(self.target_structure.shape, self.pid_particle_pair[1].position)

##############
class testPair(Others):

    def __init__(self, single1, single2):
        # Note: for the Others superclass nothing is to be initialized.

        # assert both reset
        assert single1.is_reset()
        assert single2.is_reset()

        self.single1 = single1
        self.single2 = single2
        self.pid_particle_pair1 = single1.pid_particle_pair
        self.pid_particle_pair2 = single2.pid_particle_pair
        self.com, self.iv = self.do_transform()         # NOTE The IV always points from particle1 to particle2
        self.r0 = length(self.iv)
        if self.r0 < self.sigma:
            raise protoDomainError('distance_from_sigma (pair gap) between %s and %s = %s < 0' % \
                (self.single1, self.single2, FORMAT_DOUBLE % (self.r0 - self.sigma)))


    def get_D_tot(self):
        return self.pid_particle_pair1[1].D + \
               self.pid_particle_pair2[1].D
    D_tot = property(get_D_tot)
    D_r   = property(get_D_tot)

    def get_D_R(self):
        return (self.pid_particle_pair1[1].D *
                self.pid_particle_pair2[1].D) / self.D_tot
    D_R = property(get_D_R)

    def do_transform(self):
        # transform the pos1 and pos2 of particles1 and 2 to the CoM and IV vectors
        pass        # overridden in subclasses

class testSimplePair(testPair):

    def __init__(self, single1, single2):

        # Simple pairs are pairs of particles that are on the same structure
        assert single1.structure == single2.structure
        self.structure = single1.structure

        testPair.__init__(self, single1, single2)

    def get_sigma(self):
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        return radius1 + radius2
    sigma = property(get_sigma)

    def do_transform(self):
        # transform the pos1 and pos2 of particles1 and 2 to the CoM and IV vectors

        pos1 = self.pid_particle_pair1[1].position
        pos2 = self.pid_particle_pair2[1].position
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        com = self.world.calculate_pair_CoM(pos1, pos2, D_1, D_2)
        com = self.world.apply_boundary(com)
        # make sure that the com is in the structure of the particle (assume that this is done correctly in
        # calculate_pair_CoM)

        pos2t = self.world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        return com, iv

    def get_min_pair_size(self):
        # This calculates the minimal 'size' of a simple pair.
        # Note that the 'size' is interpreted differently dependent on the type of pair. When the pair is a
        # SphericalPair          -> size = radius
        # PlanarSurfacePair      -> size = radius
        # CylindricalSurfacePair -> size = half_length

        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        dist_from_com1 = self.r0 * D_1 / self.D_tot      # particle distance from CoM
        dist_from_com2 = self.r0 * D_2 / self.D_tot
        iv_shell_size1 = dist_from_com1 + radius1       # the shell should surround the particles
        iv_shell_size2 = dist_from_com2 + radius2

        # also fix the minimum shellsize for the CoM domain
        com_shell_size = max(radius1, radius2)

        # calculate total radii including the margin for the burst volume for the particles
        shell_size = max(iv_shell_size1 + com_shell_size + radius1 * SINGLE_SHELL_FACTOR,
                         iv_shell_size2 + com_shell_size + radius2 * SINGLE_SHELL_FACTOR)

        return shell_size

class testMixedPair2D3D(testPair):

    def __init__(self, single2D, single3D):

        if __debug__: 
            assert isinstance(single2D.structure, PlanarSurface) 
            assert isinstance(single3D.structure, CuboidalRegion) 
        self.surface   = single2D.structure         # note that these need to be initialized before calling testPair.__init__
        self.structure = single3D.structure         # since then these are needed for calculating the transform

        testPair.__init__(self, single2D, single3D) # note: this makes self.single1/self.pid_particle_pair1 the 2D particle
                                                    #              and self.single2/self.pid_particle_pair2 the 3D particle

    def get_sigma(self):
        # rescale sigma to correct for the rescaling of the coordinate system
        # This is the sigma that is used for the evaluation of the Green's function and is in this case slightly
        # different than the sums of the radii of the particleso
        xi = self.get_scaling_factor()
        xi_inv = 1.0/xi
        alpha = math.acos(xi_inv)
        sigma = self.pid_particle_pair1[1].radius + self.pid_particle_pair2[1].radius
        rho = abs(sigma * math.sqrt(0.5 + (alpha * xi/(2.0*math.sin(alpha)))))
        return rho
    sigma = property(get_sigma)

    def do_transform(self):
        # transform the pos1 and pos2 of particles1 and 2 to the CoM and IV vectors
        pos1 = self.pid_particle_pair1[1].position
        pos2 = self.pid_particle_pair2[1].position
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        # the CoM is calculated in a similar way to a normal 3D pair
        com = self.world.calculate_pair_CoM(pos1, pos2, D_1, D_2)

        # and then projected onto the plane to make sure the CoM is in the surface
        com = self.world.cyclic_transpose(com, self.surface.shape.position)
        com, _ = self.surface.projected_point (com)
        com = self.world.apply_boundary(com)


        # calculate the interparticle vector
        pos2t = self.world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        # rescale the iv in the axis normal to the plane
        iv_z_component = numpy.dot (iv, self.surface.shape.unit_z)
        iv_z     = self.surface.shape.unit_z * iv_z_component
        new_iv_z = self.surface.shape.unit_z * (iv_z_component * self.get_scaling_factor())

        iv = iv - iv_z + new_iv_z

        return com, iv

    def get_scaling_factor(self):
        # This returns the scaling factor needed to make the anisotropic diffusion problem of the IV
        # into a isotropic one.
        D_2 = self.pid_particle_pair2[1].D      # particle 2 is in 3D and is the only contributor to diffusion
                                                # normal to the plane
        return math.sqrt(self.D_r/D_2)


##################################
# These classes define the geometric aspect of the existing domains. When making a domain, the domain
# has a Shell (hence the hasShell class name).
# They are initialized using the testShell of the appropriate type.
class hasShell(object):

    def __init__(self, testShell):
        pass

#######################
class hasSphericalShell(hasShell):
    # This class is inherited by domains that have a Spherical shell such as the SphericalSingle and a SphericalPair.
    # The class makes the appropriate shell for the domain and allows a testShell to ask what the maximum dimensions
    # of the testShell can be. This means, how far can the testShell be scaled according to its scaling rules before
    # it hits THIS Spherical shell.

    def __init__(self, sphericaltestShell, domain_id):
        # Note1: the Multi should not call __init__ because we don't want to initialize the shell
        # Note2: the hasShell.__init__ is not called because nothing happens there

        assert isinstance(sphericaltestShell, SphericaltestShell)
        self.testShell = sphericaltestShell

        # Make the appropriate shell.
        self.shell_center = sphericaltestShell.center
        self.shell_radius = sphericaltestShell.radius
        self.shell = self.create_new_shell(self.shell_center,
                                           self.shell_radius,
                                           domain_id)

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

#########################
class hasCylindricalShell(hasShell):
    # This class is inherited by domains that have a Cylindrical shell such as the PlanarSurfaceNonInteractionSingle
    # and a MixedPair2D3D.
    # The class makes the appropriate shell for the domain and allows a testShell to ask what the maximum dimensions
    # of the testShell can be. This means, how far can the testShell be scaled according to its scaling rules before
    # it hits THIS Cylindrical shell.

    def __init__(self, cylindricaltestShell, domain_id):
        # Note: the hasShell.__init__ is not called because nothing happens there

        assert isinstance(cylindricaltestShell, CylindricaltestShell)
        self.testShell = cylindricaltestShell

        # Compute origin, radius and half_length of cylinder.
        self.shell_center, self.shell_radius, self.shell_half_length = \
                    self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                       self.testShell.get_orientation_vector(),
                                                       self.testShell.dr,
                                                       self.testShell.dz_right,
                                                       self.testShell.dz_left)
        self.shell_center = self.testShell.world.apply_boundary(self.shell_center)
        # Make the shell.
        self.shell = self.create_new_shell(self.shell_center,
                                           self.shell_radius,
                                           self.shell_half_length,
                                           domain_id)

    @ classmethod
    def r_zright_zleft_to_r_center_hl(cls, referencepoint, orientation_vector, dr, dz_right, dz_left):
        # Converts the dr, dz_right and dz_left to the center, radius and half_length of the cylinder
        # The dr, dz_right, dz_left are used to scale the cylinder, whereas the center, radius and
        # half_length are used to use the cylinderical shell as an actual shell (and no longer as a testShell).

        center = referencepoint + ((dz_right - dz_left)/2.0) * orientation_vector
        radius = dr
        half_length = (dz_left + dz_right) / 2.0
        return center, radius, half_length

    def create_new_shell(self, position, radius, half_length, domain_id):
        orientation = self.testShell.get_orientation_vector()
        return CylindricalShell(domain_id, Cylinder(position, radius, orientation, half_length))


###########################
# functions to size up a testShell to a spherical shape object.
def get_dr_dzright_dzleft_to_SphericalShape(shape, testShell, r, z_right, z_left):
    # This function returns the dr, dz_right, dz_left parameters for the cylindrical 'testShell'
    # using the spherical 'shell' as its closest neighbor. Note that it returns the minimum the newly calculated
    # dr, dz_right, dz_left and the old r, z_right, z_left. The 'testShell' can therefore only become
    # smaller.

    # Note that the 'shell' is querried from this domain earlier by the testShell.
    # -> the shell MUST be a Sphere
    assert (type(shape) is Sphere)

    # Laurens' algorithm (part1)
    shell_position = shape.position
    shell_radius = shape.radius 

    # get the reference point and orientation of the domain to scale
    reference_point = testShell.get_referencepoint()
    orientation_vector = testShell.get_orientation_vector()

    # FIXME DIRTY HACK TO MAKE DOMAINS LESS GREEDY
#    r = min(r, testShell.world.distance(shape.position, reference_point)/2)

    # determine on what side the midpoint of the shell is relative to the reference_point
    shell_position_t = testShell.world.cyclic_transpose (shell_position, reference_point)
    ref_to_shell_vec = shell_position_t - reference_point
    ref_to_shell_z = numpy.dot(ref_to_shell_vec, orientation_vector)

    # if the shell is on the side of orientation_vector -> use z_right
    # also use this one if the shell is in plane with the reference_point
    if ref_to_shell_z >= 0:
        direction = 1           # improve direction specification such that following calculations still work
                                # Otherwise problems when direction = 0
        scale_angle = testShell.right_scalingangle
        scale_center_r, scale_center_z = testShell.right_scalingcenter 
        r1_function =  testShell.r_right
        z1_function = testShell.z_right
        z2_function = testShell.z_left
        z1 = z_right
        z2 = z_left
    # else -> use z_left
    else:
        direction = -1
        scale_angle = testShell.left_scalingangle
        scale_center_r, scale_center_z = testShell.left_scalingcenter 
        r1_function =  testShell.r_left
        z1_function = testShell.z_left
        z2_function = testShell.z_right
        z1 = z_left
        z2 = z_right

    # calculate ref_to_shell_r/z in the cylindrical coordinate system on the right/left side
    ref_to_shell_z_vec = ref_to_shell_z * orientation_vector
    ref_to_shell_r_vec = ref_to_shell_vec - ref_to_shell_z_vec
    ref_to_shell_r = length(ref_to_shell_r_vec)     # FIXME length does work well when ref_to_shell_r_vec = null vector
    ref_to_shell_z = ref_to_shell_z * direction

    # calculate the distances in r/z from the scaling center to the shell
    scale_center_to_shell_r = ref_to_shell_r - scale_center_r
    scale_center_to_shell_z = ref_to_shell_z - scale_center_z

    # calculate the angle 'shell_angle' of the vector from the 'scale_center' to the shell with the vector
    # from the 'reference_point' to the 'scale_center' (which is +- the orientation_vector).
    shell_angle = math.atan(scale_center_to_shell_r/scale_center_to_shell_z)

    if scale_center_to_shell_z < 0.0:
        shell_angle += Pi           # if the shell was too much to the side we correct the angle to be positive
    # elif: a negative angle was caused by a negative scale_center_to_shell we want a negative angle -> do nothing
    # otherwise: shell_angle is ok -> do nothing

    ### check what situation arrises
    # The shell can hit the cylinder on the flat side (situation 1),
    #                                on the edge (situation 2),
    #                                or on the round side (situation 3).
    # I think this also works for scale_angle == 0 and scale_angle == Pi/2
    if shell_angle <= 0.0:
    # The midpoint of the shell lies on the axis of the cylinder
        situation = 1

    elif 0.0 < shell_angle and shell_angle < scale_angle:
    # The spherical shell can touch the cylinder on its flat side or its edge

        if scale_angle == Pi/2.0:
            scale_center_to_shell_z_thres = 0
        else:
            # scale_angle == 0 should not get here
            scale_center_to_shell_z_thres = abs(scale_center_to_shell_r)/math.tan(scale_angle)

        shell_radius_thres = scale_center_to_shell_z - scale_center_to_shell_z_thres
        if shell_radius < shell_radius_thres:
            situation = 1
        else:
            situation = 2

    elif scale_angle <= shell_angle and shell_angle < Pi/2.0:
    # The (spherical) shell can touch the cylinder on its edge or its radial side
        tan_scale_angle = math.tan(scale_angle)                             # scale_angle == Pi/2 should not get here
        shell_radius_thres = scale_center_to_shell_r - scale_center_to_shell_z * tan_scale_angle  # TODO same here

        if shell_radius > shell_radius_thres:
            situation = 2
        else:
            situation = 3

    elif Pi/2.0 <= shell_angle:
    # The shell is always only on the radial side of the cylinder
        situation = 3

    else:
    # Don't know how we would get here, but it shouldn't happen
        raise RuntimeError('Error: shell_angle was not in valid range. shell_angle = %s, scale_angle = %s' %
                           (FORMAT_DOUBLE % shell_angle, FORMAT_DOUBLE % scale_angle))

    ### Get the right values for z and r for the given situation
    if situation == 1:      # the spherical shell hits cylinder on the flat side
        z1_new = min(z1, (ref_to_shell_z - shell_radius))
        r_new  = min(r,  r1_function(z1_new))

    elif situation == 2:    # shell hits sphere on the edge
        shell_radius_sq = shell_radius*shell_radius

        ss_sq = (scale_center_to_shell_z**2 + scale_center_to_shell_r**2)
        scale_center_to_shell = math.sqrt(ss_sq)

        angle_diff = shell_angle - scale_angle
        sin_angle_diff = math.sin(angle_diff)
        cos_angle_diff = math.cos(angle_diff)

        scale_center_shell_dist = (scale_center_to_shell * cos_angle_diff - 
                                   math.sqrt(shell_radius_sq - (ss_sq * sin_angle_diff * sin_angle_diff) ))

        if scale_angle <= Pi/4.0:
            cos_scale_angle = math.cos(scale_angle)
            z1_new = min(z1, (scale_center_z + cos_scale_angle * scale_center_shell_dist))
            r_new  = min(r, r1_function(z1_new))
        else:
            sin_scale_angle = math.sin(scale_angle)
            r_new  = min(r,  (scale_center_r + sin_scale_angle * scale_center_shell_dist))
            z1_new = min(z1, z1_function(r_new))

    elif situation == 3:    # shell hits cylinder on the round side
        r_new = min(r, (ref_to_shell_r - shell_radius))
        z1_new = min(z1, z1_function(r_new))
    else:
        raise RuntimeError('Bad situation for cylindrical shell making to sphere.')

    z2_new = min(z2, z2_function(r_new))

    # switch the z values in case it's necessary. r doesn't have to be switched.
    r = r_new
    if direction >= 0.0:
        z_right = z1_new
        z_left  = z2_new
    else:
        z_right = z2_new
        z_left  = z1_new

    return r, z_right, z_left


def get_radius_to_SphericalShape(shape, testShell, r):
    # This function returns the radius for the spherical 'testShell' using the spherical 'shell' as its closest
    # neighbor. Note that it returns the minimum the newly calculated radius and the old radius. The 'testShell' can
    # therefore only become smaller.

    # Note that the 'shell' is querried from this domain earlier by the testShell.
    # -> the shell MUST be a Sphere
    assert (type(shape) is Sphere)

    # Do simple distance calculation to sphere
    scale_point = testShell.center
    r_new = testShell.world.distance(shape, scale_point)

    # FIXME DIRTY HACK TO MAKE DOMAINS LESS GREEDY
#    r_new = min(r_new, testShell.world.distance(shape.position, scale_point)/2)

    # if the newly calculated dimensions are smaller than the current one, use them
    return min(r, r_new)

# functions to size up a testShell to a Cylindrical shape object
def get_dr_dzright_dzleft_to_CylindricalShape(shape, testShell, r, z_right, z_left):
    # This function returns the dr, dz_right, dz_left parameters for the cylindrical 'testShell'
    # using the cylindrical 'shell' as its closest neighbor. Note that it returns the minimum the newly calculated
    # dr, dz_right, dz_left and the old r, z_right, z_left. The 'testShell' can therefore only become
    # smaller.

    # Note that the 'shell' is querried from this domain earlier by the testShell.
    # -> the shell MUST be a Cylinder.
    assert (type(shape) is Cylinder)

    # Laurens' algorithm (part2)
    shell_position = shape.position
    shell_radius = shape.radius 
    shell_half_length = shape.half_length 


    # get the reference point and orientation of the domain to scale
    reference_point = testShell.get_referencepoint()
    orientation_vector = testShell.get_orientation_vector()

    # FIXME DIRTY HACK TO MAKE DOMAINS LESS GREEDY
#    r = min(r, testShell.world.distance(shape.position, reference_point)/2)

    # determine on what side the midpoint of the shell is relative to the reference_point
    shell_position_t = testShell.world.cyclic_transpose (shell_position, reference_point)
    ref_to_shell_vec = shell_position_t - reference_point
    ref_to_shell_z = numpy.dot(ref_to_shell_vec, orientation_vector)

    # if the shell is on the side of orientation_vector -> use z_right
    # also use this one if the shell is in plane with the reference_point
    if ref_to_shell_z >= 0:
        direction = 1           # improve direction specification such that following calculations still work
                                    # Otherwise problems when direction = 0
        scale_angle = testShell.right_scalingangle
        scale_center_r, scale_center_z = testShell.right_scalingcenter 
        r1_function =  testShell.r_right
        z1_function = testShell.z_right
        z2_function = testShell.z_left
        z1 = z_right
        z2 = z_left
    # else -> use z_left
    else:
        direction = -1
        scale_angle = testShell.left_scalingangle
        scale_center_r, scale_center_z = testShell.left_scalingcenter 
        r1_function =  testShell.r_left
        z1_function = testShell.z_left
        z2_function = testShell.z_right
        z1 = z_left
        z2 = z_right

    relative_orientation = abs(numpy.dot(orientation_vector, shape.unit_z))
    if feq(relative_orientation, 1.0):
    ### If the cylinders are oriented parallelly

        # calculate ref_to_shell_r/z in the cylindrical coordinate system on the right/left side
        ref_to_shell_z_vec = ref_to_shell_z * orientation_vector
        ref_to_shell_r_vec = ref_to_shell_vec - ref_to_shell_z_vec
        ref_to_shell_r = length(ref_to_shell_r_vec)         # the radius is always positive
        ref_to_shell_z = ref_to_shell_z * direction         # ref_to_shell_z is positive on the scaling side (right/left)


        # calculate the distances in r/z from the scaling center to the shell
        scale_center_to_shell_r = ref_to_shell_r - scale_center_r
        scale_center_to_shell_z = ref_to_shell_z - scale_center_z

        # get angles
        to_edge_angle = math.atan( (scale_center_to_shell_r - shell_radius) / (scale_center_to_shell_z - shell_half_length) )

        if (scale_center_to_shell_z - shell_half_length) < 0.0:
            to_edge_angle += Pi           # if the shell was too much to the side we correct the angle to be positive
        # elif: a negative angle was caused by a negative scale_center_to_shell we want a negative angle -> do nothing
        # otherwise: shell_angle is ok -> do nothing

        if to_edge_angle <= scale_angle:
            # shell hits the scaling cylinder on top
            z1_new = min(z1, (ref_to_shell_z - shell_half_length))
            r_new  = min(r,  r1_function(z1_new))           # TODO if z1 hasn't changed we also don't have to recalculate this
        else:
            # shell hits the scaling cylinder on the radial side
            r_new = min(r, (ref_to_shell_r - shell_radius))
            z1_new = min(z1, z1_function(r_new))

    elif feq(relative_orientation, 0):
    ### If the cylinders are oriented perpendicularly

        # The local x-axis is defined as the vector parallel to the axis of the shell in the direction from the
        # reference_point of the scaling cylinder to the cylindrical shell
        # The local y-axis is defined as the vector parallel to the flat side of the shell, perpendicular to the
        # axis of the scaling cylinder. It also point from the reference_point of the scaling cylinder to the shell.
        # The local z-axis is defined as the vector parallel to the axis of the scaling cylinder in the direction
        # from the reference point to the shell.
        if numpy.dot(ref_to_shell_vec, shape.unit_z) >= 0:
            local_x = shape.unit_z
        else:
            local_x = -shape.unit_z
        local_z = orientation_vector * direction
        cross = numpy.cross(local_x, local_z)
        if numpy.dot(ref_to_shell_vec, cross) >= 0:
            local_y = cross
        else:
            local_y = -cross

        ref_to_shell_x = numpy.dot(ref_to_shell_vec, local_x)
        ref_to_shell_y = numpy.dot(ref_to_shell_vec, local_y)
        ref_to_shell_z = numpy.dot(ref_to_shell_vec, local_z)
        assert (ref_to_shell_x >= 0.0) and (ref_to_shell_y >= 0.0) and (ref_to_shell_z >= 0.0)  # By the definition of the local coordinate system

        scale_center_to_shell_x = ref_to_shell_x
        scale_center_to_shell_y = ref_to_shell_y
        scale_center_to_shell_z = ref_to_shell_z - scale_center_z
        
        ref_to_shell_x2 = ref_to_shell_x - shell_half_length
        ref_to_shell_y2 = ref_to_shell_y - shell_radius

        if (ref_to_shell_x2 < 0) and (ref_to_shell_y2 < 0):
            # quadrant 1
#            print "quadrant 1"
            situation = 4

        elif (ref_to_shell_x2 >= 0) and (ref_to_shell_y2 < 0):
            # quadrant 2
#            print "quadrant 2"
            scale_center_to_shell_x -= scale_center_r
            scale_center_to_flatend_x = scale_center_to_shell_x - shell_half_length
            scale_center_to_critpoint_z = scale_center_to_shell_z - math.sqrt(shell_radius**2 - scale_center_to_shell_y**2) 
            scale_center_to_lowpoint_z = scale_center_to_shell_z - shell_radius

            # angle1
            if scale_center_to_critpoint_z == 0:
                scale_center_to_shell_crit_angle_y = Pi/2
            else:
                scale_center_to_shell_crit_angle_y = math.atan(scale_center_to_flatend_x/ scale_center_to_critpoint_z)
                if scale_center_to_critpoint_z < 0.0:
                    scale_center_to_shell_crit_angle_y += Pi
            # angle2
            if scale_center_to_lowpoint_z == 0:
                scale_center_to_shell_low_angle_y = Pi/2
            else:
                scale_center_to_shell_low_angle_y  = math.atan(scale_center_to_flatend_x/ scale_center_to_lowpoint_z )
                if scale_center_to_lowpoint_z < 0.0:
                    scale_center_to_shell_low_angle_y += Pi


            if scale_angle <= scale_center_to_shell_crit_angle_y:
                situation = 1
            elif scale_center_to_shell_low_angle_y <= scale_angle:
                situation = 4
            else:
                assert scale_center_to_shell_crit_angle_y < scale_angle and \
                       scale_angle < scale_center_to_shell_low_angle_y
                situation = 2
                r1_min = (scale_center_to_shell_x - shell_half_length)*(1.0+TOLERANCE)
                h1_min = r1_min/math.tan(scale_angle)

        elif (ref_to_shell_x2 < 0) and (ref_to_shell_y2 >= 0):
            # quadrant 3
#            print "quadrant 3"
            scale_center_to_shell_y -= scale_center_r
            scale_center_to_critpoint_y = scale_center_to_shell_y - shell_radius
            scale_center_to_lowpoint_z = scale_center_to_shell_z - shell_radius

            # angle1
            if scale_center_to_shell_z == 0.0:
                scale_center_to_shell_crit_angle_x = Pi/2.0
            else:
                scale_center_to_shell_crit_angle_x = math.atan(scale_center_to_critpoint_y / scale_center_to_shell_z)
                if scale_center_to_shell_z < 0.0:
                    scale_center_to_shell_crit_angle_x += Pi
            # angle2
            if scale_center_to_lowpoint_z == 0.0:
                scale_center_to_shell_low_angle_x = Pi/2.0
            else:
                scale_center_to_shell_low_angle_x  = math.atan(scale_center_to_shell_y/scale_center_to_lowpoint_z)
                if scale_center_to_lowpoint_z < 0.0:
                    scale_center_to_shell_low_angle_x += Pi


            if scale_angle <= scale_center_to_shell_crit_angle_x:
                situation = 6
            elif scale_center_to_shell_low_angle_x <= scale_angle:
                situation = 4
            else:
                assert scale_center_to_shell_crit_angle_x < scale_angle and \
                       scale_angle < scale_center_to_shell_low_angle_x
                situation = 5

        else:
            # quadrant 4
#            print "quadrant 4"
            assert (ref_to_shell_x2 >= 0) and (ref_to_shell_y2 >= 0)

            scale_center_to_critpoint_r = math.sqrt((scale_center_to_shell_y - shell_radius)**2 + (scale_center_to_shell_x - shell_half_length)**2) - scale_center_r
            scale_center_to_lowpoint_r = math.sqrt(scale_center_to_shell_y**2 + (scale_center_to_shell_x - shell_half_length)**2) - scale_center_r
            scale_center_to_lowpoint_z = scale_center_to_shell_z - shell_radius

            # angle1
            if scale_center_to_shell_z == 0.0:
                scale_center_to_shell_crit_angle_xy = Pi/2.0
            else:
                scale_center_to_shell_crit_angle_xy = math.atan(scale_center_to_critpoint_r / scale_center_to_shell_z)
                if scale_center_to_shell_z < 0.0:
                    scale_center_to_shell_crit_angle_xy += Pi
            # angle2
            if scale_center_to_lowpoint_z == 0.0:
                scale_center_to_shell_low_angle_xy = Pi/2.0
            else:
                scale_center_to_shell_low_angle_xy  = math.atan(scale_center_to_lowpoint_r / scale_center_to_lowpoint_z )
                if scale_center_to_lowpoint_z < 0.0:
                    scale_center_to_shell_low_angle_xy += Pi


            if scale_angle <= scale_center_to_shell_crit_angle_xy:
                situation = 3
            elif scale_center_to_shell_low_angle_xy <= scale_angle:
                situation = 4
            else:
                assert scale_center_to_shell_crit_angle_xy < scale_angle and \
                       scale_angle < scale_center_to_shell_low_angle_xy
                situation = 2
                r1_min = math.sqrt((scale_center_to_shell_x-shell_half_length)**2 + (scale_center_to_shell_y-shell_radius)**2)*(1.0+TOLERANCE)
                h1_min = r1_min/math.tan(scale_angle)

#        print "situation= ", situation
        #################
        if situation == 1:
            # shell hits the scaling cylinder with its flat surface on the radial side
            r_new = min(r, (ref_to_shell_x - shell_half_length))
            z1_new = min(z1, z1_function(r_new))

        elif situation == 2:
            # shell hits the scaling cylinder with its edge on the edge
            # TODO we have a solution but it can only be found with a root finder -> slow
            tan_scale_angle = math.tan(scale_angle)

            if scale_angle <= Pi/4.0:
                def h1(x):
                    return x - scale_center_to_shell_z + \
                           math.sqrt(shell_radius**2 - (scale_center_to_shell_y - math.sqrt((x*tan_scale_angle)**2 - (scale_center_to_shell_x - shell_half_length)**2) )**2 )

                h_touch = scale_center_z + findroot(h1, h1_min, scale_center_to_shell_z*(1.0-TOLERANCE))
                z1_new = min(z1, h_touch)
                r_new  = min(r,  r1_function(z1_new))
            else:
                def r1(x):
                    shell_radius_sq = shell_radius**2
                    return x/tan_scale_angle - scale_center_to_shell_z + \
                           math.sqrt(shell_radius_sq - (scale_center_to_shell_y - math.sqrt(x**2 - (scale_center_to_shell_x - shell_half_length)**2))**2)

                r_touch = findroot(r1, r1_min, math.sqrt((scale_center_to_shell_x-shell_half_length)**2 + (scale_center_to_shell_y)**2)*(1.0-TOLERANCE))
                r_new  = min(r, r_touch)
                z1_new = min(z1, z1_function(r_new))

        elif situation == 3:
            # shell hits the scaling cylinder with its edge on the radial side
            r_new = min(r, math.sqrt( (ref_to_shell_x-shell_half_length)**2 + (ref_to_shell_y - shell_radius)**2))
            z1_new = min(z1, z1_function(r_new))

        elif situation == 4:
            # shell hits the scaling cylinder with its round side on the top flat side
            z1_new = min(z1, (ref_to_shell_z - shell_radius))
            r_new  = min(r,  r1_function(z1_new))

        elif situation == 5:
            # shell hist the scaling cylinder with its round side on the edge of the scaling cylinder.
            shell_radius_sq = shell_radius*shell_radius

            ss_sq = (scale_center_to_shell_z**2 + scale_center_to_shell_y**2)
            scale_center_to_shell = math.sqrt(ss_sq)

            shell_angle_x = math.atan(scale_center_to_shell_y/scale_center_to_shell_z)
            angle_diff = shell_angle_x - scale_angle
            sin_angle_diff = math.sin(angle_diff)
            cos_angle_diff = math.cos(angle_diff)

            scale_center_shell_dist = (scale_center_to_shell * cos_angle_diff -
                                       math.sqrt(shell_radius_sq - (ss_sq * sin_angle_diff * sin_angle_diff) ))

            if scale_angle <= Pi/4.0:
                cos_scale_angle = math.cos(scale_angle)
                z1_new = min(z1, (scale_center_z + cos_scale_angle * scale_center_shell_dist))
                r_new  = min(r, r1_function(z1_new))
            else:
                sin_scale_angle = math.sin(scale_angle)
                r_new  = min(r,  (scale_center_r + sin_scale_angle * scale_center_shell_dist))
                z1_new = min(z1, z1_function(r_new))

        elif situation == 6:
            # shell hits the scaling cylinder with its round side on the radial side
            r_new = min(r, (ref_to_shell_y - shell_radius))
            z1_new = min(z1, z1_function(r_new))

        else:
            raise RuntimeError('get_dr_dzright_dzleft_to_CylindricalShape: Bad situation for making cylinders against cylinders.')

    else:
        raise NotImplementedError('get_dr_dzright_dzleft_to_CylindricalShape: Cylinders should be oriented parallel or perpendicular.')

    z2_new = min(z2, z2_function(r_new))

    # switch the z values in case it's necessary. r doesn't have to be switched.
    r = r_new
    if direction == 1:
        z_right = z1_new
        z_left  = z2_new
    else:
        z_right = z2_new
        z_left  = z1_new

    return r, z_right, z_left

def get_radius_to_CylindricalShape(shape, testShell, r):
    # This function returns the radius for the spherical 'testShell' using the cylindrical 'shell' as its closest
    # neighbor. Note that it returns the minimum the newly calculated radius and the old radius. The 'testShell' can
    # therefore only become smaller.

    # Note that the 'shell' is querried from this domain earlier by the testShell.
    # -> the shell MUST be a Cylinder
    assert (type(shape) is Cylinder)

    # Do simple distance calculation to cylinder
    scale_point = testShell.center
    r_new = testShell.world.distance(shape, scale_point)

    # if the newly calculated dimensions are smaller than the current one, use them
    return min(r, r_new)

def get_dr_dzright_dzleft_to_PlanarShape(shape, testShell, r, z_right, z_left):
    # This function returns the dr, dz_right, dz_left for the cylindrical 'testShell' using the planar 'shape' as its closest
    # neighbor. Note that it returns the minimum of the newly calculated parameters and the old parameters. The 'testShell' can
    # therefore only become smaller.

    assert (type(shape) is Plane)

    plane_position = shape.position

    # get the reference point and orientation of the domain to scale
    reference_point = testShell.get_referencepoint()
    orientation_vector = testShell.get_orientation_vector()

    # determine what parameter of the cylinder (dr of dz) to scale.
    relative_orientation = numpy.dot(orientation_vector, shape.unit_z)
    if feq(relative_orientation, 0):

        distance = testShell.world.distance(shape, reference_point)
        # FIXME ugly hack to make planar surface singles overlap with membranes
        if isinstance(testShell, PlanarSurfaceSingletestShell):
            distance += testShell.pid_particle_pair[1].radius

        # The unit_z of the plane is perpendicular to the cylinder -> scale r
        r = min(r, distance)
        z_right = min(z_right, testShell.z_right(r))
        z_left  = min(z_left,  testShell.z_left (r))

    elif feq(abs(relative_orientation), 1):
        # The unit_z of the plane is parallel to the axis of the cylinder -> scale z_right/z_left

        # determine on what side the midpoint of the plane is relative to the reference_point
        plane_position_t = testShell.world.cyclic_transpose (plane_position, reference_point)
        ref_to_plane_vec = plane_position_t - reference_point
        ref_to_plane_z = numpy.dot(ref_to_plane_vec, orientation_vector)

        # if the shell is on the side of orientation_vector -> use z_right
        if ref_to_plane_z >= 0:
            z_right = min(z_right, testShell.world.distance(shape, reference_point))
            r       = min(r, testShell.r_right(z_right))
            z_left  = min(z_left, testShell.z_left(r))
        else:
            z_left  = min(z_left, testShell.world.distance(shape, reference_point))
            r       = min(r, testShell.r_left(z_left))
            z_right = min(z_right, testShell.z_right(r))
            
    else:
        raise NotImplementedError('Only perpendicular planes are supported.')

    return r, z_right, z_left


def get_radius_to_PlanarShape(shape, testShell, r):
    # This function returns the radius for the spherical 'testShell' using the planar 'shape' as its closest
    # neighbor. Note that it returns the minimum of the newly calculated radius and the old radius. The 'testShell' can
    # therefore only become smaller.

    assert (type(shape) is Plane)

    # Do simple distance calculation to sphere
    scale_point = testShell.center
    r_new = testShell.world.distance(shape, scale_point)

    # if the newly calculated dimensions are smaller than the current one, use them
    return min(r, r_new)


#######################################################
#######################################################
class testShell(object):
    # The testShell represents the geometric aspect of the domain that we are trying to make
    # A testShell can be Spherical (SphericaltestShell) or Cylindrical (CylindricaltestShell)
    # Both have different rules for doing the nearest neighbor search and geometrically sizing up
    #  the prospective domain

    def __init__(self, geometrycontainer, domains):
        self.geometrycontainer = geometrycontainer
        self.world = geometrycontainer.world
        self.domains = domains

    def get_neighbors(self, base_structure_id, ignore, ignores):
        # The get_neighbors is general to all types of testShells.
        # It finds the neighboring domains and surfaces.
        searchpoint = self.get_searchpoint()
        radius = self.get_searchradius()

        neighbor_domains  = self.geometrycontainer.get_neighbor_domains(searchpoint, self.domains, ignore)
        neighbor_surfaces = self.geometrycontainer.get_neighbor_surfaces(searchpoint, base_structure_id, ignores)
        return neighbor_domains, neighbor_surfaces

    def get_orientation_vector(self):
        # The orientation of the shell
        pass    # overridden later

########################
class SphericaltestShell(testShell):
    # The SphericaltestShell is the geometric aspect of testDomains with a spherical shell.

    def __init__(self, geometrycontainer, domains):
        testShell.__init__(self, geometrycontainer, domains)

        self.center = None          # These are the parameters that define the dimensions of
        self.radius = 0             # the spherical shell -> will be determined later.

    def determine_possible_shell(self, base_structure_id, ignore, ignores):
        # This determines the largest possible radius of the spherical testShell or throws an
        # exception if the domain could not be made due to geometrical constraints.

        neighbor_domains, neighbor_surfaces = self.get_neighbors (base_structure_id, ignore, ignores)
        min_radius = self.apply_safety(self.get_min_radius())
        max_radius = self.get_max_radius()

        # initialize the radius to start the scaling of the sphere
        radius = max_radius

        # first check the maximum radius against the surfaces
        for surface, distance in neighbor_surfaces:
            # FIXME ugly hack to make spherical singles overlap with membranes
            if isinstance(surface, PlanarSurface) and isinstance(self, SphericalSingletestShell):
                distance += self.pid_particle_pair[1].radius
            radius = min(radius, distance)
            if radius < min_radius:
                raise ShellmakingError('Surface too close to make spherical testshell, '
                                       'surface = %s, distance = %s, testShell = %s' %
                                       (surface, distance, self))

        # TODO first sort the neighbors -> faster to find if we fail

        # then check against all the shells of the neighboring domains (remember that domains can have
        # multiple shells).
        for neighbor, _ in neighbor_domains:

            shell_list = self.get_neighbor_shell_list(neighbor)
            for _, shell_appearance in shell_list:
                if isinstance(shell_appearance.shape, Sphere):
                    radius = get_radius_to_SphericalShape(shell_appearance.shape, self, radius)
                elif isinstance(shell_appearance.shape, Cylinder):
                    radius = get_radius_to_CylindricalShape(shell_appearance.shape, self, radius)

                if radius < min_radius:
                    raise ShellmakingError('Domain too close to make spherical testshell, '
                                           'domain = %s, radius = %s, min_radius = %s, testShell = %s' %
                                           (neighbor, radius, min_radius, self))

        # we calculated a valid radius -> suscess!
        assert radius >= min_radius, 'SphericaltestShell radius smaller than the minimum, radius = %s, min_radius = %s.' % \
                                     (radius, min_radius)

        radius = self.apply_safety(radius)
        return radius

    def get_searchpoint(self):
        # The search point from where to search for neighboring domains.
        return self.center
    def get_searchradius(self):
        # the search radius is the radius in which to find neighboring domains. Usually just the
        # max shell size.
        return self.geometrycontainer.get_max_shell_size()
    def get_orientation_vector(self):
        # The orientation of the sphere is irrelevant, but needs to be defined.
        return self.structure.shape.unit_z    # NOTE: Is self.structure always well defined??

    def get_max_radius(self):
        # The maximum radius of the spherical testShell.
        return self.get_searchradius()

    def apply_safety(self, r):
        return r/SAFETY
#####
class SphericalSingletestShell(SphericaltestShell, testNonInteractionSingle):

    def __init__(self, pid_particle_pair, structure, geometrycontainer, domains):
        SphericaltestShell.__init__(self, geometrycontainer, domains)
        testNonInteractionSingle.__init__(self, pid_particle_pair, structure)

        self.center = self.pid_particle_pair[1].position
        self.radius = self.pid_particle_pair[1].radius
        # Here determine_possible_shell is not called since the NonInteractionSingle should never fail

    def get_min_radius(self):
        return self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR   # the minimum radius of a NonInteractionSingle

#####
class SphericalPairtestShell(SphericaltestShell, testSimplePair):

    def __init__(self, single1, single2, geometrycontainer, domains):
        SphericaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testSimplePair.__init__(self, single1, single2)

        self.center = self.com
        try:
            self.radius = self.determine_possible_shell(self.structure.id, [self.single1.domain_id, self.single2.domain_id],
                                                        [])
        except ShellmakingError as e:
            raise testShellError('(SphericalPair). %s' %
                                 (str(e)))

    def get_min_radius(self):
        return self.get_min_pair_size()

#####
class PlanarSurfaceTransitionSingletestShell(SphericaltestShell, testTransitionSingle):

    def __init__(self, single, target_structure, geometrycontainer, domains):
        SphericaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testTransitionSingle.__init__(self, single, target_structure)

        assert isinstance(self.origin_structure, PlanarSurface)
        assert isinstance(self.target_structure, PlanarSurface)

        self.center = self.pid_particle_pair[1].position
        try:
            self.radius = self.determine_possible_shell(self.structure.id, [self.single.domain_id],
                                                        [self.target_structure.id])
        except ShellmakingError as e:
            raise testShellError('(PlanarSurfaceTransitionSingle). %s' %
                                 (str(e)))

    def get_min_radius(self):
        # The minimal shell size also deterimines the minimal radius of the later 2D domain.
        # min_offset here is the minimal length that the circular 2D domain protrudes into
        # the target surface. Make sure that min_offset > 0.0, because only then it is
        # guaranteed that there is always a certain probability to cross the edge even
        # if the shell with minimal radius is constructed
        return self.distance_to_target_structure + self.pid_particle_pair[1].radius


##########################
class CylindricaltestShell(testShell):
    # The CylindricaltestShell is the geometric aspect of testDomains with a cylindrical shell.
    # Note that the dimensions of the cylindrical test shell are expressed in dr, dz_right and
    # dz_left, whereas the dimensions of a 'real' cylindrical shell are expressed in center,
    # radius and half_length.

    def __init__(self, geometrycontainer, domains):
        testShell.__init__(self, geometrycontainer, domains)

        self.dr       = 0           # These are the parameters that define the dimensions of
        self.dz_right = 0           # the cylindrical shell -> will be determined later.
        self.dz_left  = 0
        # dr is the radius of the cylinder.
        # dz_right is the distance between the 'reference point' and the face of the
        #   cylinder on the side of the 'orientation vector'.
        # dz_left is the distance from the 'reference point' to the face of the cylinder
        #   on the side oposite of the 'orientation vector'.


    def determine_possible_shell(self, base_structure_id, ignore, ignores):
        # This determines the maximum dr, dz_right, dz_left of the cylindrical testShell or
        # throws an exception if the domain could not be made due to geometrical constraints.

        neighbor_domains, neighbor_surfaces = self.get_neighbors (base_structure_id, ignore, ignores)
        min_dr, min_dz_right, min_dz_left = self.get_min_dr_dzright_dzleft()
        max_dr, max_dz_right, max_dz_left = self.get_max_dr_dzright_dzleft()

        # initialize the parameters to start the scaling of the cylinder
        min_dr, min_dz_right, min_dz_left = self.apply_safety(min_dr, min_dz_right, min_dz_left)
        dr, dz_right, dz_left             = max_dr, max_dz_right, max_dz_left
 
        # first check the maximum dr, dz_right, dz_left against the surfaces
        # NOTE: we assume that all relevant surfaces are cylindrical and parallel to the testCylinder
        # or planar surfaces and parallel to the testCylinder axis
        for surface, distance in neighbor_surfaces:
            # TODO
            if isinstance(surface, CylindricalSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_CylindricalShape(surface.shape, self,
                                                                                              dr, dz_right, dz_left)
            elif isinstance(surface, PlanarSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_PlanarShape(surface.shape, self,
                                                                                         dr, dz_right, dz_left)
            else:
                assert False, "Wrong type of surface used"

            if (dr < min_dr) or (dz_right < min_dz_right) or (dz_left < min_dz_left):
                raise ShellmakingError('Surface too close to make cylindrical testshell, '
                                       'surface = %s, dr = %s, dz_right = %s, dz_left = %s, testShell = %s, min_dr = %s' %
                                       (surface, dr, dz_right, dz_left, self, min_dr))
        # TODO first sort the neighbors to distance -> faster to find if we fail

        # then check against all the shells of the neighboring domains (remember that domains can have
        # multiple shells).
        for neighbor, distance in neighbor_domains:

            shell_list = self.get_neighbor_shell_list(neighbor)
            for _, shell_appearance in shell_list:
                if isinstance(shell_appearance.shape, Sphere):
                    dr, dz_right, dz_left = get_dr_dzright_dzleft_to_SphericalShape(shell_appearance.shape, self,
                                                                                                dr, dz_right, dz_left)
                elif isinstance(shell_appearance.shape, Cylinder):
                    dr, dz_right, dz_left = get_dr_dzright_dzleft_to_CylindricalShape(shell_appearance.shape, self,
                                                                                                  dr, dz_right, dz_left)
                else:
                    assert False, "Wrong type of shell shape"

                if (dr < min_dr) or (dz_right < min_dz_right) or (dz_left < min_dz_left):
                    raise ShellmakingError('Domain too close to make cylindrical testshell, '
                                           'domain = %s, dr = %s, dz_right = %s, dz_left = %s, min_dr = %s, min_dz_right = %s, min_dz_left = %s.' % \
                                           (neighbor, dr, dz_right, dz_left, min_dr, min_dz_right, min_dz_left))

        # we calculated valid dimensions -> success!
        assert dr >= min_dr and dz_right >= min_dz_right and dz_left >= min_dz_left, \
               'CylindricaltestShell dimensions smaller than the minimum, \
                dr = %s, dz_right = %s, dz_left = %s, min_dr = %s, min_dz_right = %s, min_dz_left = %s.' % \
               (dr, dz_right, dz_left, min_dr, min_dz_right, min_dz_left)

        # Note that dr, dz_right, dz_left can now actually be smaller than the minimum
        dr, dz_right, dz_left = self.apply_safety(dr, dz_right, dz_left)

        return dr, dz_right, dz_left

    def get_searchpoint(self):
        # The search point from where to search for neighboring domains.
        pass

    def get_searchradius(self):
        # the search radius is the radius in which to find neighboring domains. Usually just the
        # max shell size.
        return self.geometrycontainer.get_max_shell_size()

    def get_orientation_vector(self):
        # The orientation_vector just defines the direction of the shell. It also defines what
        # SIDE is 'right' and what side is 'left'.
        pass

    def get_referencepoint(self):
        # The reference point is where the 'orientation_vector' starts and defines what PART of the
        # cylinder is 'right' and what part if 'left'.
        pass

    # These methods are used to calculate the new r/z_right/z_left after one of the parameters 
    # r/z_right/z_left has changed. The scaling centers (r0, z0) and scaling directions drdz
    # are defined differently for every testShell.
    def r_right(self, z_right):
        return self.drdz_right * (z_right - self.z0_right) + self.r0_right
    def z_right(self, r_right):
        return self.dzdr_right * (r_right - self.r0_right) + self.z0_right
    def r_left(self, z_left):
        return self.drdz_left * (z_left - self.z0_left) + self.r0_left
    def z_left(self, r_left):
        return self.dzdr_left * (r_left - self.r0_left) + self.z0_left

    def get_right_scalingcenter(self):
        # Returns the location of the scaling center on the 'right side' of the cylinder in the
        # coordinates r, z defined by the plane cutting through the axis of the cylindrical
        # testShell and the center of the shell that is being scaled to.
        # Note0: the 'right side' is the side of the cylinder were the 'orientation_vector' points.
        # Note1: we assume cylindrical symmetry around the axis.
        # Note2: the coordinates are relative to the 'side' (right/left) of the cylinder.
        return self.r0_right, self.z0_right
    right_scalingcenter = property(get_right_scalingcenter)

    def get_left_scalingcenter(self):
        # Returns the location of the scaling center on the 'left side' of the cylinder.
        # Note: the 'left side' is the side of the cylinder opposite of were the 'orientation_vector' points.
        return self.r0_left, self.z0_left
    left_scalingcenter = property(get_left_scalingcenter)

    def get_right_scalingangle(self):
        # The scaling angle is the angle with the axis of the cylinder of the cone along which the
        # edge of the cylinder moves while scaling.
        # Note again that 'right' is the side where the 'orientation_vector' points.
        if self.drdz_right == numpy.inf:
            right_scaling_angle = Pi/2.0       # atan(infinity) == Pi/2 -> angle is 90 degrees
        elif self.drdz_right == 0.0:
            right_scaling_angle = 0.0          # atan(0.0) = 0.0 -> faster
        else:
            right_scaling_angle = math.atan(self.drdz_right)
        return right_scaling_angle
    right_scalingangle = property(get_right_scalingangle)

    def get_left_scalingangle(self):
        if self.drdz_left == numpy.inf:
            left_scaling_angle = Pi/2.0       # atan(infinity) == Pi/2 -> angle is 90 degrees
        elif self.drdz_left == 0.0:
            left_scaling_angle = 0.0          # atan(0.0) = 0.0 -> faster
        else:
            left_scaling_angle = math.atan(self.drdz_left)
        return left_scaling_angle
    left_scalingangle = property(get_left_scalingangle)

    def get_min_dr_dzright_dzleft(self):
        # This defines the minimal dimensions in terms of r, z_right, z_left of the cylindrical domain
        # that we're trying to make. The exact values vary between type of domain.
        # If a nearby domain would make it smaller than this, an exception will be thrown.

        # NOTE the minimum domain also needs to obey the scaling rules or it might be scaled down to a size
        # smaller than the minimum
        pass

    def get_max_dr_dzright_dzleft(self):
        # This defines the maximal dimensions in terms of r, z_right, z_left of the cylindrical domain.
        # The cylindrical domain needs to fit inside the search_radius (centered around the search_point)
        # because 'it doesn't know what is outside of the searchradius'.

        # NOTE also that the maximum cylinder needs to obey the scaling rules or it might be scaled up to
        # a size larger than the search_radius.
        pass

#####
class PlanarSurfaceSingletestShell(CylindricaltestShell, testNonInteractionSingle):

    def __init__(self, pid_particle_pair, structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)
        testNonInteractionSingle.__init__(self, pid_particle_pair, structure)

        # initialize the scaling parameters
        self.dzdr_right = 0.0
        self.drdz_right = numpy.inf
        self.r0_right   = 0.0
        self.z0_right   = self.pid_particle_pair[1].radius
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.pid_particle_pair[1].radius

        # sizing up the shell to a zero shell
        self.dz_right = self.pid_particle_pair[1].radius
        self.dz_left  = self.pid_particle_pair[1].radius
        self.dr       = self.pid_particle_pair[1].radius
        # Here determine_possible_shell is not called since the making of a NonInteractionSingle
        # should never fail

    def get_orientation_vector(self):
        return self.structure.shape.unit_z   # just copy from structure

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.pid_particle_pair[1].position

    def get_min_dr_dzright_dzleft(self):
        # we need to make sure that the minimal cylinder for the Single fits inside the putative
        # spherical Multi shell.
        # NOTE the MULTI_SHELL_FACTOR should therefore be at least sqrt(2)!!
        dr       = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        dz_right = self.pid_particle_pair[1].radius
        dz_left  = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dz_right = self.pid_particle_pair[1].radius
        dz_left = dz_right
        dr = math.sqrt((self.get_searchradius())**2 - dz_right**2) # stay within the searchradius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r/SAFETY, z_right, z_left

#####
class PlanarSurfacePairtestShell(CylindricaltestShell, testSimplePair):

    def __init__(self, single1, single2, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testSimplePair.__init__(self, single1, single2)

        # initialize the scaling parameters
        self.dzdr_right = 0.0
        self.drdz_right = numpy.inf
        self.r0_right   = 0.0
        self.z0_right   = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.z0_right

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single1.domain_id, self.single2.domain_id],
                                                          [])
        except ShellmakingError as e:
            raise testShellError('(PlanarSurfacePair). %s' %
                                 (str(e)))

    def get_orientation_vector(self):
        return self.structure.shape.unit_z

    def get_searchpoint(self):
        return self.com

    def get_referencepoint(self):
        return self.com

    def get_min_dr_dzright_dzleft(self):
        dr = self.get_min_pair_size()
        dz_right = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        dz_left = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dz_right = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        dz_left = dz_right
        dr = math.sqrt((self.get_searchradius())**2 - dz_right**2) # stay within the searchradius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r/SAFETY, z_right, z_left

#####
class CylindricalSurfaceSingletestShell(CylindricaltestShell, testNonInteractionSingle):

    def __init__(self, pid_particle_pair, structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)
        testNonInteractionSingle.__init__(self, pid_particle_pair, structure)

        # initialize the scaling parameters
        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.pid_particle_pair[1].radius
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.pid_particle_pair[1].radius
        self.z0_left    = 0.0


        # sizing up the shell to a zero shell
        self.dz_right = self.pid_particle_pair[1].radius
        self.dz_left  = self.pid_particle_pair[1].radius
        self.dr       = self.pid_particle_pair[1].radius

    def get_orientation_vector(self):
        return self.structure.shape.unit_z   # just copy from structure

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.pid_particle_pair[1].position

    def get_min_dr_dzright_dzleft(self):
        # we need to make sure that the minimal cylinder for the Single fits inside the putative
        # spherical Multi shell.
        # NOTE the MULTI_SHELL_FACTOR should therefore be at least sqrt(2)!!
        dr       = self.pid_particle_pair[1].radius
        dz_right = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        dz_left  = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dr = self.pid_particle_pair[1].radius
        dz_right = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
        dz_left = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class CylindricalSurfacePairtestShell(CylindricaltestShell, testSimplePair):

    def __init__(self, single1, single2, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testSimplePair.__init__(self, single1, single2)

        # initialize scaling parameters
        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        self.z0_left    = 0.0

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single1.domain_id, self.single2.domain_id],
                                                          [])
        except ShellmakingError as e:
            raise testShellError('(CylindricalSurfacePair). %s' %
                                 (str(e)))

        # make sure the domain is symmetric around the CoM
        self.dz_right = min(self.dz_right, self.dz_left)
        self.dz_left  = self.dz_right

    def get_orientation_vector(self):
        return self.structure.shape.unit_z

    def get_searchpoint(self):
        return self.com

    def get_referencepoint(self):
        return self.com

    def get_min_dr_dzright_dzleft(self):
        dr = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        dz_right = self.get_min_pair_size()
        dz_left = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dr = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius)
        dz_right = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
        dz_left = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class PlanarSurfaceInteractiontestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, surface, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, surface)

        # initialize scaling parameters
        self.dzdr_right = 1.0
        self.drdz_right = 1.0
        self.r0_right   = 0.0
        self.z0_right   = self.particle_surface_distance # calculated in the __init__ of testInteractionSingle
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.pid_particle_pair[1].radius

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single.domain_id], [self.surface.id])
        except ShellmakingError as e:
            raise testShellError('(PlanarSurfaceInteration). %s' %
                                 (str(e)))


    def get_orientation_vector(self):
        return self.reference_vector      # calculated in the __init__ of testInteractionSingle

    def get_searchpoint(self):
#        return self.pid_particle_pair[1].position   # since the shell is scaled equally in all directions
        # TODO this can be improved
        search_distance = self.get_searchradius()/math.sqrt(2) - self.z0_left
        return self.get_referencepoint() + search_distance * self.get_orientation_vector()

    def get_referencepoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_min_dr_dzright_dzleft(self):
        dr       = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        dz_right = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR + self.particle_surface_distance
        dz_left  = self.pid_particle_pair[1].radius
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        # TODO this can be improved
        dr       = (self.get_searchradius()/math.sqrt(2))
        dz_right = dr + self.particle_surface_distance          # Note that we assume that dr >> particle_surface_distance
        dz_left  = self.pid_particle_pair[1].radius             # so that the cylinder always fits inside the search_radius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return self.r_right(z_right/SAFETY), z_right/SAFETY, z_left

#####
class CylindricalSurfaceInteractiontestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, surface, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, surface)

        # initialize scaling parameters
        self.dzdr_right = 1.0
        self.drdz_right = 1.0
        self.r0_right   = 0.0
        self.z0_right   = -(self.particle_surface_distance + self.surface.shape.radius)   # note that this is because dzdr_right = 1.0
        self.dzdr_left  = self.dzdr_right
        self.drdz_left  = self.drdz_right
        self.r0_left    = self.r0_right
        self.z0_left    = self.z0_right

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single.domain_id], [self.surface.id])
        except ShellmakingError as e:
            raise testShellError('(CylindricalSurfaceInteraction). %s' %
                                 (str(e)))


    def get_orientation_vector(self):
        return self.surface.shape.unit_z

    def get_searchpoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_referencepoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_min_dr_dzright_dzleft(self):
        dr       = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR + (self.particle_surface_distance + self.surface.shape.radius)
        dz_right = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        dz_left  = dz_right
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        # TODO clarify
        dz_right = ((-(self.particle_surface_distance + self.surface.shape.radius) + \
                     math.sqrt(2*(self.get_searchradius()**2) - (self.particle_surface_distance + self.surface.shape.radius)**2)) / 2.0)
        dz_left  = dz_right
        dr       = ((self.particle_surface_distance + self.surface.shape.radius) + dz_right)
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r/SAFETY, self.z_right(r/SAFETY), self.z_left(r/SAFETY)

class CylindricalSurfaceSinktestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, surface, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, surface)

        # initialize scaling parameters
        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.pid_particle_pair[1].radius
        self.z0_right   = 0.0
        self.dzdr_left  = self.dzdr_right
        self.drdz_left  = self.drdz_right
        self.r0_left    = self.r0_right
        self.z0_left    = self.z0_right

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single.domain_id], [self.surface.id])
        except ShellmakingError as e:
            raise testShellError('(CylindricalSurfaceSink). %s' %
                                 (str(e)))

        # make sure the domain is symmetric around the particle
        self.dz_right = min(self.dz_right, self.dz_left)
        self.dz_left  = self.dz_right
        # check if the domain is still large enough for the interaction
        _, _, min_dz_left = self.get_min_dr_dzright_dzleft()
        if self.dz_left < min_dz_left:
            raise testShellError('(CylindricalSurfaceSink). Domain too small after symmetricalizing.')


    def get_orientation_vector(self):
        pos = self.world.cyclic_transpose(self.pid_particle_pair[1].position, self.surface.shape.position)
        orientation_vector = pos - self.surface.shape.position
        # if the dot product is negative, then take the -unit_z, otherwise take unit_z
        return self.surface.shape.unit_z * cmp(numpy.dot(orientation_vector, self.surface.shape.unit_z), 0)

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.pid_particle_pair[1].position

    def get_min_dr_dzright_dzleft(self):
        dr       = self.pid_particle_pair[1].radius
        dz_left  = self.particle_surface_distance + self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        dz_right = dz_left              # we make the domain to be symmetrical around the particle
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        dr       = self.pid_particle_pair[1].radius
        dz_left  = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
        dz_right = dz_left
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class MixedPair2D3DtestShell(CylindricaltestShell, testMixedPair2D3D):

    def __init__(self, single2D, single3D, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        # The initialization of r0 can fail in testPair.__init__

#        raise testShellError('(MixedPair2D3D).')

        try:
            testMixedPair2D3D.__init__(self, single2D, single3D)
        except protoDomainError as e:
            raise testShellError('(MixedPair2D3D). %s' %
                                 (str(e)))

        # initialize commonly used constants
        self.sqrt_DRDr = math.sqrt((2*self.D_R)/(3*self.D_r))

        # initialize the scaling parameters
        bla = self.get_scaling_factor() * ( self.sqrt_DRDr + max( self.pid_particle_pair1[1].D/self.D_tot,
                                                                  self.pid_particle_pair2[1].D/self.D_tot ))
        minimum = self.r0+self.pid_particle_pair1[1].radius
        self.drdz_right = self.r_right(minimum + 1.0) - self.r_right(minimum)
        assert feq(bla, self.drdz_right), 'bla= %s, drdz_right=%s' % (bla, self.drdz_right)

        bla2 = self.z_right(minimum + 1.0) - self.z_right(minimum)

        minr, minz_right, _ = self.get_min_dr_dzright_dzleft()
        self.dzdr_right = 1.0/self.drdz_right
        assert feq(self.dzdr_right, bla2), 'bla2= %s, dzdr_right= %s' % (bla2, self.dzdr_right)
        self.r0_right   = 0.0
        self.z0_right   = self.z_right(self.r0_right)
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.pid_particle_pair1[1].radius

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure.id, [self.single1.domain_id, self.single2.domain_id],
                                                          [self.surface.id])
        except ShellmakingError as e:
            raise testShellError('(MixedPair2D3D). %s' %
                                 (str(e)))


    def get_orientation_vector(self):
        # if the dot product is negative, then take the -unit_z, otherwise take unit_z
        return self.surface.shape.unit_z * cmp(numpy.dot(self.iv, self.surface.shape.unit_z), 0)

    def get_searchpoint(self):
        # TODO this can be improved, now same as PlanarSurfaceInteraction
        search_distance = self.get_searchradius()/math.sqrt(2) - self.z0_left
        return self.get_referencepoint() + search_distance * self.get_orientation_vector()

    def get_referencepoint(self):
        return self.com

    def get_min_dr_dzright_dzleft(self):
        # defines the minimum parameters of the testShell.
        # Note that also the minimum testShell must obey the scaling rules

        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D


        ### calculate the minimal height z_right1 of the shell including burst radius
        iv_z = self.r0      # approximation, actually the height of the domain can be slightly lower.
        dz_right1 = iv_z + radius2 * SINGLE_SHELL_FACTOR
        # with the accompanying radius r1
        dr_1 = self.r_right(dz_right1)


#        ### calculate the minimal radius r2 of the shell including the burst radius
#        # First calculate the minimal size of the iv shell (r0 is at outer boundary)
#        iv_shell_radius1 = self.r0 * D_1 / self.D_tot
#        iv_shell_radius2 = self.r0 * D_2 / self.D_tot
#
#        # fix the minimum shell size for the CoM domain
#        com_shell_radius = max(radius1, radius2)
#
#        # calculate the minimal dimensions of the protective domain including space for the
#        # burst volumes of the particles
#        dr_2 = max(iv_shell_radius1 + com_shell_radius + radius1 * SINGLE_SHELL_FACTOR,
#                   iv_shell_radius2 + com_shell_radius + radius2 * SINGLE_SHELL_FACTOR)
#        dz_right2 = self.z_right(dr_2)
#
#        # of both alternatives pick the largest one. Here we just compare the height, but the
#        # radius scales accordingly.
#        dz_right, dr = max((dz_right1, dr_1),(dz_right2, dr_2))

        dr = dr_1
        dz_right = dz_right1

        dz_left  = self.pid_particle_pair1[1].radius
#        print "min= ", dr, dz_right, dz_left
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        # TODO this can be improved
        dr       = self.get_searchradius()/math.sqrt(2)
        dz_left  = self.pid_particle_pair1[1].radius

        # now adjust dz_right such that they obey the scaling laws and are still with the search_radius
        # We assume that dz_right needs to be adjusted to dr and not the other way around since the cylinder is
        # always a bit flatter.
        dz_right = self.z_right(dr)
#        print "max= ", dr, dz_right, dz_left
        assert dr > dz_right

        return dr, dz_right, dz_left

    def z_right(self, r_right):
        # if the radius is known and we want to determine the height z_right
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        # calculate a_r such that the expected first-passage for the CoM and IV are equal
        a_r1 = (r_right - radius1 + self.r0*self.sqrt_DRDr ) / (self.sqrt_DRDr + (D_1/self.D_tot) )
        a_r2 = (r_right - radius2 + self.r0*self.sqrt_DRDr ) / (self.sqrt_DRDr + (D_2/self.D_tot) )
        # take the smallest that, if entered in the function for r below, would lead to this z_right
        a_r = min (a_r1, a_r2)

        z_right = (a_r/self.get_scaling_factor()) + radius2

#        assert feq(self.r_right(z_right), r_right), 'inconsistent z_right function \
#                                                     r_right= %s, r_right(z_right)=%s ' % \
#                                                     (r_right, self.r_right(z_right))

        return z_right

    def r_right(self, z_right):
        # if the z_right is known and we want to know the radius r
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        # we first calculate the a_r, since it is the only radius that depends on z_right only.
        a_r = (z_right - radius2) * self.get_scaling_factor()

        # We equalize the estimated first passage time for the CoM (2D) and IV (3D) for a given a_r
        # Note that for the IV we only take the distance to the outer boundary into account.
        a_R = (a_r - self.r0)*self.sqrt_DRDr

        # We calculate the maximum space needed for particles A and B based on maximum IV displacement
        iv_max = max( (D_1/self.D_tot * a_r + radius1),
                      (D_2/self.D_tot * a_r + radius2))

        r_right = a_R + iv_max

        assert feq(self.z_right(r_right), z_right), 'inconsistent r_right function \
                                                  z_right= %s, z_right(r_right)=%s ' % \
                                                  (z_right, self.z_right(r_right))

        return r_right

    def apply_safety(self, r, z_right, z_left):
        SAFETY = 1.1
        return r/SAFETY, self.z_right(r/SAFETY), z_left

#####
#class MixedPair3D1DtestShell(CylindricaltestShell, Others):
