#!/usr/bin/python

import math

from _gfrd import (
    Sphere,
    SphericalShell,
    SphericalSurface,
    Cylinder,
    CylindricalShell,    
    CylindricalSurface,
    Disk,
    DiskSurface,
    Plane,
    PlanarSurface,
    CuboidalRegion,
    )

from utils import *
from gfrdbase import (
    get_neighbor_surfaces,
    )

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
    'testPlanarSurfaceTransitionPair',
    'testMixedPair2D3D',
    'testMixedPair1DCap',
    'hasSphericalShell',
    'hasCylindricalShell',   
    'SphericalSingletestShell',
    'SphericalPairtestShell',
    'PlanarSurfaceSingletestShell',
    'PlanarSurfacePairtestShell',
    'PlanarSurfaceTransitionSingletestShell',
    'PlanarSurfaceTransitionPairtestShell',
    'CylindricalSurfaceSingletestShell',
    'DiskSurfaceSingletestShell',
    'CylindricalSurfacePairtestShell',
    'PlanarSurfaceInteractiontestShell',
    'CylindricalSurfaceInteractiontestShell',
    'CylindricalSurfaceCapInteractiontestShell',
    'CylindricalSurfaceSinktestShell',
    'MixedPair2D3DtestShell',
    'MixedPair1DCaptestShell',
    ]

import logging
log = logging.getLogger('ecell')


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
    # NOTE the two methods below have quite a different purpose!
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
        # This means that more bursting will take place than necessary because cylindrical domains do not keep the 
        # domains that are out of the 'diffusion dimension' at the proper distance -> TODO

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

    def __init__(self, single, origin_structure, target_structure):

        # assert that the single is reset
        assert single.is_reset()

        testSingle.__init__(self, single.pid_particle_pair, origin_structure)
        # Note: for the Others superclass nothing is to be initialized.

        self.single           = single
        self.origin_structure = origin_structure
        self.target_structure = target_structure

        ### initialize some constants that are common for the testInteractionSingle domains.

        # calculate the Projected_point and distance to the target_structure
        # Cyclic transpose needed when calling target_structure.project_point!
        pos_transposed = self.world.cyclic_transpose(self.pid_particle_pair[1].position,
                                                     self.target_structure.shape.position)
        self.reference_point, (_, dist_to_surface_edge) = self.target_structure.project_point(pos_transposed)
        if dist_to_surface_edge >= 0:
            raise testShellError('(testInteractionSingle). Projected point of particle is not in surface.')

        self.min_dist_proj_to_edge = abs(dist_to_surface_edge)
        # projection_distance can be negative in case of plane
        # note that the distance is to the center of the cylinder in case of cylinders
        self.particle_surface_distance = self.world.distance(self.target_structure.shape, self.pid_particle_pair[1].position)

        # The reference_vector is the normalized vector from the reference_point to the particle
        if any(pos_transposed != self.reference_point) :
            diff_vector = pos_transposed - self.reference_point

            # Compare to the largest component of the difference vector
            # If a component is TOLERANCE times smaller, regard it as zero
            # This is to avoid calculation errors when the distance to the
            # surface is very small / has very small nonzero components
            max_component = max(diff_vector, key = lambda x : abs(x))
            for c in range(0,3):
                if abs(diff_vector[c]) <= 1.0e-6 * abs(max_component):  # FIXME using TOLERANCE does not always work here                    
                    log.warning('Setting component %s of difference vector %s to zero.' % (c, diff_vector))
                    diff_vector[c] = 0.0

            self.reference_vector = normalize(diff_vector)
        else:
            self.reference_vector = 0
            raise testShellError('(testInteractionSingle). reference vector = 0.')

class testTransitionSingle(testSingle, Others):

    def __init__(self, single, target_structure):

        testSingle.__init__(self, single.pid_particle_pair, single.structure)
        assert (self.structure.sid == target_structure.sid)     # Both structures are of same structure_type

        self.single = single

        self.origin_structure = self.structure
        self.target_structure = target_structure
        self.structure1       = self.origin_structure
        self.structure2       = self.target_structure

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

        # In general we allow the two particles to be on different structures
        self.structure1 = single1.structure
        self.structure2 = single2.structure

        self.com, self.iv = self.do_transform()   # NOTE The IV always points from particle1 to particle2
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
        testPair.__init__(self, single1, single2)
        self.structure = single1.structure # equal to self.structure1 and self.structure2,
                                           # needed by some class methods downstream

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

        # First define for clarity:
        self.single2D = single2D
        self.single3D = single3D
        self.ppp2D      = single2D.pid_particle_pair
        self.ppp3D      = single3D.pid_particle_pair
        self.particle2D = single2D.pid_particle_pair[1]
        self.particle3D = single3D.pid_particle_pair[1]

        self.structure2D = single2D.structure   # equal to self.structure1
        self.structure3D = single3D.structure   # equal to self.structure2

        # These methods need above definitions.
        if __debug__: 
                assert isinstance(single2D.structure, PlanarSurface) 
                assert isinstance(single3D.structure, CuboidalRegion)
        testPair.__init__(self, single2D, single3D)
          # note: this makes self.single1/self.pid_particle_pair1/self.structure1 for the 2D particle
          #              and self.single2/self.pid_particle_pair2/self.structure2 for the 3D particle              

        assert self.single2D == self.single1           and self.single3D == self.single2
        assert self.ppp2D == self.pid_particle_pair1   and self.ppp3D == self.pid_particle_pair2
        assert self.structure2D == self.structure1     and self.structure3D == self.structure2

        # The following will be needed by some domain overlap checking functions
        # TODO: modification by Martijn Wehrens (jintram@gmail.com)
        #       -- not absolutely sure this is correct place
        self.origin_structure = self.structure3D
        self.target_structure = self.structure2D        

    def get_sigma(self):
        # rescale sigma to correct for the rescaling of the coordinate system.
        # This is the sigma that is used for the evaluation of the Green's function and is in this case 
        # slightly different from the sums of the radii of the particles. We rescale sigma in a way
        # that the total flux over the reactive surface of the prolate spheroidal boundary with greater radius
        # sigma is equal to the total flux over a spherical reactive surface with radius = rescaled sigma.
        xi = self.get_scaling_factor()
        xi_inv = 1.0/xi
        alpha = math.acos(xi_inv)

        sigma = self.particle2D.radius + self.particle3D.radius
        rescaled_sigma = abs(sigma * math.sqrt(0.5 + (alpha * xi/(2.0*math.sin(alpha)))))

        return rescaled_sigma
    sigma = property(get_sigma)

    def do_transform(self):
        # transform the pos1 and pos2 of particles1 and 2 to the CoM and IV vectors
        pos2D = self.particle2D.position
        pos3D = self.particle3D.position
        D_2D  = self.particle2D.D
        D_3D  = self.particle3D.D

        # CoM calculation
        # the CoM is calculated in a similar way to a normal 3D pair
        com = self.world.calculate_pair_CoM(pos2D, pos3D, D_2D, D_3D)

        # and then projected onto the plane to make sure the CoM is in the surface
        com = self.world.cyclic_transpose(com, self.structure2D.shape.position)
        com, (_, dist) = self.structure2D.project_point(com)
        if dist >= 0:
            raise testShellError('(testMixedPair2D3D). Projected point of CoM is not in surface.')
        else:
            self.min_dist_proj_to_edge = -dist      # FIXME this is a bit of a hack, find better place to put this.

        com = self.world.apply_boundary(com)

        # IV calculation
        pos3D_t = self.world.cyclic_transpose(pos3D, pos2D)
        iv      = pos3D_t - pos2D

        # rescale the iv in the axis normal to the plane
        iv_z_component = numpy.dot (iv, self.structure2D.shape.unit_z)
        iv_z     = self.structure2D.shape.unit_z * iv_z_component
        new_iv_z = self.structure2D.shape.unit_z * (iv_z_component * self.get_scaling_factor())

        iv = iv - iv_z + new_iv_z

        return com, iv

    def get_scaling_factor(self):
        # This returns the scaling factor needed to make the anisotropic diffusion problem of the IV
        # into a isotropic one.
        D_3D = self.particle3D.D    # 3D particle is the only contributor to diffusion
                                    # normal to the plane
        return math.sqrt(self.D_r/D_3D)

class testMixedPair1DCap(testPair):

    def __init__(self, single1D, cap_single):

        if __debug__: 
                assert isinstance(single1D.structure, CylindricalSurface) 
                assert isinstance(cap_single.structure, DiskSurface)
        testPair.__init__(self, single1D, cap_single)
          # note: this makes self.single1/self.pid_particle_pair1/self.structure1 for the 1D particle
          #              and self.single2/self.pid_particle_pair2/self.structure2 for the cap particle

        # for clarity define:
        self.single1D       = single1D
        self.particle1D     = self.pid_particle_pair1[1]
        self.structure1D    = single1D.structure    # equal to self.structure1  
        self.cap_single     = cap_single
        self.cap_particle   = self.pid_particle_pair2[1]        
        self.cap_structure  = cap_single.structure  # equal to self.structure2
        # This will be inherited by MixedPair1DCaptestShell and MixedPair1DCap.
        # Note that we have to use pid_particle_pair1(2) in the methods defined below
        # because they are called at the initialization of testPair().

        # Make sure cap particle is immobile; this is required for testSimplePair routines
        # to give the proper results for this case.
        assert self.cap_particle.D == 0
        assert self.cap_particle.v == 0

    def get_sigma(self):
        # Copied from SimplePair
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        return radius1 + radius2
    sigma = property(get_sigma)

    def do_transform(self):

        pos1 = self.pid_particle_pair1[1].position
        pos2 = self.pid_particle_pair2[1].position
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        # In the MixedPair1DCap only the 1D particle moves. We therefore use
        # iv to sample its position and com for the position of the static
        # cap particle; this is initialized here:
        com = pos2
        com = self.world.apply_boundary(pos2)
        
        pos2t = self.world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        return com, iv

class testPlanarSurfaceTransitionPair(testPair):

    # This is essentially a copy of testSimplePair; however, the latter does
    # not allow for particles being on different structures. Ergo we have to
    # define a new class for that purpose.
    # The main difference is in do_transform and get_min_pair_size which use
    # the CoM calculated after projection of single2 into the plane of single1
    # via deflect_back.

    def __init__(self, single1, single2):

        if __debug__: 
                assert isinstance(single1.structure, PlanarSurface) 
                assert isinstance(single2.structure, PlanarSurface) 

        testPair.__init__(self, single1, single2) # note: this makes self.single1/self.pid_particle_pair1 and the same for particle2
        self.structure  = self.structure1         # some pair methods need this to be defined

    def get_sigma(self):
        # Nothing changed as compared to the simple pair here
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        return radius1 + radius2
    sigma = property(get_sigma)

    def do_transform(self):
        # Transform the pos1 and pos2 of particles1 and 2 to the CoM and IV vectors
        # As an indermediate step project position of single2 into plane of single1 and add its
        # orthogonal component rotated by 90 deg. to transform the whole problem into a single 2D plane
        pos1 = self.pid_particle_pair1[1].position
#        pos2 = self.structure1.deflect_back(self.pid_particle_pair2[1].position, self.structure2.shape.unit_z)
        pos2, _ = self.world.cyclic_transpose((self.pid_particle_pair2[1].position, self.pid_particle_pair2[1].structure_id), self.structure1)
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        com = self.world.calculate_pair_CoM(pos1, pos2, D_1, D_2)
        com = self.world.apply_boundary(com)
        # TODO make sure that the com is in the structure of the particle
        # (assume that this is done correctly in calculate_pair_CoM)

        pos2t = self.world.cyclic_transpose(pos2, pos1)
        iv = pos2t - pos1

        return com, iv

    def get_min_pair_size(self):
        # This calculates the minimal 'size' of a simple pair.
        # Note that the 'size' is interpreted differently dependent on the type of pair. When the pair is a
        # SphericalPair          -> size = radius
        # PlanarSurfacePair      -> size = radius
        # CylindricalSurfacePair -> size = half_length

        # TODO: Make sure this also works correctly when IP-vectors are enlarged to remove
        # overlaps at deflection.

        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        dist_from_com1 = self.r0 * D_1 / self.D_tot     # particle distance from CoM
        dist_from_com2 = self.r0 * D_2 / self.D_tot
        iv_shell_size1 = dist_from_com1 + radius1       # the shell should surround the particles
        iv_shell_size2 = dist_from_com2 + radius2

        # also fix the minimum shellsize for the CoM domain
        com_shell_size = max(radius1, radius2)

        # calculate total radii including the margin for the burst volume for the particles
        shell_size = max(iv_shell_size1 + com_shell_size + radius1 * SINGLE_SHELL_FACTOR,
                         iv_shell_size2 + com_shell_size + radius2 * SINGLE_SHELL_FACTOR)

        return shell_size


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
    # of the testShell can be. This means, how far can (the other) testShell be scaled according to its scaling rules
    # before it hits THIS Cylindrical shell.

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
        scale_angle     = testShell.right_scalingangle
        tan_scale_angle = testShell.tan_right_scalingangle
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
        tan_scale_angle = testShell.tan_left_scalingangle
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
            scale_center_to_shell_z_thres = abs(scale_center_to_shell_r)/tan_scale_angle

        shell_radius_thres = scale_center_to_shell_z - scale_center_to_shell_z_thres
        if shell_radius < shell_radius_thres:
            situation = 1
        else:
            situation = 2

    elif scale_angle <= shell_angle and shell_angle < Pi/2.0:
    # The (spherical) shell can touch the cylinder on its edge or its radial side
        # scale_angle == Pi/2 should not get here
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
    # using a cylindrical 'shell' as its closest neighbor. Note that it returns the minimum of 
    # the newly calculated dr, dz_right, dz_left and the old r, z_right, z_left. The 'testShell'
    # therefore can only become smaller.

    # Note that the neighboring 'shell' is querried from this domain earlier by the testShell.
    # -> the shell MUST be a Cylinder.
    assert(type(shape) is Cylinder)

    # Laurens' algorithm (part2)
    shell_position = shape.position
    shell_radius = shape.radius 
    shell_half_length = shape.half_length 

    # Get the reference point and orientation of the domain to scale
    reference_point = testShell.get_referencepoint()
    orientation_vector = testShell.get_orientation_vector()

    # FIXME DIRTY HACK TO MAKE DOMAINS LESS GREEDY
#    r = min(r, testShell.world.distance(shape.position, reference_point)/2)

    # Determine on what side the midpoint of the shell is relative to the reference_point
    shell_position_t = testShell.world.cyclic_transpose (shell_position, reference_point)
    ref_to_shell_vec = shell_position_t - reference_point
    ref_to_shell_z = numpy.dot(ref_to_shell_vec, orientation_vector)
    # The length ref_to_shell_z is positive when the shell is on the right side of the
    # cylinder and negative when it is on the left side

    # If the shell is on the side of orientation_vector -> use z_right
    # Also use this one if the shell is in plane with the reference_point
    if ref_to_shell_z >= 0:
        direction = 1           # improve direction specification such that following calculations still work
                                # Otherwise problems when direction = 0
        scale_angle     = testShell.right_scalingangle
        tan_scale_angle = testShell.tan_right_scalingangle
        scale_center_r, scale_center_z = testShell.right_scalingcenter 
        r1_function = testShell.r_right
        z1_function = testShell.z_right
        z2_function = testShell.z_left
        z1 = z_right
        z2 = z_left
    # else -> use z_left
    else:
        direction = -1
        scale_angle     = testShell.left_scalingangle
        tan_scale_angle = testShell.tan_left_scalingangle
        scale_center_r, scale_center_z = testShell.left_scalingcenter 
        r1_function = testShell.r_left
        z1_function = testShell.z_left
        z2_function = testShell.z_right
        z1 = z_left
        z2 = z_right
    
    # This must always be the case:
    assert scale_angle >= 0.0

    # Check how the cylinders are oriented with respect to each other
    relative_orientation = abs(numpy.dot(orientation_vector, shape.unit_z))

    if feq(relative_orientation, 1.0):
    #### If the cylinders are oriented parallelly ####

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
    #### If the cylinders are oriented perpendicularly ####

        # First construct a local coordinate system:
        #   - The local x-axis is defined as the vector parallel to the axis of the neighboring shell pointing
        #   from the reference_point of the scaling cylinder towards the neighboring shell.
        #   - The local y-axis is defined as the vector parallel to the flat side of the shell, perpendicular to the
        #   axis of the scaling cylinder. It also points from the reference_point of the scaling cylinder to the shell.
        #   - The local z-axis is defined as the vector parallel to the axis of the scaling cylinder pointing
        #   from the reference point towards the neighboring shell.
        if numpy.dot(ref_to_shell_vec, shape.unit_z) >= 0:
            local_x = shape.unit_z
        else:
            local_x = -shape.unit_z
        local_z = orientation_vector * direction
        cross = numpy.cross(local_x, local_z)
        if numpy.dot(ref_to_shell_vec, cross) >= 0:
            # cross already points towards the shell, we are fine
            local_y = cross
        else:
            # it points away from the shell, thus negate it
            local_y = -cross

        # Transform the vector pointing from the ref. point towards the shell
        # into the local coordinate system
        ref_to_shell_x = numpy.dot(ref_to_shell_vec, local_x)
        ref_to_shell_y = numpy.dot(ref_to_shell_vec, local_y)
        ref_to_shell_z = numpy.dot(ref_to_shell_vec, local_z)
        assert (ref_to_shell_x >= 0.0) and (ref_to_shell_y >= 0.0) and (ref_to_shell_z >= 0.0)  # by the definition of the local coordinate system
        # Also calculate the coordinates of the scale center
        # FIXME Attention: The scale center is not always on the z-axis!
        scale_center_to_shell_x = ref_to_shell_x
        scale_center_to_shell_y = ref_to_shell_y
        scale_center_to_shell_z = ref_to_shell_z - scale_center_z
        # And the coordinates of the closest corner of a box surrounding the shape
        ref_to_shell_x2 = ref_to_shell_x - shell_half_length
        ref_to_shell_y2 = ref_to_shell_y - shell_radius
        ref_to_shell_z2 = ref_to_shell_z - shell_radius # not really needed, but defined here for completeness

        # Define situation (collision) specifiers
        # Depending on the geometric constellation the scaling will eventually lead to different collision behaviour
        # which has to be treated in different ways; the following collision specifiers are from the point of view
        # of the scaled cylinder:
        BARREL_HITS_FLAT, EDGE_HITS_EDGE, BARREL_HITS_EDGE, FLAT_HITS_BARREL, EDGE_HITS_BARREL, BARREL_HITS_BARREL = range(6)
        situation_string = ["BARREL_HITS_FLAT", "EDGE_HITS_EDGE", "BARREL_HITS_EDGE", "FLAT_HITS_BARREL", "EDGE_HITS_BARREL", "BARREL_HITS_BARREL"]

        #### (1) Now first determine which type of collision happens:

        if (ref_to_shell_x2 < 0) and (ref_to_shell_y2 < 0):
            # Quadrant 1
            quadrant = 1

            # Scale cylinder one until its height is equal to the z-distance minus the radius of the cyl. shell 2;
            # Then check whether any point of the axis of cyl. 2 projected on the top side of cyl. 1 is inside the radius of (scaled) cyl. 1.
            # If this is the case, we have the flat side of cyl. 1 hitting the barrel of cyl. 2. Else we have the edge hitting the barrel.
            r1_touch = (scale_center_to_shell_z - shell_radius)*tan_scale_angle

            if r1_touch >= scale_center_to_shell_y:
                # At least a part of the projected axis is within r_touch => flat side of scaled cyl. hits the barrel of cyl. 2
                # Note that we know that at least a part of the shell is above the top side of the scaled cylinder (quadrant 1 condition)
                situation = FLAT_HITS_BARREL
            else:
                situation = EDGE_HITS_BARREL

        elif (ref_to_shell_x2 >= 0) and (ref_to_shell_y2 < 0):
            # Quadrant 2
            quadrant = 2

            scale_center_to_shell_x_minus_half_length = scale_center_to_shell_x - shell_half_length

            # The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            # because in this case the mathematics in the standard case misdetect the collision situation
            if scale_angle == 0:
                
                if math.sqrt( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length + \
                              scale_center_to_shell_y*scale_center_to_shell_y) < r :
                # The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder.
                # Therefore, when the height is scaled, the flat side must hit the barrel of the static cylinder at the lowpoint.
                    situation = FLAT_HITS_BARREL

                else:
                    situation = EDGE_HITS_EDGE
                    # This will also properly treat the case in which there is no collision because the projection of
                    # the static cylinder does not overlap with the flat side circle of the scaled cylinder
            else:
                # case scale_angle >0

                # Two points on the edge of the static cylinder are of interest here:
                # - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
                # - the "critpoint" (critical point), which is the point on the edge with the shortest distance to the scale center                  
                scale_center_to_flatend_x   = scale_center_to_shell_x_minus_half_length - scale_center_r
                scale_center_to_lowpoint_x  = math.sqrt( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length +\
                                                         scale_center_to_shell_y*scale_center_to_shell_y ) - scale_center_r
                    # TODO Is there no better way to include scale_center_r here?
                scale_center_to_critpoint_z = scale_center_to_shell_z - math.sqrt(shell_radius*shell_radius - scale_center_to_shell_y*scale_center_to_shell_y) 
                scale_center_to_lowpoint_z  = scale_center_to_shell_z - shell_radius

                # Strategy:
                # - When the scale angle is smaller than the angle between the scaled cylinder's axis and the line
                #   that links the critpoint then we have a BARREL_HITS_FLAT situation.
                # - When the scale angle is bigger than the angle between the scaled cylinder's axis and the line
                #   that links the lowpoint then we are in a FLAT_HITS_BARREL situation.                  
                # - Everything else results in EDGE_HITS_EDGE.

                # Calculate the angle between the scaled cylinder's axis and the line that links the scale center and lowpoint
                if scale_center_to_critpoint_z == 0:
                    scale_center_to_shell_crit_angle_y = Pi/2
                else:
                    scale_center_to_shell_crit_angle_y = math.atan(scale_center_to_flatend_x/ scale_center_to_critpoint_z)
                    if scale_center_to_critpoint_z < 0.0:
                        scale_center_to_shell_crit_angle_y += Pi

                # Calculate the angle between the scaled cylinder's axis and the line that links the scale center and critpoint
                if scale_center_to_lowpoint_z == 0:
                    scale_center_to_shell_low_angle_y = Pi/2
                else:
                    scale_center_to_shell_low_angle_y  = math.atan(scale_center_to_lowpoint_x/ scale_center_to_lowpoint_z )
                    if scale_center_to_lowpoint_z < 0.0:
                        scale_center_to_shell_low_angle_y += Pi

                # Compare the angles to determine the situation
                if scale_angle <= scale_center_to_shell_crit_angle_y:
                    situation = BARREL_HITS_FLAT
                elif scale_angle >= scale_center_to_shell_low_angle_y:
                    situation = FLAT_HITS_BARREL
                else:
                    assert scale_center_to_shell_crit_angle_y < scale_angle and \
                            scale_angle < scale_center_to_shell_low_angle_y
                    situation = EDGE_HITS_EDGE
                    r1_min = (scale_center_to_shell_x - shell_half_length)*(1.0+TOLERANCE)
                    h1_min = r1_min/tan_scale_angle

        elif (ref_to_shell_x2 < 0) and (ref_to_shell_y2 >= 0):
            # Quadrant 3
            quadrant = 3

            # The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            # because in this case the mathematics in the standard case misdetect the collision situation
            if scale_angle == 0:

                if scale_center_to_shell_y < r :
                # The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder.
                # Therefore, when the height is scaled, the flat side must hit the barrel of the static cylinder at the lowpoint.
                    situation = FLAT_HITS_BARREL

                elif scale_center_to_shell_y - shell_radius <= r:
                    situation = EDGE_HITS_BARREL

                else:
                    # In this case the cylinders do not hit; this is treated properly by the BARREL_HITS_BARREL routine
                    situation = BARREL_HITS_BARREL

            else:

                # Two points on the edge of the static cylinder are of interest here:
                # - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
                # - the "critpoint" (critical point), which is the point on the edge with the shortest distance to the scale center   
                scale_center_to_shell_y -= scale_center_r  # TODO Is there no better way to include scale_center_r here?
                scale_center_to_critpoint_y = scale_center_to_shell_y - shell_radius
                scale_center_to_lowpoint_z  = scale_center_to_shell_z - shell_radius

                # Calculate the angle between the scaled cylinder's axis and the line that links the scale center and critpoint
                if scale_center_to_shell_z == 0.0:
                    scale_center_to_shell_crit_angle_x = Pi/2.0
                else:
                    scale_center_to_shell_crit_angle_x = math.atan(scale_center_to_critpoint_y / scale_center_to_shell_z)
                    if scale_center_to_shell_z < 0.0:
                        scale_center_to_shell_crit_angle_x += Pi
                # Calculate the angle between the scaled cylinder's axis and the line that links the scale center and lowpoint
                if scale_center_to_lowpoint_z == 0.0:
                    scale_center_to_shell_low_angle_x = Pi/2.0
                else:
                    scale_center_to_shell_low_angle_x  = math.atan(scale_center_to_shell_y/scale_center_to_lowpoint_z)
                    if scale_center_to_lowpoint_z < 0.0:
                        scale_center_to_shell_low_angle_x += Pi

                # Determine which collision will happen:
                # If the scale angle is bigger than the angle between the scaled cylinder's axis and the lowpoint
                # then we hit the barrel of the static cylinder when scaling.
                # If the scale angle is smaller than the angle between the scaled cylinder's axis and the critpoint
                # the barrel of the scaled cylinder will hit the barrel of the static cylinder.
                # Otherwise we have the edge of the scaled cylinder hitting the barrel of the static cylinder.
                if scale_angle <= scale_center_to_shell_crit_angle_x:
                    situation = BARREL_HITS_BARREL
                elif scale_center_to_shell_low_angle_x <= scale_angle:
                    situation = FLAT_HITS_BARREL
                else:
                    assert scale_center_to_shell_crit_angle_x < scale_angle and \
                          scale_angle < scale_center_to_shell_low_angle_x
                    situation = EDGE_HITS_BARREL

        else:
            # Quadrant 4
            quadrant = 4

            assert (ref_to_shell_x2 >= 0) and (ref_to_shell_y2 >= 0)

            scale_center_to_shell_x_minus_half_length = scale_center_to_shell_x - shell_half_length
            scale_center_to_shell_y_minus_radius      = scale_center_to_shell_y - shell_radius

            # Two points on the edge of the static cylinder are of interest here:
            # - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
            # - the "critpoint" (critical point), which is the point on the edge with the minimal y-distance to
            #   the scale center in the zy-plane
            scale_center_to_critpoint_r = math.sqrt(scale_center_to_shell_y_minus_radius*scale_center_to_shell_y_minus_radius + \
                                                    scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length) - scale_center_r  # a_r
            scale_center_to_lowpoint_r = math.sqrt(scale_center_to_shell_y*scale_center_to_shell_y + \
                                                   scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length) - scale_center_r   # b_r
            scale_center_to_lowpoint_z = scale_center_to_shell_z - shell_radius                                                                            # b_z

            # Calculate the angle between the z-axis of the scaled cylinder and the line that
            # connects this axis with the "critpoint" on the static cylinder's edge
            if scale_center_to_shell_z == 0.0:
                scale_center_to_shell_crit_angle_xy = Pi/2.0
            else:
                scale_center_to_shell_crit_angle_xy = math.atan(scale_center_to_critpoint_r / scale_center_to_shell_z)
                if scale_center_to_shell_z < 0.0:
                    scale_center_to_shell_crit_angle_xy += Pi

            # Calculate the angle between the z-axis of the scaled cylinder and the line that
            # connects this axis with the "lowpoint" on the static cylinder's edge
            if scale_center_to_lowpoint_z == 0.0:
                scale_center_to_shell_low_angle_xy = Pi/2.0
            else:
                scale_center_to_shell_low_angle_xy  = math.atan(scale_center_to_lowpoint_r / scale_center_to_lowpoint_z )
                if scale_center_to_lowpoint_z < 0.0:
                    scale_center_to_shell_low_angle_xy += Pi

            # First treat the special case: when the scale center is further away from the z-axis than the "critpoint" (in the xy-plane).
            # This also treats the case in which scale_angle = 0 and the scaled cylinder is directly below the static shell.
            if scale_center_to_critpoint_r <= 0.0:
                situation = EDGE_HITS_EDGE

                if scale_center_to_lowpoint_r <= 0.0:
                    situation = FLAT_HITS_BARREL

            # Now the standard case:
            else:

                if scale_angle <= scale_center_to_shell_crit_angle_xy:
                    situation = BARREL_HITS_EDGE

                elif scale_center_to_shell_low_angle_xy <= scale_angle:
                    situation = FLAT_HITS_BARREL

                else:
                    assert scale_center_to_shell_crit_angle_xy < scale_angle and scale_angle < scale_center_to_shell_low_angle_xy and scale_angle > 0.0
                    situation = EDGE_HITS_EDGE
                    r1_min = math.sqrt(scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length + \
                                       scale_center_to_shell_y_minus_radius*scale_center_to_shell_y_minus_radius) * (1.0+TOLERANCE)
                    h1_min = r1_min/tan_scale_angle

        ##TODO TESTING Debug info
        #log.info("  *** QUADRANT = %s ***" % str(quadrant) )
        #log.info("  *** Situation = %s" % str(situation_string[situation]) )
        #log.info("  *** testShell=%s, r=%s, z1=%s" % (str(testShell), r, z1) )
        #log.info("  *** scale_angle=%s, tan_scale_angle=%s, scale_center_z=%s, scale_center_r=%s" % (scale_angle, tan_scale_angle, scale_center_z, scale_center_r) )
        #log.info("  *** shell_radius=%s, shell_half_length=%s" % (shell_radius, shell_half_length) )
        #log.info("  *** ref_to_shell_x=%s, ref_to_shell_y=%s, ref_to_shell_z=%s" % (ref_to_shell_x, ref_to_shell_y, ref_to_shell_z) )
        #log.info("  *** scale_center_to_shell_x=%s, scale_center_to_shell_y=%s, scale_center_to_shell_z=%s" % (scale_center_to_shell_x, scale_center_to_shell_y, scale_center_to_shell_z) )

        ##### (2) Treat the situation accordingly        
        if situation == BARREL_HITS_FLAT:
            # the scaling cylinder hits the flat side of 'shell' with its barrel side
            r_new = min(r, (ref_to_shell_x - shell_half_length))
            z1_new = min(z1, z1_function(r_new))

        elif situation == EDGE_HITS_EDGE:             
            # the scaling cylinder hits the edge of 'shell' with its edge
            # TODO we have a solution but it can only be found with a root finder -> slow
            shell_radius_sq = shell_radius*shell_radius
            scale_center_to_shell_edge_x = scale_center_to_shell_x - shell_half_length
            scale_center_to_shell_edge_y = scale_center_to_shell_y - shell_radius
            scale_center_to_shell_edge_z = scale_center_to_shell_z - shell_radius

            if scale_angle == 0.0:

                # scale_angle = 0 means that the radius r will stay constant at scaling
                # Therefore, if the projected edge of the static shell lies outside of the 
                # circle defined by r there will be no collision:
                if r < math.sqrt( scale_center_to_shell_edge_x*scale_center_to_shell_edge_x + \
                                  scale_center_to_shell_edge_y*scale_center_to_shell_edge_y ) :

                    # The height and radius of the scaled cylinder remain unaltered
                    r_new = r
                    z1_new = z1         # Note that z1_function would give infinity,
                                        # but we would take the min(z1_function(r_new), z1)
                else:
                # There is a collision possible when z1 is scaled
                    
                    y1_collide = math.sqrt(r*r - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x)
                    y2_collide = scale_center_to_shell_y - y1_collide
                    h_collide = math.sqrt(shell_radius*shell_radius - y2_collide*y2_collide)

                    h_touch = scale_center_to_shell_z - h_collide
                    z1_new = min(z1, h_touch)
                    r_new  = min(r,  r1_function(z1_new))

            else:

                # TESTING Code snippet to check whether rootfinder run is necessary
                # We only want to run it if the radius increase at intersection of the cylinder edges is "significant".
                # This is sub-optimal but it avoids rootfinder errors caused by the fact that the radius search interval
                # becomes to small, which is equivalent to the potential r-increase becoming small.
                #
                # We define the increase as "significant" by comparing the potential radius increase with the radial (xy-) distance
                # from the scale center to the static shell; if increase / distance >= 10% the gain is "significant".
                #radial_distance = math.sqrt( scale_center_to_shell_x**2 + scale_center_to_shell_y**2 )
                #potential_r_increase = math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2 ) - scale_center_to_shell_edge_x  
                #assert radial_distance >= potential_r_increase and potential_r_increase >= 0

                #if( potential_r_increase / radial_distance >= 0.10 ):

                    #run_rootfinder = True

                #else:

                    #run_rootfinder = False
        
                    #if scale_center_to_shell_x >= scale_center_to_shell_y:
                        #r_touch = scale_center_to_shell_edge_x
                    #else:
                        #r_touch = scale_center_to_shell_edge_y
                    
                    #r_new  = min(r, r_touch)
                    #z1_new = min(z1, z1_function(r_new))


            #if run_rootfinder:

                if scale_angle <= Pi/4.0:

                    def h1_eq(x):

                        #rhs = (scale_center_to_shell_z - \
                              #math.sqrt(shell_radius**2 - (scale_center_to_shell_y - math.sqrt( (x*tan_scale_angle)**2 - (scale_center_to_shell_x - shell_half_length)**2 ) )**2 ) )

                        # TODO TESTING REMOVE THIS WHEN DONE
                        #log.debug( "***** ROOTFINDER CALL: Calling h1_eq() with: *****" )
                        #log.debug( "  r1 = %s" % (r1_function(x+scale_center_z)) )
                        #log.debug( "  sqrt_arg = %s" % ((x*tan_scale_angle)**2 - scale_center_to_shell_edge_x**2) )
                        #log.debug( "    tan_scale_angle = %s" % (tan_scale_angle) )
                        #log.debug( "    x*tan_scale_angle = %s" % (x*tan_scale_angle) )
                        #log.debug( "    scale_center_to_shell_edge_x = %s" % scale_center_to_shell_edge_x )

                        # We will take the square root of the following below
                        sqrt_arg = (x*tan_scale_angle)*(x*tan_scale_angle) - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x

                        if sqrt_arg < 0.0 and abs(sqrt_arg) <= TOLERANCE*scale_center_to_shell_edge_x*scale_center_to_shell_edge_x:

                            sqrt_arg = 0.0      # This safety check is to prevent math domain errors
                                                # in case sqrt_arg is close to zero and taking the
                                                # difference results in very small negative numbers
                            log.warn('Orthogonal cylinder scaling, EDGE_HITS_EDGE case: Setting negative sqrt argument to zero within TOLERANCE.')

                        eq_value = (scale_center_to_shell_z - x)*(scale_center_to_shell_z - x) - shell_radius_sq +\
                                      (scale_center_to_shell_y - math.sqrt( sqrt_arg ))*(scale_center_to_shell_y - math.sqrt( sqrt_arg ))
                        
                        # TODO TESTING REMOVE THIS WHEN DONE                    
                        #log.debug( "  h1 = %s, value = %s" % (x, eq_value ) )

                        return eq_value

                    # Self-adaptive initial value guessing for the rootfinder procedure:
                    # Start with prefactors (1+TOLERANCE), (1-TOLERANCE) and check whether the h1_eq values for these
                    # interval boundaries have different signs; if not, slightly enlarge the interval.
                    n=1
                    nmax=100
                    h1_eq_product = 1
                    assert TOLERANCE>0.0 and TOLERANCE<1.0
                    while h1_eq_product > 0.0:

                        h1_interval_start = (1.0+TOLERANCE**n) * (z1_function(scale_center_to_shell_edge_x) - scale_center_z)
                        h1_interval_end   = (1.0-TOLERANCE**n) * (z1_function( math.sqrt( scale_center_to_shell_y*scale_center_to_shell_y + \
                                                                                          scale_center_to_shell_edge_x*scale_center_to_shell_edge_x ) ) - scale_center_z)
                        h1_eq_product     = h1_eq(h1_interval_start) * h1_eq(h1_interval_end)

                        n = n+1

                        if n>nmax:
                            raise testShellError('get_dr_dzright_dzleft_to_CylindricalShape: Could not find suitable rootfinder boundaries in EDGE_HITS_EDGE cylinder scaling. nmax=%s' % str(nmax) )

                    assert(h1_interval_start >= 0)
                    assert(h1_interval_end   >= h1_interval_end)

                    # TODO TESTING REMOVE THIS WHEN DONE
                    #log.debug( "***** NEW ROOTFINDER ITERATION *****" )
                    #log.debug( "  Dy-r2 = %s" % (scale_center_to_shell_y-shell_radius) )
                    #log.debug( "  Dx-h2 = %s" % scale_center_to_shell_edge_x )
                    #log.debug( "  h1_interval_start = %s" % h1_interval_start )
                    #log.debug( "  h1_interval_end = %s"   % h1_interval_end )
                    #log.debug( "  h1(i_start) = %s"       % h1_eq(h1_interval_start) )
                    #log.debug( "  h1(i_end) = %s"         % h1_eq(h1_interval_end) )

                    h_touch = scale_center_z + findroot(h1_eq, h1_interval_start, h1_interval_end)
                    z1_new = min(z1, h_touch)
                    r_new  = min(r,  r1_function(z1_new))

                else:
                    # This uses the same rootfinding equation and interval as above with r1=h1/tan_scale_angle
                    # instead of h1; notice that scale_angle > 0 here, so we can divide by it
                    def r1_eq(x):

                        # We will take the square root of the following below
                        sqrt_arg = x*x - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x

                        if sqrt_arg < 0.0 and abs(sqrt_arg) <= TOLERANCE*scale_center_to_shell_edge_x*scale_center_to_shell_edge_x:

                            sqrt_arg = 0.0      # This safety check is to prevent math domain errors
                                                # in case sqrt_arg is close to zero and taking the
                                                # difference results in very small negative numbers
                            log.warn('Orthogonal cylinder scaling, EDGE_HITS_EDGE case: Setting negative sqrt argument to zero within TOLERANCE.')

                        eq_value = (scale_center_to_shell_z - x/tan_scale_angle)*(scale_center_to_shell_z - x/tan_scale_angle) - shell_radius_sq +\
                                      (scale_center_to_shell_y - math.sqrt( sqrt_arg ) )*(scale_center_to_shell_y - math.sqrt( sqrt_arg ) )

                        # TODO TESTING REMOVE THIS WHEN DONE
                        #log.debug( "***** ROOTFINDER CALL: Calling r1_eq() with: *****" )
                        #log.debug( "  r1 = %s, value = %s" % (x, eq_value ) )
                        #log.debug( "  h1 = %s" % (z1_function(x)) )

                        return eq_value

                    n=1
                    nmax=100
                    r1_eq_product = 1
                    assert TOLERANCE>0.0 and TOLERANCE<1.0
                    while r1_eq_product > 0.0:

                        r1_interval_start = (1.0+TOLERANCE**n) * scale_center_to_shell_edge_x
                        r1_interval_end   = (1.0-TOLERANCE**n) * math.sqrt( scale_center_to_shell_y*scale_center_to_shell_y + \
                                                                            scale_center_to_shell_edge_x*scale_center_to_shell_edge_x )
                        r1_eq_product     = r1_eq(r1_interval_start) * r1_eq(r1_interval_end)
                        
                        n = n+1

                        if n>nmax:
                            raise testShellError('get_dr_dzright_dzleft_to_CylindricalShape: Could not find suitable rootfinder boundaries in EDGE_HITS_EDGE cylinder scaling. nmax=%s' % str(nmax) )

                    # Non-adaptive version; TODO Remove this after TESTING
                    #r1_interval_start = scale_center_to_shell_edge_x
                    #r1_interval_end   = math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2 )
                    assert(r1_interval_start >= 0)
                    assert(r1_interval_end   >= r1_interval_start)

                    # TODO TESTING REMOVE THIS WHEN DONE
                    #print "***** NEW ROOTFINDER ITERATION *****"
                    #print "  Dy-r2 = %s" % (scale_center_to_shell_y-shell_radius)
                    #print "  Dx-h2 = %s" % scale_center_to_shell_edge_x
                    #print "  r1_interval_start = %s" % r1_interval_start
                    #print "  r1_interval_end = %s"   % r1_interval_end
                    #print "  r1(i_start) = %s"       % r1_eq(scale_center_to_shell_edge_x)
                    #print "  r1(i_end) = %s"         % r1_eq(math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2))

                    r_touch = findroot(r1_eq, r1_interval_start, r1_interval_end)
                              
                    r_new  = min(r, r_touch)
                    z1_new = min(z1, z1_function(r_new))


        elif situation == BARREL_HITS_EDGE:
            # the scaling cylinder hits the edge of 'shell' with its barrel side
            r_new = min(r, math.sqrt( (ref_to_shell_x-shell_half_length)*(ref_to_shell_x-shell_half_length) \
                                      + (ref_to_shell_y - shell_radius)*(ref_to_shell_y - shell_radius)))
            z1_new = min(z1, z1_function(r_new))

        elif situation == FLAT_HITS_BARREL:
            # the scaling cylinder hits the barrel of 'shell' with its top flat side
            z1_new = min(z1, (ref_to_shell_z - shell_radius))
            r_new  = min(r,  r1_function(z1_new))

        elif situation == EDGE_HITS_BARREL:
            # the scaling cylinder hits the barrel of 'shell' with its edge
            shell_radius_sq = shell_radius*shell_radius

            # If scale_angle == 0, the scale center is offset from the reference point,
            # located at the barrel of the scaled cylinder. In order for the below calculation
            # to succeed we have to take this into account. In general scale_center_r == 0
            # so that this step only affects the scale_angle == 0 cases:
            scale_center_to_shell_y -= scale_center_r

            ss_sq = (scale_center_to_shell_z*scale_center_to_shell_z + scale_center_to_shell_y*scale_center_to_shell_y)
            scale_center_to_shell = math.sqrt(ss_sq)

            shell_angle_yz = math.atan(scale_center_to_shell_y/scale_center_to_shell_z)
            angle_diff = abs(shell_angle_yz - scale_angle)
            sin_angle_diff = math.sin(angle_diff)
            cos_angle_diff = math.cos(angle_diff)
            
            assert scale_center_to_shell >= shell_radius # should never fail
            #ss_angle = math.asin(math.sin(angle_diff)*scale_center_to_shell/shell_radius) ## TESTING Old version
            ss_angle = math.pi - math.asin(math.sin(angle_diff)*scale_center_to_shell/shell_radius)
            assert(ss_angle >= math.pi/2.0), 'EDGE_HITS_BARREL: ss_angle must be >= Pi/2.0; ss_angle = %s' % str(ss_angle)
                    # We know that this angle must be larger than Pi/2 here by construction of the problem

            scale_center_shell_dist = shell_radius * math.sin(math.pi-(angle_diff+ss_angle)) / math.sin(angle_diff)
            assert scale_center_shell_dist>0, 'EDGE_HITS_BARREL: scale_center_shell_dist not positive, scale_center_shell_dist = %s'\
                                              % str(scale_center_shell_dist)

            # Check whether the calculated distance is within the forseen bounds and warn if not
            if scale_center_shell_dist > math.sqrt(scale_center_to_shell_y*scale_center_to_shell_y \
                                                   + scale_center_to_shell_z*scale_center_to_shell_z) * (1.0+TOLERANCE):
                log.warning('Orthogonal cylinder scaling, EDGE_HITS_BARREL case: scale-center-to-shell distance is out of foreseen bounds:')
                log.warning('   distance=%s, scale_center_to_shell_y=%s, scale_center_to_shell_z=%s' \
                             % (scale_center_shell_dist, scale_center_to_shell_y, scale_center_to_shell_z) )

            if scale_angle <= Pi/4.0:
                cos_scale_angle = math.cos(scale_angle)
                z1_new = min(z1, (scale_center_z + cos_scale_angle * scale_center_shell_dist))
                r_new  = min(r, r1_function(z1_new))
            else:
                sin_scale_angle = math.sin(scale_angle)
                r_new  = min(r,  (scale_center_r + sin_scale_angle * scale_center_shell_dist))
                z1_new = min(z1, z1_function(r_new))

        elif situation == BARREL_HITS_BARREL:
            # The scaling cylinder hits the barrel of 'shell' with its barrel
            if scale_angle == 0.0:
                # In this case the cylinders can never hit; just leave the lenghts as they are
                r_new = r
                z1_new = z1
            else:
                # case scale_angle > 0
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

    if( r < 0 or z_left < 0 or z_right < 0 ):
        raise testShellError('get_dr_dzright_dzleft_to_CylindricalShape: Negative length in shell scaling: r=%s, z_left=%s, z_right=%s' % (r, z_left, z_right) )

    return r, z_right, z_left

def get_radius_to_CylindricalShape(shape, testShell, r):
    # This function returns the radius for the spherical 'testShell' using the cylindrical 'shell' as its closest
    # neighbor. Note that it returns the minimum of the newly calculated radius and the old radius. The 'testShell' can
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

def get_dr_dzright_dzleft_to_DiskShape(shape, testShell, r, z_right, z_left):
    # This function returns the dr, dz_right, dz_left for the cylindrical 'testShell' using the disk 'shape' as its closest
    # neighbor. Note that it returns the minimum of the newly calculated parameters and the old parameters. The 'testShell' can
    # therefore only become smaller.

    assert (type(shape) is Disk)

    # For now just consider the Disk as a special case of a cylinder
    # TODO implement properly, I'm sure things can be more optimal
    return get_dr_dzright_dzleft_to_CylindricalShape(Cylinder(shape.position, shape.radius, shape.unit_z, 0.0), testShell, r, z_right, z_left)

def get_radius_to_DiskShape(shape, testShell, r):
    # This function returns the radius for the spherical 'testShell' using the disk 'shape' as its closest
    # neighbor. Note that it returns the minimum of the newly calculated radius and the old radius. The 'testShell' can
    # therefore only become smaller.

    assert (type(shape) is Disk)

    # Do simple distance calculation to disk
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
        neighbor_surfaces = get_neighbor_surfaces(self.world, searchpoint, base_structure_id, ignores)
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
            self.radius = self.determine_possible_shell(self.origin_structure.id, [self.single.domain_id],
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
        return self.distance_to_target_structure + self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR

#####
class PlanarSurfaceTransitionPairtestShell(SphericaltestShell, testPlanarSurfaceTransitionPair):

    def __init__(self, single1, single2, geometrycontainer, domains):
        SphericaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition

        try:
            testPlanarSurfaceTransitionPair.__init__(self, single1, single2)
        except protoDomainError as e:
            raise testShellError('(PlanarSurfaceTransitionPair). %s' %
                                 (str(e)))

        self.center = self.com
        try:
            self.radius = self.determine_possible_shell(self.structure1.id, [self.single1.domain_id, self.single2.domain_id],
                                                        [self.structure2.id])
        except ShellmakingError as e:
            raise testShellError('(PlanarSurfaceTransitionPair). %s' %
                                 (str(e)))

    def get_min_radius(self):
        return self.get_min_pair_size()


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

        ## TESTING
        #assert dr >= min_dr and dz_right >= min_dz_right and dz_left >= min_dz_left, \
               #'CylindricaltestShell dimensions smaller than the minimum, \
                #dr = %s, dz_right = %s, dz_left = %s, min_dr = %s, min_dz_right = %s, min_dz_left = %s.' % \
               #(dr, dz_right, dz_left, min_dr, min_dz_right, min_dz_left)
 
        # first check the maximum dr, dz_right, dz_left against the surfaces
        # NOTE: we assume that all relevant surfaces are cylindrical and parallel to the testCylinder
        # or planar surfaces and parallel to the testCylinder axis
        for surface, distance in neighbor_surfaces:
            # TODO
            if isinstance(surface, SphericalSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_SphericalShape(surface.shape, self, dr, dz_right, dz_left)
            elif isinstance(surface, CylindricalSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_CylindricalShape(surface.shape, self, dr, dz_right, dz_left)            
            elif isinstance(surface, PlanarSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_PlanarShape(surface.shape, self, dr, dz_right, dz_left)
            elif isinstance(surface, DiskSurface):
                dr, dz_right, dz_left = get_dr_dzright_dzleft_to_DiskShape(surface.shape, self, dr, dz_right, dz_left)            
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
                    dr, dz_right, dz_left = get_dr_dzright_dzleft_to_SphericalShape(shell_appearance.shape, self, dr, dz_right, dz_left)
                elif isinstance(shell_appearance.shape, Cylinder):
                    dr, dz_right, dz_left = get_dr_dzright_dzleft_to_CylindricalShape(shell_appearance.shape, self, dr, dz_right, dz_left)
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
    # property self.right_scalingangle set by constructor of derived classes

    def get_left_scalingangle(self):
        if self.drdz_left == numpy.inf:
            left_scaling_angle = Pi/2.0       # atan(infinity) == Pi/2 -> angle is 90 degrees
        elif self.drdz_left == 0.0:
            left_scaling_angle = 0.0          # atan(0.0) = 0.0 -> faster
        else:
            left_scaling_angle = math.atan(self.drdz_left)
        return left_scaling_angle
    # property self.left_scalingangle set by constructor of derived classes

    def get_tan_right_scalingangle(self):
        return self.tan_right_scalingangle

    def get_tan_left_scalingangle(self):
        return self.tan_left_scalingangle

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

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

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
        dz_left  = dz_right
        dr_sr    = math.sqrt((self.get_searchradius())**2 - dz_right**2) # stay within the searchradius
#        dr_edge  = self.structure.min_dist_proj_to_edge(self.get_referencepoint()) # TODO rethink this
#        dr       = min(dr_sr, dr_edge)
        dr = dr_sr
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

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

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
        dz_left  = dz_right
        dr_sr    = math.sqrt((self.get_searchradius())**2 - dz_right**2) # stay within the searchradius
#        dr_edge  = self.structure.min_dist_proj_to_edge(self.get_referencepoint()) # TODO rethink this
#        dr       = min(dr_sr, dr_edge)
        dr = dr_sr
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r/SAFETY, z_right, z_left

#####
class CylindricalSurfaceSingletestShell(CylindricaltestShell, testNonInteractionSingle):

    def __init__(self, pid_particle_pair, structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)
        testNonInteractionSingle.__init__(self, pid_particle_pair, structure)

        # Initialize the scaling parameters
        self.dr_const   = self.pid_particle_pair[1].radius * CYLINDER_R_FACTOR

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.dr_const
        self.z0_left    = 0.0

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # Sizing up the shell to a zero shell
        # Important: leave this as it is, because we always want to construct
        # a zero-single by default!
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
        dr       = self.dr_const
        dz_right = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        dz_left  = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dr = self.dr_const
        #dz_right_edge = self.structure.min_dist_proj_to_edge(self.get_referencepoint())    # TODO rethink this
        ##_, (_, pp_dist) = self.structure.project_point(self.get_referencepoint())
        ##assert(pp_dist <= 0.0)
        ##dz_right_edge = abs(pp_dist)
        dz_right_sr   = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
        #dz_right      = min(dz_right_sr, dz_right_edge)
        dz_right      = dz_right_sr
        dz_left       = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

class DiskSurfaceSingletestShell(CylindricaltestShell, testNonInteractionSingle):

    def __init__(self, pid_particle_pair, structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)
        testNonInteractionSingle.__init__(self, pid_particle_pair, structure)

        # Initialize scaling parameters
        # - the reference point is at the cap disk position
        # - the orientation vector points from the cap towards the rod center, i.e.:
        #       - the right scaling point is for the scaling "inwards", i.e. onto the rod
        #       - the left scaling point is for the scaling "outwards", i.e. away from the rod (and cap)
        #       - dz_left stays constant and is equal to the cap-bound particle radius
        #       - dz_right can be scaled, minimum is set accordingly
        # - dr stays constant and is determined by the maximal radius involved
        self.dr_const   = self.pid_particle_pair[1].radius*CYLINDER_R_FACTOR

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.dr_const
        self.z0_left    = 0.0

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # Sizing up the shell to a zero shell
        # Important: leave this as it is, because we always want to construct
        # a zero-single by default!
        self.dz_right = self.pid_particle_pair[1].radius
        self.dz_left  = self.pid_particle_pair[1].radius
        self.dr       = self.pid_particle_pair[1].radius

    def get_orientation_vector(self):
        return self.structure.shape.unit_z
        # just copy from disk structure

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.pid_particle_pair[1].position

    def get_min_dr_dzright_dzleft(self):
        # TODO This will never be called, right? Why do dz_right/dz_left have value larger than particle_radius?
        dr       = self.dr_const
        dz_right = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        dz_left  = dz_right
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        # Radius is not scaled here so we do not to check for distance to shape edge
        dr       = self.dr_const
        dz_right = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0) # same as the minimum, i.e. no scaling of this length
        dz_left  = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class CylindricalSurfacePairtestShell(CylindricaltestShell, testSimplePair):

    def __init__(self, single1, single2, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testSimplePair.__init__(self, single1, single2)

        # initialize scaling parameters
        self.dr_const   = max(self.pid_particle_pair1[1].radius, self.pid_particle_pair2[1].radius) * CYLINDER_R_FACTOR

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.dr_const
        self.z0_left    = 0.0

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

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
        dr       = self.dr_const
        dz_right = self.get_min_pair_size()
        dz_left  = dz_right
        return dr, dz_right, dz_left

    def get_max_dr_dzright_dzleft(self):
        dr          = self.dr_const
#        dz_right_edge = self.structure.min_dist_proj_to_edge(self.get_referencepoint())    # TODO rethink
        dz_right_sr = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
#        dz_right    = min(dz_right_sr, dz_right_edge)
        dz_right    = dz_right_sr
        dz_left     = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class PlanarSurfaceInteractiontestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, target_structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, target_structure)

        # initialize scaling parameters
        self.dzdr_right = 1.0
        self.drdz_right = 1.0
        self.r0_right   = 0.0
        self.z0_right   = self.particle_surface_distance # calculated in the __init__ of testInteractionSingle
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.pid_particle_pair[1].radius

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.origin_structure.id, [self.single.domain_id], [self.target_structure.id])
        except ShellmakingError as e:
            raise testShellError('(PlanarSurfaceInteraction). %s' %
                                 (str(e)))

    def get_orientation_vector(self):
        return self.reference_vector     # from testInteractionSingle; reference vector always points from the
                                         # particle position projected onto the target surface towards particle pos.

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
        dr_sr    = (self.get_searchradius()/math.sqrt(2))
        dr_edge  = self.min_dist_proj_to_edge
        dr       = min(dr_sr, dr_edge)
        dz_right = dr + self.particle_surface_distance          # Note that we assume that dr >> particle_surface_distance        
        dz_left  = self.pid_particle_pair[1].radius             # so that the cylinder always fits inside the search_radius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return self.r_right(z_right/SAFETY)/SAFETY, z_right/SAFETY, z_left

#####
class CylindricalSurfaceInteractiontestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, target_structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, target_structure)

        # initialize scaling parameters
        self.dzdr_right = 1.0
        self.drdz_right = 1.0
        self.r0_right   = 0.0
        self.z0_right   = -(self.particle_surface_distance + self.target_structure.shape.radius)   # note that this is because dzdr_right = 1.0
        self.dzdr_left  = self.dzdr_right
        self.drdz_left  = self.drdz_right
        self.r0_left    = self.r0_right
        self.z0_left    = self.z0_right

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.origin_structure.id, [self.single.domain_id], [self.target_structure.id])
        except ShellmakingError as e:
            raise testShellError('(CylindricalSurfaceInteraction). %s' %
                                 (str(e)))


    def get_orientation_vector(self):
        return self.target_structure.shape.unit_z

    def get_searchpoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_referencepoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_min_dr_dzright_dzleft(self):
        dr       = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR + (self.particle_surface_distance + self.target_structure.shape.radius)
        dz_right = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        dz_left  = dz_right
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        # TODO clarify
        dz_right_edge = self.min_dist_proj_to_edge
        dz_right_sr   = ((-(self.particle_surface_distance + self.target_structure.shape.radius) + \
                       math.sqrt(2*(self.get_searchradius()**2) - (self.particle_surface_distance + self.target_structure.shape.radius)**2)) / 2.0)
        dz_right      = min(dz_right_sr, dz_right_edge)
        dz_left       = dz_right
        dr            = ((self.particle_surface_distance + self.target_structure.shape.radius) + dz_right)
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r/SAFETY, self.z_right(r/SAFETY)/SAFETY, self.z_left(r/SAFETY)/SAFETY

class CylindricalSurfaceCapInteractiontestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, target_structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, target_structure)

        # Initialize scaling parameters
        # - the reference point is at the cap disk position
        # - the orientation vector points from the cap towards the rod particle, i.e.:
        #       - the right scaling point is for the scaling "inwards", i.e. onto the rod
        #       - the left scaling point is for the scaling "outwards", i.e. away from the rod (and cap)
        #       - dz_left stays constant and is equal to the cap-bound particle radius * a safety factor
        #       - dz_right can be scaled, minimum is set accordingly
        # - dr stays constant and is determined by the maximal radius involved
        self.dr_const   = self.pid_particle_pair[1].radius * CYLINDER_R_FACTOR

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.dr_const
        self.z0_left    = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.origin_structure.id, [self.single.domain_id], [self.target_structure.id])
        except ShellmakingError as e:
            raise testShellError('(CylindricalSurfaceCapInteraction). %s' %
                                 (str(e)))

    def get_orientation_vector(self):
        # The orientation vector is the (normalized) vector that points from the cap
        # towards the particle on the rod
        return normalize(self.pid_particle_pair[1].position - self.target_structure.shape.position)

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.reference_point         # calculated in the __init__ of testInteractionSingle

    def get_min_dr_dzright_dzleft(self):
        dr       = self.dr_const
        dz_left  = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        dz_right = self.particle_surface_distance + self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        dr       = self.dr_const
        dz_left  = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0) # same as the minimum, i.e. no scaling of this length
        # TODO include max distance to other side of rod?
        dz_right = math.sqrt((self.get_searchradius())**2 - dr**2) + self.particle_surface_distance  # stay within the searchradius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left

class CylindricalSurfaceSinktestShell(CylindricaltestShell, testInteractionSingle):

    def __init__(self, single, target_structure, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testInteractionSingle.__init__(self, single, single.structure, target_structure)

        # initialize scaling parameters
        self.dr_const   = self.pid_particle_pair[1].radius * CYLINDER_R_FACTOR

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = self.dzdr_right
        self.drdz_left  = self.drdz_right
        self.r0_left    = self.dr_const
        self.z0_left    = self.z0_right

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.origin_structure.id, [self.single.domain_id], [self.target_structure.id])
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
        pos = self.world.cyclic_transpose(self.pid_particle_pair[1].position, self.target_structure.shape.position)
        orientation_vector = pos - self.target_structure.shape.position
        # if the dot product is negative, then take the -unit_z, otherwise take unit_z
        return self.target_structure.shape.unit_z * cmp(numpy.dot(orientation_vector, self.target_structure.shape.unit_z), 0)

    def get_searchpoint(self):
        return self.pid_particle_pair[1].position

    def get_referencepoint(self):
        return self.pid_particle_pair[1].position

    def get_min_dr_dzright_dzleft(self):
        dr       = self.dr_const
        dz_left  = self.particle_surface_distance + self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        dz_right = dz_left              # we make the domain to be symmetrical around the particle
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        dr       = self.dr_const
#        dz_right_edge = self.origin_structure.min_dist_proj_to_edge(self.get_referencepoint()) # TODO rethink
        dz_right_sr   = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
#        dz_right      = min(dz_right_sr, dz_right_edge)
        dz_right      = dz_right_sr
        dz_left       = dz_right
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left/SAFETY

#####
class MixedPair2D3DtestShell(CylindricaltestShell, testMixedPair2D3D):

    def __init__(self, single2D, single3D, geometrycontainer, domains):
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        # The initialization of r0 can fail in testPair.__init__

        try:
            testMixedPair2D3D.__init__(self, single2D, single3D)
        except protoDomainError as e:
            raise testShellError('(MixedPair2D3D). %s' %
                                 (str(e)))

        # Initialize commonly used constants
        self.sqrt_DRDr = math.sqrt((2*self.D_R)/(3*self.D_r))

        # Initialize the scaling parameters
        # To determine the scale angle we calculate the radius r corresponding to the minimal z_right
        # of the cylinder in a way that the first passage time of the CoM (towards r) and of the
        # interparticle vector (towards z_right) are equal (see self.r_right(), self.z_right() below).
        self.iv_z_minimum = self.r0 + self.particle3D.radius * SINGLE_SHELL_FACTOR
            # this is an upper bound, actually the height of the domain could be slightly lower than this
        self.drdz_right = self.r_right(self.iv_z_minimum) / self.iv_z_minimum
        self.dzdr_right = 1.0/self.drdz_right        

        self.r0_right   = 0.0
        self.z0_right   = self.z_right(self.r0_right)
        self.dzdr_left  = 0.0
        self.drdz_left  = numpy.inf
        self.r0_left    = 0.0
        self.z0_left    = self.particle2D.radius

        # TODO This is probably not needed any more:
        min_r, min_dz_right, _ = self.get_min_dr_dzright_dzleft()

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        ### TESTING TESTING TESTING remove when done
        log.info(' MixedPair2D3D: Initalized shellscaling: scale_angle = %s, drdz  = %s' % (self.right_scalingangle, self.drdz_right))
        log.info('                                         r0_right = %s, z0_right = %s' % (self.r0_right, self.z0_right))
        log.info('                                         r0_left  = %s, z0_left  = %s' % (self.r0_left,  self.z0_left))

        log.info('      min_dr_dzright_dzleft = %s' % str(self.get_min_dr_dzright_dzleft() ))
        log.info('      max_dr_dzright_dzleft = %s' % str(self.get_max_dr_dzright_dzleft() ))


        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure3D.id, [self.single2D.domain_id, self.single3D.domain_id],
                                                          [self.structure2D.id])
        except ShellmakingError as e:
            raise testShellError('(MixedPair2D3D). %s' %
                                 (str(e)))


    def get_orientation_vector(self):
        # if the dot product is negative, then take the -unit_z, otherwise take unit_z
        return self.structure2D.shape.unit_z * cmp(numpy.dot(self.iv, self.structure2D.shape.unit_z), 0)

    def get_searchpoint(self):
        # TODO this can be improved, now same as PlanarSurfaceInteraction
        search_distance = self.get_searchradius()/math.sqrt(2.0) - self.z0_left
        return self.get_referencepoint() + search_distance * self.get_orientation_vector()

    def get_referencepoint(self):
        return self.com

    def get_min_dr_dzright_dzleft(self):
        # defines the minimum parameters of the testShell.
        # Note that also the minimum testShell must obey the scaling rules

        radius2D = self.particle2D.radius
        radius3D = self.particle3D.radius
        D_2D = self.particle2D.D
        D_3D = self.particle3D.D


        ### calculate the minimal height z_right1 of the shell including burst radius
        dz_right1 = self.iv_z_minimum
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

        dz_left  = self.particle2D.radius

        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        # TODO this can be improved
        dr_sr    = self.get_searchradius()/math.sqrt(2)
        dr_edge  = self.min_dist_proj_to_edge
        dr       = min(dr_sr, dr_edge)
        dz_left  = self.particle2D.radius

        # now adjust dz_right such that they obey the scaling laws and are still within the search_radius
        # We assume that dz_right needs to be adjusted to dr and not the other way around since the cylinder is
        # always a bit flatter.
        dz_right = self.z_right(dr)        
        ## Debug output for TESTING
        #log.debug( "max_dr=%s, max_dz_right=%s, max_dz_left=%s" % (dr, dz_right, dz_left) )
        #log.debug( "min_dist_proj_to_edge=%s" % str(self.min_dist_proj_to_edge) )
        #assert dr > dz_right   ## Removed for TESTING

        return dr, dz_right, dz_left

    def z_right(self, r_right):
        # This calculates the height z_right of the newly constructed domain
        # (above the plane) if the radius r is known

        # The trivial case
        if r_right == 0.0:
            return 0.0

        # Some constants that we need
        radius2D = self.particle2D.radius
        radius3D = self.particle3D.radius
        D_2D = self.particle2D.D
        D_3D = self.particle3D.D

        # calculate a_r such that the expected first-passage for the CoM and IV are equal        
        #sqrt_DRDr_2D = math.sqrt( 2.0/3.0 * self.D_R/D_2D )  # TESTING
        #sqrt_DRDr_3D = math.sqrt( 2.0/3.0 * self.D_R/D_3D )
        a_r_2D = (r_right - radius2D + self.r0*self.sqrt_DRDr ) / (self.sqrt_DRDr + (D_2D/self.D_tot) )
        a_r_3D = (r_right - radius3D + self.r0*self.sqrt_DRDr ) / (self.sqrt_DRDr + (D_3D/self.D_tot) )        
        # take the smaller a_r that, if entered in the function for r below, would lead to this r_right
        a_r = min(a_r_2D, a_r_3D)

        z_right = (a_r/self.get_scaling_factor()) + radius3D

        ## Debug output for TESTING
        #log.debug( " a_r_2D=%s, a_r_3D=%s, z_right=%s, r_right=%s" % (a_r_2D, a_r_3D, z_right, r_right))
        #log.debug('r_right= %s, r_right(z_right)=%s, z_right=%s' % (r_right, self.r_right(z_right), z_right) )

        #assert feq(self.r_right(z_right), r_right), 'Inconsistent z_right function in MixedPair2D3DtestShell: \
                                                     #r_right= %s, r_right(z_right)=%s, z_right=%s' % \
                                                     #(r_right, self.r_right(z_right), z_right)

        if r_right < 0.0 or z_right <0.0:
            log.warning('Negative cylinder dimensions in MixedPair2D3DtestShell: z_right= %s, r_right=%s' % (z_right, r_right) )

        if abs(z_right) < abs(r_right) * TOLERANCE:
            log.warning('Setting z_right to zero within TOLERANCE in MixedPair2D3DtestShell, z_right=%s' % str(z_right))
            z_right = 0.0

        return z_right

    def r_right(self, z_right):
        # This calculates the radius r of the newly constructed domain
        # if the height z_right is known

        # The trivial case
        if z_right == 0.0:
            return 0.0

        # Some constants that we need
        radius2D = self.particle2D.radius
        radius3D = self.particle3D.radius
        D_2D = self.particle2D.D
        D_3D = self.particle3D.D

        # We first calculate the a_r, since it is the only radius that depends on z_right only.
        # The dependence is given by the scaling factor that converts spherical coordinates into
        # prolate spheroidal coordinates
        a_r = (z_right - radius3D) * self.get_scaling_factor()

        # We equalize the estimated first passage time for the CoM (2D) and IV (3D) for a given a_r
        # via a_R/(a_r-self.r0) == sqrt(4 D_R t)/sqrt(6 D_r t) (for an arbitrary t).
        # This gives us a CoM domain radius a_R as a function of a_r at equal expected first passage time.
        # Note that for the IV we only take the distance to the outer boundary into account.
        a_R = (a_r - self.r0)*self.sqrt_DRDr

        # We calculate the maximum space needed for particles A and B based on maximum IV displacement
        # taking into account the particle radii
        iv_max = max( (D_2D/self.D_tot * a_r + radius2D),
                      (D_3D/self.D_tot * a_r + radius3D))

        # The total domain radius is then given by the sum of the max. IV displacement and
        # the estimated a_R value for the same first passage time
        r_right = a_R + iv_max

        ### TESTING
        log.info('MixedPair2D3DtestShell: z_right= %s, z_right(r_right)=%s, r_right=%s' % (z_right, self.z_right(r_right), r_right) )

        if r_right < 0.0 or z_right <0.0:
            log.warning('Negative cylinder dimensions in MixedPair2D3DtestShell: z_right= %s, r_right=%s' % (z_right, r_right) )

        assert feq(self.z_right(r_right), z_right), 'Inconsistent r_right function in MixedPair2D3DtestShell: \
                                                     z_right= %s, z_right(r_right)=%s, r_right=%s' % \
                                                     (z_right, self.z_right(r_right), r_right)

        if abs(r_right) < abs(z_right) * TOLERANCE:
            log.warning('Setting r_right to zero within TOLERANCE in MixedPair2D3DtestShell, r_right=%s' % str(r_right))
            r_right = 0.0

        return r_right

    def apply_safety(self, r, z_right, z_left):
        SAFETY = 1.1 # FIXME What the hell is that? Why is that not the standard?
        # Make sure to keep the right scaling ratio when applying SAFETY
        # because we determine the math. domains for IV and CoM from the shell dimensions
        return r/SAFETY, self.z_right(r/SAFETY), z_left


class MixedPair1DCaptestShell(CylindricaltestShell, testMixedPair1DCap):

    def __init__(self, single1D, cap_single, geometrycontainer, domains): # TODO Do we need target_structure (for the ignore list)?
        CylindricaltestShell.__init__(self, geometrycontainer, domains)  # this must be first because of world definition
        testMixedPair1DCap.__init__(self, single1D, cap_single)

        # NOTE: The following is already defined by testMixedPair1DCap:
        # - self.single1D, self.structure1D, self.particle1D
        # - self.cap_single, self.cap_structure, self.cap_particle

        # Initialize scaling parameters
        # - the reference point is at the cap disk position
        # - the orientation vector points from the cap towards the rod particle, i.e.:
        #       - the right scaling point is for the scaling "inwards", i.e. onto the rod
        #       - the left scaling point is for the scaling "outwards", i.e. away from the rod (and cap)
        #       - dz_left stays constant and is equal to the cap-bound particle radius * a safety factor
        #       - dz_right can be scaled, minimum is set accordingly
        # - dr stays constant and is determined by the maximal radius involved
        CRF             = CYLINDER_R_FACTOR
        self.dr_const   = max([CRF * self.cap_particle.radius, CRF * self.particle1D.radius])

        self.dzdr_right = numpy.inf
        self.drdz_right = 0.0
        self.r0_right   = self.dr_const
        self.z0_right   = 0.0
        self.dzdr_left  = numpy.inf
        self.drdz_left  = 0.0
        self.r0_left    = self.dr_const
        self.z0_left    = self.cap_particle.radius * SINGLE_SHELL_FACTOR

        # Now we can also define the scaling angle
        self.right_scalingangle = self.get_right_scalingangle()
        self.left_scalingangle  = self.get_left_scalingangle()
        # In particular store the tangent, because math.tan is expensive!
        self.tan_right_scalingangle = math.tan(self.right_scalingangle)
        self.tan_left_scalingangle  = math.tan(self.left_scalingangle)

        # This will determine if the shell is possible.
        # If possible, it will write the dr, dz_right, dz_left defining the dimensions of the cylindrical shell.
        # If not possible, it throws an exception and the construction of the testShell IS ABORTED!
        try:
            self.dr, self.dz_right, self.dz_left = \
                            self.determine_possible_shell(self.structure1D.id, [self.single1D.domain_id, self.cap_single.domain_id], [self.cap_structure.id])
        except ShellmakingError as e:
            raise testShellError('(MixedPair1DCap). %s' %
                                 (str(e)))

    def get_orientation_vector(self):        
        # The orientation vector is the (normalized) vector that points from the cap
        # towards the particle on the rod
        return normalize(self.particle1D.position - self.cap_structure.shape.position)

    def get_searchpoint(self):
        return self.particle1D.position

    def get_referencepoint(self):        
        return self.cap_structure.shape.position

    def get_min_dr_dzright_dzleft(self):
        dr       = self.dr_const
        dz_left  = self.cap_particle.radius * SINGLE_SHELL_FACTOR
        dz_right = self.get_min_pair_size() # TODO make sure this does the right thing!
        return dr, dz_right, dz_left
        
    def get_max_dr_dzright_dzleft(self):
        dr       = self.dr_const
        dz_left  = self.cap_particle.radius * SINGLE_SHELL_FACTOR # same as the minimum, i.e. no scaling of this length
        dz_right = math.sqrt((self.get_searchradius())**2 - dr**2) # stay within the searchradius
        return dr, dz_right, dz_left

    def apply_safety(self, r, z_right, z_left):
        return r, z_right/SAFETY, z_left

    def get_min_pair_size(self):
        # This calculates the minimal 'size' of a simple pair.
        # Note that the 'size' is interpreted differently dependent on the type of pair. When the pair is a
        # SphericalPair          -> size = radius
        # PlanarSurfacePair      -> size = radius
        # CylindricalSurfacePair -> size = half_length

        D_1 = self.particle1D.D
        D_2 = self.cap_particle.D
        radius1 = self.particle1D.radius
        radius2 = self.cap_particle.radius

        dist_from_com1 = self.r0 * D_1 / self.D_tot     # particle distance from CoM
        dist_from_com2 = self.r0 * D_2 / self.D_tot
        iv_shell_size1 = dist_from_com1 + radius1       # the shell should surround the particles
        iv_shell_size2 = dist_from_com2 + radius2

        # also fix the minimum shellsize for the CoM domain
        com_shell_size = max(radius1, radius2)

        # calculate total radii including the margin for the burst volume for the particles
        shell_size = max(iv_shell_size1 + com_shell_size + radius1 * SINGLE_SHELL_FACTOR,
                         iv_shell_size2 + com_shell_size + radius2 * SINGLE_SHELL_FACTOR)

        return shell_size

#####
#class MixedPair3D1DtestShell(CylindricaltestShell, Others):
