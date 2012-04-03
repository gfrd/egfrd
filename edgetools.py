from _gfrd import (
    SphericalShell,
    CylindricalShell,
    Sphere,
    Cylinder,
    CuboidalRegion,
    PlanarSurface,
    CylindricalSurface,
    )

from utils import *

from domain import (
    Domain,
    ProtectiveDomain)

from shells import *

__all__ = [
    'EdgeTools'
    ]


class EdgeTools :

    def __init__(self, testShell, start_position):

        assert isinstance(testShell, testTransitionSingle);

        # Might be used internally by classes inheriting from here
        self.origin_structure  = testShell.origin_structure
        self.target_structure  = testShell.target_structure

        self.origin_center = self.origin_structure.shape.position
        self.target_center = self.target_structure.shape.position
        self.origin_half_extent = self.origin_structure.shape.half_extent
        self.target_half_extent = self.target_structure.shape.half_extent

        self.start_position  = start_position
        self.target_distance = testShell.distance_to_target_structure

        self.edge_point = self.get_edge_point()
        self.origin_unit_perp, self.origin_unit_par, self.origin_half_extent_perp = self.get_origin_unit_vectors()
        self.target_unit_perp, self.target_unit_par, self.target_half_extent_perp = self.get_target_unit_vectors()        
        

    def process_new_position_vector(self, oldpos, displacement):
        # This routine checks whether the new position in the surface of origin
        # reaches out of it and if necessary calculates the position in the
        # target surface.
        # First calculate the coordinate perpendicular to the edge
        # in the old surface
        newpos_in_origin = (oldpos - self.origin_center) + displacement
        d_origin_perp = numpy.dot(newpos_in_origin, self.origin_unit_perp)
        assert (d_origin_perp >= 0.0 )
        # Calculate the part of it that reaches out of the surface of origin
        d_out = d_origin_perp - self.origin_half_extent_perp
        if d_out <= 0.0 :
            # The new position is still in the old surface;
            # no transform required
            newpos = oldpos + displacement  
            new_structure_id = self.origin_structure.id  # i.e. no change here
        else :
            # The new position is out of the old surface;
            # transform it to the target surface.
            # First calculate the coordinate in the target surface
            # that is perpendicular to the edge:
            d_target_perp = self.target_half_extent_perp - d_out
            # Now the same for the parallel component
            # Take into account that the parallel unit vector in the target surface might be
            # either parallel or antiparallel to the parallel unit vector in the origin surface
            d_origin_par = numpy.dot(newpos_in_origin, self.origin_unit_par)
            d_target_par = d_origin_par * numpy.dot(self.origin_unit_par, self.target_unit_par)
            # Construct the new position using the predefined unit vectors of the target surface
            newpos = self.target_center + d_target_perp*self.target_unit_perp + d_target_par*self.target_unit_par
            new_structure_id = self.target_structure.id

        return newpos, new_structure_id

    def get_origin_unit_vectors(self):
        # calculate the unit vector perpendicular and parallel to the edge
        # in the surface of origin and the half_extent of this surface in the
        # perpendicular vector direction
        # ATTENTION! THIS DOES NOT YET WORK WITH PERIODIC BC!
        w_z = self.target_structure.shape.unit_z
        u_perp = (1.0/self.target_distance * numpy.dot(self.target_center-self.start_position, w_z )) * w_z
        assert feq( length(u_perp), 1.0) 

        u_x = self.origin_structure.shape.unit_x
        u_y = self.origin_structure.shape.unit_y
        if feq(abs(numpy.dot(u_x, u_perp)), 1.0):
             u_par  = u_y
             h_perp = self.origin_half_extent[0] # half_extent in u_x direction
        else:
             u_par  = u_x
             h_perp = self.origin_half_extent[1] # half_extent in u_y direction

        return u_perp, u_par, h_perp

    def get_target_unit_vectors(self):
        # calculate the unit vector perpendicular and parallel to the edge
        # in the target surface and the half_extent of this surface in the
        # perpendicular vector direction
        # ATTENTION! THIS DOES NOT YET WORK WITH PERIODIC BC!
        w_perp_unnorm = (self.origin_center - self.target_center) + self.origin_half_extent_perp * self.origin_unit_perp
        w_perp = ( 1.0/numpy.sqrt(numpy.dot(w_perp_unnorm,w_perp_unnorm)) ) * w_perp_unnorm
        assert feq(numpy.sqrt(numpy.dot(w_perp, w_perp)), 1.0)

        w_x = self.target_structure.shape.unit_x
        w_y = self.target_structure.shape.unit_y
        if feq(abs(numpy.dot(w_x, w_perp)), 1.0):
             w_par = w_y
             h_perp = self.target_half_extent[0] # half_extent in w_x direction
        else:
             w_par = w_x
             h_perp = self.target_half_extent[1] # half_extent in w_y direction

        return w_perp, w_par, h_perp

    def get_edge_point(self):
        self.origin_unit_perp, _, _ = self.get_origin_unit_vectors()
        ep = self.start_position + (self.target_distance * self.origin_unit_perp)
        return ep

