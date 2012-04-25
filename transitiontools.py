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
    'TransitionTools'
    ]


class TransitionTools :

    def __init__(self, testShell, start_position):

        assert isinstance(testShell, testTransitionSingle);

        # Might be used internally by classes inheriting from here
        self.origin_structure  = testShell.origin_structure
        self.target_structure  = testShell.target_structure
        self.structure1        = testShell.structure1
        self.structure2        = testShell.structure2

        self.origin_center = self.origin_structure.shape.position
        self.target_center = self.target_structure.shape.position
        self.origin_half_extent = self.origin_structure.shape.half_extent
        self.target_half_extent = self.target_structure.shape.half_extent

        self.start_position  = start_position
        self.target_distance = testShell.distance_to_target_structure
        

    def process_new_position_vector(self, old_pos, displacement):
                
        # Construct the new position using the deflection function of the target surface            
        new_pos, changeflag = self.target_structure.deflect(old_pos, displacement)        
                
        if changeflag==0 :
            # The new position is still in the old surface;
            new_structure_id = self.origin_structure.id  # i.e. no change here
        else :
            # The new position is out of the old surface;            
            new_structure_id = self.target_structure.id

        return new_pos, new_structure_id


