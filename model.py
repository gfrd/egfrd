import _gfrd
from _gfrd import create_cuboidal_region, create_cylindrical_surface, \
        create_disk_surface, create_planar_surface, create_double_sided_planar_surface
import numpy

__all__ = [
    'Species',
    'ParticleModel',
    'create_unimolecular_reaction_rule',
    'create_decay_reaction_rule',
    'create_annihilation_reaction_rule',
    'create_binding_reaction_rule',
    'create_unbinding_reaction_rule',
    'create_surface_absorption_reaction_rule',
    'create_surface_binding_reaction_rule',
    'create_surface_unbinding_reaction_rule',

    # From _gfrd. Should be part of the world class.
    'create_cuboidal_region',
    'create_cylindrical_surface',
    'create_disk_surface',
    'create_planar_surface',
    'create_double_sided_planar_surface',
    ]


# Define _gfrd docstrings here, much easier to format than in C++.
# TODO These functions should be moved to gfrdbase.py because structures are no longer a part of the model,
# but part of the world instead.
_gfrd.create_cuboidal_region.__doc__ = \
"""create_cuboidal_region(sid, id, corner, diagonal)

Create and return a new cuboidal Region.

Arguments:
    - sid
        the structure type of the cuboidal region.
    - id
        the name of the structure
    - corner
        the point [x, y, z] of the cuboidal Region closest to
        [0, 0, 0]. Units: [meters, meters, meters]
    - diagonal
        the vector [x, y, z] from the corner closest to [0, 0, 0], to 
        the corner furthest away from [0, 0, 0]. Units:
        [meters, meters, meters]

"""

_gfrd.create_cylindrical_surface.__doc__ = \
"""create_cylindrical_surface(sid, id, corner, radius, orientation, length)

Create and return a new cylindrical Surface.

Arguments:
    - sid
        the structure type of the cylindrical surface.
    - id
        the name of the structure
    - corner
        the point [x, y, z] on the axis of the cylinder closest to 
        [0, 0, 0]. Units: [meters, meters, meters]
    - radius
        the radius of the cylinder. Units: meters.
    - orientation
        the unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the 
        axis of the cylinder.
    - length
        the length of the cylinder. Should be equal to the world_size. 
        Units: meters.

Surfaces are not allowed to touch or overlap.

"""

_gfrd.create_disk_surface.__doc__ = \
"""create_disk_surface(sid, id, center, radius, orientation)

Create and return a new disk Surface.

Arguments:
    - sid
        the structure type of the cylindrical surface.
    - id
        the name of the structure
    - center
        The center point of the circle defining the disk.
        Units: [meters, meters, meters]
    - radius
        the radius of the disk. Units: meters.
    - orientation
        the unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1]
        perpendicular to the disk. This vector typically defines
        the direction of dissociation from the disk if it is a
        cap.

Surfaces are not allowed to touch or overlap.

"""

_gfrd.create_planar_surface.__doc__ = \
"""create_planar_surface(sid, id, corner, unit_x, unit_y, length_x, length_y)

Create and return a new planar surface with one-sided particle unbinding,
i.e. particles on this surface will unbind in the direction of the normal
vector unit_z (cross product of unit_x and unit_y).

Arguments:
    - sid
        the structure type of the planar surface.
    - id
        the name of the structure
    - corner
        the point [x, y, z] on the plane closest to [0, 0, 0]. Units: 
        [meters, meters, meters]
    - unit_x
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the 
        plane.
    - unit_y
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the plane 
        and perpendicular to unit_x.
    - length_x
        the length of the plane along the unit vector unit_x. Should be 
        equal to the world_size. Units: meters.
    - length_y
        the length of the plane along the unit vector unit_y. Should be 
        equal to the world_size. Units: meters.

Surfaces are not allowed to overlap.

Todo: allow the user to specify the position of a planar surface 
relative to a region.

"""

_gfrd.create_double_sided_planar_surface.__doc__ = \
"""create_double_sided_planar_surface(sid, id, corner, unit_x, unit_y, length_x, length_y)

Create and return a new planar surface with double-sided particle unbinding,
i.e. particles on this surface will unbind randomly in either unit_z or -unit_z
direction (unit_z = cross product of unit_x and unit_y).

Arguments:
    - sid
        the structure type of the planar surface.
    - id
        the name of the structure
    - corner
        the point [x, y, z] on the plane closest to [0, 0, 0]. Units: 
        [meters, meters, meters]
    - unit_x
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the 
        plane.
    - unit_y
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the plane 
        and perpendicular to unit_x.
    - length_x
        the length of the plane along the unit vector unit_x. Should be 
        equal to the world_size. Units: meters.
    - length_y
        the length of the plane along the unit vector unit_y. Should be 
        equal to the world_size. Units: meters.

Surfaces are not allowed to overlap.

Todo: allow the user to specify the position of a planar surface 
relative to a region.

"""

_gfrd.Model.add_species_type.im_func.__doc__ = \
"""add_species_type(self, species)

Add a Species to the ParticleModel.

Arguments:
    - species
        a Species created with the function model.Species.

"""

# Species is just a function, not a class. It seems to be a frontend to the
# constructor of _gfrd.SpeciesType
def Species(name, D, radius=0, structure_type=None, drift=0):
    """Define a new Species (in/on a specific Region or Surface).

    Arguments:
        - name
            the name of this Species.
        - D
            the diffusion constant for this Species in/on this 
            Region or Surface. Units: meters^2/second.
        - radius
            the radius for this Species in/on this Region or Surface. 
            Units: meters.
        - structure_type
            the Region or Surface in/on which this Species can exist.  
            Optional. If you do not specify a Structure the Species is 
            added to the "world".
        - drift
            the drift term for this ParticleType on a 
            CylindricalSurface (1D drift). Units: meters/second. 
            Optional.

    If a certain Species should be able to exist in the "world" as 
    well as in/on one of the previously created Regions or Surfaces, 
    then two distinct Species should be created. One with and one 
    without an explicit Structure argument.

    """
    # The SpeciesType actually only holds the information needed for a general species
    # The information needed for spatial simulations (D, r, v etc) is stored in SpeciesInfo
    # Here it is stored temporarily in string fields and only when the world in create
    # put in the SpeciesInfo -> TODO SpeciesInfo should be in the ParticleModel
    st = _gfrd.SpeciesType()
    st["name"] = str(name)
    st["D"] = str(D)
    st["v"] = str(drift)
    st["radius"] = str(radius)
#    st["structure_type"] = str(id(structure_type))  # The most ugly way of getting the reference in there.

    # Particles of a Species whose Surface is not specified will be 
    # added to the "world". 
    if structure_type:
        st["structure_type"] = structure_type['name']  # The most ugly way of getting the reference in there.
    else:
        st["structure_type"] = "world"

#    # create the associated SpeciesInfo
#    si = _gfrd.SpeciesInfo(st.id, structure_type.id,
#                           D, radius, v)

    return st


class ParticleModel(_gfrd.ParticleModel):
    """
    """
    def __init__(self, world_size):
        """Create a new ParticleModel.

        Arguments:
            - world_size
                the size of one side of the simulation "world". Units: 
                meters.

        The simulation "world" is always assumed to be a cube with 
        *periodic boundary conditions*, with 1 corner at [0, 0, 0] and 
        the corner furthest away from [0, 0, 0] being at
        [world_size, world_size, world_size].

        """
        _gfrd.ParticleModel.__init__(self)
        # Note that in _gfrd.ParticleModel the default structure_type is
        # created and added to the model.

        self.world_size = world_size
        # Dimensions don't matter, except for 
        # visualization.

#        self.structures = {}


### TODO change to add_structure_type
#    def add_structure(self, structure):
#        """Add a Structure (Region or Surface) to the ParticleModel.
#
#        Arguments:
#            - structure
#              a Region or Surface created with one of the functions
#              model.create_<>_region or model.create_<>_surface.
#
#        """
#        assert isinstance(structure, _gfrd.Structure)
#        self.structures[structure.id] = structure
#        return structure
#
#    def get_structure(self, id): 
#        return self.structures[id]

    def add_reaction_rule(self, reaction_rule):
        """Add a ReactionRule to the ParticleModel.

        Argument:
            - reaction rule
                a ReactionRule created by one of the functions
                model.create_<>_reaction_rule.

        """
        self.network_rules.add_reaction_rule(reaction_rule)

    def set_all_repulsive(self):
        """Set all 'other' possible ReactionRules to be repulsive.

        By default an EGFRDSimulator will assume:
            - a repulsive bimolecular reaction rule (k=0) for each 
              possible combination of reactants for which no 
              bimolecular reaction rule is specified. 
          
        This method explicitly adds these ReactionRules to the 
        ParticleModel.

        """
        nr = self.network_rules
        # Maybe the user has defined a reaction rule for any 2 species since a 
        # previous call to this method, so remove *all* repulsive reaction 
        # rules first.
        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is not None:
                    for reaction_rule in gen:
                        if float(reaction_rule['k']) == 0.0:
                            nr.remove_reaction_rule(reaction_rule)

        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is None or len(set(gen)) == 0:
                    rr = _gfrd.ReactionRule([species1, species2], [])
                    rr['k'] = '0.0'
                    nr.add_reaction_rule(rr)

def create_unimolecular_reaction_rule(reactant, product, k):
    """Example: A -> B.

    Arguments:
        - reactant
            a Species.
        - product 
            a Species.
        - k
            reaction rate. Units: per second. (Rough order of magnitude: 
            1e-2 /s to 1e2 /s).

    The reactant and the product should be in/on the same 
    Region or Surface.

    There is no distinction between an intrinsic and an overall reaction 
    rate for a unimolecular ReactionRule.

    A unimolecular reaction rule defines a Poissonian process.

    """
    rr = _gfrd.ReactionRule([reactant], [product])
    rr['k'] = '%.16g' % k
    return rr

def create_decay_reaction_rule(reactant, k):
    """Example: A -> 0.

    Arguments:
        - reactant
            a Species.
        - k
            reaction rate. Units: per second. (Rough order of magnitude: 
            1e-2 /s to 1e2 /s).

    There is no distinction between an intrinsic and an overall reaction 
    rate for a decay ReactionRule.

    A decay reaction rule defines a Poissonian process.

    """
    rr = _gfrd.ReactionRule([reactant], [])
    rr['k'] = '%.16g' % k
    return rr

def create_creation_reaction_rule(product, k):
    """Example: 0 -> A.

    Arguments:
        - product
            a Species
        - k
            reaction rate. Units: per second.

    A creation reaction rule defines a Poissonian process.

    """

    rr = _gfrd.ReactionRule([], [product])
    rr['k'] = '%.16g' % k
    return rr

def create_annihilation_reaction_rule(reactant1, reactant2, ka):
    """Example: A + B -> 0.

    Arguments:
        - reactant1
            a Species.
        - reactant2
            a Species.
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s).

    The reactants should be in/on the same Region or Surface.

    ka should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (kon) to an intrinsic reaction rate (ka) with 
    the function utils.k_a(kon, kD), but only for reaction rules in 3D.

    By default an EGFRDSimulator will assume a repulsive 
    bimolecular reaction rule (ka=0) for each possible combination of 
    reactants for which no bimolecular reaction rule is specified. 
    You can explicitly add these reaction rules to the model with the 
    method model.ParticleModel.set_all_repulsive.

    """
    rr = _gfrd.ReactionRule([reactant1, reactant2], [])
    rr['k'] = '%.16g' % ka
    return rr

def create_binding_reaction_rule(reactant1, reactant2, product, ka):
    """Example: A + B -> C.

    Arguments:
        - reactant1
            a Species.
        - reactant2
            a Species.
        - product
            a Species.
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s)

    The reactants and the product should be in/on the same 
    Region or Surface.

    A binding reaction rule always has exactly one product.

    ka should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (kon) to an intrinsic reaction rate (ka) with 
    the function utils.k_a(kon, kD), but only for reaction rules in 3D.

    By default an EGFRDSimulator will assume a repulsive 
    bimolecular reaction rule (ka=0) for each possible combination of 
    reactants for which no bimolecular reaction rule is specified. 
    You can explicitly add these reaction rules to the model with the 
    method model.ParticleModel.set_all_repulsive.

    """
    rr = _gfrd.ReactionRule([reactant1, reactant2], [product])
    rr['k'] = '%.16g' % ka
    return rr

def create_unbinding_reaction_rule(reactant, product1, product2, kd):
    """Example: A -> B + C.

    Arguments:
        - reactant
            a Species.
        - product1
            a Species.
        - product2
            a Species.
        - kd
            intrinsic reaction rate. Units: per second. (Rough order of 
            magnitude: 1e-2 /s to 1e2 /s).

    The reactant and the products should be in/on the same 
    Region or Surface.

    An unbinding reaction rule always has exactly two products.

    kd should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (koff) for this reaction rule to an intrinsic 
    reaction rate (kd) with the function utils.k_d(koff, kon, kD) or 
    utils.k_d_using_ka(koff, ka, kD).

    An unbinding reaction rule defines a Poissonian process.

    """
    rr = _gfrd.ReactionRule([reactant], [product1, product2])
    rr['k'] = '%.16g' % kd
    return rr

def create_surface_absorption_reaction_rule(reactant, surface, ka):
    """Example: A + some_surface -> 0 + some_surface

    Arguments:
        - reactant
            a Species in the "world" or in a Region.
        - surface
            a Surface created with one of the functions 
            model.create_<>_surface.
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s)

    ka should be an *intrinsic* reaction rate. No analytical expression 
    is currently known to convert an overall reaction rate (kon) to 
    an intrinsic reaction rate (ka) for surface binding or absorption 
    reaction rules.

    By default an EGFRDSimulator will assume a repulsive 
    surface binding reaction rule (ka=0) for each possible 
    combination of reactant and Surface for which no surface 
    binding or absorption reaction rule is specified. You can 
    explicitly add these reaction rules to the model by calling 
    model.ParticleModel.set_all_repulsive().

    """
    assert isinstance(reactant['structure'], _gfrd.BoxShapedRegion)

def create_surface_binding_reaction_rule(reactant, surface, product, ka):
    """Example: A + some_surface -> A_on_surface + some_surface

    Arguments:
        - reactant
            a Species in the "world" or in a Region.
        - surface
            a Surface created with one of the functions 
            model.create_<>_surface.
        - product
            a Species that can exist on the Surface mentioned in the 
            previous argument. For example: if you assigned that 
            Surface to the variable some_surface, then a valid product 
            Species would be:
            model.Species('A_on_surface', some_D, some_r, some_surface)
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s)

    A surface binding reaction rule always has exactly one product.

    ka should be an *intrinsic* reaction rate. No analytical expression 
    is currently known to convert an overall reaction rate (kon) to 
    an intrinsic reaction rate (ka) for surface binding or absorption 
    reaction rules.

    By default an EGFRDSimulator will assume a repulsive 
    surface binding reaction rule (ka=0) for each possible 
    combination of reactant and Surface for which no surface 
    binding or absorption reaction rule is specified. You can 
    explicitly add these reaction rules to the model by calling 
    model.ParticleModel.set_all_repulsive().

    """
    assert product['structure'] == surface

def create_surface_unbinding_reaction_rule(reactant, surface, product, kd):
    """Example: A_on_surface + some_surface -> A + some_surface

    Arguments:
        - reactant
            a Species that can exist on the Surface mentioned in the 
            next argument. For example: if you assign that 
            Surface to the variable some_surface, then a valid reactant 
            Species would be:
            model.Species('A_on_surface', some_D, some_r, some_surface)
        - surface
            a Surface created with one of the functions 
            model.create_<>_surface.
        - product
            a Species in the "world" or in a Region.
        - kd
            intrinsic reaction rate. Units: per second. (Rough order of 
            magnitude: 1e-2 /s to 1e2 /s).

    A surface unbinding reaction rule always has exactly one product.

    kd should be an *intrinsic* reaction rate. No analytical expression 
    is currently known to convert an overall reaction rate (koff) to 
    an intrinsic reaction rate (kd) for surface unbinding reaction rules.

    A surface unbinding reaction rule defines a Poissonian process.

    """
    assert reactant['structure'] == surface

def create_membrane_traversal_reaction_rule(reactant, surface1, 
                                            product, surface2, ka):
    """

    """
    # Implementation: also flip orientation of cylinder.
    pass
