#ifndef STRUCTURE_UTILS_HPP
#define STRUCTURE_UTILS_HPP

#include <string>
#include <typeinfo>
#include "linear_algebra.hpp"
#include "geometry.hpp"
#include "exceptions.hpp"
#include "Surface.hpp"
#include "Region.hpp"

template<typename Tsim_>
struct StructureUtils
{
    typedef Tsim_                                   simulator_type;
    typedef typename simulator_type::traits_type    traits_type;

    // some typedefs
    typedef typename traits_type::world_type::position_type             position_type;
    typedef typename traits_type::world_type::length_type               length_type;
    typedef typename traits_type::world_type::structure_name_type       structure_name_type;
    typedef typename traits_type::world_type::structure_id_type         structure_id_type;
    typedef typename traits_type::world_type::structure_type            structure_type;
    typedef typename traits_type::world_type::structure_type_id_type    structure_type_id_type;

    typedef typename simulator_type::surface_type               surface_type;
    typedef typename simulator_type::region_type                region_type;
    typedef typename simulator_type::sphere_type                sphere_type;
    typedef typename simulator_type::cylinder_type              cylinder_type;
    typedef typename simulator_type::disk_type                  disk_type;
    typedef typename simulator_type::box_type                   box_type;
    typedef typename simulator_type::plane_type                 plane_type;
    typedef typename simulator_type::spherical_surface_type     spherical_surface_type;
    typedef typename simulator_type::cylindrical_surface_type   cylindrical_surface_type;
    typedef typename simulator_type::disk_surface_type          disk_surface_type;
    typedef typename simulator_type::planar_surface_type        planar_surface_type;
    typedef typename simulator_type::cuboidal_region_type       cuboidal_region_type;
    typedef typename simulator_type::world_type::traits_type::rng_type rng_type;   
 
    static planar_surface_type* create_planar_surface(
            structure_type_id_type const& sid,          // This refers to the structure type of the planar surface
            structure_name_type   const& name,
            position_type const& corner,
            position_type const& unit_x,
            position_type const& unit_y,
            length_type const& lx,
            length_type const& ly,            
            structure_id_type const& parent_struct_id)
    {
        BOOST_ASSERT(is_cartesian_versor(unit_x));
        BOOST_ASSERT(is_cartesian_versor(unit_y));
        BOOST_ASSERT(is_cartesian_versor(cross_product(unit_x, unit_y)));

	// Note that when calling the function the origin is in the corner and the length
	// is the whole length of the plane, whereas for 'plane_type' the origin is in the
	// center of the plane (pos) and that the length is only half lengths 'half_lx' and 'half_ly'.
        const length_type half_lx(lx / 2);
        const length_type half_ly(ly / 2);

        const position_type pos(add(add(corner, multiply(unit_x, half_lx)),
                                                multiply(unit_y, half_ly)));
        const bool is_one_sided(true);

        return new planar_surface_type(name, sid, parent_struct_id,
                                       plane_type(pos, unit_x, unit_y,
                                                  half_lx, half_ly, is_one_sided));
    }
    
    static planar_surface_type* create_double_sided_planar_surface(
            structure_type_id_type const& sid,          // This refers to the structure type of the planar surface
            structure_name_type   const& name,
            position_type const& corner,
            position_type const& unit_x,
            position_type const& unit_y,
            length_type const& lx,
            length_type const& ly,            
            structure_id_type const& parent_struct_id)
    {
        BOOST_ASSERT(is_cartesian_versor(unit_x));
        BOOST_ASSERT(is_cartesian_versor(unit_y));
        BOOST_ASSERT(is_cartesian_versor(cross_product(unit_x, unit_y)));

        // Note that when calling the function the origin is in the corner and the length
        // is the whole length of the plane, whereas for 'plane_type' the origin is in the
        // center of the plane (pos) and that the length is only half lengths 'half_lx' and 'half_ly'.
        const length_type half_lx(lx / 2);
        const length_type half_ly(ly / 2);

        const position_type pos(add(add(corner, multiply(unit_x, half_lx)),
                                                multiply(unit_y, half_ly)));
        const bool is_one_sided(false);

        return new planar_surface_type(name, sid, parent_struct_id,
                                       plane_type(pos, unit_x, unit_y,
                                                  half_lx, half_ly, is_one_sided));
    }

    static spherical_surface_type* create_spherical_surface(
            structure_type_id_type const& sid,
            structure_name_type const& name,
            position_type const& pos,
            length_type const& radius,
            structure_id_type const& parent_struct_id)
    {
        return new spherical_surface_type(name, sid, parent_struct_id,
                                          sphere_type(pos, radius));
    }

    static cylindrical_surface_type* create_cylindrical_surface(
            structure_type_id_type const& sid,
            structure_name_type const& name,
            position_type const& corner,
            length_type const& radius,
            position_type const& unit_z,
            length_type const& length,
            structure_id_type const& parent_struct_id)
    {
        BOOST_ASSERT(is_cartesian_versor(unit_z));

        const length_type half_length(length / 2);
        const position_type pos(add(corner, multiply(unit_z, half_length)));

        return new cylindrical_surface_type(name, sid, parent_struct_id,
                                            cylinder_type(pos, radius, unit_z, half_length));
    }
    
    static disk_surface_type* create_disk_surface(
            structure_type_id_type const& sid,
            structure_name_type const& name,
            position_type const& center,
            length_type const& radius,
            position_type const& unit_z,
            structure_id_type const& parent_struct_id)
    {
        BOOST_ASSERT(is_cartesian_versor(unit_z));

        return new disk_surface_type(name, sid, parent_struct_id,
                disk_type(center, radius, unit_z));
    }

    static cuboidal_region_type* create_cuboidal_region(
            structure_type_id_type const& sid,
            structure_name_type const& name,
            position_type const& corner,
            boost::array<length_type, 3> const& extent,
            structure_id_type const& parent_struct_id)
    {
        const boost::array<length_type, 3> half_extent(divide(extent, 2));
        return new cuboidal_region_type(name, sid, parent_struct_id,
                                        box_type(add(corner, half_extent),
                                                 create_vector<position_type>(1, 0, 0),
                                                 create_vector<position_type>(0, 1, 0),
                                                 create_vector<position_type>(0, 0, 1),
                                                 half_extent));
    }

    static position_type random_vector(structure_type const& structure,
            length_type const& r, rng_type& rng)
    {
        return structure.random_vector(r, rng);
    }

    static position_type random_position(structure_type const& structure, rng_type& rng)
    {
        return structure.random_position(rng);
    }
/*
    static length_type minimal_distance_from_surface(surface_type const& surface, length_type const& radius)
    {
        return surface.minimal_distance(radius);
    }
*/
};

#endif /* STRUCTURE_UTILS_HPP */
