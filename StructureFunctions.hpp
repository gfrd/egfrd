#ifndef STRUCTFUNC_HPP
#define STRUCTFUNC_HPP

//#include "CuboidalRegion.hpp"
//#include "SphericalSurface.hpp"
//#include "CylindricalSurface.hpp"
//#include "DiskSurface.hpp"
//#include "PlanarSurface.hpp"
#include "linear_algebra.hpp"
#include "exceptions.hpp"
#include "Logger.hpp"

// Forward declaration of structures
template <typename Ttraits_>
class Structure;

template <typename Ttraits_>
class CuboidalRegion;

template <typename Ttraits_>
class SphericalSurface;

template <typename Ttraits_>
class CylindricalSurface;

template <typename Ttraits_>
class DiskSurface;

template <typename Ttraits_>
class PlanarSurface;


/******************************************************************************************/

// Define the structure functions and the actions performed in each case.
// There are two types: 1) yielding one new position and structure id
//                      2) yielding two new positions and structure ids
//                         (for single reactions with two products)

/******************************************************************************************/

/************************/
/*** ONE NEW POSITION ***/
/************************/

/******************************/
/* Coming from CuboidalRegion */
/******************************/
// CuboidalRegion -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;
    
    // Currently only species change and decay are supported
    if(origin_structure.id() == target_structure.id()){
     
          return std::make_pair( old_pos, origin_structure.id() );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (CuboidalRegion->CuboidalRegion).");
};

// CuboidalRegion -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (CuboidalRegion->Sphere).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CuboidalRegion -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id() );
    position_type       new_pos(   target_structure.project_point(old_pos).first );
    length_type         proj_dist( target_structure.project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition (CuboidalRegion->Cylinder).");
    
};

// CuboidalRegion -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id() );
    position_type       new_pos(   target_structure.project_point(old_pos).first );
    length_type         proj_dist( target_structure.project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition (CuboidalRegion->Disk).");
};

// CuboidalRegion -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id() );
    position_type       new_pos(   target_structure.project_point(old_pos).first );
    length_type         proj_dist( target_structure.project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition (CuboidalRegion->Plane).");
};

/********************************/
/* Coming from SphericalSurface */
/********************************/
// SphericalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->CuboidalRegion).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Cylinder).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Disk).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Plane).");
    
    return std::make_pair(position_type(), structure_id_type());
};

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
// CylindricalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;
    
    position_type displacement( origin_structure.surface_dissociation_vector(rng, offset, reaction_length) );
    position_type new_pos( add(old_pos, displacement) );    
    // TODO assert that new_pos is in target_structure
    
    return std::make_pair(new_pos, target_structure.id());
};

// CylindricalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently only species change and decay are supported
    if(origin_structure.id() == target_structure.id()){
     
          return std::make_pair( old_pos, origin_structure.id() );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (Cylinder->Cylinder).");
};

// CylindricalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Cylinder->Sphere).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CylindricalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id() );
    position_type       new_pos(   target_structure.project_point(old_pos).first );
    length_type         proj_dist( target_structure.project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    // Note that this function also correctly handles the pair reaction between a particle on a disk
    // and a particle on a cylinder. Since we assume that the product will end up on the disk, new_pos
    // (which is the projection of the center of mass of both particles on the disk) will be the
    // correct product position.
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition (Cylinder->Disk).");
};

// CylindricalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Cylinder->Plane).");
    
    return std::make_pair(position_type(), structure_id_type());
};

/***************************/
/* Coming from DiskSurface */
/***************************/
// DiskSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    position_type displacement( origin_structure.surface_dissociation_vector(rng, offset, reaction_length) );
    position_type new_pos( add(old_pos, displacement) );
    // TODO assert that new_pos is in target_structure
    
    return std::make_pair(new_pos, target_structure.id());
};

// DiskSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Sphere).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// DiskSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    position_type u( origin_structure.surface_dissociation_unit_vector(rng) );
    
    // When a particle is dissociation from a disk to the cylinder, the direction is randomized
    Real r(    rng.uniform(0.,1.) );
    Real r_rl( rng.uniform(0.,reaction_length) );
    Real r_displace( r < 0.5 ? -1.0*r_rl : +1.0*r_rl );
    
    position_type new_pos( add(old_pos, multiply(u, r_displace)) );
    // TODO assert that new_pos is in target_structure
    
    return std::make_pair(new_pos, target_structure.id());
};

// DiskSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently only species change and decay are supported
    if(origin_structure.id() == target_structure.id()){
     
          return std::make_pair( old_pos, origin_structure.id() );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (Disk->Disk).");    
};

// DiskSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Plane).");
    
    return std::make_pair(position_type(), structure_id_type());
};

/*****************************/
/* Coming from PlanarSurface */
/*****************************/
// PlanarSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    position_type displacement( origin_structure.surface_dissociation_vector(rng, offset, reaction_length) );
    position_type new_pos( add(old_pos, displacement) );
    // TODO assert that new_pos is in target_structure
    
    return std::make_pair(new_pos, target_structure.id());
};

// PlanarSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Sphere).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Cylinder).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Disk).");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& typical_length,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           &rng              )
{
    typedef typename Ttraits_::structure_id_type          structure_id_type;
    typedef typename Ttraits_::position_type              position_type;
    typedef typename Ttraits_::length_type                length_type;
    typedef std::pair<length_type, length_type>           length_pair_type;
    typedef std::pair<position_type, length_pair_type>    projected_type;    
    
    // Species change and decay: no structure change
    if( origin_structure.id() == target_structure.id() )
    {     
          return std::make_pair( old_pos, origin_structure.id() );
    }
    // Function is called with two different structure of the same type: Pair reaction "over the edge"
    else if( origin_structure.sid() == target_structure.sid() )
    {        
          // In this case old_pos = CoM of the two reactants, 
          // origin_structure = the origin structure of the first reactant,
          // target_structure = the origin structure of the second reactant.
          // We assume here that the problem has been correctly transformed
          // into the structure passed as target_structure to this function
          // before; however, to be sure, we check again and return the 
          // ID of the structure that the CoM (old_pos) lies in.
          // As a standard, this should be target_structure.
          
          // Calculate the projections of the CoM (old_pos) into both involved planes
          // This will also give the normal components of old_pos in the CS of the planes
          projected_type proj_o( origin_structure.project_point(old_pos) );
          projected_type proj_t( target_structure.project_point(old_pos) );
          
          // Check in which plane old_pos lies by comparing the normal components of the projection
          // Note: pair entries proj_X.second.first are the normal components
          // Then return the right position-structure_id pair; for safety we return the projected CoM
          // Note: apply_boundary will modify the product position accordingly if it lies outside of the (bounded) plane
          assert(typical_length > 0.0); // typical_length should contain the particle radius if this function was called correctly
          if( feq(proj_t.second.first, 0.0, typical_length) )
              // CoM is in target_structure
              return std::make_pair( proj_t.first, target_structure.id() );
          
          else if( feq(proj_o.second.first, 0.0, typical_length) )
              // CoM is in origin_structure
              return std::make_pair( proj_o.first, origin_structure.id() );
          else
              // CoM is in neither plane; something is terribly wrong...
              throw illegal_propagation_attempt("Center of mass is not in plane of either origin_structure or target_structure in pair reaction involving two PlanarSurfaces.");                    
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure type must be equal to target structure type for this type of structure transition (Plane->Plane).");        
};


/******************************************************************************************/


/*************************/
/*** TWO NEW POSITIONS ***/
/*************************/

/******************************/
/* Coming from CuboidalRegion */
/******************************/
// CuboidalRegion -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently we do not allow for pair forward/backward reactions between two different cubes
    if(origin_structure.id() == target_structure.id()){
     
        structure_id_type new_sid( origin_structure.id() );
        std::pair<position_type, position_type> new_positions( origin_structure.geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
        // geminate_dissociation_positions will produce two new positions close to old_pos taking into account
        // the type of origin_structure and the properties of the two product species
        // (the displacements from old_pos are weighted by the diffusion constants)
            
        return std::make_pair(      std::make_pair(new_positions.first,  new_sid),
                                    std::make_pair(new_positions.second, new_sid)    );          
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (CuboidalRegion->CuboidalRegion/CuboidalRegion).");        
};

// CuboidalRegion -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (CuboidalRegion->CuboidalRegion/Sphere).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CuboidalRegion -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (CuboidalRegion->CuboidalRegion/Cylinder).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
    
};

// CuboidalRegion -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (CuboidalRegion->CuboidalRegion/Disk).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CuboidalRegion -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (CuboidalRegion->CuboidalRegion/Plane).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/********************************/
/* Coming from SphericalSurface */
/********************************/
// SphericalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere/CuboidalRegion).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere/Sphere).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere/Cylinder).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere/Disk).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Sphere->Sphere/Plane).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
// CylindricalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    std::pair<position_type, position_type> new_positions( origin_structure.special_geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
    // special_geminate_dissociation_positions will produce two new positions close to old_pos taking into account
    // the types of origin_structure and target_structure and the properties of the two product species
    // (the displacements from old_pos are weighted by the diffusion constants);
    // the first pair entry is the particle staying on the surface, the second entry the particle
    // that goes into the bulk
            
    return std::make_pair(      std::make_pair(new_positions.first,  origin_structure.id()),
                                std::make_pair(new_positions.second, target_structure.id())    );  
};

// CylindricalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently we do not allow for pair forward/backward reactions between two different cylinders
    if(origin_structure.id() == target_structure.id()){
     
        structure_id_type new_sid( origin_structure.id() );
        std::pair<position_type, position_type> new_positions( origin_structure.geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
        // geminate_dissociation_positions will produce two new positions close to old_pos taking into account
        // the type of origin_structure and the properties of the two product species
        // (the displacements from old_pos are weighted by the diffusion constants)
            
        return std::make_pair(      std::make_pair(new_positions.first,  new_sid),
                                    std::make_pair(new_positions.second, new_sid)    );          
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (Cylinder->Cylinder/Cylinder).");        
};

// CylindricalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Cylinder->Cylinder/Sphere).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Cylinder->Cylinder/Disk).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Cylinder->Cylinder/Plane).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/***************************/
/* Coming from DiskSurface */
/***************************/
// DiskSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    std::pair<position_type, position_type> new_positions( origin_structure.special_geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
    // special_geminate_dissociation_positions will produce two new positions close to old_pos taking into account
    // the types of origin_structure and target_structure and the properties of the two product species
    // (the displacements from old_pos are weighted by the diffusion constants);
    // the first pair entry is the particle staying on the surface, the second entry the particle
    // that goes into the bulk
            
    return std::make_pair(      std::make_pair(new_positions.first,  origin_structure.id()),
                                std::make_pair(new_positions.second, target_structure.id())    );
};

// DiskSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Disk/Sphere).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // TODO TODO TODO
    // IMPLEMENT THE UNBINDING FROM THE DISK ACTING AS A SINK!
    
    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Disk/Cylinder).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Disk/Disk).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Disk->Disk/Plane).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/*****************************/
/* Coming from PlanarSurface */
/*****************************/
// PlanarSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;
    
    std::pair<position_type, position_type> new_positions;
    length_type dist_to_edge_1;
    length_type dist_to_edge_2;
    
    bool positions_legal(false);
    
    do{
        // Produce two new positions
        new_positions = origin_structure.special_geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length);
          // special_geminate_dissociation_positions will produce two new positions close to old_pos taking into account
          // the types of origin_structure and target_structure and the properties of the two product species
          // (the displacements from old_pos are weighted by the diffusion constants);
          // the first pair entry is the particle staying on the surface, the second entry the particle
          // that goes into the bulk
        
        // Check whether the newly produced positions are above the origin plane
        dist_to_edge_1 = origin_structure.project_point(new_positions.first).second.second;
        dist_to_edge_2 = origin_structure.project_point(new_positions.second).second.second;
        
        positions_legal = (dist_to_edge_1 < 0.0) and (dist_to_edge_2 < 0.0);    
    }
    while( not(positions_legal) );
            
    return std::make_pair(      std::make_pair(new_positions.first,  origin_structure.id()),
                                std::make_pair(new_positions.second, target_structure.id())    );    
};

// PlanarSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Plane/Sphere).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Plane/Cylinder).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported (Plane->Plane/Disk).");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s_orig,
                       typename Ttraits_::species_type        const& s_targ,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            &rng                )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // As a default we produce two new positions on the same plane; in principle the particles can end up
    // on different planes, but this should be treated afterwards via apply_boundary.
    // TODO: Move apply_boundary into structure functions?
    if(origin_structure.id() == target_structure.id()){
     
        structure_id_type new_sid( origin_structure.id() );
        std::pair<position_type, position_type> new_positions( origin_structure.geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
        // geminate_dissociation_positions will produce two new positions close to old_pos taking into account
        // the type of origin_structure and the properties of the two product species
        // (the displacements from old_pos are weighted by the diffusion constants)
            
        return std::make_pair(      std::make_pair(new_positions.first,  new_sid),
                                    std::make_pair(new_positions.second, new_sid)    );          
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (Plane->Plane/Plane).");        
};

/******************************************************************************************************/


#endif /* STRUCTUREFUNCTIONS_HPP */

