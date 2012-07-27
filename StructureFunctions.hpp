#ifndef STRUCTFUNC_HPP
#define STRUCTFUNC_HPP

//#include "CuboidalRegion.hpp"
//#include "SphericalSurface.hpp"
//#include "CylindricalSurface.hpp"
//#include "DiskSurface.hpp"
//#include "PlanarSurface.hpp"
#include "linear_algebra.hpp"
#include "exceptions.hpp"

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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
    
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
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
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
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");        
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");        
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");        
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    
    std::pair<position_type, position_type> new_positions( origin_structure.special_geminate_dissociation_positions(rng, s_orig, s_targ, old_pos, reaction_length) );
    // special_geminate_dissociation_positions will produce two new positions close to old_pos taking into account
    // the types of origin_structure and target_structure and the properties of the two product species
    // (the displacements from old_pos are weighted by the diffusion constants);
    // the first pair entry is the particle staying on the surface, the second entry the particle
    // that goes into the bulk
            
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
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
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");        
};

/******************************************************************************************************/


#endif /* STRUCTUREFUNCTIONS_HPP */

