#ifndef STRUCTFUNC_HPP
#define STRUCTFUNC_HPP


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

// Now declare the structure functions.
// There are two types: 1) yielding one new position and structure id
//                      2) yielding two new positions and structure ids
//                         (for single reactions with two products)

/******************************************************************************************/

/************************/
/*** ONE NEW POSITION ***/
/************************/

/******************************/
/* Some Defaults declarations */
/******************************/
// These are actually never called but have to be there for type consistency reasons
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  Structure<Ttraits_>                   const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  Structure<Ttraits_>                   const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  Structure<Ttraits_>                   const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  Structure<Ttraits_>                   const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>              const& origin_structure,
                  Structure<Ttraits_>                   const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( Structure<Ttraits_>                   const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( Structure<Ttraits_>                   const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( Structure<Ttraits_>                   const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( Structure<Ttraits_>                   const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( Structure<Ttraits_>                   const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

/******************************/
/* Coming from CuboidalRegion */
/******************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

/********************************/
/* Coming from SphericalSurface */
/********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

/***************************/
/* Coming from DiskSurface */
/***************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
/*****************************/
/* Coming from PlanarSurface */
/*****************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos,
                  typename Ttraits_::length_type        const& offset,
                  typename Ttraits_::length_type        const& reaction_length,
                  typename Ttraits_::rng_type           const& rng              );

                  
/******************************************************************************************/


/*************************/
/*** TWO NEW POSITIONS ***/
/*************************/

/******************************/
/* Some Defaults declarations */
/******************************/
// These are actually never called but have to be there for type consistency reasons
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       Structure<Ttraits_>                    const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       Structure<Ttraits_>                    const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       Structure<Ttraits_>                    const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       Structure<Ttraits_>                    const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       Structure<Ttraits_>                    const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 )
{
    throw illegal_propagation_attempt("Calling structure function without type specification!");
}

/******************************/
/* Coming from CuboidalRegion */
/******************************/
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>               const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

/********************************/
/* Coming from SphericalSurface */
/********************************/
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>             const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>           const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure,
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

/***************************/
/* Coming from DiskSurface */
/***************************/
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                  const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
/*****************************/
/* Coming from PlanarSurface */
/*****************************/
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       CuboidalRegion<Ttraits_>               const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       SphericalSurface<Ttraits_>             const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       CylindricalSurface<Ttraits_>           const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       DiskSurface<Ttraits_>                  const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );
                  
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>                const& origin_structure,
                       PlanarSurface<Ttraits_>                const& target_structure, 
                       typename Ttraits_::position_type       const& old_pos,
                       typename Ttraits_::species_type        const& s1,
                       typename Ttraits_::species_type        const& s2,
                       typename Ttraits_::length_type         const& reaction_length,
                       typename Ttraits_::rng_type            const& rng                 );

                       
/******************************************************************************************************/

                  
#endif /* STRUCTUREFUNCTIONS_HPP */

