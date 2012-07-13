#ifndef SPHERICAL_SURFACE_HPP
#define SPHERICAL_SURFACE_HPP

#include "Surface.hpp"
#include "Sphere.hpp"
#include "freeFunctions.hpp"
#include "StructureFunctions.hpp"

template <typename Tobj_, typename Tid_, typename Ttraits_>
class StructureContainer;

template<typename Ttraits_>
class SphericalSurface
    : public BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::length_type> >
{
public:
    typedef BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::length_type> > base_type;
    typedef Ttraits_ traits_type;
    typedef typename base_type::structure_name_type     structure_name_type;        // THis is just the name of the structure
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::shape_type              shape_type;
    typedef typename base_type::rng_type                rng_type;
    typedef typename base_type::position_type           position_type;
    typedef typename base_type::length_type             length_type;
    typedef typename base_type::side_enum_type          side_enum_type;
    typedef typename traits_type::species_type          species_type;
    typedef typename traits_type::structure_type        structure_type;

    typedef StructureContainer<typename traits_type::structure_type, structure_id_type, traits_type>    structure_container_type;

    typedef std::pair<position_type, position_type>                               position_pair_type;
    typedef std::pair<position_type, structure_id_type>                           position_structid_pair_type;
    typedef std::pair<position_structid_pair_type, position_structid_pair_type>   position_structid_pair_pair_type;


    /*** Simple structure-specific sampling functions ***/
    virtual position_type random_position(rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    /*** New BD scheme functions ***/
    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const
    {
        return Real(); //TODO
    }
    
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const
    {
        return Real(); //TODO
    }

    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const
    {
        return Real(); //TODO
    }
    
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const
    {
        return Real(); //TODO
    }

    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        return position_type(); //TODO  
    }
        
    virtual position_type surface_dissociation_unit_vector( rng_type& rng ) const
    {
        return position_type(); //TODO  
    }
    
    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, length_type const& rl ) const
    {
        return position_pair_type(); //TODO
    }
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, position_type const& op_surf, length_type const& rl ) const
    {
        return position_pair_type(); //TODO
    }
    
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const
    {
        return base_type::distance(new_pos);
    }
/*
    virtual length_type minimal_distance(length_type const& radius) const
    {
        return 0.; // TODO
    }
*/

    /*** Boundary condition handling ***/
    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id,
                                                       structure_container_type const& structure_container) const
    {
        // The apply_boundary for a spherical surface is trivial (since there is no boundary!)
        return pos_struct_id;
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_container_type const& structure_container) const
    {
        return pos_struct_id;       // Two spherical surface cannot be connected (there is no boundary!)
    }


    // *** Dynamic dispatch for the structure functions *** //
    // *** 1 *** - One new position
    // This requires a double dynamic dispatch.
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& target_structure, position_type const& position,
                                                         length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {
        return target_structure.get_pos_sid_pair_helper(*this, position, offset, reaction_length, rng);
    }
//     // the associated helper function
//     template <typename Tstruct_>
//     position_structid_pair_type get_pos_sid_pair_helper(Tstruct_ const& origin_structure, position_type const& position,
//                                                         length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//         // Dispatch to function with well-defined typing
//         return ::get_pos_sid_pair(origin_structure, *this, position, offset, reaction_length, rng);
//     }
    
    // *** 2 *** - Two new positions
    // Same principle as above, but different return type
    virtual position_structid_pair_pair_type get_pos_sid_pair_pair(structure_type const& target_structure, position_type const& position,
                                                                   species_type const& s1, species_type const& s2, length_type const& reaction_length, rng_type const& rng) const
    {
        return target_structure.get_pos_sid_pair_pair_helper(*this, position, s1, s2, reaction_length, rng);
    }
//     // the associated helper function
//     template <typename Tstruct_>
//     position_structid_pair_pair_type get_pos_sid_pair_pair_helper(Tstruct_ const& origin_structure, position_type const& position,
//                                                                   species_type const& s1, species_type const& s2, length_type const& reaction_length, rng_type const& rng) const
//     {
//         // Dispatch to function with well-defined typing
//         return ::get_pos_sid_pair_pair(origin_structure, *this, position, s1, s2, reaction_length, rng);
//     }
    
    // *** 3 *** - Pair reactions => two origin structures => triple switchbox
    // Overloading get_pos_sid_pair with signature (origin_structure2, target_structure_type_id, ...) and
    // introducing another double dynamic dispatch.
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type_id_type const& target_sid,
                                                         position_type const& CoM, length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {
        return origin_structure2.get_pos_sid_pair_helper_two_origins(*this, target_sid, CoM, offset, reaction_length, rng);
    }
//     // two associated helper function
//     template <typename Tstruct_>
//     position_structid_pair_type get_pos_sid_pair_helper_two_origins(Tstruct_ const& origin_structure1, structure_type_id_type const& target_sid,
//                                                                     position_type const& CoM, length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//           // The types of both origin structures now are determined.
//           // As a next step, determine which one is the target structure and call the right
//           // structure function with the other one as origin structure.
//           
//           if( this->structure_type_id() != origin_structure1->structure_type_id() )
//           // if the two pair reactants come from different types of structure
//           {
// 
//               if ( (origin_structure1->structure_id() == this->id) && (origin_structure1->structure_type_id() == target_sid) )
//               // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
//               {
//                   // Call transition function with *this as origin_structure and origin_structure1 as target_structure
//                   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
//               }
//               else if ( (this->structure_id() == origin_structure1->id) && (this->structure_type_id() == target_sid) )
//               // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
//               {
//                   // Call transition function with origin_structure1 as origin_structure and *this as target_structure
//                   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
//               }
//               else
//               {
//                   throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
//                   // TODO In principle this constraint could be dropped; the target structure should be determined
//                   // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
//               }  
//         }
//         else
//         // the reactants live on the same structure type, i.e. the product will also end up on 
//         {
//             // Call transition function with origin_structure = target_structure = *this
//             return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
//             // Note that apply_boundary should place the product particle on the right one of the two structures
//             // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
//             // structure_type_id's are the same!).
//         }
//         
//     }
    
//     // *** 4 *** - Generalized functions for pair reactions with two origin structures and one target structure
//     // NOTE: This is yet unused, but possibly useful in the future.
//     // Overloading get_pos_sid_pair again with signature (origin_structure2, target_structure, ...) and introducing
//     // a triple dynamic dispatch.
//     virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type const& target_structure, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//         return origin_structure2.get_pos_sid_pair_helper1(*this, target_structure, position, offset, reaction_length, rng);
//     }
//     // two associated helper functions this time
//     template <typename Tstruct1_>
//     position_structid_pair_type get_pos_sid_pair_helper1(Tstruct1_ const& origin_structure1, structure_type const& target_structure, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//         return target_structure.get_pos_sid_pair_helper2(origin_structure1, *this, position, offset, reaction_length, rng);
//     }
//     template <typename Tstruct1_, typename Tstruct2_>
//     position_structid_pair_type get_pos_sid_pair_helper2(Tstruct1_ const& origin_structure1, Tstruct2_ const& origin_structure2, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//         if(origin_structure1.id == *this->id)
//             // Dispatch to function with well-defined typing
//             return ::get_pos_sid_pair(origin_structure2, *this, position, offset, reaction_length, rng);
//         
//         else if(origin_structure2.id == *this->id)
//             // Dispatch to function with well-defined typing
//             return ::get_pos_sid_pair(origin_structure1, *this, position, offset, reaction_length, rng);
//         
//         else
//             throw propagation_error("Target structure must be one of the origin structures for pair reaction.");
//     }
    
    
    /*** Formerly used functions of the Morelli scheme ***/
    // DEPRECATED
    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {    
        return length_type(); // TODO
    }
    // DEPRECATED
    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {    
        return Real(); //TODO
    }
    // DEPRECATED
    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return position_type(); //TODO
    }
    
    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    // The Constructor
    SphericalSurface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id, shape) {}
};

#endif /* SPHERICAL_SURFACE_HPP */
