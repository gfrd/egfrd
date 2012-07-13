#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <sstream>
#include "Vector3.hpp"
#include "exceptions.hpp"
#include "freeFunctions.hpp"
#include "SpeciesTypeID.hpp"
#include "StructureFunctions.hpp"

// Forward declarations
template <typename Tobj_, typename Tid, typename Ttraits_>
class StructureContainer;

template<typename T_>
class CuboidalRegion;

template<typename T_>
class CylindricalSurface;

template<typename T_>
class SphericalSurface;

template<typename T_>
class DiskSurface;

template<typename T_>
class PlanarSurface;



template<typename Ttraits_>
class Structure
{
public:
    typedef Ttraits_ traits_type;
    // shorthands for types that we use
    typedef typename traits_type::rng_type                                        rng_type;
    typedef typename traits_type::structure_name_type                             structure_name_type;
    typedef typename traits_type::structure_id_type                               structure_id_type;
    typedef typename traits_type::structure_type                                  structure_type;
    typedef typename traits_type::length_type                                     length_type;
    typedef typename traits_type::position_type                                   position_type;
    typedef typename traits_type::species_type                                    species_type;
    typedef typename traits_type::structure_type_id_type                          structure_type_id_type;
    typedef std::pair<length_type, length_type>                                   components_pair_type;
    typedef std::pair<position_type, components_pair_type>                        projected_type;
    typedef std::pair<position_type, position_type>                               position_pair_type;
    typedef std::pair<position_type, bool>                                        position_flag_pair_type;
    typedef std::pair<position_type, structure_id_type>                           position_structid_pair_type;
    typedef std::pair<position_structid_pair_type, position_structid_pair_type>   position_structid_pair_pair_type;
    typedef StructureContainer<typename traits_type::structure_type, structure_id_type, traits_type>   structure_container_type;

public:
    virtual ~Structure() {}

    structure_id_type const& id() const
    {
        if (!id_)
        {
            throw illegal_state("ID for structure not defined");
        }
        return id_;
    }

    void set_id(structure_id_type const& id)
    {
        id_ = id;
    }

    structure_name_type const& name()
    {
        return name_;
    }

    // Get the StructureType of the structure
    structure_type_id_type const& sid() const
    {
        if (!sid_)
        {
            throw illegal_state("not bound to StructureType");
        }
        return sid_;
    }

    structure_type_id_type& sid()
    {
        return sid_;
    }

    structure_id_type const& structure_id() const
    {
        return parent_struct_id_;
    }

    virtual bool operator==(Structure const& rhs) const
    {
        return id_ == rhs.id() && sid_ == rhs.sid();
    }

    bool operator!=(Structure const& rhs) const
    {
        return !operator==(rhs);
    }

    virtual position_type random_position(rng_type& rng) const = 0;

    virtual position_type random_vector(length_type const& r, rng_type& rng) const = 0;


    // Methods used in the 'old' BDPropagator // DEPRECATED
    virtual position_type dissociation_vector(rng_type& rng, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const = 0;
    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const = 0;
    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, Real const& D0, Real const& D1, Real const& v0, Real const& v1) const = 0;

    // Methods used in the 'new' BDPropagator
    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const = 0;
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const = 0;

    // TODO this are just functions->move somewhere else
    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const = 0;
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const = 0;
    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const = 0;
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const = 0;     // This does contain a surface dependent component.
    
    // Methods used to calculate dissociation positions
    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const = 0;
    virtual position_type surface_dissociation_unit_vector( rng_type& rng ) const = 0;
    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, length_type const& rl ) const = 0;
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, position_type const& op_surf, length_type const& rl ) const = 0;
    
    // General method for getting some measures/info
    virtual projected_type project_point(position_type const& pos) const = 0;
    virtual projected_type project_point_on_surface(position_type const& pos) const = 0;
    virtual length_type distance(position_type const& pos) const = 0;
    virtual position_type const& position() const = 0;

    // Methods used for edge crossing (only for the planes so far)
    virtual position_flag_pair_type deflect(position_type const& pos0, position_type const& displacement) const = 0;
//    virtual position_type deflect_back(position_type const& pos, position_type const& u_z) const = 0;

    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id,
                                                       structure_container_type const& structure_container) const = 0;
    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_container_type const& structure_container) const = 0;
    // Structure functions dynamic dispatch
    // 1 - Producing one new position
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& target_structure, position_type const& position,
                                                         length_type const& offset, length_type const& rl, rng_type const& rng) const = 0;                                                         
    position_structid_pair_type get_pos_sid_pair_helper(CuboidalRegion<traits_type> const& origin_structure, position_type const& position,
                                                        length_type const& offset, length_type const& rl, rng_type const& rng) const
                                                        { return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng); };    
    position_structid_pair_type get_pos_sid_pair_helper(SphericalSurface<traits_type> const& origin_structure, position_type const& position,
                                                        length_type const& offset, length_type const& rl, rng_type const& rng) const
                                                        { return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng); };    
    position_structid_pair_type get_pos_sid_pair_helper(CylindricalSurface<traits_type> const& origin_structure, position_type const& position,
                                                        length_type const& offset, length_type const& rl, rng_type const& rng) const
                                                        { return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng); };    
    position_structid_pair_type get_pos_sid_pair_helper(DiskSurface<traits_type> const& origin_structure, position_type const& position,
                                                        length_type const& offset, length_type const& rl, rng_type const& rng) const
                                                        { return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng); };    
    position_structid_pair_type get_pos_sid_pair_helper(PlanarSurface<traits_type> const& origin_structure, position_type const& position,
                                                        length_type const& offset, length_type const& rl, rng_type const& rng) const
                                                        { return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng); };
    
    // 2 - Producing two new positions
    virtual position_structid_pair_pair_type get_pos_sid_pair_pair(structure_type const& target_structure, position_type const& position,
                                                                   species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const = 0;
        
    position_structid_pair_pair_type get_pos_sid_pair_pair_helper(CuboidalRegion<traits_type> const& origin_structure, position_type const& position,
                                                                  species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const
                                                                  { return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng); };
    position_structid_pair_pair_type get_pos_sid_pair_pair_helper(SphericalSurface<traits_type> const& origin_structure, position_type const& position,
                                                                  species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const
                                                                  { return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng); };
    position_structid_pair_pair_type get_pos_sid_pair_pair_helper(CylindricalSurface<traits_type> const& origin_structure, position_type const& position,
                                                                  species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const
                                                                  { return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng); };
    position_structid_pair_pair_type get_pos_sid_pair_pair_helper(DiskSurface<traits_type> const& origin_structure, position_type const& position,
                                                                  species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const
                                                                  { return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng); };
    position_structid_pair_pair_type get_pos_sid_pair_pair_helper(PlanarSurface<traits_type> const& origin_structure, position_type const& position,
                                                                  species_type const& s_orig, species_type const& s_targ, length_type const& rl, rng_type const& rng) const
                                                                  { return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng); };
        
    // 3 - Pair reactions => two origin structures
    // Overloading method call structure.get_pos_sid_pair
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type_id_type const& target_sid, position_type const& CoM,
                                                         length_type const& offset, length_type const& reaction_length, rng_type const& rng) const = 0;
    
    // TODO This is just so ugly, outsource this somehow!
    position_structid_pair_type get_pos_sid_pair_helper_two_origins(CuboidalRegion<traits_type> const& origin_structure1, structure_type_id_type const& target_sid, position_type const& CoM,
                                                                    length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {          
          // The types of both origin structures now are determined.
          // As a next step, determine which one is the target structure and call the right
          // structure function with the other one as origin structure.          
          
          structure_id_type      os1_id( origin_structure1.id() );
          structure_id_type      os1_parent_id( origin_structure1.structure_id() );
          structure_type_id_type os1_sid( origin_structure1.sid() );
          
          structure_id_type      this_id( this->id() );
          structure_id_type      this_parent_id( this->structure_id() );
          structure_type_id_type this_sid( this->sid() );
          
          if( this->sid() != os1_sid )
          // if the two pair reactants come from different types of structure
          {

              if ( (os1_parent_id == this_id) && (os1_sid == target_sid) )
              // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
              {
                  // Call transition function with *this as origin_structure and origin_structure1 as target_structure
                  return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
              }
              else if ( (this_parent_id == os1_id) && (this_sid == target_sid) )
              // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
              {
                  // Call transition function with origin_structure1 as origin_structure and *this as target_structure
                  return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
              }
              else
              {
                  throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
                  // TODO In principle this constraint could be dropped; the target structure should be determined
                  // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
              }  
        }
        else
        // the reactants live on the same structure type, i.e. the product will also end up on 
        {
            // Call transition function with origin_structure = target_structure = *this
            return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
            // Note that apply_boundary should place the product particle on the right one of the two structures
            // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
            // structure_type_id's are the same!).
        }
    } 
    position_structid_pair_type get_pos_sid_pair_helper_two_origins(SphericalSurface<traits_type> const& origin_structure1, structure_type_id_type const& target_sid, position_type const& CoM,
                                                                    length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {          
          // The types of both origin structures now are determined.
          // As a next step, determine which one is the target structure and call the right
          // structure function with the other one as origin structure.          
          
          structure_id_type      os1_id( origin_structure1.id() );
          structure_id_type      os1_parent_id( origin_structure1.structure_id() );
          structure_type_id_type os1_sid( origin_structure1.sid() );
          
          structure_id_type      this_id( this->id() );
          structure_id_type      this_parent_id( this->structure_id() );
          structure_type_id_type this_sid( this->sid() );
          
          if( this->sid() != os1_sid )
          // if the two pair reactants come from different types of structure
          {

              if ( (os1_parent_id == this_id) && (os1_sid == target_sid) )
              // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
              {
                  // Call transition function with *this as origin_structure and origin_structure1 as target_structure
                  return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
              }
              else if ( (this_parent_id == os1_id) && (this_sid == target_sid) )
              // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
              {
                  // Call transition function with origin_structure1 as origin_structure and *this as target_structure
                  return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
              }
              else
              {
                  throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
                  // TODO In principle this constraint could be dropped; the target structure should be determined
                  // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
              }  
        }
        else
        // the reactants live on the same structure type, i.e. the product will also end up on 
        {
            // Call transition function with origin_structure = target_structure = *this
            return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
            // Note that apply_boundary should place the product particle on the right one of the two structures
            // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
            // structure_type_id's are the same!).
        }
    } 
    position_structid_pair_type get_pos_sid_pair_helper_two_origins(CylindricalSurface<traits_type> const& origin_structure1, structure_type_id_type const& target_sid, position_type const& CoM,
                                                                    length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {          
          // The types of both origin structures now are determined.
          // As a next step, determine which one is the target structure and call the right
          // structure function with the other one as origin structure.          
          
          structure_id_type      os1_id( origin_structure1.id() );
          structure_id_type      os1_parent_id( origin_structure1.structure_id() );
          structure_type_id_type os1_sid( origin_structure1.sid() );
          
          structure_id_type      this_id( this->id() );
          structure_id_type      this_parent_id( this->structure_id() );
          structure_type_id_type this_sid( this->sid() );
          
          if( this->sid() != os1_sid )
          // if the two pair reactants come from different types of structure
          {

              if ( (os1_parent_id == this_id) && (os1_sid == target_sid) )
              // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
              {
                  // Call transition function with *this as origin_structure and origin_structure1 as target_structure
                  return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
              }
              else if ( (this_parent_id == os1_id) && (this_sid == target_sid) )
              // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
              {
                  // Call transition function with origin_structure1 as origin_structure and *this as target_structure
                  return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
              }
              else
              {
                  throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
                  // TODO In principle this constraint could be dropped; the target structure should be determined
                  // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
              }  
        }
        else
        // the reactants live on the same structure type, i.e. the product will also end up on 
        {
            // Call transition function with origin_structure = target_structure = *this
            return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
            // Note that apply_boundary should place the product particle on the right one of the two structures
            // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
            // structure_type_id's are the same!).
        }
    } 
    position_structid_pair_type get_pos_sid_pair_helper_two_origins(DiskSurface<traits_type> const& origin_structure1, structure_type_id_type const& target_sid, position_type const& CoM,
                                                                    length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {          
          // The types of both origin structures now are determined.
          // As a next step, determine which one is the target structure and call the right
          // structure function with the other one as origin structure.          
          
          structure_id_type      os1_id( origin_structure1.id() );
          structure_id_type      os1_parent_id( origin_structure1.structure_id() );
          structure_type_id_type os1_sid( origin_structure1.sid() );
          
          structure_id_type      this_id( this->id() );
          structure_id_type      this_parent_id( this->structure_id() );
          structure_type_id_type this_sid( this->sid() );
          
          if( this->sid() != os1_sid )
          // if the two pair reactants come from different types of structure
          {

              if ( (os1_parent_id == this_id) && (os1_sid == target_sid) )
              // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
              {
                  // Call transition function with *this as origin_structure and origin_structure1 as target_structure
                  return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
              }
              else if ( (this_parent_id == os1_id) && (this_sid == target_sid) )
              // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
              {
                  // Call transition function with origin_structure1 as origin_structure and *this as target_structure
                  return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
              }
              else
              {
                  throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
                  // TODO In principle this constraint could be dropped; the target structure should be determined
                  // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
              }  
        }
        else
        // the reactants live on the same structure type, i.e. the product will also end up on 
        {
            // Call transition function with origin_structure = target_structure = *this
            return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
            // Note that apply_boundary should place the product particle on the right one of the two structures
            // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
            // structure_type_id's are the same!).
        }
    } 
    position_structid_pair_type get_pos_sid_pair_helper_two_origins(PlanarSurface<traits_type> const& origin_structure1, structure_type_id_type const& target_sid, position_type const& CoM,
                                                                    length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {          
          // The types of both origin structures now are determined.
          // As a next step, determine which one is the target structure and call the right
          // structure function with the other one as origin structure.          
          
          structure_id_type      os1_id( origin_structure1.id() );
          structure_id_type      os1_parent_id( origin_structure1.structure_id() );
          structure_type_id_type os1_sid( origin_structure1.sid() );
          
          structure_id_type      this_id( this->id() );
          structure_id_type      this_parent_id( this->structure_id() );
          structure_type_id_type this_sid( this->sid() );
          
          if( this->sid() != os1_sid )
          // if the two pair reactants come from different types of structure
          {

              if ( (os1_parent_id == this_id) && (os1_sid == target_sid) )
              // *this is the parent of origin_structure1, i.e. target_structure is origin_structure1
              {
                  // Call transition function with *this as origin_structure and origin_structure1 as target_structure
                  return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);
              }
              else if ( (this_parent_id == os1_id) && (this_sid == target_sid) )
              // origin_structure1 is the parent of *this, i.e. target_structure is *this (= origin_structure2)
              {
                  // Call transition function with origin_structure1 as origin_structure and *this as target_structure
                  return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
              }
              else
              {
                  throw propagation_error("Particles can be at most one hierarchical level apart for a pair reaction.");
                  // TODO In principle this constraint could be dropped; the target structure should be determined
                  // already by comparing structure_type_id's with target_sid (which is fed in from the reaction rules).
              }  
        }
        else
        // the reactants live on the same structure type, i.e. the product will also end up on 
        {
            // Call transition function with origin_structure = target_structure = *this
            return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);
            // Note that apply_boundary should place the product particle on the right one of the two structures
            // afterwards in case that *this->id != origin_structure1->id (right now we postulate only that the
            // structure_type_id's are the same!).
        }
    } 
    
    
    // TODO Outsource the ugly stuff above at least partly into this!
    bool is_parent_of_or_same_as(structure_type s)
    {    
        return true;
    }
                                                                    
    // 4 - Generalized functions for pair reactions => two origin structures and one target_structure
    // This introduces a triple dynamic dispatch, overloading method call structure.get_pos_sid_pair once more.
    // NOTE: As yet these methods are unused but might prove useful in the future.
//     virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type const& target_structure, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const = 0;
//     template <typename Tstruct1_>
//     position_structid_pair_type get_pos_sid_pair_helper1(Tstruct1_ const& origin_structure1, structure_type const& target_structure, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const;
//     template <typename Tstruct1_, typename Tstruct2_>
//     position_structid_pair_type get_pos_sid_pair_helper2(Tstruct1_ const& origin_structure1, Tstruct2_ const& origin_structure2, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const;




    virtual std::size_t hash() const
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<structure_name_type>()(name_) ^
               hash<structure_type_id_type>()(sid_);
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Structure(" << id() << ", " << sid() << ")";
        return out.str();
    }

    // Constructor
    Structure(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id)
        : name_(name), sid_(sid), parent_struct_id_(parent_struct_id) {}

////// Member variables
protected:
    structure_name_type     name_;                  // just the name
    structure_type_id_type  sid_;                   // id of the structure_type of the structure
    structure_id_type       parent_struct_id_;

    structure_id_type       id_;        // id of the structure (filled in later)
};


//////// Inline functions
template<typename Tstrm, typename Ttraits, typename T_traits>
inline std::basic_ostream<Tstrm, Ttraits>& operator<<(std::basic_ostream<Tstrm, Ttraits>& strm, const Structure<T_traits>& v)
{
    strm << v.as_string(); 
    return strm;
}

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename Ttraits>
struct hash<Structure<Ttraits> >
{
    typedef Structure<Ttraits> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return val.hash();
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* STRUCTURE_HPP */
