#ifndef SPECIES_INFO_HPP
#define SPECIES_INFO_HPP

#include <set>
#include <string>
#include <ostream>
#include "Defs.hpp"

template<typename Tid_, typename TD_, typename Tlen_, typename Tstructure_type_id_>
struct SpeciesInfo
// SpeciesInfo seems to extend the SpeciesType with extra information which is relevant for spatial simulations.
// A Species itself holds very little information.
// Note that the id is the same id as the SpeciesType and can be used to refer to the SpeciesType.
{
    // defining shorthands for used types
    typedef Tid_                identifier_type;        // This is the SpeciesTypeID
    typedef TD_                 D_type;
    typedef TD_                 v_type;
    typedef Tlen_               length_type;
    typedef Tstructure_type_id_ structure_type_id_type; // The structure type that the species lives on

    // constructors
    SpeciesInfo() {}

    SpeciesInfo(identifier_type const& id, structure_type_id_type const& s, D_type const& D = 0., 
                length_type const& r = 0., v_type const& v = 0.) 
        : id_(id), diffusion_coef_(D), drift_velocity_(v), radius_(r), structure_type_id_(s) {}



    // Get the id
    identifier_type const& id() const
    {
        return id_;
    }

    // Get the particle radius
    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    // Get the id of the structure type the species lives on
    structure_type_id_type const& structure_type_id() const
    {
        return structure_type_id_;
    }

    structure_type_id_type& structure_type_id()
    {
        return structure_type_id_;
    }

    // Get the diffusion constant
    D_type const& D() const
    {
        return diffusion_coef_;
    }

    D_type& D()
    {
        return diffusion_coef_;
    }

    // Get the drift
    v_type const& v() const
    {
        return drift_velocity_;
    }

    v_type& v()
    {
        return drift_velocity_;
    }

    // Check equality/inequality
    bool operator==(SpeciesInfo const& rhs) const
    {
        return id_ == rhs.id() && 
               diffusion_coef_ == rhs.D() &&
               drift_velocity_ == rhs.v() &&
               radius_ == rhs.radius() &&
               structure_type_id_ == rhs.structure_type_id();
    }

    bool operator!=(SpeciesInfo const& rhs) const
    {
        return !operator==(rhs);
    }

/////// Member variables
private:
    identifier_type         id_;
    D_type                  diffusion_coef_;
    v_type                  drift_velocity_;
    length_type             radius_;
    structure_type_id_type  structure_type_id_;      // The structure type that the particle lives on
};


//////// Inline functions
template<typename Tchar_, typename Ttraits_, typename Tid_, typename TD_, typename Tlen_, typename Tstructure_type_id_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& strm, const SpeciesInfo<Tid_, TD_, Tlen_, Tstructure_type_id_>& s)
// Allows for printing of the object
// -> provides a stream of 'attributename : attribute_value' pairs.
{
    strm << "SpeciesInfo(id=" << s.id() << ", D=" << s.D() << ", v=" << s.v() << ", radius=" << s.radius() << ", structure_type=" << s.structure_type_id() << ")";
    return strm;
}

#endif /* SPECIES_INFO_HPP */
