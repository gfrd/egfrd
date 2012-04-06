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
#include "freeFunctions.hpp"
#include "SpeciesTypeID.hpp"

template<typename Ttraits_>
class Structure
{
public:
    typedef Ttraits_ traits_type;
    // shorthands for types that we use
    typedef typename traits_type::rng_type                  rng_type;
    typedef typename traits_type::structure_name_type       identifier_type;
    typedef typename traits_type::structure_id_type         structure_id_type;
    typedef typename traits_type::length_type               length_type;
    typedef typename traits_type::position_type             position_type;
    typedef typename traits_type::base_type::species_type   species_type;
    typedef typename traits_type::structure_type_id_type    structure_type_id_type;
    typedef std::pair<position_type, length_type>           projected_type;
    typedef std::pair<position_type, position_type>         position_pair_type;
    typedef std::pair<position_type, bool>                  position_flag_pair_type;

public:
    virtual ~Structure() {}

    identifier_type const& name()
    {
        return id_;
    }

    structure_id_type const& real_id() const
    {
        if (!real_id_)
        {
            throw illegal_state("ID for structure not defined");
        }
        return real_id_;
    }

    void set_id(structure_id_type const& id)
    {
        real_id_ = id;
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

    virtual bool operator==(Structure const& rhs) const
    {
        return real_id_ == rhs.real_id() && sid_ == rhs.sid();
    }

    bool operator!=(Structure const& rhs) const
    {
        return !operator==(rhs);
    }

    virtual position_type random_position(rng_type& rng) const = 0;

    virtual position_type random_vector(length_type const& r, rng_type& rng) const = 0;

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const = 0;

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const = 0;

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, Real const& D0, Real const& D1, Real const& v0, Real const& v1) const = 0;

    virtual position_type dissociation_vector(rng_type& rng, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const = 0;
    
    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const = 0;
    
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const = 0;

    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const = 0;
    
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const = 0;
    
    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const = 0;
    
    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, length_type const& rl ) const = 0;
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, position_type const& op_surf, length_type const& rl ) const = 0;
    
    virtual bool bounced(position_type const& old_pos, position_type const& new_pos, length_type const& dist_to_surface, length_type const& particle_radius) const = 0;
    
    virtual bool in_reaction_volume( length_type const& dist_to_surface, length_type const& particle_radius, length_type const& rl ) const = 0;

    virtual projected_type projected_point(position_type const& pos) const = 0;
    
    virtual projected_type projected_point_on_surface(position_type const& pos) const = 0;
    
    virtual length_type distance(position_type const& pos) const = 0;
    
    virtual position_flag_pair_type deflect(position_type const& pos0, position_type const& displacement) const = 0;
    
    virtual position_type const& structure_position() const = 0;

    virtual std::size_t hash() const
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<identifier_type>()(id_) ^
               hash<structure_type_id_type>()(sid_);
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Structure(" << real_id() << ", " << sid() << ")";
        return out.str();
    }

    // Constructor
    Structure(identifier_type const& id, structure_type_id_type const& sid)
        : id_(id), sid_(sid) {}

////// Member variables
protected:
    structure_id_type       real_id_;   // id of the structure
    identifier_type         id_;        // This is now just the name
    structure_type_id_type  sid_;       // id of the structure_type of the structure
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
