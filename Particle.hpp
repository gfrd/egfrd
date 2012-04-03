#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include "Sphere.hpp"
#include "Shape.hpp"

template<typename T_, typename Td_, typename Tsid_, typename Tstructure_id_>
struct Particle
// T_ is a length_type
{
    // defining shorthands for used types
    typedef Sphere<T_>                          shape_type;
    typedef Td_                                 D_type;
    // the drift v has the same type as diffusion constant D for now, may be generalized at a later stage
    typedef Td_                                 v_type;
    typedef Tsid_                               species_id_type;
    typedef typename shape_type::position_type  position_type;
    typedef typename shape_type::length_type    length_type;
    typedef Tstructure_id_                      structure_id_type;    // The structure (not structure_type) on which the particle currently lives.

    // constructors
    // Note that the structure_id is mandatory (?)
    Particle(): shape_(), species_id_(), structure_id_(), D_(0.), v_(0.) {}

/*    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D)
        : shape_(shape), species_id_(species_id), structure_id_(), D_(D), v_(0.){}

    Particle(species_id_type const& species_id, shape_type const& shape,
             D_type const& D, v_type const& v)
        : shape_(shape), species_id_(species_id), structure_id_(), D_(D), v_(v) {}
*/
    Particle(species_id_type const& species_id, shape_type const& shape,
             structure_id_type const& structure_id, D_type const& D)
        : shape_(shape), species_id_(species_id), structure_id_(structure_id), D_(D), v_(0.){}

    Particle(species_id_type const& species_id, shape_type const& shape,
             structure_id_type const& structure_id, D_type const& D, v_type const& v)
        : shape_(shape), species_id_(species_id), structure_id_(structure_id), D_(D), v_(v) {}


    // Get basic properties such as position, radius, D, v (drift)
    position_type& position()
    {
        return shape_.position();
    }

    position_type const& position() const
    {
        return shape_.position();
    }

    length_type& radius()
    {
        return shape_.radius();
    }

    length_type const& radius() const
    {
        return shape_.radius();
    }

    D_type& D()
    {
        return D_;
    }

    D_type const& D() const
    {
        return D_;
    }
    
    v_type& v()
    {
        return v_;
    }

    v_type const& v() const
    {
        return v_;
    }

    // Get the shape of the particle (Sphere)
    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
    }

    // Get the species_id
    species_id_type const& sid() const
    {
        return species_id_;
    }

    species_id_type& sid()
    {
        return species_id_;
    }

    // Get the id of the structure that the particle lives on.
    structure_id_type const& structure_id() const
    {
        return structure_id_;
    }

    structure_id_type& structure_id()
    {
        return structure_id_;
    }

    // Check equality/inequality
    bool operator==(Particle const& rhs) const
    {
        return species_id_ == rhs.sid() &&
               shape_ == rhs.shape() &&
               structure_id_ == rhs.structure_id();
    }

    bool operator!=(Particle const& rhs) const
    {
        return !operator==(rhs);
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

////// Member variables
private:
    shape_type          shape_;
    species_id_type     species_id_;
    structure_id_type   structure_id_;
    D_type              D_;
    v_type              v_;
};


///// Inline functions
template<typename Tstrm_, typename Ttraits_, typename T_, typename Td_, typename Tsid_, typename Tstructure_id_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const Particle<T_, Td_, Tsid_, Tstructure_id_>& p)
{
    strm << "Particle(" << p.shape() << ", D=" << p.D() << ", v=" << p.v() << ", " << p.sid() << ", " << p.structure_id() << ")";
    return strm;
}



#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_, typename Td_, typename Tsid_, typename Tstructure_id_>
struct hash<Particle<T_, Td_, Tsid_, Tstructure_id_> >
{
    typedef Particle<T_, Td_, Tsid_, Tstructure_id_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
               hash<typename argument_type::length_type>()(val.radius()) ^
               hash<typename argument_type::D_type>()(val.D()) ^
               hash<typename argument_type::species_id_type>()(val.sid()) ^
               hash<typename argument_type::structure_id_type>()(val.structure_id());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* PARTICLE_HPP */
