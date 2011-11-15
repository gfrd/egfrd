#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include <cmath>
#include "Vector3.hpp"
#include "Shape.hpp"

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
template<typename T_>
class Cylinder
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Cylinder()
        : position_(), radius_(0), unit_z_(), half_length_(0) {}

    Cylinder(position_type const& position, length_type const& radius,
             position_type const& unit_z, length_type const& half_length )
        : position_(position), radius_(radius), unit_z_(unit_z),
          half_length_(half_length) {}

    bool operator==(const Cylinder& rhs) const
    {
        return position_ == rhs.position() && radius_ == rhs.radius() && unit_z_ == rhs.unit_z() && half_length_ == rhs.half_length();
    }

    bool operator!=(const Cylinder& rhs) const
    {
        return !operator==(rhs);
    }

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    position_type const& unit_z() const
    {
        return unit_z_;
    }

    position_type& unit_z()
    {
        return unit_z_;
    }

    length_type const& half_length() const
    {
        return half_length_;
    }

    length_type& half_length()
    {
        return half_length_;
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

private:
    position_type position_; // centre.
    length_type radius_;
    position_type unit_z_; // Z-unit_z. should be normalized.
    length_type half_length_;
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << ", " << v.unit_z() << ", " << v.half_length() << "}";
    return strm;
}

template<typename T_>
inline std::pair<typename Cylinder<T_>::length_type,
                 typename Cylinder<T_>::length_type>
to_internal(Cylinder<T_> const& obj, typename Cylinder<T_>::position_type const& pos)
{
    // Return pos relative to position of cylinder. 
    typedef typename Cylinder<T_>::position_type position_type;
    typedef typename Cylinder<T_>::length_type length_type;

    const position_type pos_vector(subtract(pos, obj.position()));
    // z can be < 0
    const length_type z(dot_product(pos_vector, obj.unit_z()));
    // r is always >= 0
    const length_type r(length(pos_vector - multiply(obj.unit_z(), z)));

    return std::make_pair(r, z);
}

template<typename T_>
inline std::pair<typename Cylinder<T_>::position_type,
                 typename Cylinder<T_>::length_type>
projected_point(Cylinder<T_> const& obj,
                typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::length_type length_type;

    // The projection lies on the z-axis.
    std::pair<length_type, length_type> r_z(to_internal(obj, pos));
    return std::make_pair(
        add(obj.position(), multiply(obj.unit_z(), r_z.second)),
        r_z.first);
}

//Allmost equal to projected point method, but for the substraction of the cylinder radius from the radial distance r.
//And projected point now lies on the surface, not on the central axis.
template<typename T_>
inline std::pair<typename Cylinder<T_>::position_type,
                 typename Cylinder<T_>::length_type>
projected_point_on_surface(Cylinder<T_> const& obj,
                typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::length_type length_type;
    typedef typename Cylinder<T_>::position_type position_type;

    std::pair<length_type, length_type> r_z(to_internal(obj, pos));
    position_type within_surface( add(obj.position(), multiply(obj.unit_z(), r_z.second)) );

    return std::make_pair( add(within_surface, multiply( normalize( subtract( pos, within_surface ) ), obj.radius() )),
        r_z.first - obj.radius() );
}

template<typename T_>
inline typename Cylinder<T_>::length_type
distance(Cylinder<T_> const& obj,
                typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::position_type position_type;
    typedef typename Cylinder<T_>::length_type length_type;

    /* First compute the (r,z) components of pos in a coordinate system 
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const std::pair<length_type, length_type> r_z(to_internal(obj, pos));

    /* Then compute distance to cylinder. */
    const length_type dz(std::fabs(r_z.second) - obj.half_length());
    const length_type dr(r_z.first - obj.radius());
    length_type distance;
    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r_z.first > obj.radius())
        {
            // Compute distance to edge.
            distance = std::sqrt( dz * dz + dr * dr );
        }
        else
        {
            distance = dz;
        }
    }
    else
    {
        if (dr > obj.radius())
        {
            // pos is somewhere 'parallel' to the cylinder.
            distance = dr;
        }
        else
        {
            // Inside cylinder. 
            distance = std::max(dr, dz);
        }
    }
    return distance;
}

template<typename T, typename Trng>
inline typename Cylinder<T>::position_type
random_position(Cylinder<T> const& shape, Trng& rng)
{
    // -1 < rng() < 1. See for example CylindricalSurface.hpp.
    return add(shape.position(),
               multiply(shape.unit_z(), rng() * shape.half_length()));
}

template<typename T_>
inline Cylinder<T_> const& shape(Cylinder<T_> const& shape)
{
    return shape;
}

template<typename T_>
inline Cylinder<T_>& shape(Cylinder<T_>& shape)
{
    return shape;
}

template<typename T_>
struct is_shape<Cylinder<T_> >: public boost::mpl::true_ {};

template<typename T_>
struct shape_position_type<Cylinder<T_> >
{
    typedef typename Cylinder<T_>::position_type type;
};

template<typename T_>
struct shape_length_type<Cylinder<T_> > {
    typedef typename Cylinder<T_>::length_type type;
};

template<typename T>
inline typename shape_length_type<Cylinder<T> >::type const& shape_size(Cylinder<T> const& shape)
{
    return shape.radius();
} 

template<typename T>
inline typename shape_length_type<Cylinder<T> >::type& shape_size(Cylinder<T>& shape)
{
    return shape.radius();
} 

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_>
struct hash<Cylinder<T_> >
{
    typedef Cylinder<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::position_type>()(val.unit_z()) ^
            hash<typename argument_type::length_type>()(val.half_length());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* CYLINDER_HPP */
