#ifndef PLANE_HPP
#define PLANE_HPP

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <utility>
#include <algorithm>
#include "utils/array_helper.hpp"
#include "Shape.hpp"
#include "linear_algebra.hpp"

template<typename T_>
class Plane
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Plane(position_type const& position = position_type())
        : position_(position),
          units_(array_gen(
            create_vector<position_type>(1., 0., 0.),
            create_vector<position_type>(0., 1., 0.),
            create_vector<position_type>(0., 0., 1.))),
          half_extent_(array_gen<length_type>(0.5, 0.5)) {}

    template<typename Tarray_>
    Plane(position_type const& position, Tarray_ const& half_extent)
        : position_(position),
          units_(array_gen(
            create_vector<position_type>(1., 0., 0.),
            create_vector<position_type>(0., 1., 0.),
            create_vector<position_type>(0., 0., 1.)))
    {
        std::copy(boost::begin(half_extent), boost::end(half_extent),
                  boost::begin(half_extent_));
    }

    template<typename Tarray1, typename Tarray2>
    Plane(position_type const& position,
        Tarray1 const& units, Tarray2 const& half_extent)
        : position_(position)
    {
        std::copy(boost::begin(units), boost::end(units),
                  boost::begin(units_));
        std::copy(boost::begin(half_extent), boost::end(half_extent),
                  boost::begin(half_extent_));
    }

    template<typename Tarray_>
    Plane(position_type const& position,
        position_type const& vx,
        position_type const& vy,
        Tarray_ const& half_extent = array_gen<length_type>(0.5, 0.5))
        : position_(position), units_(array_gen(vx, vy, cross_product(vx, vy)))
    {
        std::copy(boost::begin(half_extent), boost::end(half_extent),
                  boost::begin(half_extent_));
    }

    Plane(position_type const& position,
        position_type const& vx,
        position_type const& vy,
        length_type const& half_lx,
        length_type const& half_ly)
        : position_(position), units_(array_gen(vx, vy, cross_product(vx, vy))),
          half_extent_(array_gen<length_type>(half_lx, half_ly)) {}

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    position_type const& unit_x() const
    {
        return units_[0];
    }

    position_type& unit_x()
    {
        return units_[0];
    }

    position_type const& unit_y() const
    {
        return units_[1];
    }

    position_type& unit_y()
    {
        return units_[1];
    }

    position_type const& unit_z() const
    {
        return units_[2];
    }

    position_type& unit_z()
    {
        return units_[2];
    }

    boost::array<position_type, 3> const& units() const
    {
        return units_;
    }

    boost::array<position_type, 3>& units()
    {
        return units_;
    }

    length_type const Lx() const
    { 
        return 2 * half_extent_[0];
    }

    length_type Lx()
    {
        return 2 * half_extent_[0];
    }

    length_type const Ly() const
    {
        return 2 * half_extent_[1];
    }

    length_type Ly()
    {
        return 2 * half_extent_[1];
    }

    boost::array<length_type, 2> const& half_extent() const
    {
        return half_extent_;
    }

    boost::array<length_type, 2>& half_extent()
    {
        return half_extent_;
    }

    bool operator==(const Plane& rhs) const
    {
        return position_ == rhs.position_ && units_ == rhs.units_ &&
               half_extent_ == rhs.half_extent_;
    }

    bool operator!=(const Plane& rhs) const
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

protected:
    position_type position_;
    boost::array<position_type, 3> units_;
    boost::array<length_type, 2> half_extent_;
};

template<typename T_>
inline boost::array<typename Plane<T_>::length_type, 3>
to_internal(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
// The function calculates the coefficients to express 'pos' into the base of the plane 'obj'
{
    typedef typename Plane<T_>::position_type position_type;
    position_type pos_vector(subtract(pos, obj.position()));

    return array_gen<typename Plane<T_>::length_type>(
        dot_product(pos_vector, obj.unit_x()),
        dot_product(pos_vector, obj.unit_y()),
        dot_product(pos_vector, obj.unit_z()));
}

template<typename T_>
inline std::pair<typename Plane<T_>::position_type,
                 typename Plane<T_>::length_type>
projected_point(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
// Calculates the projection of 'pos' onto the plane 'obj' and also returns the coefficient
// for the normal component (z) of 'pos' in the basis of the plane
{
    boost::array<typename Plane<T_>::length_type, 3> x_y_z(to_internal(obj, pos));
    return std::make_pair(
        add(add(obj.position(), multiply(obj.unit_x(), x_y_z[0])),
                                multiply(obj.unit_y(), x_y_z[1])),
        x_y_z[2]);
}

template<typename T_>
inline std::pair<typename Plane<T_>::position_type,
                 typename Plane<T_>::length_type>
projected_point_on_surface(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
// Since the projected point on the plane, is already on the surface of the plane,
// this function is just a wrapper of projected point.
{
    return projected_point(obj, pos);
}

template<typename T_>
inline typename Plane<T_>::length_type
distance(Plane<T_> const& obj, typename Plane<T_>::position_type const& pos)
// Calculates the distance from 'pos' to plane 'obj' Note that when the plane is finite,
// and also calculates the distance to the edge of the plane if necessary
{
    typedef typename Plane<T_>::length_type length_type;
    boost::array<length_type, 3> const x_y_z(to_internal(obj, pos));

    length_type const dx(subtract( abs(x_y_z[0]), obj.half_extent()[0]));
    length_type const dy(subtract( abs(x_y_z[1]), obj.half_extent()[1]));

    if (dx < 0 && dy < 0) {
        // pos is positioned over the plane (projected point is in the plane and
	    // not next to it).
        return abs(x_y_z[2]);
    }

    if (dx > 0) // outside the plane in the x direction
    {
        if (dy > 0)
        {
            // outside the plane in both x and y direction
            return std::sqrt(gsl_pow_2(dx) + gsl_pow_2(dy) +
                             gsl_pow_2(x_y_z[2]));
        }
        else
        {
	    // outside the plane in x, but inside in y direction
            return std::sqrt(gsl_pow_2(dx) + gsl_pow_2(x_y_z[2]));
        }
    }
    else   // inside the plane in x direction
    {
        if (dy > 0)
        {
	    // outside the plane in y, but inside in x direction
            return std::sqrt(gsl_pow_2(dy) + gsl_pow_2(x_y_z[2]));
        }
        else
        {
            // inside the plane in both x and y direction (see above)
            return abs(x_y_z[2]);
        }
    }
}

template<typename T, typename Trng>
inline typename Plane<T>::position_type
random_position(Plane<T> const& shape, Trng& rng)
{
    typedef typename Plane<T>::length_type length_type;

    // -1 < rng() < 1. See for example PlanarSurface.hpp.
    return add(
        shape.position(),
        add(multiply(shape.units()[0], shape.half_extent()[0] * rng()),
            multiply(shape.units()[1], shape.half_extent()[1] * rng())));
}

template<typename T>
inline Plane<T> const& shape(Plane<T> const& shape)
{
    return shape;
}

template<typename T>
inline Plane<T>& shape(Plane<T>& shape)
{
    return shape;
}

template<typename T_>
struct is_shape<Plane<T_> >: public boost::mpl::true_ {};

template<typename T_>
struct shape_position_type<Plane<T_> >
{
    typedef typename Plane<T_>::position_type type;
};

template<typename Tstrm_, typename Ttraits_, typename T_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm,
        const Plane<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.unit_x() << ", " << v.unit_y() << "," << v.Lx() << ", " << v.Ly() << "}";
    return strm;
}


#if defined(HAVE_TR1_FUNCTIONAL)
namespace std { namespace tr1 {
#elif defined(HAVE_STD_HASH)
namespace std {
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost {
#endif

template<typename T_>
struct hash<Plane<T_> >
{
    typedef Plane<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::position_type>()(val.unit_x()) ^
            hash<typename argument_type::position_type>()(val.unit_y()) ^
            hash<typename argument_type::length_type>()(val.half_extent()[0]) ^
            hash<typename argument_type::length_type>()(val.half_extent()[1]);
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} } // namespace std::tr1
#elif defined(HAVE_STD_HASH)
} // namespace std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // namespace boost
#endif

#endif /* PLANE_HPP */
