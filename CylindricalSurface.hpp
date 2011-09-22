#ifndef CYLINDRICAL_SURFACE_HPP
#define CYLINDRICAL_SURFACE_HPP

#include <boost/bind.hpp>
#include <iostream>
#include "Surface.hpp"
#include "Cylinder.hpp"
//#include "freeFunctions.hpp"

template<typename Ttraits_>
class CylindricalSurface
    : public BasicSurfaceImpl<Ttraits_, Cylinder<typename Ttraits_::world_type::traits_type::length_type> >
{
public:
    typedef BasicSurfaceImpl<Ttraits_, Cylinder<typename Ttraits_::world_type::traits_type::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(),
                (rng.uniform_int(0, 1) * 2 - 1) * r);
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(), rng.normal(mean, r));
    }

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {
        return drawR_gbd_1D(rnd, r01, dt, D01, v);
    }

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {
        /*
            The I_bd factors used for calculating the acceptance probability are dependent on the direction 
            of the overlap step (r = r_1 - r_0), compared to the direction of the drift. 
            The I_bd factors are defined for a particle creating an overlap comming from the right (r < 0).
            Since the I_bd terms calulated here are for the backward move, we have to invert their drifts.
            When the particle comes from the left (r > 0) we have to invert its drift again.

            ---Code below is used for drift dependent backstep.

            Real numerator = g_bd_1D(ipv, r01, dt, D0, -v0);
            Real denominator = g_bd_1D(ipv, r01, dt, D0, v0)*exp( ipv/abs_ipv*(abs_ipv - r01)*v/D01 );
            Real correction = numerator/denominator;
            if( ipv < 0 )
                return correction*( k_a * dt / ( I_bd_1D(r01, dt, D0, -v0) + I_bd_1D(r01, dt, D1, v1) ) );
            else
                return correction*( k_a * dt / ( I_bd_1D(r01, dt, D0, v0) + I_bd_1D(r01, dt, D1, -v1) ) );

            Also change v -> -v in drawR for the dissociation move.
        */

        return 0.5*( k_a * dt / ( I_bd_1D(r01, dt, D0, v0) + I_bd_1D(r01, dt, D1, v1) ) );
  
    }

    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return random_vector(drawR_gbd(rng.uniform(0., 1.), r01, dt, D01, v), rng);
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        length_type cylinder_radius = base_type::shape().radius();
        // Return minimal distance *to* surface.
        return (cylinder_radius + radius) * traits_type::MINIMAL_SEPARATION_FACTOR - cylinder_radius;
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    CylindricalSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* CYLINDRICAL_SURFACE_HPP */
