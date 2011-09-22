#ifndef SPHERICAL_SURFACE_HPP
#define SPHERICAL_SURFACE_HPP

#include "Surface.hpp"
#include "Sphere.hpp"

template<typename Ttraits_>
class SphericalSurface
    : public BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::world_type::traits_type::length_type> >
{
public:
    typedef BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::world_type::traits_type::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

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

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {    
        return length_type(); // TODO
    }

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {    
        return Real(); //TODO
    }

    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return position_type(); //TODO
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        return 0.; // TODO
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    SphericalSurface(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* SPHERICAL_SURFACE_HPP */
