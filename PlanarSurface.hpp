#ifndef PLANAR_SURFACE_HPP
#define PLANAR_SURFACE_HPP

#include <boost/bind.hpp>
#include "Surface.hpp"
#include "Plane.hpp"
#include "freeFunctions.hpp"

template <typename Tobj_, typename Tid_, typename Ttraits_>
class StructureContainer;


template<typename Ttraits_>
class PlanarSurface
    : public BasicSurfaceImpl<Ttraits_, Plane<typename Ttraits_::length_type> >
{
    // The planar surface is an implementation of a Basic surface parameterized with the plane

public:
    typedef BasicSurfaceImpl<Ttraits_, Plane<typename Ttraits_::length_type> > base_type;
    typedef Ttraits_ traits_type;
    typedef typename base_type::structure_name_type     structure_name_type;        // This is just the name of the structure
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::shape_type              shape_type;
    typedef typename base_type::rng_type                rng_type;
    typedef typename base_type::position_type           position_type;
    typedef typename base_type::length_type             length_type;
    typedef typename base_type::side_enum_type          side_enum_type;
    typedef typename traits_type::species_type          species_type;

    typedef StructureContainer<Structure<traits_type>, structure_id_type, traits_type>    structure_container_type;

    typedef std::pair<position_type, position_type>         position_pair_type;
    typedef std::pair<position_type, structure_id_type>     position_structid_pair_type;


    virtual position_type random_position(rng_type& rng) const
    // Selects a random position in the plane
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    // Provides a random vector in the plane
    {
        return multiply(
            normalize(
                add(
                    multiply(
                        base_type::shape().units()[0], rng.uniform(-1., 1.)),
                    multiply(
                        base_type::shape().units()[1], rng.uniform(-1., 1.)))), r);
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    // Calculated a bd displacement for particles diffusing in the plane?
    {
        length_type const x(rng.normal(0., r)), y(rng.normal(0., r));
        return add(
            multiply(base_type::shape().unit_x(), x),
            multiply(base_type::shape().unit_y(), y));
    }

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {
        //TODO: use the 2D BD function instead of the 3D one - failed on very hard integral.
        return drawR_gbd_3D(rnd, r01, dt, D01);
    }

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {
        //TODO: use the 2D BD function instead of the 3D one. - Solution known
        return k_a * dt / ((I_bd_3D(r01, dt, D0) + I_bd_3D(r01, dt, D1)) * 4.0 * M_PI);
    }

    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return random_vector( drawR_gbd(rng(), r01, dt, D01, v), rng ); 
    }

    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const
    {
        return k / (2 * M_PI * r01);    
    }
    
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const
    {
        return k;    
    }

    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const
    {
        length_type r01l( r01 + rl );
        length_type r01l_sq( r01l * r01l );
        length_type r01_sq( r01 * r01 );

        return M_PI * ( r01l_sq - r01_sq );
    }
    
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const
    {
        //factor 2 because surface has two possible reaction sides due to cyclic BCn.
        return 2*rl;
    }
    
    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        Real X( rng.uniform(0.,1.) );
        length_type diss_vec_length( X*rl );

//        position_type unit_z( cross_product( base_type::shape().unit_x(), base_type::shape().unit_y() ) );
//
//        unit_z = normalize ( unit_z );

        return multiply( base_type::shape().unit_z(), (rng.uniform_int(0, 1) * 2 - 1) * diss_vec_length );
    }

    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, 
        length_type const& rl ) const
    {
        length_type const r01( s0.radius() + s1.radius() );
        Real const D01( s0.D() + s1.D() );
    
        Real X( rng.uniform(0.,1.) );

        length_type const r01l( r01 + rl );
        length_type const r01l_sq( r01l * r01l );
        length_type const r01_sq( r01 * r01 );
        
        length_type const diss_vec_length( sqrt( X * (r01l_sq - r01_sq) + r01_sq ) );

        position_type const m( random_vector( diss_vec_length, rng ) );

        return position_pair_type( op - m * s0.D() / D01,
                                    op + m * s1.D() / D01 );
    }
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, 
        position_type const& op_surf, length_type const& rl ) const
    {
        length_type const r01( s_bulk.radius() + s_surf.radius() );
        Real const D01( s_bulk.D() + s_surf.D() );
        Real const D_bulk_D01( s_bulk.D() / D01 );
        Real const D_surf_D01( s_surf.D() / D01 );
        
        //Real theta_max( M_PI/2 - asin(s_bulk.radius() / ( r01 )) );
        //theta_max = theta_max < 0 ? 0 : theta_max;
        
        Real const theta( rng.uniform(0.,1.) * M_PI );        
        Real const phi( rng.uniform(0.,1.) * 2 * M_PI );
        
        Real const X( rng.uniform(0.,1.) );
        length_type const r01l( r01 + rl );
        length_type const r01l_cb( r01l * r01l * r01l );
        length_type const r01_cb( r01 * r01 * r01 );
        
        length_type const diss_vec_length( cbrt( X * (r01l_cb - r01_cb ) + r01_cb ) );   
        
        position_type unit_z( cross_product( base_type::shape().unit_x(), base_type::shape().unit_y() ) );
        unit_z = normalize ( unit_z );
        
        //unit_z = multiply(unit_z, rng.uniform_int(0, 1) * 2 - 1);
        
        length_type const x( diss_vec_length * sin( theta ) * cos( phi ) );
        length_type const y( diss_vec_length * sin( theta ) * sin( phi ) );
        length_type const z( diss_vec_length * cos( theta ) );
        
        position_pair_type pp01;
        
        pp01.first = subtract( op_surf, 
                       add( base_type::shape().unit_x() * (x * D_surf_D01),
                              base_type::shape().unit_y() * (y * D_surf_D01) ) );
        
        pp01.second = add( op_surf, 
                       add( base_type::shape().unit_x() * (x * D_bulk_D01),
                         add( base_type::shape().unit_y() * (y * D_bulk_D01),
                                unit_z * z ) ) );

        return pp01;
        
    }
    
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const
    {
        const boost::array<length_type, 2> half_lengths(base_type::shape().half_extent());
        const boost::array<length_type, 3> new_pos_xyz(::to_internal(base_type::shape(), new_pos));
        const boost::array<length_type, 3> old_pos_xyz(::to_internal(base_type::shape(), old_pos));
        if (new_pos_xyz[2] * old_pos_xyz[2] < 0 &&
            ((abs(new_pos_xyz[0]) < half_lengths[0] && abs(new_pos_xyz[1]) < half_lengths[1]) ||
             (abs(old_pos_xyz[0]) < half_lengths[0] && abs(old_pos_xyz[1]) < half_lengths[1])))
        {
            return -1.0 * base_type::distance(new_pos) + sigma;
        }
        else
        {
            return base_type::distance(new_pos) + sigma;
        }
    }
/*
    virtual length_type minimal_distance(length_type const& radius) const
    {
        // PlanarSurface has thickness of 0.
        return radius * traits_type::MINIMAL_SEPARATION_FACTOR;
    }
*/
    // FIXME This is a mess but it works. See ParticleContainerBase.hpp for explanation.
    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id,
                                                       structure_container_type const& structure_container) const
    {
        return structure_container.apply_boundary(*this, pos_struct_id);
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_container_type const& structure_container) const
    {
        return structure_container.cyclic_transpose(*this, pos_struct_id);
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    PlanarSurface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id, shape) {}
};


#endif /* PLANAR_SURFACE_HPP */
