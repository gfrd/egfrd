#ifndef CYLINDRICAL_SURFACE_HPP
#define CYLINDRICAL_SURFACE_HPP

#include <boost/bind.hpp>
#include "Surface.hpp"
#include "Cylinder.hpp"
#include "freeFunctions.hpp"
#include "geometry.hpp"

template<typename Ttraits_>
class CylindricalSurface
    : public BasicSurfaceImpl<Ttraits_, Cylinder<typename Ttraits_::world_type::traits_type::length_type> >
{
    // The CylindricalSurface is the implementation of a Basic surface parameterized with a Cylinder

public:
    typedef BasicSurfaceImpl<Ttraits_, Cylinder<typename Ttraits_::world_type::traits_type::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::structure_name_type     structure_name_type;        // This is just the name of the structure
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::shape_type              shape_type;
    typedef typename base_type::rng_type                rng_type;
    typedef typename base_type::position_type           position_type;
    typedef typename base_type::length_type             length_type;
    typedef typename Ttraits_::world_type::species_type species_type;
    typedef std::pair<position_type, position_type>     position_pair_type;

    virtual position_type random_position(rng_type& rng) const
    // Gives a random position in the cylinder (1D)
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    // Gives a random vector in the axis of the cylinder
    {
        return multiply(base_type::shape().unit_z(),
                (rng.uniform_int(0, 1) * 2 - 1) * r);
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    // BD displacement for BD on the line of the cylinder?
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
    
    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const
    {
        return k;
    }
    
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const
    {
        return k / ( 2 * M_PI * (base_type::shape().radius() + r0 ) );
    }

    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const
    {
        return rl;
    }

    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const
    {
        length_type rc( base_type::shape().radius() + r0 ); 
        length_type rcl( rc + rl );
        length_type rcl_sq( rcl * rcl );
        length_type rc_sq( rc * rc );

        return M_PI * ( rcl_sq - rc_sq );
    }
    
    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        Real X( rng.uniform(0.,1.) );
        length_type const rod_radius = base_type::shape().radius();
        position_type const unit_z = base_type::shape().unit_z();

        length_type const rrl( rod_radius + r0 + rl );
        length_type const rrl_sq( gsl_pow_2(rrl) );
        length_type const rr_sq( gsl_pow_2(rod_radius + r0) );
        
        length_type const diss_vec_length( sqrt( X * (rrl_sq - rr_sq) + rr_sq ) );

        position_type v(rng.uniform(0.,1.) - .5, rng.uniform(0.,1.) - .5, rng.uniform(0.,1.) - .5);        
        v = normalize( subtract(v, multiply( unit_z, dot_product( unit_z, v ) ) ) );
         
        return multiply( v, diss_vec_length); 
    }

    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, 
        length_type const& rl ) const
    {
        length_type const r01( s0.radius() + s1.radius() );
        Real const D01( s0.D() + s1.D() );
        
        Real const X( rng.uniform(0.,1.) );
        
        length_type const diss_vec_length( X*rl + r01 );

        position_type const m( random_vector( diss_vec_length, rng ) );
        
        return position_pair_type( op - m * s0.D() / D01,
                                    op + m * s1.D() / D01 );
    }
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, 
        position_type const& op_surf, length_type const& rl ) const
    {
        Real const rod_radius( base_type::shape().radius() );
        
        //Species living on the rod should have a larger radius than the rod.
        assert( rod_radius < s_surf.radius() );
    
        length_type const r01( s_bulk.radius() + s_surf.radius() );
        Real const D01( s_bulk.D() + s_surf.D() );
        Real const D_bulk_D01( s_bulk.D() / D01 );
        Real const D_surf_D01( s_surf.D() / D01 );
        
        //Commented code for direct binding case with c.o.m. reaction.
        //Real const theta_min( asin(rod_radius / r01) );       
        Real const theta_min( asin( (rod_radius + s_bulk.radius()) / r01) );       
        Real const theta( theta_min + rng.uniform(0.,1.) * (M_PI - 2 * theta_min) );  
        Real const phi( rng.uniform(0.,1.) * 2 * M_PI );
        
        Real const X( rng.uniform(0.,1.) );
        length_type const r01l( r01 + rl );
        length_type const r01l_cb( r01l * r01l * r01l );
        length_type const r01_cb( r01 * r01 * r01 );
        
        length_type const diss_vec_length( cbrt( X * (r01l_cb - r01_cb) + r01_cb ) );
        
        position_type v;
        v[0] = 1.; v[1] = 1.; v[2] = 1.;
        
        position_type const unit_z( base_type::shape().unit_z() );      
        position_type const unit_x( normalize( subtract( v, 
                                        multiply( unit_z, dot_product(v, unit_z) ) ) ) );
        position_type const unit_y( normalize( cross_product( unit_x, unit_z ) ) );
                
        length_type const x( diss_vec_length * sin( theta ) * cos( phi ) );
        length_type const y( diss_vec_length * sin( theta ) * sin( phi ) );
        length_type const z( diss_vec_length * cos( theta ) );
                                
        position_pair_type pp01;
        
        pp01.first = subtract( op_surf, unit_z * (z * D_surf_D01) );
        
        pp01.second = add( op_surf, 
                        add( unit_x * x,
                          add( unit_y * y, 
                                 unit_z * (z * D_bulk_D01) ) ) );
        
        return pp01;
    }

    /* Determine if particle has crossed the 'surface' of the structure. */
    virtual bool bounced(position_type const& old_pos, position_type const& new_pos, 
        length_type const& dist_to_surface, length_type const& particle_radius) const
    {       
        if( dist_to_surface < particle_radius )
            return true;
        else
            return false;
    }
    
    virtual bool in_reaction_volume( length_type const& dist_to_surface, length_type const& particle_radius, length_type const& rl ) const
    {
        if( dist_to_surface - particle_radius <= rl )
            return true;
        else
            return false;
    }

    // This should replace above two methods.
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const
    {
        return base_type::distance(new_pos);
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

    CylindricalSurface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id, shape) {}
};

#endif /* CYLINDRICAL_SURFACE_HPP */
