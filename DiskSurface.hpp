#ifndef DISK_SURFACE_HPP
#define DISK_SURFACE_HPP

#include <boost/bind.hpp>
#include "Surface.hpp"
#include "Disk.hpp"
#include "freeFunctions.hpp"
#include "StructureFunctions.hpp"
#include "geometry.hpp"

template <typename Tobj_, typename Tid_, typename Ttraits_>
class StructureContainer;


template<typename Ttraits_>
class DiskSurface
    : public BasicSurfaceImpl<Ttraits_, Disk<typename Ttraits_::length_type> >
{
    // The DiskSurface is the implementation of a Basic surface parameterized with a Disk

public:
    typedef BasicSurfaceImpl<Ttraits_, Disk<typename Ttraits_::length_type> > base_type;
    typedef Ttraits_    traits_type;
    typedef typename base_type::structure_name_type     structure_name_type;
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::shape_type              shape_type;
    typedef typename base_type::rng_type                rng_type;
    typedef typename base_type::position_type           position_type;
    typedef typename base_type::length_type             length_type;
    typedef typename base_type::side_enum_type          side_enum_type;
    typedef typename traits_type::species_type          species_type;
    typedef typename traits_type::structure_type        structure_type;

    typedef StructureContainer<Structure<traits_type>, structure_id_type, traits_type>    structure_container_type;

    typedef std::pair<position_type, position_type>                               position_pair_type;
    typedef std::pair<position_type, structure_id_type>                           position_structid_pair_type;
    typedef std::pair<position_structid_pair_type, position_structid_pair_type>   position_structid_pair_pair_type;

    
    /*** Simple structure-specific sampling functions ***/
    // Produce a "random position" in the disk, which is always its center (the only legal pos.)
    virtual position_type random_position(rng_type& rng) const        
    {
        return ::random_position(base_type::shape(), boost::bind(&rng_type::uniform, rng, -1., 1.));
        // return base_type::shape().position(); // TODO maybe this variant is better
    }

    // Procude a "random vector" on the disk; returns the same as random_position()
    virtual position_type random_vector(length_type const& r, rng_type& rng) const    
    {
        return base_type::shape().position();
    }

    // BD displacement = zero vector because there is only one legal position on the disk
    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    {
        return multiply(base_type::shape().unit_z(), 0.0);  // TODO is there not cheaper way to pass a zero vector?
    }
    
    /*** New BD scheme functions ***/
    // Rate for binding to particle on the structure
    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const
    {
        // This rate does not make sense on the ("zero-dimensional") disk object
        // but is used to determine the maximal time step in a Multi object and
        // therefore should be defined in a proper way.
        return 0.0;
    }
    
    // Rate for binding to the structure
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const
    {
        return k;
    }

    // Reaction volume for binding to particle in the structure
    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const
    {
        // The disk can only hold 1 particle; this function therefore never should be called.
        return rl;
    }

    // Reaction volume for binding to the structure
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const
    {        
        // The reaction volume for a particle on the rod interacting with a disk;
        // should be the same as for two particles interacting with each other on the rod.
        return rl;
    }
    
    // Vector of dissociation from the structure into the bulk
    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        // This function produces a position for a particle unbinding from the disk onto the rod.
        // It should lie within the reaction volume around the disk.
        // TODO: THIS IS THE CODE COPIED FROM CylindricalSurface
        // FIXME WTF is all this calculation!! The dissociation vector is supposed to be super simple!
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
    
    // Normed direction of dissociation from the structure to parent structure
    virtual position_type surface_dissociation_unit_vector( rng_type& rng ) const
    {
        return base_type::shape().unit_z();
    }

    // Positions created at dissociation of one particle on the structure into two particles on the structure
    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, 
        length_type const& rl ) const
    {
        // The positions of a particle dissociating into two new ones on the disk; should never happen,
        // therefore this function just returns a dummy positions pair (2x the disk center)        
        
        return position_pair_type( base_type::shape().position(), base_type::shape().position() );
    }
    
    // Positions created at dissociation of one particle on the structure into two particles, one of which ends up in the bulk
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, 
        position_type const& op_surf, length_type const& rl ) const
    {
        // TODO 
        // This function produces two new positions for a dissociating particle in the case
        // that one stays on the surface and the other one changes to the parent structure.
        // We have to distinguish between sink and cap here and between dissociation onto
        // the rod and into the bulk for the cap!
        
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
   
    // Used by newBDPropagator
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const
    {
        const length_type disk_radius (base_type::shape().radius());
        const boost::array<length_type, 2> new_pos_rz(::to_internal(base_type::shape(), new_pos));
        const boost::array<length_type, 2> old_pos_rz(::to_internal(base_type::shape(), old_pos));
        if (new_pos_rz[1] * old_pos_rz[1] < 0 &&
            ( (new_pos_rz[0] < disk_radius ) || (old_pos_rz[0] < disk_radius) ) )
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
        // TODO
        length_type cylinder_radius = base_type::shape().radius();
        // Return minimal distance *to* surface.
        return (cylinder_radius + radius) * traits_type::MINIMAL_SEPARATION_FACTOR - cylinder_radius;
    }
*/
    /*** Boundary condition handling ***/
    // FIXME This is a mess but it works. See ParticleContainerBase.hpp for explanation.
    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id,
                                                       structure_container_type const& structure_container) const
    {
        return pos_struct_id;   // This seems a little strange, but we assume that particles are immobile on the disk
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_container_type const& structure_container) const
    {
        return pos_struct_id;   // Disks can also not be connected, so no cyclic transpose.
    }


    // *** Dynamic dispatch for the structure functions *** //
    // *** 1 *** - One new position
    // This requires a double dynamic dispatch.
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& target_structure, position_type const& position,
                                                         length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {
        return target_structure.get_pos_sid_pair_helper(*this, position, offset, reaction_length, rng);
    }
    // *** 2 *** - Two new positions
    // Same principle as above, but different return type
    virtual position_structid_pair_pair_type get_pos_sid_pair_pair(structure_type const& target_structure, position_type const& position,
                                                                   species_type const& s1, species_type const& s2, length_type const& reaction_length, rng_type const& rng) const
    {
        return target_structure.get_pos_sid_pair_pair_helper(*this, position, s1, s2, reaction_length, rng);
    }    
    // *** 3 *** - Pair reactions => two origin structures
    // Overloading get_pos_sid_pair with signature (origin_structure2, target_structure_type_id, ...)
    virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type_id_type const& target_sid,
                                                         position_type const& CoM, length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
    {
        return origin_structure2.get_pos_sid_pair_helper_two_origins(*this, target_sid, CoM, offset, reaction_length, rng);
    }    
//     // *** 4 *** - Generalized functions for pair reactions with two origin structures and one target structure
//     // NOTE: This is yet unused, but possibly useful in the future.
//     // Overloading get_pos_sid_pair again with signature (origin_structure2, target_structure, ...) and introducing
//     // a triple dynamic dispatch.
//     virtual position_structid_pair_type get_pos_sid_pair(structure_type const& origin_structure2, structure_type const& target_structure, position_type const& position,
//                                                          length_type const& offset, length_type const& reaction_length, rng_type const& rng) const
//     {
//         return origin_structure2.get_pos_sid_pair_helper1(*this, target_structure, position, offset, reaction_length, rng);
//     }

    
    /*** Formerly used functions of the Morelli scheme ***/
    // DEPRECATED
    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {
         // TODO: This is part of the old BD scheme and should be removed at some point
        return drawR_gbd_1D(rnd, r01, dt, D01, v);
    }
    // DEPRECATED
    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {
        // TODO: This is part of the old BD scheme and should be removed at some point
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
    // DEPRECATED
    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        // TODO: This is part of the old BD scheme and should be removed at some point
        return random_vector(drawR_gbd(rng.uniform(0., 1.), r01, dt, D01, v), rng);
    }
    
    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }
        
    DiskSurface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id, shape) {}
};

#endif /* DISK_SURFACE_HPP */
