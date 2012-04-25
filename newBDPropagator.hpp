#ifndef NEW_BD_PROPAGATOR_HPP
#define NEW_BD_PROPAGATOR_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include "Defs.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "freeFunctions.hpp"
#include "utils.hpp"
#include "utils/random.hpp"
#include "utils/get_default_impl.hpp"
#include "Logger.hpp"

#include <iostream>

template<typename Ttraits_>
class newBDPropagator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename Ttraits_::world_type::particle_container_type  particle_container_type;

    // shorthand typedef that we use
    typedef typename particle_container_type::length_type           length_type;
    typedef typename particle_container_type::position_type         position_type;
    typedef typename particle_container_type::species_type          species_type;
    typedef typename particle_container_type::species_id_type       species_id_type;
    typedef typename particle_container_type::particle_id_type      particle_id_type;
    typedef typename particle_container_type::particle_type         particle_type;
    typedef typename particle_container_type::particle_id_pair      particle_id_pair;
    typedef typename particle_container_type::particle_shape_type   particle_shape_type;
    typedef typename particle_container_type::structure_type        structure_type;
    typedef typename particle_container_type::structure_id_type     structure_id_type;    

    typedef typename particle_container_type::structure_id_and_distance_pair     structure_id_and_distance_pair;
    typedef typename particle_container_type::particle_id_pair_generator         particle_id_pair_generator;
    typedef typename particle_container_type::particle_id_pair_and_distance      particle_id_pair_and_distance;
    typedef typename particle_container_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
    typedef typename traits_type::world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type                         time_type;
    typedef typename traits_type::network_rules_type                network_rules_type;
    typedef typename network_rules_type::reaction_rules             reaction_rules;
    typedef typename network_rules_type::reaction_rule_type         reaction_rule_type;
    typedef typename traits_type::reaction_record_type              reaction_record_type;
    typedef typename traits_type::reaction_recorder_type            reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type               volume_clearer_type;

    typedef std::vector<particle_id_type>                           particle_id_vector_type;
    typedef std::pair<position_type, position_type>                 position_pair_type;
    typedef std::pair<position_type, length_type>                   projected_type;

public:
    // The constructor
    template<typename Trange_>
    newBDPropagator(
        particle_container_type& tx, network_rules_type const& rules,
        rng_type& rng, time_type dt, int max_retry_count, length_type reaction_length,
        reaction_recorder_type* rrec, volume_clearer_type* vc,
        Trange_ const& particle_ids)
        : tx_(tx), rules_(rules), rng_(rng), dt_(dt),
          max_retry_count_(max_retry_count), rrec_(rrec), vc_(vc),
          queue_(), rejected_move_count_(0), reaction_length_( reaction_length )
    {
        call_with_size_if_randomly_accessible(
            boost::bind(&particle_id_vector_type::reserve, &queue_, _1),
            particle_ids);

        // initialize the queue with the particle_ids of the particles to propagate.
        for (typename boost::range_const_iterator<Trange_>::type
                i(boost::begin(particle_ids)),
                e(boost::end(particle_ids)); i != e; ++i)
        {
            queue_.push_back(*i);
        }
        // randomize the queue
        shuffle(rng, queue_);
    }


    //// MAIN 'step' method 
    // This method is called to process the change of position and putative reaction of every particle.
    // The method only treats one particle at a time and returns 'True' as long as there are still particles
    //    in the container that need to be processed.
    bool operator()()
    {
        // if there are no more particle to treat -> end
        if (queue_.empty())
            return false;        

        //// 0.
        // Get the next particle from the queue
        particle_id_type pid(queue_.back());
        queue_.pop_back();
        particle_id_pair pp(tx_.get_particle(pid));
        // Log
        LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));
        


        /* Consider decay reactions first. Particles can move after decay, but this is done
           inside the 'attempt_single_reaction' function. */
        try
        {
            if (attempt_single_reaction(pp))
                return true;
        }
        catch (propagation_error const& reason)
        {
            log_.info("first-order reaction rejected (reason: %s)", reason.what());
            ++rejected_move_count_;
            return true;
        }
        

        //// 1.?
        // Get info of the particle
        const species_type                            pp_species(tx_.get_species(pp.second.sid()));
        const length_type                             r0( pp_species.radius() );
        const position_type                           old_pos( pp.second.position() );
        const structure_id_type                       old_struct_id(pp.second.structure_id());
        const boost::shared_ptr<const structure_type> pp_structure( tx_.get_structure( old_struct_id ) );
        // declare variables to use
        structure_id_and_distance_pair  struct_id_distance_pair;
        length_type                     particle_surface_distance (0);
        position_type                   new_pos;
        structure_id_type               new_structure_id( old_struct_id );

        //// 2. New position
        /* Sample a potential move, and check if the particle _core_ has overlapped with another 
           particle or surface. If this is the case the particle bounces, and is returned
           to it's original position. For a particle with D = 0, new_pos = old_pos */        
        if (pp_species.D() != 0.)
        {
            // get a new position, dependent on the structure the particle lives on.
            const position_type displacement( pp_structure->bd_displacement(pp_species.v() * dt_,
                                                                            std::sqrt(2.0 * pp_species.D() * dt_),
                                                                            rng_) );

            // TODO check if the new position is actually still in the structure
            // else
            // get 'connecting' structure and apply deflection. Should we do this recursively until the coordinate stays in the currents structure?
            // This actually kinda similar to the boundary condition application.
            new_pos = tx_.apply_boundary( add( old_pos, displacement ) );
/*            const structure_id_type target_struct_id(pp_structure->get_target_structure(new_pos));
            const structure_id_type target_struct_id(pp_structure->get_target_structure(new_pos));

            if (target_struct_id != pp_structure.id())
            // The particle jumped to a new position outside of the current surface.
            {
                const structure_type target_structure(tx_.get_structure(target_struct_id));
                new_structure_id = target_struct_id;
                const position_flag_pair_type newpos_flag (target_structure->deflect(old_pos, displacement));
                new_pos = newpos_flag.first;
            }
            // Else, just use the current information.
*/
        }
        else
        {
            new_pos = old_pos;
        }
        
        // Get all the particles that are in the reaction volume at the new position (and that may consequently also be in the core)
        /* Use a spherical shape with radius = particle_radius + reaction_length.
	       We use it to check for overlaps with other particles. */
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                tx_.check_overlap(particle_shape_type( new_pos, r0 + reaction_length_ ), pp.first));
        int particles_in_overlap(overlapped ? overlapped->size(): 0);
        bool bounced( false );


        //// 3.1 Check for core overlaps with particles
        /* Check if the particle at new_pos overlaps with any particle cores. */        
        int j( 0 );
        while(!bounced && j < particles_in_overlap)
            bounced = overlapped->at(j++).second < r0;
        
        //// 3.2 Check for core overlaps structures
        /* If the particle has not bounced with another particle and lives in the bulk: 
           check for core overlap with a surface. */
        // TODO always check for overlaps with surface ignoring its own surface (regardless of in bulk)
        if(!bounced && pp.second.structure_id() == tx_.get_def_structure_id())
        {
            // TODO Get all the structures within the reaction volume (and that may consequently also be in the core)
            struct_id_distance_pair = tx_.get_closest_surface( new_pos, new_structure_id);

            /* Check for overlap with nearest surface. */
            // TODO check all the surfaces while bounced != true.
            particle_surface_distance = struct_id_distance_pair.second;
            boost::shared_ptr<structure_type> closest_struct( tx_.get_structure( struct_id_distance_pair.first ) );
            bounced = closest_struct->bounced( tx_.cyclic_transpose( old_pos, closest_struct->position() ), 
                                               tx_.cyclic_transpose( new_pos, closest_struct->position() ), 
                                               particle_surface_distance, r0);
        }
         
        //// 4. Finalize the position
        /* If particle is bounced, restore old position and check reaction_volume for reaction partners at old_pos. */
        if(bounced)
        {
            // restore old position and structure_id
            new_pos = old_pos;
            new_structure_id = pp.second.structure_id();
            
            // re-get the reaction partners (particles), now on old position.
            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_after_bounce( 
                    tx_.check_overlap( particle_shape_type( new_pos, r0 + reaction_length_ ), pp.first) );
            overlapped.swap( overlapped_after_bounce );     // FIXME is there no better way?
                                            
            // NOTE that it is asserted that the particle overlap criterium for the particle with other particles
            // and surfaces is False!
            particles_in_overlap = overlapped ? overlapped->size(): 0; 
            
            // re-get the reaction partners (structures), at old position. TODO don't just do this when in the bulk.
            if(pp.second.structure_id() == tx_.get_def_structure_id())
            {
                // TODO get all near surfaces instead of just the closest.
                struct_id_distance_pair = tx_.get_closest_surface( new_pos, pp.second.structure_id() );
                // TODO get the number of surfaces.
                particle_surface_distance = struct_id_distance_pair.second;  
            }
        }
        

        //// 5.
        /* Attempt a reaction (and/or interaction) with all the particles (and/or a surface) that 
           are/is inside the reaction volume. */
        Real accumulated_prob = 0;
        Real const rnd( rng_() );
        
        /* First, if a surface is inside the reaction volume, and the particle is in the 3D attempt an interaction. */
        // TODO also interaction should be allowed when a particle is on a cylinder.
        // TODO don't check only the closest but check all overlapping surfaces.
        if( pp.second.structure_id() == tx_.get_def_structure_id())
        {
            // TODO get all the surfaces that overlap with the reaction volume (depends on current structure)
            boost::shared_ptr<structure_type> const closest_surf( tx_.get_structure( struct_id_distance_pair.first) );

            // TODO check that the projected point of the position of the particle is in the surface. If not, than no
            // reaction is possible.
            if( closest_surf->in_reaction_volume( particle_surface_distance, r0, reaction_length_ ) )
            {        
                accumulated_prob += k_total( pp.second.sid(), closest_surf->sid() ) * dt_ / 
                                    closest_surf->surface_reaction_volume( r0, reaction_length_ );
        
                if(accumulated_prob >= 1.)
                {
                    LOG_WARNING(("the acceptance probability of an interaction exeeded one; %f.",
                                 accumulated_prob));
                } 
                
                if( accumulated_prob > rnd )
                {            
                    try
                    {
                        LOG_DEBUG(("fire surface interaction"));
                        if(attempt_interaction(pp, closest_surf ))
                            return true;
                    }   
                    catch (propagation_error const& reason)
                    {
                        log_.info("surface interaction rejected (reason: %s)", reason.what());
                        ++rejected_move_count_;
                    }
                }
                else                
                {
                    LOG_DEBUG(("Particle attempted an interaction with the non-interactive surface %s.", 
                               boost::lexical_cast<std::string>(closest_surf->id()).c_str()));     // TODO is this the correct interpretation?
                }
            }

        }
        
        ////
        /* Now attempt a reaction with all particles inside the reaction volume. */
        j = 0;
        while(j < particles_in_overlap)
        {
            const particle_id_pair_and_distance & closest( overlapped->at(j) );

            species_type s0(pp_species);
            species_type s1(tx_.get_species(closest.first.second.sid()));
                
            /* If the structures of the reactants do not equal, one of the reactants has to come from the bulk,
               and we let this be s1, the particle from the surface is named s0. */
            if(s0.structure_type_id() != s1.structure_type_id())
            {
                if( !(s0.structure_type_id() == tx_.get_def_structure_type_id() || s1.structure_type_id() == tx_.get_def_structure_type_id())  )
                {
                    log_.warn("A surface particle overlapped with a particle living on a different surface. No reaction implemented.");
                    continue;
                }
                if(s0.structure_type_id() == tx_.get_def_structure_type_id())
                    std::swap(s0,s1);
            }
            
            const boost::shared_ptr<const structure_type> s1_struct( tx_.get_structure( closest.first.second.structure_id()) );
            const Real reaction_volume( s1_struct->particle_reaction_volume( s0.radius() + s1.radius(), reaction_length_ ) );
            accumulated_prob += k_total(s0.id(), s1.id()) * dt_ / ( 2. * reaction_volume ); 
            
            if (accumulated_prob >= 1.)
            {
                LOG_WARNING((
                    "the accumulated acceptance probability inside a reaction volume exeeded one; %f.",
                    accumulated_prob));
            } 
            
            if(accumulated_prob > rnd)
            {
                LOG_DEBUG(("fire reaction"));
    		    try
                {
                    if( attempt_pair_reaction(pp, closest.first, s0, s1) )
                        return true;
                }
                catch (propagation_error const& reason)
                {
                    log_.info("second-order reaction rejected (reason: %s)", reason.what());
                    ++rejected_move_count_;
                }
            }
            else                
            {
                LOG_DEBUG(("Particle attempted a reaction with particle %s and failed.", 
                    boost::lexical_cast<std::string>(closest.first.first).c_str()));    // TODO is this really a rejected move?
            }
            
            j++;
        }
                               
        /* If the particle did neither react, interact or bounce, update it to it's new position. */
        if(!bounced)
        {   
            const particle_id_pair particle_to_update( pp.first, 
                                                       particle_type(pp_species.id(),
                                                                     particle_shape_type(new_pos, r0),
                                                                     new_structure_id,
                                                                     pp_species.D(),
                                                                     pp_species.v()) );
                
            if (vc_)
            {
                if (!(*vc_)(particle_to_update.second.shape(), particle_to_update.first))
                {
                    log_.info("propagation move rejected.");
                    return true;
                }
            }
                
            tx_.update_particle(particle_to_update);
        }

        return true;
    }


    // Getter of rejected_move
    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

////// Private methods
private:
    bool attempt_single_reaction(particle_id_pair const& pp)
    // Handles all monomolecular reactions
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const Real rnd(rng_() / dt_);
        Real k_cumm = 0.;

        // select one of the available reaction rules that led to the monomolecular reaction
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            k_cumm += r.k();

            if (k_cumm > rnd)
            // We have found the reaction rule that is in effect
            {
                // get the products
                typename reaction_rule_type::species_id_range products(r.get_products());

                switch (::size(products))
                {
                    // When no product particles
                    case 0:
                    {
                        tx_.remove_particle(pp.first);  // pp was just popped of the local queue so only has to be removed from the upstream particlecontainer.
                        break;
                    }

                    // When one product particle
                    case 1:
                    {
                        const boost::shared_ptr<const structure_type> reactant_structure( tx_.get_structure(pp.second.structure_id()) );
                        const species_type product_species(tx_.get_species(products[0]));

                        // default values for the position and structure of the product
                        position_type     new_pos          (pp.second.position());
                        structure_id_type new_structure_id (pp.second.structure_id());
                        
                        /* If the product particle does NOT live on the same structure type as the reactant => surface -> bulk dissociation. */
                        if( product_species.structure_type_id() != reactant_structure->sid() )
                        {
                            // dissociation from a structure_type into another must always be TO the bulk (default structure/structure_type)
                            if( product_species.structure_type_id() != tx_.get_def_structure_type_id() )
                                throw not_implemented("No surface -> surface dissociation allowed");
                            new_structure_id = tx_.get_def_structure_id();
                        
                            // get dissociation position.
                            const position_type displacement( reactant_structure->surface_dissociation_vector(rng_, product_species.radius(), reaction_length_ ) );
                            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
                            
                            //Particle is allowed to move after dissociation from surface.
                            new_pos = make_move(product_species, new_pos, pp.first);
                        }
                      

                        const particle_shape_type new_shape(new_pos, product_species.radius());
                        // check for overlap with particles that are inside the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped(tx_.check_overlap(new_shape, pp.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space due to other particle");
                        }
                        // check for overlap with particles outside of the propagator (for example the eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(new_shape, pp.first))
                            {
                                throw propagation_error("no space due to other particle (vc)");
                            }
                        }

                        // check overlap with surfaces
                        if( product_species.structure_type_id() == tx_.get_def_structure_type_id() && pp.second.structure_id() == tx_.get_def_structure_id() )
                        {
                            const structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( new_pos, pp.second.structure_id()) );
                            const boost::shared_ptr<const structure_type> closest_struct( tx_.get_structure( struct_id_distance_pair.first ) );

                            if( closest_struct->bounced( tx_.cyclic_transpose( new_pos, closest_struct->position() ), 
                                                         tx_.cyclic_transpose( new_pos, closest_struct->position() ), 
                                                         struct_id_distance_pair.second, product_species.radius() ) )
                                throw propagation_error("no space due to near surface");
                        }

                        // Process changes
                        tx_.remove_particle(pp.first);
                        // Make new particle on right surface
                        const particle_id_pair product_particle( tx_.new_particle(product_species.id(), new_structure_id, new_pos) );

                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(product_particle.first), pp.first));
                        }
                        break;
                    }

                    // When two product particles
                    case 2:
                    {
                        const position_type reactant_pos(pp.second.position());
                        const boost::shared_ptr<const structure_type> reactant_structure( tx_.get_structure( pp.second.structure_id() ) );
                        
                        species_type product0_species(tx_.get_species(products[0]));
                        species_type product1_species(tx_.get_species(products[1]));
                        position_pair_type pos0pos1_pair;
                        structure_id_type prod0_struct_id;
                        structure_id_type prod1_struct_id;

                        /* Create positions (np0, np1) for the reaction products. 
                        The ipv between the products lies within the reation volume. */
                        int i (max_retry_count_);
                        while (true)
                        {
                            // Escape clause #1: Too many tries
                            if (i-- == 0)
                            {
                                throw propagation_error("no space due to particle or surface");
                            }
                            

                            //// Generate positions and structure_id for product particles after reaction

                            /* If the product particles do NOT live on the same structure => surface -> surface + bulk dissociation.
                               Else, the products live on the same structure AND on the same structure as the reactant. 
                               => standard geminate dissociation. 
                               
                            We set: species 0 lives on the surface; species 1 lives in the bulk. */
                            if( product0_species.structure_type_id() != product1_species.structure_type_id() )
                            {
                                if( !(product0_species.structure_type_id() == tx_.get_def_structure_type_id() || 
                                      product1_species.structure_type_id() == tx_.get_def_structure_type_id())  )
                                    throw not_implemented("No surface -> surface + surface - dissociation allowed");
                                                 
                                if(product0_species.structure_type_id() == tx_.get_def_structure_type_id())    
                                     std::swap(product0_species,product1_species);

                                pos0pos1_pair = reactant_structure-> special_geminate_dissociation_positions( rng_, product0_species, product1_species,
                                                                                                              reactant_pos, reaction_length_ );
                                prod0_struct_id = pp.second.structure_id();         // The 'original' structure
                                prod1_struct_id = tx_.get_def_structure_id();       // The bulk
                            }
                            else
                            {
                                pos0pos1_pair = reactant_structure-> geminate_dissociation_positions( rng_, product0_species, product1_species,
                                                                                                      reactant_pos, reaction_length_ );
                                prod0_struct_id = pp.second.structure_id();         // same structure as the reactant
                                prod1_struct_id = pp.second.structure_id();         // ditto
                            }
                            pos0pos1_pair.first  = tx_.apply_boundary( pos0pos1_pair.first );
                            pos0pos1_pair.second = tx_.apply_boundary( pos0pos1_pair.second );


                            // check overlap with particle in the propagator
                            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped_product0(
                                tx_.check_overlap( particle_shape_type(pos0pos1_pair.first, product0_species.radius()),  pp.first));
                            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped_product1(
                                tx_.check_overlap( particle_shape_type(pos0pos1_pair.second, product1_species.radius()), pp.first));
                                    
                            /* Only when both products live in the bulk, check for possible overlap with a surface. */
                            bool surface_bounce0( false );
                            bool surface_bounce1( false );
                            if( product0_species.structure_type_id() == tx_.get_def_structure_type_id() && product1_species.structure_type_id() == tx_.get_def_structure_type_id() )
                            {
                                structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( pos0pos1_pair.first, tx_.get_def_structure_id() ) );
                                boost::shared_ptr<structure_type> closest_struct( tx_.get_structure( struct_id_distance_pair.first ) );

                                surface_bounce0 = closest_struct->bounced( tx_.cyclic_transpose( pp.second.position(), closest_struct->position() ), 
                                                                           tx_.cyclic_transpose( pos0pos1_pair.first, closest_struct->position() ), 
                                                                           struct_id_distance_pair.second, product0_species.radius() );

                                struct_id_distance_pair = tx_.get_closest_surface( pos0pos1_pair.second, tx_.get_def_structure_id() );
                                closest_struct = tx_.get_structure( struct_id_distance_pair.first );

                                surface_bounce1 = closest_struct->bounced( tx_.cyclic_transpose( pp.second.position(), closest_struct->position() ), 
                                                                           tx_.cyclic_transpose( pos0pos1_pair.second, closest_struct->position() ), 
                                                                           struct_id_distance_pair.second, product1_species.radius() );
                            }
                            
                            // If the positions were allowed -> SUCCESS and leave the loop
                            if ( !(overlapped_product0 && overlapped_product0->size() > 0) && 
                                 !(overlapped_product1 && overlapped_product1->size() > 0) &&
                                 !surface_bounce0 && !surface_bounce1 )
                                break;
                        }

                        /* After dissociation both particles are allowed to move. */
                        pos0pos1_pair = make_pair_move(product0_species, product1_species, pos0pos1_pair, pp.first);

                        // check overlap with particle outside of propagator (for example in eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(pos0pos1_pair.first,  product0_species.radius()), pp.first) ||
                                !(*vc_)(particle_shape_type(pos0pos1_pair.second, product1_species.radius()), pp.first))
                            {
                                throw propagation_error("no space due to near particle (vc)");
                            }
                        }

                        // Process changes
                        tx_.remove_particle(pp.first);
                        const particle_id_pair product0_particle(tx_.new_particle(product0_species.id(), prod0_struct_id, pos0pos1_pair.first));
                        const particle_id_pair product1_particle(tx_.new_particle(product1_species.id(), prod1_struct_id, pos0pos1_pair.second));

                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(),
                                    array_gen(product0_particle.first, product1_particle.first),
                                    pp.first));
                        }
                        break;
                    }
                    default:
                        throw not_implemented("monomolecular reactions that produce more than two products are not supported");
                }
                // After we processed the reaction, we return -> no other reaction should be done.
                return true;
            }
        }
        // No monomolecular reaction has taken place. 
        return false;
    }


    const bool attempt_pair_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1,
                                     species_type const& s0, species_type const& s1)
    // This handles the reaction between a pair of particles.
    // pp0 is the initiating particle (also no longer in the queue) and pp1 is the target particle (still in the queue).
    // pp0, pp1 and s0,s1 can either be in 3D/2D/1D or mixed.
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.sid(), pp1.second.sid()));
        if(::size(rules) == 0)
            throw propagation_error("trying to fire a reaction between particles without a reaction rule.");
        
        const Real k_tot( k_total(pp0.second.sid(), pp1.second.sid()) );
        const Real rnd( k_tot * rng_() );
        Real k_cumm = 0;    
        
        // select one of the available reaction rules that led to the pair reaction
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& reaction_rule(*i);  
            k_cumm += reaction_rule.k();

            if (k_cumm >= rnd)
            // We have found the reaction rule that is in effect
            {
                // get the products
                const typename reaction_rule_type::species_id_range products(reaction_rule.get_products());

                switch (::size(products))
                {
                    // When no products
                    case 0:
                    {
                        tx_.remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        break;
                    }
                    // When one product particle
                    case 1:
                    {
//                        const species_id_type product(products[0]);
                        const species_type product_species(tx_.get_species(products[0]));
                        
                        position_type new_pos( tx_.apply_boundary(
                                        divide(
                                        add(multiply(pp0.second.position(), s1.D()),
                                            multiply(tx_.cyclic_transpose(
                                                pp1.second.position(),
                                                pp0.second.position()), s0.D())),
                                        (s0.D() + s1.D()))) );

//                        position_type new_pos (tx_.calculate_pair_CoM(pp0.second.position(), pp1.second.position()), s0.D(), s1.D());
                        structure_id_type new_structure_id;

                        //For unequal structures, project new_pos on surface.
                        if(product_species.structure_type_id() != s1.structure_type_id())
                        {
                            // make sure that pp1 is the particle that is in 3D.
                            // FIXME Ugly!!
                            boost::shared_ptr<structure_type> sp_struct;
                            if(pp0.second.structure_id() == tx_.get_def_structure_id())
                            {
                                new_structure_id = pp1.second.structure_id();
                            }
                            else
                            {
                                new_structure_id = pp0.second.structure_id();
                            }
                            sp_struct = tx_.get_structure( new_structure_id );
                            projected_type const position_on_surface( sp_struct->
                                                                      projected_point( tx_.cyclic_transpose( new_pos, 
                                                                                                             sp_struct->position() ) ) );
                            // check that the projected point is not outside of the product surface.
                            if (sp_struct->distance( position_on_surface.first) < 0)
                                throw propagation_error("position product particle was not in surface.");
                            new_pos = tx_.apply_boundary( position_on_surface.first );
                        }


                        // check for overlap with particles that are inside the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, product_species.radius()),
                                              pp0.first, pp1.first));
                        if( overlapped && overlapped->size() > 0 )
                        {
                            throw propagation_error("no space due to particle");
                        }
                        // check for overlap with particles outside of the propagator (for example the eGFRD simulator)
                        if( vc_ )
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, product_species.radius()), 
                                    pp0.first, pp1.first))
                            {
                                throw propagation_error("no space due to particle (vc)");
                            }
                        }
                        
                        // Check overlap with surfaces
                        if( product_species.structure_type_id() == tx_.get_def_structure_type_id() )
                        {
                            const structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( new_pos, tx_.get_def_structure_id() ) );
                            const boost::shared_ptr<const structure_type> closest_struct( tx_.get_structure( struct_id_distance_pair.first ) );

                            if( closest_struct->bounced( tx_.cyclic_transpose( pp0.second.position(), closest_struct->position() ), 
                                                         tx_.cyclic_transpose( pp1.second.position(), closest_struct->position() ), 
                                                         struct_id_distance_pair.second, product_species.radius() ) )
                                throw propagation_error("no space due to near surface");
                        }


                        // Process changes
                        tx_.remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        // Make new particle on right structure
                        particle_id_pair product_particle(tx_.new_particle(product_species.id(), new_structure_id, new_pos));
                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    reaction_rule.id(), array_gen(product_particle.first), pp0.first, pp1.first));
                        }
                        break;
                    }
                    // When number of products was not '0' or '1'
                    default:
                        throw not_implemented("bimolecular reactions that produce more than one product are not supported");
                }
                // After we processed the reaction, we return -> no other reaction should be done.
                return true;
            }
        }
        //Should not reach here.
        throw illegal_state("Pair reaction failed, no reaction selected when reaction was supposed to happen.");
    }

        
    const bool attempt_interaction(particle_id_pair const& pp, boost::shared_ptr<structure_type> const& surface)
    // This handles the 'reaction' (interaction) between a particle and a surface.
    // Returns True if the interaction was succesfull
    {
        // Get the reaction rules for the interaction with this surface.
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid(), surface->sid() ));
        if (::size(rules) == 0)
            throw propagation_error("trying to fire an interaction with a surface without a reaction rule.");
        
                
        const Real k_tot( k_total(pp.second.sid(), surface->sid()) );
        const Real rnd( k_tot * rng_() );
        Real k_cumm = 0;  

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& reaction_rule(*i); 
            k_cumm += reaction_rule.k();

            if (k_cumm >= rnd)
            // We have found the reaction rule that is in effect
            {
                // get the products
                const typename reaction_rule_type::species_id_range products(reaction_rule.get_products());

                switch (::size(products))
                {
                    // When no products
                    case 0:
                    {

                        tx_.remove_particle(pp.first);  // pp was just popped of the local queue so only has to be removed from the upstream particlecontainer.
                        break;
                    }
                    // When one product particle
                    case 1:
                    {
                        const species_type product_species(tx_.get_species(products[0]));

                        // Get the new position of the particle on the surface
                        const position_type new_pos( tx_.apply_boundary( surface->projected_point( pp.second.position() ).first ) );
                            
                        // check for overlap with particles that are in the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, product_species.radius()), pp.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }
                        // check for overlap with particles outside of the progator (e.g. the eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, product_species.radius()), 
                                    pp.first))
                            {
                                throw propagation_error("no space (vc)");
                            }
                        }
                        // TODO check overlap with surfaces.

                        // Process changes
                        tx_.remove_particle(pp.first);
                        // Make new particle in interaction surface 
                        const particle_id_pair product_particle( tx_.new_particle(product_species.id(), surface->id(), new_pos) );
                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    reaction_rule.id(), array_gen(product_particle.first), pp.first));
                        }
                        break;
                    }
                    // When number of products was not '0' or '1'
                    default:
                        throw not_implemented("surface interactions that produce more than one product are not supported");
                }
                // After we processed the reaction, we return -> no other reaction should be done.
                return true;
            }
        }
        //should not reach here.
        throw illegal_state("Surface interaction failed, no interaction selected when interaction was supposed to happen.");
    }
    

    const position_type make_move(species_type const& species, position_type const& old_pos, particle_id_type const& ignore) const
    // generates a move for a particle and checks if the move was allowed.
    {

        // FIXME this is just redoing what we are doing above when moving a particle.
        position_type     new_pos(old_pos);
        structure_id_type new_structure_id();

        if(species.D() == 0)
        {
            new_pos = old_pos;
//            new_structure_id = old_structure_id;
        }
        else
        {
            const boost::shared_ptr<const structure_type> s0_struct( tx_.get_structure( tx_.get_particle(ignore).second.structure_id() ) );

            const position_type displacement( s0_struct->bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_) );
            new_pos = tx_.apply_boundary( add( old_pos, displacement ) );
            // TODO apply deflection
                        
            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlapped( 
                                        tx_.check_overlap(
                                        particle_shape_type(new_pos, species.radius()),
                                        ignore));
                                                                
            // See if it overlaps with any particles
            if ( (overlapped && overlapped->size() > 0) )
                return old_pos;
            // See if it overlaps with any surfaces.
            if(species.structure_type_id() == tx_.get_def_structure_type_id())
            {
                const structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( new_pos, tx_.get_particle(ignore).second.structure_id() ) );
                const boost::shared_ptr<const structure_type> closest_struct( tx_.get_structure( struct_id_distance_pair.first ) );

                if( closest_struct->bounced( tx_.cyclic_transpose( old_pos, closest_struct->position() ), 
                                             tx_.cyclic_transpose( new_pos, closest_struct->position() ), 
                                             struct_id_distance_pair.second, species.radius() ) )    
                    return old_pos;        
            }
        }    
        return new_pos;
    }


    /* Given a position pair and two species, the function generates two new 
       positions based on two moves drawn from the free propagator. */
    const position_pair_type make_pair_move(species_type const& s0, species_type const& s1,
                                            position_pair_type const& pos0_pos1, particle_id_type const& ignore) const
    {
        const length_type min_distance_01( s0.radius() + s1.radius() );
        
        position_type temp_pos0(pos0_pos1.first),
                      temp_pos1(pos0_pos1.second);
        
        //Randomize which species moves first, and if there will be one or two moves. (50% of the cases)
        int  move_count( 1 + rng_.uniform_int(0, 1) );  // determine number of particles to move (1 or 2)
        bool move_first( rng_.uniform_int(0, 1) );      // Set to true if the first particle is moved, false means second particle.
        while( move_count-- )
        {
            if( move_first )
            // Move the first particle
            {
                temp_pos0 = make_move(s0, pos0_pos1.first, ignore);
                // check for overlap with the second particle (and potentially revert move)
                temp_pos0 = tx_.distance(temp_pos0, temp_pos1) - min_distance_01 < 0 ? pos0_pos1.first : temp_pos0;
            }                       
            else
            // Move the second particle
            {
                temp_pos1 = make_move(s1, pos0_pos1.second, ignore);
                // check for overlap with the first particle (and potentially revert move)
                temp_pos1 = tx_.distance(temp_pos1, temp_pos0) - min_distance_01 < 0 ? pos0_pos1.second : temp_pos1;
            }

            // Swap what particle to move
            move_first = !move_first;
        }  
            
        return std::make_pair(temp_pos0, temp_pos1);
    }


    /* Function returns the accumulated intrinsic reaction rate between to species. */
    /* Note that the species_id_type is equal to the structure_type_id_type for structure_types */
    const Real k_total(species_id_type const& s0_id, species_id_type const& s1_id) const
    {   
        Real k_tot (0);
        reaction_rules const& rules(rules_.query_reaction_rule(s0_id, s1_id));

//        if(::size(rules) != 0)
//        {
//            return 0;
//        }
//        else
//        {
        // This should also work when 'rules' is empty.
            for (typename boost::range_const_iterator<reaction_rules>::type
                    i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
            {
                k_tot += (*i).k();
            }
//        }
        return k_tot;
    }


    void remove_particle(particle_id_type const& pid)
    // Removes a particle from the propagator and particle container.
    // Note that this should only be called for the particles that are NOT processed currently (i.e. the particle that are just popped of the queue)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));

        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
        {
            queue_.erase(i);
        }
        else
        {
            throw not_found(std::string("Unknown particle (id=") + boost::lexical_cast<std::string>(pid) + ")");
        }
    }


////// Member variables
private:
    particle_container_type&    tx_;            // The 'ParticleContainer' object that contains the particles.
    network_rules_type const&   rules_;
    rng_type&                   rng_;
    const Real                  dt_;
    const int                   max_retry_count_;
    reaction_recorder_type* const rrec_;
    volume_clearer_type* const  vc_;
    particle_id_vector_type     queue_;
    int                         rejected_move_count_;
    const Real                  reaction_length_;
    static Logger&              log_;
};

template<typename Ttraits_>
Logger& newBDPropagator<Ttraits_>::log_(Logger::get_logger("ecell.newBDPropagator"));

#endif /* NEW_BD_PROPAGATOR_HPP */

