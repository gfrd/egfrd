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

    typedef typename particle_container_type::particle_id_pair_generator         particle_id_pair_generator;
    typedef typename particle_container_type::particle_id_pair_and_distance      particle_id_pair_and_distance;
    typedef typename particle_container_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
    typedef typename particle_container_type::structure_id_pair_and_distance      structure_id_pair_and_distance;
    typedef typename particle_container_type::structure_id_pair_and_distance_list structure_id_pair_and_distance_list;
    typedef typename traits_type::world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type                         time_type;
    typedef typename traits_type::network_rules_type                network_rules_type;
    typedef typename network_rules_type::reaction_rules             reaction_rules;
    typedef typename network_rules_type::reaction_rule_type         reaction_rule_type;
    typedef typename traits_type::reaction_record_type              reaction_record_type;
    typedef typename traits_type::reaction_recorder_type            reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type               volume_clearer_type;

    typedef std::vector<particle_id_type>                                       particle_id_vector_type;
    typedef std::pair<position_type, position_type>                             position_pair_type;
    typedef std::pair<position_type, structure_id_type>                         position_structid_pair_type;
    typedef std::pair<position_structid_pair_type, position_structid_pair_type> posstructid_posstructid_pair_type;
    typedef std::pair<length_type, length_type>                                 component_pair_type;
    typedef std::pair<position_type, component_pair_type>                       projected_type;

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
        const species_type                      pp_species(tx_.get_species(pp.second.sid()));
        const length_type                       r0( pp_species.radius() );
        const position_type                     old_pos( pp.second.position() );
        const structure_id_type                 old_struct_id(pp.second.structure_id());
        const boost::shared_ptr<structure_type> pp_structure( tx_.get_structure( old_struct_id ) );
        // declare variables to use
        position_type                   new_pos(old_pos);
        structure_id_type               new_structure_id( old_struct_id );

        //// 2. New position
        /* Sample a potential move, and check if the particle _core_ has overlapped with another 
           particle or surface. If this is the case the particle bounces, and is returned
           to it's original position. For a particle with D = 0, new_pos = old_pos */        
        if (pp_species.D() != 0.0)
        {
            // get a new position, dependent on the structure the particle lives on.
            const position_type displacement( pp_structure->bd_displacement(pp_species.v() * dt_,
                                                                            std::sqrt(2.0 * pp_species.D() * dt_),
                                                                            rng_) );

            const position_structid_pair_type pos_structid(tx_.apply_boundary( std::make_pair( add(old_pos, displacement), old_struct_id )));
            new_pos = pos_structid.first;
            new_structure_id = pos_structid.second;
        }

        bool bounced( false );
        // Get all the particles that are in the reaction volume at the new position (and that may consequently also be in the core)
        /* Use a spherical shape with radius = particle_radius + reaction_length.
	       We use it to check for overlaps with other particles. */
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlap_particles(
                tx_.check_overlap(particle_shape_type( new_pos, r0 + reaction_length_ ), pp.first));
        int particles_in_overlap(overlap_particles ? overlap_particles->size(): 0);

        //// 3.1 Check for core overlaps with particles
        /* Check if the particle at new_pos overlaps with any particle cores. */        
        int j( 0 );
        while(!bounced && j < particles_in_overlap)
            bounced = overlap_particles->at(j++).second < r0;
        
        //// 3.2 Check for core overlaps structures
        /* If the particle has not bounced with another particle check for overlap with a surface. */
        boost::scoped_ptr<structure_id_pair_and_distance_list> overlap_structures;
        int structures_in_overlap(0);
        if (!bounced)
        {
            // Get all the structures within the reaction volume (and that may consequently also be in the core)
            boost::scoped_ptr<structure_id_pair_and_distance_list> overlap_structures_tmp(
                tx_.check_surface_overlap(particle_shape_type( new_pos, r0 + reaction_length_ ), old_pos, new_structure_id, r0, old_struct_id));
            overlap_structures.swap(overlap_structures_tmp);
            int structures_in_overlap(overlap_structures ? overlap_structures->size(): 0);
 
            j = 0;
            while(!bounced && j < structures_in_overlap)
                bounced = overlap_structures->at(j++).second < r0;
        }
         
        //// 4. Finalize the position
        /* If particle is bounced, restore old position and check reaction_volume for reaction partners at old_pos. */
        if(bounced)
        {
            // restore old position and structure_id
            new_pos = old_pos;
            new_structure_id = old_struct_id;
            
            // re-get the reaction partners (particles), now on old position.
            boost::scoped_ptr<particle_id_pair_and_distance_list> overlap_particles_after_bounce( 
                    tx_.check_overlap( particle_shape_type( old_pos, r0 + reaction_length_ ), pp.first) );
            overlap_particles.swap( overlap_particles_after_bounce );     // FIXME is there no better way?
            // NOTE that it is asserted that the particle overlap criterium for the particle with other particles
            // and surfaces is False!
            particles_in_overlap = overlap_particles ? overlap_particles->size(): 0; 
            
            // re-get the reaction partners (structures), now on old position.
            boost::scoped_ptr<structure_id_pair_and_distance_list> overlap_structures_after_bounce( 
                    tx_.check_surface_overlap( particle_shape_type( old_pos, r0 + reaction_length_ ), old_pos, old_struct_id, r0) );
            overlap_structures.swap( overlap_structures_after_bounce );     // FIXME is there no better way?
            // NOTE that it is asserted that the particle overlap criterium for the particle with other particles
            // and surfaces is False!
            structures_in_overlap = overlap_structures ? overlap_structures->size(): 0; 
        }
        

        //// 5.
        /* Attempt a reaction (and/or interaction) with all the particles (and/or a surface) that 
           are/is inside the reaction volume. */
        Real accumulated_prob (0);
        Real rnd( rng_() );
        
        const boost::shared_ptr<const structure_type> current_struct( tx_.get_structure( new_structure_id) );
        /* First, if a surface is inside the reaction volume, and the particle is in the 3D attempt an interaction. */
        // TODO also interaction should be allowed when a particle is on a cylinder.
        // TODO don't check only the closest but check all overlapping surfaces.
        j = 0;
        while(j < structures_in_overlap)
        {
            const structure_id_pair_and_distance & overlap_struct( overlap_structures->at(j) );
            const projected_type new_pos_projected (overlap_struct.first.second->project_point(new_pos));

            // For an interaction the structure must be:
            if( (overlap_struct.first.second->sid() != current_struct->sid()) &&    // - of a different structure_type
                (overlap_struct.second < r0 + reaction_length_) &&                  // - within the reaction volume
                (new_pos_projected.second.second < 0) )                             // - 'alongside' of the particle
            {
                accumulated_prob += k_total( pp.second.sid(), overlap_struct.first.second->sid() ) * dt_ / 
                                    overlap_struct.first.second->surface_reaction_volume( r0, reaction_length_ );
        
                if(accumulated_prob >= 1.)
                {
                    LOG_WARNING(("the acceptance probability of an interaction/reaction exeeded one; %f.",
                                 accumulated_prob));
                } 
                
                if( accumulated_prob > rnd )
                {            
                    try
                    {
                        LOG_DEBUG(("fire surface interaction"));
                        if(attempt_interaction(pp, new_pos_projected.first, overlap_struct.first.second ))
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
                               boost::lexical_cast<std::string>(overlap_struct.first.first).c_str()));
                }
            }

            j++;
        }
        
        ////
        /* Now attempt a reaction with all particles inside the reaction volume. */
        j = 0;
        while(j < particles_in_overlap)
        {
            const particle_id_pair_and_distance & overlap_particle( overlap_particles->at(j) );

            species_type s0(pp_species);
            species_type s1(tx_.get_species(overlap_particle.first.second.sid()));
                
/*            // If the structure_types of the reactants are not equal, one of the reactants has to come from the bulk,
            // and we let this be s1, the particle from the surface is named s0.
            if(s0.structure_type_id() != s1.structure_type_id())
            {
                if( !(s0.structure_type_id() == tx_.get_def_structure_type_id() || s1.structure_type_id() == tx_.get_def_structure_type_id())  )
                {
                    log_.warn("A surface particle overlapped with a particle living on a different structure_type. No reaction implemented.");
                    continue;           // FIXME this now blocks the cap dissociation.
                }
                if(s0.structure_type_id() == tx_.get_def_structure_type_id())
                    std::swap(s0,s1);       // ??????
            }
*/            
            const boost::shared_ptr<const structure_type> s1_struct( tx_.get_structure( overlap_particle.first.second.structure_id()) );
            accumulated_prob += k_total(s0.id(), s1.id()) * dt_ /
                                ( 2. * s1_struct->particle_reaction_volume( s0.radius() + s1.radius(), reaction_length_ ) ); 
            
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
                    if( attempt_pair_reaction(pp, overlap_particle.first, s0, s1) )
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
                    boost::lexical_cast<std::string>(overlap_particle.first.first).c_str()));
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
            reaction_rule_type const& reaction_rule(*i);
            k_cumm += reaction_rule.k();

            if (k_cumm > rnd)
            // We have found the reaction rule that is in effect
            {
                // get the products
                typename reaction_rule_type::species_id_range products(reaction_rule.get_products());

                switch (::size(products))
                {
                    // When no product particles
                    case 0:
                    {
                        remove_particle(pp.first);  // pp was just popped of the local queue so only has to be removed from the upstream particlecontainer.
                        break;
                    }

                    // When one product particle
                    case 1:
                    {
                        const position_type                             reactant_pos(pp.second.position());
                        const boost::shared_ptr<const structure_type>   reactant_structure( tx_.get_structure(pp.second.structure_id()) );
                        const species_type                              product_species(tx_.get_species(products[0]));

                        // default values for the position and structure of the product
                        position_structid_pair_type     product_pos_struct_id (std::make_pair(reactant_pos, pp.second.structure_id()));
                        
                        /* If the product particle does NOT live on the same structure type as the reactant => structure -> superstructure dissociation. */
                        if( product_species.structure_type_id() != reactant_structure->sid() )
                        {
                            // get dissociation position and structure
                            // TODO merge the structure_id thing into the 'surface_dissociation_vector' method.
                            // TODO when dissociating from the 'cap' structure, we go up two levels instead of one.
                            const position_type displacement( reactant_structure->surface_dissociation_vector(rng_, product_species.radius(), reaction_length_ ) );
                            const position_type product_pos (tx_.apply_boundary( add( reactant_pos, displacement ) ));
                            // When dissociating, the new structure is the parent structure (is this always true?)
                            const structure_id_type product_structure_id (reactant_structure->structure_id());

                            product_pos_struct_id = std::make_pair(product_pos, product_structure_id);
                            // Particle is allowed to move after dissociation from surface. TODO Isn't it allways allowed to move?
                            product_pos_struct_id = make_move(product_species, product_pos_struct_id, pp.first);
                        }


                        ////// Check for overlaps
                        const particle_shape_type new_shape(product_pos_struct_id.first, product_species.radius());
                        // check for overlap with particles that are inside the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles(
                                tx_.check_overlap(new_shape, pp.first));
                        if (overlap_particles && overlap_particles->size() > 0)
                        {
                            throw propagation_error("no space due to other particle");
                        }
                        // check overlap with surfaces
                        const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_surfaces(
                                tx_.check_surface_overlap(new_shape, reactant_pos, product_pos_struct_id.second, product_species.radius()));
                        if (overlap_surfaces && overlap_surfaces->size() > 0)
                        {
                            throw propagation_error("no space due to near surface");
                        }
                        // check for overlap with particles outside of the propagator (for example the eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(new_shape, pp.first))
                            {
                                throw propagation_error("no space due to other particle (vc)");
                            }
                        }


                        // Process changes
                        remove_particle(pp.first);
                        // Make new particle on right surface
                        const particle_id_pair product_particle( tx_.new_particle(product_species.id(), product_pos_struct_id.second, product_pos_struct_id.first) );
                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    reaction_rule.id(), array_gen(product_particle.first), pp.first));
                        }
                        break;
                    }

                    // When two product particles
                    case 2:
                    {
                        const position_type                             reactant_pos(pp.second.position());
                        const boost::shared_ptr<const structure_type>   reactant_structure( tx_.get_structure( pp.second.structure_id() ) );
                        const boost::shared_ptr<const structure_type>   parent_structure( tx_.get_structure( reactant_structure->structure_id() ) );
                        
                        // We set: species 0 lives on the current structure; species 1 lives in the parent structure (if applicable).
                        species_type product0_species(tx_.get_species(products[0]));
                        species_type product1_species(tx_.get_species(products[1]));
                        if( product0_species.structure_type_id() != product1_species.structure_type_id() )
                        {
                            if( !(product0_species.structure_type_id() == parent_structure->sid() || 
                                  product1_species.structure_type_id() == parent_structure->sid())  )
                                throw not_implemented("One of product species should live on parent structure.");
                                             
                            if ( !(product0_species.structure_type_id() == reactant_structure->sid()) )
                                std::swap(product0_species, product1_species);
                        }

                        position_pair_type pos0pos1;
                        structure_id_type prod0_struct_id;
                        structure_id_type prod1_struct_id;
                        posstructid_posstructid_pair_type pos0pos1_pair;

                        /* Create positions (np0, np1) for the reaction products. 
                        The iv between the products lies within the reaction volume. */
                        int i (max_retry_count_);
                        while (true)
                        {
                            // Escape clause #1: Too many tries
                            if (i-- == 0)
                            {
                                throw propagation_error("no space due to particle or surface");
                            }
                            

                            //// Generate positions and structure_id for product particles after reaction

                            /* If the product particles do NOT live on the same structure => structure -> structure + parent_structure dissociation.
                               Else, the products live on the same structure AND on the same structure as the reactant. 
                               => standard geminate dissociation. */
                            if( product0_species.structure_type_id() != product1_species.structure_type_id() )
                            {
                                // TODO include structure_id in return_type of dissociation_positions
                                pos0pos1 = reactant_structure-> special_geminate_dissociation_positions( rng_, product0_species, product1_species,
                                                                                                         reactant_pos, reaction_length_ );
                                prod0_struct_id = pp.second.structure_id();             // The 'original' structure
                                prod1_struct_id = reactant_structure->structure_id();   // The parent structure
                            }
                            else
                            {
                                // TODO include structure_id in return_type of dissociation_positions
                                pos0pos1 = reactant_structure-> geminate_dissociation_positions( rng_, product0_species, product1_species,
                                                                                                 reactant_pos, reaction_length_ );
                                prod0_struct_id = pp.second.structure_id();         // same structure as the reactant
                                prod1_struct_id = pp.second.structure_id();         // ditto
                            }
                            // TODO clean this up
                            const boost::shared_ptr<const structure_type> prod0_structure(reactant_structure);
                            const boost::shared_ptr<const structure_type> prod1_structure(tx_.get_structure(prod1_struct_id));
                            const position_structid_pair_type pos_structid0 (tx_.apply_boundary( std::make_pair(pos0pos1.second, prod1_struct_id) ));
                            const position_structid_pair_type pos_structid1 (tx_.apply_boundary( std::make_pair(pos0pos1.first,  prod0_struct_id) ));
                            pos0pos1_pair = std::make_pair(pos_structid0, pos_structid1);


                            ////// Check for overlaps
                            const particle_shape_type new_shape0(pos0pos1_pair.first.first, product0_species.radius());
                            // check overlap with particle in the propagator
                            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles_product0(
                                    tx_.check_overlap( new_shape0,  pp.first));
                            if (overlap_particles_product0 && overlap_particles_product0->size() > 0)
                            {
                                continue;   // try other positions pair
                            }
                            /* Only when both products live in the bulk, check for possible overlap with a surface. */
                            const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_structures_product0(
                                    tx_.check_surface_overlap(new_shape0, reactant_pos, pos0pos1_pair.first.second, product0_species.radius()));
                            if (overlap_structures_product0 && overlap_structures_product0->size() > 0)
                            {
                                continue;   // try other positions pair
                            }
                            const particle_shape_type new_shape1(pos0pos1_pair.second.first, product1_species.radius());
                            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles_product1(
                                    tx_.check_overlap( new_shape1, pp.first));
                            if (overlap_particles_product1 && overlap_particles_product1->size() > 0)
                            {
                                continue;   // try other positions pair
                            }
                            const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_structures_product1(
                                    tx_.check_surface_overlap(new_shape1, reactant_pos, pos0pos1_pair.second.second, product1_species.radius()));
                            if (overlap_structures_product1 && overlap_structures_product1->size() > 0)
                            {
                                continue;   // try other positions pair
                            }
                            // If no overlap occured -> positions are ok
                            break;
                        }

                        /* After dissociation both particles are allowed to move. */
                        pos0pos1_pair = make_pair_move(product0_species, product1_species, pos0pos1_pair, pp.first);

                        // check overlap with particle outside of propagator (for example in eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(pos0pos1_pair.first.first,  product0_species.radius()), pp.first) ||
                                !(*vc_)(particle_shape_type(pos0pos1_pair.second.first, product1_species.radius()), pp.first))
                            {
                                throw propagation_error("no space due to near particle (vc)");
                            }
                        }

                        // Process changes
                        remove_particle(pp.first);
                        const particle_id_pair product0_particle(tx_.new_particle(product0_species.id(), pos0pos1_pair.first.second,  pos0pos1_pair.first.first));
                        const particle_id_pair product1_particle(tx_.new_particle(product1_species.id(), pos0pos1_pair.second.second, pos0pos1_pair.second.first));

                        // Record changes
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    reaction_rule.id(),
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
                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        break;
                    }
                    // When one product particle
                    case 1:
                    {
                        const structure_id_type             reactant0_structure_id (pp0.second.structure_id());
                        const structure_id_type             reactant1_structure_id (pp1.second.structure_id());
                        boost::shared_ptr<structure_type>   reactant0_structure (tx_.get_structure( reactant0_structure_id ));
                        boost::shared_ptr<structure_type>   reactant1_structure (tx_.get_structure( reactant1_structure_id ));
                        position_type                       reactant0_pos (pp0.second.position());
                        position_type                       reactant1_pos (pp1.second.position());

                        const species_type                  product_species(tx_.get_species(products[0]));
                        structure_id_type                   product_structure_id(reactant0_structure_id);

                        // If the particles lived on adjacent structures of the same type
                        if ((s0.structure_type_id() == s1.structure_type_id()) &&
                            (reactant0_structure_id != reactant1_structure_id))
                        {
                            const position_structid_pair_type pos_cyclic1 (tx_.cyclic_transpose(std::make_pair(reactant1_pos, reactant1_structure_id),
                                                                                                (*reactant0_structure) ));
                            reactant1_pos = pos_cyclic1.first;
                        }

                        // calculate the product position (CoM of the two particles)
                        position_type product_pos( tx_.apply_boundary(
                                                    divide(
                                                    add(multiply(reactant0_pos, s1.D()),
                                                        multiply(tx_.cyclic_transpose(
                                                            reactant1_pos,
                                                            reactant0_pos), s0.D())),
                                                    (s0.D() + s1.D()))) );
//                        position_type product_pos (tx_.calculate_pair_CoM(pp0.second.position(), pp1.second.position()), s0.D(), s1.D());

                        // For unequal structure types, project product_pos on the 'lesser' structure.
                        // The 'lesser' structure is the structure that is the lowest in the hierarchy of the two.
                        // Note that this means that a reaction can only cross ONE hierarchical level of structures.
                        if( s0.structure_type_id() != s1.structure_type_id() )
                        {

                            boost::shared_ptr<structure_type> product_structure;
                            if ( reactant0_structure->structure_id() == reactant1_structure_id )        // if pp1 lives on the parent structure of pp0
                            {
                                product_structure_id = reactant0_structure_id;                          // -> lesser structure is structure0
                                product_structure = reactant0_structure;
                            }
                            else if ( reactant1_structure->structure_id() == reactant0_structure_id )   // pp0 lives on the parent structure of pp1
                            {
                                product_structure_id = reactant1_structure_id;                          // -> lesser structure is structure1
                                product_structure = reactant1_structure;
                            }
                            else
                            {
                                throw propagation_error("Particles can only be one hierarchical level apart.");
                            }

                            projected_type const position_on_surface( product_structure-> project_point(
                                    tx_.cyclic_transpose( product_pos, product_structure->position() )));       // TODO check of cyclic transpose is still needed.

                            // check that the projected point is not outside of the product surface.
                            if ( position_on_surface.second.second > 0 )
                                throw propagation_error("position product particle was not in surface.");

                            product_pos = tx_.apply_boundary( position_on_surface.first );
                        }
                        else
                        // The structure_types are the same
                        {
                            if (reactant0_structure_id != reactant1_structure_id)
                            {
                                const position_structid_pair_type bla (tx_.apply_boundary(std::make_pair(product_pos, reactant0_structure_id)));
                                // deflect the product coordinate and get right structure_id.
                                product_pos = bla.first;
                                product_structure_id = bla.second;

                            }
                        }


                        //// Check for overlaps
                        const particle_shape_type new_shape(product_pos, product_species.radius());
                        // check for overlap with particles that are inside the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles(
                                tx_.check_overlap(new_shape, pp0.first, pp1.first));
                        if( overlap_particles && overlap_particles->size() > 0 )
                        {
                            throw propagation_error("no space due to particle");
                        }
                        // Check overlap with surfaces
                        const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_structures(
                                tx_.check_surface_overlap(new_shape, product_pos, product_structure_id, product_species.radius()));
                        if (overlap_structures && overlap_structures->size() > 0)
                        {
                            throw propagation_error("no space due to near surface");
                        }
                        // check for overlap with particles outside of the propagator (for example the eGFRD simulator)
                        if( vc_ )
                        {
                            if (!(*vc_)(new_shape, pp0.first, pp1.first))
                            {
                                throw propagation_error("no space due to particle (vc)");
                            }
                        }


                        // Process changes
                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        // Make new particle on right structure
                        particle_id_pair product_particle(tx_.new_particle(product_species.id(), product_structure_id, product_pos));
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
        // Should not reach here.
        throw illegal_state("Pair reaction failed, no reaction selected when reaction was supposed to happen.");
    }

        
    const bool attempt_interaction(particle_id_pair const& pp, position_type const& pos_in_struct, boost::shared_ptr<structure_type> const& structure)
    // This handles the 'reaction' (interaction) between a particle and a structure.
    // Returns True if the interaction was succesfull
    {
        // Get the reaction rules for the interaction with this structure.
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid(), structure->sid() ));
        if (::size(rules) == 0)
            throw propagation_error("trying to fire an interaction with a structure without a reaction rule.");
        
                
        const Real k_tot( k_total(pp.second.sid(), structure->sid()) );
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

                        remove_particle(pp.first);
                        break;
                    }
                    // When one product particle
                    case 1:
                    {
                        const species_type product_species(tx_.get_species(products[0]));

                        // Get the new position of the particle on the structure
                        const position_type product_pos( tx_.apply_boundary(pos_in_struct) );


                        ///// Check for overlaps   
                        const particle_shape_type new_shape(product_pos, product_species.radius());
                        // check for overlap with particles that are in the propagator
                        const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles(
                            tx_.check_overlap(new_shape, pp.first));
                        if (overlap_particles && overlap_particles->size() > 0)
                        {
                            throw propagation_error("no space");
                        }
                        // check overlap with structures.
                        // No need to include old structure as ignore because going down the structure hierarchy anyway.
                        const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_structures(
                                tx_.check_surface_overlap(new_shape, product_pos, structure->id(), product_species.radius()));
                        if (overlap_structures && overlap_structures->size() > 0)
                        {
                            throw propagation_error("no space due to near surface");
                        }
                        // check for overlap with particles outside of the propagator (e.g. the eGFRD simulator)
                        if (vc_)
                        {
                            if (!(*vc_)(new_shape, pp.first))
                            {
                                throw propagation_error("no space (vc)");
                            }
                        }

                        // Process changes
                        remove_particle(pp.first);
                        // Make new particle in interaction structure 
                        const particle_id_pair product_particle( tx_.new_particle(product_species.id(), structure->id(), product_pos) );
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
                        throw not_implemented("structure interactions that produce more than one product are not supported");
                }
                // After we processed the reaction, we return -> no other reaction should be done.
                return true;
            }
        }
        //should not reach here.
        throw illegal_state("Surface interaction failed, no interaction selected when interaction was supposed to happen.");
    }
    

    const position_structid_pair_type make_move(species_type const& species, position_structid_pair_type const& old_pos_struct_id,
                                                particle_id_type const& ignore) const
    // generates a move for a particle and checks if the move was allowed.
    {

        // FIXME this is just redoing what we are doing above when moving a particle.
        position_type     new_pos(old_pos_struct_id.first);
        structure_id_type new_structure_id(old_pos_struct_id.second);

        if(species.D() != 0)
        {
            const boost::shared_ptr<structure_type> old_structure( tx_.get_structure( new_structure_id ) );

            const position_type displacement( old_structure->bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_) );

            const position_structid_pair_type pos_structid(tx_.apply_boundary( std::make_pair( add(old_pos_struct_id.first, displacement), old_pos_struct_id.second) ));
            new_pos = pos_structid.first;
            new_structure_id = pos_structid.second;

            // See if it overlaps with any particles
            const particle_shape_type new_shape(new_pos, species.radius());
            const boost::scoped_ptr<const particle_id_pair_and_distance_list> overlap_particles( 
                    tx_.check_overlap( new_shape, ignore));
            if (overlap_particles && overlap_particles->size() > 0)
                return old_pos_struct_id;

            // See if it overlaps with any surfaces.
            const boost::scoped_ptr<const structure_id_pair_and_distance_list> overlap_surfaces(
                    tx_.check_surface_overlap(new_shape, old_pos_struct_id.first, new_structure_id, species.radius())); // TODO include to ignore old_structure
            if (overlap_surfaces && overlap_surfaces->size() > 0)
                return old_pos_struct_id;
        }    

        return std::make_pair(new_pos, new_structure_id);
    }


    /* Given a position pair and two species, the function generates two new 
       positions based on two moves drawn from the free propagator. */
    const posstructid_posstructid_pair_type make_pair_move(species_type const& s0, species_type const& s1,
                                                           posstructid_posstructid_pair_type const& posstruct0_posstruct1,
                                                           particle_id_type const& ignore) const
    {
        const length_type min_distance_01( s0.radius() + s1.radius() );
        
        position_structid_pair_type temp_posstruct0(posstruct0_posstruct1.first),
                                    temp_posstruct1(posstruct0_posstruct1.second);
        
        //Randomize which species moves first, and if there will be one or two moves. (50% of the cases)
        int  move_count( 1 + rng_.uniform_int(0, 1) );  // determine number of particles to move (1 or 2)
        bool move_first( rng_.uniform_int(0, 1) );      // Set to true if the first particle is moved, false means second particle.
        while( move_count-- )
        {
            if( move_first )
            // Move the first particle
            {
                temp_posstruct0 = make_move(s0, posstruct0_posstruct1.first, ignore);
                // check for overlap with the second particle and potentially revert move
                if (tx_.distance(temp_posstruct0.first, temp_posstruct1.first) < min_distance_01)
                    temp_posstruct0 = posstruct0_posstruct1.first;
            }                       
            else
            // Move the second particle
            {
                temp_posstruct1 = make_move(s1, posstruct0_posstruct1.second, ignore);
                // check for overlap with the first particle and potentially revert move
                if (tx_.distance(temp_posstruct1.first, temp_posstruct0.first) < min_distance_01)
                    temp_posstruct1 = posstruct0_posstruct1.second;
            }

            // Swap what particle to move
            move_first = !move_first;
        }  
            
        return std::make_pair(temp_posstruct0, temp_posstruct1);
    }


    /* Function returns the accumulated intrinsic reaction rate between to species. */
    /* Note that the species_id_type is equal to the structure_type_id_type for structure_types */
    const Real k_total(species_id_type const& s0_id, species_id_type const& s1_id) const
    {   
        Real k_tot (0);
        reaction_rules const& rules(rules_.query_reaction_rule(s0_id, s1_id));

        // This should also work when 'rules' is empty.
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            k_tot += (*i).k();
        }
        return k_tot;
    }


    void remove_particle(particle_id_type const& pid)
    // Removes a particle from the propagator and particle container.
    // Note that if the particle is no longer in the local list (we just popped it off and treating it currently), it will
    // only be removed in the particle container (tx_)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));

        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
        {
            queue_.erase(i);
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

