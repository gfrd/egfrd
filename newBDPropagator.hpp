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
    typedef typename Ttraits_::world_type::particle_container_type particle_container_type;
    typedef typename particle_container_type::species_id_type species_id_type;
    typedef typename particle_container_type::position_type position_type;
    typedef typename particle_container_type::particle_shape_type particle_shape_type;
    typedef typename particle_container_type::species_type species_type;
    typedef typename particle_container_type::length_type length_type;
    typedef typename particle_container_type::particle_id_type particle_id_type;
    typedef typename particle_container_type::particle_type particle_type;
    typedef typename particle_container_type::particle_id_pair particle_id_pair;
    typedef typename particle_container_type::structure_id_and_distance_pair structure_id_and_distance_pair;
    typedef std::vector<particle_id_type> particle_id_vector_type;
    typedef std::pair<position_type, position_type> position_pair_type;
    typedef std::pair<position_type, length_type> projected_type;
    typedef typename particle_container_type::particle_id_pair_generator particle_id_pair_generator;
    typedef typename particle_container_type::particle_id_pair_and_distance particle_id_pair_and_distance;
    typedef typename particle_container_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
    typedef typename particle_container_type::structure_type structure_type;
    typedef typename particle_container_type::structure_id_type structure_id_type;    
    typedef typename traits_type::world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type volume_clearer_type;


public:
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
        for (typename boost::range_const_iterator<Trange_>::type
                i(boost::begin(particle_ids)),
                e(boost::end(particle_ids)); i != e; ++i)
        {
            queue_.push_back(*i);
        }
        shuffle(rng, queue_);
    }

    bool operator()()
    {
        if (queue_.empty())
            return false;        

        particle_id_type pid(queue_.back());
        queue_.pop_back();
        particle_id_pair pp(tx_.get_particle(pid));

        LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));
        
        /* Consider decay reactions first. Particles can move after decay, but this is done
           inside the 'attempt_reaction' function. */
        try
        {
            if (attempt_reaction(pp))
                return true;
        }
        catch (propagation_error const& reason)
        {
            log_.info("first-order reaction rejected (reason: %s)", reason.what());
            ++rejected_move_count_;
            return true;
        }
        
        species_type const pp_species(tx_.get_species(pp.second.sid()));
        length_type const r0( pp_species.radius() );
        boost::shared_ptr<structure_type> const pp_structure( tx_.get_structure( pp_species.structure_id() ) );
        
        structure_id_and_distance_pair struct_id_distance_pair;
        length_type particle_surface_distance = 0;
        position_type new_pos;
        
        /* Sample a potential move, and check if the particle _core_ has overlapped with another 
           particle or surface. If this is the case the particlem bounces, and is returned
           to it's original position. For a particle with D = 0, new_pos = old_pos */        
        if (pp_species.D() != 0.)
        {
            // Make a step, dependent on the surface the particle lives on.
            position_type const displacement( pp_structure->
                                            bd_displacement(pp_species.v() * dt_, std::sqrt(2.0 * pp_species.D() * dt_), rng_) );

            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
        }
        else
        {
            new_pos = pp.second.position();
        }
        
        /* Use a spherical shape with radius = particle_radius + reaction_length.
	       We use it to check for overlaps with other particels. */
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                tx_.check_overlap(particle_shape_type( new_pos, r0 + reaction_length_ ), pp.first));
        int particles_in_overlap(overlapped ? overlapped->size(): 0);

        /* Check if the particle at new_pos overlaps with any particle cores. */        
        bool bounced = false;
        int j = 0;
        while(!bounced && j < particles_in_overlap)
            bounced = overlapped->at(j++).second < r0;
        
        /* If the particle has not bounced with another particle and lives in the bulk: 
           check for core overlap with a surface. */
        if(!bounced && pp_structure->id() == "world")
        {
            /* Check for overlap with nearest surface. */
            struct_id_distance_pair = tx_.get_closest_surface( new_pos, pp_structure->id() );
            particle_surface_distance = struct_id_distance_pair.second;
        
            bounced = determine_surface_passage(new_pos, pp.second.position(), struct_id_distance_pair, r0);
        }
         
        /* If particle is bounced, check reaction_volume for reaction partners at old_pos. */     
        if(bounced)
        {
            new_pos = pp.second.position();
            
            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_after_bounce( 
                    tx_.check_overlap( particle_shape_type( new_pos, r0 + reaction_length_ ), 
                            pp.first) );
                            
            overlapped.swap( overlapped_after_bounce );
                                            
            particles_in_overlap = overlapped ? overlapped->size(): 0; 
            
            if(pp_structure->id() == "world")
            {
                struct_id_distance_pair = tx_.get_closest_surface( new_pos, pp_structure->id() );
                particle_surface_distance = struct_id_distance_pair.second;  
            }
        }
        
        /* Attempt a reaction (and/or interaction) with all the particles (and/or a surface) that 
           are/is inside the reaction volume. */
        Real accumulated_prob = 0;
        Real const rnd( rng_() );
        
        /* First, if a surface is inside the reaction volume, attempt an interaction. */
        if( particle_surface_distance < reaction_length_ && pp_structure->id() == "world")
        {        
            boost::shared_ptr<structure_type> const closest_surf( tx_.get_structure(struct_id_distance_pair.first) );
            
            accumulated_prob += k_total( pp.second.sid(), closest_surf->sid() ) * dt_ / 
                closest_surf->surface_reaction_volume( r0, reaction_length_ );
        
            if(accumulated_prob >= 1.)
            {
                LOG_WARNING((
                    "the acceptance probability of an interaction exeeded one; %f.",
                    accumulated_prob));
            } 
        
            if( accumulated_prob > rnd )
            {            
                try
                {
                    LOG_DEBUG(("fire surface interaction"));
                    if(fire_interaction(pp, closest_surf ))
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
                    boost::lexical_cast<std::string>(closest_surf->id()).c_str()));
                ++rejected_move_count_;
            }
        }
        
        /* Now attempt a reaction with all particles inside the reaction volume. */
        j = 0;
        while(j < particles_in_overlap)
        {
            particle_id_pair_and_distance const& closest( overlapped->at(j) );

            species_type s0(pp_species),s1(tx_.get_species(closest.first.second.sid()));
                
            /* If the structures of the reactants do not equal, one of the reactants has to come from the bulk,
               and we let this be s1, the particle from the surface is named s0. */
            if(s0.structure_id() != s1.structure_id())
            {
                if( !(s0.structure_id() == "world" || s1.structure_id() == "world")  )
                {
                    log_.warn("A surface particle overlapped with a particle living on a different surface. No reaction implemented.");
                    continue;
                }
                if(s0.structure_id() == "world")
                    std::swap(s0,s1);
            }
            
            boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id()) );
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
                    if( fire_reaction(pp, closest.first) )
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
                LOG_DEBUG(("Particle attempted a reaction with a non-reactive particle %s.", 
                    boost::lexical_cast<std::string>(closest.first.first).c_str()));
                ++rejected_move_count_;
            }
            
            j++;
        }
                       
        /* If the particle did neither react, interact or bounce, update it to it's new position. */
        if(!bounced)
        {   
            particle_id_pair particle_to_update( pp.first, 
                        particle_type(pp_species.id(), particle_shape_type(new_pos, r0), pp_species.D()) );
                
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

    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

private:

    bool attempt_reaction(particle_id_pair const& pp)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const Real rnd(rng_() / dt_);
        Real prob = 0.;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            prob += r.k();
            if (prob > rnd)
            {
                typename reaction_rule_type::species_id_range products(
                        r.get_products());
                switch (::size(products))
                {
                case 0:
                    tx_.remove_particle(pp.first);
                    break;

                case 1:
                    {
                        const species_type s0(tx_.get_species(products[0]));
                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species(pp.second.sid()).structure_id() ) );
                        position_type new_pos = pp.second.position();
                        
                        /* If the product particle does NOT live of the same structure as the reactant => surface -> bulk dissociation. */
                        if( s0.structure_id() != pp_struct->id() )
                        {
                            if( s0.structure_id() != "world" )
                                throw not_implemented("No surface -> surface dissociation allowed");
                        
                            position_type const displacement( pp_struct->surface_dissociation_vector(rng_, s0.radius(), reaction_length_ ) );

                            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
                            
                            //Particle is allowed to move after dissociation from surface.
                            new_pos = make_move(s0, new_pos, pp.first);
                        }
                      
                        const particle_id_pair new_p(
                            pp.first, particle_type(products[0],
                                particle_shape_type(new_pos,
                                                    s0.radius()),
                                                    s0.D()));
                                                    
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(tx_.check_overlap(new_p.second.shape(), new_p.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space due to other particle");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(new_p.second.shape(), pp.first))
                            {
                                throw propagation_error("no space due to other particle (vc)");
                            }
                        }

                        if( s0.structure_id() == "world" && pp_struct->id() == "world" )
                        {
                            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos, s0.structure_id()) );
                            
                            if( determine_surface_passage(new_pos, pp.second.position(), struct_id_distance_pair, s0.radius()) )
                                throw propagation_error("no space due to near surface");
                        }

                        tx_.update_particle(new_p);

                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(new_p.first), pp.first));
                        }
                    }
                    break;

                case 2:
                    {
                        species_type s0(tx_.get_species(products[0])),
                                s1(tx_.get_species(products[1]));
                        int i = max_retry_count_;
                        position_pair_type pp01;

                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species( pp.second.sid() ).structure_id() ) );
                        
                        /* Create positions (np0, np1) for the reaction products. 
                        The ipv between the products lies within the reation volume. */
                        for (;;)
                        {
                            if (--i < 0)
                            {
                                throw propagation_error("no space due to particle or surface");
                            }
                            
                            /* If the product particles do NOT live on the same structure => surface -> surface + bulk dissociation.
                               Else, the products live on the same structure AND on the same structure as the reactant. 
                               => standard geminate dissociation. 
                               
                               We set: species 0 lives on the surface; species 1 lives in the bulk. */
                            if( s0.structure_id() != s1.structure_id() )
                            {
                                if( !(s0.structure_id() == "world" || s1.structure_id() == "world")  )
                                    throw not_implemented("No surface -> surface + surface - dissociation allowed");
                                                 
                                if(s0.structure_id() == "world")    
                                     std::swap(s0,s1);

                                pp01 = pp_struct->
                                        special_geminate_dissociation_positions( rng_, s0, s1, pp.second.position(), reaction_length_ );
                            }
                            else
                            {
                                pp01 = pp_struct->
                                        geminate_dissociation_positions( rng_, s0, s1, pp.second.position(), reaction_length_ );
                            }
                            
                            pp01.first = tx_.apply_boundary( pp01.first );
                            pp01.second = tx_.apply_boundary( pp01.second );
                                                                                   
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s0(
                                tx_.check_overlap(
                                    particle_shape_type(pp01.first, s0.radius()),
                                    pp.first));
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s1(
                                tx_.check_overlap(
                                    particle_shape_type(pp01.second, s1.radius()),
                                    pp.first));
                                    
                            /* Only when both products live in the bulk, check for possible overlap with a surface. */
                            bool surface_overlap = false;
                            if( s0.structure_id() == "world" && s1.structure_id() == "world" )
                            {
                                structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( pp01.first, s0.structure_id() ) );
                                surface_overlap = determine_surface_passage(pp01.first, pp.second.position(), struct_id_distance_pair, s0.radius());

                                struct_id_distance_pair = tx_.get_closest_surface( pp01.second, s1.structure_id() );
                                surface_overlap = determine_surface_passage(pp01.second, pp.second.position(), struct_id_distance_pair, s1.radius());                          
                            }
                                    
                            if (!(overlapped_s0 && overlapped_s0->size() > 0) && !(overlapped_s1 && overlapped_s1->size() > 0) && !surface_overlap)
                                break;
                        }

                        /* After dissociation both particles are allowed to move. */
                        pp01 = make_pair_move(s0, s1, pp01, pp.first);

                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(pp01.first, s0.radius()), pp.first) || !(*vc_)(particle_shape_type(pp01.second, s1.radius()), pp.first))
                            {
                                throw propagation_error("no space due to near particle (vc)");
                            }
                        }

                        tx_.remove_particle(pp.first);
                        const particle_id_pair
                            npp0(tx_.new_particle(s0.id(), pp01.first)),
                            npp1(tx_.new_particle(s1.id(), pp01.second));

                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(),
                                    array_gen(npp0.first, npp1.first),
                                    pp.first));
                        }
                    }
                    break;
                default:
                    throw not_implemented("monomolecular reactions that produce more than two products are not supported");
                }
                return true;
            }
        }
        return false;
    }
        
    /* Function returns the accumulated intrinsic reaction rate between to species. */
    Real k_total(species_id_type const& s0_id, species_id_type const& s1_id)
    {   
        reaction_rules const& rules(rules_.query_reaction_rule(s0_id, s1_id));
        if(::size(rules) == 0)
        {
            return 0;
        }

        Real k_tot = 0;
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            k_tot += (*i).k();
        }
        
        return k_tot;
    }
    
    
    bool fire_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.sid(), pp1.second.sid()));
        if(::size(rules) == 0)
            throw propagation_error("trying to fire a reaction between non-reactive particles");

        species_type s0(tx_.get_species(pp0.second.sid())),
                s1(tx_.get_species(pp1.second.sid()));
        
        boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id() ) );
        
        const Real k_tot( k_total(pp0.second.sid(), pp1.second.sid()) );
        const Real rnd( k_tot * rng_() );
        Real prob = 0;    
        
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);  
            prob += r.k();
            if (prob > rnd)
            {
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const species_type sp(tx_.get_species(product));
                        
                        position_type new_pos( tx_.apply_boundary(
                                        divide(
                                        add(multiply(pp0.second.position(), s1.D()),
                                            multiply(tx_.cyclic_transpose(
                                                pp1.second.position(),
                                                pp0.second.position()), s0.D())),
                                        (s0.D() + s1.D()))) );

                        //For unequal structures, project new_pos on surface.
                        if(sp.structure_id() != s1.structure_id())
                        {
                            boost::shared_ptr<structure_type> sp_struct( tx_.get_structure( sp.structure_id() ) );
                            projected_type const position_on_surface( sp_struct->projected_point(new_pos) );
                            
                            new_pos = tx_.apply_boundary( position_on_surface.first );
                        }

                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, sp.radius()),
                                              pp0.first, pp1.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space due to particle");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius()), 
                                    pp0.first, pp1.first))
                            {
                                throw propagation_error("no space due to particle (vc)");
                            }
                        }
                        
                        if( sp.structure_id() == "world" )
                        {
                            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos, sp.structure_id() ) );
                            
                            if( determine_surface_passage(new_pos, pp0.second.position(), struct_id_distance_pair, sp.radius()) )
                                throw propagation_error("no space due to near surface");
                        }

                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        particle_id_pair npp(tx_.new_particle(product, new_pos));
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(npp.first), pp0.first, pp1.first));
                        }
                        break;
                    }
                case 0:
                    {
                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        break;
                    }
                
                default:
                    throw not_implemented("bimolecular reactions that produce more than one product are not supported");
                }
                
                return true;
            }
        }
        //Should not reach here.
        return false;
    }

        
    bool fire_interaction(particle_id_pair const& pp, boost::shared_ptr<structure_type> const& surface)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid(), surface->sid() ));
        if (::size(rules) == 0)
            throw propagation_error("trying to fire an interaction with a non-interactive surface");
        
        const species_type s0(tx_.get_species(pp.second.sid()));
                
        const Real k_tot( k_total(pp.second.sid(), surface->sid()) );
        const Real rnd( k_tot * rng_() );
        Real prob = 0;  

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i); 
            prob += r.k();

            if (prob > rnd)
            {
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const species_type sp(tx_.get_species(product));

                        const position_type new_pos( tx_.apply_boundary( surface->projected_point( pp.second.position() ).first ) );
                            
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, sp.radius()),
                                              pp.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius()), 
                                    pp.first))
                            {
                                throw propagation_error("no space (vc)");
                            }
                        }

			            const particle_id_pair new_p(
			              pp.first, particle_type(sp.id(),
				            particle_shape_type(new_pos,
					          sp.radius()),
					            sp.D()));

                        tx_.update_particle(new_p);
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(new_p.first), pp.first));
                        }
                        break;
                    }
                case 0:
                    {
                        remove_particle(pp.first);
                        break;
                    }
                
                default:
                    throw not_implemented("surface interactions that produce more than one product are not supported");
                }

                return true;
            }
        }
        //should not reach here.
        return false;
    }
    
    /* Function determines wether old_pos lies on the other side of the surface than new_pos. */
    bool determine_surface_passage(position_type const& new_pos, position_type const& old_pos, 
                                    structure_id_and_distance_pair const& struct_id_distance_pair, length_type const& radius)
    {
        if( struct_id_distance_pair.second > radius )
            return false;
    
        boost::shared_ptr<structure_type> surface( tx_.get_structure( struct_id_distance_pair.first ) );
        
        projected_type const old_projection( surface->projected_point_on_surface( old_pos ) );
        projected_type const new_projection( surface->projected_point_on_surface( new_pos ) );
        
        if( old_projection.second * new_projection.second <= 0 )
            return true;
        else
            return false;
    }
    
    /* Given a position pair and two species, the function generates two new 
       positions based on two moves drawn from the free propegator. */
    position_pair_type const make_pair_move(species_type const& s0, species_type const& s1, position_pair_type const& pp01, 
                                                particle_id_type const& ignore)
    {
        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );
        boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id() ) );

        length_type const r01( s0.radius() + s1.radius() );
        
        position_type temp_pos0, temp_pos1;
        position_pair_type new_pp01( pp01 );
        
        //Randomize which species moves first, and if there will be one or two moves. (50% of the cases)
        bool move_prt0_first( rng_.uniform_int(0, 1) );
        int move_count( 1 + rng_.uniform_int(0, 1) );
        while( move_count-- )
        {
            if( move_prt0_first )
            {
                temp_pos0 = make_move(s0, pp01.first, ignore);
                new_pp01.first = tx_.distance(temp_pos0, new_pp01.second) - r01 < 0 ? pp01.first : temp_pos0;
            }                       
            else
            {
                temp_pos1 = make_move(s1, pp01.second, ignore);
                new_pp01.second = tx_.distance(temp_pos1, new_pp01.first) - r01 < 0 ? pp01.second : temp_pos1;
            }
            move_prt0_first = !move_prt0_first;
        }  
            
        return new_pp01;
    }

    position_type const make_move(species_type const& s0, position_type const& old_pos, particle_id_type const& ignore)
    {
        if(s0.D() == 0)
            return old_pos;
    
        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );

        position_type displacement( s0_struct->
                                        bd_displacement(s0.v() * dt_, std::sqrt(2.0 * s0.D() * dt_), rng_) );

        position_type new_pos(tx_.apply_boundary( add( old_pos, displacement ) ));
                        
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped( 
                                    tx_.check_overlap(
                                    particle_shape_type(new_pos, s0.radius()),
                                    ignore));
                                                                
        if ( (overlapped && overlapped->size() > 0) )
            return old_pos;
        
        if(s0.structure_id() == "world")
        {
            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos, s0.structure_id() ) );
                             
            if( determine_surface_passage(new_pos, old_pos, struct_id_distance_pair, s0.radius()) )
                return old_pos;        
        }
            
        return new_pos;
    }
    
    void remove_particle(particle_id_type const& pid)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));
        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
            queue_.erase(i);
    }

private:
    position_type random_unit_vector()
    {
        position_type v(rng_() - .5, rng_() - .5, rng_() - .5);
        return v / length(v);
    }

private:
    particle_container_type& tx_;
    network_rules_type const& rules_;
    rng_type& rng_;
    Real const dt_;
    int const max_retry_count_;
    reaction_recorder_type* const rrec_;
    volume_clearer_type* const vc_;
    particle_id_vector_type queue_;
    int rejected_move_count_;
    Real const reaction_length_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& newBDPropagator<Ttraits_>::log_(Logger::get_logger("ecell.newBDPropagator"));

#endif /* NEW_BD_PROPAGATOR_HPP */

