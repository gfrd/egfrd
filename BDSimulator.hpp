#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
#include "NetworkRules.hpp"
#include "newBDPropagator.hpp"
#include "World.hpp"
#include "ParticleSimulator.hpp"
#include "utils/pair.hpp"

template<typename Tworld_>
struct BDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
};

template<typename Ttraits_>
class BDSimulator: public ParticleSimulator<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;
    typedef typename traits_type::world_type world_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::species_id_type species_id_type;
    typedef typename world_type::species_type species_type;
    typedef typename world_type::length_type length_type;
    typedef typename world_type::structure_type structure_type;
    typedef typename traits_type::time_type time_type;  
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename traits_type::rate_type rate_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef std::pair<const Real, const Real> real_pair;

public:
    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(boost::shared_ptr<world_type> world, 
                boost::shared_ptr<network_rules_type const> network_rules,
                rng_type& rng, Real dt_factor = .1,
                int dissociation_retry_moves = 1, length_type reaction_length_factor = traits_type::MULTI_SHELL_FACTOR)
        : base_type(world, network_rules, rng),
          dt_factor_(dt_factor), num_retries_(dissociation_retry_moves), 
          reaction_length_factor_(reaction_length_factor)
    {
        calculate_dt_and_reaction_length();
    }

    void calculate_dt_and_reaction_length()
    {
        //base_type::dt_ = dt_factor_ * determine_dt(*base_type::world_);
        base_type::dt_ = determine_dt();
        reaction_length_ = determine_reaction_length();
        log_.info( "dt=%e, reaction length=%e", base_type::dt_, reaction_length_ );
    }

    virtual void step()
    {
        _step(base_type::dt_);
    }

    virtual bool step(time_type upto)
    {
        time_type const lt(upto - base_type::t_);
        if (lt <= 0.)
            return false;
        if (base_type::dt_ < lt)
        {
            _step(base_type::dt_);
        }
        else
        {
            _step(lt);
            base_type::t_ = upto;
        }
        
        return true;
    }

    static Real determine_dt(world_type const& world)
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, world.get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        return gsl_pow_2(radius_min * 2) / (D_max * 2);
    }
    
    Real determine_dt()
    {
        std::pair<const Real, const Real> maxD_minr( maxD_minsigma() );
        const Real k_max( get_max_rate() );
        const Real D_max( maxD_minr.first );
        const Real r_min( maxD_minr.second );
        const Real Pacc_max( 0.01 ); //Maximum value of the acceptance probability.
        const Real tau_D( 2 * r_min * r_min / D_max ); //~step over particle - time.
        Real dt;
        
        if( k_max > 0)
        {
            Real dt_temp( 2 * Pacc_max * reaction_length_factor_ * r_min / k_max );
            dt = std::min( dt_temp, dt_factor_ * tau_D ); // dt_factor * tau_D is upper limit of dt.
        }
        else
            dt = dt_factor_ * tau_D;

        return dt;    
    }
    
    length_type determine_reaction_length()
    {
        real_pair maxD_minr( maxD_minsigma() );
        return reaction_length_factor_ * maxD_minr.second;
    }
    
    /* Returns largest diffusion constant and smallest particle radius in the multi. */
    real_pair maxD_minsigma()
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, base_type::world_->get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        
        return real_pair(D_max, radius_min);
    }
    
    /* Functions returns largest _1D_ intrinsic reaction rate in the world. */
    Real get_max_rate()
    {    
        Real k_max(0.);
        int i = 0, j= 0;

        BOOST_FOREACH(species_type s, base_type::world_->get_species())
        {
            if( s.structure_type_id() != base_type::world_->get_def_structure_type_id() )
                continue;

            BOOST_FOREACH(boost::shared_ptr<structure_type> structure, base_type::world_->get_structures())
            {            
                species_id_type const strc_sid( structure->sid() );
                reaction_rules const& rrules((*base_type::network_rules_).query_reaction_rule( s.id(), strc_sid ));
                if (::size(rrules) == 0)
                   continue;

                for (typename boost::range_const_iterator<reaction_rules>::type
                    it(boost::begin(rrules)), e(boost::end(rrules)); it != e; ++it)
                {
                    Real const k( structure->get_1D_rate_surface( (*it).k(), s.radius() ) ); 
    
                    if ( k_max < k )
                        k_max = k;
                }
            }
        }
        
        //Since surface rates are not devided by 2 to compensate for double reaction attempts.
        k_max *= 2.0;
        
        BOOST_FOREACH(species_type s0, base_type::world_->get_species())
        {
            j = 0;
            BOOST_FOREACH(species_type s1, base_type::world_->get_species())
            {
                if( j++ < i )
                    continue;
                
                reaction_rules const& rrules((*base_type::network_rules_).query_reaction_rule( s0.id(), s1.id() ));
                if (::size(rrules) == 0)
                {
                    continue;
                }
                                
                for (typename boost::range_const_iterator<reaction_rules>::type
                it(boost::begin(rrules)), e(boost::end(rrules)); it != e; ++it)
                {
//                    const length_type r01( s0.radius() + s1.radius() ); // TODO
                    Real k;
                    
                    if(s0.structure_type_id() != s1.structure_type_id())
                    {
                        if(s0.structure_type_id() == base_type::world_->get_def_structure_type_id())
                            k = 0.001; //TODO k = base_type::world_->get_structure( s0.structure_type_id() )->get_1D_rate_geminate( (*it).k(), r01 );
                        else
                            k = 0.001; //TODO k = base_type::world_->get_structure( s1.structure_type_id() )->get_1D_rate_geminate( (*it).k(), r01 );
                    }
                    else
                    {
                        k = 0.001;     //TODO k = base_type::world_->get_structure( s0.structure_id() )->get_1D_rate_geminate( (*it).k(), r01 ); 
                    }
                
                    if ( k_max < k )
                        k_max = k;
                }          
            }
            i++;
        }
        
        return k_max;
    }
    
    Real get_reaction_length() const
    {
        return reaction_length_;    
    }
    
    void set_reaction_length_factor(length_type new_reaction_length_factor, time_type new_dt_factor)
    {        
        dt_factor_ = new_dt_factor;
        reaction_length_factor_ = new_reaction_length_factor;
        
        calculate_dt_and_reaction_length();     
    }

protected:
    void _step(time_type dt)
    {
        {
            newBDPropagator<traits_type> propagator(
                *base_type::world_,
                *base_type::network_rules_,
                base_type::rng_,
                dt, num_retries_, reaction_length_,
                base_type::rrec_.get(), 0,
                make_select_first_range(base_type::world_->
                                        get_particles_range()));
            while (propagator());
            LOG_DEBUG(("%d: t=%lg, dt=%lg", base_type::num_steps_, 
                       base_type::t_, dt));
        }
        ++base_type::num_steps_;
        base_type::t_ += dt;
    }

private:
    Real            dt_factor_;
    int const       num_retries_;
    length_type     reaction_length_factor_;
    length_type     reaction_length_;
    static Logger&  log_;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
