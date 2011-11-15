#ifndef TRANSACTION_HPP
#define TRANSACTION_HPP

#include <vector>
#include <map>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include "utils.hpp"
#include "ParticleContainer.hpp"
#include "sorted_list.hpp"
#include "generator.hpp"
#include "utils/unassignable_adapter.hpp"
#include "utils/stringizer.hpp"

template<typename Ttraits_>
class Transaction: public ParticleContainer<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;

    virtual ~Transaction() {}

    virtual particle_id_pair_generator* get_added_particles() const = 0;

    virtual particle_id_pair_generator* get_removed_particles() const = 0;

    virtual particle_id_pair_generator* get_modified_particles() const = 0;

    virtual void rollback() = 0;
};

template<typename Tpc_>
class TransactionImpl: public Transaction<typename Tpc_::traits_type>
{
public:
    typedef Tpc_ particle_container_type;
    typedef typename particle_container_type::traits_type traits_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef std::pair<structure_id_type, length_type> structure_id_and_distance_pair;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;

private:
    typedef std::map<typename particle_id_pair::first_type,
            typename particle_id_pair::second_type> particle_id_pair_set_type;
    typedef sorted_list<std::vector<particle_id_type> > particle_id_list_type;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef select_second<typename structure_map::value_type> surface_second_selector_type;

public:    
    typedef boost::transform_iterator<surface_second_selector_type,
            typename structure_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<surface_iterator> structures_range;


    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        particle_id_pair retval(pc_.new_particle(sid, pos));
        const bool result(added_particles_.push_no_duplicate(retval.first));
        BOOST_ASSERT(result);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        BOOST_ASSERT(removed_particles_.end() ==
                removed_particles_.find(pi_pair.first));
        std::pair<typename particle_id_pair_set_type::iterator, bool> r(
                orig_particles_.insert(particle_id_pair(
                    pi_pair.first, particle_type())));
        if (r.second &&
            added_particles_.end() == added_particles_.find(pi_pair.first))
        {
            modified_particles_.push_no_duplicate(pi_pair.first);
            particle_type _v(pc_.get_particle(pi_pair.first).second);
            std::swap((*r.first).second, _v);
        }
        return pc_.update_particle(pi_pair);
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        std::pair<typename particle_id_pair_set_type::iterator, bool> r(
                orig_particles_.insert(particle_id_pair(
                    id, particle_type())));
        if (r.second)
        {
            particle_type _v(pc_.get_particle(id).second);
            std::swap((*r.first).second, _v);
        }

        if (added_particles_.erase(id) == 0)
        {
            modified_particles_.erase(id);
            const bool result(removed_particles_.push_no_duplicate(id));
            BOOST_ASSERT(result);
        }
        else
        {
            orig_particles_.erase(id);
        }
        return pc_.remove_particle(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        return pc_.get_particle(id);
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return pc_.has_particle(id);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return pc_.check_overlap(s);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return pc_.check_overlap(s, ignore);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return pc_.check_overlap(s, ignore1, ignore2);
    }

    virtual Transaction<traits_type>* create_transaction()
    {
        return new TransactionImpl<particle_container_type>(*this);
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return pc_.get_structure(id);
    }
    
    virtual structures_range get_structures() const
    {
        return pc_.get_structures();
    }    
        
    virtual structure_id_and_distance_pair get_closest_surface(position_type const& pos, structure_id_type const& ignore) const
    {
        return pc_.get_closest_surface( pos, ignore );
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        return pc_.get_species(id);
    }

    virtual size_type num_particles() const
    {
        return pc_.num_particles();
    }

    virtual length_type world_size() const
    {
        return pc_.world_size();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return pc_.get_particles();
    }

    virtual particle_id_pair_generator* get_added_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(added_particles_,
                boost::bind(&TransactionImpl::get_particle, this, _1)));
            
    }

    virtual particle_id_pair_generator* get_removed_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(removed_particles_,
                boost::bind(&TransactionImpl::get_original_particle, this, _1)));
    }

    virtual particle_id_pair_generator* get_modified_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(modified_particles_,
                boost::bind(&TransactionImpl::get_particle, this, _1)));
    }

    virtual void rollback()
    {
        for (typename particle_id_pair_set_type::iterator
                i(orig_particles_.begin()), e(orig_particles_.end());
                i != e; ++i)
        {
            pc_.update_particle(*i);
        }

        for (typename particle_id_list_type::iterator
                i(added_particles_.begin()), e(added_particles_.end());
                i != e; ++i)
        {
            pc_.remove_particle(*i);
        }
        added_particles_.clear();
        modified_particles_.clear();
        removed_particles_.clear();
        orig_particles_.clear();
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return pc_.distance(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return pc_.apply_boundary(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return pc_.apply_boundary(v);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return pc_.cyclic_transpose(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return pc_.cyclic_transpose(p0, p1);
    }

    virtual ~TransactionImpl() {}

    TransactionImpl(particle_container_type& pc): pc_(pc) {}

private:
    particle_id_pair get_original_particle(particle_id_type const& id) const
    {
        typename particle_id_pair_set_type::const_iterator i(orig_particles_.find(id));
        if (orig_particles_.end() == i)
        {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

private:
    particle_container_type& pc_;
    particle_id_list_type added_particles_;
    particle_id_list_type modified_particles_;
    particle_id_pair_set_type orig_particles_;
    particle_id_list_type removed_particles_;
};

#endif /* TRANSACTION_HPP */
