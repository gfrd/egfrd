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
// A transaction 'is a' ParticleContainer that makes it possible to revert the moves of the particles
// The Transaction is actually just an abstract datatype, not a class from which can be instantiated.
{
public:
    typedef Ttraits_ traits_type;

    // useful shorthands for types that we are using.
    typedef typename traits_type::length_type                   length_type;
    typedef typename traits_type::position_type                 position_type;
    typedef typename traits_type::size_type                     size_type;
    typedef typename traits_type::particle_type                 particle_type;
    typedef typename traits_type::particle_id_type              particle_id_type;
    typedef typename particle_type::shape_type                  particle_shape_type;
    typedef typename traits_type::species_type                  species_type;
    typedef typename traits_type::species_id_type               species_id_type;

    typedef std::pair<const particle_id_type, particle_type>    particle_id_pair;
    typedef abstract_limited_generator<particle_id_pair>        particle_id_pair_generator;

    virtual ~Transaction() {}

    virtual particle_id_pair_generator* get_added_particles() const = 0;

    virtual particle_id_pair_generator* get_removed_particles() const = 0;

    virtual particle_id_pair_generator* get_modified_particles() const = 0;

    virtual void rollback() = 0;
};

template<typename Tpc_>
class TransactionImpl: public Transaction<typename Tpc_::traits_type>
// Tpc_ is the 'particlecontainer_type'
{
public:
    typedef Tpc_ particle_container_type;
    typedef typename particle_container_type::traits_type   traits_type;        // get the traits that are passed on from the particlecontainer.

    // define shorthands to the types that we are using here.
    typedef typename traits_type::length_type               length_type;
    typedef typename traits_type::size_type                 size_type;
    typedef typename traits_type::position_type             position_type;
    typedef typename traits_type::species_type              species_type;
    typedef typename traits_type::species_id_type           species_id_type;
    typedef typename traits_type::particle_type             particle_type;
    typedef typename traits_type::particle_id_type          particle_id_type;
    typedef typename particle_type::shape_type              particle_shape_type;
    typedef typename traits_type::structure_type_type       structure_type_type;
    typedef typename traits_type::structure_type_id_type    structure_type_id_type;
    typedef typename traits_type::structure_type            structure_type;
    typedef typename traits_type::structure_id_type         structure_id_type;

    typedef typename particle_container_type::particle_id_pair_and_distance_list    particle_id_pair_and_distance_list;
    typedef typename particle_container_type::structure_id_pair_and_distance_list   structure_id_pair_and_distance_list;

    typedef std::pair<const particle_id_type, particle_type>    particle_id_pair;
    typedef abstract_limited_generator<particle_id_pair>        particle_id_pair_generator;
    typedef std::pair<structure_id_type, length_type>           structure_id_and_distance_pair;

private:
    typedef std::map<typename particle_id_pair::first_type,
            typename particle_id_pair::second_type>                         particle_id_pair_set_type;
    typedef sorted_list<std::vector<particle_id_type> >                     particle_id_list_type;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef select_second<typename structure_map::value_type>               structure_second_selector_type;
    typedef std::map<structure_type_id_type, structure_type_type >          structure_type_map;
    typedef select_second<typename structure_type_map::value_type>          structure_type_second_selector_type;

public:    
    typedef boost::transform_iterator<structure_second_selector_type,
            typename structure_map::const_iterator>                         structure_iterator;
    typedef sized_iterator_range<structure_iterator>                        structures_range;
    typedef boost::transform_iterator<structure_type_second_selector_type,
            typename structure_type_map::const_iterator>                    structure_type_iterator;
    typedef typename particle_container_type::structure_types_range         structure_types_range;
    typedef typename particle_container_type::structure_id_set              structure_id_set;
    typedef typename particle_container_type::structure_id_pair             structure_id_pair;
    typedef typename particle_container_type::position_structid_pair_type   position_structid_pair_type;

    // Particle Stuff
    virtual particle_id_pair new_particle(species_id_type const& sid, structure_id_type const& structure_id,
            position_type const& pos)
    {
        particle_id_pair retval(pc_.new_particle(sid, structure_id, pos));
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

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma) const
    {
        return pc_.check_surface_overlap(s, old_pos, current, sigma);
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma, structure_id_type const& ignore) const
    {
        return pc_.check_surface_overlap(s, old_pos, current, sigma, ignore);
    }

    virtual Transaction<traits_type>* create_transaction()
    {
        return new TransactionImpl<particle_container_type>(*this);
    }

    ///// All the other methods are passed on to the associated ParticleContainer.
    // Species stuff
    virtual species_type const& get_species(species_id_type const& id) const
    {
        return pc_.get_species(id);
    }

    // StructureType stuff
//    virtual bool add_structure_type(structure_type_type const& structure_type);   // TODO

    virtual structure_type_type get_structure_type(structure_type_id_type const& sid) const
    {
        return pc_.get_structure_type(sid);
    }
    virtual structure_types_range get_structure_types() const
    {
        return pc_.get_structure_types();
    }
    virtual structure_type_id_type get_def_structure_type_id() const
    {
        return pc_.get_def_structure_type_id();
    }    

    // Start Structure stuff
//    virtual structure_id_type add_structure(structure_type const& structure);   // TODO

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return pc_.get_structure(id);
    }
    virtual structures_range get_structures() const
    {
        return pc_.get_structures();
    }    
    template <typename Tstructid_pair_>
    bool update_structure(Tstructid_pair_ const& structid_pair)
    {
        return pc_.update_structure(structid_pair);
    }
    virtual bool remove_structure(structure_id_type const& id)
    {
        return pc_.remove_structure(id);
    }
    virtual structure_id_set get_structure_ids(structure_type_id_type const& sid) const
    {
        return pc_.get_structure_ids(sid);
    }
    virtual structure_id_type get_def_structure_id() const
    {
        return pc_.get_def_structure_id();
    }    
    virtual structure_id_and_distance_pair get_closest_surface(position_type const& pos, structure_id_type const& ignore) const
    {
        return pc_.get_closest_surface( pos, ignore );
    }
    virtual structure_id_pair_and_distance_list* get_close_structures(position_type const& pos, structure_id_type const& current_struct_id,
                                                                      structure_id_type const& ignore) const
    {
        return pc_.get_close_structures(pos, current_struct_id, ignore);
    }
    // End structure stuff


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

    // Begin Transaction methods. The following methods are prototyped in the base class Transaction (see above)
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
    // end Transaction methods.


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

    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id,
                                                       const boost::shared_ptr<const structure_type> structure) const
    {
        return pc_.apply_boundary(pos_struct_id, structure);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return pc_.cyclic_transpose(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return pc_.cyclic_transpose(p0, p1);
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         const boost::shared_ptr<const structure_type> structure) const
    {
        return pc_.cyclic_transpose(pos_struct_id, structure);
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

//////// Member variables.
private:
    particle_container_type&    pc_;                    // the associated particle container
    particle_id_list_type       added_particles_;
    particle_id_list_type       modified_particles_;
    particle_id_pair_set_type   orig_particles_;
    particle_id_list_type       removed_particles_;
};

#endif /* TRANSACTION_HPP */
