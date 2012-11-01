#ifndef PARTICLE_CONTAINER_HPP
#define PARTICLE_CONTAINER_HPP

#include <utility>
#include <boost/shared_ptr.hpp>
#include "generator.hpp"
#include "utils/get_default_impl.hpp"
#include "utils/unassignable_adapter.hpp"
#include "CuboidalRegion.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "SphericalSurface.hpp"
#include "DiskSurface.hpp"

template<typename Ttraits_>
class Transaction;

template<typename Ttraits_>
class ParticleContainer
{
// This class just defines the properties that the ParticleContainer should have.
// The ParticleContainerBase implements some of them.
public:
    typedef Ttraits_ traits_type;

    // some shorthand typename inherited from the traits
    typedef typename traits_type::length_type                   length_type;
    typedef typename traits_type::size_type                     size_type;
    typedef typename traits_type::position_type                 position_type;
    typedef typename traits_type::species_type                  species_type;
    typedef typename traits_type::species_id_type               species_id_type;
    typedef typename traits_type::particle_type                 particle_type;
    typedef typename traits_type::particle_id_type              particle_id_type;
    typedef typename traits_type::structure_type_type           structure_type_type;
    typedef typename traits_type::structure_type_id_type        structure_type_id_type;
    typedef typename traits_type::structure_type                structure_type;
    typedef typename traits_type::structure_id_type             structure_id_type;

    typedef typename particle_type::shape_type                  particle_shape_type;
    typedef Transaction<traits_type>                            transaction_type;

    typedef CuboidalRegion<traits_type>                         cuboidal_region_type;
    typedef PlanarSurface<traits_type>                          planar_surface_type;
    typedef CylindricalSurface<traits_type>                     cylindrical_surface_type;
    typedef DiskSurface<traits_type>                            disk_surface_type;
    typedef SphericalSurface<traits_type>                       spherical_surface_type;

//    typedef std::pair<structure_id_type, boost::shared_ptr<cuboidal_region_type> >       cuboidal_region_id_pair_type;
//    typedef std::pair<structure_id_type, boost::shared_ptr<planar_surface_type> >        planar_surface_id_pair_type;
//    typedef std::pair<structure_id_type, boost::shared_ptr<cylindrical_surface_type> >   cylindrsurf_id_pair_type;
//    typedef std::pair<structure_id_type, boost::shared_ptr<disk_surface_type> >          disk_surface_id_pair_type;
//    typedef std::pair<structure_id_type, boost::shared_ptr<spherical_surface_type> >     spherical_surface_id_pair_type;

    typedef std::set<particle_id_type>                                              particle_id_set;
    typedef std::pair<const particle_id_type, particle_type>                        particle_id_pair;
    typedef abstract_limited_generator<particle_id_pair>                            particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type>                                particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance,
                                 get_default_impl::std::vector>                     particle_id_pair_and_distance_list;
    typedef std::set<structure_id_type>                                             structure_id_set;
    typedef std::pair<const structure_id_type, boost::shared_ptr<structure_type> >  structure_id_pair;
    typedef abstract_limited_generator<structure_id_pair>                           structure_id_pair_generator;
    typedef std::pair<structure_id_pair, length_type>                               structure_id_pair_and_distance;
    typedef unassignable_adapter<structure_id_pair_and_distance,
                                 get_default_impl::std::vector>                     structure_id_pair_and_distance_list;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> >         structure_map;

private:
    typedef select_second<typename structure_map::value_type>                   structure_second_selector_type;
    typedef boost::transform_iterator<structure_second_selector_type,
            typename structure_map::const_iterator>                             structure_iterator;
    typedef std::map<structure_type_id_type, structure_type_type >              structure_type_map;
    typedef select_second<typename structure_type_map::value_type>              structure_type_second_selector_type;
    typedef boost::transform_iterator<structure_type_second_selector_type,
            typename structure_type_map::const_iterator>                        structure_type_iterator;

public:    
    typedef sized_iterator_range<structure_iterator>                            structures_range;
    typedef sized_iterator_range<structure_type_iterator>                       structure_types_range;
    typedef std::pair<position_type, structure_id_type>                         position_structid_pair_type;


    virtual ~ParticleContainer() {};

    virtual size_type num_particles() const = 0;

    virtual length_type world_size() const = 0;

    // Species stuff
    virtual species_type const& get_species(species_id_type const& id) const = 0;

    // StructureType stuff
//    virtual bool add_structure_type(structure_type_type const& structure_type) = 0;   // TODO

    virtual structure_type_type get_structure_type(structure_type_id_type const& sid) const = 0;

    virtual structure_types_range get_structure_types() const = 0;

    virtual structure_type_id_type get_def_structure_type_id() const = 0;

    // Structure stuff
//    virtual structure_id_type add_structure(structure_type const& structure) = 0;   // TODO

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const = 0;
    
    virtual structures_range get_structures() const = 0;
    
    virtual boost::shared_ptr<structure_type> get_some_structure_of_type(structure_type_id_type const& sid) const = 0;

//    virtual bool update_structure(structure_id_pair const& structid_pair) = 0;    

    virtual bool remove_structure(structure_id_type const& id) = 0;

    virtual structure_id_set get_structure_ids(structure_type_id_type const& sid) const = 0;

    virtual structure_id_type get_def_structure_id() const = 0;

    virtual structure_id_pair_and_distance_list* get_close_structures(position_type const& pos, structure_id_type const& current_struct_id,
                                                                      structure_id_type const& ignore) const = 0;

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma) const = 0;

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma, structure_id_type const& ignore) const = 0;

    // Particle stuff
    virtual particle_id_pair new_particle(species_id_type const& sid, structure_id_type const& structure_id,
            position_type const& pos) = 0;

    virtual bool update_particle(particle_id_pair const& pi_pair) = 0;

    virtual bool remove_particle(particle_id_type const& id) = 0;

    virtual particle_id_pair get_particle(particle_id_type const& id) const = 0;
    
    virtual bool has_particle(particle_id_type const& id) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const = 0;

    virtual particle_id_pair_generator* get_particles() const = 0;

    virtual transaction_type* create_transaction() = 0;

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const = 0;

    virtual position_type apply_boundary(position_type const& v) const = 0;

    virtual length_type apply_boundary(length_type const& v) const = 0;

    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id) const = 0;

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const = 0;

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const = 0;

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_type const& structure) const = 0;

};


#endif /* PARTICLE_CONTAINER_HPP */
