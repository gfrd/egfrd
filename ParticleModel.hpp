#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "Defs.hpp"
#include "Model.hpp"
#include "StructureType.hpp"

class ParticleModel: public Model
// The ParticleModel class add spatial properties (e.g. particles) to the model. The model by itself could
// namely also just model well mixed systems.
{
public:
    // defining shorthands for the used types
    typedef Model                                   base_type;

    typedef StructureType                                           structure_type_type;
    typedef structure_type_type::identifier_type                    structure_type_id_type;

    typedef std::map<structure_type_id_type, boost::shared_ptr<structure_type_type> >   structure_type_map;
    typedef select_second<structure_type_map::value_type>                               structure_type_second_selector_type;

public:
    typedef boost::transform_iterator<structure_type_second_selector_type,
            structure_type_map::const_iterator>                                         structure_type_iterator;
    typedef boost::iterator_range<structure_type_iterator>                               structure_types_range;

public:
    // Constructor
    ParticleModel();

    virtual ~ParticleModel();

    // Gets a structure type
    boost::shared_ptr<structure_type_type> get_structure_type_by_id(structure_type_id_type const& id) const;

    // Add a structure type
    void add_structure_type(boost::shared_ptr<structure_type_type> const& structure_type);

    // Get all structure types in the model
    structure_types_range get_structure_types() const;

    // Gets and sets the default structure_type
    structure_type_id_type get_def_structure_type_id() const;

/////// Member variables
public:
    structure_type_map      structure_type_map_;            // mapping: structure_type_id -> structure_type
    structure_type_id_type  default_structure_type_id_;     // The id of the default structure_type ("world")
};


#endif /* MODEL_HPP */
