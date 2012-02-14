#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/lexical_cast.hpp>

#include "utils/fun_wrappers.hpp"

#include "SpeciesType.hpp"
#include "NetworkRules.hpp"
#include "BasicNetworkRulesImpl.hpp"
#include "Model.hpp"

Model::Model(): network_rules_(new BasicNetworkRulesImpl())
// The model class is more or less an implementation of the network rules
// The network rules define the possible reactions between species.
{
}

Model::~Model()
{
}

void Model::add_species_type(boost::shared_ptr<species_type_type> const& species)
{
    species->bind_to_model(this, species_type_id_generator_());
    species_type_map_.insert(std::make_pair(species->id(), species));
    // the species_type_map_ is a mapping (container) between the species_id and the species object
}

boost::shared_ptr<Model::species_type_type> Model::get_species_type_by_id(species_id_type const& id) const
// Returns the species if you give it an species 'id'
{
    species_type_map_type::const_iterator i(species_type_map_.find(id));
    if (species_type_map_.end() == i)
    {
        throw not_found(boost::lexical_cast<std::string>(id));
    }

    return (*i).second;
}

Model::species_type_range Model::get_species_types() const
// Just gets all the species in the model
{
    return species_type_range(
        species_type_iterator(species_type_map_.begin(), species_second_selector_type()),
        species_type_iterator(species_type_map_.end(), species_second_selector_type()));
}

std::string const& Model::operator[](std::string const& name) const
// Gets the name of the species? Or gets the value of an attribute by name?
{
    string_map_type::const_iterator i(attrs_.find(name));
    if (i == attrs_.end())
        throw not_found((boost::format("key %s not found") % name).str());
    return (*i).second;
}

std::string& Model::operator[](std::string const& name)
{
    return attrs_[name];
}
