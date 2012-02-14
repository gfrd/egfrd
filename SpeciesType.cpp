#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "SpeciesType.hpp"

SpeciesType::identifier_type const& SpeciesType::id() const
// Gets the id of the SpeciesType (ONLY DEFINED WHEN SpeciesType is bound to model)
{
    if (!model_)
    {
        throw illegal_state("not bound to Model");
    }
    return id_;
}
    
std::string const& SpeciesType::operator[](std::string const& name) const
// Get an attribute value by name. Usage for example: a = 'species_type['name']'
{
    string_map_type::const_iterator i(attrs_.find(name));
    if (i == attrs_.end())
        throw not_found((boost::format("key %s not found") % name).str());
    return (*i).second;
}

std::string& SpeciesType::operator[](std::string const& name)
// Same thing as above but messier?
{
    return attrs_[name];
}

SpeciesType::attributes_range SpeciesType::attributes() const
// Gets all the attributes
{
    return attributes_range(attrs_.begin(), attrs_.end());
}
