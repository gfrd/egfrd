#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "StructureType.hpp"

StructureType::identifier_type const& StructureType::id() const
{
    if (!model_)
    {
        throw illegal_state("not bound to Model");
    }
    return id_;
}
    
StructureType::structure_type_id_type const& StructureType::structure_type_id() const
{
    if (!structure_type_id_)
    {
        throw illegal_state("no structure_type defined");
    }
    return structure_type_id_;
}

// Check equality/inequality
bool StructureType::operator==(StructureType const& rhs) const
{
    return id_ == rhs.id() &&
           structure_type_id_ == rhs.structure_type_id() &&
           model_ == rhs.model();
}

bool StructureType::operator!=(StructureType const& rhs) const
{
    return !operator==(rhs);
}


std::string const& StructureType::operator[](std::string const& name) const
{
    string_map_type::const_iterator i(attrs_.find(name));
    if (i == attrs_.end())
        throw not_found((boost::format("key %s not found") % name).str());
    return (*i).second;
}

std::string& StructureType::operator[](std::string const& name)
{
    return attrs_[name];
}

StructureType::attributes_range StructureType::attributes() const
{
    return attributes_range(attrs_.begin(), attrs_.end());
}
