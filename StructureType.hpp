#ifndef STRUCTURE_TYPE_HPP
#define STRUCTURE_TYPE_HPP

#include <ostream>
#include <string>
#include <boost/format.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "exceptions.hpp"
#include "utils/get_mapper_mf.hpp"
#include "SpeciesTypeID.hpp"

class ParticleModel;

class StructureType
{
    friend class ParticleModel;
private:
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef SpeciesTypeID                               identifier_type;    // NOTE: we use the same identifier as for the species!
    typedef string_map_type::const_iterator             string_map_iterator;
    typedef boost::iterator_range<string_map_iterator>  attributes_range;
    typedef identifier_type                             structure_type_id_type;

public:
    // Constructor
    // TODO should add StructureType as property of StructureType so that we can
    //      build a hyarchie of StructureTypes living in Structuretypes just like structure living in structures.
    StructureType(): model_(0) {}
 
    // Get the id
    identifier_type const& id() const;

    // Get the id of the structure type the structure_type lives on/in
    structure_type_id_type const& structure_type_id() const;

    // Check equality/inequality
    bool operator==(StructureType const& rhs) const;

    bool operator!=(StructureType const& rhs) const;

    std::string const& operator[](std::string const& name) const;

    std::string& operator[](std::string const& name);

    // Get all the attributes
    attributes_range attributes() const;

    // Get the particle model to which the structuretype is associated
    ParticleModel* model() const
    {
        return model_;
    }

protected:
    void bind_to_model(ParticleModel* model, identifier_type const& id)
    {
        model_ = model; 
        id_ = id;
    }


//////// Member variables
private:
    ParticleModel*          model_;             // to what model is it associated
    identifier_type         id_;                // identifier Note that these are only defined if the object is bound to a model.
    structure_type_id_type  structure_type_id_; // The structure type that the structure_type lives on/in
    string_map_type         attrs_;
};

/////// Inline functions
template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, const StructureType& v)
{
    bool first = true;
    out << "StructureType(id=" << v.id() << ", attributes={";

    typename StructureType::attributes_range attributes(v.attributes());
    for (typename boost::range_const_iterator<
        typename StructureType::attributes_range>::type
            i(attributes.begin()), e(attributes.end()); i != e; ++i)
    {
        typename boost::range_value<typename StructureType::attributes_range>::type
                const& pair(*i);
        if (!first)
            out << ", ";
        out << pair.first << ":" << pair.second;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* STRUCTURE_TYPE_HPP */
