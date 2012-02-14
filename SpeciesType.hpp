#ifndef SPECIES_TYPE_HPP
#define SPECIES_TYPE_HPP

#include <ostream>
#include <string>
#include <boost/format.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "SpeciesTypeID.hpp"

#include "exceptions.hpp"
#include "utils/get_mapper_mf.hpp"

class Model;

class SpeciesType
{
    friend class Model;
private:
    typedef get_mapper_mf<std::string, std::string>::type string_map_type;

public:
    typedef SpeciesTypeID identifier_type;                          // The identifier
    typedef string_map_type::const_iterator string_map_iterator;
    typedef boost::iterator_range<string_map_iterator> attributes_range;

public:
    identifier_type const& id() const;

    std::string const& operator[](std::string const& name) const;

    std::string& operator[](std::string const& name);

    attributes_range attributes() const;

    Model* model() const
    // The model with which this species type is associated?
    {
        return model_;
    }

    SpeciesType(): model_(0) {}
    // The constructor initializes the model with '0'
 
protected:
    void bind_to_model(Model* model, identifier_type const& id)
    // Later the species type can be bound to the model.
    // Note that also here the id is fixed -> Only when bound to a model is the id known
    {
        model_ = model; 
        id_ = id;
    }


////////// Member variables
private:
    Model* model_;              // The model to which the species type is bound
    identifier_type id_;        // The id of the species type
    string_map_type attrs_;     // all the attributes?
};


template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, const SpeciesType& v)
// Allows for printing of the object
// -> provides a stream of 'attributename : attribute_value' pairs.
{
    bool first = true;
    out << "SpeciesType(id=" << v.id() << ", attributes={";

    typename SpeciesType::attributes_range attributes(v.attributes());
    for (typename boost::range_const_iterator<
        typename SpeciesType::attributes_range>::type
            i(attributes.begin()), e(attributes.end()); i != e; ++i)
    {
        typename boost::range_value<typename SpeciesType::attributes_range>::type
                const& pair(*i);
        if (!first)
            out << ", ";
        out << pair.first << ":" << pair.second;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* SPECIES_TYPE_HPP */
