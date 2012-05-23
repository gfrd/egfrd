#ifndef CONNECTIVITY_CONTAINER_HPP
#define CONNECTIVITY_CONTAINER_HPP

#include <boost/lexical_cast.hpp>
#include "exceptions.hpp"


template<typename Tobj_, typename Tdata_, std::size_t num_neighbors>
class ConnectivityContainer
// implements the datastructure of objects that each have exactly 'num_neighbors' neighbors and have Tdata_ data associated
// with each edge from the one of the object to the other.
// NOTE the edges to neighbors are MANDATORY and have a specific ORDER.

// It allows you to check what structure you enter into when you leave the current one.
{
public:
    typedef Tobj_          obj_type;               // is in our case the StructureID
    typedef Tdata_         data_type;              // is the normal vector into the neighboring structure
    typedef std::size_t    size_type;              // the type for the index in the array

    typedef std::pair<obj_type, data_type>                                          obj_data_pair_type;

protected:
    // the structureID of the neighbor and the normal vector pointing from the
    // intersection in the direction of the center of the neighbor
    typedef boost::array<obj_data_pair_type, num_neighbors>                         obj_data_pair_array_type;
    // the mapping from the structure and all its neighbors and how they are connected.
    typedef std::map<obj_type, obj_data_pair_array_type>                            obj_neighbor_objdata_array_map;
    typedef select_second<typename obj_neighbor_objdata_array_map::value_type>      obj_neighbor_objdata_second_selector_type;
    typedef boost::transform_iterator<obj_neighbor_objdata_second_selector_type,
            typename obj_neighbor_objdata_array_map::const_iterator>                obj_neighbor_objdata_array_iterator;

public:
    static const size_type max = num_neighbors;

public:
    void set_neighbor_info(obj_type const& obj, size_type const n, obj_data_pair_type const& obj_data_pair)
    {
        // check that the index is not too large
        if ( n >= max )
            throw illegal_argument(std::string("Index out of range for neighbor (n= ") + boost::lexical_cast<std::string>(n) + "), max= " + boost::lexical_cast<std::string>(max));

        neighbor_mapping_[obj][n] = obj_data_pair;
    }

    obj_data_pair_type get_neighbor_info(obj_type const& id, size_type const n) const
    {
        // check that the index is not too large
        if ( n >= max )
            throw illegal_argument(std::string("Index out of range for neighbor (n= ") + boost::lexical_cast<std::string>(n) + "), max= " + boost::lexical_cast<std::string>(max));

        typename obj_neighbor_objdata_array_map::const_iterator i(neighbor_mapping_.find(id));
        if (neighbor_mapping_.end() == i)
        {
            throw not_found(std::string("Object not found(id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second[n];
    }

    ConnectivityContainer() {}
    // do we start with an empty container or do we construct it as a static object with n structures already connected.

private:
    // mapping StructureID -> array<(StructureID, vector3D), num_neighbors>
    obj_neighbor_objdata_array_map  neighbor_mapping_;
};


////// Inline functions related to ConnectivityContainer

#endif /* CONNECTIVITY_CONTAINER_HPP */

