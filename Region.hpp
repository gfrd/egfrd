#ifndef REGION_HPP
#define REGION_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <sstream>
#include "ParticleSimulationStructure.hpp"
#include "Box.hpp"

template<typename Ttraits_>
class Region: public ParticleSimulationStructure<Ttraits_>
{
public:
    typedef ParticleSimulationStructure<Ttraits_>       base_type;
    typedef typename base_type::identifier_type         identifier_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;

public:
    virtual ~Region() {}

    // Constructor
    Region(identifier_type const& id, structure_type_id_type const& sid): base_type(id, sid) {}
};



template<typename Ttraits_, typename Tshape_>
class BasicRegionImpl: public Region<Ttraits_>
{
public:
    typedef Region<Ttraits_>                            base_type;
    typedef Tshape_                                     shape_type;
    typedef typename base_type::identifier_type         identifier_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::length_type             length_type;
    typedef typename base_type::position_type           position_type;
    typedef std::pair<position_type, bool>              position_flag_pair_type;

public:
    virtual ~BasicRegionImpl() {}

    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
    }

    virtual bool operator==(Structure<typename Ttraits_::world_type::traits_type> const& rhs) const
    {
        BasicRegionImpl const* _rhs(dynamic_cast<BasicRegionImpl const*>(&rhs));
        return _rhs &&
               base_type::real_id_ == rhs.real_id() &&
               base_type::sid_ == rhs.sid() &&
               shape_ == _rhs->shape();
    }

    virtual std::size_t hash() const
    {    // Displacements are not deflected on cylinders (yet),
    // but this function has to be defined for every shape to be used in structure
    // For now it just returns the new position
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<identifier_type>()(base_type::id_) ^
               hash<structure_type_id_type>()(base_type::sid_) ^
               hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Region(" << base_type::real_id_ << ":" << shape() << ")";
        return out.str();
    }

    virtual std::pair<position_type, length_type> projected_point(position_type const& pos) const
    {
        return ::projected_point(shape(), pos);
    }
    
    virtual std::pair<position_type, length_type> projected_point_on_surface(position_type const& pos) const
    {
        return ::projected_point(shape(), pos);
    }
    
    virtual length_type distance(position_type const& pos) const
    {
        return ::distance(shape(), pos);
    }
    
    virtual position_type const& position() const
    {
        return shape_.position();
    }

    virtual position_flag_pair_type deflect(position_type const& pos0, position_type const& displacement) const
    {
        return ::deflect(shape(), pos0, displacement);
    }
    
    // The constructor
    BasicRegionImpl(identifier_type const& id, structure_type_id_type const& sid, shape_type const& shape)
        : base_type(id, sid), shape_(shape) {}

protected:
    shape_type shape_;
};

#endif /* REGION_HPP */
