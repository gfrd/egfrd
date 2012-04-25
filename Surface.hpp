#ifndef SURFACE_HPP
#define SURFACE_HPP

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
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Plane.hpp"

template<typename Ttraits_>
class Surface: public ParticleSimulationStructure<Ttraits_>
{
public:
    typedef ParticleSimulationStructure<Ttraits_>       base_type;
    typedef typename base_type::structure_name_type     structure_name_type;    // This is just a string
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::length_type             length_type;

public:
    virtual ~Surface() {}

    // Constructor
    Surface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id):
         base_type(name, sid, parent_struct_id) {}

    virtual length_type minimal_distance(length_type const& radius) const = 0;
};





template<typename Ttraits_, typename Tshape_>
class BasicSurfaceImpl: public Surface<Ttraits_>
{
public:
    typedef Surface<Ttraits_> base_type;
    typedef Tshape_ shape_type;

    typedef typename base_type::structure_name_type     structure_name_type;
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::length_type             length_type;
    typedef typename base_type::position_type           position_type;
    typedef std::pair<position_type, length_type>       projected_type;
    typedef std::pair<position_type, bool>              position_flag_pair_type;

public:
    virtual ~BasicSurfaceImpl() {}

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
        BasicSurfaceImpl const* _rhs(dynamic_cast<BasicSurfaceImpl const*>(&rhs));
        return _rhs &&
               base_type::id_ == rhs.id() &&
               base_type::sid_ == rhs.sid() &&
               shape_ == _rhs->shape();
    }
    
    virtual std::size_t hash() const
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<structure_name_type>()(base_type::name_) ^
               hash<structure_type_id_type>()(base_type::sid_) ^
               hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Surface(" << base_type::id_ << ":" << shape() << ")";
        return out.str();
    }

    virtual projected_type projected_point(position_type const& pos) const
    {
        return ::projected_point(shape(), pos);
    }
    
    virtual projected_type projected_point_on_surface(position_type const& pos) const
    {
        return ::projected_point_on_surface(shape(), pos);
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
    
    virtual position_type deflect_back(position_type const& pos, position_type const& u_z) const
    {
        return ::deflect_back(shape(), pos, u_z);
    }
    
    // Constructor
    BasicSurfaceImpl(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id), shape_(shape) {}

protected:
    shape_type shape_;
};

#endif /* SURFACE_HPP */
