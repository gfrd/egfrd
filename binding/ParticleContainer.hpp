#ifndef BINDING_PARTICLE_CONTAINER_HPP
#define BINDING_PARTICLE_CONTAINER_HPP

#include <utility>
#include <map>

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/override.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include "../exceptions.hpp"
#include "../utils/get_default_impl.hpp"
#include "../generator.hpp"
#include "../utils/pair.hpp"

#include "peer/compat.h"
#include "peer/wrappers/generator/pyiterator_generator.hpp"
#include "peer/wrappers/generator/generator_wrapper.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/converters/generator/to_python.hpp"

namespace binding {

//// Converters. These convert C++ types to something that Python can handle.

// 
template<typename Timpl_>
struct particle_id_pair_generator_converter
{
    typedef Timpl_ native_type;

    struct to_native_converter
    {
        static void* convert(PyObject* ptr)
        {
            boost::python::handle<> iter(
                boost::python::allow_null(PyObject_GetIter(ptr)));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }
            return new peer::wrappers::pyiterator_generator<
                    typename native_type::result_type>(iter);
        }

        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }
    };

    // Also need to register the converter.
    static void __register()
    {
        peer::wrappers::generator_wrapper<
                ptr_generator<native_type, std::auto_ptr<native_type> > >
                ::__register_class("ParticleIDPairGenerator");
        boost::python::to_python_converter<
                native_type*,
                peer::converters::ptr_generator_to_pyiterator_converter<
                    native_type, std::auto_ptr<native_type> > >();
                    

        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Timpl_>
struct structure_id_pair_generator_converter
{
    typedef Timpl_ native_type;

    struct to_native_converter
    {
        static void* convert(PyObject* ptr)
        {
            boost::python::handle<> iter(
                boost::python::allow_null(PyObject_GetIter(ptr)));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }
            return new peer::wrappers::pyiterator_generator<
                    typename native_type::result_type>(iter);
        }

        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }
    };

    // Also need to register the converter.
    static void __register()
    {
        peer::wrappers::generator_wrapper<
                ptr_generator<native_type, std::auto_ptr<native_type> > >
                ::__register_class("StructureIDPairGenerator");
        boost::python::to_python_converter<
                native_type*,
                peer::converters::ptr_generator_to_pyiterator_converter<
                    native_type, std::auto_ptr<native_type> > >();
                    

        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Timpl_>
struct particle_id_pair_and_distance_list_converter
{
    typedef Timpl_ native_type;

    struct to_python_converter
    {
        template<typename T_>
        struct policy
        {
            typedef typename boost::range_size<T_>::type size_type;
            typedef typename boost::range_value<T_>::type value_type;
            typedef value_type const& reference;
            typedef value_type const& const_reference;
            typedef typename boost::range_const_iterator<T_>::type iterator;
            typedef typename boost::range_const_iterator<T_>::type const_iterator;

            static size_type size(T_ const& c)
            {
                return ::size(c);
            }

            static void set(T_& c, size_type i, const_reference v)
            {
                c.set(i, v);
            }

            static reference get(T_ const& c, size_type i)
            {
                return c.at(i);
            }

            static iterator begin(T_ const& c)
            {
                return boost::begin(c);
            }

            static iterator end(T_ const& c)
            {
                return boost::end(c);
            }
        };

        typedef peer::wrappers::stl_container_wrapper<native_type, std::auto_ptr<native_type>, peer::wrappers::default_policy_generator<policy> > wrapper_type;
        static PyObject* convert(native_type* v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(std::auto_ptr<native_type>(v ? v: new native_type())));
        }
    };

    struct to_native_converter
    {
        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }

        static void* convert(PyObject* pyo)
        {
            if (pyo == Py_None)
            {
                return 0;
            }

            Py_ssize_t hint(-1);
            hint = PyObject_Size(pyo);
            if (hint < 0)
            {
                PyErr_Clear();
            }

            boost::python::handle<> iter(PyObject_GetIter(pyo));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }

            if (hint < 0)
            {
                hint = compat_PyObject_LengthHint(iter.get());
                if (hint < 0)
                {
                    hint = 0;
                }
            }

            native_type* obj(new native_type);
            obj->reserve(hint);
            while (PyObject* item = PyIter_Next(iter.get()))
            {
                obj->push_back(boost::python::extract<typename native_type::value_type>(item)());
            }

            return obj;
        }
    };

    // Also need to register the converter.
    static void __register()
    {
        to_python_converter::wrapper_type::__class_init__("ParticleIDAndDistanceVector", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type*, to_python_converter>();
        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Timpl_>
struct structure_id_pair_and_distance_list_converter
{
    typedef Timpl_ native_type;

    struct to_python_converter
    {
        template<typename T_>
        struct policy
        {
            typedef typename boost::range_size<T_>::type size_type;
            typedef typename boost::range_value<T_>::type value_type;
            typedef value_type const& reference;
            typedef value_type const& const_reference;
            typedef typename boost::range_const_iterator<T_>::type iterator;
            typedef typename boost::range_const_iterator<T_>::type const_iterator;

            static size_type size(T_ const& c)
            {
                return ::size(c);
            }

            static void set(T_& c, size_type i, const_reference v)
            {
                c.set(i, v);
            }

            static reference get(T_ const& c, size_type i)
            {
                return c.at(i);
            }

            static iterator begin(T_ const& c)
            {
                return boost::begin(c);
            }

            static iterator end(T_ const& c)
            {
                return boost::end(c);
            }
        };

        typedef peer::wrappers::stl_container_wrapper<native_type, std::auto_ptr<native_type>, peer::wrappers::default_policy_generator<policy> > wrapper_type;
        static PyObject* convert(native_type* v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(std::auto_ptr<native_type>(v ? v: new native_type())));
        }
    };

    struct to_native_converter
    {
        static PyTypeObject const* expected_pytype()
        {
            return &PyBaseObject_Type;
        }

        static void* convert(PyObject* pyo)
        {
            if (pyo == Py_None)
            {
                return 0;
            }

            Py_ssize_t hint(-1);
            hint = PyObject_Size(pyo);
            if (hint < 0)
            {
                PyErr_Clear();
            }

            boost::python::handle<> iter(PyObject_GetIter(pyo));
            if (!iter)
            {
                boost::python::throw_error_already_set();
            }

            if (hint < 0)
            {
                hint = compat_PyObject_LengthHint(iter.get());
                if (hint < 0)
                {
                    hint = 0;
                }
            }

            native_type* obj(new native_type);
            obj->reserve(hint);
            while (PyObject* item = PyIter_Next(iter.get()))
            {
                obj->push_back(boost::python::extract<typename native_type::value_type>(item)());
            }

            return obj;
        }
    };

    // Also need to register the converter.
    static void __register()
    {
        to_python_converter::wrapper_type::__class_init__("StructureIDAndDistanceVector", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type*, to_python_converter>();
        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};


// ParticleContainer wrapper class
template<typename Tbase_>
class ParticleContainerWrapper
    : public Tbase_, public boost::python::wrapper<Tbase_>
{
public:
    typedef boost::python::wrapper<Tbase_> py_wrapper_type;
    typedef Tbase_ wrapped_type;
    // Get types that are defined in the ParticleContainer class that we are wrapping.
    typedef typename wrapped_type::size_type                            size_type;
    typedef typename wrapped_type::length_type                          length_type;
    typedef typename wrapped_type::particle_id_type                     particle_id_type;
    typedef typename wrapped_type::particle_shape_type                  particle_shape_type;
    typedef typename wrapped_type::position_type                        position_type;
    typedef typename wrapped_type::species_type                         species_type;
    typedef typename wrapped_type::species_id_type                      species_id_type;
    typedef typename wrapped_type::structure_type                       structure_type;
    typedef typename wrapped_type::structure_id_type                    structure_id_type;
    typedef typename wrapped_type::structure_type_type                  structure_type_type;
    typedef typename wrapped_type::structure_type_id_type               structure_type_id_type;

    typedef typename wrapped_type::particle_id_pair                     particle_id_pair;
    typedef typename wrapped_type::structure_id_pair                    structure_id_pair;
    typedef typename wrapped_type::particle_id_pair_generator           particle_id_pair_generator;
    typedef typename wrapped_type::structure_id_pair_generator          structure_id_pair_generator;
    typedef typename wrapped_type::particle_id_pair_and_distance_list   particle_id_pair_and_distance_list;
    typedef typename wrapped_type::structure_id_pair_and_distance_list  structure_id_pair_and_distance_list;
    typedef typename wrapped_type::position_structid_pair_type          position_structid_pair_type;
    typedef typename wrapped_type::structure_id_set                     structure_id_set;
    typedef typename wrapped_type::structures_range                     structures_range;
    typedef typename wrapped_type::structure_types_range                structure_types_range;

//    typedef typename wrapped_type::cuboidal_region_id_pair_type         cuboidal_region_id_pair_type;
//    typedef typename wrapped_type::planar_surface_id_pair_type          planar_surface_id_pair_type;
//    typedef typename wrapped_type::cylindrsurf_id_pair_type             cylindrsurf_id_pair_type;
//    typedef typename wrapped_type::disk_surface_id_pair_type            disk_surface_id_pair_type;
//    typedef typename wrapped_type::spherical_surface_id_pair_type       spherical_surface_id_pair_type;


    typedef typename wrapped_type::transaction_type                     transaction_type;

public:    

    //
    virtual ~ParticleContainerWrapper() {}

    virtual size_type num_particles() const
    {
        boost::python::handle<> retval(
            boost::python::allow_null(
                PyObject_GetAttrString(
                    boost::python::detail::wrapper_base_::get_owner(*this),
                    "num_particles")));
        return boost::python::extract<size_type>(retval.get())();
    }

    virtual length_type world_size() const
    {
        boost::python::handle<> retval(
            boost::python::allow_null(
                PyObject_GetAttrString(
                    boost::python::detail::wrapper_base_::get_owner(*this),
                    "world_size")));
        return boost::python::extract<length_type>(retval.get())();
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_species")(id).template unchecked<species_type const&>();
    }

    // structure stuff
    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_structure")(id);
    }
    
    virtual structures_range get_structures() const
    {
        return py_wrapper_type::get_override("get_structures")();
    }
/*
    virtual bool update_structure(cuboidal_region_id_pair_type const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }

    virtual bool update_structure(planar_surface_id_pair_type const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }

    virtual bool update_structure(cylindrsurf_id_pair_type const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }

    virtual bool update_structure(disk_surface_id_pair_type const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }

    virtual bool update_structure(spherical_surface_id_pair_type const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }

    template<typename Tstructid_pair_>
    bool update_structure(Tstructid_pair_ const& structid_pair)
    {
        return py_wrapper_type::get_override("update_structure")(structid_pair);
    }
*/
    
    virtual bool remove_structure(structure_id_type const& id)
    {
        return py_wrapper_type::get_override("remove_structure")(id);
    }

    virtual structure_id_set get_structure_ids(structure_type_id_type const& sid) const
    {
        return py_wrapper_type::get_override("get_structure_ids")(sid);
    }

    virtual structure_id_type get_def_structure_id() const
    {
        return py_wrapper_type::get_override("get_def_structure_id")();
    }

    virtual structure_id_pair_and_distance_list* get_close_structures(position_type const& pos, structure_id_type const& current_struct_id,
                                                                      structure_id_type const& ignore) const
    {
        return py_wrapper_type::get_override("get_close_structures")(
                pos, current_struct_id, boost::python::make_tuple(ignore))
               .template unchecked<structure_id_pair_and_distance_list*>();
    }

    // end structure stuff

    // Begin StructureType stuff
    virtual structure_type_type get_structure_type(structure_type_id_type const& sid) const
    {
        return py_wrapper_type::get_override("get_structure_type")(sid);
    }

    virtual structure_types_range get_structure_types() const
    {
        return py_wrapper_type::get_override("get_structure_types")();
    }

    virtual structure_type_id_type get_def_structure_type_id() const
    {
        return py_wrapper_type::get_override("get_def_structure_type_id")();
    }
    // end StructureType stuff

    virtual particle_id_pair new_particle(species_id_type const& sid,
            structure_id_type const& structure_id, position_type const& pos)
    {
        return py_wrapper_type::get_override("new_particle")(sid, structure_id, pos);
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return py_wrapper_type::get_override("update_particle")(pi_pair);
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return py_wrapper_type::get_override("remove_particle")(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        return py_wrapper_type::get_override("get_particle")(id);
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return py_wrapper_type::get_override("has_particle")(id);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_id_pair const& s) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s.second.shape(), boost::python::make_tuple(s.first))
                .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return py_wrapper_type::get_override("check_overlap")(s, boost::python::tuple())
                .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s, boost::python::make_tuple(ignore))
               .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return py_wrapper_type::get_override("check_overlap")(
                s, boost::python::make_tuple(ignore1, ignore2))
               .template unchecked<particle_id_pair_and_distance_list*>();
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma) const
    {
        return py_wrapper_type::get_override("check_surface_overlap")(
                s, old_pos, current, sigma)
               .template unchecked<structure_id_pair_and_distance_list*>();
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma, structure_id_type const& ignore) const
    {
        return py_wrapper_type::get_override("check_surface_overlap")(
                s, old_pos, current, sigma, ignore)
               .template unchecked<structure_id_pair_and_distance_list*>();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return py_wrapper_type::get_override("__iter__")()
                .template unchecked<particle_id_pair_generator*>();
    }

    virtual transaction_type* create_transaction()
    {
        return py_wrapper_type::get_override("create_transaction")()
                .template unchecked<transaction_type*>();
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return py_wrapper_type::get_override("distance")(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return py_wrapper_type::get_override("apply_boundary")(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return py_wrapper_type::get_override("apply_boundary")(v);
    }

    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id) const
    {
        return py_wrapper_type::get_override("apply_boundary")(pos_struct_id);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return py_wrapper_type::get_override("cyclic_transpose")(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return py_wrapper_type::get_override("cyclic_transpose")(p0, p1);
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
                                                         structure_type const& structure) const
    {
        return py_wrapper_type::get_override("cyclic_transpose")(pos_struct_id, structure);
    }

};


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_particle_container_class(
        char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // registering converters from standard boost templates.
    peer::converters::register_tuple_converter<typename impl_type::particle_id_pair>();
    peer::converters::register_tuple_converter<typename impl_type::structure_id_pair>();
    peer::converters::register_tuple_converter<typename impl_type::particle_id_pair_and_distance>();
    peer::converters::register_tuple_converter<typename impl_type::structure_id_pair_and_distance>();
    peer::converters::register_tuple_converter<typename impl_type::position_structid_pair_type>();

    // register the converters that are defined above (not sure what this does)
    particle_id_pair_and_distance_list_converter<typename impl_type::particle_id_pair_and_distance_list>::__register();
    particle_id_pair_generator_converter<typename impl_type::particle_id_pair_generator>::__register();
    structure_id_pair_and_distance_list_converter<typename impl_type::structure_id_pair_and_distance_list>::__register();
    structure_id_pair_generator_converter<typename impl_type::structure_id_pair_generator>::__register();

    // defining the python class
    return class_<ParticleContainerWrapper<impl_type>, boost::noncopyable>(name)
        .add_property("num_particles", &impl_type::num_particles)
        .add_property("world_size", &impl_type::world_size)
        // Species stuff
        .def("get_species",
            pure_virtual((typename impl_type::species_type const&(impl_type::*)(typename impl_type::species_id_type const&) const)&impl_type::get_species),
            return_internal_reference<>())
        // StructureType stuff
        .def("get_structure_type", pure_virtual(&impl_type::get_structure_type))
        .def("get_structure_types", pure_virtual(&impl_type::get_structure_types))
        .def("get_def_structure_type_id", pure_virtual(&impl_type::get_def_structure_type_id))
        // Structure stuff
        .def("get_structure", pure_virtual(&impl_type::get_structure))
        .def("get_structures", pure_virtual(&impl_type::get_structures))
//        .def("update_structure", pure_virtual(&impl_type::update_structure))
//        .def("update_structure", &impl_type::template update_structure<typename impl_type::cuboidal_region_id_pair_type>)
//        .def("update_structure", &impl_type::template update_structure<typename impl_type::planar_surface_id_pair_type>)
//        .def("update_structure", &impl_type::template update_structure<typename impl_type::cylindrsurf_id_pair_type>)
//        .def("update_structure", &impl_type::template update_structure<typename impl_type::disk_surface_id_pair_type>)
//        .def("update_structure", &impl_type::template update_structure<typename impl_type::spherical_surface_id_pair_type>)
        .def("remove_structure", pure_virtual(&impl_type::remove_structure))
        .def("get_structure_ids", pure_virtual(&impl_type::get_structure_ids))
        .def("get_def_structure_id", pure_virtual(&impl_type::get_def_structure_id))
        .def("get_close_structures", pure_virtual((typename impl_type::structure_id_pair_and_distance_list*(impl_type::*)(typename impl_type::position_type const& pos, typename impl_type::structure_id_type const&, typename impl_type::structure_id_type const&) const)&impl_type::get_close_structures), return_value_policy<return_by_value>())
        // Particle stuff
        .def("new_particle", pure_virtual(&impl_type::new_particle))
        .def("update_particle", pure_virtual(&impl_type::update_particle))
        .def("remove_particle", pure_virtual(&impl_type::remove_particle))
        .def("get_particle", pure_virtual(&impl_type::get_particle))
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::particle_id_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::particle_id_type const&, typename impl_type::particle_id_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_overlap", pure_virtual((typename impl_type::particle_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&) const)&impl_type::check_overlap), return_value_policy<return_by_value>())
        .def("check_surface_overlap", pure_virtual((typename impl_type::structure_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::position_type const&, typename impl_type::structure_id_type const&, typename impl_type::length_type const&, typename impl_type::structure_id_type const&) const)&impl_type::check_surface_overlap), return_value_policy<return_by_value>())
        .def("check_surface_overlap", pure_virtual((typename impl_type::structure_id_pair_and_distance_list*(impl_type::*)(typename impl_type::particle_type::shape_type const&, typename impl_type::position_type const&, typename impl_type::structure_id_type const&, typename impl_type::length_type const&) const)&impl_type::check_surface_overlap), return_value_policy<return_by_value>())
        .def("create_transaction", pure_virtual(&impl_type::create_transaction),
                return_value_policy<manage_new_object>())
        .def("distance", pure_virtual(&impl_type::distance))
        .def("apply_boundary", pure_virtual((typename impl_type::position_type(impl_type::*)(typename impl_type::position_type const&) const)&impl_type::apply_boundary))
        .def("apply_boundary", pure_virtual((typename impl_type::length_type(impl_type::*)(typename impl_type::length_type const&) const)&impl_type::apply_boundary))
        .def("apply_boundary", pure_virtual((typename impl_type::position_structid_pair_type(impl_type::*)(typename impl_type::position_structid_pair_type const&) const)&impl_type::apply_boundary))
        .def("cyclic_transpose", pure_virtual((typename impl_type::position_type(impl_type::*)(typename impl_type::position_type const&, typename impl_type::position_type const&) const)&impl_type::cyclic_transpose))
        .def("cyclic_transpose", pure_virtual((typename impl_type::length_type(impl_type::*)(typename impl_type::length_type const&, typename impl_type::length_type const&) const)&impl_type::cyclic_transpose))
        .def("cyclic_transpose", pure_virtual((typename impl_type::position_structid_pair_type(impl_type::*)(typename impl_type::position_structid_pair_type const&, typename impl_type::structure_type const&) const)&impl_type::cyclic_transpose))
        .def("__contains__", pure_virtual(&impl_type::has_particle))
        .def("__iter__", pure_virtual(&impl_type::get_particles),
                return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_PARTICLE_CONTAINER_HPP */
