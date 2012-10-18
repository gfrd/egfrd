#ifndef BINDING_WORLD_HPP
#define BINDING_WORLD_HPP

#include <set>
#include <boost/python.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/const_iterator.hpp>

#include "peer/utils.hpp"
#include "peer/set_indexing_suite.hpp"
#include "utils/range.hpp"
#include "utils/pair.hpp"

namespace binding {

//// Converters. These convert C++ types to something that Python can handle.

// Species_range converter
template<typename Timpl_>
struct species_range_converter: public boost::python::default_call_policies
{
    typedef Timpl_ native_type;

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
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static reference get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0;
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

    struct instance_holder
    {
        instance_holder(native_type const& instance,
                        boost::python::handle<> owner)
            : instance_(instance), owner_(owner) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type instance_;
        boost::python::handle<> owner_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct result_converter
    {
        template<typename T_>
        struct apply
        {
            struct type {
                PyObject* operator()(native_type const& val) const
                {
                    return wrapper_type::create(instance_holder(val, boost::python::handle<>()));
                }

                PyTypeObject const* get_pytype() const
                {
                    return &wrapper_type::__class__;
                }
            };
        };
    };

    template<typename Targs>
    static PyObject* postcall(Targs const& arg, PyObject* result)
    {
        reinterpret_cast<wrapper_type*>(result)->ptr().owner_ = boost::python::handle<>(boost::python::borrowed(PyTuple_GET_ITEM(arg, 0)));
        return result;
    }

    // Register the converter.
    static void __register()
    {
        wrapper_type::__class_init__("SpeciesRange", boost::python::scope().ptr());
    }
};

// Structure_types_range converter
template<typename Timpl_>
struct structure_types_range_converter: public boost::python::default_call_policies
{
    typedef Timpl_ native_type;

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
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static reference get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0;
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

    struct instance_holder
    {
        instance_holder(native_type const& instance,
                        boost::python::handle<> owner)
            : instance_(instance), owner_(owner) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type instance_;
        boost::python::handle<> owner_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct result_converter
    {
        template<typename T_>
        struct apply
        {
            struct type {
                PyObject* operator()(native_type const& val) const
                {
                    return wrapper_type::create(instance_holder(val, boost::python::handle<>()));
                }

                PyTypeObject const* get_pytype() const
                {
                    return &wrapper_type::__class__;
                }
            };
        };
    };

    template<typename Targs>
    static PyObject* postcall(Targs const& arg, PyObject* result)
    {
        reinterpret_cast<wrapper_type*>(result)->ptr().owner_ = boost::python::handle<>(boost::python::borrowed(PyTuple_GET_ITEM(arg, 0)));
        return result;
    }

    // Register the converter.
    static void __register()
    {
        wrapper_type::__class_init__("StructureTypesRange", boost::python::scope().ptr());
    }
};

// Structures_range converter
template<typename Timpl_>
struct structures_range_converter: public boost::python::default_call_policies
{
    typedef Timpl_ native_type;

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
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
        }

        static reference get(T_ const& c, size_type i)
        {
            PyErr_SetString(PyExc_TypeError, "object does not support indexing");
            boost::python::throw_error_already_set();
            throw 0;
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

    struct instance_holder
    {
        instance_holder(native_type const& instance,
                        boost::python::handle<> owner)
            : instance_(instance), owner_(owner) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type instance_;
        boost::python::handle<> owner_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder, peer::wrappers::default_policy_generator<policy> > wrapper_type;

    struct result_converter
    {
        template<typename T_>
        struct apply
        {
            struct type {
                PyObject* operator()(native_type const& val) const
                {
                    return wrapper_type::create(instance_holder(val, boost::python::handle<>()));
                }

                PyTypeObject const* get_pytype() const
                {
                    return &wrapper_type::__class__;
                }
            };
        };
    };

    template<typename Targs>
    static PyObject* postcall(Targs const& arg, PyObject* result)
    {
        reinterpret_cast<wrapper_type*>(result)->ptr().owner_ = boost::python::handle<>(boost::python::borrowed(PyTuple_GET_ITEM(arg, 0)));
        return result;
    }

    // Register the converter.
    static void __register()
    {
        wrapper_type::__class_init__("StructuresRange", boost::python::scope().ptr());
    }
};

template<typename T>
static typename get_select_first_range<typename T::particle_id_pair_range>::type
World_get_particle_ids(T const& world)
{
    return make_select_first_range(world.get_particles_range());
}



////// Registering master function
template<typename Timpl_, typename Tbase_, typename Tsim>
inline boost::python::objects::class_base register_world_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;       // The class World that is being wrapped.

    typedef species_range_converter<typename impl_type::species_range> species_range_converter_type;
    typedef structure_types_range_converter<typename impl_type::structure_types_range> structure_types_range_converter_type;
    typedef structures_range_converter<typename impl_type::structures_range> structures_range_converter_type;

    // registering the converters defined above
    species_range_converter_type::__register();
    structure_types_range_converter_type::__register();
    structures_range_converter_type::__register();

    // registering converters from standard boost templates
    class_<std::set<typename impl_type::particle_id_type> >("ParticleIDSet")
        .def(peer::util::set_indexing_suite<std::set<typename impl_type::particle_id_type> >())
        ;

    class_<std::set<typename impl_type::structure_id_type> >("StructureIDSet")
        .def(peer::util::set_indexing_suite<std::set<typename impl_type::structure_id_type> >())
        ;

    // defining the python class
    return class_<impl_type, bases<Tbase_>,
                  boost::shared_ptr<impl_type> >(
        "World", init<typename impl_type::length_type, typename impl_type::size_type>())
        .add_property("cell_size", &impl_type::cell_size)
        .add_property("matrix_size", &impl_type::matrix_size)
        .add_property("species",
            make_function(
                (typename impl_type::species_range(impl_type::*)() const)&impl_type::get_species, species_range_converter_type()))
        .add_property("structures",
            make_function(
                (typename impl_type::structures_range(impl_type::*)() const)&impl_type::get_structures, structures_range_converter_type()))
        .add_property("structure_types",
            make_function(
                (typename impl_type::structure_types_range(impl_type::*)() const)&impl_type::get_structure_types,
                 structure_types_range_converter_type()))
        .add_property("particle_ids", &World_get_particle_ids<impl_type>)
        .def("get_particle_ids", &impl_type::get_particle_ids)
        .def("get_particle_ids_on_struct", &impl_type::get_particle_ids_on_struct)
        .def("add_species", &impl_type::add_species)
        // Structure stuff
        .def("add_structure", &impl_type::template add_structure<typename impl_type::planar_surface_type>)
        .def("add_structure", &impl_type::template add_structure<typename impl_type::cuboidal_region_type>)
        .def("add_structure", &impl_type::template add_structure<typename impl_type::cylindrical_surface_type>)
        .def("add_structure", &impl_type::template add_structure<typename impl_type::disk_surface_type>)
        .def("add_structure", &impl_type::template add_structure<typename impl_type::spherical_surface_type>)
        .def("connect_structures", &impl_type::template connect_structures<typename impl_type::planar_surface_type,      int>)
        .def("connect_structures", &impl_type::template connect_structures<typename impl_type::cylindrical_surface_type, int>)
        .def("get_neighbor_info",  &impl_type::template get_neighbor_info<typename impl_type::planar_surface_type,      int>)
        .def("get_neighbor_info",  &impl_type::template get_neighbor_info<typename impl_type::cylindrical_surface_type, int>)
              // TODO get_neighbor_info still needs a type converter to work properly!
        .def("get_neighbor_id",    &impl_type::template get_neighbor_id<typename impl_type::planar_surface_type,      int>)
        .def("get_neighbor_id",    &impl_type::template get_neighbor_id<typename impl_type::cylindrical_surface_type, int>)
        .def("set_def_structure",  &impl_type::set_def_structure)
        // StructureType stuff
        .def("add_structure_type", &impl_type::add_structure_type)
        .def("set_def_structure_type_id", &impl_type::set_def_structure_type_id)
        .def("distance", &impl_type::template distance<typename impl_type::position_type>)
        .def("distance", &impl_type::template distance<typename Tsim::sphere_type>)
        .def("distance", &impl_type::template distance<typename Tsim::cylinder_type>)
        .def("distance", &impl_type::template distance<typename Tsim::disk_type>)
        .def("distance", &impl_type::template distance<typename Tsim::box_type>)
        .def("distance", &impl_type::template distance<typename Tsim::plane_type>)
        .def("calculate_pair_CoM", &impl_type::template calculate_pair_CoM<typename impl_type::position_type>)
        ;
}

} // namespace binding

#endif /* BINDING_WORLD_HPP */
