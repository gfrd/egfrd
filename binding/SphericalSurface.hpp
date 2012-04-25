#ifndef BINDING_SPHERICAL_SURFACE_HPP
#define BINDING_SPHERICAL_SURFACE_HPP

#include <boost/python.hpp>

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_spherical_surface_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // defining the python class
    return class_<impl_type, bases<typename impl_type::base_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(
            name, init<typename impl_type::structure_name_type,
                       typename impl_type::structure_type_id_type,
                       typename impl_type::structure_id_type,
                       typename impl_type::shape_type>())
        .add_property("shape",
            make_function((typename impl_type::shape_type const&(impl_type::*)()const)&impl_type::shape, return_value_policy<return_by_value>()))
        ;
}

} // namespace binding

#endif /* BINDING_SPHERICAL_SURFACE_HPP */
