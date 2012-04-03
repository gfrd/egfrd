#ifndef BINDING_SURFACE_HPP
#define BINDING_SURFACE_HPP

#include <boost/python.hpp>

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_surface_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // defining the python class
    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(
            name, no_init)
        .def("minimal_distance", &impl_type::minimal_distance)
        ;
}

} // namespace binding

#endif /* BINDING_SURFACE_HPP */
