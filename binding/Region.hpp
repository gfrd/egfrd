#ifndef BINDING_REGION_HPP
#define BINDING_REGION_HPP

#include <boost/python.hpp>

namespace binding {

////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_region_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // defining the python class (kinda empty)
    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(
            name, no_init);
}

} // namespace binding

#endif /* BINDING_REGION_HPP */
