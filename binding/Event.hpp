#ifndef BINDING_PYEVENT_HPP
#define BINDING_PYEVENT_HPP

#include <boost/python.hpp>

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::object register_event_class(char const* name)
{
    using namespace boost::python;

    return class_<Timpl, bases<>, boost::shared_ptr<Timpl>, boost::noncopyable>(name, init<typename Timpl::time_type>())
        .add_property("time",
            make_function(
                &Timpl::time, return_value_policy<copy_const_reference>()))
        ;
}

} // namespace binding

#endif /* BINDING_PYEVENT_HPP */
