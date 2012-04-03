#ifndef BINDING_BOX_HPP
#define BINDING_BOX_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"
#include "peer/numpy/ndarray_converters.hpp"
#include "peer/converters/sequence.hpp"

namespace binding {

// Not sure
template<typename Timpl_>
static std::string Box___str__(Timpl_* impl)
{
    return boost::lexical_cast<std::string>(*impl);
}


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_box_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // registering converters from standard boost templates
    to_python_converter<boost::array<typename impl_type::length_type, 3>,
            peer::util::detail::to_ndarray_converter<
                boost::array<typename impl_type::length_type, 3> > >();
    peer::converters::register_iterable_to_ra_container_converter<
            boost::array<typename impl_type::length_type, 3>, 3>();

    // defining the python class
    return class_<impl_type>(name)
        .def(init<typename impl_type::position_type, 
                  typename impl_type::position_type, 
                  typename impl_type::position_type, 
                  typename impl_type::position_type, 
                  typename impl_type::length_type,
                  typename impl_type::length_type, 
                  typename impl_type::length_type>())
        .def(init<typename impl_type::position_type, 
                  boost::array<typename impl_type::length_type, 3> >())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position,
                    &impl_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::position,
                    &impl_type::position>::set))
        .add_property("half_extent",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    boost::array<typename impl_type::length_type, 3>,
                    &impl_type::half_extent,
                    &impl_type::half_extent>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    boost::array<typename impl_type::length_type, 3>,
                    &impl_type::half_extent,
                    &impl_type::half_extent>::set))
        .add_property("unit_x",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_x,
                    &impl_type::unit_x>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_x,
                    &impl_type::unit_x>::set))
        .add_property("unit_y",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_y,
                    &impl_type::unit_y>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_y,
                    &impl_type::unit_y>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_z,
                    &impl_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    typename impl_type::position_type,
                    &impl_type::unit_z,
                    &impl_type::unit_z>::set))
        .def("__str__", &Box___str__<impl_type>)
        .def("show", &impl_type::show);
}

} // namespace binding

#endif /* BINDING_BOX_HPP */
