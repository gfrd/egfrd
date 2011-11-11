#ifndef BINDING_STRUCTURE_HPP
#define BINDING_STRUCTURE_HPP

#include <boost/python.hpp>
#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl>
inline boost::python::objects::class_base register_structure_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    peer::converters::register_tuple_converter<
            typename impl_type::projected_type>();

    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, no_init)
        .add_property("id", 
            make_function(&impl_type::id,
                          return_value_policy<return_by_value>()))
        .add_property("sid",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::species_id_type,
                    &impl_type::sid,
                    &impl_type::sid>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::species_id_type,
                &impl_type::sid,
                &impl_type::sid>::set)
        .def("random_position", &impl_type::random_position)
        .def("random_vector", &impl_type::random_vector)
        .def("bd_displacement", &impl_type::bd_displacement)
        .def("projected_point", &impl_type::projected_point)
        ;
}

} // namespace binding

#endif /* BINDING_STRUCTURE_HPP */
