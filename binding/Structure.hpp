#ifndef BINDING_STRUCTURE_HPP
#define BINDING_STRUCTURE_HPP

#include <boost/python.hpp>
#include "peer/converters/tuple.hpp"

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_structure_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    // registering converters from standard boost templates
    // registers the projected_point/(projected_distance, dist_to_edge) tuple defined in ../Structure.hpp
    peer::converters::register_tuple_converter<typename impl_type::projected_type>();
    peer::converters::register_tuple_converter<typename impl_type::components_pair_type>();
    peer::converters::register_tuple_converter<typename impl_type::position_flag_pair_type>();

    // defining the python class
    return class_<impl_type, boost::shared_ptr<impl_type>,
                  boost::noncopyable>(name, no_init)
        .add_property("id", 
            make_function(&impl_type::id,
                          return_value_policy<return_by_value>()))
        .add_property("sid",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::structure_type_id_type,
                    &impl_type::sid,
                    &impl_type::sid>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::structure_type_id_type,
                &impl_type::sid,
                &impl_type::sid>::set)
        .add_property("structure_id", 
            make_function(&impl_type::structure_id,
                          return_value_policy<return_by_value>()))
        .add_property("name", 
            make_function(&impl_type::name,
                          return_value_policy<return_by_value>()))
        .def("random_position", &impl_type::random_position)
        .def("random_vector", &impl_type::random_vector)
        .def("bd_displacement", &impl_type::bd_displacement)
        .def("project_point", &impl_type::project_point)
        .def("deflect", &impl_type::deflect)
//        .def("deflect_back", &impl_type::deflect_back)
        .def("get_pos_sid_pair", &impl_type::get_pos_sid_pair)
        .def("get_pos_sid_pair_pair", &impl_type::get_pos_sid_pair_pair)
        .def("get_pos_sid_pair_2o", &impl_type::get_pos_sid_pair_2o)
        ;
}

} // namespace binding

#endif /* BINDING_STRUCTURE_HPP */
