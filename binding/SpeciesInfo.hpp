#ifndef BINDING_SPECIES_INFO_HPP
#define BINDING_SPECIES_INFO_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {


////// Registering master function
template<typename Timpl_>
inline boost::python::objects::class_base register_species_info_class(
        char const* name)
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    // defining the python class
    return class_<impl_type>("SpeciesInfo",
            init<>() )
        .def(init<typename impl_type::identifier_type, typename impl_type::structure_type_id_type,
                  typename impl_type::length_type, typename impl_type::D_type, typename impl_type::v_type>())
        .add_property("id",
            make_function(&impl_type::id,
                return_value_policy<return_by_value>()))
        .add_property("structure_type_id",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::structure_type_id_type,
                    &impl_type::structure_type_id,
                    &impl_type::structure_type_id>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::structure_type_id_type,
                &impl_type::structure_type_id,
                &impl_type::structure_type_id>::set)
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::length_type,
                    &impl_type::radius,
                    &impl_type::radius>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::length_type,
                &impl_type::radius,
                &impl_type::radius>::set)
        .add_property("D",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::D_type,
                    &impl_type::D,
                    &impl_type::D>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::D_type,
                &impl_type::D,
                &impl_type::D>::set)
         .add_property("v",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type, typename impl_type::v_type,
                    &impl_type::v,
                    &impl_type::v>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                impl_type, typename impl_type::v_type,
                &impl_type::v,
                &impl_type::v>::set)
        ;

}

} // namespace peer

#endif /* BINDING_SPECIES_INFO_HPP */
