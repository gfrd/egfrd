#ifndef BINDING_CYLINDRICAL_PAIR_HPP
#define BINDING_CYLINDRICAL_PAIR_HPP

#include <boost/python.hpp>
#include "fcpair_reference_accessor_wrapper.hpp"

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_cylindrical_pair_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(name,
        init<typename impl_type::identifier_type,
             typename impl_type::particle_id_pair,
             typename impl_type::particle_id_pair,
             typename impl_type::shell_id_pair,
             typename impl_type::position_type,
             typename impl_type::reaction_rule_vector>())
        .add_property("shell",
            make_function(
                 &fcpair_reference_accessor_wrapper<
                    impl_type, typename impl_type::shell_id_pair,
                    &impl_type::shell, &impl_type::shell>::get,
                 return_value_policy<return_by_value>()),
            make_function(
                 &fcpair_reference_accessor_wrapper<
                    impl_type, typename impl_type::shell_id_pair,
                    &impl_type::shell, &impl_type::shell>::set))
        .add_property("iv",
            make_function(&impl_type::iv,
                return_value_policy<return_by_value>()))
        .add_property("r0",
            make_function(&impl_type::r0,
                return_value_policy<return_by_value>()))
        .add_property("reactions",
            make_function(&impl_type::reactions,
                return_value_policy<return_by_value>()))
        .add_property("a_R",
            make_function(&impl_type::a_R,
                return_value_policy<return_by_value>()))
        .add_property("a_r",
            make_function(&impl_type::a_r,
                return_value_policy<return_by_value>()))
        .add_property("sigma",
            make_function(&impl_type::sigma,
                return_value_policy<return_by_value>()))
        .add_property("D_tot",
            make_function(&impl_type::D_tot,
                return_value_policy<return_by_value>()))
        .add_property("D_geom",
            make_function(&impl_type::D_geom,
                return_value_policy<return_by_value>()))
        .add_property("D_R",
            make_function(&impl_type::D_R,
                return_value_policy<return_by_value>()));
}

} // namespace binding

#endif /* BINDING_CYLINDRICAL_PAIR_HPP */
