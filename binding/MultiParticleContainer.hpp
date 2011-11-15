#ifndef BINDING_MULTI_PARTICLE_CONTAINER_HPP
#define BINDING_MULTI_PARTICLE_CONTAINER_HPP

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

#include "peer/converters/tuple.hpp"

namespace binding {

template<typename Timpl, typename Tworld>
boost::python::objects::class_base register_multi_particle_container_class(char const* name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    peer::converters::register_tuple_converter<typename impl_type::real_pair>();

    return class_<impl_type, bases<typename Tworld::particle_container_type>,
           boost::noncopyable>(name, init<Tworld&>())
        .def("determine_dt_and_reaction_length",&impl_type::determine_dt_and_reaction_length)
        ;
}

} // namespace binding

#endif /* BINDING_MULTI_PARTICLE_CONTAINER_HPP */
