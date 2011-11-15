#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "newBDPropagator.hpp"
#include "binding_common.hpp"

namespace binding {

void register_new_bd_propagator_class()
{
    register_new_bd_propagator_class<newBDPropagator>("newBDPropagator");
}

} // namespace binding
