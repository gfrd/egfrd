#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Identifier.hpp"
#include "SerialIDGenerator.hpp"
#include "binding_common.hpp"
#include "peer/utils.hpp"

namespace binding {

void register_structure_id_class()
{
    IdentifierWrapper<StructureID>::__register_class("StructureID");
    register_serial_id_generator_class<StructureID>("StructureIDGenerator");
}

} // namespace binding
