#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Disk.hpp"
#include "binding_common.hpp"

namespace binding {

void register_disk_class()
{
    register_disk_class<Disk>("Disk");
}

} // namespace binding
