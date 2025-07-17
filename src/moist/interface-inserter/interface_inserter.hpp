#ifndef MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
#define MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/metrics.hpp"

#include "slice.hpp"

namespace moist
{
    void insert_constraints(moist::MeshSlice& a, moist::MeshSlice& b, moist::Interface& interface, moist::metrics::Metrics_ptr metrics);
}

#endif // MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
