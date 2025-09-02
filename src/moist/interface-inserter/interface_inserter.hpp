#ifndef MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
#define MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/metrics.hpp"

#include "exact_types.hpp"

namespace moist
{
    struct CrossedEdgeExact
    {
        std::size_t v0;
        std::size_t v1;

        moist::exact::Point p;
        std::size_t vp;
    };

    struct RemeshingParameters
    {
        double hmin;
        double hmax;
    };

    void create_interface_mesh(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane, const moist::RemeshingParameters remeshing, moist::metrics::Metrics_ptr metrics = nullptr);
}

#endif // MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
