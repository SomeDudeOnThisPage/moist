#ifndef MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
#define MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/metrics.hpp"

#include "slice.hpp"

namespace moist
{
    bool insert_constraints(moist::MeshSlice& a, moist::MeshSlice& b, moist::Interface& interface, moist::metrics::Metrics_ptr metrics);

    struct CrossedEdge
    {
        g_index e_v0;
        g_index e_v1;
        g_index p;

        g_index e_interface; // (custom) interface edge global index, used in decimation to find a point to decimate onto

        bool operator==(const CrossedEdge& other) const
        {
            return p == other.p && ((e_v0 == other.e_v0 && e_v1 == other.e_v1) || (e_v0 == other.e_v1 && e_v1 == other.e_v0));
        }
    };

    typedef struct
    {
        g_index v0; // alias v_line
        g_index v1;
        g_index v2;
        g_index v3;
    } CreatedTetrahedon;
}


namespace std
{
    template <>
    struct hash<moist::CrossedEdge>
    {
        std::size_t operator()(const moist::CrossedEdge& edge) const
        {
            return std::hash<geo::index_t>()(std::min(edge.e_v0, edge.e_v1)) ^ (std::hash<geo::index_t>()(std::max(edge.e_v0, edge.e_v1)) << 1);
        }
    };
}

#endif // MOIST_INTERFACE_INSERTER_INTERFACE_INSERTER_HPP_
