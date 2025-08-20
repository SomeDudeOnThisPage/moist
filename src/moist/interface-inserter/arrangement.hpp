#ifndef MOIST_INTERFACE_INSERTER_ARRANGEMENT_HPP_
#define MOIST_INTERFACE_INSERTER_ARRANGEMENT_HPP_

#include <vector>
#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

#include "exact_types.hpp"

namespace moist::exact
{
    using Triangulation = std::vector<moist::exact::Triangle>;

    enum class ArrangementConfiguration
    {
        OVERLAP,
        NO_OVERLAP,
        // extra condition, in this case we merge the two closest vertices instead of creating new tets
        // this is, because it happens where envelopes overlap, and a tet from B with vertices completely outside
        // mesh A intersects with a tet from mesh A.
        // in this case it's better to move vertices instead of creating almost degenerate tetrahedra, which will
        // likely become degenerate when discretizing them to double (or worse, float, when saving as .mesh).
        OVERLAP_NO_SHARED_VERTEX
    };

    bool arrangeable(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1);
    Triangulation arrange(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1);
}

#endif // MOIST_INTERFACE_INSERTER_ARRANGEMENT_HPP_
