#ifndef __OOC_GEOMETRY_CPP
#define __OOC_GEOMETRY_CPP

#include <geogram/mesh/mesh.h>

#include "../core.hpp"

#include "ooc_mesh.hpp"
#include "predicates.inl"


namespace incremental_meshing::geometry
{
    constexpr g_index NO_ELEMENT = -1;

    PURE INLINE g_index non_coplanar_opposite(const g_index cell, const g_index a, const g_index b, const SubMesh& mesh, const AxisAlignedInterfacePlane& plane)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (v != a && v != b && !incremental_meshing::predicates::point_on_plane(mesh.vertices.point(v), plane))
            {
                return v;
            }
        }

        return NO_ELEMENT;
    }

    PURE INLINE g_index opposite(const g_index cell, const g_index a, const g_index b, const g_index c, const SubMesh& mesh)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (v != a && v != b && v != c)
            {
                return v;
            }
        }

        return NO_ELEMENT;
    }
}

#endif // __OOC_GEOMETRY_CPP
