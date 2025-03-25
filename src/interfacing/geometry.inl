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

    PURE INLINE std::tuple<g_index, g_index, g_index> other(const g_index cell, g_index opposite, const SubMesh& mesh)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (opposite == v)
            {
                return std::make_tuple(
                    mesh.cells.vertex(cell, (lv + 1) % 4),
                    mesh.cells.vertex(cell, (lv + 2) % 4),
                    mesh.cells.vertex(cell, (lv + 3) % 4)
                );
            }
        }

        return std::make_tuple(NO_ELEMENT, NO_ELEMENT, NO_ELEMENT);
    }

    PURE INLINE g_index other(const g_index cell, const g_index a, const g_index b, const g_index c, const SubMesh& mesh)
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

    PURE INLINE g_index non_interface_vertex(const g_index cell, const SubMesh& mesh, const AxisAlignedInterfacePlane& plane)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (!incremental_meshing::predicates::point_on_plane(mesh.vertices.point(v), plane))
            {
                return v;
            }
        }

        return NO_ELEMENT;
    }
}

#endif // __OOC_GEOMETRY_CPP
