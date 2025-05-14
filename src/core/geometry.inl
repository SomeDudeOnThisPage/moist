#ifndef __OOC_GEOMETRY_CPP
#define __OOC_GEOMETRY_CPP

#include <array>

#include <geogram/mesh/mesh.h>

#include "core.hpp"

#include "ooc_mesh.hpp"
#include "predicates.inl"

namespace moist::geometry
{
    constexpr g_index NO_ELEMENT = -1;
    struct Edge
    {
        g_index v0;
        g_index v1;

        Edge(const g_index v0, const g_index v1) : v0(v0), v1(v1) {};

        bool operator==(const Edge& other) const
        {
            return (v0 == other.v0 && v1 == other.v1) || (v0 == other.v1 && v1 == other.v0);
        }
    };

    struct EdgeHash
    {
        EdgeHash() = default;
        size_t operator()(const Edge& e) const
        {
            return std::hash<g_index>()(e.v0) ^ std::hash<g_index>()(e.v1);
        }
    };

    PURE INLINE bool point_of_cell(const MeshSlice &mesh, const g_index cell, const vec3 &point)
    {
#ifdef OPTION_UNROLL_LOOPS
#pragma unroll 4
#endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            if (mesh.vertices.point(mesh.cells.vertex(cell, lv)) == point)
            {
                return true;
            }
        }

        return false;
    }

    PURE INLINE g_index non_coplanar_opposite(const g_index cell, const g_index a, const g_index b, const MeshSlice &mesh, const AxisAlignedInterfacePlane &plane)
    {
#ifdef OPTION_UNROLL_LOOPS
#pragma unroll 4
#endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (v != a && v != b && !moist::predicates::point_on_plane(mesh.vertices.point(v), plane))
            {
                return v;
            }
        }

        return NO_ELEMENT;
    }

    PURE INLINE bool has_duplicate_vertex(const g_index cell, const MeshSlice &mesh)
    {
        const vec3 a = mesh.vertices.point(mesh.cells.vertex(cell, 0));
        const vec3 b = mesh.vertices.point(mesh.cells.vertex(cell, 1));
        const vec3 c = mesh.vertices.point(mesh.cells.vertex(cell, 2));
        const vec3 d = mesh.vertices.point(mesh.cells.vertex(cell, 3));

        return a == b ||
               a == c ||
               a == d ||
               b == c ||
               b == d ||
               c == d;
    }

    PURE INLINE std::array<g_index, 4> cell_vertices(const g_index cell, const MeshSlice &mesh)
    {
        return {
            mesh.cells.vertex(cell, 0),
            mesh.cells.vertex(cell, 1),
            mesh.cells.vertex(cell, 2),
            mesh.cells.vertex(cell, 3)
        };
    }

    PURE INLINE std::tuple<g_index, g_index, g_index> other(const g_index cell, g_index opposite, const MeshSlice &mesh)
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
                    mesh.cells.vertex(cell, (lv + 3) % 4));
            }
        }

        return std::make_tuple(NO_ELEMENT, NO_ELEMENT, NO_ELEMENT);
    }

    PURE INLINE g_index other(const g_index cell, const g_index a, const g_index b, const g_index c, const MeshSlice &mesh)
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

    PURE INLINE std::tuple<g_index, g_index, g_index> interface_vertices(const g_index cell, const MeshSlice &mesh, const AxisAlignedInterfacePlane &plane)
    {
#ifdef OPTION_UNROLL_LOOPS
#pragma unroll 4
#endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (!moist::predicates::point_on_plane(mesh.vertices.point(v), plane))
            {
                return std::make_tuple(
                    mesh.cells.vertex(cell, (lv + 1) % 4),
                    mesh.cells.vertex(cell, (lv + 2) % 4),
                    mesh.cells.vertex(cell, (lv + 3) % 4));
            }
        }

        return std::make_tuple(NO_ELEMENT, NO_ELEMENT, NO_ELEMENT);
    }

    PURE INLINE g_index non_interface_vertex(const g_index cell, const MeshSlice &mesh, const AxisAlignedInterfacePlane &plane)
    {
#ifdef OPTION_UNROLL_LOOPS
#pragma unroll 4
#endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (!moist::predicates::point_on_plane(mesh.vertices.point(v), plane))
            {
                return v;
            }
        }

        return NO_ELEMENT;
    }

    PURE INLINE std::unordered_set<Edge, EdgeHash> collect_edges(const geogram::Mesh& mesh)
    {
        std::unordered_set<Edge, EdgeHash> edges;
        edges.reserve(mesh.vertices.nb());
        for (const g_index f : mesh.facets)
        {
        #ifdef OPTION_UNROLL_LOOPS
            #pragma unroll 3
        #endif
            for (size_t i = 0; i < 3; i++)
            {
                edges.emplace(
                    mesh.facets.vertex(f, i),
                    mesh.facets.vertex(f, (i + 1) % 3)
                );
            }
        }
        return edges; // return by value not ideal should probably be passed as reference parameter
    }
}

#endif // __OOC_GEOMETRY_CPP
