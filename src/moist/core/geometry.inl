#ifndef MOIST_CORE_GEOMETRY_INL_
#define MOIST_CORE_GEOMETRY_INL_

#include <array>
#include <unordered_set>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/predicates.inl"
#include "moist/core/attributes.inl"

/**
 * @brief Geometric-Predicates.
 */
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

    PURE INLINE bool vertex_of_cell(const geo::Mesh& mesh, const g_index cell, const g_index v)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            if (mesh.cells.vertex(cell, lv) == v)
            {
                return true;
            }
        }

        return false;
    }

    PURE INLINE bool point_of_cell(const geo::Mesh &mesh, const g_index cell, const vec3 &point)
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

    PURE INLINE g_index non_coplanar_opposite(const g_index cell, const g_index a, const g_index b, const geo::Mesh &mesh, const AxisAlignedPlane &plane)
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

    PURE INLINE bool has_duplicate_vertex(const g_index cell, const geo::Mesh &mesh)
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

    PURE INLINE std::array<g_index, 4> cell_vertices(const g_index cell, const geo::Mesh &mesh)
    {
        return {
            mesh.cells.vertex(cell, 0),
            mesh.cells.vertex(cell, 1),
            mesh.cells.vertex(cell, 2),
            mesh.cells.vertex(cell, 3)
        };
    }

    PURE INLINE std::array<g_index, 3> other(const g_index cell, g_index opposite, const geo::Mesh &mesh)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (opposite == v)
            {
                return std::array{
                    mesh.cells.vertex(cell, (lv + 1) % 4),
                    mesh.cells.vertex(cell, (lv + 2) % 4),
                    mesh.cells.vertex(cell, (lv + 3) % 4)
                };
            }
        }

        return std::array{NO_ELEMENT, NO_ELEMENT, NO_ELEMENT};
    }

    PURE INLINE g_index other(const g_index cell, const g_index a, const g_index b, const g_index c, const geo::Mesh &mesh)
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

    PURE INLINE std::vector<geo::index_t> interface_vertices(const geo::index_t& cell, const geo::Mesh &mesh)
    {
        std::vector<geo::index_t> vertices;
        const auto v_interface = geo::Attribute<bool>(mesh.vertices.attributes(), moist::attributes::V_INTERFACE);

        #pragma unroll 4
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(cell, lv);
            if (v_interface[v])
            {
                vertices.push_back(v);
            }
        }

        return vertices;
    }

    PURE INLINE std::tuple<g_index, g_index, g_index> interface_vertices(const g_index cell, const geo::Mesh &mesh, const AxisAlignedPlane &plane)
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

    PURE INLINE g_index non_interface_vertex(const g_index c, const geo::Mesh& mesh)
    {
        const auto v_interface = geo::Attribute<bool>(mesh.vertices.attributes(), moist::attributes::V_INTERFACE);

        #ifdef OPTION_UNROLL_LOOPS
            #pragma unroll 4
        #endif
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = mesh.cells.vertex(c, lv);
            LOCK_ATTRIBUTES;
            if (!v_interface[v])
            {
                return v;
            }
        }
    };

    PURE INLINE g_index non_interface_vertex(const g_index cell, const geo::Mesh &mesh, const AxisAlignedPlane &plane)
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

    PURE INLINE std::vector<g_index> non_interface_vertices(const g_index c, const geo::Mesh &mesh)
    {
        std::vector<g_index> indices;

        const auto v_interface = geo::Attribute<bool>(mesh.vertices.attributes(), moist::attributes::V_INTERFACE);
        for (const g_index v : moist::geometry::cell_vertices(c, mesh))
        {
            if (!v_interface[v])
            {
                indices.push_back(v);
            }
        }

        return indices;
    }

    PURE INLINE std::vector<g_index> other(const g_index c, const std::vector<g_index>& v_other, const geo::Mesh &mesh)
    {
        std::vector<g_index> indices;

        for (const g_index v : moist::geometry::cell_vertices(c, mesh))
        {
            if (std::find(v_other.begin(), v_other.end(), v) == v_other.end())
            {
                indices.push_back(v);
            }
        }

        return indices;
    }

    PURE INLINE std::unordered_set<Edge, EdgeHash> collect_edges(const geo::Mesh& mesh, const AxisAlignedPlane& plane)
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
                g_index v0 = mesh.facets.vertex(f, i);
                g_index v1 = mesh.facets.vertex(f, (i + 1) % 3);

                if (v0 > v1)
                {
                    std::swap(v0, v1);
                }

                const vec3 p0 = mesh.vertices.point(v0);
                const vec3 p1 = mesh.vertices.point(v1);
                if (moist::predicates::point_on_plane(p0, plane) && moist::predicates::point_on_plane(p1, plane))
                {
                    edges.emplace(
                        v0,
                        v1
                    );
                }
            }
        }
        return edges; // return by value not ideal should probably be passed as reference parameter
    }

    PURE INLINE std::unordered_set<Edge, EdgeHash> collect_interface_edges(const geo::Mesh& mesh, const moist::AxisAlignedPlane& plane)
    {
        const auto nv = mesh.cells.nb();
        std::unordered_set<Edge, EdgeHash> edges;
        for (const g_index c : mesh.cells)
        {
            for (g_index le = 0; le < mesh.cells.nb_edges(c); le++)
            {
                const g_index v0 = mesh.cells.edge_vertex(c, le, 0);
                const g_index v1 = mesh.cells.edge_vertex(c, le, 1);
                const vec3 cp0 = mesh.vertices.point(v0);
                const vec3 cp1 = mesh.vertices.point(v1);

                if (!moist::predicates::edge_on_plane(cp0, cp1, plane))
                {
                    continue;
                }

                edges.emplace(
                    mesh.cells.edge_vertex(c, le, 0),
                    mesh.cells.edge_vertex(c, le, 1)
                );
            }
        }

        return edges;
    }

    PURE INLINE std::unordered_set<Edge, EdgeHash> collect_edges(const geo::Mesh& mesh)
    {
        std::unordered_set<Edge, EdgeHash> edges;
        if (mesh.facets.nb() > 0)
        {
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
        }
        else if (mesh.edges.nb() > 0)
        {
            for (const geo::index_t e : mesh.edges)
            {
                edges.emplace(mesh.edges.vertex(e, 0), mesh.edges.vertex(e, 1));
            }
        }

        return edges; // return by value not ideal should probably be passed as reference parameter
    }
}

#endif // MOIST_CORE_GEOMETRY_INL_
