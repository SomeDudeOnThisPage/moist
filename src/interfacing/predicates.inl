#ifndef __OOC_PREDICATES_CPP
#define __OOC_PREDICATES_CPP

#include <cmath>
#include <optional>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/vecg.h>

#include "../core.hpp"
#include "interface.hpp"

// TODO: THIS IS WRONG!
#define in_range(x, p, eps) (std::abs(x - p) <= eps) || (std::abs(x + p) <= eps)

#define orient_3d geogram::PCK::orient_3d
#define orient_2d geogram::PCK::orient_2d
#define EPS 1e-8

namespace incremental_meshing::predicates
{
    PURE INLINE bool vec_eq_2d(const geogram::vec3 &v0, const geogram::vec3 &v1, const incremental_meshing::AxisAlignedInterfacePlane& plane)
    {
        switch (plane.axis)
        {
            case incremental_meshing::Axis::X:
                return v0.y == v1.y && v0.z == v1.z;
            case incremental_meshing::Axis::Y:
                return v0.x == v1.x && v0.z == v1.z;
            case incremental_meshing::Axis::Z:
                return v0.x == v1.x && v0.y == v1.y;
        }

        return false;
    }

    PURE INLINE bool point_on_plane(const geogram::vec3& point, const incremental_meshing::AxisAlignedInterfacePlane& plane)
    {
        switch (plane.axis)
        {
            case incremental_meshing::Axis::X:
                return in_range(point.x, plane.extent, plane.epsilon);
            case incremental_meshing::Axis::Y:
                return in_range(point.y, plane.extent, plane.epsilon);
            case incremental_meshing::Axis::Z:
                return point.z == plane.extent; // TODO
                //return in_range(point.z, plane.extent, plane.epsilon);
        }

        return false;
    }

    PURE INLINE bool facet_on_plane(const geogram::index_t facet, const geogram::Mesh& mesh, const incremental_meshing::AxisAlignedInterfacePlane& plane)
    {
        for (geogram::index_t lv = 0; lv < /* mesh.facets.nb_vertices(facet) */ 3; lv++)
        {
            const auto point = mesh.vertices.point(mesh.facets.vertex(facet, lv));
            if (!point_on_plane(point, plane))
            {
                return false;
            }
        }

        return true;
    }

    PURE INLINE bool cell_on_plane(const geogram::index_t cell, const geogram::Mesh& mesh, const incremental_meshing::AxisAlignedInterfacePlane& plane)
    {
        unsigned char points_on_plane = 0;
        for (unsigned char i = 0; i < mesh.cells.nb_vertices(cell); i++)
        {
            if (point_on_plane(mesh.vertices.point(mesh.cells.vertex(cell, i)), plane))
            {
                points_on_plane++;
            }
        }

        return points_on_plane == 3;
    }

    PURE INLINE bool cell_facet_on_plane(const geogram::index_t cell, const geogram::index_t lf, const geogram::Mesh& mesh, const incremental_meshing::AxisAlignedInterfacePlane& plane)
    {
        for (geogram::index_t lv = 0; lv < mesh.cells.facet_nb_vertices(cell, lf); lv++)
        {
            const auto point = mesh.vertices.point(mesh.cells.facet_vertex(cell, lf, lv));
            if (!point_on_plane(point, plane))
            {
                return false;
            }
        }

        return true;
    }

    // Source: geogram/mesh/mesh_AABB.cpp#175
    // TODO: this can likely be replaced with less instructions because we work in 2d...
    PURE INLINE bool point_in_tet(const geogram::Mesh& mesh, geogram::index_t t, const geogram::vec3& p)
    {
        const auto p0 = mesh.vertices.point(mesh.cells.vertex(t, 0));
        const auto p1 = mesh.vertices.point(mesh.cells.vertex(t, 1));
        const auto p2 = mesh.vertices.point(mesh.cells.vertex(t, 2));
        const auto p3 = mesh.vertices.point(mesh.cells.vertex(t, 3));

        geogram::Sign s[4];
        s[0] = orient_3d(p, p1, p2, p3);
        s[1] = orient_3d(p0, p, p2, p3);
        s[2] = orient_3d(p0, p1, p, p3);
        s[3] = orient_3d(p0, p1, p2, p);

        return (
            (s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0) ||
            (s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0)
        );
    }

    /*inline bool lines_intersect(const geogram::vec3& p0, const geogram::vec3& p1, const geogram::vec3& p2, const geogram::vec3& p3)
    {
        const auto d0 = p1 - p0;
        const auto d1 = p3 - p2;

        const auto n = geogram::cross(d0, d1);
        if (n.length() < EPS)
        {
            return false;
        }

        const auto dot = geogram::dot(n, p3 - p1);
        return EPS > dot && dot > -EPS;
    }*/

    PURE /* INLINE */ inline std::optional<geogram::vec2> get_line_intersection(const geogram::vec2& p0, const geogram::vec2& p1, const geogram::vec2& p2, const geogram::vec2& p3)
    {
        geogram::Sign s[4];
        s[0] = orient_2d(p0, p1, p2);
        s[1] = orient_2d(p0, p1, p3);
        s[2] = orient_2d(p2, p3, p0);
        s[3] = orient_2d(p2, p3, p1);

        if (s[0] != s[1] && s[2] != s[3])
        {
            double denom = (p0.x - p1.x) * (p2.y - p3.y) - (p0.y - p1.y) * (p2.x - p3.x);
            if (denom == 0) return std::nullopt;

            double t = ((p0.x - p2.x) * (p2.y - p3.y) - (p0.y - p2.y) * (p2.x - p3.x)) / denom;
            double u = ((p0.x - p2.x) * (p0.y - p1.y) - (p0.y - p2.y) * (p0.x - p1.x)) / denom;

            if ((t > 0 && t < 1) && (u > 0 && u < 1))
            {
                geogram::vec2 intersection = {
                    p0.x + t * (p1.x - p0.x),
                    p0.y + t * (p1.y - p0.y)
                };

                return intersection;
            }
        }

        return std::nullopt;
    }

    PURE INLINE bool lines_intersect(const geogram::vec2& p0, const geogram::vec2& p1, const geogram::vec2& p2, const geogram::vec2& p3)
    {
        geogram::Sign s[4];
        s[0] = orient_2d(p0, p1, p2);
        s[1] = orient_2d(p0, p1, p3);
        s[2] = orient_2d(p2, p3, p0);
        s[3] = orient_2d(p2, p3, p1);

        if (s[0] != s[1] && s[2] != s[3])
        {
            double denom = (p0.x - p1.x) * (p2.y - p3.y) - (p0.y - p1.y) * (p2.x - p3.x);
            if (denom == 0) return false;

            double t = ((p0.x - p2.x) * (p2.y - p3.y) - (p0.y - p2.y) * (p2.x - p3.x)) / denom;
            double u = ((p0.x - p2.x) * (p0.y - p1.y) - (p0.y - p2.y) * (p0.x - p1.x)) / denom;

            return (t > 0 && t < 1) && (u > 0 && u < 1);
        }

        return false;
    }
}

#undef orient_3d
#undef orient_2d

#endif // __OOC_PREDICATES_CPP
