#ifndef MOIST_CORE_PREDICATES_INL_
#define MOIST_CORE_PREDICATES_INL_

#include <cmath>
#include <optional>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/vecg.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"

// TODO: THIS IS WRONG!
#define in_range(x, p, eps) (std::abs(x - p) <= eps) || (std::abs(x + p) <= eps)

/**
 * @brief Numeric-Predicates.
 */
namespace moist::predicates
{

    static double round16(double value)
    {
        return value;
        //const double threshold = 1e-3;
        //value = std::round(value * 1e3) / 1e3;
        //return (std::abs(value) < threshold) ? 0.0 : value;
    }

    /**
     * @brief Epsilon-tolerant implementation of orient3d to treat rounded steiner points...
     * Very roughly adapter from geogram::PCK::orient3d().
     *
     * @param p0
     * @param p1
     * @param p2
     * @param p3
     * @return Sign
     */
    PURE INLINE geogram::Sign orient3d(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        double abx = round16(b[0] - a[0]);
        double aby = round16(b[1] - a[1]);
        double abz = round16(b[2] - a[2]);

        double acx = round16(c[0] - a[0]);
        double acy = round16(c[1] - a[1]);
        double acz = round16(c[2] - a[2]);

        double adx = round16(d[0] - a[0]);
        double ady = round16(d[1] - a[1]);
        double adz = round16(d[2] - a[2]);

        // scalar triple product?
        const double volume = abx * (acy * adz - acz * ady)
                    - aby * (acx * adz - acz * adx)
                    + abz * (acx * ady - acy * adx);

        if (std::abs(volume) < 1e-12)
        {
            return geogram::Sign::ZERO;
        }

        return (volume > 0)
            ? geogram::Sign::POSITIVE
            : geogram::Sign::NEGATIVE;
    }

    PURE INLINE bool facet_matches_cell(const g_index cell, const g_index facet, const geogram::Mesh& mesh, const geogram::Mesh& interface)
    {
        vec3 f0 = interface.vertices.point(interface.facets.vertex(facet, 0));
        vec3 f1 = interface.vertices.point(interface.facets.vertex(facet, 1));
        vec3 f2 = interface.vertices.point(interface.facets.vertex(facet, 2));

        // check if three vertices of the cell match the facet
        size_t matches = 0;
        for (g_index lv = 0; lv < mesh.cells.nb_vertices(cell); lv++)
        {
            const vec3 p = mesh.vertices.point(mesh.cells.vertex(cell, lv));
            if (p == f0 || p == f1 || p == f2)
            {
                matches++;
            }
        }

        return matches >= 3;
    }

    PURE INLINE bool vec_eq_2d(const geogram::vec3 &v0, const geogram::vec3 &v1, const moist::AxisAlignedInterfacePlane& plane)
    {
        // TODO [Axis-Support]: Add support for all axis.
        switch (plane.axis)
        {
            case moist::Axis::X:
                return v0.y == v1.y && v0.z == v1.z;
            case moist::Axis::Y:
                return v0.x == v1.x && v0.z == v1.z;
            case moist::Axis::Z:
                return v0.x == v1.x && v0.y == v1.y;
        }

        return false;
    }

    PURE INLINE bool point_on_plane(const geogram::vec3& point, const moist::AxisAlignedInterfacePlane& plane)
    {
        // TODO [Axis-Support]: Add support for all axis.
        switch (plane.axis)
        {
            case moist::Axis::X:
                return in_range(point.x, plane.extent, plane.epsilon);
            case moist::Axis::Y:
                return in_range(point.y, plane.extent, plane.epsilon);
            case moist::Axis::Z:
                return std::abs(point.z - plane.extent) <= plane.epsilon;
        }

        return false;
    }

    PURE INLINE bool edge_on_plane(const geogram::vec3& p0, const geogram::vec3& p1, const moist::AxisAlignedInterfacePlane& plane)
    {
        return point_on_plane(p0, plane) && point_on_plane(p1, plane);
    }

    PURE INLINE bool facet_on_plane(const geogram::index_t facet, const geogram::Mesh& mesh, const moist::AxisAlignedInterfacePlane& plane)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 3
    #endif
        for (l_index lv = 0; lv < /* mesh.facets.nb_vertices(facet) */ 3; lv++)
        {
            const auto point = mesh.vertices.point(mesh.facets.vertex(facet, lv));
            if (!point_on_plane(point, plane))
            {
                return false;
            }
        }

        return true;
    }

    PURE INLINE bool cell_on_plane(const geogram::index_t cell, const geogram::Mesh& mesh, const moist::AxisAlignedInterfacePlane& plane)
    {
    #ifdef OPTION_UNROLL_LOOPS
        #pragma unroll 4
    #endif
        for (unsigned char i = 0; i < /*mesh.cells.nb_vertices(cell)*/ 4; i++)
        {
            if (point_on_plane(mesh.vertices.point(mesh.cells.vertex(cell, i)), plane))
            {
                return true;
            }
        }

        return false;
    }

    PURE INLINE bool cell_facet_on_plane(const geogram::index_t cell, const geogram::index_t lf, const geogram::Mesh& mesh, const moist::AxisAlignedInterfacePlane& plane)
    {
        for (geogram::index_t lv = 0; lv < /*mesh.cells.facet_nb_vertices(cell, lf)*/ 3; lv++)
        {
            const auto point = mesh.vertices.point(mesh.cells.facet_vertex(cell, lf, lv));
            if (!point_on_plane(point, plane))
            {
                return false;
            }
        }

        return true;
    }

    PURE INLINE bool point_of_tet(const geogram::Mesh& mesh, const geogram::index_t t, const geogram::vec3& p)
    {
        const vec3 p0 = mesh.vertices.point(mesh.cells.vertex(t, 0));
        const vec3 p1 = mesh.vertices.point(mesh.cells.vertex(t, 1));
        const vec3 p2 = mesh.vertices.point(mesh.cells.vertex(t, 2));
        const vec3 p3 = mesh.vertices.point(mesh.cells.vertex(t, 3));

        return p == p0 || p == p1 || p == p2 || p == p3;
    }

    enum class PointInTet
    {
        NONE,
        VERTEX,
        EDGE,
        FACET,
        INSIDE
    };

    // Source: geogram/mesh/mesh_AABB.cpp#175
    PURE INLINE PointInTet point_in_tet(const geogram::Mesh& mesh, const geogram::index_t t, const geogram::vec3& p, bool exclude_existing_points = false)
    {
        const auto p0 = mesh.vertices.point(mesh.cells.vertex(t, 0));
        const auto p1 = mesh.vertices.point(mesh.cells.vertex(t, 1));
        const auto p2 = mesh.vertices.point(mesh.cells.vertex(t, 2));
        const auto p3 = mesh.vertices.point(mesh.cells.vertex(t, 3));

        geogram::Sign s[4];
        //s[0] = moist::predicates::orient3d(p, p1, p2, p3);
        //s[1] = moist::predicates::orient3d(p0, p, p2, p3);
        //s[2] = moist::predicates::orient3d(p0, p1, p, p3);
        //s[3] = moist::predicates::orient3d(p0, p1, p2, p);
        s[0] = geogram::PCK::orient_3d(p, p1, p2, p3);
        s[1] = geogram::PCK::orient_3d(p0, p, p2, p3);
        s[2] = geogram::PCK::orient_3d(p0, p1, p, p3);
        s[3] = geogram::PCK::orient_3d(p0, p1, p2, p);
        const bool inside_or_on_boundary = (
            (s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0) ||
            (s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0)
        );

        int nz = 0;
        for(int i = 0; i < 4; i++)
        {
            if(s[i] == 0)
            {
                nz++;
            }
        }

        if (inside_or_on_boundary)
        {
            if (nz == 0)
            {
                return PointInTet::INSIDE;
            }
            if (nz == 1)
            {
                return PointInTet::FACET;
            }
            else if (nz == 2)
            {
                return PointInTet::EDGE;
            }
            else if (nz == 3)
            {
                return PointInTet::VERTEX;
            }
        }

        return PointInTet::NONE;
    }

    namespace xy
    {
        PURE INLINE bool check_lines_aabb(const geogram::vec2& p0, const geogram::vec2& p1, const geogram::vec2& p2, const geogram::vec2& p3)
        {
            return (std::min(p0.x, p1.x) < std::max(p2.x, p3.x)) && (std::max(p0.x, p1.x) > std::min(p2.x, p3.x))
                && (std::min(p0.y, p1.y) < std::max(p2.y, p3.y)) && (std::max(p0.y, p1.y) > std::min(p2.y, p3.y));
        }

        PURE INLINE bool point_in_facet(const geogram::vec2& point, const g_index facet, const geogram::Mesh& triangulation)
        {
            const auto p0 = reinterpret_cast<const geogram::vec2&>(triangulation.vertices.point(triangulation.facets.vertex(facet, 0)));
            const auto p1 = reinterpret_cast<const geogram::vec2&>(triangulation.vertices.point(triangulation.facets.vertex(facet, 1)));
            const auto p2 = reinterpret_cast<const geogram::vec2&>(triangulation.vertices.point(triangulation.facets.vertex(facet, 2)));

            geogram::Sign s[3];
            s[0] = geogram::PCK::orient_2d(p0, p1, point);
            s[1] = geogram::PCK::orient_2d(p1, p2, point);
            s[2] = geogram::PCK::orient_2d(p2, p0, point);

            return (
                (s[0] >= 0 && s[1] >= 0 && s[2] >= 0) ||
                (s[0] <= 0 && s[1] <= 0 && s[2] <= 0)
            );
        }

        PURE /* INLINE */ inline std::optional<geogram::vec2> get_line_intersection(const geogram::vec2& p0, const geogram::vec2& p1, const geogram::vec2& p2, const geogram::vec2& p3, const double rnd = 1e12)
        {
            // prevent line from cutting itself by checking if a point of the line is also a point checked...
            if (p0 == p2 || p1 == p3 || p0 == p3 || p1 == p2)
            {
                return std::nullopt;
            }

            geogram::Sign s[4];
            s[0] = geogram::PCK::orient_2d(p0, p1, p2);
            s[1] = geogram::PCK::orient_2d(p0, p1, p3);
            s[2] = geogram::PCK::orient_2d(p2, p3, p0);
            s[3] = geogram::PCK::orient_2d(p2, p3, p1);

            if (s[0] != s[1] && s[2] != s[3])
            {
                double denom = (p0.x - p1.x) * (p2.y - p3.y) - (p0.y - p1.y) * (p2.x - p3.x);
                if (denom == 0) return std::nullopt;

                double t = ((p0.x - p2.x) * (p2.y - p3.y) - (p0.y - p2.y) * (p2.x - p3.x)) / denom;
                double u = ((p0.x - p2.x) * (p0.y - p1.y) - (p0.y - p2.y) * (p0.x - p1.x)) / denom;

                const auto eps = 2 * moist::__DOUBLE_EPSILON;
                if ((t >= eps && t <= 1.0 - eps) && (u >= eps && u <= 1.0 - eps))
                {
                    return vec2
                    {
                        std::round((p0.x + t * (p1.x - p0.x)) * rnd) / rnd,
                        std::round((p0.y + t * (p1.y - p0.y)) * rnd) / rnd
                    };
                }
            }

            return std::nullopt;
        }
    }

    PURE INLINE bool lines_intersect(const geogram::vec2& p0, const geogram::vec2& p1, const geogram::vec2& p2, const geogram::vec2& p3)
    {
        geogram::Sign s[4];
        s[0] = geogram::PCK::orient_2d(p0, p1, p2);
        s[1] = geogram::PCK::orient_2d(p0, p1, p3);
        s[2] = geogram::PCK::orient_2d(p2, p3, p0);
        s[3] = geogram::PCK::orient_2d(p2, p3, p1);

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

#endif // MOIST_CORE_PREDICATES_INL_
