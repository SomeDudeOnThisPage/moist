#ifndef MOIST_CORE_GEOMETRY_EXACT_INL_
#define MOIST_CORE_GEOMETRY_EXACT_INL_

#include <array>
#include <algorithm>
#include <unordered_set>
#include <optional>

#include <CGAL/intersections.h>
#include <boost/variant.hpp>

#include "moist/core/defines.hpp"

#include "exact_types.hpp"
#include "exact_mesh.hpp"

namespace moist::geometry::exact
{
    PURE INLINE moist::exact::CellType get_cell_type(std::array<std::size_t, 4>& vertices, const moist::ExactMesh& mesh)
    {
        std::size_t nv = 0;
        for (const auto v : vertices)
        {
            nv += (mesh.Point(v)._interface) ? 1 : 0;
        }

        switch(nv)
        {
        case 1:
            return moist::exact::CellType::I;
        case 2:
            return moist::exact::CellType::II;
        case 3:
            return moist::exact::CellType::III;
        default:
            OOC_ERROR("invalid interface cell type");
        }
        return moist::exact::CellType::I; // stop compiler complain
    };

    PURE INLINE bool cell_has_interface_facet(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
    {
        std::size_t nv = 0;
        for (const auto v : cell._points)
        {
            if (mesh.Point(v)._interface)
            {
                nv++;
            }
        }

        return nv == 3;
    }

    PURE INLINE std::vector<std::size_t> interface_vertices(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
    {
        std::vector<std::size_t> vertices;
        for (const auto v : cell._points)
        {
            if (mesh.Point(v)._interface)
            {
                vertices.push_back(v);
            }
        }
        return vertices;
    }

    PURE INLINE std::array<std::size_t, 3> other(const std::size_t& c, const std::size_t& v_opposite, const moist::ExactMesh& mesh)
    {
        const auto& cell = mesh.Cell(c);
        for (std::size_t lv = 0; lv < 4; lv++)
        {
            const auto& v = cell._points[lv];
            if (v_opposite == v)
            {
                return std::array
                {
                    cell._points[(lv + 1) % 4],
                    cell._points[(lv + 2) % 4],
                    cell._points[(lv + 3) % 4]
                };
            }
        }

        return std::array
        {
            std::size_t(-1U), std::size_t(-1U), std::size_t(-1U)
        };
    }

    PURE INLINE bool is_degenerate(const moist::exact::Cell& cell)
    {
        std::unordered_set<std::size_t> points;
        for (const auto v : cell._points)
        {
            if (points.contains(v))
            {
                return true;
            }
            points.insert(v);
        }

        return false;
    }

    static bool is_strictly_inside(const moist::exact::Kernel::Segment_2& seg, const moist::exact::Kernel::Point_2& p)
    {
        return CGAL::collinear_are_ordered_along_line(seg.source(), p, seg.target()) &&
            p != seg.source() &&
            p != seg.target();
    }

    PURE INLINE std::optional<std::vector<moist::exact::Triangle>> arrange()
    {
        std::vector<moist::exact::Triangle> polygon_points;
        return std::nullopt;
    }

    PURE INLINE std::optional<std::vector<moist::exact::Point>> get_intersection(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1)
    {
        std::vector<moist::exact::Point> polygon_points;

        auto result = CGAL::intersection(t0._t, t1._t);
        if (const std::vector<moist::exact::Kernel::Point_3>* polygon = std::get_if<std::vector<moist::exact::Kernel::Point_3>>(&*result))
        {
            polygon_points.reserve(polygon->size());
            std::transform(polygon->begin(), polygon->end(), polygon_points.begin(), [](const moist::exact::Kernel::Point_3 point)
            {
                return moist::exact::Point(static_cast<moist::exact::Kernel::Point_3>(point));
            });
        }
        else if (const moist::exact::Kernel::Triangle_3* triangle = std::get_if<moist::exact::Kernel::Triangle_3>(&*result))
        {
            polygon_points.reserve(3);
            polygon_points[0] = moist::exact::Point(static_cast<moist::exact::Kernel::Point_3>(triangle->vertex(0)));
            polygon_points[1] = moist::exact::Point(static_cast<moist::exact::Kernel::Point_3>(triangle->vertex(1)));
            polygon_points[2] = moist::exact::Point(static_cast<moist::exact::Kernel::Point_3>(triangle->vertex(2)));
        }
        else
        {
            // We don't care about no intersection, intersection in a single point, or intersection on one edge.
            return std::nullopt;
        }

        return polygon_points;
    }

    /**
     * @brief Since we are matching points from meshes with limited accuracy, in order to sort out intersections in equal lines, apply an epsilon to points here...
     *
     * @param p0
     * @param p1
     * @param epsilon
     * @return true
     * @return false
     */
    PURE INLINE bool points_are_equal(const moist::exact::Kernel::Point_3& p0, const moist::exact::Kernel::Point_3& p1, double epsilon = 1e-8)
    {
        double dx = std::abs(CGAL::to_double(p0.x() - p1.x()));
        double dy = std::abs(CGAL::to_double(p0.y() - p1.y()));
        return dx <= epsilon && dy <= epsilon;
    }

    PURE INLINE std::optional<moist::exact::Point> intersection(const moist::exact::EdgePoints& e0, const moist::exact::EdgePoints& e1)
    {
        const std::array<moist::exact::Point, 4> points
        {{
            e0.p0, e0.p1, e1.p0, e1.p1
        }};

        const moist::exact::Kernel::Segment_2 s0(
            moist::exact::Kernel::Point_2(e0.p0._p.x(), e0.p0._p.y()),
            moist::exact::Kernel::Point_2(e0.p1._p.x(), e0.p1._p.y())
        );

        const moist::exact::Kernel::Segment_2 s1(
            moist::exact::Kernel::Point_2(e1.p0._p.x(), e1.p0._p.y()),
            moist::exact::Kernel::Point_2(e1.p1._p.x(), e1.p1._p.y())
        );

        auto result = CGAL::intersection(s0, s1);
        if (result)
        {
            if (const moist::exact::Kernel::Point_2* p = std::get_if<moist::exact::Kernel::Point_2>(&*result))
            {
                if (is_strictly_inside(s0, *p) && is_strictly_inside(s1, *p))
                {
                    const double z = (e0.p0.z() + e0.p1.z()) / 2.0;
                    const auto point = moist::exact::Point(moist::exact::Kernel::Point_3(p->x(), p->y(), z));
                    for (const auto edge_point : points)
                    {
                        if (points_are_equal(edge_point._p, point._p))
                        {
                            return std::nullopt;
                        }
                    }
                    return point;
                }
            }
        }

        return std::nullopt;
    }

    PURE INLINE bool edge_exists(const moist::exact::EdgePoints& edge, const moist::ExactMesh& mesh)
    {
        for (const auto cell : mesh.Cells())
        {
            for (const auto& [i, j] : TET_EDGE_DESCRIPTOR)
            {
                const auto p0 = mesh.Point(i);
                const auto p1 = mesh.Point(j);
                bool has_ep0 = moist::geometry::exact::points_are_equal(p0._p, edge.p0._p) || moist::geometry::exact::points_are_equal(p1._p, edge.p0._p);
                bool has_ep1 = moist::geometry::exact::points_are_equal(p0._p, edge.p1._p) || moist::geometry::exact::points_are_equal(p1._p, edge.p1._p);

                if (has_ep0 && has_ep1)
                {
                    return true;
                }
            }
        }
        return false;
    }
}

#endif // MOIST_CORE_GEOMETRY_EXACT_INL_
