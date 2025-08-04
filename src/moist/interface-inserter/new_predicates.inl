#ifndef MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_
#define MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>
#include <CGAL/enum.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <cmath>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/predicates.inl"
#include "exact_mesh.hpp"

namespace moist::new_predicates
{
    typedef CGAL::Simple_cartesian<double> K_Approx;
    typedef K_Approx::Point_3 ApproxPoint;
    typedef K_Approx::Segment_3 ApproxSegment;

    static inline bool approx_collinear_helper(const ApproxPoint& a, const ApproxPoint& b, const ApproxPoint& c, double epsilon)
    {
        double area = std::abs((b.x() - a.x()) * (c.y() - a.y()) -
                            (c.x() - a.x()) * (b.y() - a.y()));
        return area <= epsilon;
    }
    static inline bool is_point_on_tetrahedron_edge(const ApproxPoint& p, const std::array<ApproxPoint, 4>& tetra, double epsilon)
    {
        std::array<std::pair<int, int>, 6> edges =
        {{
            {0, 1}, {0, 2}, {0, 3},
            {1, 2}, {1, 3}, {2, 3}
        }};

        for (const auto& [i, j] : edges)
        {
            ApproxSegment edge(tetra[i], tetra[j]);

            if (/*CGAL::collinear(tetra[i], tetra[j], p)*/ approx_collinear_helper(tetra[i], tetra[j], p, epsilon))
            {
                double d = std::abs(std::sqrt(CGAL::to_double(CGAL::squared_distance(p, edge))));
                if (d <= epsilon)
                {
                    return true;
                }
            }
        }
        return false;
    }

    static inline ApproxPoint convert_point(const moist::exact::Point& point)
    {
        return ApproxPoint(
            CGAL::to_double(point._p.x()),
            CGAL::to_double(point._p.y()),
            CGAL::to_double(point._p.z())
        );
    }

    PURE INLINE bool approx_collinear(const moist::exact::Point& p0, const moist::exact::Point& p1, const moist::exact::Point& p2, const double eps)
    {
        static CGAL::Cartesian_converter<moist::exact::Kernel, K_Approx> converter;
        return approx_collinear_helper(converter(p0._p), converter(p1._p), converter(p2._p), eps);
    }

    PURE INLINE moist::predicates::PointInTet point_in_tet_exact(const moist::ExactMesh& mesh, const std::size_t c, const moist::exact::Point& p, bool exclude_existing_points = false)
    {
        const auto p0 = mesh.Point(mesh.Cell(c)._points[0]);
        const auto p1 = mesh.Point(mesh.Cell(c)._points[1]);
        const auto p2 = mesh.Point(mesh.Cell(c)._points[2]);
        const auto p3 = mesh.Point(mesh.Cell(c)._points[3]);

        const auto ap = vec3(p.x(), p.y(), p.z());
        const auto ap0 = vec3(p0.x(), p0.y(), p0.z());
        const auto ap1 = vec3(p1.x(), p1.y(), p1.z());
        const auto ap2 = vec3(p2.x(), p2.y(), p2.z());
        const auto ap3 = vec3(p3.x(), p3.y(), p3.z());

        const std::array<ApproxPoint, 4> approx_tet = {{
            convert_point(mesh.Point(mesh.Cell(c)._points[0])),
            convert_point(mesh.Point(mesh.Cell(c)._points[1])),
            convert_point(mesh.Point(mesh.Cell(c)._points[2])),
            convert_point(mesh.Point(mesh.Cell(c)._points[3]))
        }};

        CGAL::Sign s[4];
        s[0] = CGAL::orientation(p._p, p1._p, p2._p, p3._p);
        s[1] = CGAL::orientation(p0._p, p._p, p2._p, p3._p);
        s[2] = CGAL::orientation(p0._p, p1._p, p._p, p3._p);
        s[3] = CGAL::orientation(p0._p, p1._p, p2._p, p._p);

        const bool inside_or_on_boundary = (
            (s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0) ||
            (s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0)
        );

        int nz = 0;
        for(int i = 0; i < 4; i++)
        {
            if(s[i] == CGAL::Sign::ZERO)
            {
                nz++;
            }
        }

        if (inside_or_on_boundary)
        {
            if (nz == 0)
            {
                return moist::predicates::PointInTet::INSIDE;
            }
            if (nz == 1)
            {
                return moist::predicates::PointInTet::FACET;
            }
            else if (nz == 2)
            {
                return moist::predicates::PointInTet::EDGE;
            }
            else if (nz == 3)
            {
                return moist::predicates::PointInTet::VERTEX;
            }
        }

        if (is_point_on_tetrahedron_edge(convert_point(p), approx_tet, 1e-14))
        {
            OOC_DEBUG("Approximating point " << convert_point(p) << " on edge");
            return moist::predicates::PointInTet::EDGE_APPROXIMATED;
        }

        return moist::predicates::PointInTet::NONE;
    }
}

#endif // MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_
