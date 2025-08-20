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
#include <CGAL/squared_distance_3.h>
#include <cmath>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/predicates.inl"
#include "exact_mesh.hpp"

namespace moist::new_predicates
{
    typedef CGAL::Simple_cartesian<double> K_Approx;
    typedef K_Approx::Point_3 ApproxPoint;
    typedef K_Approx::Vector_3 ApproxVector;
    typedef K_Approx::Segment_3 ApproxSegment;

    static inline bool approx_collinear_helper(const ApproxPoint& P, const ApproxPoint& A, const ApproxPoint& B, double epsilon)
    {
        ApproxVector AB = B - A;
        ApproxVector AP = P - A;

        // Check collinearity via cross product
        ApproxVector cross = CGAL::cross_product(AB, AP);
        double cross_norm_sq = cross.squared_length();
        double AB_norm_sq = AB.squared_length();
        if (cross_norm_sq > epsilon * epsilon * AB_norm_sq)
        {
            return false;
        }

        double dot = AP * AB;
        double t = dot / AB.squared_length();

        return t >= -epsilon && t <= 1.0 + epsilon;
    }

    static inline ApproxPoint convert_point(const moist::exact::Point& point)
    {
        return ApproxPoint(
            CGAL::to_double(point._p.x()),
            CGAL::to_double(point._p.y()),
            CGAL::to_double(point._p.z())
        );
    }

    PURE INLINE double approx_volume(const moist::exact::Point& p0, const moist::exact::Point& p1, const moist::exact::Point& p2, const moist::exact::Point& p3)
    {
        const vec3 ap1 = vec3(p0.x(), p0.y(), p0.z());
        const vec3 ap2 = vec3(p1.x(), p1.y(), p1.z());
        const vec3 ap3 = vec3(p2.x(), p2.y(), p2.z());
        const vec3 ap4 = vec3(p3.x(), p3.y(), p3.z());

        return geo::dot(ap2 - ap1, geo::cross(ap3 - ap1, ap4 - ap1)) / 6.0;
    }

    PURE INLINE bool approx_collinear(const moist::exact::Point& p0, const moist::exact::Point& p1, const moist::exact::Point& p2, const double eps)
    {
        static CGAL::Cartesian_converter<moist::exact::Kernel, K_Approx> converter;
        return approx_collinear_helper(converter(p0._p), converter(p1._p), converter(p2._p), eps);
    }

    static inline bool is_equal_test(const ApproxPoint& a, const ApproxPoint& b, double epsilon = 1e-12)
    {
        return std::fabs(a.x() - b.x()) <= epsilon && std::fabs(a.y() - b.y()) <= epsilon && std::fabs(a.z() - b.z()) <= epsilon;
    }

    static inline bool is_point_on_tetrahedron_edge(const moist::exact::Point& point, const std::size_t& c, const moist::ExactMesh& mesh, double epsilon)
    {
        const std::array<ApproxPoint, 4> approx_tet = {{
            convert_point(mesh.Point(mesh.Cell(c)._points[0])),
            convert_point(mesh.Point(mesh.Cell(c)._points[1])),
            convert_point(mesh.Point(mesh.Cell(c)._points[2])),
            convert_point(mesh.Point(mesh.Cell(c)._points[3]))
        }};

        std::array<std::pair<int, int>, 6> edges =
        {{
            {0, 1}, {0, 2}, {0, 3},
            {1, 2}, {1, 3}, {2, 3}
        }};

        for (const auto& [i, j] : edges)
        {
            if (!mesh.Point(mesh.Cell(c)._points[i])._interface || !mesh.Point(mesh.Cell(c)._points[j])._interface)
            {
                continue;
            }

            ApproxPoint ap = convert_point(point);
            ApproxSegment edge(approx_tet[i], approx_tet[j]);

            if (approx_collinear_helper(ap, approx_tet[i], approx_tet[j], epsilon))
            {
                //double d = std::abs(std::sqrt(CGAL::to_double(CGAL::squared_distance(ap, edge))));
                //if (d <= epsilon)
                //{
                    if (!is_equal_test(ap, approx_tet[0]) && !is_equal_test(ap, approx_tet[1]) && !is_equal_test(ap, approx_tet[2]) && !is_equal_test(ap, approx_tet[3]))
                    {
                        return true;
                    }
                //}
            }
        }
        return false;
    }

    PURE INLINE CGAL::Sign orient3d(const ApproxPoint& a, const ApproxPoint& b, const ApproxPoint& c, const ApproxPoint& d)
    {
        double abx = b.x() - a.x();
        double aby = b.y() - a.y();
        double abz = b.z() - a.z();

        double acx = c.x() - a.x();
        double acy = c.y() - a.y();
        double acz = c.z() - a.z();

        double adx = d.x() - a.x();
        double ady = d.y() - a.y();
        double adz = d.z() - a.z();

        // scalar triple product?
        const double volume = abx * (acy * adz - acz * ady)
                    - aby * (acx * adz - acz * adx)
                    + abz * (acx * ady - acy * adx);

        if (std::abs(volume) < 1e-12)
        {
            return CGAL::Sign::ZERO;
        }

        return (volume > 0)
            ? CGAL::Sign::POSITIVE
            : CGAL::Sign::NEGATIVE;
    }

    PURE INLINE moist::predicates::PointInTet point_in_tet_exact(const moist::ExactMesh& mesh, const std::size_t c, const moist::exact::Point& exact_p, bool exclude_existing_points = false)
    {
        // TODO: calc orient3d with exact triple scalar, and apply epsilon only onto that!
        const auto exact_p0 = mesh.Point(mesh.Cell(c)._points[0]);
        const auto exact_p1 = mesh.Point(mesh.Cell(c)._points[1]);
        const auto exact_p2 = mesh.Point(mesh.Cell(c)._points[2]);
        const auto exact_p3 = mesh.Point(mesh.Cell(c)._points[3]);

        const auto p = convert_point(exact_p._p);
        const auto p0 = convert_point(exact_p0._p);
        const auto p1 = convert_point(exact_p1._p);
        const auto p2 = convert_point(exact_p2._p);
        const auto p3 = convert_point(exact_p3._p);

        CGAL::Sign s[4];
        s[0] = CGAL::orientation(exact_p._p, exact_p1._p, exact_p2._p, exact_p3._p);
        s[1] = CGAL::orientation(exact_p0._p, exact_p._p, exact_p2._p, exact_p3._p);
        s[2] = CGAL::orientation(exact_p0._p, exact_p1._p, exact_p._p, exact_p3._p);
        s[3] = CGAL::orientation(exact_p0._p, exact_p1._p, exact_p2._p, exact_p._p);
        //s[0] = orient3d(p, p1, p2, p3);
        //s[1] = orient3d(p0, p, p2, p3);
        //s[2] = orient3d(p0, p1, p, p3);
        //s[3] = orient3d(p0, p1, p2, p);

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
                // shared edges of the faces evaluated by orient3d above
                // e.g. the faces which determine s[0] and s[1] have p2 and p3 in common
                if (s[0] == CGAL::ZERO && s[1] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE23;
                }
                if (s[0] == CGAL::ZERO && s[2] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE13;
                }
                if (s[0] == CGAL::ZERO && s[3] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE12;
                }
                if (s[1] == CGAL::ZERO && s[2] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE03;
                }
                if (s[1] == CGAL::ZERO && s[3] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE02;
                }
                if (s[2] == CGAL::ZERO && s[3] == CGAL::ZERO)
                {
                    return moist::predicates::PointInTet::EDGE01;
                }
                OOC_WARNING("could not identify edge!");
                return moist::predicates::PointInTet::EDGE;
            }
            else if (nz == 3)
            {
                return moist::predicates::PointInTet::VERTEX;
            }
        }

        return moist::predicates::PointInTet::NONE;
    }
}

#endif // MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_
