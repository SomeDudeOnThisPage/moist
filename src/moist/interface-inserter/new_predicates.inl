#ifndef MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_
#define MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>
#include <CGAL/enum.h>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/predicates.inl"
#include "exact_mesh.hpp"

namespace moist::new_predicates
{
    PURE INLINE moist::predicates::PointInTet point_in_tet_exact(const moist::ExactMesh& mesh, const std::size_t c, const moist::ExactMesh::ExactPoint& p, bool exclude_existing_points = false)
    {
        const auto p0 = mesh.Point(mesh.Cell(c)._points[0]);
        const auto p1 = mesh.Point(mesh.Cell(c)._points[1]);
        const auto p2 = mesh.Point(mesh.Cell(c)._points[2]);
        const auto p3 = mesh.Point(mesh.Cell(c)._points[3]);

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

        return moist::predicates::PointInTet::NONE;
    }
}

#endif // MOIST_INTERFACE_INSERTER_NEW_PREDICATES_INL_
