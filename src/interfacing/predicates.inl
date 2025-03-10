#ifndef __OOC_PREDICATES_CPP
#define __OOC_PREDICATES_CPP

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/vecg.h>

#include "../core.hpp"
#include "interface.hpp"

namespace ooc::predicates
{
    // TODO: Needs correct projection to/from x/y-plane. For now this only works if the interface plane is parallel to the x/y axis.
    // TODO: Needs to account for +- epsilon along the normal of the plane.
    inline bool point_on_plane(const geogram::vec3& point, const ooc::InterfacePlane& plane)
    {
        return geogram::dot(point - plane.v, geogram::cross(plane.a, plane.b)) == 0.0;
    }

    inline bool facet_on_plane(const geogram::index_t facet, const geogram::Mesh& mesh, const ooc::InterfacePlane& plane)
    {
        for (geogram::index_t local_vertex = 0; local_vertex < mesh.facets.nb_vertices(facet); local_vertex++)
        {
            auto point = mesh.vertices.point(mesh.facets.vertex(facet, local_vertex));
            if (!point_on_plane(point, plane))
            {
                return false;
            }
        }

        return true;
    }

    // Source: geogram/mesh/mesh_AABB.cpp#175
    inline bool point_in_tet(const geogram::Mesh& M, geogram::index_t t, const geogram::vec3& p)
    {
        const geogram::vec3& p0 = geogram::Geom::mesh_vertex(M, M.cells.vertex(t, 0));
        const geogram::vec3& p1 = geogram::Geom::mesh_vertex(M, M.cells.vertex(t, 1));
        const geogram::vec3& p2 = geogram::Geom::mesh_vertex(M, M.cells.vertex(t, 2));
        const geogram::vec3& p3 = geogram::Geom::mesh_vertex(M, M.cells.vertex(t, 3));

        geogram::Sign s[4];
        s[0] = geogram::PCK::orient_3d(p, p1, p2, p3);
        s[1] = geogram::PCK::orient_3d(p0, p, p2, p3);
        s[2] = geogram::PCK::orient_3d(p0, p1, p, p3);
        s[3] = geogram::PCK::orient_3d(p0, p1, p2, p);

        return (
            (s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0) ||
            (s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0)
        );
    }
}

#endif // __OOC_PREDICATES_CPP
