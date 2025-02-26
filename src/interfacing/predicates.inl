#ifndef __OOC_PREDICATES_CPP
#define __OOC_PREDICATES_CPP

#include <geogram/mesh/mesh.h>
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
}

#endif // __OOC_PREDICATES_CPP
