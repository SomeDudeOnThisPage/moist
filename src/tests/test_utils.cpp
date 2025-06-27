#include "test_utils.hpp"

#include <string>

#include <geogram/basic/process.h>

#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"
#include "moist/core/attributes.inl"

// all of this code is shit, and only used for testing...
bool moist::test::contains_overlapping_constraints(geogram::Mesh& a, geogram::Mesh& b, moist::Interface& interface)
{
    // check if a and b have any overlapping edges on interface
    const auto edges_a = moist::geometry::collect_interface_edges(a, *interface.Plane());
    const auto edges_b = moist::geometry::collect_interface_edges(b, *interface.Plane());

    for (const auto e0 : edges_a)
    {
        const vec3 e0p0 = a.vertices.point(e0.v0);
        const vec3 e0p1 = a.vertices.point(e0.v1);

        for (const auto e1 : edges_b)
        {
            const vec3 e1p0 = b.vertices.point(e1.v0);
            const vec3 e1p1 = b.vertices.point(e1.v1);

            if (e0p0 == e1p0 && e0p1 == e1p1 || e0p0 == e1p1 && e0p1 == e1p0)
            {
                continue;
            }

            if (!moist::predicates::xy::check_lines_aabb(
                reinterpret_cast<const vec2&>(e0p0),
                reinterpret_cast<const vec2&>(e0p1),
                reinterpret_cast<const vec2&>(e1p0),
                reinterpret_cast<const vec2&>(e1p1)
            ))
            {
                continue;
            }

            const auto intersection_opt = moist::predicates::xy::get_line_intersection(
                reinterpret_cast<const vec2&>(e0p0),
                reinterpret_cast<const vec2&>(e0p1),
                reinterpret_cast<const vec2&>(e1p0),
                reinterpret_cast<const vec2&>(e1p1)
            );

            if (!intersection_opt.has_value())
            {
                continue;
            }

            if (intersection_opt.value() == reinterpret_cast<const vec2&>(e0p0) || intersection_opt.value() == reinterpret_cast<const vec2&>(e0p1))
            {
                // we don't handle the "intersecting on a vertex case" here
                continue;
            }

            OOC_WARNING("intersecting edges [" << e0p0.x << ", " << e0p0.y << "] -> [" << e0p1.x << ", " << e0p1.y << "], ["
                << e1p0.x << ", " << e1p0.y << "] -> [" << e1p1.x << ", " << e1p1.y << "] " << " at [" << intersection_opt.value().x << ", " << intersection_opt.value().y << "]");
        }
    }
    return true;
}

// this code is shit
bool moist::test::contains_constraints(geogram::Mesh& mesh, geogram::Mesh& constraints, moist::Interface& interface, size_t steiner_points)
{
    constexpr std::string_view attribute = "__contains_constraints";
    constraints.facets.compute_borders();

    const auto constraint_edges = moist::geometry::collect_edges(constraints);
    const auto mesh_edges = moist::geometry::collect_interface_edges(mesh, *interface.Plane());
    std::vector<bool> found(constraint_edges.size(), false);

    size_t i = 0;
    for (const auto e0 : constraint_edges)
    {
        for (const auto e1 : mesh_edges)
        {
            const vec3 e0p0 = constraints.vertices.point(e0.v0);
            const vec3 e0p1 = constraints.vertices.point(e0.v1);
            const vec3 e1p0 = mesh.vertices.point(e1.v0);
            const vec3 e1p1 = mesh.vertices.point(e1.v1);

            if (e0p0 == e1p0 && e0p1 == e1p1 || e0p0 == e1p1 && e0p1 == e1p0)
            {
                found[i] = true;
                break;
            }
        }
        i++;
    }


    /*size_t i = 0;
    for (const auto e : constraint_edges)
    {
        //const g_index v0 = constraints.edges.vertex(e, 0);
        //const g_index v1 = constraints.edges.vertex(e, 1);
        const vec3 p0 = constraints.vertices.point(e.v0);
        const vec3 p1 = constraints.vertices.point(e.v1);

        // precheck if both vertices are part of the mesh
        bool has_v0, has_v1 = false;
        for (const g_index v : mesh.vertices)
        {
            if (mesh.vertices.point(v) == p0)
            {
                has_v0 = true;
            }

            if (mesh.vertices.point(v) == p1)
            {
                has_v1 = true;
            }
        }

        // ignore edges that are not in the mesh due to being outside...
        if (!has_v0 || !has_v1)
        {
            continue;
        }

        //geogram::parallel_for(0, mesh.cells.nb(), [&mesh, &constraints, &attribute, &p0, &p1, &e] (const g_index c)
        //{
        for (const g_index c : mesh.cells)
        {
            for (l_index le = 0; le < mesh.cells.nb_edges(c); le++)
            {
                const g_index lev0 = mesh.cells.edge_vertex(c, le, 0);
                const g_index lev1 = mesh.cells.edge_vertex(c, le, 1);
                const vec3 lep0 = mesh.vertices.point(lev0);
                const vec3 lep1 = mesh.vertices.point(lev1);

                if (p0 == lep0 && p1 == lep1 || p1 == lep0 && p0 == lep1)
                {
                    found[i] = true;
                    break;
                }
            }
        }

        //});

        i++;
    }*/

    size_t missing = 0;
    i = 0;
    for (const auto e : constraint_edges)
    {
        if (!found[i])
        {
            const vec3 p0 = constraints.vertices.point(e.v0);
            const vec3 p1 = constraints.vertices.point(e.v1);

            // attempt to find a cell which the middle between p0 and p1 belong to... if none exist, we can be reasonably sure that this is a false missing due to depressed cells from the interface
            // if a cell exists, we can be reasonably sure that there is some form of fuckery going on
            const vec3 middle = ((p0 + p1) / 2.0);
            for (const g_index c : mesh.cells)
            {
                const auto pit = moist::predicates::point_in_tet(mesh, c, middle);
                if (pit == moist::predicates::PointInTet::FACET || pit == moist::predicates::PointInTet::EDGE || pit == moist::predicates::PointInTet::INSIDE || pit == moist::predicates::PointInTet::VERTEX)
                {
                    OOC_WARNING("missing edge " << i << ": [" << p0.x << ", " << p0.y << "] -> [" << p1.x << ", " << p1.y << "]");
                    missing++;
                    break;
                }
            }

        }
        i++;
    }

    // TODO: Subtract edges split by steiner points...
    OOC_WARNING("missing " << missing << " edges with " << steiner_points << " steiner points, totaling " << (missing - steiner_points) << " out of " << constraints.edges.nb() << " missing edges...");

    return !(missing - steiner_points > 0);
}
