#include "local_operations.hpp"

#include "moist/core/geometry.inl"
#include "moist/core/attributes.inl"
#include "moist/core/utils.hpp"

#ifndef NDEBUG
#define VECE(a, b, c, d) ((a) == (b) || (a) == (c) || (a) == (d) || \
                          (b) == (c) || (b) == (d) || \
                          (c) == (d))
#endif // NDEBUG

void moist::operation::edge_split_1to2(MeshSlice &mesh, const g_index cell, const CrossedEdge &edge, const AxisAlignedInterfacePlane &plane)
{
    // bandaid-fix, must fix predicates...
    const auto p0 = mesh.vertices.point(edge.p);
    const auto p1 = mesh.vertices.point(edge.e_v0);
    const auto p2 = mesh.vertices.point(edge.e_v1);
    if (p0 == p1 || p0 == p2)
    {
        return;
    }

    const g_index v_opposite = moist::geometry::non_coplanar_opposite(cell, edge.e_v0, edge.e_v1, mesh, plane);
    const g_index v_coplanar_opposite = moist::geometry::other(cell, edge.e_v0, edge.e_v1, v_opposite, mesh);

    {
        // points must be marked as discardable, as they must ultimately be collapsed onto other vertices.
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_discard(mesh.vertices.attributes(), ATTRIBUTE_DISCARD);

        v_discard[edge.e_v0] = false;
        v_discard[edge.e_v1] = false;
        v_discard[edge.p] = true;
    }

#ifndef NDEBUG
    // check to make sure bandaid-fix worked
    const auto p0t = mesh.vertices.point(edge.e_v0);
    const auto p1t = mesh.vertices.point(v_coplanar_opposite);
    const auto p2t = mesh.vertices.point(v_opposite);
    const auto p3t = mesh.vertices.point(edge.p);
    const auto p4t = mesh.vertices.point(edge.e_v1);

    const auto c0 = mesh.vertices.point(mesh.cells.vertex(cell, 0));
    const auto c1 = mesh.vertices.point(mesh.cells.vertex(cell, 1));
    const auto c2 = mesh.vertices.point(mesh.cells.vertex(cell, 2));
    const auto c3 = mesh.vertices.point(mesh.cells.vertex(cell, 3));
    if (VECE(p0t, p1t, p2t, p3t) || VECE(p4t, p1t, p2t, p3t))
    {
        OOC_WARNING("zero volume tet in 1->2 split");
    }
#endif // NDEBUG

    mesh.CreateTetrahedra({
        {edge.e_v0, v_coplanar_opposite, v_opposite, edge.p},
        {edge.e_v1, v_coplanar_opposite, v_opposite, edge.p}
    });
    mesh.DeleteTetrahedra(cell);
}

void moist::operation::edge_split_1to3(MeshSlice &mesh, const g_index cell, const CrossedEdge &edge0, const CrossedEdge &edge1, const AxisAlignedInterfacePlane &plane)
{
    const g_index shared = (edge0.e_v0 == edge1.e_v0 || edge0.e_v0 == edge1.e_v1) ? edge0.e_v0 : edge0.e_v1;
    const g_index v_opposite = moist::geometry::other(
        cell,
        edge0.e_v0,
        edge0.e_v1,
        (edge1.e_v0 == shared) ? edge1.e_v1 : edge1.e_v0,
        mesh
    );

    const g_index v_coplanar_opposite_p0 = (edge1.e_v0 != shared) ? edge1.e_v0 : edge1.e_v1;
    const g_index v_coplanar_opposite_p1 = (edge0.e_v0 != shared) ? edge0.e_v0 : edge0.e_v1;

    if (moist::geometry::point_of_cell(mesh, cell, mesh.vertices.point(edge0.p)) ||  moist::geometry::point_of_cell(mesh, cell, mesh.vertices.point(edge1.p)))
    {
        OOC_DEBUG("prevented 1 -> 3 split due to exising point");
        return;
    }

    {
        // points must be marked as discardable, as they must ultimately be collapsed onto other vertices.
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_discard(mesh.vertices.attributes(), ATTRIBUTE_DISCARD);

        v_discard[edge0.p] = true;
        v_discard[edge1.p] = true;
        v_discard[v_opposite] = false;
        v_discard[shared] = false;
        v_discard[v_coplanar_opposite_p0] = false;
        v_discard[v_coplanar_opposite_p1] = false;
    }

#ifndef NDEBUG
{
    // check to make sure bandaid-fix worked
    const auto p0t = mesh.vertices.point(v_opposite);
    const auto p1t = mesh.vertices.point(shared);
    const auto p2t = mesh.vertices.point(edge0.p);
    const auto p3t = mesh.vertices.point(edge1.p);
    if (p2t.x == -0.5954272150993347 || p3t.x == -0.5954272150993347)
    {
        int i  = 222;
    }


    const auto e0v0 = mesh.vertices.point(edge0.e_v0);
    const auto e0v1 = mesh.vertices.point(edge0.e_v1);
    const auto e1v0 = mesh.vertices.point(edge1.e_v0);
    const auto e1v1 = mesh.vertices.point(edge1.e_v1);
    if (VECE(p0t, p1t, p2t, p3t))
    {
        OOC_WARNING("zero volume tet in 1->3 split");
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_dbg(mesh.vertices.attributes(), "DebugAttribute");
        v_dbg[edge0.p] = 10*10;
        v_dbg[shared] = 10*10*10;
        v_dbg[v_coplanar_opposite_p0] = 10*10*10*10;
        v_dbg[v_coplanar_opposite_p1] = 10*10*10*10;
        moist::utils::geo::save("debug_attribute_mesh_zero_volume.geogram", mesh);
    }
}

{
    // check to make sure bandaid-fix worked
    const auto p0t = mesh.vertices.point(v_opposite);
    const auto p1t = mesh.vertices.point(v_coplanar_opposite_p0);
    const auto p2t = mesh.vertices.point(edge0.p);
    const auto p3t = mesh.vertices.point(edge1.p);
    if (VECE(p0t, p1t, p2t, p3t))
    {
        OOC_WARNING("zero volume tet in 1->3 split");
    }
}

{
    // check to make sure bandaid-fix worked
    const auto p0t = mesh.vertices.point(v_opposite);
    const auto p1t = mesh.vertices.point(v_coplanar_opposite_p0);
    const auto p2t = mesh.vertices.point(edge0.p);
    const auto p3t = mesh.vertices.point(v_coplanar_opposite_p1);
    if (VECE(p0t, p1t, p2t, p3t))
    {
        OOC_WARNING("zero volume tet in 1->3 split");
    }
}

#endif // NDEBUG

    mesh.CreateTetrahedra({
        {v_opposite, shared, edge0.p, edge1.p},
        {v_opposite, edge0.p, edge1.p, v_coplanar_opposite_p0},
        {v_opposite, edge0.p, v_coplanar_opposite_p0, v_coplanar_opposite_p1}
    });
    mesh.DeleteTetrahedra(cell);
}

static bool are_colinear(const vec3& a, const vec3& b, const vec3& c)
{
    return geogram::length(geogram::cross(b - a, c - a)) < moist::__DOUBLE_EPSILON;
}

void moist::operation::vertex_insert_1to2(MeshSlice &mesh, const g_index cell, const vec3& point, const moist::AxisAlignedInterfacePlane &plane)
{
    // honestly this is so niche I don't have the time right now...
    if (moist::predicates::point_of_tet(mesh, cell, point))
    {
        // TODO: recreate vertex, as to "reorder" it to the back of the vertices array.
        return;
    }

    const g_index p = mesh.vertices.create_vertex(point.data());
    const g_index v_opposite = moist::geometry::non_interface_vertex(cell, mesh, plane);
    const auto [v0, v1, v2] = moist::geometry::other(cell, v_opposite, mesh);
    const vec3 p0 = mesh.vertices.point(v0);
    const vec3 p1 = mesh.vertices.point(v1);
    const vec3 p2 = mesh.vertices.point(v2);
    const vec3 p3 = mesh.vertices.point(v_opposite);
    if (are_colinear(p0, p1, point))
    {
        mesh.CreateTetrahedra({
            {v_opposite, v0, v2, p},
            {v_opposite, v1, v2, p}
        });
    }
    else if (are_colinear(p1, p2, point))
    {
        mesh.CreateTetrahedra({
            {v_opposite, v0, v1, p},
            {v_opposite, v0, v2, p}
        });
    }
    else if (are_colinear(p0, p2, point))
    {
        mesh.CreateTetrahedra({
            {v_opposite, v1, v0, p},
            {v_opposite, v1, v2, p}
        });
    }

    mesh.DeleteTetrahedra(cell);
}

void moist::operation::vertex_insert_1to3(MeshSlice &mesh, const g_index cell, const vec3 &point, const moist::AxisAlignedInterfacePlane &plane)
{
    if (moist::predicates::point_of_tet(mesh, cell, point))
    {
        // TODO: recreate vertex, as to "reorder" it to the back of the vertices array.
        return;
    }

    const g_index p = mesh.vertices.create_vertex(point.data());
    const g_index v_opposite = moist::geometry::non_interface_vertex(cell, mesh, plane);
    const auto [v0, v1, v2] = moist::geometry::other(cell, v_opposite, mesh);

#ifndef NDEBUG
    // check to make sure bandaid-fix worked
    const auto p0t = mesh.vertices.point(v0);
    const auto p1t = mesh.vertices.point(v1);
    const auto p2t = mesh.vertices.point(v2);
    const auto p3t = mesh.vertices.point(p);
    const auto p4t = mesh.vertices.point(v_opposite);
    if (VECE(v_opposite, v0, v1, p) || VECE(v_opposite, v1, v2, p) || VECE(v_opposite, v2, v0, p))
    {
        OOC_ERROR("zero volume tet");
    }
#endif // NDEBUG

    mesh.CreateTetrahedra({
        {v_opposite, v0, v1, p},
        {v_opposite, v1, v2, p},
        {v_opposite, v2, v0, p}
    });
    mesh.DeleteTetrahedra(cell);
}
