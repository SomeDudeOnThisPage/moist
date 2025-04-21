#include "local_operations.hpp"

#include "geometry.inl"
#include "attributes.inl"

void incremental_meshing::operation::edge_split_1to2(MeshSlice &mesh, const g_index cell, const CrossedEdge &edge, const AxisAlignedInterfacePlane &plane)
{
    const g_index v_opposite = incremental_meshing::geometry::non_coplanar_opposite(cell, edge.e_v0, edge.e_v1, mesh, plane);
    const g_index v_coplanar_opposite = incremental_meshing::geometry::other(cell, edge.e_v0, edge.e_v1, v_opposite, mesh);

    {
        std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR);
        geogram::Attribute<incremental_meshing::attributes::VertexDescriptorFlags> v_descriptor(mesh.vertices.attributes(), ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS);
        v_descriptor[edge.p] |= incremental_meshing::attributes::VertexDescriptorFlags::DISCARD;
    }

    mesh.CreateTetrahedra({
        {edge.e_v0, v_coplanar_opposite, v_opposite, edge.p},
        {edge.e_v1, v_coplanar_opposite, v_opposite, edge.p}
    });
    mesh.DeleteTetrahedra(cell);
}

void incremental_meshing::operation::edge_split_1to3(MeshSlice &mesh, const g_index cell, const CrossedEdge &edge0, const CrossedEdge &edge1, const AxisAlignedInterfacePlane &plane)
{
    const g_index shared = (edge0.e_v0 == edge1.e_v0 || edge0.e_v0 == edge1.e_v1) ? edge0.e_v0 : edge0.e_v1;
    const g_index v_opposite = incremental_meshing::geometry::other(
        cell,
        edge0.e_v0,
        edge0.e_v1,
        (edge1.e_v0 == shared) ? edge1.e_v1 : edge1.e_v0,
        mesh
    );

    const g_index v_coplanar_opposite_p0 = (edge1.e_v0 != shared) ? edge1.e_v0 : edge1.e_v1;
    const g_index v_coplanar_opposite_p1 = (edge0.e_v0 != shared) ? edge0.e_v0 : edge0.e_v1;

    if (incremental_meshing::geometry::point_of_cell(mesh, cell, mesh.vertices.point(edge0.p)) ||  incremental_meshing::geometry::point_of_cell(mesh, cell, mesh.vertices.point(edge1.p)))
    {
        OOC_DEBUG("prevented 1 -> 3 split due to exising point");
        return;
    }

    {
        std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR);
        geogram::Attribute<incremental_meshing::attributes::VertexDescriptorFlags> v_descriptor(mesh.vertices.attributes(), ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS);
        v_descriptor[edge0.p] |= incremental_meshing::attributes::VertexDescriptorFlags::DISCARD;
        v_descriptor[edge1.p] |= incremental_meshing::attributes::VertexDescriptorFlags::DISCARD;
    }

    mesh.CreateTetrahedra({
        {v_opposite, shared, edge0.p, edge1.p},
        {v_opposite, edge0.p, edge1.p, v_coplanar_opposite_p0},
        {v_opposite, edge0.p, v_coplanar_opposite_p0, v_coplanar_opposite_p1}
    });
    mesh.DeleteTetrahedra(cell);
}

void incremental_meshing::operation::vertex_insert_1to2()
{
    // honestly this is so niche I don't have the time right now...
}

// why did I name this 2to4??
void incremental_meshing::operation::vertex_insert_2to4(MeshSlice &mesh, const g_index cell, const vec3 &point, const incremental_meshing::AxisAlignedInterfacePlane &plane)
{
    // if we already have the vertex, only mark it as "interface"
    for (l_index lv = 0; lv < 4; lv++)
    {
        const g_index v = mesh.cells.vertex(cell, lv);
        const vec3 p = mesh.vertices.point(v);
        if (incremental_meshing::predicates::vec_eq_2d(point, p, plane))
        {
            std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR);
            geogram::Attribute<incremental_meshing::attributes::VertexDescriptorFlags> v_descriptor(mesh.vertices.attributes(), ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS);
            v_descriptor[v] |= incremental_meshing::attributes::VertexDescriptorFlags::INTERFACE;
            return;
        }
    }

    const g_index v_opposite = incremental_meshing::geometry::non_interface_vertex(cell, mesh, plane);
    const auto [v0, v1, v2] = incremental_meshing::geometry::other(cell, v_opposite, mesh);

    const g_index p = mesh.vertices.create_vertex(point.data());
    {
        geogram::Attribute<incremental_meshing::attributes::VertexDescriptorFlags> v_descriptor(mesh.vertices.attributes(), ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS);
        v_descriptor[p] |= incremental_meshing::attributes::VertexDescriptorFlags::INTERFACE;
    }

    mesh.CreateTetrahedra({
        {v_opposite, v0, v1, p},
        {v_opposite, v1, v2, p},
        {v_opposite, v2, v0, p}
    });
    mesh.DeleteTetrahedra(cell);
}
