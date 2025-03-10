#include <vector>
#include <geogram/mesh/mesh_repair.h>

#include "interface.hpp"
#include "predicates.inl"

ooc::Interface::Interface(const InterfacePlane plane) : _constraints(geogram::Mesh(2)), _triangulation(geogram::Mesh(3)), _plane(plane)
{
    if (geogram::length(geogram::cross(_plane.a, _plane.b)) <= 1e-6)
    {
        OOC_ERROR("Vectors do not span a plane!");
    }
}

void ooc::Interface::AddConstraints(std::string name, geogram::Mesh& mesh)
{
    std::map<GEO::index_t, GEO::index_t> mesh_to_interface;

    // Add vertices of current mesh to union of interface vertices.
    for (auto vertex_id : mesh.vertices)
    {
        auto point = mesh.vertices.point(vertex_id);
        if (predicates::point_on_plane(point, this->_plane))
        {
            GEO::index_t interface_vertex_id = this->_constraints.vertices.create_vertex(point.data());
            this->_indices[{point[0], point[1]}] = interface_vertex_id;
            mesh_to_interface[vertex_id] = interface_vertex_id;

            this->_interface_vertices[name][interface_vertex_id] = vertex_id;
        }
    }

    // Collect edges not shared between more than 2 Facets...
    std::map<std::pair<GEO::index_t, GEO::index_t>, unsigned int> edge_count_map;
    for (auto facet : mesh.facets)
    {
        // We only care about facets on the interface plane! TODO: Make this inline with other loop, can we add vertices at the same time?
        if (!predicates::facet_on_plane(facet, mesh, this->_plane))
        {
            continue;
        }

        GEO::index_t nb_local_vertices = mesh.facets.nb_vertices(facet);
        for (GEO::index_t local_vertex = 0; local_vertex < nb_local_vertices; ++local_vertex)
        {
            GEO::index_t v1 = mesh.facets.vertex(facet, local_vertex);
            GEO::index_t v2 = mesh.facets.vertex(facet, (local_vertex + 1) % nb_local_vertices);

            if (v1 > v2)
            {
                std::swap(v1, v2);
            }

            edge_count_map[{v1, v2}]++;
        }
    }

    for (const auto& entry : edge_count_map)
    {
        // One adjacent interface facet -> edge must be a boundary edge.
        if (entry.second == 1)
        {
            // Translate Edge-Vertices from Mesh-Global index to Interface-Local index.
            auto global = entry.first;
            this->_constraints.edges.create_edge(mesh_to_interface[global.first], mesh_to_interface[global.second]);
        }
    }
}

void ooc::Interface::Triangulate()
{
    geogram::mesh_repair(this->_constraints, geogram::MESH_REPAIR_COLOCATE);

    // Only the constraint-mesh is required, must pass 0 and nullptr here to set_vertices when using triangle.
    // TODO: REALLY annoying thing about triangle... if it creates a bbox ("master-triangles") where a vertex of the bbox corresponds with a vertex of the
    //       input, a double vertex is found and "ignored", which leads to a segfault when geogram tries to read the data back into it's
    //       own data structure... somehow create own bbox for triangle?
    auto delaunay = geogram::Delaunay::create(2, "triangle");
    delaunay->set_constraints(&this->_constraints);
    delaunay->set_vertices(0, nullptr);
    OOC_DEBUG("Triangulated #" << delaunay->nb_vertices() << " interface vertices");

    geogram::vector<double> vertices(delaunay->nb_vertices() * 3);
    for(auto v = 0; v < delaunay->nb_vertices(); v++)
    {
        vertices[3 * v] = delaunay->vertex_ptr(v)[0];
        vertices[3 * v + 1] = delaunay->vertex_ptr(v)[1];
        vertices[3 * v + 2] = 0.0;
    }

    geogram::vector<geogram::index_t> triangles(delaunay->nb_cells() * 3);
    for(auto t = 0; t < delaunay->nb_cells(); t++)
    {
        triangles[3 * t] = geogram::index_t(delaunay->cell_vertex(t, 0));
        triangles[3 * t + 1] = geogram::index_t(delaunay->cell_vertex(t, 1));
        triangles[3 * t + 2] = geogram::index_t(delaunay->cell_vertex(t, 2));
    }

    this->_triangulation.facets.assign_triangle_mesh((GEO::coord_index_t) 3, vertices, triangles, true);
}

const bool ooc::Interface::HasMeshConstraints(std::string mesh)
{
    return _interface_vertices.find(mesh) != _interface_vertices.end();
}

const geogram::index_t ooc::Interface::GetMappedVertex(std::string mesh, geogram::index_t v_id)
{
    if (this->_interface_vertices[mesh].contains(v_id))
    {
        return this->_interface_vertices[mesh][v_id];
    }

    return false;
}

const geogram::Mesh* ooc::Interface::Triangulation()
{
    return &this->_triangulation;
}

#ifndef NDEBUG
void ooc::export_delaunay(const std::string filename, geogram::Mesh& mesh, int dimension)
{
    GEO::MeshIOFlags flags;
    flags.set_dimension(dimension);
    flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    if (!GEO::mesh_save(mesh, filename, flags))
    {
        OOC_ERROR("Error: Could not save mesh to " << filename);
    }
}
#endif // NDEBUG
