#include "interface_generator.hpp"

#include <unordered_set>

#include <geogram/delaunay/delaunay.h>
#include <geogram/mesh/mesh_repair.h>
#ifndef OPTION_PARALLEL_LOCAL_OPERATIONS
#include <geogram/basic/process.h>
#endif

#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/mesh_quality.inl"

moist::InterfaceGenerator::InterfaceGenerator(const AxisAlignedInterfacePlane plane) : _constraints(geogram::Mesh(2, false))
{
    this->_triangulation = std::make_shared<geogram::Mesh>(3);
    this->_plane = std::make_shared<AxisAlignedInterfacePlane>(plane);
    this->_unique_vertices = 0;
    this->WriteMetadata();
}

void moist::InterfaceGenerator::AddConstraints(const geogram::Mesh &mesh)
{
    std::map<g_index, g_index> mesh_to_interface;
    std::map<std::pair<g_index, g_index>, uint8_t> edge_count_map; // if we have more than 255 incident edges we have a waaaay bigger problem anyway...

#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    std::mutex mtx;
    geogram::parallel_for(0, mesh.vertices.nb(), [this, &mesh, &mtx, &mesh_to_interface](const g_index v) {
#else
    for (const g_index v : mesh.vertices) {
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS
        const vec3 point = mesh.vertices.point(v);
        if (predicates::point_on_plane(point, *this->_plane))
        {
        #ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
            std::lock_guard<std::mutex> lock(mtx);
        #endif
            g_index interface_vertex_id = this->_constraints.vertices.create_vertex(point.data());
            this->_unique_vertices++;
            mesh_to_interface[v] = interface_vertex_id;
        #ifndef NDEBUG
            this->_required_vertices.push_back(point);
        #endif // NDEBUG
        }
    }
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    );
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

//#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
//    geogram::parallel_for(0, mesh.facets.nb(), [this, &mesh, &edge_count_map, &mtx](const g_index f) {
//#else
    for (const g_index f : mesh.facets) {
//#endif // OPTION_PARALLEL_LOCAL_OPERATIONS
    //#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS // even when locking the entire parallel lambda, inside triangle occurs a segfault...
    //    std::lock_guard<std::mutex> lock(mtx);
    //#endif
        if (predicates::facet_on_plane(f, mesh, *this->_plane))
        {
            const geogram::index_t nb_local_vertices = mesh.facets.nb_vertices(f);
            for (l_index v = 0; v < nb_local_vertices; v++)
            {
                g_index v1 = mesh.facets.vertex(f, v);
                g_index v2 = mesh.facets.vertex(f, (v + 1) % nb_local_vertices);

                if (v1 > v2)
                {
                    std::swap(v1, v2);
                }

                {
                    edge_count_map[{v1, v2}]++;
                }
            }
        }
    }
//#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
//    );
//#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    for (const auto& entry : edge_count_map)
    {
        // one adjacent interface facet -> edge must be a boundary edge.
        if (entry.second == 1)
        {
            // translate edge-vertices from mesh-global index to interface-local index.
            auto global = entry.first;
            this->_constraints.edges.create_edge(mesh_to_interface[global.first], mesh_to_interface[global.second]);
        }
    }
}

void moist::InterfaceGenerator::Triangulate()
{
    // Only the constraint-mesh is required, must pass 0 and nullptr here to set_vertices when using triangle.
    // TODO: REALLY annoying thing about triangle... it seems to create a bbox ("master-triangles") where a vertex of the bbox corresponds with a vertex of the
    //       input, a double vertex is found and "ignored", which leads to a segfault when geogram tries to read the data back into it's
    //       own data structure... somehow create own bbox for triangle?
    // TODO: this doesn't actually help fix things...
    geogram::mesh_repair(this->_constraints, geogram::MESH_REPAIR_COLOCATE);

    auto delaunay = geogram::Delaunay::create(2, "triangle");
    delaunay->set_constraints(&this->_constraints);
    delaunay->set_vertices(0, nullptr);

    OOC_DEBUG("Triangulated #" << delaunay->nb_vertices() << " interface vertices");

    geogram::vector<double> vertices(delaunay->nb_vertices() * 3);
    for(g_index v = 0; v < delaunay->nb_vertices(); v++)
    {
        vertices[3 * v] = delaunay->vertex_ptr(v)[0];
        vertices[3 * v + 1] = delaunay->vertex_ptr(v)[1];
        vertices[3 * v + 2] = 0.0;
    }

    geogram::vector<geogram::index_t> triangles(delaunay->nb_cells() * 3);
    for(g_index t = 0; t < delaunay->nb_cells(); t++)
    {
        triangles[3 * t] = geogram::index_t(delaunay->cell_vertex(t, 0));
        triangles[3 * t + 1] = geogram::index_t(delaunay->cell_vertex(t, 1));
        triangles[3 * t + 2] = geogram::index_t(delaunay->cell_vertex(t, 2));
    }

    this->_triangulation->facets.assign_triangle_mesh((geogram::coord_index_t) 3, vertices, triangles, true);
    this->_triangulation->facets.compute_borders();

    // assign attribute constraint edges to edges in triangulation
    geogram::Attribute<int> e_constrained(this->_triangulation->edges.attributes(), ATTRIBUTE_CONSTRAINT_EDGE);
    for (const g_index constraint_edge : this->_constraints.edges)
    {
        // TODO: This also only works with z-growing...
        const auto cep0p = this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(constraint_edge, 0));
        const auto cep1p = this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(constraint_edge, 1));
        const vec2 cep0 = vec2(cep0p[0], cep0p[1]);
        const vec2 cep1 = vec2(cep1p[0], cep1p[1]);

        // find edge in interface mesh, and mark it
        for (const g_index edge : this->_triangulation->edges)
        {
            const vec2 ep0 = reinterpret_cast<const vec2&>(_triangulation->vertices.point(this->_triangulation->edges.vertex(edge, 0)));
            const vec2 ep1 = reinterpret_cast<const vec2&>(_triangulation->vertices.point(this->_triangulation->edges.vertex(edge, 1)));
            if (ep0 == cep0 && ep1 == cep1 || ep0 == cep1 && ep1 == cep0)
            {
                e_constrained[edge] = 1;
            }
        }
    }

    geogram::Attribute<int> m_target_vertices(this->_triangulation->cells.attributes(), ATTRIBUTE_INTERFACE_TARGET_VERTICES);
    m_target_vertices[0] = this->_unique_vertices / 2.0;

#ifndef NDEBUG
    OOC_DEBUG("checking for missing vertices...");
    uint32_t missing = 0;
    for (const g_index v : this->_triangulation->vertices)
    {
        const vec3 point = this->_triangulation->vertices.point(v);
        bool found = false;
        for (const vec3 other : this->_required_vertices) // too lazy to add operator==...
        {
            if (reinterpret_cast<const vec2&>(point) == reinterpret_cast<const vec2&>(other))
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_DEBUG("missing point " << point);
            missing++;
        }
    }
    OOC_DEBUG("finished checking vertices, missing " << missing << " points...");
#endif // NDEBUG
}

static vec3 middle(const vec3 v0, const vec3 v1)
{
    return vec3(
        (v0.x + v1.x) / 2.0,
        (v0.y + v1.y) / 2.0,
        (v0.z + v1.z) / 2.0
    );
}

void moist::InterfaceGenerator::Decimate()
{
    const size_t n = (this->_triangulation->facets.nb() - this->_constraints.edges.nb()) / 2;

    // enumerate the n worst triangles
    std::vector<g_index> triangles(this->_triangulation->facets.nb());
    for (const g_index facet : this->_triangulation->facets)
    {
        triangles.push_back(facet);
    }

    const auto triangulation = this->_triangulation;
    std::sort(triangles.begin(), triangles.end(), [triangulation](const g_index a, const g_index b)
    {
        return moist::mesh_quality::tri::aspect_ratio(a, *triangulation) > moist::mesh_quality::tri::aspect_ratio(b, *triangulation);
    });

    // decimate shortest edge of n triangles
    for (size_t i = 0; i < n; i++)
    {
        const g_index f = triangles[i];
        l_index shortest_edge = -1;
        double shortest_edge_length = std::numeric_limits<double>::max();
        for (l_index le = 0; le < 3; le++)
        {
            const vec3 p0 = triangulation->vertices.point(triangulation->facets.vertex(f, le));
            const vec3 p1 = triangulation->vertices.point(triangulation->facets.vertex(f, (le + 1) % 3));
            const double distance = geogram::distance(p0, p1);
            if (distance < shortest_edge_length)
            {
                shortest_edge_length = distance;
                shortest_edge = le;
            }
        }

        const vec3 p0 = triangulation->vertices.point(triangulation->facets.vertex(f, shortest_edge));
        const vec3 p1 = triangulation->vertices.point(triangulation->facets.vertex(f, (shortest_edge + 1) % 3));
        const vec3 mid = middle(p0, p1);
        triangulation->vertices.point(triangulation->facets.vertex(f, shortest_edge)).x = mid.x;
        triangulation->vertices.point(triangulation->facets.vertex(f, shortest_edge)).y = mid.y;
        triangulation->vertices.point(triangulation->facets.vertex(f, shortest_edge)).z = mid.z;
        triangulation->vertices.point(triangulation->facets.vertex(f, (shortest_edge + 1) % 3)).x = mid.x;
        triangulation->vertices.point(triangulation->facets.vertex(f, (shortest_edge + 1) % 3)).y = mid.y;
        triangulation->vertices.point(triangulation->facets.vertex(f, (shortest_edge + 1) % 3)).z = mid.z;
    }

    // delete all resulting 0-volume triangles
    geogram::vector<g_index> deleted_facets(triangulation->facets.nb());
    for (g_index i = 0; i < triangulation->facets.nb(); i++)
    {
        const vec3 p0 = triangulation->vertices.point(triangulation->facets.vertex(i, 0));
        const vec3 p1 = triangulation->vertices.point(triangulation->facets.vertex(i, 1));
        const vec3 p2 = triangulation->vertices.point(triangulation->facets.vertex(i, 2));

        if (p0 == p1 || p0 == p2 || p1 == p2)
        {
            deleted_facets[i] = 1;
        }
    }

    triangulation->facets.delete_elements(deleted_facets, true);
}
