#include "interface_generator.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include <geogram/delaunay/delaunay.h>
#include <geogram/mesh/mesh_repair.h>
#ifndef OPTION_PARALLEL_LOCAL_OPERATIONS
#include <geogram/basic/process.h>
#endif

#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/mesh_quality.inl"
#include "moist/core/utils.hpp"

moist::InterfaceGenerator::InterfaceGenerator(const AxisAlignedInterfacePlane plane) : _constraints(geogram::Mesh(2, false))
{
    this->_triangulation = std::make_shared<geogram::Mesh>(3);
    this->_plane = std::make_shared<AxisAlignedInterfacePlane>(plane);
    this->_unique_vertices = 0;
    this->WriteMetadata();
}

static double round16(double value)
{
    const double threshold = 1e-15;
    value = std::round(value * 1e16) / 1e16;
    return (std::abs(value) < threshold) ? 0.0 : value;
}

void moist::InterfaceGenerator::AddConstraints(const geogram::Mesh &mesh)
{
    std::map<g_index, g_index> mesh_to_interface;
    std::map<std::pair<g_index, g_index>, uint8_t> edge_count_map; // if we have more than 255 incident edges we have a waaaay bigger problem anyway...

    for (const g_index v : mesh.vertices)
    {
        // geogram uses precision values up to 16 decimal places...
        vec3 point = mesh.vertices.point(v);
        //point.x = round16(point.x);
        //point.y = round16(point.y);
        //point.z = round16(point.z);

        if (predicates::point_on_plane(point, *this->_plane))
        {
            bool is_duplicate = false;
            for (const auto& p : _inserted_points)
            {
                const vec3 o = p.second;
                if (o.x == point.x && o.y == point.y)
                {
                    mesh_to_interface[v] = p.first;
                    is_duplicate = true;
                    break;
                }
            }

            if (!is_duplicate)
            {
                g_index interface_vertex_id = this->_constraints.vertices.create_vertex(point.data());
                this->_unique_vertices++;
                mesh_to_interface[v] = interface_vertex_id;
                _inserted_points[interface_vertex_id] = vec3(point.x, point.y, this->_plane->extent);
            }
        }
    }

    for (const g_index v0 : this->_constraints.vertices)
    {
        const auto p0 =_constraints.vertices.point_ptr(v0);
        for (const g_index v1 : this->_constraints.vertices)
        {
            if (v0 == v1)
            {
                continue;
            }

            const auto p1 =_constraints.vertices.point_ptr(v1);
            if (p0[0] == p1[0] && p0[1] == p1[1])
            {
                OOC_DEBUG("equal points " << v0 << ", " << v1 << " at " << p0[0] << ", " << p0[1]);
                break;
            }

        }
    }

    for (const g_index f : mesh.facets)
    {
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

    //if (this->_constraints.edges.nb() != 0)
    //{
    //    return;
    //}

    for (const auto& entry : edge_count_map)
    {
        std::unordered_set<g_index> to_delete;
        // one adjacent interface facet -> edge must be a boundary edge.
        if (entry.second == 1)
        {
            // translate edge-vertices from mesh-global index to interface-local index.
            auto global = entry.first;
            if (mesh_to_interface.contains(global.first) && mesh_to_interface.contains(global.second))
            {
                bool intersected = false;
                for (const g_index e : this->_constraints.edges)
                {
                    if (to_delete.contains(e))
                    {
                        continue;
                    }

                    const vec2 p0 = vec2(this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(e, 0)));
                    const vec2 p1 = vec2(this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(e, 1)));

                    const vec2 cp0 = vec2(this->_constraints.vertices.point_ptr(mesh_to_interface[global.first]));
                    const vec2 cp1 = vec2(this->_constraints.vertices.point_ptr(mesh_to_interface[global.second]));

                    const auto intersection = moist::predicates::xy::get_line_intersection(p0, p1, cp0, cp1, 1e16);

                    if (intersection != std::nullopt)
                    {
                        // create vertex and re-create edges around it
                        const g_index v = this->_constraints.vertices.create_vertex(intersection.value());

                        // split the original edge...
                        //this->_constraints.edges.create_edge(this->_constraints.edges.vertex(e, 0), v);
                        //this->_constraints.edges.create_edge(v, this->_constraints.edges.vertex(e, 1));

                        // create the new edge as split edge
                        //this->_constraints.edges.create_edge(mesh_to_interface[global.first], v);
                        //this->_constraints.edges.create_edge(v, mesh_to_interface[global.second]);

                        //to_delete.insert(e);
                        intersected = true;
                    }
                }

                // if no intersecting edge has been found, just create the edge in the constraints...
                if (!intersected)
                {
                    this->_constraints.edges.create_edge(mesh_to_interface[global.first], mesh_to_interface[global.second]);
                }

                geogram::vector<g_index> deleted(this->_constraints.edges.nb());
                for (const g_index d : to_delete)
                {
                    deleted[d] = true;
                }
                this->_constraints.edges.delete_elements(deleted, false);
            }
            else
            {
                OOC_WARNING("attempt to constrain edge of nonexistent interface vertices...");
            }
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
    geogram::mesh_repair(this->_constraints, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    g_index nb_edges = this->_constraints.edges.nb();
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

#ifndef NDEBUG
    // assign attribute constraint edges to edges in triangulation
    geogram::Attribute<double> e_constrained(this->_triangulation->edges.attributes(), ATTRIBUTE_CONSTRAINT_EDGE);
    for (const g_index edge : this->_triangulation->edges)
    {
        e_constrained[edge] = 0.5;
    }

    OOC_DEBUG("total of " << this->_triangulation->edges.nb() << " edges...");
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
                e_constrained[edge] = 1.0;
            }
        }
    }

    moist::utils::geo::save("interface.constraints.geogram", this->_constraints);
#endif // NDEBUG
    //geogram::Attribute<int> m_target_vertices(this->_triangulation->cells.attributes(), ATTRIBUTE_INTERFACE_TARGET_VERTICES);
    //m_target_vertices[0] = this->_unique_vertices / 2.0;
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
