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

moist::InterfaceGenerator::InterfaceGenerator(const AxisAlignedPlane plane) :
    _constraints(geo::Mesh(2, false)),
    _cdt(CDT::Triangulation<double>(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::TryResolve, 1e-12))
{
    this->_triangulation = std::make_shared<geo::Mesh>(3);
    this->_plane = std::make_shared<AxisAlignedPlane>(plane);
    this->_unique_vertices = 0;
    this->WriteMetadata();
}

static double round16(double value)
{
    const double threshold = 1e-15;
    value = std::round(value * 1e16) / 1e16;
    return (std::abs(value) < threshold) ? 0.0 : value;
}

#ifndef NDEBUG
static geo::Mesh edge_mesh(3);
#endif // NDEBUG

void moist::InterfaceGenerator::AddConstraints(const geo::Mesh &mesh)
{
    std::map<g_index, g_index> mesh_to_interface;
    std::map<std::pair<g_index, g_index>, uint8_t> edge_count_map; // if we have more than 255 incident edges we have a waaaay bigger problem anyway...

    g_index i = _interface_vertices.size();
    for (const g_index v : mesh.vertices)
    {
        // geogram uses precision values up to 16 decimal places...
        vec3 point = mesh.vertices.point(v);

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
                _interface_vertices.push_back(CDT::V2d<double>(point.x, point.y));
                g_index interface_vertex_id = i;
                _inserted_points[interface_vertex_id] = vec3(point.x, point.y, this->_plane->extent);
                mesh_to_interface[v] = interface_vertex_id;
                i++;
                //g_index interface_vertex_id = this->_constraints.vertices.create_vertex(point.data());
                //this->_unique_vertices++;
                //mesh_to_interface[v] = interface_vertex_id;
                //_inserted_points[interface_vertex_id] = vec3(point.x, point.y, this->_plane->extent);
            }
        }
    }

    for (const geo::index_t c : mesh.cells)
    {
        if (!predicates::cell_full_facet_on_plane(c, mesh, *_plane))
        {
            continue;
        }

        for (geo::index_t le = 0; le < mesh.cells.nb_edges(c); le++)
        {
            geo::index_t v1 = mesh.cells.edge_vertex(c, le, 0);
            geo::index_t v2 = mesh.cells.edge_vertex(c, le, 1);

            if (!moist::predicates::edge_on_plane(mesh.vertices.point(v1), mesh.vertices.point(v2), *_plane))
            {
                continue;
            }

            if (v1 > v2)
            {
                std::swap(v1, v2);
            }

            edge_count_map[{v1, v2}]++;
        }
    }

    std::unordered_set<g_index> to_delete;
    for (const auto& entry : edge_count_map)
    {
        // one adjacent interface facet -> edge must be a boundary edge.
        if (entry.second == 1)
        {
            // translate edge-vertices from mesh-global index to interface-local index.
            auto global = entry.first;
            if (mesh_to_interface.contains(global.first) && mesh_to_interface.contains(global.second))
            {
                // "New" Ogre CDT can handle intersections in constrained edges more robustly than triangle, so we can just dump every edge here
                // Update: Nope it cannot ;(
                _interface_edges.push_back(CDT::Edge(mesh_to_interface[global.first], mesh_to_interface[global.second]));
                edge_mesh.edges.create_edge(
                    edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.first]]),
                    edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.second]])
                );

                // Check if our new edge intersects already with an interface edge... if so, create four new edges and one point
                /*bool intersected = false;
                std::size_t nb_edges = _interface_edges.size();
                for (std::size_t e = 0; e < nb_edges; e++)
                {
                    if (to_delete.contains(e))
                    {
                        continue;
                    }

                    const auto edge = _interface_edges.at(e);

                    // Points of the original edge which we want to insert
                    const auto p0 = _interface_vertices.at(mesh_to_interface[global.first]);
                    const auto p1 = _interface_vertices.at(mesh_to_interface[global.second]);

                    const auto existing_p0 = _interface_vertices.at(edge.v1());
                    const auto existing_p1 = _interface_vertices.at(edge.v2());

                    // Edges cannot intersect in a single point if they share an endpoint (besides that one endpoint, which we don't care about)
                    if (p0 == existing_p0 || p0 == existing_p1 || p1 == existing_p0 || p1 == existing_p1)
                    {
                        continue;
                    }

                    const auto geo_p0 = vec2(p0.x, p0.y);
                    const auto geo_p1 = vec2(p1.x, p1.y);
                    const auto geo_p2 = vec2(existing_p0.x, existing_p0.y);
                    const auto geo_p3 = vec2(existing_p1.x, existing_p1.y);
                    const auto intersection = moist::predicates::xy::get_line_intersection(geo_p0, geo_p1, geo_p2, geo_p3, 1e16);

                    if (intersection != std::nullopt)
                    {
                        const auto point = intersection.value();
                        if (geo_p0 == point || geo_p1 == point || geo_p2 == point || geo_p3 == point)
                        {
                            continue;
                        }

                        _interface_vertices.push_back(CDT::V2d<double>(point.x, point.y));
                        g_index v = _interface_vertices.size() - 1;
                        _inserted_points[v] = vec3(point.x, point.y, this->_plane->extent);

                        // split the original edge...
                        _interface_edges.push_back(CDT::Edge(edge.v1(), v));
                        _interface_edges.push_back(CDT::Edge(v, edge.v2()));
                        // create the new edge as split edge
                        this->_constraints.edges.create_edge(mesh_to_interface[global.first], v);
                        this->_constraints.edges.create_edge(v, mesh_to_interface[global.second]);
                        to_delete.insert(e);
                        intersected = true;

                    #ifndef NDEBUG
                        edge_mesh.edges.create_edge(
                            edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.first]]),
                            edge_mesh.vertices.create_vertex(_inserted_points[v])
                        );

                        edge_mesh.edges.create_edge(
                            edge_mesh.vertices.create_vertex(_inserted_points[v]),
                            edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.second]])
                        );
                        edge_mesh.edges.create_edge(
                            edge_mesh.vertices.create_vertex(_inserted_points[edge.v1()]),
                            edge_mesh.vertices.create_vertex(_inserted_points[v])
                        );
                        edge_mesh.edges.create_edge(
                            edge_mesh.vertices.create_vertex(_inserted_points[v]),
                            edge_mesh.vertices.create_vertex(_inserted_points[edge.v2()])
                        );
                    #endif // NDEBUG
                    }
                }

                if (!intersected)
                {
                    _interface_edges.push_back(CDT::Edge(mesh_to_interface[global.first], mesh_to_interface[global.second]));
                #ifndef NDEBUG
                    edge_mesh.edges.create_edge(
                        edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.first]]),
                        edge_mesh.vertices.create_vertex(_inserted_points[mesh_to_interface[global.second]])
                    );
                #endif // NDEBUG
                }


                // for (const g_index e : this->_constraints.edges)
                // {
                //     if (to_delete.contains(e))
                //     {
                //         continue;
                //     }
                //     const vec2 p0 = vec2(this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(e, 0)));
                //     const vec2 p1 = vec2(this->_constraints.vertices.point_ptr(this->_constraints.edges.vertex(e, 1)));
                //     const vec2 cp0 = vec2(this->_constraints.vertices.point_ptr(mesh_to_interface[global.first]));
                //     const vec2 cp1 = vec2(this->_constraints.vertices.point_ptr(mesh_to_interface[global.second]));
                //     if (p0 == cp0 || p0 == cp1 || p1 == cp0 || p1 == cp1)
                //     {
                //         continue;
                //     }
                //     const auto intersection = moist::predicates::xy::get_line_intersection(p0, p1, cp0, cp1, 1e16);
                //     if (intersection != std::nullopt)
                //     {
                //         if (p0 == intersection.value() || p1 == intersection.value() || cp0 == intersection.value() || cp1 == intersection.value())
                //         {
                //             continue;
                //         }
                //         // create vertex and re-create edges around it
                //         const g_index v = this->_constraints.vertices.create_vertex(intersection.value());
                //         // split the original edge...
                //         this->_constraints.edges.create_edge(this->_constraints.edges.vertex(e, 0), v);
                //         this->_constraints.edges.create_edge(v, this->_constraints.edges.vertex(e, 1));
                //         // create the new edge as split edge
                //         this->_constraints.edges.create_edge(mesh_to_interface[global.first], v);
                //         this->_constraints.edges.create_edge(v, mesh_to_interface[global.second]);
                //         //to_delete.insert(e);
                //         intersected = true;
                //     }
                // }

                // if no intersecting edge has been found, just create the edge in the constraints...
                // if (!intersected)
                // {
                // this->_constraints.edges.create_edge(mesh_to_interface[global.first], mesh_to_interface[global.second]);
                // }

                // geo::vector<g_index> deleted(this->_constraints.edges.nb());
                // for (const g_index d : to_delete)
                // {
                //     deleted[d] = true;
                // }
                // this->_constraints.edges.delete_elements(deleted, false);*/
            }
            else
            {
                OOC_WARNING("attempt to constrain edge of nonexistent interface vertices...");
            }
        }
    }

    // CDT does not allow for copying / moving edges in a vector, so do some shenanigans where we delete from the back
    /*std::vector<std::size_t> indices(to_delete.begin(), to_delete.end());
    std::sort(indices.rbegin(), indices.rend());

    CDT::EdgeVec new_edges;
    new_edges.reserve(_interface_edges.size() - to_delete.size()); // upper bound

    for (int idx : indices)
    {
        if (idx >= 0 && static_cast<size_t>(idx) < _interface_edges.size())
        {
            _interface_edges.erase(_interface_edges.begin() + idx);
        }
    }*/

#ifndef NDEBUG
    /*geo::vector<geo::index_t> to_delete_geo(edge_mesh.edges.nb());
    for (geo::index_t e = 0; e < edge_mesh.edges.nb(); e++)
    {
        if (to_delete.contains(e))
        {
            to_delete_geo[e] = true;
        }
    }
    edge_mesh.edges.delete_elements(to_delete_geo);*/
    moist::utils::geogram::save("edge_mesh.obj", edge_mesh);
#endif // NDEBUG
}

static bool vec2_eq_float_precision(const vec2& a, const vec2& b)
{
    constexpr float epsilon = 1e-6f;

    return std::fabs(static_cast<float>(a.x) - static_cast<float>(b.x)) < epsilon &&
           std::fabs(static_cast<float>(a.y) - static_cast<float>(b.y)) < epsilon;
}

void moist::InterfaceGenerator::Triangulate()
{
    _cdt.insertVertices(_interface_vertices);
    _cdt.insertEdges(_interface_edges);
    _cdt.eraseSuperTriangle();
    geo::vector<double> vertices(_cdt.vertices.size() * 3);
    geo::vector<geo::index_t> triangles(_cdt.triangles.size() * 3);

    for (int v = 0; v < _cdt.vertices.size(); v++)
    {
        vertices[3 * v] = _cdt.vertices[v].x;
        vertices[3 * v + 1] = _cdt.vertices[v].y;
        vertices[3 * v + 2] = 0.0;
    }

    for (int f = 0; f < _cdt.triangles.size(); f++)
    {
        triangles[3 * f] = geo::index_t(_cdt.triangles[f].vertices[0]);
        triangles[3 * f + 1] = geo::index_t(_cdt.triangles[f].vertices[1]);
        triangles[3 * f + 2] = geo::index_t(_cdt.triangles[f].vertices[2]);
    }

    // Only the constraint-mesh is required, must pass 0 and nullptr here to set_vertices when using triangle.
    // TODO: REALLY annoying thing about triangle... it seems to create a bbox ("master-triangles") where a vertex of the bbox corresponds with a vertex of the
    //       input, a double vertex is found and "ignored", which leads to a segfault when geogram tries to read the data back into it's
    //       own data structure... somehow create own bbox for triangle?
    // TODO: this doesn't actually help fix things...
    // geo::mesh_repair(this->_constraints, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);

    // g_index nb_edges = this->_constraints.edges.nb();
    // auto delaunay = geo::Delaunay::create(2, "triangle");
    // delaunay->set_constraints(&this->_constraints);
    // delaunay->set_vertices(0, nullptr);

    // OOC_DEBUG("Triangulated #" << delaunay->nb_vertices() << " interface vertices");

    // geo::vector<double> vertices(delaunay->nb_vertices() * 3);
    // for(g_index v = 0; v < delaunay->nb_vertices(); v++)
    // {
    //     vertices[3 * v] = delaunay->vertex_ptr(v)[0];
    //     vertices[3 * v + 1] = delaunay->vertex_ptr(v)[1];
    //     vertices[3 * v + 2] = 0.0;
    // }

    // geo::vector<geo::index_t> triangles(delaunay->nb_cells() * 3);
    // for(g_index t = 0; t < delaunay->nb_cells(); t++)
    // {
    //     triangles[3 * t] = geo::index_t(delaunay->cell_vertex(t, 0));
    //     triangles[3 * t + 1] = geo::index_t(delaunay->cell_vertex(t, 1));
    //     triangles[3 * t + 2] = geo::index_t(delaunay->cell_vertex(t, 2));
    // }

    this->_triangulation->facets.assign_triangle_mesh((geo::coord_index_t) 3, vertices, triangles, true);
    this->_triangulation->facets.compute_borders();

#ifndef NDEBUG
    // assign attribute constraint edges to edges in triangulation
    geo::Attribute<bool> v_constrained(this->_triangulation->vertices.attributes(), ATTRIBUTE_CONSTRAINT_EDGE);
    for (const g_index edge : this->_triangulation->edges)
    {
        v_constrained[_triangulation->edges.vertex(edge, 0)] = false;
        v_constrained[_triangulation->edges.vertex(edge, 1)] = false;
    }

    OOC_DEBUG("total of " << this->_triangulation->edges.nb() << " edges...");
    for (int e = 0; e < _interface_edges.size(); e++)
    {
        // TODO: This also only works with z-growing...
        const CDT::Edge constraint_edge = _interface_edges[e];

        const vec2 cep0 = vec2(_interface_vertices[constraint_edge.v1()].x, _interface_vertices[constraint_edge.v1()].y);
        const vec2 cep1 = vec2(_interface_vertices[constraint_edge.v2()].x, _interface_vertices[constraint_edge.v2()].y);

        // find edge in interface mesh, and mark it
        for (const g_index edge : this->_triangulation->edges)
        {
            const vec2 ep0 = reinterpret_cast<const vec2&>(_triangulation->vertices.point(this->_triangulation->edges.vertex(edge, 0)));
            const vec2 ep1 = reinterpret_cast<const vec2&>(_triangulation->vertices.point(this->_triangulation->edges.vertex(edge, 1)));

            //if (ep0 == cep0 && ep1 == cep1 || ep0 == cep1 && ep1 == cep0)
            if ((vec2_eq_float_precision(ep0, cep0) && vec2_eq_float_precision(ep1, cep1)) || (vec2_eq_float_precision(ep0, cep1) && vec2_eq_float_precision(ep1, cep0)))
            {
                v_constrained[_triangulation->edges.vertex(edge, 0)] = true;
                v_constrained[_triangulation->edges.vertex(edge, 1)] = true;
            }
        }
    }

    moist::utils::geogram::save("interface.constraints.geogram", this->_constraints);
#endif // NDEBUG
    //geogram::Attribute<int> m_target_vertices(this->_triangulation->cells.attributes(), ATTRIBUTE_INTERFACE_TARGET_VERTICES);
    //m_target_vertices[0] = this->_unique_vertices / 2.0;

    // this->Decimate();
}

static vec3 middle(const vec3 v0, const vec3 v1)
{
    return vec3(
        (v0.x + v1.x) / 2.0,
        (v0.y + v1.y) / 2.0,
        (v0.z + v1.z) / 2.0
    );
}

static std::pair<geo::index_t, geo::index_t> find_shortest_edge(const geo::index_t& f, const geo::Mesh& mesh)
{
    double shortest = std::numeric_limits<double>::max();
    std::pair<geo::index_t, geo::index_t> e_shortest;
    for (geo::index_t le = 0; le < 3; le++)
    {
        const geo::index_t v0 = mesh.facets.vertex(f, le);
        const geo::index_t v1 = mesh.facets.vertex(f, (le + 1) % 3);

        const double d = geo::distance(mesh.vertices.point(v0), mesh.vertices.point(v1));
        if (d <= shortest)
        {
            shortest = d;
            e_shortest = std::make_pair(v0, v1);
        }
    }

    return e_shortest;
}

void moist::InterfaceGenerator::Decimate()
{
    size_t n = (this->_triangulation->facets.nb() - this->_constraints.edges.nb()) / 2;
    OOC_DEBUG("decimating " << n << " worst unconstrained facets");

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

    geo::Attribute<bool> v_constrained(this->_triangulation->vertices.attributes(), ATTRIBUTE_CONSTRAINT_EDGE);
    for (const geo::index_t f : triangles)
    {
        if (n == 0)
        {
            break;
        }

        const auto e = find_shortest_edge(f, *triangulation);

        // Check if the edge is decimatable
        if (v_constrained[e.first] && v_constrained[e.second])
        {
            continue;
        }

        if (v_constrained[e.first])
        {
            triangulation->vertices.point(e.second) = triangulation->vertices.point(e.first);
            v_constrained[e.second] = true;
        }
        else if (v_constrained[e.second])
        {
            triangulation->vertices.point(e.first) = triangulation->vertices.point(e.second);
            v_constrained[e.first] = true;
        }
        else
        {
            // idk just merge in middle
            const vec3 mid = middle(triangulation->vertices.point(e.first), triangulation->vertices.point(e.second));
            triangulation->vertices.point(e.first) = mid;
            triangulation->vertices.point(e.second) = mid;
            v_constrained[e.first] = true;
            v_constrained[e.second] = true;
        }

        n--;
    }

    // decimate shortest edge of n triangles
    /*for (size_t i = 0; i < n; i++)
    {
        const g_index f = triangles[i];
        l_index shortest_edge = -1;
        double shortest_edge_length = std::numeric_limits<double>::max();
        for (l_index le = 0; le < 3; le++)
        {
            const vec3 p0 = triangulation->vertices.point(triangulation->facets.vertex(f, le));
            const vec3 p1 = triangulation->vertices.point(triangulation->facets.vertex(f, (le + 1) % 3));
            const double distance = geo::distance(p0, p1);
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
    }*/

    // delete all resulting 0-volume triangles
    geo::vector<g_index> deleted_facets(triangulation->facets.nb());
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
