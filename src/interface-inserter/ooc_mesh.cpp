#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>
#include <limits>
#include <mutex>

#include "local_operations.hpp"
#include "ooc_mesh.hpp"
#include "attributes.inl"
#include "predicates.inl"
#include "geometry.inl"
#include "utils.hpp"

// COMPLETELY DIFFERENT IDEA:
// FOR EACH INTERSECTION BETWEEN AN EXISTING EDGE AND AN INTERFACE EDGE, INSERT ONE TET THAT
// HAS ONE OF THE INTERFACE EDGE-VERTICES AND THE OTHERS OF THE "ORIGINAL" TET THE CROSSING LIES IN
// THAT WAY WE CAN INSTANTLY CREATE THE NECCESSARY TETRAHEDRA???
// Downside: Not really easily parallelizable
// Also not difficultly parallelizable

typedef struct EDGE_STRUCT
{
    g_index v0;
    g_index v1;

    g_index e;

    bool operator==(const EDGE_STRUCT& other) const
    {
        return (v0 == other.v0 && v1 == other.v1) || (v0 == other.v1 && v1 == other.v0);
    }
} Edge;

typedef struct EDGE_STRUCT_POSITION
{
    vec2 p0;
    vec2 p1;

    bool operator==(const EDGE_STRUCT_POSITION& other) const
    {
        return (p0 == other.p0 && p1 == other.p1) || (p0 == other.p1 && p1 == other.p0);
    }
} EdgePosition;

struct Vec3HashOperator
{
    std::size_t operator()(const vec3& v) const
    {
        std::hash<double> hasher;
        std::size_t hx = hasher(v.x);
        std::size_t hy = hasher(v.y);
        std::size_t hz = hasher(v.z);

        return hx ^ (hy << 1) ^ (hz << 2);
    }
};

struct Vec3EqualOperator
{
    bool operator()(const GEO::vecng<3, double>& a, const GEO::vecng<3, double>& b) const
    {
        return std::fabs(a.x - b.x) < incremental_meshing::__DOUBLE_EPSILON &&
            std::fabs(a.y - b.y) < incremental_meshing::__DOUBLE_EPSILON &&
            std::fabs(a.z - b.z) < incremental_meshing::__DOUBLE_EPSILON;
    }
};

namespace std
{
    template <>
    struct hash<vec2>
    {
        std::size_t operator()(const vec2& v) const
        {
            std::hash<float> hasher;
            std::size_t h1 = hasher(v.x);
            std::size_t h2 = hasher(v.y);
            return h1 ^ (h2 << 1);
        }
    };

    template <>
    struct hash<incremental_meshing::CROSSED_EDGE_STRUCT>
    {
        std::size_t operator()(const incremental_meshing::CROSSED_EDGE_STRUCT& edge) const
        {
            return std::hash<geogram::index_t>()(std::min(edge.e_v0, edge.e_v1)) ^ (std::hash<geogram::index_t>()(std::max(edge.e_v0, edge.e_v1)) << 1);
        }
    };

    template <>
    struct hash<EDGE_STRUCT>
    {
        std::size_t operator()(const EDGE_STRUCT& edge) const
        {
            return std::hash<geogram::index_t>()(std::min(edge.v0, edge.v1)) ^ (std::hash<geogram::index_t>()(std::max(edge.v0, edge.v1)) << 1);
        }
    };

    template <>
    struct hash<EDGE_STRUCT_POSITION>
    {
        std::size_t operator()(const EDGE_STRUCT_POSITION& e) const
        {
            return std::hash<vec2>()(e.p0) ^ std::hash<vec2>()(e.p1);
        }
    };

}

incremental_meshing::MeshSlice::MeshSlice(geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision)
{
    _deleted_tets = std::unordered_set<geogram::index_t>();
}

void incremental_meshing::MeshSlice::CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra)
{
    for (const CreatedTetrahedon tet : tetrahedra)
    {
        _created_tets.push_back(tet);
    }
}

void incremental_meshing::MeshSlice::DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra)
{
    for (const g_index tet : tetrahedra)
    {
        _deleted_tets.insert(tet);
    }
}

void incremental_meshing::MeshSlice::InsertInterface(incremental_meshing::Interface& interface)
{
    this->InsertInterfaceVertices(interface);
    this->InsertInterfaceEdges(interface);
    this->DecimateNonInterfaceEdges(interface);
}

// incremental approach... much less parallelizable... need to lock affected tetrahedra... somehow sort into edge-pools that have no common overlap?
// maybe only parallelize local operations... like searching tetrahedra, etc...
//#define EDGE_VEC3(t, v) t->vertices.point(v)
/*void incremental_meshing::MeshSlice::InsertInterface(incremental_meshing::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();

    std::unordered_set<Edge> edges; // set to avoid duplicates... TODO: make this inline loop and just check set exists...
    for (const g_index f_id : triangulation->facets)
    {
        for (int i = 0; i < 3; i++)
        {
            edges.emplace(Edge {
                triangulation->facets.vertex(f_id, i),
                triangulation->facets.vertex(f_id, (i + 1) % 3)
            });
        }
    }

    uint32_t i = 0;
    for (const auto edge : edges)
    {
        this->_created_tets_idx.clear();

        vec3 p0 = EDGE_VEC3(triangulation, edge.v0);
        vec3 p1 = EDGE_VEC3(triangulation, edge.v0);
        p0.z = plane->extent;
        p1.z = plane->extent;

        this->InsertVertex(p0, *plane);
        this->InsertVertex(p1, *plane);

        // find all tetrahedra that lie on the line between the two points
        const auto cells = this->cells.nb();
        std::vector<g_index> crossed_cells;
        for (g_index cell = 0; cell < cells; cell++)
        {
            if (_deleted_tets.contains(cell))
            {
                PARALLEL_CONTINUE;
            }

            // insert vertices where the line crosses edges of the tetrahedra
            std::vector<CrossedEdge> crossed_edges;
            std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices;
            for (g_index e_id = 0; e_id < this->cells.nb_edges(cell); e_id++)
            {
                const g_index v0 = this->cells.edge_vertex(cell, e_id, 0);
                const g_index v1 = this->cells.edge_vertex(cell, e_id, 1);
                const vec3 cp0 = this->vertices.point(v0);
                const vec3 cp1 = this->vertices.point(v1);

                if (!incremental_meshing::predicates::edge_on_plane(p0, p1, *plane))
                {
                    continue;
                }

                // internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3
                if (!incremental_meshing::predicates::xy::check_lines_aabb(
                    reinterpret_cast<const vec2&>(cp0),
                    reinterpret_cast<const vec2&>(cp1),
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1)
                ))
                {
                    continue;
                }

                const auto intersection_opt = incremental_meshing::predicates::xy::get_line_intersection(
                    reinterpret_cast<const vec2&>(cp0),
                    reinterpret_cast<const vec2&>(cp1),
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1)
                );

                if (!intersection_opt.has_value())
                {
                    continue;
                }

                const vec3 new_vertex = vec3(intersection_opt.value().x, intersection_opt.value().y, plane->extent);
                const g_index p = created_vertices.find(new_vertex) != created_vertices.end()
                    ? created_vertices.at(new_vertex)
                    : this->vertices.create_vertex(new_vertex.data());

                created_vertices.emplace(new_vertex, p);
                crossed_edges.push_back({ v0, v1, p });
            }

            switch (crossed_edges.size())
            {
                case 1:
                    incremental_meshing::operation::edge_split_1to2(*this, cell, crossed_edges[0], *plane);
                    // TODO: directly add the new shitty cells to crossed_cells since they need to be collapsed anyway
                    break;
                case 2:
                    // TODO: directly add the new shitty cells to crossed_cells since they need to be collapsed anyway
                    incremental_meshing::operation::edge_split_1to3(*this, cell, crossed_edges[0], crossed_edges[1], *plane);
                    break;
            default:
                    break;
            }

            this->CreateTetrahedra();
        }

        // decimate unwanted tetrahedra by collapsing along the line
        for (const g_index created_cell : this->_created_tets_idx)
        {
            // find vertex that is on the line
            for (l_index lv = 0; lv < 4; lv++)
            {
                LOCK_ATTRIBUTES;
                geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
                const g_index v = this->cells.vertex(created_cell, lv);
                if (!v_discard[v])
                {
                    continue;
                }

                // move onto the closer edge vertex
                const auto d_e0 = geogram::Geom::distance(p0, this->vertices.point(v));
                const auto d_e1 = geogram::Geom::distance(p1, this->vertices.point(v));
                if (d_e0 < d_e1)
                {
                    this->vertices.point(v).x = p0.x;
                    this->vertices.point(v).y = p0.y;
                    this->vertices.point(v).z = p0.z;
                }
                else
                {
                    this->vertices.point(v).x = p1.x;
                    this->vertices.point(v).y = p1.y;
                    this->vertices.point(v).z = p1.z;
                }
            }
        }

        for (const auto cell : this->cells)
        {
            if (geogram::mesh_cell_volume(*this, cell) == 0 || incremental_meshing::geometry::has_duplicate_vertex(cell, *this))
            {
                _deleted_tets.insert(cell);
            }
        }

        this->FlushTetrahedra();

    #ifndef NDEBUG
        if (i == 0 || i == 5 || i == 10)
        {
            incremental_meshing::utils::dump_mesh(*this, std::format("test/target/_after_edge_insert_{}.geogram", i++));
        }
    #endif
    }
}*/

void incremental_meshing::MeshSlice::InsertInterfaceVertices(incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface vertices");

    // mark vertices as discardable if they lie on the interface plane (later, during vertex insertion, vertices that match the interface will
    // be marked as interface, and thus non-discardable...)
    for (const g_index v : this->vertices)
    {
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
        geogram::Attribute<int> v_interface(this->vertices.attributes(), ATTRIBUTE_INTERFACE);
        geogram::Attribute<int> v_cluster_direction(this->vertices.attributes(), ATTRIBUTE_CLUSTER_ONTO);

        v_discard[v] = true;
        v_interface[v] = incremental_meshing::predicates::point_on_plane(this->vertices.point(v), *interface.Plane());
        v_cluster_direction[v] = -1;
    }

    const auto triangulation = interface.Triangulation();
    for (auto v_id : triangulation->vertices)
    {
        this->InsertVertex(
            geogram::vec3(triangulation->vertices.point(v_id).x, triangulation->vertices.point(v_id).y, interface.Plane()->extent),
            *interface.Plane()
        );
    }

    this->FlushTetrahedra();
    TIMER_END;

#ifndef NDEBUG
    OOC_DEBUG("Validating point insertion...");
    for (const auto interface_vertex : triangulation->vertices)
    {
        auto interface_point = triangulation->vertices.point(interface_vertex);
        bool found = false;
        for (const auto mesh_vertex : this->vertices)
        {
            if (incremental_meshing::predicates::vec_eq_2d(this->vertices.point(mesh_vertex), interface_point, *interface.Plane()))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            OOC_WARNING("Missing point " << interface_point << " in mesh " << this->_identifier);
        }
    }
    OOC_DEBUG("Done validating point insertion...");
    incremental_meshing::utils::dump_mesh(*this, "after_point_insertion.geogram");
#endif // NDEBUG
}

void incremental_meshing::MeshSlice::InsertInterfaceEdges(incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface edges into mesh '" + _identifier + "'");

    auto triangulation = interface.Triangulation();
    auto plane = interface.Plane();

    // sometimes edges are just... missing... even when calling compute_borders, so do it manually here.
    // note: this needs to be backreferenced and transformed like this again when decimating edges!!!
    std::unordered_set<Edge> edges;
    g_index e = 0;
    for (const g_index f : triangulation->facets)
    {
        for (l_index le = 0; le < 3; le++)
        {
            edges.emplace(Edge {
                triangulation->facets.vertex(f, le),
                triangulation->facets.vertex(f, (le + 1) % 3),
                e++ // custom global edge indexing, MUST be represented by iterating facets and then le from 0->3 with p0 = f.v(f,le); p1 ? f.v(f,(le + 1) % 3)
            });
        }
    }

    for (const auto interface_edge : edges)
    {
        const auto e_p0 = triangulation->vertices.point(interface_edge.v0);
        const auto e_p1 = triangulation->vertices.point(interface_edge.v1);

        const auto cells = this->cells.nb();
        std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices;
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        geogram::parallel_for(0, cells, [&](const g_index c_id)
        {
#else
        for (g_index c_id = 0; c_id < cells; c_id++)
        {
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS
            if (_deleted_tets.contains(c_id))
            {
                PARALLEL_CONTINUE;
            }

            std::vector<CrossedEdge> crossed_edges;
            for (g_index e_id = 0; e_id < this->cells.nb_edges(c_id); e_id++)
            {
                const g_index v0 = this->cells.edge_vertex(c_id, e_id, 0);
                const g_index v1 = this->cells.edge_vertex(c_id, e_id, 1);
                const vec3 p0 = this->vertices.point(v0);
                const vec3 p1 = this->vertices.point(v1);

                if (!incremental_meshing::predicates::edge_on_plane(p0, p1, *plane))
                {
                    continue;
                }

                // internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3
                if (!incremental_meshing::predicates::xy::check_lines_aabb(
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1),
                    reinterpret_cast<const vec2&>(e_p0),
                    reinterpret_cast<const vec2&>(e_p1)
                ))
                {
                    continue;
                }

                const auto intersection_opt = incremental_meshing::predicates::xy::get_line_intersection(
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1),
                    reinterpret_cast<const vec2&>(e_p0),
                    reinterpret_cast<const vec2&>(e_p1)
                );

                if (!intersection_opt.has_value())
                {
                    continue;
                }

                const vec3 new_vertex = vec3(intersection_opt.value().x, intersection_opt.value().y, plane->extent);
                const g_index p = created_vertices.find(new_vertex) != created_vertices.end()
                    ? created_vertices.at(new_vertex)
                    : this->vertices.create_vertex(new_vertex.data());

                created_vertices.emplace(new_vertex, p);
                crossed_edges.push_back({ v0, v1, p, interface_edge.e });
            }

            switch (crossed_edges.size())
            {
                case 1:
                    incremental_meshing::operation::edge_split_1to2(*this, c_id, crossed_edges[0], *plane);
                    break;
                case 2:
                    incremental_meshing::operation::edge_split_1to3(*this, c_id, crossed_edges[0], crossed_edges[1], *plane);
                    break;
               default:
                    break;
            }
        }
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        );
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

        this->CreateTetrahedra();
    }

    this->FlushTetrahedra();

    TIMER_END;

#ifndef NDEBUG
    incremental_meshing::utils::dump_mesh(*this, "after_edge_insertion.geogram");
#endif // NDEBUG
}

void incremental_meshing::MeshSlice::CreateTetrahedra()
{
    //OOC_DEBUG("creating #" << this->_created_tets.size() << " tets");
    for (const auto tet : this->_created_tets)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
        _created_tets_idx.push_back(t);
    #ifndef NDEBUG
        const auto volume = geogram::mesh_cell_volume(*this, t);
        if (volume <= 0.00000) // this happens because I insert existing vertices from the triangulation into the tetmesh...
        {
            const auto p0 = this->vertices.point(tet.v0);
            const auto p1 = this->vertices.point(tet.v1);
            const auto p2 = this->vertices.point(tet.v2);
            const auto p3 = this->vertices.point(tet.v3);
            OOC_WARNING("cell t0 " << t << " has zero volume: " << volume);
            _deleted_tets.insert(t);
        }
    #endif // NDEBUG
    }

    this->_created_tets.clear();
}

void incremental_meshing::MeshSlice::FlushTetrahedra()
{
    OOC_DEBUG("flushing " << _deleted_tets.size() << " out of " << this->cells.nb() << " tets");
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, false);
    _deleted_tets.clear();
}

void incremental_meshing::MeshSlice::DecimateNonInterfaceEdges(incremental_meshing::Interface& interface)
{

#ifndef NDEBUG
{
    int discardable = 0;
    LOCK_ATTRIBUTES;
    geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
    for (const g_index v : this->vertices)
    {
        if (v_discard[v])
        {
            discardable++;
        }
    }
    OOC_DEBUG(discardable << " discardable vertices");
}
#endif // NDEBUG

    // TODO: build this once... need it at multiple places
    // sometimes edges are just... missing... even when calling compute_borders, so do it manually here.
    std::unordered_set<EdgePosition> edges;
    for (const g_index f_id : interface.Triangulation()->facets)
    {
        for (int i = 0; i < 3; i++)
        {
            edges.emplace(EdgePosition {
                reinterpret_cast<const vec2&>(interface.Triangulation()->vertices.point(interface.Triangulation()->facets.vertex(f_id, i))),
                reinterpret_cast<const vec2&>(interface.Triangulation()->vertices.point(interface.Triangulation()->facets.vertex(f_id, (i + 1) % 3)))
            });
        }
    }

    // create buckets of all collapsable vertices containing their index
    std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator> map;
    for (const g_index v : this->vertices)
    {
#ifndef NDEBUG // sanity check... this should not happen as we map created vertices when splitting edges...
        for (const g_index vo : this->vertices)
        {
            if (vo != v && this->vertices.point(v) == this->vertices.point(vo))
            {
                OOC_DEBUG("colocated vertices " << v << ", " << vo << " at " << this->vertices.point(v));
            }
        }
#endif // NDEBUG
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);

        if (v_discard[v])
        {
            const vec3 point = this->vertices.point(v);
            map[point].insert(v);
        }
    }

    std::unordered_map<g_index, std::unordered_set<g_index>> v_to_c; // vertex to cell incidence array only containing interface relevant elements
    for (g_index cell = this->cells.nb() - 1; cell > 0; cell--)
    {
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_interface(this->vertices.attributes(), ATTRIBUTE_INTERFACE);
        geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
        for (const g_index v : incremental_meshing::geometry::cell_vertices(cell, *this))
        {
            if (v == -1)
            {
                OOC_DEBUG("invalid vertex...");
            }
            // TODO: replace this all with bitflags...
            if (!incremental_meshing::predicates::point_on_plane(this->vertices.point(v), *interface.Plane()))
            {
                continue;
            }

            v_to_c[v].insert(cell);
        }
    }
    OOC_DEBUG("#v_to_c = " << v_to_c.size());

    // iterate all vertices in our incidence array
    // if the vertex needs to be collapsed, search incident cells for a interface-nondiscardable vertex.
    for (const auto& [v, _] : v_to_c)
    {
        LOCK_ATTRIBUTES;
        geogram::Attribute<int> v_interface(this->vertices.attributes(), ATTRIBUTE_INTERFACE);
        geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
        if (!v_discard[v])
        {
            continue;
        }

        g_index actual_cluster_onto = -1; // if we have no interface vertices in any cells around the current vertex, cluster onto any...
        g_index cluster_onto = -1;
        for (const g_index cell : v_to_c[v])
        {
            for (const g_index v_to : incremental_meshing::geometry::cell_vertices(cell, *this))
            {
                if (v_to == -1 || !incremental_meshing::predicates::point_on_plane(this->vertices.point(v_to), *interface.Plane()))
                {
                    // OOC_DEBUG("invalid vertex...");
                    continue;
                }

                if (actual_cluster_onto == -1)
                {
                    actual_cluster_onto = v_to;
                }

                // prioritize collapse along interface edge if possible
                if (edges.contains(EdgePosition {
                    reinterpret_cast<const vec2&>(this->vertices.point(v)),
                    reinterpret_cast<const vec2&>(this->vertices.point(v_to))
                }))
                {
                    cluster_onto = v_to;
                }
            }
        }

        actual_cluster_onto = (cluster_onto != -1) ? cluster_onto : actual_cluster_onto;

        if (actual_cluster_onto == -1)
        {
            OOC_DEBUG("invalid vertex");
            continue;
        }

        const vec3 point_to = this->vertices.point(actual_cluster_onto);
        const vec3 point_from = this->vertices.point(v);

        for (const g_index v_movable : map[point_from])
        {
            this->vertices.point(v_movable).x = point_to.x;
            this->vertices.point(v_movable).y = point_to.y;
            this->vertices.point(v_movable).z = point_to.z;

            map[point_to].insert(v_movable);

            if (actual_cluster_onto == cluster_onto)
            {
                LOCK_ATTRIBUTES;
                geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
                geogram::Attribute<int> v_interface(this->vertices.attributes(), ATTRIBUTE_INTERFACE);
                v_discard[v_movable] = false;
                v_interface[v_movable] = true;
            }
        }
    }

    for (const auto cell : this->cells)
    {
        if (geogram::mesh_cell_volume(*this, cell) == 0 || incremental_meshing::geometry::has_duplicate_vertex(cell, *this))
        {
            _deleted_tets.insert(cell);
        }
    }
    this->FlushTetrahedra();

#ifndef NDEBUG
    incremental_meshing::utils::dump_mesh(*this, "after_clustering.geogram");
#endif // NDEBUG
}

static std::mutex _m_deleted_tets;
void incremental_meshing::MeshSlice::InsertVertex(const geogram::vec3& point, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    geogram::parallel_for(0, this->cells.nb(), [this, point, plane](const g_index cell)
    {
#else
    for (const g_index cell : this->cells)
    {
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

        {
            std::lock_guard<std::mutex> lock(_m_deleted_tets);
            if (_deleted_tets.contains(cell))
            {
                PARALLEL_CONTINUE;
            }
        }

        if (incremental_meshing::predicates::point_in_tet(*this, cell, point))
        {
            incremental_meshing::operation::vertex_insert_1to3(*this, cell, point, plane);
            PARALLEL_BREAK;
        }
    }
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    );
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    this->CreateTetrahedra();
}
