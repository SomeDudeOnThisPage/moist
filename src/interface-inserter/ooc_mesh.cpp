#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>
#include <mutex>
#include <limits>

#include "local_operations.hpp"
#include "ooc_mesh.hpp"
#include "attributes.inl"
#include "predicates.inl"
#include "geometry.inl"
#include "utils.hpp"

typedef struct EDGE_STRUCT
{
    g_index v0;
    g_index v1;

    bool operator==(const EDGE_STRUCT& other) const
    {
        return (v0 == other.v0 && v1 == other.v1) || (v0 == other.v1 && v1 == other.v0);
    }
} Edge;

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

void incremental_meshing::MeshSlice::InsertInterfaceVertices(incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface vertices");

    geogram::Attribute<incremental_meshing::InterfaceVertexStrategy> v_strategy(this->vertices.attributes(), incremental_meshing::INTERFACE_VERTEX_STRATEGY_ATTRIBUTE);
    geogram::Attribute<double> v_test(this->vertices.attributes(), incremental_meshing::INTERFACE);

    // mark vertices as discardable if they lie on the interface plane (later, during vertex insertion, vertices that match the interface will
    // be marked as interface, and thus non-discardable...)
    for (const g_index v : this->vertices)
    {
        std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR);
        geogram::Attribute<incremental_meshing::attributes::VertexDescriptorFlags> v_descriptor(this->vertices.attributes(), ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS);
        v_descriptor[v] |= incremental_meshing::predicates::point_on_plane(this->vertices.point(v), *interface.Plane())
            ? incremental_meshing::attributes::VertexDescriptorFlags::DISCARD
            : incremental_meshing::attributes::VertexDescriptorFlags::NONE;
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
    incremental_meshing::utils::dump_mesh(*this, "_after_point_insertion.msh");
#endif // NDEBUG
}

void incremental_meshing::MeshSlice::InsertInterfaceEdges(incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface edges into mesh '" + _identifier + "'");

    auto triangulation = interface.Triangulation();
    auto plane = interface.Plane();

    // sometimes edges are just... missing... even when calling compute_borders, so do it manually here.
    std::unordered_set<Edge> edges;
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
                crossed_edges.push_back({ v0, v1, p });
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
    incremental_meshing::utils::dump_mesh(*this, "_after_edge_insertion.msh");
#endif // NDEBUG
}

void incremental_meshing::MeshSlice::CreateTetrahedra()
{
    //OOC_DEBUG("creating #" << this->_created_tets.size() << " tets");
    for (const auto tet : this->_created_tets)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
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
    this->cells.delete_elements(tets_to_delete, true);
    _deleted_tets.clear();
}

void incremental_meshing::MeshSlice::DecimateNonInterfaceEdges(incremental_meshing::Interface& interface)
{
    geogram::Attribute<incremental_meshing::InterfaceVertexStrategy> v_strategy(this->vertices.attributes(), incremental_meshing::INTERFACE_VERTEX_STRATEGY_ATTRIBUTE);
    // map positions onto a list of vertices
    std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator> map;

    for (g_index v : this->vertices)
    {
        const auto strat = v_strategy[v];
        if (!incremental_meshing::predicates::point_on_plane(this->vertices.point(v), *interface.Plane()))
        {
            continue;
        }

        map[this->vertices.point(v)] = std::unordered_set<g_index>();
        map[this->vertices.point(v)].insert(v);
    }

    int collapsed = 0;
    const g_index size = this->cells.nb();
    for (g_index c = 0; c < size; c++)
    {
        for (l_index lv = 0; lv < 4; lv++)
        {
            const g_index v = this->cells.vertex(c, lv);
            if (v_strategy[v] == incremental_meshing::InterfaceVertexStrategy::KEEP || !incremental_meshing::predicates::point_on_plane(this->vertices.point(v), *interface.Plane()))
            {
                continue;
            }

            // find vertex to move onto
            g_index v_collapse_onto = -1;
            bool keep;
            for (l_index to_lv = (lv + 1) % 4; to_lv != lv; to_lv = (to_lv + 1) % 4)
            {
                const g_index to_v = this->cells.vertex(c, to_lv);
                if (!incremental_meshing::predicates::point_on_plane(this->vertices.point(to_v), *interface.Plane()))
                {
                    continue;
                }

                v_collapse_onto = to_v;

                if (v_strategy[to_v] == incremental_meshing::InterfaceVertexStrategy::KEEP)
                {
                    keep = true;
                    break;
                }
            }

            if (v_collapse_onto == -1)
            {
                continue;
            }

            const vec3 v_keep_point = this->vertices.point(v_collapse_onto);
            for (const auto move_v : map[this->vertices.point(v)])
            {
                const auto ll = map[this->vertices.point(v)];
                if (keep)
                {
                    v_strategy[move_v] = incremental_meshing::InterfaceVertexStrategy::KEEP;
                }

                // OOC_DEBUG("moving " << this->vertices.point(move_v) << " -> " << v_keep_point);
                // actually move all of the points
                this->vertices.point_ptr(move_v)[0] = v_keep_point.x;
                this->vertices.point_ptr(move_v)[1] = v_keep_point.y;
                this->vertices.point_ptr(move_v)[2] = v_keep_point.z;
                map[v_keep_point].insert(move_v);
                map[this->vertices.point(move_v)].clear();
            }

        }

        collapsed++;
    }

    OOC_DEBUG("collapsed " << collapsed << " tets");
    for (const auto cell : this->cells)
    {
        if (geogram::mesh_cell_volume(*this, cell) == 0 || incremental_meshing::geometry::has_duplicate_vertex(cell, *this))
        {
            _deleted_tets.insert(cell);
        }
    }
    this->FlushTetrahedra();

#ifndef NDEBUG
    incremental_meshing::utils::dump_mesh(*this, "_after_clustering.msh");
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
            incremental_meshing::operation::vertex_insert_2to4(*this, cell, point, plane);
            PARALLEL_BREAK;
        }
    }
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    );
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    this->CreateTetrahedra();
}
