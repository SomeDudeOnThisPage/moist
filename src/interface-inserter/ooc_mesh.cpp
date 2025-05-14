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
        return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
            std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON &&
            std::fabs(a.z - b.z) < moist::__DOUBLE_EPSILON;
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
    struct hash<moist::CROSSED_EDGE_STRUCT>
    {
        std::size_t operator()(const moist::CROSSED_EDGE_STRUCT& edge) const
        {
            return std::hash<geogram::index_t>()(std::min(edge.e_v0, edge.e_v1)) ^ (std::hash<geogram::index_t>()(std::max(edge.e_v0, edge.e_v1)) << 1);
        }
    };

}

moist::MeshSlice::MeshSlice(geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision)
{
    _deleted_tets = std::unordered_set<geogram::index_t>();
}

void moist::MeshSlice::CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra)
{
    for (const CreatedTetrahedon tet : tetrahedra)
    {
        _created_tets.push_back(tet);
    }
}

void moist::MeshSlice::DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra)
{
    for (const g_index tet : tetrahedra)
    {
        _deleted_tets.insert(tet);
    }
}

// incremental approach... much less parallelizable... need to lock affected tetrahedra... somehow sort into edge-pools that have no common overlap?
// maybe only parallelize local operations... like searching tetrahedra, etc...
#define EDGE_VEC3(t, v) t->vertices.point(v)
void moist::MeshSlice::InsertInterface(moist::Interface& interface)
{
    this->InsertInterfaceVertices(interface);
    this->InsertInterfaceEdges(interface);
}

void moist::MeshSlice::InsertInterfaceVertices(moist::Interface& interface)
{
    TIMER_START("insert interface vertices");
    const auto triangulation = interface.Triangulation();
    for (const g_index v : triangulation->vertices)
    {
        const vec3 p = this->vertices.point(v);
        std::mutex m_deleted_tets;
    #ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        geogram::parallel_for(0, this->cells.nb(), [this, point, plane](const g_index cell)
        {
    #else
        for (const g_index cell : this->cells)
        {
    #endif // OPTION_PARALLEL_LOCAL_OPERATIONS

            {
                std::lock_guard<std::mutex> lock(m_deleted_tets);
                if (_deleted_tets.contains(cell))
                {
                    PARALLEL_CONTINUE;
                }
            }

            if (moist::predicates::point_in_tet(*this, cell, p, true))
            {
                moist::operation::vertex_insert_1to3(*this, cell, p, *interface.Plane());
                PARALLEL_BREAK;
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
    OOC_DEBUG("Validating point insertion...");
    for (const auto interface_vertex : triangulation->vertices)
    {
        auto interface_point = triangulation->vertices.point(interface_vertex);
        bool found = false;
        for (const auto mesh_vertex : this->vertices)
        {
            if (moist::predicates::vec_eq_2d(this->vertices.point(mesh_vertex), interface_point, *interface.Plane()))
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
    moist::utils::dump_mesh(*this, "after_point_insertion.geogram");
#endif // NDEBUG
}

void moist::MeshSlice::InsertInterfaceEdges(moist::Interface& interface)
{

    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();
    const auto edges = moist::geometry::collect_edges(*triangulation);

    uint32_t i = 0;
    for (const auto edge : edges)
    {
        this->_created_tets_idx.clear();

        vec3 p0 = EDGE_VEC3(triangulation, edge.v0);
        vec3 p1 = EDGE_VEC3(triangulation, edge.v1);
        p0.z = plane->extent;
        p1.z = plane->extent;

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

                if (!moist::predicates::edge_on_plane(cp0, cp1, *plane))
                {
                    continue;
                }

                // internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3
                if (!moist::predicates::xy::check_lines_aabb(
                    reinterpret_cast<const vec2&>(cp0),
                    reinterpret_cast<const vec2&>(cp1),
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1)
                ))
                {
                    continue;
                }

                const auto intersection_opt = moist::predicates::xy::get_line_intersection(
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
                    moist::operation::edge_split_1to2(*this, cell, crossed_edges[0], *plane);
                    // TODO: directly add the new shitty cells to crossed_cells since they need to be collapsed anyway
                    break;
                case 2:
                    // TODO: directly add the new shitty cells to crossed_cells since they need to be collapsed anyway
                    moist::operation::edge_split_1to3(*this, cell, crossed_edges[0], crossed_edges[1], *plane);
                    break;
            default:
                    break;
            }

            // TODO: Flushing this often is not ideal... better mark as deleted!
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
            if (geogram::mesh_cell_volume(*this, cell) == 0 || moist::geometry::has_duplicate_vertex(cell, *this))
            {
                _deleted_tets.insert(cell);
            }
        }

        this->FlushTetrahedra();

    #ifndef NDEBUG
        // check if the edge actually exists...
        bool exists = false;
        for (const auto debug_cell : this->cells)
        {
            for (g_index debug_edge = 0; debug_edge < this->cells.nb_edges(debug_cell); debug_edge++)
            {
                const g_index debug_v0 = this->cells.edge_vertex(debug_cell, debug_edge, 0);
                const g_index debug_v1 = this->cells.edge_vertex(debug_cell, debug_edge, 1);
                const vec3 debug_p0 = this->vertices.point(debug_v0);
                const vec3 debug_p1 = this->vertices.point(debug_v1);
                if (debug_p0 == p0 && debug_p1 == p1 || debug_p0 == p1 && debug_p1 == p0)
                {
                    exists = true;
                }
            }
        }

        if (!exists)
        {
            OOC_DEBUG("Edge insertion failed!");
            moist::utils::dump_mesh(*this, std::format("test/target/_after_edge_insert_fail{}.geogram", i));
        }
    #endif

    }

#ifndef NDEBUG
    moist::utils::dump_mesh(*this, "test/target/__after_edge_insert_final.geogram");
#endif
}

void moist::MeshSlice::CreateTetrahedra()
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

void moist::MeshSlice::FlushTetrahedra()
{
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, false);
    _deleted_tets.clear();
}
