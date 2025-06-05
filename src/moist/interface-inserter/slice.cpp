#include "slice.hpp"

#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>
#include <limits>
#include <mutex>

#include "moist/core/metrics.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/mesh_quality.inl"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

#include "local_operations.hpp"

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

void moist::MeshSlice::InsertInterface(moist::Interface& interface, moist::metrics::TimeMetrics_ptr metrics)
{
    moist::descriptor::LocalInterface descriptor;
    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceVertices", metrics);
        this->InsertInterfaceVertices(interface);
    }

    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceEdges", metrics);
        this->InsertInterfaceEdges(interface);
    }

    {
        auto timer = moist::Timer("MeshSlice::InsertTetQuality", metrics);
        this->InsertTetQuality(interface);
    }

#ifndef NDEBUG
    // this->Validate(interface);
#endif
}

void moist::MeshSlice::InsertTetQuality(moist::Interface& interface)
{
    for (const g_index cell : this->cells)
    {
        if (!moist::predicates::cell_on_plane(cell, *this, *interface.Plane()))
        {
            continue;
        }

        // find corresponding interface facet to attach quality to
        // TODO: make a reducing list like "unmatched_facets" so we don't need to iterate all for each cell?
        // TODO: or simply parallelize this...
        for (const g_index facet : interface.Triangulation()->facets)
        {
            if (!moist::predicates::facet_matches_cell(cell, facet, *this, *interface.Triangulation()))
            {
                continue;
            }
        }
    }
}


void moist::MeshSlice::InsertInterfaceVertices(moist::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    for (const g_index v : triangulation->vertices)
    {
        const vec3 p = triangulation->vertices.point(v);
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
    // moist::utils::save(*this, "after_point_insertion.geogram");
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

        vec3 p0 = triangulation->vertices.point(edge.v0);
        vec3 p1 = triangulation->vertices.point(edge.v1);
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
            if (moist::geometry::has_duplicate_vertex(cell, *this))
            {
                _deleted_tets.insert(cell);
            }
        }

    /*#ifndef NDEBUG
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
            OOC_DEBUG("Edge insertion from " << p0 << " to " << p1 << " failed!");
            for (const g_index vv : this->vertices)
            {
                if (this->vertices.point(vv) == p0 || this->vertices.point(vv) == p1)
                {
                    LOCK_ATTRIBUTES;
                    geogram::Attribute<int> v_debug(this->vertices.attributes(), ATTRIBUTE_DEBUG_V);
                    v_debug[vv] = 1;
                }
            }

            // moist::utils::dump_mesh(*this, std::format("test/target/_after_edge_insert_fail{}.geogram", i));
        }
    #endif*/

    }

    this->FlushTetrahedra();
}

void moist::MeshSlice::CreateTetrahedra()
{
    for (const auto tet : this->_created_tets)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
        _created_tets_idx.push_back(t);

        // TODO [Bugs]: Figure out why this happens. Also, figure out why this leads to a validatably correct output...
        const auto volume = geogram::mesh_cell_volume(*this, t);
        if (volume <= 0.00000)
        {
        #ifndef NDEBUG
            const auto p0 = this->vertices.point(tet.v0);
            const auto p1 = this->vertices.point(tet.v1);
            const auto p2 = this->vertices.point(tet.v2);
            const auto p3 = this->vertices.point(tet.v3);
            OOC_WARNING("cell t0 " << t << " has zero volume: " << volume);
        #endif // NDEBUG
            _deleted_tets.insert(t);
        }
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

// TODO [Testing]: Move this code into a gtest module...
void moist::MeshSlice::Validate(moist::Interface& interface)
{
    // validate vertices mesh -> interface
    for (const g_index v : this->vertices)
    {
        const vec3 p = this->vertices.point(v);
        if (moist::predicates::point_on_plane(p, *interface.Plane()))
        {
            bool found = false;
            for (const g_index v_i : interface.Triangulation()->vertices)
            {
                const vec3 p_i = interface.Triangulation()->vertices.point(v_i);
                if (p == p_i)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                OOC_DEBUG("validation: invalid point " << p << " does not exist in the interface");
            }
        }
    }

    // validate vertices interface -> mesh
    for (const g_index v_i : interface.Triangulation()->vertices)
    {
        const vec3 p_i = interface.Triangulation()->vertices.point(v_i);

        bool found = false;
        for (const g_index v : this->vertices)
        {
            const vec3 p = this->vertices.point(v);
            if (p == p_i)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            //OOC_DEBUG("validation: missing point " << p_i << " does not exist in the mesh");
        }
    }

    const auto interface_edges = moist::geometry::collect_edges(*interface.Triangulation());
    const auto mesh_interface_edges = moist::geometry::collect_edges(*interface.Triangulation(), *interface.Plane());

    for (const auto interface_edge : interface_edges)
    {
        const auto p0 = interface.Triangulation()->vertices.point(interface_edge.v0);
        const auto p1 = interface.Triangulation()->vertices.point(interface_edge.v1);

        bool found = false;
        for (const auto mesh_interface_edge : mesh_interface_edges)
        {
            const auto ep0 = this->vertices.point(mesh_interface_edge.v0);
            const auto ep1 = this->vertices.point(mesh_interface_edge.v1);

            if (ep0 == p0 && ep1 == p1 || ep0 == p1 && ep1 == p0)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_DEBUG("validation: missing edge (" << p0 << " -> " << p1 << ") in slice");
        }
    }

    for (const auto mesh_interface_edge : mesh_interface_edges)
    {
        const auto p0 = this->vertices.point(mesh_interface_edge.v0);
        const auto p1 = this->vertices.point(mesh_interface_edge.v1);

        bool found = false;
        for (const auto interface_edge : interface_edges)
        {


            const auto ep0 = interface.Triangulation()->vertices.point(interface_edge.v0);
            const auto ep1 = interface.Triangulation()->vertices.point(interface_edge.v1);

            if (ep0 == p0 && ep1 == p1 || ep0 == p1 && ep1 == p0)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_DEBUG("validation: additional edge (" << p0 << " -> " << p1 << ") in slice");
        }
    }
}
