#include "slice.hpp"

#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>
#include <limits>
#include <format>
#include <cmath>
#include <mutex>

#include "moist/core/metrics.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/mesh_quality.inl"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

#include "local_operations.hpp"

/* public */ moist::MeshSlice::MeshSlice(const geo::index_t dimension, const bool single_precision) : geo::Mesh(dimension, single_precision)
{
}

/* public */ void moist::MeshSlice::InsertEdges(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane)
{
    // Clamp vertices that are within eps of the plane, and set the v_interface attribute for later use.
    auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
    geo::parallel_for(0, this->vertices.nb(), [&](const geo::index_t v)
    {
        if (!moist::predicates::point_on_plane(this->vertices.point(v), plane))
        {
            { LOCK_ATTRIBUTES; v_interface[v] = false; }
            return;
        }

        { LOCK_ATTRIBUTES; v_interface[v] = true; }

        switch (plane.axis)
        {
            case moist::Axis::X:
                this->vertices.point(v).x = plane.extent;
                break;
            case moist::Axis::Y:
                this->vertices.point(v).y = plane.extent;
                break;
            case moist::Axis::Z:
                this->vertices.point(v).z = plane.extent;
                break;
            default:
                throw std::runtime_error("unknown axis definition");
        }
    });

    // Reorder all cells, so the cells touching the plane are packed tight at the end of the cell array.
    // Vertices don't really need this treatment, as most vertex accesses are indexed through cells anyway.
    std::mutex m;
    geo::parallel_for(0, this->cells.nb(), [&](const geo::index_t c)
    {
        if (!moist::predicates::cell_on_plane(c, *this, plane))
        {
            return;
        }

        std::lock_guard<std::mutex> lock(m);
        this->CreateTetrahedra({
            this->cells.vertex(c, 0),
            this->cells.vertex(c, 1),
            this->cells.vertex(c, 2),
            this->cells.vertex(c, 3)
        });

        this->DeleteTetrahedra(c);
    });

    this->FlushTetrahedra();
    this->_start_interface_cell = this->CreateTetrahedra();

    {
        //auto timer = moist::Timer("MeshSlice::InsertInterfaceVertices", metrics);
        this->InsertVertices(edge_mesh, plane);
    }

    {
        //auto timer = moist::Timer("MeshSlice::InsertInterfaceEdges", metrics);
        this->InsertEdges2(edge_mesh, plane);
    }
}


/* private */ void moist::MeshSlice::InsertVertices(const geo::Mesh &edge_mesh, const moist::AxisAlignedPlane &plane)
{
    for (const g_index v : edge_mesh.vertices)
    {
        this->InsertVertex2(edge_mesh.vertices.point(v));
    }

    this->FlushTetrahedra();
    OOC_DEBUG("mesh after vertex insertion: " << this->vertices.nb() << " vertices, " << this->cells.nb() << " cells");
}

/* public */ vec3& moist::MeshSlice::Point(const geo::index_t& v)
{
    return this->vertices.point(v);
}

/* private */ void moist::MeshSlice::InsertVertex2(const geo::vec3 &p)
{
    std::mutex m_deleted_tets;
    std::mutex m_1to2;

    this->_created_cell_ids.clear();
    auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
    auto c_deleted = geo::Attribute<bool>(this->cells.attributes(), "c_deleted");

    geo::parallel_for(this->_start_interface_cell, this->cells.nb(), [&](const g_index c)
    {
        g_index v_1to2 = geo::NO_VERTEX;

        if (c_deleted[c] || !moist::predicates::is_interface_cell(c, *this))
        {
            return;
        }

        const auto location = moist::predicates::point_in_tet(*this, c, p, true);
        if (location == moist::predicates::PointInTet::FACET)
        {
            const g_index v = this->vertices.create_vertex(p);
            /*LOCK_ATTRIBUTES;*/ v_interface[v] = true;

            moist::operation::InsertVertexOnCellBoundaryFacet(c, v, *this);
        }
        else if (location == moist::predicates::PointInTet::EDGE)
        {
            std::lock_guard<std::mutex> lock(m_1to2);
            if (v_1to2 == geo::NO_VERTEX)
            {
                v_1to2 = this->vertices.create_vertex(p);
                /*LOCK_ATTRIBUTES;*/ v_interface[v_1to2] = true;
            }

            moist::operation::InsertVertexOnCellBoundaryEdge(c, v_1to2, *this);
            // moist::operation::vertex_insert_1to2(*this, c, v_1to2, p);
        }
    });

    this->CreateTetrahedra();
}

void moist::MeshSlice::IsDeleted(const geo::index_t c)
{
}

void moist::MeshSlice::CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra)
{
    for (const CreatedTetrahedon tet : tetrahedra)
    {
        _created_cells.push_back(tet);
    }
}

void moist::MeshSlice::DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra)
{
    auto c_deleted = geo::Attribute<bool>(this->cells.attributes(), "c_deleted");
    for (const geo::index_t c : tetrahedra)
    {
        c_deleted[c] = true;
        _deleted_cells.insert(c);
    }
}

moist::SteinerPoints moist::MeshSlice::InsertInterface(moist::Interface& interface, moist::metrics::Metrics_ptr metrics)
{
    SteinerPoints steiner_points;
    moist::descriptor::LocalInterface descriptor;
    {
        //auto timer = moist::Timer("MeshSlice::InsertInterfaceVertices", metrics);
        this->InsertInterfaceVertices(interface);
    }

    {
        //auto timer = moist::Timer("MeshSlice::InsertInterfaceEdges", metrics);
        this->InsertInterfaceEdges(interface, steiner_points);
    }

    return steiner_points;
}

void moist::MeshSlice::GetFixedGeometry(geo::Mesh& mesh)
{
    // If a cell contains a fixed vertex v, add all edges of this cell adjacent to v to fixed geometry.
    std::unordered_map<g_index, g_index> vertices;
    std::unordered_set<moist::geometry::Edge, moist::geometry::EdgeHash> edges;

    geo::Attribute<bool> v_fixed(this->vertices.attributes(), "v_fixed");
    geo::Attribute<bool> v_interface(this->vertices.attributes(), "v_interface");

    for (g_index c = this->_start_interface_cell; c < this->cells.nb(); c++)
    {
        if (!moist::predicates::is_interface_cell(c, *this))
        {
            continue;
        }

        const auto c_vertices = moist::geometry::cell_vertices(c, *this);
        if (std::none_of(c_vertices.begin(), c_vertices.end(), [&](const g_index v) { return v_fixed[v]; }))
        {
            continue;
        }

        for (l_index le = 0; le < this->cells.nb_edges(c); le++)
        {
            const g_index v0 = this->cells.edge_vertex(c, le, 0);
            const g_index v1 = this->cells.edge_vertex(c, le, 1);

            if (!v_interface[v0] || !v_interface[v1])
            {
                continue;
            }

            const g_index nv0 = vertices.contains(v0) ? vertices[v0] : mesh.vertices.create_vertex(this->vertices.point(v0));
            const g_index nv1 = vertices.contains(v1) ? vertices[v1] : mesh.vertices.create_vertex(this->vertices.point(v1));
            vertices[v0] = nv0;
            vertices[v1] = nv1;

            const moist::geometry::Edge edge = {nv0, nv1};
            if (!edges.contains(edge))
            {
                edges.insert(edge);
                mesh.edges.create_edge(edge.v0, edge.v1);
            }
        }
    }

#ifndef NDEBUG
    moist::utils::geogram::save("additional_edges.obj", mesh, false);
#endif // NDEBUG
}

void moist::MeshSlice::InsertVertex(const geo::vec3 &point, const moist::AxisAlignedPlane &plane)
{
    this->_created_cell_ids.clear();

    std::mutex m_deleted_tets;
    std::mutex m_1to2;

    auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");

    geo::parallel_for(this->_start_interface_cell, this->cells.nb(), [&](const g_index cell)
    {
        g_index v_1to2 = geo::NO_VERTEX;

        {
            std::lock_guard<std::mutex> lock(m_deleted_tets);
            if (this->_deleted_cells.contains(cell))
            {
                return;
            }
        }

        const auto location = moist::predicates::point_in_tet(*this, cell, point, true);
        if (location == moist::predicates::PointInTet::FACET)
        {
            const g_index v = this->vertices.create_vertex(point);
            LOCK_ATTRIBUTES; v_interface[v] = true;

            moist::operation::vertex_insert_1to3(*this, cell, v, plane);
        }
        else if (location == moist::predicates::PointInTet::EDGE)
        {
            std::lock_guard<std::mutex> lock(m_1to2);
            if (v_1to2 == geo::NO_VERTEX)
            {
                v_1to2 = this->vertices.create_vertex(point);
                LOCK_ATTRIBUTES; v_interface[v_1to2] = true;
            }

            moist::operation::vertex_insert_1to2(*this, cell, v_1to2, plane);
        }
    });

    this->CreateTetrahedra();
}

void moist::MeshSlice::InsertInterfaceVertices(moist::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();


    for (const g_index v : triangulation->vertices)
    {
        this->InsertVertex2(triangulation->vertices.point(v));
        // this->InsertVertex(triangulation->vertices.point(v), *interface.Plane());
    }

    this->FlushTetrahedra();
    OOC_DEBUG("mesh after vertex insertion: " << this->vertices.nb() << " vertices, " << this->cells.nb() << " cells");

#ifndef NDEBUG
    moist::utils::geogram::save("after_insert.msh", *this);
#endif
}

/**
 * @brief Checks if vertex v can be moved onto point to without tetrahedral inversions of incident cells.
 *
 * @param v
 * @param to
 * @param cluster
 * @return true
 * @return false
 */
bool moist::MeshSlice::CanMoveVertex(const g_index v, const vec3& to, const std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator>& incident)
{
    auto set_vertices = [&](const vec3 from, const vec3 to)
    {
        for (const g_index cv : incident.at(from))
        {
            this->vertices.point(cv) = to;
        }
    };

    auto sign = [&](const g_index c)
    {
        return geo::geo_sgn(std::round(
            geo::Geom::tetra_signed_volume(
                this->vertices.point(this->cells.vertex(c, 0)),
                this->vertices.point(this->cells.vertex(c, 1)),
                this->vertices.point(this->cells.vertex(c, 2)),
                this->vertices.point(this->cells.vertex(c, 3))
            ) * 1e12
        ) / 1e12);
    };

    std::vector<geo::Sign> signs(0);
    const vec3 from = this->vertices.point(v);

    for (size_t ci = 0; ci < this->_created_cell_ids.size(); ci++)
    {
        if (_deleted_cells.contains(_created_cell_ids[ci])) continue;

        signs.push_back(sign(this->_created_cell_ids[ci]));
    }

    set_vertices(from, to);

    for (size_t ci = 0; ci < this->_created_cell_ids.size(); ci++)
    {
        if (_deleted_cells.contains(_created_cell_ids[ci])) continue;

        const geo::Sign s = sign(this->_created_cell_ids[ci]);
        if (/* signs[ci] != geogram::Sign::ZERO && */ s != geo::Sign::ZERO && s != signs[ci])
        {
            set_vertices(from, from);
            return false;
        }
    }

    set_vertices(from, from);
    return true;
}

void moist::MeshSlice::DecimateCreatedTetrahedra(const vec3& p_e0, const vec3& p_e1, const std::unordered_set<g_index>& vertices, SteinerPoints& steiner_points)
{
    if (vertices.empty() || this->_created_cell_ids.empty())
    {
        return;
    }

    std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator> cluster;
    std::unordered_map<g_index, std::unordered_set<g_index>> adjacencies;

    geo::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);

    // Preprocess: Initialize cluster and adjacency data.
    for (const g_index v : vertices)
    {
        const vec3 p = this->vertices.point(v);
        if (!cluster.contains(p))
        {
            cluster.emplace(p, std::unordered_set<g_index>());
        }
        cluster[p].insert(v);
    }

    // we can have cases where edges lead out of the boundary with envelope boundaries e.g. with tetwild...
    // in this case, e0 or e1 may not be part of created vertices...
    // need to take the closest vertex in that case and make it "fake" e0 or "fake" e1...
    double nearest_e0 = std::numeric_limits<double>::max();
    double nearest_e1 = std::numeric_limits<double>::max();;
    g_index nearest_v_e0 = geo::NO_VERTEX;
    g_index nearest_v_e1 = geo::NO_VERTEX;
    vec3 e0, e1;
    for (const g_index c : this->_created_cell_ids)
    {
        if (nearest_v_e0 == -1.0 && nearest_v_e1 == -1.0)
        {
            break;
        }
        for (const g_index lv : moist::geometry::cell_vertices(c, *this))
        {
            if (nearest_v_e0 == -1.0 && nearest_v_e1 == -1.0)
            {
                break;
            }
            const vec3 p = this->vertices.point(lv);

            if (!vertices.contains(lv) && p != p_e0 && p != p_e1)
            {
                continue;
            }

            if (p == p_e0)
            {
                nearest_e0 = -1.0;
                continue;
            }
            else if (p == p_e1)
            {
                nearest_e1 = -1.0;
                continue;
            }

            if (geo::distance(p, p_e0) < nearest_e0)
            {
                nearest_e0 = geo::distance(p, p_e0);
                nearest_v_e0 = lv;
            }

            if (geo::distance(p, p_e1) < nearest_e1)
            {
                nearest_e1 = geo::distance(p, p_e1);
                nearest_v_e1 = lv;
            }
        }
    }

    e0 = nearest_e0 != -1 ? this->vertices.point(nearest_v_e0) : p_e0;
    e1 = nearest_e1 != -1 ? this->vertices.point(nearest_v_e1) : p_e1;

    cluster.emplace(e0, std::unordered_set<g_index>());
    cluster.emplace(e1, std::unordered_set<g_index>());

    bool has_steiner_points = false;
    for (const g_index v : vertices)
    {
        bool must_be_steiner = true;
        const vec3 p = this->vertices.point(v);
        // Find neighbouring vertices in tetrahedra...
        std::unordered_set<vec3, Vec3HashOperator, Vec3EqualOperator> neighbours;
        for (const g_index c : this->_created_cell_ids)
        {
            // If vertex is part of cell, we can search this cell for other points on the line...
            if (!moist::geometry::point_of_cell(*this, c, p))
            {
                continue;
            }

            if (moist::geometry::has_duplicate_vertex(c, *this))
            {
                continue;
            }

            for (const g_index cv : moist::geometry::cell_vertices(c, *this))
            {
                const vec3 cvp = this->vertices.point(cv);
                if (cluster.contains(cvp) && cvp != p)
                {
                    neighbours.insert(cvp);
                }
            }
        }

        // Attempt to first cluster onto line endpoints...
        if (neighbours.contains(e0) || neighbours.contains(e1))
        {
            // Can contain both (last decimation)...
            if (neighbours.contains(e0) && this->CanMoveVertex(v, e0, cluster))
            {
                for (const g_index cluster_v : cluster[p])
                {
                    this->vertices.point(cluster_v) = e0;
                    cluster[e0].insert(cluster_v);
                }
                cluster.erase(p);
                must_be_steiner = false;
                continue;
            }

            if (neighbours.contains(e1) && this->CanMoveVertex(v, e1, cluster))
            {
                for (const g_index cluster_v : cluster[p])
                {
                    this->vertices.point(cluster_v) = e1;
                    cluster[e1].insert(cluster_v);
                }
                cluster.erase(p);
                must_be_steiner = false;
                continue;
            }
        }

        // Find a neighbour to cluster onto...
        bool moved = false;
        for (const vec3 neighbour : neighbours)
        {
            if (this->CanMoveVertex(v, neighbour, cluster))
            {
                for (const g_index cluster_v : cluster[p])
                {
                    this->vertices.point(cluster_v) = neighbour;
                    cluster[neighbour].insert(cluster_v);
                }
                cluster.erase(p);
                must_be_steiner = false;
                break;
            }
        }

        if (must_be_steiner)
        {
            OOC_DEBUG("inserting interface steiner point " << v << " at [" << p.x << ", " << p.y << ", " << p.z << "]");
            steiner_points.insert(p);

            // Mark the vertex for later extraction in re-hole-filling.
            geo::Attribute<bool> v_fixed(this->vertices.attributes(), "v_fixed");
            v_fixed[v] = true;
        }
    }

    // Cleanup any zero volume tetrahedra.
    for (const g_index c : this->_created_cell_ids)
    {
        if (geo::mesh_cell_volume(*this, c) == 0.0 || moist::geometry::has_duplicate_vertex(c, *this))
        {
            _deleted_cells.insert(c);
        }
    }

}

/* private */ void moist::MeshSlice::InsertEdges2(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane)
{
    const auto edges = moist::geometry::collect_edges(edge_mesh);

    for (const auto edge : edges)
    {
        this->_created_cell_ids.clear();

        const vec3 p0 = edge_mesh.vertices.point(edge.v0);
        const vec3 p1 = edge_mesh.vertices.point(edge.v1);

        // find all tetrahedra that lie on the line between the two points
        const auto nb_cells = this->cells.nb();
        std::vector<CreatedTetrahedon> crossed_cells;
        std::unordered_set<g_index> edge_vertices(0);
        std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices(0);

        for (g_index c = this->_start_interface_cell; c < nb_cells; c++)
        {
            if (this->_deleted_cells.contains(c))
            {
                continue;
            }

            if (geo::mesh_cell_volume(*this, c) == 0.0)
            {
                this->_deleted_cells.insert(c);
                continue;
            }

            // insert vertices where the line crosses edges of the tetrahedra
            std::vector<CrossedEdge> crossed_edges;
            for (l_index le = 0; le < this->cells.nb_edges(c); le++)
            {
                const g_index v0 = this->cells.edge_vertex(c, le, 0);
                const g_index v1 = this->cells.edge_vertex(c, le, 1);
                const vec3 cp0 = this->vertices.point(v0);
                const vec3 cp1 = this->vertices.point(v1);

                if (!moist::predicates::edge_on_plane(cp0, cp1, plane))
                {
                    continue;
                }

                // replace this at all with the aabb grid
                if (!moist::predicates::xy::check_lines_aabb(reinterpret_cast<const vec2&>(cp0), reinterpret_cast<const vec2&>(cp1), reinterpret_cast<const vec2&>(p0), reinterpret_cast<const vec2&>(p1)))
                {
                    continue;
                }

                // Internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3.
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

                const vec3 p = vec3(intersection_opt.value().x, intersection_opt.value().y, plane.extent);
                const g_index v = created_vertices.contains(p) ? created_vertices.at(p) : this->vertices.create_vertex(p.data());

                created_vertices.emplace(p, v);
                crossed_edges.push_back({ v0, v1, v });
                edge_vertices.insert(v);

                auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
                v_interface[v] = true;
            }

            switch (crossed_edges.size())
            {
                case 0:
                    break;
                case 1:
                    moist::operation::edge_split_1to2(*this, c, crossed_edges[0], plane);
                    break;
                case 2:
                    if (crossed_edges[0].p == crossed_edges[1].p)
                    {
                        OOC_WARNING("invalid edge split configuration - possible near-zero volume cell");
                    }
                    else
                    {
                        moist::operation::edge_split_1to3(*this, c, crossed_edges[0], crossed_edges[1], plane);
                    }
                    break;
                default:
                    OOC_WARNING("invalid edge split configuration - possible near-zero volume cell");
                    break;
            }

            this->CreateTetrahedra();
        }

        // if (_is) this->DebugMesh("dbg-created-cells.msh", this->_created_cell_ids);
        SteinerPoints sp;
        this->DecimateCreatedTetrahedra(p0, p1, edge_vertices, sp);
        // if (_is) this->DebugMesh("dbg-after-decimation.msh", this->_created_cell_ids);
    }

    this->FlushTetrahedra(true);
}

void moist::MeshSlice::InsertInterfaceEdges(moist::Interface &interface, moist::SteinerPoints &steiner_points)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();
    const auto edges = moist::geometry::collect_edges(*triangulation);

    for (const auto edge : edges)
    {
        _created_cell_ids.clear();

        const vec3 p0 = triangulation->vertices.point(edge.v0);
        const vec3 p1 = triangulation->vertices.point(edge.v1);

        // find all tetrahedra that lie on the line between the two points
        const auto nb_cells = this->cells.nb();
        std::vector<CreatedTetrahedon> crossed_cells;
        std::unordered_set<g_index> edge_vertices(0);
        std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices(0);

        for (g_index c = _start_interface_cell; c < nb_cells; c++)
        {
            if (_deleted_cells.contains(c))
            {
                continue;
            }

            if (geo::mesh_cell_volume(*this, c) == 0.0)
            {
                _deleted_cells.insert(c);
                continue;
            }

            // insert vertices where the line crosses edges of the tetrahedra
            std::vector<CrossedEdge> crossed_edges;
            for (l_index le = 0; le < this->cells.nb_edges(c); le++)
            {
                const g_index v0 = this->cells.edge_vertex(c, le, 0);
                const g_index v1 = this->cells.edge_vertex(c, le, 1);
                const vec3 cp0 = this->vertices.point(v0);
                const vec3 cp1 = this->vertices.point(v1);

                if (!moist::predicates::edge_on_plane(cp0, cp1, *plane))
                {
                    continue;
                }

                // replace this at all with the aabb grid
                if (!moist::predicates::xy::check_lines_aabb(reinterpret_cast<const vec2&>(cp0), reinterpret_cast<const vec2&>(cp1), reinterpret_cast<const vec2&>(p0), reinterpret_cast<const vec2&>(p1)))
                {
                    continue;
                }

                // Internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3.
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

                const vec3 p = vec3(intersection_opt.value().x, intersection_opt.value().y, plane->extent);
                const g_index v = created_vertices.contains(p) ? created_vertices.at(p) : this->vertices.create_vertex(p.data());

                created_vertices.emplace(p, v);
                crossed_edges.push_back({ v0, v1, v });
                edge_vertices.insert(v);

                auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
                v_interface[v] = true;
            }

            switch (crossed_edges.size())
            {
                case 1:
                    moist::operation::edge_split_1to2(*this, c, crossed_edges[0], *plane);
                    break;
                case 2:
                    if (crossed_edges[0].p == crossed_edges[1].p)
                    {
                        OOC_WARNING("invalid edge split configuration - possible near-zero volume cell");
                    }
                    else
                    {
                        moist::operation::edge_split_1to3(*this, c, crossed_edges[0], crossed_edges[1], *plane);
                    }
                    break;
            default:
                    break;
            }

            this->CreateTetrahedra();
        }

        // if (_is) this->DebugMesh("dbg-created-cells.msh", this->_created_cell_ids);
        const size_t num_steiner_points = steiner_points.size();
        this->DecimateCreatedTetrahedra(p0, p1, edge_vertices, steiner_points);
        // if (_is) this->DebugMesh("dbg-after-decimation.msh", this->_created_cell_ids);
    }

    this->FlushTetrahedra(true);
}

g_index moist::MeshSlice::ReorderCells(const moist::AxisAlignedPlane &plane)
{
    return 0;
}

g_index moist::MeshSlice::CreateTetrahedra()
{
    g_index first = geo::NO_CELL;
    for (const auto tet : this->_created_cells)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
        if (first == geo::NO_CELL)
        {
            first = t;
        }

        _created_cell_ids.push_back(t);

    #ifndef NDEBUG
        if (geo::mesh_cell_volume(*this, t) <= 0.0)
        {
            const auto p0 = this->vertices.point(tet.v0);
            const auto p1 = this->vertices.point(tet.v1);
            const auto p2 = this->vertices.point(tet.v2);
            const auto p3 = this->vertices.point(tet.v3);
            OOC_WARNING("cell " << t << " has zero volume");
        }
     #endif // NDEBUG
    }

    this->_created_cells.clear();

    return first;
}

void moist::MeshSlice::FlushTetrahedra(bool delete_zero_volume)
{
    geo::vector<g_index> flushed_elements(this->cells.nb());

    if (delete_zero_volume)
    {
        for (const g_index c : this->cells)
        {
            if (geo::mesh_cell_volume(*this, c) == 0.0)
            {
                flushed_elements[c] = true;
            }
        }
    }

    for (const g_index c : _deleted_cells)
    {
        flushed_elements[c] = true;
    }

    this->cells.delete_elements(flushed_elements, false);
    _deleted_cells.clear();
}

#ifndef NDEBUG
void moist::MeshSlice::DebugMesh(std::string file, std::vector<g_index>& tetrahedra)
{
    geo::Mesh dbg(3);
    geo::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
    geo::Attribute<int> v_discard_dbg(dbg.vertices.attributes(), ATTRIBUTE_DISCARD);

    for (const g_index c : tetrahedra)
    {
        if (!_deleted_cells.contains(c) && c < this->cells.nb())
        {
            g_index vertices[4];
            for (l_index lv = 0; lv < 4; lv++)
            {
                const g_index v = this->cells.vertex(c, lv);
                vertices[lv] = dbg.vertices.create_vertex(this->vertices.point(v));
                v_discard_dbg[vertices[lv]] = v_discard[v];
            }
            dbg.cells.create_tet(vertices[0], vertices[1], vertices[2], vertices[3]);
        }
    }

    moist::utils::geogram::save(file, dbg);
}
#endif
