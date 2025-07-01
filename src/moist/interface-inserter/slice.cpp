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

moist::MeshSlice::MeshSlice(geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision)
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
    for (const g_index tet : tetrahedra)
    {
        _deleted_cells.insert(tet);
    }
}

moist::SteinerPoints moist::MeshSlice::InsertInterface(moist::Interface& interface, moist::metrics::Metrics_ptr metrics)
{
    this->_start_interface_cell = this->ReorderCells(*interface.Plane());

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

    {
        //auto timer = moist::Timer("MeshSlice::InsertTetQuality", metrics);
        this->InsertTetQuality(interface);
    }

    return steiner_points;
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

void moist::MeshSlice::InsertVertex(const geogram::vec3 &point, const moist::AxisAlignedInterfacePlane &plane)
{
    this->_created_cell_ids.clear();
    std::mutex m_deleted_tets;
    size_t edge_insert_cells = 0;
//#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
//    geogram::parallel_for(_start_interface_cell, this->cells.nb(), [this, &point, &plane, &m_deleted_tets, &edge_insert_cells](const g_index cell)
//    {
//#else
    const size_t nb = this->cells.nb();
    for (g_index cell = 0; cell < nb; cell++)
    {
//#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

        {
            std::lock_guard<std::mutex> lock(m_deleted_tets);
            if (_deleted_cells.contains(cell))
            {
                PARALLEL_CONTINUE;
            }
        }

        // 0, -0.583800, -1
        moist::predicates::PointInTet pit = moist::predicates::point_in_tet(*this, cell, point, true);
        if (pit == moist::predicates::PointInTet::FACET)
        {
            std::lock_guard<std::mutex> lock(m_deleted_tets);
            geogram::Mesh dbg(3);
            const g_index t = dbg.cells.create_tet(
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 0))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 1))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 2))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 3)))
            );
            moist::operation::vertex_insert_1to3(*this, cell, point, plane);

            //if (cell == 898)
            //{
            //    moist::utils::geo::save("debug/dbg-crossed-cell.msh", dbg);
            //    this->CreateTetrahedra();
            //    this->DebugMesh("debug/dbg-created_cells.msh", this->_created_cell_ids);
            //}
            //PARALLEL_BREAK;
            continue;
        }

        if (pit == moist::predicates::PointInTet::EDGE)
        {
            std::lock_guard<std::mutex> lock(m_deleted_tets);
            geogram::Mesh dbg(3);
            const g_index t = dbg.cells.create_tet(
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 0))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 1))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 2))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 3)))
            );
            //moist::utils::geo::save("debug/dbg-crossed-cell.msh", dbg);
            moist::operation::vertex_insert_1to2(*this, cell, point, plane);
            edge_insert_cells++;
            //PARALLEL_CONTINUE;
            continue;
        }

        if (pit == moist::predicates::PointInTet::VERTEX)
        {
            int isd = 0;
            // PARALLEL_BREAK;
            continue;
        }
    }
//ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
//    );
//#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    // if (edge_insert_cells > 0 && edge_insert_cells == 1)
    // {
        // OOC_WARNING("inserted vertex on edge - only one cell was split... possible boundary cell, otherwise invalid geometry...");
    // }

    this->CreateTetrahedra();
    if (!this->_created_cell_ids.empty())
    {
        //this->DebugMesh("debug/dbg-created_cells.msh", this->_created_cell_ids);
        this->_created_cell_ids.clear();
    }
}

void moist::MeshSlice::InsertInterfaceVertices(moist::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();
    int i = 0;
    for (const g_index v : this->vertices)
    {
        i++;
        if (moist::predicates::point_on_plane(this->vertices.point(v), *plane))
        {
            this->vertices.point(v).z = plane->extent;
        }
    }

    this->_start_interface_cell = this->ReorderCells(*plane);
    for (const g_index v : triangulation->vertices)
    {
        this->InsertVertex(triangulation->vertices.point(v), *interface.Plane());
    }

    this->FlushTetrahedra();
    OOC_DEBUG("mesh after vertex insertion: " << this->vertices.nb() << " vertices, " << this->cells.nb() << " cells");

    // geogram::mesh_repair(*this, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);
    for (const g_index cell : this->cells)
    {
        if (geogram::mesh_cell_volume(*this, cell) == 0.0 || moist::geometry::has_duplicate_vertex(cell, *this))
        {
            OOC_DEBUG("0 volume cell after insert...");
        }
    }

#ifndef NDEBUG
    for (const g_index iv : triangulation->vertices)
    {
        bool found = false;
        const vec3 ip = triangulation->vertices.point(iv);
        for (const g_index v : this->vertices)
        {
            if (this->vertices.point(v) == ip)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_WARNING("missing point [" << ip.x << ", " << ip.y << ", " << ip.z << "]");
        }
    }
    moist::utils::geo::save("after_insert.msh", *this);
#endif
}

bool moist::MeshSlice::CanMoveVertex(const g_index v, const vec3& p, const std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator>& cluster)
{
    std::vector<geogram::Sign> signs(0);
    for (size_t ci = 0; ci < this->_created_cell_ids.size(); ci++)
    {
        const g_index c = this->_created_cell_ids[ci];

        const vec3 p0 = this->vertices.point(this->cells.vertex(c, 0));
        const vec3 p1 = this->vertices.point(this->cells.vertex(c, 1));
        const vec3 p2 = this->vertices.point(this->cells.vertex(c, 2));
        const vec3 p3 = this->vertices.point(this->cells.vertex(c, 3));
        signs.push_back(geogram::geo_sgn(geogram::Geom::tetra_signed_volume(p0, p1, p2, p3)));
    }

    const auto original_pos = this->vertices.point(v);

    for (const g_index cv : cluster.at(original_pos))
    {
        this->vertices.point(cv) = p;
    }

    for (size_t ci = 0; ci < this->_created_cell_ids.size(); ci++)
    {
        if (_deleted_cells.contains(_created_cell_ids[ci])) continue;

        const g_index c = this->_created_cell_ids[ci];
        const vec3 p0 = this->vertices.point(this->cells.vertex(c, 0));
        const vec3 p1 = this->vertices.point(this->cells.vertex(c, 1));
        const vec3 p2 = this->vertices.point(this->cells.vertex(c, 2));
        const vec3 p3 = this->vertices.point(this->cells.vertex(c, 3));
        const auto volume = geogram::Geom::tetra_signed_volume(p0, p1, p2, p3);
        const auto sign = geogram::geo_sgn(std::round(volume * 1e12) / 1e12);
        // TODO: add _deleted_cells during decimation, so I can remove signs[ci] != geogram::Sign::ZERO
        if (signs[ci] != geogram::Sign::ZERO && sign != geogram::Sign::ZERO && sign != signs[ci])
        {
            for (const g_index cv : cluster.at(original_pos))
            {
                this->vertices.point(cv) = original_pos;
            }
            return false;
        }
    }

    for (const g_index cv : cluster.at(original_pos))
    {
        this->vertices.point(cv) = original_pos;
    }
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

    geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);

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
    double nearest_e0 = 9999.99;
    double nearest_e1 = 9999.99;
    g_index nearest_v_e0 = geogram::NO_VERTEX;
    g_index nearest_v_e1 = geogram::NO_VERTEX;
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

            if (geogram::distance(p, p_e0) < nearest_e0)
            {
                nearest_e0 = geogram::distance(p, p_e0);
                nearest_v_e0 = lv;
            }

            if (geogram::distance(p, p_e1) < nearest_e1)
            {
                nearest_e1 = geogram::distance(p, p_e1);
                nearest_v_e1 = lv;
            }
        }
    }

    e0 = nearest_e0 != -1 ? this->vertices.point(nearest_v_e0) : p_e0;
    e1 = nearest_e1 != -1 ? this->vertices.point(nearest_v_e1) : p_e1;

    cluster.emplace(e0, std::unordered_set<g_index>());
    cluster.emplace(e1, std::unordered_set<g_index>());

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
        }
    }

    for (const g_index c : this->_created_cell_ids)
    {
        if (geogram::mesh_cell_volume(*this, c) == 0.0 || moist::geometry::has_duplicate_vertex(c, *this))
        {
            _deleted_cells.insert(c);
        }
    }
}

void moist::MeshSlice::InsertInterfaceEdges(moist::Interface& interface, moist::SteinerPoints& steiner_points)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();
    const auto edges = moist::geometry::collect_edges(*triangulation);

    uint32_t i = 0;
    uint32_t highest = 0;

    for (const auto edge : edges)
    {
        i++;
        this->_created_cell_ids.clear();

        vec3 p0 = triangulation->vertices.point(edge.v0);
        vec3 p1 = triangulation->vertices.point(edge.v1);
        p0.z = plane->extent;
        p1.z = plane->extent;

        // find all tetrahedra that lie on the line between the two points
        const auto cells = this->cells.nb();
        std::vector<CreatedTetrahedon> crossed_cells;
        std::unordered_set<g_index> edge_vertices(0);

        size_t nb_created_cells = 0;
        size_t nb_crossed_cells = 0;
        geogram::Mesh dbg_crossed_cells(3);
        geogram::Mesh dbg_created_cells(3);

        std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices(0);
        for (g_index cell = _start_interface_cell; cell < cells; cell++)
        {
            nb_created_cells = 0;
            if (_deleted_cells.contains(cell))
            {
                continue;
            }

            if (geogram::mesh_cell_volume(*this, cell) == 0.0)
            {
                _deleted_cells.insert(cell);
                continue;
            }

            // insert vertices where the line crosses edges of the tetrahedra
            std::vector<CrossedEdge> crossed_edges;
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
                const g_index v = created_vertices.contains(new_vertex)
                    ? created_vertices.at(new_vertex)
                    : this->vertices.create_vertex(new_vertex.data());

                created_vertices.emplace(new_vertex, v);
                crossed_edges.push_back({ v0, v1, v });
                edge_vertices.insert(v);

#ifndef NDEBUG
                dbg_crossed_cells.cells.create_tet(
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 0))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 1))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 2))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 3)))
                );
#endif // NDEBUG

            }

            if (crossed_edges.size() > 0)
            {
                nb_crossed_cells++;
            }

            switch (crossed_edges.size())
            {
                case 1:
                    moist::operation::edge_split_1to2(*this, cell, crossed_edges[0], *plane);
                    nb_created_cells += 2;
                    break;
                case 2:
                    if (crossed_edges[0].p == crossed_edges[1].p)
                    {
                        OOC_WARNING("invalid edge split configuration...");
                    }
                    else
                    {
                        moist::operation::edge_split_1to3(*this, cell, crossed_edges[0], crossed_edges[1], *plane);
                        nb_created_cells += 3;
                    }
                    break;
            default:
                    break;
            }

            this->CreateTetrahedra();

            if (nb_created_cells > highest)
            {
                highest = nb_created_cells;
                OOC_DEBUG("new highest: " << i);
            }
        }

        // if (_is) this->DebugMesh("dbg-created-cells.msh", this->_created_cell_ids);
        const size_t num_steiner_points = steiner_points.size();
        this->DecimateCreatedTetrahedra(p0, p1, edge_vertices, steiner_points);
        // if (_is) this->DebugMesh("dbg-after-decimation.msh", this->_created_cell_ids);

        this->_created_cell_ids.clear();
    }

    this->FlushTetrahedra(true);
}

g_index moist::MeshSlice::ReorderCells(const moist::AxisAlignedInterfacePlane &plane)
{
    for (const g_index c : this->cells)
    {
        if (moist::predicates::cell_on_plane(c, *this, plane))
        {
            this->CreateTetrahedra({
                this->cells.vertex(c, 0),
                this->cells.vertex(c, 1),
                this->cells.vertex(c, 2),
                this->cells.vertex(c, 3)
            });
            this->DeleteTetrahedra(c);
        }
    }

    this->FlushTetrahedra();
    return this->CreateTetrahedra();
}

g_index moist::MeshSlice::CreateTetrahedra()
{
    g_index first = geogram::NO_CELL;
    for (const auto tet : this->_created_cells)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
        if (first == geogram::NO_CELL)
        {
            first = t;
        }

        _created_cell_ids.push_back(t);

    #ifndef NDEBUG
        if (geogram::mesh_cell_volume(*this, t) <= 0.0)
        {
            const auto p0 = this->vertices.point(tet.v0);
            const auto p1 = this->vertices.point(tet.v1);
            const auto p2 = this->vertices.point(tet.v2);
            const auto p3 = this->vertices.point(tet.v3);
            OOC_WARNING("cell " << t << " has zero volume");
            geogram::Mesh dbg(3);
            dbg.cells.create_tet(
                dbg.vertices.create_vertex(p0.data()),
                dbg.vertices.create_vertex(p1.data()),
                dbg.vertices.create_vertex(p2.data()),
                dbg.vertices.create_vertex(p3.data())
            );
            moist::utils::geo::save(std::format("debug/dbg-zerotet-{}.msh", t), dbg);
            //_deleted_cells.insert(t);
        }
     #endif // NDEBUG
    }

    this->_created_cells.clear();

    return first;
}

void moist::MeshSlice::FlushTetrahedra(bool delete_zero_volume)
{
    geogram::vector<g_index> flushed_elements(this->cells.nb());

    if (delete_zero_volume)
    {
        for (const g_index c : this->cells)
        {
            if (geogram::mesh_cell_volume(*this, c) == 0.0)
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

#ifndef NDEBUG
void moist::MeshSlice::DebugMesh(std::string file, std::vector<g_index>& tetrahedra)
{
    geogram::Mesh dbg(3);
    geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
    geogram::Attribute<int> v_discard_dbg(dbg.vertices.attributes(), ATTRIBUTE_DISCARD);

    for (const g_index c : tetrahedra)
    {
        if (!_deleted_cells.contains(c) && c < this->cells.nb())
        {
            g_index vertices[4];
            for (l_index lv = 0; lv < 4; lv++)
            {
                const g_index v = this->cells.vertex(c, lv);
                vertices[lv] = dbg.vertices.create_vertex(this->vertices.point(v));
                // if (v_discard[v])
                // {
                //     dbg.vertices.point(vertices[lv]).z -= 0.5;
                // }
                v_discard_dbg[vertices[lv]] = v_discard[v];
            }
            dbg.cells.create_tet(vertices[0], vertices[1], vertices[2], vertices[3]);
        }
    }

    //geogram::mesh_repair(dbg, geogram::MeshRepairMode::MESH_REPAIR_DEFAULT);
    moist::utils::geo::save(file, dbg);
}
#endif
