#include <unordered_set>
#include <optional>
#include <geogram/mesh/mesh_repair.h>

#include "ooc_mesh.hpp"
#include "predicates.inl"

incremental_meshing::SubMesh::SubMesh(std::string identifier, geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision), _identifier(identifier)
{

}

void incremental_meshing::SubMesh::InsertInterface(incremental_meshing::Interface& interface)
{
    if (!interface.HasMeshConstraints(_identifier))
    {
        OOC_ERROR("attempted inserting interface which does not contain constrains of mesh '" << _identifier << "'");
    }

    this->InsertInterfaceVertices(interface);
    this->InsertInterfaceEdges(interface);
    this->DecimateNonInterfaceEdges(interface);
}

void incremental_meshing::SubMesh::InsertInterfaceVertices(incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface vertices");
    _deleted_tets = std::unordered_set<geogram::index_t>();

    // std::vector<bool> mutex; // bitfield
    // mutex.resize(this->cells.nb()); // Need to re-allocate after each for-loop to accomodate new tets, similar to fTetWild parallelism
    // geogram::parallel_for /* ... */ // Lock mutex if a tet is deleted, disallow other operations on that tet in same loop, iterate until no points are left...
    const auto triangulation = interface.Triangulation();
    for (auto v_id : triangulation->vertices)
    {
        // Skip vertices that are part of this mesh.
        //if (interface.GetMappedVertex(_identifier, v_id))
        //{
        //    continue;
        //}

        /* const */ geogram::vec3 point = triangulation->vertices.point(v_id);
        point.z = -1.0; // TODO: This needs to be set in the reprojection of the interface...

        this->InsertVertex(point, interface);
    }

    // Remove all "deleted" tetrahedra from the underlying data structure.
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, true);
    _deleted_tets.clear();

    TIMER_END;

#ifndef NDEBUG
    OOC_DEBUG("Validating point insertion...");
    for (const auto interface_vertex : triangulation->vertices)
    {
        auto interface_point = triangulation->vertices.point(interface_vertex);
        bool found = false;
        for (const auto mesh_vertex : this->vertices)
        {
            if (incremental_meshing::predicates::vec_eq_2d(this->vertices.point(mesh_vertex), interface_point, interface.Plane()))
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
#endif // NDEBUG

#ifndef NDEBUG
    mesh_reorient(*this);
    incremental_meshing::export_delaunay(this->_identifier + "_after_vertex_insertion.mesh", *this);
#endif // NDEBUG
}

typedef struct
{
    geogram::index_t c_id;
    geogram::vec3 point;
} EdgeToSplit;

// for each edge in interface:
//  if edge exists in this mesh: continue
//  if edge does not exist, and does not cross any other edges: split 2 -> 4; continue
//  if edge does not exist, and crosses other edges:
//   collect all edges crossed, and where they are crossed
//   split each crossed edge by inserting a vertex at the crossing point (2 -> 4 split)
void incremental_meshing::SubMesh::InsertInterfaceEdges(const incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface edges into mesh '" + _identifier + "'");

    const auto triangulation = interface.Triangulation();

    for (const auto e_id : triangulation->edges)
    {
        const auto e_p0 = triangulation->vertices.point(triangulation->edges.vertex(e_id, 0));
        const auto e_p1 = triangulation->vertices.point(triangulation->edges.vertex(e_id, 1));
        OOC_DEBUG("e#" << e_id << ": " << e_p0 << " -> " << e_p1);

        std::vector<EdgeToSplit> edges_to_split;

        const auto cells = this->cells.nb();
        for (geogram::index_t c_id = 0; c_id < cells; c_id++)
        {
            if (_deleted_tets.contains(c_id))
            {
                continue;
            }

            // test
            std::vector<incremental_meshing::EdgeSplit1_2> crossed_edges_local;
            for (geogram::index_t le_id = 0; le_id < this->cells.nb_edges(c_id); le_id++)
            {
                const auto le_v0 = this->cells.edge_vertex(c_id, le_id, 0);
                const auto le_v1 = this->cells.edge_vertex(c_id, le_id, 1);

                const auto le_p0 = this->vertices.point(le_v0);
                const auto le_p1 = this->vertices.point(le_v1);

                if (!incremental_meshing::predicates::point_on_plane(le_p0, interface.Plane()) || !incremental_meshing::predicates::point_on_plane(le_p1, interface.Plane()))
                {
                    continue;
                }

                const auto opt_intersection = incremental_meshing::predicates::get_line_intersection(
                    geogram::vec2(e_p0.x, e_p0.y),
                    geogram::vec2(e_p1.x, e_p1.y),
                    geogram::vec2(le_p0.x, le_p0.y),
                    geogram::vec2(le_p1.x, le_p1.y)
                );

                if (opt_intersection.has_value())
                {
                    const geogram::vec3 intersection = geogram::vec3(opt_intersection.value().x, opt_intersection.value().y, interface.Plane().extent);
                    const incremental_meshing::EdgeSplit1_2 split
                    {
                        c_id,
                        le_v0,
                        le_v1,
                        intersection
                    };

                    crossed_edges_local.push_back(split);
                }
            }

            if (crossed_edges_local.size() == 2)
            {
                const auto p0_id = this->vertices.create_vertex(crossed_edges_local[0].point.data());
                const auto p1_id = this->vertices.create_vertex(crossed_edges_local[1].point.data());

                geogram::index_t a, b, c, d; // interface facet triangle vertices, d is the opposite vertex
                a = crossed_edges_local[0].le_v0;
                b = crossed_edges_local[0].le_v1;

                if (crossed_edges_local[0].le_v0 != crossed_edges_local[1].le_v0 && crossed_edges_local[0].le_v1 != crossed_edges_local[1].le_v0)
                {
                    c = crossed_edges_local[1].le_v0;
                }
                else
                {
                    c = crossed_edges_local[1].le_v1;
                }

                // find remaining vertex... can we make this faster / loopless?
                for (geogram::index_t lv_id = 0; lv_id < this->cells.nb_vertices(crossed_edges_local[0].c); lv_id++)
                {
                    const auto v_id = this->cells.vertex(crossed_edges_local[0].c, lv_id);
                    if (v_id != a && v_id != b && v_id != c)
                    {
                        d = v_id;
                        break;
                    }
                }

                // create three tetraheda, (a, p0, p1, d), (p0, p1, c, d), (p0, b, c, d)
                geogram::index_t t0_id, t1_id, t2_id;
                t0_id = this->cells.create_tet(p1_id, a, b, d);
                t1_id = this->cells.create_tet(p1_id, a, c, d);
                t2_id = this->cells.create_tet(p0_id, b, c, d);
                // this->cells.create_tet(a, b, c, d);

                _deleted_tets.insert(c_id);
            }
            else if (crossed_edges_local.size() == 1)
            {
                const auto p0_id = this->vertices.create_vertex(crossed_edges_local[0].point.data());
                geogram::index_t t0_id, t1_id;
                geogram::index_t a, b, c, d; // interface facet triangle vertices, d is the opposite vertex
                a = crossed_edges_local[0].le_v0;
                b = crossed_edges_local[0].le_v1;

                for (geogram::index_t lv_id = 0; lv_id < this->cells.nb_vertices(crossed_edges_local[0].c); lv_id++)
                {
                    const auto v_id = this->cells.vertex(crossed_edges_local[0].c, lv_id);
                    if (incremental_meshing::predicates::point_on_plane(this->vertices.point(v_id), interface.Plane()))
                    {
                        if (v_id == a || v_id == b)
                        {
                            continue;
                        }
                        c = v_id;
                        break;
                    }
                }

                for (geogram::index_t lv_id = 0; lv_id < this->cells.nb_vertices(crossed_edges_local[0].c); lv_id++)
                {
                    const auto v_id = this->cells.vertex(crossed_edges_local[0].c, lv_id);
                    if (v_id != a && v_id != b && v_id != c)
                    {
                        d = v_id;
                        break;
                    }
                }

                t0_id = this->cells.create_tet(p0_id, a, c, d);
                t1_id = this->cells.create_tet(p0_id, b, c, d);
                _deleted_tets.insert(c_id);
            }
        }
    }

    // Remove all "deleted" tetrahedra from the underlying data structure.
    OOC_DEBUG("deleting " << _deleted_tets.size() << " tets");
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, true);
    _deleted_tets.clear();
    TIMER_END;


#ifndef NDEBUG
    // mesh_reorient(*this);
    incremental_meshing::export_delaunay(this->_identifier + "_after_edge_insertion.mesh", *this);
#endif // NDEBUG
}

void incremental_meshing::SubMesh::DecimateNonInterfaceEdges(const incremental_meshing::Interface& interface)
{

}

void incremental_meshing::SubMesh::InsertVertex(const geogram::vec3& point, const incremental_meshing::Interface& interface)
{
    const auto plane = interface.Plane();
    // TODO: We don't need to iterate all cells here, only the subset on the interface...
    //       AABB tree is possible but probably overkill and not easily parallelizable, as it needs to be updated accordingly after each operation.
    //       This is also easily parallelizable as we only need to "find" the "right" tet.
    for (auto c_id : this->cells)
    {
        if (!_deleted_tets.contains(c_id) && incremental_meshing::predicates::point_in_tet(*this, c_id, point))
        {
            const auto p0 = this->vertices.point(this->cells.vertex(c_id, 0));
            const auto p1 = this->vertices.point(this->cells.vertex(c_id, 1));
            const auto p2 = this->vertices.point(this->cells.vertex(c_id, 2));
            const auto p3 = this->vertices.point(this->cells.vertex(c_id, 3));

            if (incremental_meshing::predicates::vec_eq_2d(point, p0, plane)
                || incremental_meshing::predicates::vec_eq_2d(point, p1, plane)
                || incremental_meshing::predicates::vec_eq_2d(point, p2, plane)
                || incremental_meshing::predicates::vec_eq_2d(point, p3, plane)
            )
            {
                continue;
            }

            this->Insert1To3(c_id, point, plane);

            _deleted_tets.insert(c_id);
            break;
        }
    }
}

void incremental_meshing::SubMesh::Insert1To2(const EdgeSplit1_2& split, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{

}

/**
 * Since this split is always on the interface-plane-facet, we don't need to care about splitting any incident faces, unlike the other splits.
 */
void incremental_meshing::SubMesh::Insert1To3(const geogram::index_t c_id, const geogram::vec3 &point, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
    const auto v0_id = this->cells.vertex(c_id, 0);
    const auto v1_id = this->cells.vertex(c_id, 1);
    const auto v2_id = this->cells.vertex(c_id, 2);
    const auto v3_id = this->cells.vertex(c_id, 3);
    const auto p0 = this->vertices.point(v0_id);
    const auto p1 = this->vertices.point(v1_id);
    const auto p2 = this->vertices.point(v2_id);
    const auto p3 = this->vertices.point(v3_id);

    // Some of these still generates tets with negative volume -> check after each insert to see which it is...
    // TODO: make this nice by indexing into an array of index_ts so we don't need 4 if statements
    geogram::index_t t0_id, t1_id, t2_id;
    const geogram::index_t centroid_index = this->vertices.create_vertex(point.data());
    if (!incremental_meshing::predicates::point_on_plane(p0, plane))
    {
        t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
        t1_id = this->cells.create_tet(v0_id, v1_id, v3_id, centroid_index);
        t2_id = this->cells.create_tet(v0_id, v2_id, v3_id, centroid_index);
    }
    else if (!incremental_meshing::predicates::point_on_plane(p1, plane))
    {
        t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
        t1_id = this->cells.create_tet(v0_id, v1_id, v3_id, centroid_index);
        t2_id = this->cells.create_tet(v1_id, v2_id, v3_id, centroid_index);
    }
    else if (!incremental_meshing::predicates::point_on_plane(p2, plane))
    {
        t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
        t1_id = this->cells.create_tet(v1_id, v2_id, v3_id, centroid_index);
        t2_id = this->cells.create_tet(v0_id, v2_id, v3_id, centroid_index);
    }
    else if (!incremental_meshing::predicates::point_on_plane(p3, plane))
    {
        t0_id = this->cells.create_tet(v3_id, v2_id, v1_id, centroid_index);
        t1_id = this->cells.create_tet(v3_id, v2_id, v0_id, centroid_index);
        t2_id = this->cells.create_tet(v3_id, v1_id, v0_id, centroid_index);
    }
    else
    {
        OOC_WARNING("attempted to insert into invalid tet " << p0 << " " << p1 << " " << p2 << " " << p3);
    }

#ifndef NDEBUG
    if (geogram::mesh_cell_volume(*this, t0_id) <= 0.0f)
    {
        OOC_WARNING("1->3 split: cell t0 " << t0_id << " has zero or negative volume: " << geogram::mesh_cell_volume(*this, t0_id));
        _deleted_tets.insert(t0_id);
    }

    if (geogram::mesh_cell_volume(*this, t1_id) <= 0.0f)
    {
        OOC_WARNING("1->3 split: cell t1 " << t1_id << " has zero or negative volume: " << geogram::mesh_cell_volume(*this, t1_id));
        _deleted_tets.insert(t1_id);
    }

    if (geogram::mesh_cell_volume(*this, t2_id) <= 0.0f)
    {
        OOC_WARNING("1->3 split: cell t2 " << t2_id << " has zero or negative volume: " << geogram::mesh_cell_volume(*this, t2_id));
        _deleted_tets.insert(t2_id);
    }
#endif // NDEBUG
}

void incremental_meshing::SubMesh::Insert2To3(const geogram::index_t c_id, const geogram::vec3 &p0, const geogram::vec3 &p1, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
}
