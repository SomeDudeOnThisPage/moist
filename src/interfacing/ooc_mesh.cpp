#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>

#include "ooc_mesh.hpp"
#include "predicates.inl"
#include "geometry.inl"

typedef struct EDGE_STRUCT
{
    g_index v0;
    g_index v1;

    bool operator==(const EDGE_STRUCT& other) const
    {
        return (v0 == other.v0 && v1 == other.v1) || (v0 == other.v1 && v1 == other.v0);
    }
} Edge;

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

incremental_meshing::SubMesh::SubMesh(std::string identifier, geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision), _identifier(identifier)
{
    _deleted_tets = std::unordered_set<geogram::index_t>();
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

    const auto triangulation = interface.Triangulation();
    for (auto v_id : triangulation->vertices)
    {
        this->InsertVertex(geogram::vec3(triangulation->vertices.point(v_id).x, triangulation->vertices.point(v_id).y, interface.Plane().extent), interface);
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
    geogram::mesh_repair(*this);
    incremental_meshing::export_delaunay(this->_identifier + "_after_vertex_insertion.mesh", *this, 3);
#endif // NDEBUG
}

void incremental_meshing::SubMesh::InsertInterfaceEdges(const incremental_meshing::Interface& interface)
{
    TIMER_START("insert interface edges into mesh '" + _identifier + "'");

    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();

    // sometimes edges are just not... added... by compute_borders, so do it manually here.
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

        const auto cells = this->cells.nb(); // parallel!
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        geogram::parallel_for(0, cells, [this, e_p0, e_p1, plane](const g_index c_id)
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
                const auto v0 = this->cells.edge_vertex(c_id, e_id, 0);
                const auto v1 = this->cells.edge_vertex(c_id, e_id, 1);
                if (!incremental_meshing::predicates::edge_on_plane(this->vertices.point(v0), this->vertices.point(v1), plane))
                {
                    continue;
                }

                const auto p0 = this->vertices.point(v0);
                const auto p1 = this->vertices.point(v1);

                // internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3
                if (!incremental_meshing::predicates::xy::check_lines_aabb(
                    reinterpret_cast<const geogram::vec2&>(p0),
                    reinterpret_cast<const geogram::vec2&>(p1),
                    reinterpret_cast<const geogram::vec2&>(e_p0),
                    reinterpret_cast<const geogram::vec2&>(e_p1)
                ))
                {
                    continue;
                }

                const auto intersection_opt = incremental_meshing::predicates::xy::get_line_intersection(
                    reinterpret_cast<const geogram::vec2&>(p0),
                    reinterpret_cast<const geogram::vec2&>(p1),
                    reinterpret_cast<const geogram::vec2&>(e_p0),
                    reinterpret_cast<const geogram::vec2&>(e_p1)
                );

                if (!intersection_opt.has_value())
                {
                    continue;
                }

                crossed_edges.push_back({ v0, v1, geogram::vec3(intersection_opt.value().x, intersection_opt.value().y, plane.extent) });
            }

            switch (crossed_edges.size())
            {
                case 1:
                    this->Split1_2(c_id, crossed_edges[0], plane);
                    break;
                case 2:
                    this->Split1_3(c_id, crossed_edges[0], crossed_edges[1], plane);
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
    incremental_meshing::export_delaunay(this->_identifier + "_after_edge_insertion.mesh", *this);
#endif // NDEBUG
}

void incremental_meshing::SubMesh::Split1_2(const g_index cell, const incremental_meshing::CrossedEdge& edge, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
    const g_index v_opposite = incremental_meshing::geometry::non_coplanar_opposite(cell, edge.e_v0, edge.e_v1, *this, plane);
    const g_index v_coplanar_opposite = incremental_meshing::geometry::opposite(cell, edge.e_v0, edge.e_v1, v_opposite, *this);
    const g_index p = this->vertices.create_vertex(edge.point.data());

    this->_created_tets.push_back({ edge.e_v0, v_coplanar_opposite, v_opposite, p });
    this->_created_tets.push_back({ edge.e_v1, v_coplanar_opposite, v_opposite, p });
    //this->cells.create_tet(edge.e_v0, v_coplanar_opposite, v_opposite, p);
    //this->cells.create_tet(edge.e_v1, v_coplanar_opposite, v_opposite, p);
    this->_deleted_tets.insert(cell);
}

void incremental_meshing::SubMesh::Split1_3(const g_index cell, const incremental_meshing::CrossedEdge& e0, const incremental_meshing::CrossedEdge& e1, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
    const g_index shared = (e0.e_v0 == e1.e_v0 || e0.e_v0 == e1.e_v1) ? e0.e_v0 : e0.e_v1;
    const g_index v_opposite = incremental_meshing::geometry::opposite(
        cell,
        e0.e_v0,
        e0.e_v1,
        (e1.e_v0 == shared) ? e1.e_v1 : e1.e_v0,
        *this
    );

    const g_index v_coplanar_opposite_p0 = (e1.e_v0 != shared) ? e1.e_v0 : e1.e_v1;
    const g_index v_coplanar_opposite_p1 = (e0.e_v0 != shared) ? e0.e_v0 : e0.e_v1;

    const g_index p0 = this->vertices.create_vertex(e0.point.data());
    const g_index p1 = this->vertices.create_vertex(e1.point.data());
    this->_created_tets.push_back({ v_opposite, shared, p0, p1 });
    this->_created_tets.push_back({ v_opposite, p0, p1, v_coplanar_opposite_p0 });
    this->_created_tets.push_back({ v_opposite, p0, v_coplanar_opposite_p0, v_coplanar_opposite_p1 });
    //this->cells.create_tet(v_opposite, p0, p1, v_coplanar_opposite_p0);
    //this->cells.create_tet(v_opposite, p0, v_coplanar_opposite_p0, v_coplanar_opposite_p1);
    //this->cells.create_tet(v_opposite, p0, p1, shared);

    this->_deleted_tets.insert(cell);
}

void incremental_meshing::SubMesh::FlushTetrahedra()
{
    OOC_DEBUG("flushing " << _deleted_tets.size() << " out of " << this->cells.nb() << " tets");
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        if (deleted_tet >= this->cells.nb())
        {
            OOC_WARNING("attempted to flush invalid tet #" << deleted_tet);
        }
        else
        {
            tets_to_delete[deleted_tet] = 1;
        }
    }
    this->cells.delete_elements(tets_to_delete, true);
    _deleted_tets.clear();
}

void incremental_meshing::SubMesh::CreateTetrahedra()
{
    //OOC_DEBUG("creating #" << this->_created_tets.size() << " tets");
    for (const auto tet : this->_created_tets)
    {
        this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
    }
    this->_created_tets.clear();
}

void incremental_meshing::SubMesh::DecimateNonInterfaceEdges(const incremental_meshing::Interface& interface)
{

}

void incremental_meshing::SubMesh::InsertVertex(const geogram::vec3& point, const incremental_meshing::Interface& interface)
{
    const auto plane = interface.Plane();
    for (const g_index c_id : this->cells)
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

// TODO: In case a vertex lies exactly on an edge, currently both tets would be split into three, creating two 0 volume tetrahedra!
void incremental_meshing::SubMesh::Insert2To3(const g_index c_id, const geogram::vec3 &p0, const geogram::vec3 &p1, const incremental_meshing::AxisAlignedInterfacePlane& plane)
{
}
