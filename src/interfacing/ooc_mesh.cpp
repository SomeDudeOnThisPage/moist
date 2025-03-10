#include "ooc_mesh.hpp"
#include "predicates.inl"

static inline bool vec_eq(const geogram::vec3 &v0, const geogram::vec3 &v1)
{
    return v0.x == v1.x && v0.y == v1.y && v0.z == v1.z;
}

ooc::OOCMesh::OOCMesh(std::string identifier, geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision), _identifier(identifier)
{

}

void ooc::OOCMesh::InsertInterface(ooc::Interface& interface)
{
    if (!interface.HasMeshConstraints(_identifier))
    {
        OOC_ERROR("attempted inserting interface which does not contain constrains of mesh '" << _identifier << "'");
    }

    _deleted_tets = std::unordered_set<geogram::index_t>();

    // std::vector<bool> mutex; // bitfield
    // mutex.resize(this->cells.nb()); // Need to re-allocate after each for-loop to accomodate new tets, similar to fTetWild parallelism
    // geogram::parallel_for /* ... */ // Lock mutex if a tet is deleted, disallow other operations on that tet in same loop, iterate until no points are left...
    const auto triangulation = interface.Triangulation();
    for (auto v_id : triangulation->vertices)
    {
        // Skip vertices that are part of this mesh.
        if (interface.GetMappedVertex(_identifier, v_id))
        {
            continue;
        }

        /* const */ geogram::vec3 point = triangulation->vertices.point(v_id);
        point.z = -1.0; // TODO: This needs to be set in the reprojection of the interface...

        // OOC_DEBUG("Inserting Vertex " << point << ", deleted.size() = " << _deleted_tets.size());
        this->InsertVertex(point);
    }

    // Remove all "deleted" tetrahedra from the underlying data structure.
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, true);

#ifndef NDEBUG
    OOC_DEBUG("Validating point insertion...");
    for (const auto interface_vertex : triangulation->vertices)
    {
        auto interface_point = triangulation->vertices.point(interface_vertex);
        interface_point.z = -1.0;
        bool found = false;
        for (const auto mesh_vertex : this->vertices)
        {
            if (vec_eq(this->vertices.point(mesh_vertex), interface_point))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            OOC_ERROR("Missing point " << interface_point << " in mesh " << this->_identifier);
        }
    }
    OOC_DEBUG("Done validating point insertion...");
#endif // NDEBUG

#ifndef NDEBUG
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    ooc::export_delaunay("mesh_insert.mesh", *this);
#endif // NDEBUG
}

void ooc::OOCMesh::InsertVertex(geogram::vec3 &point)
{
    // TODO: We don't need to iterate all cells here, only the subset on the interface...
    //       AABB tree is possible but probably overkill and not easily parallelizable, as it needs to be updated accordingly after each operation.
    //       This is also easily parallelizable as we only need to "find" the "right" tet.
    for (auto c_id : this->cells)
    {
        if (!_deleted_tets.contains(c_id) && ooc::predicates::point_in_tet(*this, c_id, point))
        {
            auto v0_id = this->cells.vertex(c_id, 0);
            auto v1_id = this->cells.vertex(c_id, 1);
            auto v2_id = this->cells.vertex(c_id, 2);
            auto v3_id = this->cells.vertex(c_id, 3);
            auto p0 = this->vertices.point(v0_id);
            auto p1 = this->vertices.point(v1_id);
            auto p2 = this->vertices.point(v2_id);
            auto p3 = this->vertices.point(v3_id);

            if (vec_eq(point, p0) || vec_eq(point, p1) || vec_eq(point, p2) || vec_eq(point, p3))
            {
                continue;
            }

            // One of these still generates tets with negative volume -> check after each insert to see which it is...
            geogram::index_t t0_id, t1_id, t2_id;
            const geogram::index_t centroid_index = this->vertices.create_vertex(point.data());
            if (p0.z > -1.0f)
            {
                t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
                t1_id = this->cells.create_tet(v0_id, v1_id, v3_id, centroid_index);
                t2_id = this->cells.create_tet(v0_id, v2_id, v3_id, centroid_index);
            }
            else if (p1.z > -1.0f)
            {
                t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
                t1_id = this->cells.create_tet(v0_id, v1_id, v3_id, centroid_index);
                t2_id = this->cells.create_tet(v1_id, v2_id, v3_id, centroid_index);
            }
            else if (p2.z > -1.0f)
            {
                t0_id = this->cells.create_tet(v0_id, v1_id, v2_id, centroid_index);
                t1_id = this->cells.create_tet(v1_id, v2_id, v3_id, centroid_index);
                t2_id = this->cells.create_tet(v0_id, v2_id, v3_id, centroid_index); // WHY DOESNT THIS WORKKKKK
                // t2_id = 1;
            }
            else if (p3.z > -1.0f)
            {
                t0_id = this->cells.create_tet(v3_id, v2_id, v1_id, centroid_index);
                t1_id = this->cells.create_tet(v3_id, v2_id, v0_id, centroid_index);
                t2_id = this->cells.create_tet(v3_id, v1_id, v0_id, centroid_index);
            }

            if (geogram::mesh_cell_volume(*this, t0_id) == 0.0f)
            {
                OOC_DEBUG("cell t0 " << t0_id << " has zero or volume. " << point << " -> " << p0 << "; " << p1 << "; " << p2 << "; " << p3);
                _deleted_tets.insert(t0_id);
            }

            if (geogram::mesh_cell_volume(*this, t1_id) <= 0.0f)
            {
                OOC_DEBUG("cell t1 " << t1_id << " has zero or negative volume.");
                _deleted_tets.insert(t1_id);
            }

            if (geogram::mesh_cell_volume(*this, t2_id) <= 0.0f)
            {
                OOC_DEBUG("cell t2 " << t2_id << " has zero or negative volume.");
                _deleted_tets.insert(t2_id);
            }

            _deleted_tets.insert(c_id);
            break;
        }
    }
}
