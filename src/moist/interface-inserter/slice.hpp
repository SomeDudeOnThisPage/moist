#ifndef MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_
#define MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_

#include <unordered_set>
#include <initializer_list>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/descriptor.hpp"
#include "moist/core/core_interface.hpp"

#ifdef OPTION_LOOKUP_GRID
    #include "lookup_grid.hpp"
#endif
#include "exact_mesh.hpp"

namespace moist
{
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
        bool operator()(const geo::vecng<3, double>& a, const geo::vecng<3, double>& b) const
        {
            return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
                std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON &&
                std::fabs(a.z - b.z) < moist::__DOUBLE_EPSILON;
        }
    };

    struct CrossedEdge
    {
        g_index e_v0;
        g_index e_v1;
        g_index p;

        g_index e_interface; // (custom) interface edge global index, used in decimation to find a point to decimate onto

        bool operator==(const CrossedEdge& other) const
        {
            return p == other.p && ((e_v0 == other.e_v0 && e_v1 == other.e_v1) || (e_v0 == other.e_v1 && e_v1 == other.e_v0));
        }
    };

    typedef struct
    {
        g_index v0; // alias v_line
        g_index v1;
        g_index v2;
        g_index v3;
    } CreatedTetrahedon;

    using SteinerPoints = std::unordered_set<vec3, Vec3HashOperator, Vec3EqualOperator>;

    class MeshSlice : public geo::Mesh
    {
    public:
        MeshSlice(const geo::index_t dimension = 3, const bool single_precision = false);
        void InsertEdges(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane);

        void InsertVertex(const geo::vec3& point, const moist::AxisAlignedPlane& plane);

        void CreateTetrahedra(const CreatedTetrahedon tet) { this->CreateTetrahedra({tet}); }
        void CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra);
        void DeleteTetrahedra(const g_index tet) { this->DeleteTetrahedra({tet}); };
        void DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra);
        void FlushTetrahedra(bool delete_zero_volume = false);

        void GetFixedGeometry(geo::Mesh& mesh);

        vec3& Point(const geo::index_t& v);
    private:
#ifdef OPTION_LOOKUP_GRID
        moist::LookupGrid _grid;
        moist::LookupGridExact _e_grid;
#endif
        moist::ExactMesh _em;

        std::vector<g_index> _created_cell_ids;
        std::vector<CreatedTetrahedon> _created_cells;
        std::unordered_set<geo::index_t> _deleted_cells;
        g_index _start_interface_cell;
        g_index _start_interface_vertex;

        void InsertVertices(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane);
        void InsertVertex2(const geo::vec3& p);
        void InsertVertexExact(const moist::ExactMesh::ExactPoint& p);

        void IsDeleted(const geo::index_t c);

        // global operations
        void InsertInterfaceVertices(moist::Interface& interface);
        void InsertEdges2(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane);

        g_index ReorderCells(const moist::AxisAlignedPlane& plane);
        g_index CreateTetrahedra();

        bool CanMoveVertex(const g_index v, const vec3& p, const std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator>& cluster);

        void DecimateCreatedTetrahedra(const vec3& e0, const vec3& e1, const std::unordered_set<g_index>& vertices, SteinerPoints& steiner_points);

#ifndef NDEBUG
        void DebugMesh(std::string file, std::vector<g_index>& tetrahedra, bool allow_deleted = false);
    #else
        void DebugMesh(std::string file, std::vector<g_index>& tetrahedra, bool allow_deleted = false) {};
    #endif
    };
}

namespace std
{
    template <>
    struct hash<moist::CrossedEdge>
    {
        std::size_t operator()(const moist::CrossedEdge& edge) const
        {
            return std::hash<geo::index_t>()(std::min(edge.e_v0, edge.e_v1)) ^ (std::hash<geo::index_t>()(std::max(edge.e_v0, edge.e_v1)) << 1);
        }
    };
}

#endif // MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_
