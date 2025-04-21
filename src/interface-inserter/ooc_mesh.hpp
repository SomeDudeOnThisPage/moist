#ifndef __OOC_MESH_HPP
#define __OOC_MESH_HPP

#include <unordered_set>
#include <initializer_list>

#include <geogram/mesh/mesh.h>

#include "core.hpp"
#include "core-interface.hpp"

namespace incremental_meshing
{
    typedef struct CROSSED_EDGE_STRUCT
    {
        g_index e_v0;
        g_index e_v1;
        g_index p;

        bool operator==(const CROSSED_EDGE_STRUCT& other) const
        {
            return p == other.p && ((e_v0 == other.e_v0 && e_v1 == other.e_v1) || (e_v0 == other.e_v1 && e_v1 == other.e_v0));
        }
    } CrossedEdge;

    typedef struct
    {
        g_index v0;
        g_index v1;
        g_index v2;
        g_index v3;
    } CreatedTetrahedon;

    class MeshSlice : public geogram::Mesh
    {
    public:
        MeshSlice(geogram::index_t dimension = 3, bool single_precision = false);
        void InsertInterface(incremental_meshing::Interface& interface);

        void CreateTetrahedra(const CreatedTetrahedon tet) { this->CreateTetrahedra({tet}); }
        void CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra);
        void DeleteTetrahedra(const g_index tet) { this->DeleteTetrahedra({tet}); };
        void DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra);
    private:
        std::string _identifier;
        std::vector<CreatedTetrahedon> _created_tets;
        std::unordered_set<geogram::index_t> _deleted_tets;

        void InsertInterfaceVertices(incremental_meshing::Interface& interface);
        void InsertInterfaceEdges(incremental_meshing::Interface& interface);
        void DecimateNonInterfaceEdges(incremental_meshing::Interface& interface);

        void InsertVertex(const geogram::vec3& point, const incremental_meshing::AxisAlignedInterfacePlane& plane);

        void CreateTetrahedra();
        void FlushTetrahedra();
    };
}

#endif // __OOC_MESH_HPP
