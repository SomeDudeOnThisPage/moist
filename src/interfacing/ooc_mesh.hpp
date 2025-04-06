#ifndef __OOC_MESH_HPP
#define __OOC_MESH_HPP

#include <unordered_set>

#include <geogram/mesh/mesh.h>

#include "../core.hpp"
#include "interface.hpp"

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

    class SubMesh : public geogram::Mesh
    {
    public:
        SubMesh(std::string identifier, geogram::index_t dimension = 3, bool single_precision = false);
        void InsertInterface(incremental_meshing::Interface& interface);

    private:
        std::string _identifier;
        std::vector<CreatedTetrahedon> _created_tets;
        std::unordered_set<geogram::index_t> _deleted_tets;

        void InsertInterfaceVertices(incremental_meshing::Interface& interface);
        void InsertInterfaceEdges(const incremental_meshing::Interface& interface);
        void DecimateNonInterfaceEdges(const incremental_meshing::Interface& interface);

        void InsertVertex(const geogram::vec3& point, const incremental_meshing::Interface& interface);

        void Split1_2(const g_index cell, const CrossedEdge& edge, const incremental_meshing::AxisAlignedInterfacePlane& plane);
        void Split1_3(const g_index cell, const incremental_meshing::CrossedEdge& e0, const incremental_meshing::CrossedEdge& e1, const incremental_meshing::AxisAlignedInterfacePlane& plane);

        void CreateTetrahedra();
        void FlushTetrahedra();

        void Insert1To3(const geogram::index_t c_id, const geogram::vec3& p0, const incremental_meshing::AxisAlignedInterfacePlane& plane);
        void Insert2To3(const geogram::index_t c_id, const geogram::vec3& p0, const geogram::vec3& p1, const incremental_meshing::AxisAlignedInterfacePlane& plane);
    };
}

#endif // __OOC_MESH_HPP
