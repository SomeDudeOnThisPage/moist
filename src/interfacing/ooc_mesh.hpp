#ifndef __OOC_MESH_HPP
#define __OOC_MESH_HPP

#include <unordered_set>

#include <geogram/mesh/mesh.h>

#include "../core.hpp"
#include "interface.hpp"

namespace incremental_meshing
{
    typedef struct
    {
        geogram::index_t c;
        geogram::index_t le_v0;
        geogram::index_t le_v1;
        geogram::vec3 point;
    } EdgeSplit1_2;

    typedef struct Edge2d
    {
        geogram::index_t le_v0;
        geogram::index_t le_v1;
        geogram::vec3 intersection;

    #ifndef NDEBUG
        geogram::vec3 v0;
        geogram::vec3 v1;
    #endif // NDEBUG

        bool operator==(Edge2d& other) {
            return le_v0 == other.le_v0 && le_v1 == other.le_v1 || le_v0 == other.le_v1 && le_v1 == other.le_v0;
        }

        bool operator==(const Edge2d& other) const {
            return le_v0 == other.le_v0 && le_v1 == other.le_v1 || le_v0 == other.le_v1 && le_v1 == other.le_v0;
        }
    } Edge2d_t;

    class SubMesh : public geogram::Mesh
    {
    public:
        SubMesh(std::string identifier, geogram::index_t dimension = 3, bool single_precision = false);
        void InsertInterface(incremental_meshing::Interface& interface);

    private:
        std::string _identifier;
        std::shared_ptr<incremental_meshing::Interface> _interface;
        std::unordered_set<geogram::index_t> _deleted_tets;

        void InsertInterfaceVertices(incremental_meshing::Interface& interface);
        void InsertInterfaceEdges(const incremental_meshing::Interface& interface);
        void DecimateNonInterfaceEdges(const incremental_meshing::Interface& interface);

        void InsertVertex(const geogram::vec3& point, const incremental_meshing::Interface& interface);

        void Insert1To2(const EdgeSplit1_2& split, const incremental_meshing::AxisAlignedInterfacePlane& plane);
        void Insert1To3(const geogram::index_t c_id, const geogram::vec3& p0, const incremental_meshing::AxisAlignedInterfacePlane& plane);
        void Insert2To3(const geogram::index_t c_id, const geogram::vec3& p0, const geogram::vec3& p1, const incremental_meshing::AxisAlignedInterfacePlane& plane);
    };
}

namespace std
{
    template <>
    struct hash<incremental_meshing::Edge2d>
    {
        size_t operator()(const incremental_meshing::Edge2d& e) const
        {
            size_t h1 = std::hash<double>{}(e.le_v0) ^ (std::hash<double>{}(e.le_v0) << 1);
            size_t h2 = std::hash<double>{}(e.le_v1) ^ (std::hash<double>{}(e.le_v1) << 1);
            return h1 ^ (h2 << 1);
        }
    };
}

#endif // __OOC_MESH_HPP
