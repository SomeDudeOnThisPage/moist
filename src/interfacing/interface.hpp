#ifndef __OOC_INTERFACE_HPP
#define __OOC_INTERFACE_HPP

#include <map>
#include <utility>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/delaunay/delaunay.h>

#include "../core.hpp"

namespace incremental_meshing {

    typedef struct
    {
        /** Defines the direction the mesh "grows" in. */
        Axis axis;
        /** Defines the position of the interface plane along the axis. */
        double extent;
        /** Defines an "envelope" in which vertices will be snapped onto the interface plane. */
        double epsilon;
    } AxisAlignedInterfacePlane;

    class Interface
    {
    public:
        Interface(const AxisAlignedInterfacePlane plane);

        void AddConstraints(std::string name, geogram::Mesh& mesh);
        void Triangulate();

        const bool HasMeshConstraints(std::string mesh) const;

        const geogram::index_t GetMappedVertex(std::string mesh, geogram::index_t v_id);
        const geogram::Mesh* Triangulation() const;

        AxisAlignedInterfacePlane Plane() const;
    private:
        geogram::Mesh _constraints;
        geogram::Mesh _triangulation;

        AxisAlignedInterfacePlane _plane;
        std::map<std::pair<double, double>, GEO::index_t> _indices;
        std::vector<std::pair<GEO::index_t, GEO::index_t>> _edges;
        std::map<std::string, std::map<geogram::index_t, geogram::index_t>> _interface_vertices;
    };

#ifndef NDEBUG
    void export_delaunay(const std::string filename, geogram::Mesh& mesh, int dimension = 3);
#endif // NDEBUG
}

#endif // __OOC_INTERFACE_HPP
