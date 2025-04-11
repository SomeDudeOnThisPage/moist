#ifndef OOC_CORE_INTERFACE_HPP
#define OOC_CORE_INTERFACE_HPP

#include <map>
#include <utility>
#include <filesystem>

#include <geogram/mesh/mesh.h>

#include "core.hpp"

namespace incremental_meshing
{
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
        Interface(const std::filesystem::path mesh, const AxisAlignedInterfacePlane plane);
        ~Interface() = default;
        std::shared_ptr<geogram::Mesh> Triangulation();
        std::shared_ptr<AxisAlignedInterfacePlane> Plane();
    protected:
        Interface() = default;
        std::shared_ptr<geogram::Mesh> _triangulation;
        std::shared_ptr<AxisAlignedInterfacePlane> _plane;
    };
}

#endif // OOC_CORE_INTERFACE_HPP
