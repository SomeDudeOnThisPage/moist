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
        /**
         * @brief Constructs a new Interface object from an existing triangulation.
         *
         * The given interface-mesh must lie on the xy-plane, as any third dimension is ignored.
         *
         * @param mesh Filesystem path to the (two-dimensional) mesh file (.geogram|.obj|.msh).
         * @param plane Plane to project the interface onto. Any local operations will use the plane as a reference.
         */
        Interface(const std::filesystem::path mesh, const AxisAlignedInterfacePlane plane);

        /**
         * @brief Destroys the Interface object
         */
        ~Interface() = default;

        /**
         * @brief Returns a pointer to the contained interface-mesh.
         *
         * @return std::shared_ptr<geogram::Mesh> The interface-triangulation.
         */
        std::shared_ptr<geogram::Mesh> Triangulation();

        /**
         * @brief Returns a pointer to the contained plane.
         *
         * @return std::shared_ptr<AxisAlignedInterfacePlane> The plane.
         */
        std::shared_ptr<AxisAlignedInterfacePlane> Plane();

    protected:
        // for subclassing in sub-projects...
        Interface() = default;
        std::shared_ptr<geogram::Mesh> _triangulation;
        std::shared_ptr<AxisAlignedInterfacePlane> _plane;
    };
}

#endif // OOC_CORE_INTERFACE_HPP
