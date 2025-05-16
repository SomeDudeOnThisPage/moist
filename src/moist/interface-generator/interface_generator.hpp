#ifndef MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_
#define MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_

#include <map>
#include <utility>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"

namespace moist
{
    class InterfaceGenerator : public Interface
    {
    public:
        /**
         * @brief Construct a new Interface Generator object.
         *
         * @param plane The plane to generate the interface from.
         */
        InterfaceGenerator(const AxisAlignedInterfacePlane plane);

        /**
         * @brief Inserts a meshes' vertices, which are coplanar to the contained plane, into this interface.
         *
         * This method also inserts boundaries, meaning coplanar edges with only one adjacent facet on the interface-plane.
         *
         * @param mesh The mesh.
         */
        void AddConstraints(const geogram::Mesh& mesh);

        /**
         * @brief Performs the triangulation of all constraint vertices.
         */
        void Triangulate();
    private:

#ifndef NDEBUG
        std::vector<geogram::vec3> _required_vertices;
#endif // NDEBUG

        geogram::Mesh _constraints;
    };
}

#endif // MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_
