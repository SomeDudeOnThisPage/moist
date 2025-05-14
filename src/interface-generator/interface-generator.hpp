#ifndef INTERFACE_GENERATOR_HPP
#define INTERFACE_GENERATOR_HPP

#include <map>
#include <utility>

#include "core-interface.hpp"

namespace moist
{
    class InterfaceGenerator : public Interface
    {
    public:
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
        std::vector<vec3> _required_vertices;
#endif // NDEBUG

        geogram::Mesh _constraints;
    };
}

#endif // INTERFACE_GENERATOR_HPP
