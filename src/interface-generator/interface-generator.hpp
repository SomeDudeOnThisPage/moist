#ifndef INTERFACE_GENERATOR_HPP
#define INTERFACE_GENERATOR_HPP

#include <map>
#include <utility>

#include "core-interface.hpp"

namespace incremental_meshing
{
    class InterfaceGenerator : public Interface
    {
    public:
        InterfaceGenerator(const AxisAlignedInterfacePlane plane);
        void AddConstraints(geogram::Mesh& mesh);
        void Triangulate();
    private:
#ifndef NDEBUG
        std::vector<vec3> _required_vertices;
#endif // NDEBUG
        geogram::Mesh _constraints;
    };
}

#endif // INTERFACE_GENERATOR_HPP
