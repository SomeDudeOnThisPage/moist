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
        geogram::Mesh _constraints;
        std::map<std::pair<double, double>, g_index> _indices;
    };
}

#endif // INTERFACE_GENERATOR_HPP
