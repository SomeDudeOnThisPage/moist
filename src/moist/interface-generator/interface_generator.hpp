#ifndef MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_
#define MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_

#include <map>
#include <unordered_map>
#include <utility>

#include <geogram/mesh/mesh.h>

#include <CDT.hpp>

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
        InterfaceGenerator(const AxisAlignedPlane plane);

        /**
         * @brief Inserts a meshes' vertices, which are coplanar to the contained plane, into this interface.
         *
         * This method also inserts boundaries, meaning coplanar edges with only one adjacent facet on the interface-plane.
         *
         * @param mesh The mesh.
         */
        void AddConstraints(const GEO::Mesh& mesh);

        /**
         * @brief Performs the triangulation of all constraint vertices.
         */
        void Triangulate();

        /**
         * @brief Decimates n = ((total - constraints) / 2) edges belonging to the worst-quality triangles.
         */
        void Decimate();
    private:

#ifndef NDEBUG
        std::vector<GEO::vec3> _required_vertices;
        std::unordered_map<g_index, vec3> _inserted_points;
#endif // NDEBUG

        GEO::Mesh _constraints;
        size_t _unique_vertices;
        CDT::Triangulation<double> _cdt;
        std::vector<CDT::V2d<double>> _interface_vertices;
        CDT::EdgeVec _interface_edges;
    };
}

#endif // MOIST_INTERFACE_GENERATOR_INTERFACE_GENERATOR_HPP_
