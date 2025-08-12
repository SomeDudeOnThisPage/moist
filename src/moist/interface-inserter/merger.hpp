#ifndef MOIST_INTERFACE_INSERTER_MERGER_HPP_
#define MOIST_INTERFACE_INSERTER_MERGER_HPP_

#include <unordered_map>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"

#include "exact_mesh.hpp"

namespace moist
{
    class Merger
    {
    public:
        Merger(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane);
        ~Merger() = default;
    private:

        enum class Slice
        {
            A, B
        };

        moist::ExactMesh _mesh_a;
        moist::ExactMesh _mesh_b;
        moist::ExactMesh _crown;
        moist::ExactMesh _crown_surface;

        std::unordered_map<geo::index_t, std::size_t> _vertices;

        void ConstructExactMesh(geo::Mesh& mesh, moist::ExactMesh& exact_mesh, const moist::AxisAlignedPlane& plane);
        void InsertAndMapPoints(const moist::ExactMesh& from, moist::ExactMesh& to);

        void Prune(moist::ExactMesh& mesh, const double min_volume = 1e-6);

        void HollowCrown();

        void InsertPoint(const std::size_t v, const moist::ExactMesh& from, moist::ExactMesh& to);
        void InsertBToA();
        void InsertCellBToA(const std::size_t c);
        void InsertCellIntoCell(const std::size_t ca, const std::size_t cb);
        void InsertPointIntoEdges(const moist::exact::Point& point, moist::ExactMesh& mesh);
    };
}


#endif // MOIST_INTERFACE_INSERTER_MERGER_HPP_
