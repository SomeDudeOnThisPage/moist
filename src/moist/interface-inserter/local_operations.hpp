#ifndef MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_
#define MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "interface_inserter.hpp"
#include "exact_types.hpp"
#include "exact_mesh.hpp"

namespace moist
{
    namespace operation
    {
        namespace exact
        {
            std::vector<moist::CrossedEdgeExact> FindIntersectedEdges(const moist::exact::EdgePoints& edge, const std::size_t& c, const moist::ExactMesh& mesh);
            void InsertVertexOnCellBoundaryFacet(const std::size_t& c, const std::size_t& v, moist::ExactMesh& mesh);
            bool InsertVertexOnCellBoundaryEdgeOld(const std::size_t& c, const std::size_t& v, moist::ExactMesh& mesh, const double eps = 0.0);
            void InsertVertexOnCellBoundaryEdge(const std::size_t& c, const std::size_t& v, const moist::predicates::PointInTet& edge, moist::ExactMesh& mesh);

            void SplitEdge1_2(const std::size_t& c, const moist::CrossedEdgeExact& edge, moist::ExactMesh& mesh, std::vector<moist::exact::Cell>& created_cells);
            void SplitEdge1_3(const std::size_t& c, const moist::CrossedEdgeExact& edge0, const moist::CrossedEdgeExact& edge1, moist::ExactMesh& mesh, std::vector<moist::exact::Cell>& created_cells);
        }
    }
}

#endif // MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_
