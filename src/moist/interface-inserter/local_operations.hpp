#ifndef MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_
#define MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "slice.hpp"

namespace moist
{
    namespace operation
    {
        /**
         * @brief Splits an edge of the given tetrahedron, resulting in two created tetrahedra.
         *
         * @param mesh Mesh reference.
         * @param cell Cell to split. Will be marked as deleted.
         * @param edge Edge information to perform the split.
         * @param plane Reference plane to detect the non-interface vertex.
         */
        void edge_split_1to2(MeshSlice& mesh, const g_index cell, const CrossedEdge& edge, const AxisAlignedInterfacePlane& plane);

        /**
         * @brief Splits two edges of the given tetrahedron, resulting in three created tetrahedra.
         *
         * @param mesh Mesh reference.
         * @param cell Cell to split. Will be marked as deleted.
         * @param edge0 Edge information to perform the split along one edge.
         * @param edge1 Edge information to perform the split along the other edge.
         * @param plane Reference plane to detect the non-interface vertex.
         */
        void edge_split_1to3(MeshSlice& mesh, const g_index cell, const CrossedEdge& edge0, const CrossedEdge& edge1, const AxisAlignedInterfacePlane& plane);

        /**
         * @brief Not implemented. Theoretically, this situation can arise when an inserted vertex lies exactly on an existing edge.
         * @todo  This should be generalized in one method, see `vertex_insert_1to3`.
         */
        void vertex_insert_1to2();

        /**
         * @brief Inserts a vertex into a given cell on a facet which lies on the given plane, resulting in three created tetrahedra.
         *
         * @param mesh Mesh reference.
         * @param cell Cell to split. Will be marked as deleted.
         * @param point Point to insert.
         * @param plane Reference plane to detect the non-interface vertex.
         *
         * @todo This should be refactored to be more "generalized", so that there needs to be no information about a plane, but rather just
         *       a point and cell. If the point lies on a facet, split into 3, if it lies on an edge, split into 2, if neither split into 4.
         */
        void vertex_insert_1to3(MeshSlice &mesh, const g_index cell, const vec3& point, const moist::AxisAlignedInterfacePlane& plane);
    }
}


#endif // MOIST_INTERFACE_INSERTER_LOCAL_OPERATIONS_HPP_
