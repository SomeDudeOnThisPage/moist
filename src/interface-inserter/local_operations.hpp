#ifndef __OOC_LOCAL_OPERATIONS_HPP
#define __OOC_LOCAL_OPERATIONS_HPP

#include "ooc_mesh.hpp"

namespace incremental_meshing
{
    namespace operation
    {
        /**
         * Splits one edge of a given tetrahedron, placing two tetrahedra into the creation buffer of a mesh.
         */
        void edge_split_1to2(MeshSlice& mesh, const g_index cell, const CrossedEdge& edge, const AxisAlignedInterfacePlane& plane);
        void edge_split_1to3(MeshSlice& mesh, const g_index cell, const CrossedEdge& edge0, const CrossedEdge& edge1, const AxisAlignedInterfacePlane& plane);
        void vertex_insert_1to2();
        void vertex_insert_2to4(MeshSlice &mesh, const g_index cell, const vec3& point, const incremental_meshing::AxisAlignedInterfacePlane& plane);
    }
}


#endif // __OOC_LOCAL_OPERATIONS_HPP
