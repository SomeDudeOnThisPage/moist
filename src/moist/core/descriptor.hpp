#ifndef MOIST_CORE_DESCRIPTOR_HPP_
#define MOIST_CORE_DESCRIPTOR_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

namespace moist
{
    /**
     * @brief Contains a collection of metadata descriptors for Interfaces, MeshSlices and GlobalMeshes.
     */
    namespace descriptor
    {
        struct BoundingBox
        {
            vec3 min;
            vec3 max;
        };

        /**
         * @brief Maps an interface triangulation cell to an index insided LocalInterface, i.e. maps a cell quality value of the local mesh to the interface.
         *
         * Used in the final decimation procedure, to determine cells that should be collapsed in both meshes.
         */
        struct InterfaceCellQuality
        {
            g_index facet_i;
            g_index cell_a;
            g_index cell_b;
            float quality_a;
            float quality_b;
        };

        /**
         * @brief Describes a local interface in a MeshSlice.
         *
         * The interface vertices must be laid out continuously in memory, starting with tthe first index equal to `v_start`, the last to `v_start + v_count - 1`.
         * All edges between these vertices are considered interface-edges when merging. Any discrepancy between the two meshes' edges will cause the merge to fail.
         */
        struct LocalInterface
        {
            /** @brief Unique identifier of this interface. */
            std::string index;
            /** @brief First interface-vertex. */
            g_index v_start;
            /** @brief Amount of interface-vertices. */
            g_index v_count;
            /** @brief Flag to check if interface-tets have already been decimated on this interface. */
            bool decimated;
        };

        /**
         * @brief Describes a MeshSlice as part of a GlobalMesh.
         *
         * This struct contains information for indexing into a GlobalMesh file, without having to read/write the entire file.
         * This enables loading parts of a msh file as a slice.
         */
        struct MeshSlice
        {
            /** @brief Unique string identifier of this slice. */
            std::string index;
            /** @brief BoundingBox encompassing this mesh slice, used for preliminary intersection checking. */
            BoundingBox bbox;
            /** @brief Byte offset of the first vertex inside the global mesh file. */
            uint64_t v_byte_start;
            /** @brief Total size of vertices in bytes inside the global mesh file. */
            uint64_t v_byte_size;
            /** @brief Byte offset of the first tetrahedron inside the global mesh file. */
            uint64_t t_byte_start;
            /** @brief Total size of tetrahedra in bytes inside the global mesh file. */
            uint64_t t_byte_size;
            /** @brief Contains a list of LocalInterfaceDescriptors that have been inserted into this MeshSlice. */
            std::vector<LocalInterface> interfaces;
        };

        /**
         * @brief Describes a GlobalMesh as a collection of MeshSlices.
         *
         * This struct is used to merge missing MeshSlices into a final GlobalMesh file.
         */
        struct Mesh
        {
            /**
             * @brief Contains a list of MeshSlices that have been already merged into a GlobalMesh
             */
            std::vector<MeshSlice> slices;
        };
    }
}

#endif // MOIST_CORE_DESCRIPTOR_HPP_
