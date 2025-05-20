#ifndef MOIST_CORE_SLICE_IO_HPP_
#define MOIST_CORE_SLICE_IO_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

namespace moist::slice_io
{
    /**
     * @brief Extracts a slice of a larger mesh.
     *
     * @param file
     * @param descriptor
     * @param slice
     */
    void extract(const std::filesystem::path file, const moist::descriptor::MeshSlice descriptor, geogram::Mesh& slice);

    /**
     * @brief Inserts a slice into a larger mesh, merging equal vertices.
     *
     * @param file
     * @param descriptor
     * @param slice
     */
    void insert(const std::filesystem::path file, const moist::descriptor::MeshSlice descriptor, geogram::Mesh& slice);
}

#endif // MOIST_CORE_SLICE_IO_HPP_
