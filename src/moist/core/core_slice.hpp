#ifndef MOIST_CORE_MESH_SLICE_HPP_
#define MOIST_CORE_MESH_SLICE_HPP_

#include <unordered_set>
#include <initializer_list>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

namespace moist
{
    class MeshSlice : public geogram::Mesh
    {
    public:
        /**
         * @brief Constructs a MeshSlice object.
         *
         * Initializes a new MeshSlice with the given dimension and precision settings.
         *
         * @param dimension The dimension of the mesh (default `3`).
         * @param single_precision If true, use `float`; otherwise, use `double` (default `false`).
         */
        MeshSlice(geogram::index_t dimension = 3, bool single_precision = false);
    };
}

#endif // MOIST_CORE_MESH_SLICE_HPP_
