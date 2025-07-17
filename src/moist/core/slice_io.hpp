#ifndef MOIST_CORE_SLICE_IO_HPP_
#define MOIST_CORE_SLICE_IO_HPP_

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

namespace moist::slice_io
{
    namespace msh
    {
        void save(const std::filesystem::path file, const GEO::Mesh& slice, const GEO::MeshIOFlags ioflags);
    }
}

#endif // MOIST_CORE_SLICE_IO_HPP_
