#ifndef MOIST_CORE_SLICE_IO_HPP_
#define MOIST_CORE_SLICE_IO_HPP_

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

#include <nlohmann/json.hpp>

namespace moist::slice_io
{
    namespace msh
    {
        void save(const std::filesystem::path file, const geogram::Mesh& slice, const geogram::MeshIOFlags ioflags);
    }
}

#endif // MOIST_CORE_SLICE_IO_HPP_
