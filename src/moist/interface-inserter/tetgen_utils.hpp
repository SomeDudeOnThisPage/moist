#ifndef MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_
#define MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_

#include <geogram/mesh/mesh.h>
#include <tetgen.h>

#include "moist/core/defines.hpp"

#include "exact_mesh.hpp"

namespace moist::tetgen
{
    void transform(moist::ExactMesh& mesh, tetgenio& tetgen_io);
    void transform(const tetgenio& from, geo::Mesh& mesh);
}

#endif // MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_
