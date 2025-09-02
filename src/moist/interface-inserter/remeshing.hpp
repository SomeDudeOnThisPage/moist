#ifndef MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_
#define MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_

#include <geogram/mesh/mesh.h>
#include <tetgen.h>
extern "C"
{
    #include "mmg/mmg3d/libmmg3d.h"
}

#include "moist/core/defines.hpp"

#include "exact_mesh.hpp"

namespace moist::tetgen
{
    void transform(moist::ExactMesh& mesh, tetgenio& tetgen_io);
    void transform(tetgenio& from, moist::ExactMesh& to, const moist::AxisAlignedPlane& plane);

    void coarsen(tetgenio& from, tetgenio& to);
    // void remesh();
}

namespace moist::mmg3d
{
    void transform(moist::ExactMesh& from, MMG5_pMesh to, MMG5_pSol solution);
    void transform(MMG5_pMesh from, MMG5_pSol solution, moist::ExactMesh& to);
    void set_solution(MMG5_pMesh mesh, MMG5_pSol solution, double hmin, double hmax, double hausd = 0.05);

    void remesh(MMG5_pMesh mesh, MMG5_pSol solution);

    namespace metrics
    {
        double avg_quality_mmg3d(const moist::ExactMesh& mesh);
    }
}

#endif // MOIST_INTERFACE_INSERTER_TETGEN_UTILS_HPP_
