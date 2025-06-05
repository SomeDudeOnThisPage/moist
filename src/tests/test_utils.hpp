#ifndef MOIST_TESTS_TEST_UTILS_HPP_
#define MOIST_TESTS_TEST_UTILS_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

namespace moist::test
{
    bool contains_constraints(geogram::Mesh& mesh, geogram::Mesh& constraints);
}

#endif // MOIST_TESTS_TEST_UTILS_HPP_
