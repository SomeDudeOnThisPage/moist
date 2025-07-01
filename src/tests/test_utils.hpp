#ifndef MOIST_TESTS_TEST_UTILS_HPP_
#define MOIST_TESTS_TEST_UTILS_HPP_

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"

namespace moist::test
{
    size_t contains_overlapping_constraints(geogram::Mesh& a, geogram::Mesh& b, moist::Interface& interface);
    bool contains_constraints(geogram::Mesh& mesh, geogram::Mesh& constraints, moist::Interface& interface, size_t steiner_points = 0);
}

#endif // MOIST_TESTS_TEST_UTILS_HPP_
