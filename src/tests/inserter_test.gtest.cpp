#include <filesystem>
#include <format>
#include <iostream>

#include <gtest/gtest.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/utils.hpp"

#include "moist/interface-inserter/slice.hpp"

#include "test_utils.hpp"

TEST(InserterTest, ContainsConstraints)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("./data/interface_inserted.geogram");

    moist::MeshSlice slice;
    moist::utils::geo::load("./data/cube.mesh", slice, 3, false);
    slice.InsertInterface(interface);
    geogram::mesh_repair(slice);

    EXPECT_TRUE(moist::test::contains_constraints(slice, *interface.Triangulation()));
}
