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

#include "moist/submesh-merger/submesh_merger.hpp"

#include "test_utils.hpp"

TEST(MergeTest, ContainsConstraints)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("./data/interface_inserted.geogram");

    moist::Merger merger("./data/cube_inserted.msh", "./data/cylinder_inserted.msh");
    merger.Merge();
    merger.CopyToOriginal("./debug/test/cube_cylinder.msh");

    geogram::Mesh merged;
    moist::utils::geo::load("./debug/test/cube_cylinder.msh", merged);
    geogram::mesh_repair(merged);

    EXPECT_TRUE(moist::test::contains_constraints(merged, *interface.Triangulation()));

    std::filesystem::remove("./debug/test/cube_cylinder.msh");
#ifndef NDEBUG
    moist::utils::geo::save("./debug/test/interface_inserted.merge.geogram", *interface.Triangulation());
#endif
}
