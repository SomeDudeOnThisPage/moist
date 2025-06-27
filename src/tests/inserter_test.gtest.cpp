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

TEST(InserterTest, ContainsConstraintsCube)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("../test/cube_cylinder/cube_cylinder.geogram");

    moist::MeshSlice slice;
    moist::utils::geo::load("../test/cube_cylinder/raw/cube.mesh", slice, 3, false);
    auto steiner_points = slice.InsertInterface(interface);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_constraints(slice, *interface.Triangulation(), interface, steiner_points.size()));
}


TEST(InserterTest, ContainsConstraintsTightCarbonate)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("../test/tight_carbonate/interface.geogram");

    moist::MeshSlice slice;
    moist::utils::geo::load("../test/tight_carbonate/raw/0-4.msh", slice, 3, false);
    auto steiner_points = slice.InsertInterface(interface);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_constraints(slice, *interface.Triangulation(), interface, steiner_points.size()));
}


TEST(InserterTest, ContainsNoOverlappingEdgesTightCarbonate)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("../test/tight_carbonate/interface.geogram");

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geo::load("../test/tight_carbonate/raw/0-4.msh", slice_a, 3, false);
    slice_a.InsertInterface(interface);

    moist::utils::geo::load("../test/tight_carbonate/raw/4-8.msh", slice_a, 3, false);
    slice_b.InsertInterface(interface);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface));
}


TEST(InserterTest, ContainsNoOverlappingEdgesCubeCylinder)
{
    moist::utils::geo::initialize();
    auto interface = moist::Interface("../test/cube_cylinder/cube_cylinder.geogram");

    moist::SteinerPoints steiner_points;
    moist::MeshSlice slice_a, slice_b;
    moist::utils::geo::load("../test/cube_cylinder/raw/cube.mesh", slice_a, 3, false);

    auto sps = slice_a.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    moist::utils::geo::load("../test/cube_cylinder/raw/cylinder.mesh", slice_b, 3, false);
    sps = slice_b.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    geogram::mesh_repair(slice_a);
    geogram::mesh_repair(slice_b);

    if (!steiner_points.empty())
    {
        slice_a.FlushTetrahedra(true);
        slice_b.FlushTetrahedra(true);
        OOC_DEBUG("inserting " << steiner_points.size() << " steiner points into both meshes...");
        for (const vec3& steiner_point : steiner_points)
        {
            slice_a.InsertVertex(steiner_point, *interface.Plane());
            slice_b.InsertVertex(steiner_point, *interface.Plane());
        }

    }

    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface));
}
