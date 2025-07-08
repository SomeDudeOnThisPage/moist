#include <filesystem>
#include <format>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include <gtest/gtest.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/utils.hpp"

#include "moist/interface-inserter/interface_inserter.hpp"
#include "moist/interface-inserter/slice.hpp"

#include "test_utils.hpp"

struct Files
{
    std::filesystem::path file_a;
    std::filesystem::path file_b;
    std::filesystem::path file_interface;
};

std::filesystem::path _csv;
std::vector<Files> parameters;

class InterfaceInserterTestFixture : public ::testing::TestWithParam<Files>
{
protected:
    static void SetUpTestSuite()
    {
        moist::utils::geogram::initialize();
        const std::time_t now = std::time(nullptr);

        std::ostringstream oss;
        oss << "data/out/" << std::put_time(std::localtime(&now), "%Y-%m-%d_%H-%M-%S") << ".csv";

        _csv = std::filesystem::path(oss.str());
    }

    void SetUp() override
    {
    }

    void TearDown() override
    {
    }
};

TEST_P(InterfaceInserterTestFixture, GenerateStatistics)
{
    const Files _files = GetParam();
}

INSTANTIATE_TEST_SUITE_P(
    InterfaceInserterTestFixtures,
    InterfaceInserterTestFixture,
    ::testing::Values(
        Files{"cube_cylinder/interface.geogram", "cube_cylinder/0.msh", "cube_cylinder/1.msh"},
        Files{"tight_carbonate/interface.geogram", "tight_carbonate/0.msh", "tight_carbonate/1.msh"},
        Files{"tight_carbonate/interface_inverted.geogram", "tight_carbonate/0_inverted.msh", "tight_carbonate/1_inverted.msh"}
    )
);

TEST(InserterTest, ContainsConstraintsCube)
{
    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/cube_cylinder/cube_cylinder.geogram");

    moist::MeshSlice slice;
    moist::utils::geogram::load("../test/cube_cylinder/raw/cube.mesh", slice, 3, false);
    auto steiner_points = slice.InsertInterface(interface);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_constraints(slice, *interface.Triangulation(), interface, steiner_points.size()));
}


TEST(InserterTest, ContainsConstraintsTightCarbonate)
{
    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/tight_carbonate/interface.geogram");

    moist::MeshSlice slice;
    moist::utils::geogram::load("../test/tight_carbonate/raw/0-4.msh", slice, 3, false);
    auto steiner_points = slice.InsertInterface(interface);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_TRUE(moist::test::contains_constraints(slice, *interface.Triangulation(), interface, steiner_points.size()));
}


TEST(InserterTest, ContainsNoOverlappingEdgesTightCarbonate)
{
    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/tight_carbonate/interface.geogram");

    moist::SteinerPoints steiner_points;
    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load("../test/tight_carbonate/raw/0-4.msh", slice_a, 3, false);
    auto sps = slice_a.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    moist::utils::geogram::load("../test/tight_carbonate/raw/4-8.msh", slice_b, 3, false);
    sps = slice_b.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    //

    moist::utils::geogram::save("test-a-before-steiner.msh", slice_a);
    moist::utils::geogram::save("test-b-before-steiner.msh", slice_b);

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

    geo::mesh_repair(slice_a, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);
    geo::mesh_repair(slice_b, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);
    moist::utils::geogram::save("test-a.msh", slice_a);
    moist::utils::geogram::save("test-b.msh", slice_b);

    EXPECT_EQ(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface), 0);
}


TEST(InserterTest, ContainsNoOverlappingEdgesCubeCylinder)
{
    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/cube_cylinder/cube_cylinder.geogram");

    moist::SteinerPoints steiner_points;
    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load("../test/cube_cylinder/raw/cube.mesh", slice_a, 3, false);

    auto sps = slice_a.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    moist::utils::geogram::load("../test/cube_cylinder/raw/cylinder.mesh", slice_b, 3, false);
    sps = slice_b.InsertInterface(interface);
    steiner_points.insert(sps.begin(), sps.end());

    geo::mesh_repair(slice_a);
    geo::mesh_repair(slice_b);

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

    moist::utils::geogram::save("test-a.msh", slice_a);
    moist::utils::geogram::save("test-b.msh", slice_b);
    // geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);

    EXPECT_EQ(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface), 0);
}
