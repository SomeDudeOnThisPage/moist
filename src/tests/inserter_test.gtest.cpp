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
//#include "moist/interface-inserter/slice.hpp"

#include "test_utils.hpp"
#include "fixtures/core_group.def.hpp"


static std::filesystem::path _csv;

class InterfaceInserterTestFixture : public ::testing::TestWithParam<moist::testing::fixture::InterfaceInserterFiles>
{
protected:
    static void SetUpTestSuite()
    {
        moist::utils::geogram::initialize();
        const std::time_t now = std::time(nullptr);

        std::ostringstream oss;
        oss << "./metrics/" << std::put_time(std::localtime(&now), "%Y-%m-%d_%H-%M-%S") << ".csv";
        _csv = std::filesystem::path(oss.str());
    }

    void SetUp() override
    {
    }

    void TearDown() override
    {
    }
};

/*TEST_P(InterfaceInserterTestFixture, GenerateStatistics)
{
    const moist::testing::fixture::InterfaceInserterFiles _files = GetParam();

    auto metrics = moist::metrics::Metrics("InterfaceInserter");
    *metrics << moist::metrics::Metric {"test_case", _files[0]};
    moist::Timer timer("InterfaceInserter::Main", metrics);

    auto interface = moist::Interface(_files[1]);
    *metrics << moist::metrics::Metric {"interface::extent", interface.Plane()->extent}
             << moist::metrics::Metric {"interface::epsilon", interface.Plane()->epsilon}
             << moist::metrics::Metric {"interface::nb_vertices", interface.Triangulation()->vertices.nb()}
             << moist::metrics::Metric {"interface::nb_edges", interface.Triangulation()->edges.nb()};

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load(_files[2], slice_a, 3, false);
    moist::utils::geogram::load(_files[3], slice_b, 3, false);

    moist::insert_constraints(slice_a, slice_b, interface, metrics);

    metrics->AppendCSV(_csv);

#ifndef NDEBUG
    moist::utils::geogram::save("test-a.msh", slice_a);
    moist::utils::geogram::save("test-b.msh", slice_b);
#endif // NDEBUG

    *metrics << moist::metrics::Metric { "overlaps", moist::test::contains_overlapping_constraints(slice_a, slice_b, interface) };

}

INSTANTIATE_TEST_SUITE_P(
    Statistics,
    InterfaceInserterTestFixture,
    ::testing::ValuesIn(moist::testing::fixture::CORE_GROUP_INTERFACE_INSERTION)
);

TEST(InserterTest, ContainsNoOverlappingEdgesTightCarbonate)
{
    auto metrics = moist::metrics::Metrics("InterfaceInserter");

    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/tight_carbonate/interface.geogram");

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load("../test/tight_carbonate/raw/0-4.msh", slice_a, 3, false);
    moist::utils::geogram::load("../test/tight_carbonate/raw/4-8.msh", slice_b, 3, false);

    moist::insert_constraints(slice_a, slice_b, interface, metrics);

    geo::mesh_repair(slice_a, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);
    geo::mesh_repair(slice_b, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);

    moist::utils::geogram::save("test-a.msh", slice_a);
    moist::utils::geogram::save("test-b.msh", slice_b);

    EXPECT_EQ(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface), 0);
}


TEST(InserterTest, ContainsNoOverlappingEdgesCubeCylinder)
{
    auto metrics = moist::metrics::Metrics("InterfaceInserter");

    moist::utils::geogram::initialize();
    auto interface = moist::Interface("../test/cube_cylinder/interface.geogram");

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load("../test/cube_cylinder/raw/cube.mesh", slice_a, 3, false);
    moist::utils::geogram::load("../test/cube_cylinder/raw/cylinder.mesh", slice_b, 3, false);

    moist::insert_constraints(slice_a, slice_b, interface, metrics);

    geo::mesh_repair(slice_a, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);
    geo::mesh_repair(slice_b, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);

    moist::utils::geogram::save("test-a.msh", slice_a);
    moist::utils::geogram::save("test-b.msh", slice_b);

    EXPECT_EQ(moist::test::contains_overlapping_constraints(slice_a, slice_b, interface), 0);
}*/
