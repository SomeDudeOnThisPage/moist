#include <array>
#include <filesystem>
#include <format>
#include <iostream>

#include <gtest/gtest.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"
#include "moist/core/utils.hpp"

#include "moist/submesh-merger/submesh_merger.hpp"

TEST(MergeTest, ContainsAllEdges)
{
    size_t missing_facets = 0;
    size_t found_facets = 0;

    // load test meshes
    moist::utils::geo::initialize();

    moist::Merger merger("./data/cube_inserted.msh", "./data/cylinder_inserted.msh");
    merger.Merge();
    merger.CopyToOriginal("./data/cube_cylinder.msh");

    geogram::Mesh merged;
    moist::utils::geo::load("./data/cube_cylinder.merged.msh", merged);
    auto interface = moist::Interface("./data/interface_inserted.geogram");

    geogram::mesh_repair(merged);
    geogram::parallel_for(0, merged.cells.nb(), [&merged, &interface, &found_facets] (const g_index c)
    {
        if (!moist::predicates::cell_on_plane(c, merged, *interface.Plane()))
        {
            return;
        }

        auto triangulation = interface.Triangulation();
        for (const g_index f : interface.Triangulation()->facets)
        {
            if (!moist::predicates::facet_matches_cell(c, f, merged, *triangulation))
            {
                continue;
            }

            LOCK_ATTRIBUTES;
            geogram::Attribute<double> f_found(interface.Triangulation()->facets.attributes(), ATTRIBUTE_INTERFACE_FACET_FOUND_AFTER_MERGE);
            if (!f_found[f])
            {
                f_found[f] = 1;
                found_facets++;
            }
            break;
        }
    });

    geogram::Attribute<double> f_found(interface.Triangulation()->facets.attributes(), ATTRIBUTE_INTERFACE_FACET_FOUND_AFTER_MERGE);
    for (const g_index f : interface.Triangulation()->facets)
    {
        EXPECT_EQ(f_found[f], 1);
    }

    std::filesystem::remove("./data/cube_cylinder.merged.msh");

#ifndef NDEBUG
    moist::utils::geo::save("./debug/test/interface_inserted.merge.geogram", *interface.Triangulation());
#endif
}
