#include <filesystem>
#include <format>
#include <iostream>

#include <gtest/gtest.h>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/process.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
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
    moist::utils::geo::load("./data/cylinder_inserted.msh", merged);
    auto interface = moist::Interface("./data/interface_inserted.geogram");

    auto triangulation = interface.Triangulation();
    std::cerr << merged.cells.nb() << std::endl;
    std::cerr << merged.facets.nb() << std::endl;
    std::cerr << interface.Triangulation()->facets.nb() << std::endl;

    //geogram::parallel_for(0, triangulation->facets.nb(), [&triangulation, &merged, &found_facets, &missing_facets](const g_index f) {
    for (g_index f = 0; f < triangulation->facets.nb(); f++)
    {
        const vec3 f0 = triangulation->vertices.point(triangulation->facets.vertex(f, 0));
        const vec3 f1 = triangulation->vertices.point(triangulation->facets.vertex(f, 1));
        const vec3 f2 = triangulation->vertices.point(triangulation->facets.vertex(f, 2));

        bool f_found = false;
        for (const g_index c : merged.cells)
        {
            for (l_index lf = 0; lf < merged.cells.nb_facets(c); lf++)
            {
                geogram::index_t matching = 0;
                for (l_index lv = 0; lv < merged.cells.facet_nb_vertices(c, lf); lv++)
                {
                    g_index lfv = merged.cells.facet_vertex(c, lf, lv);
                    vec3 point = merged.vertices.point(lfv);
                    if (point.x == f0.x && point.y == f0.y && f0.z == -1.0 || point.x == f1.x && point.y == f1.y && f1.z == -1.0 || point.x == f2.x && point.y == f2.y && f2.z == -1.0)
                    {
                        matching++;
                    }
                }
                if (matching == 3)
                {
                    LOCK_ATTRIBUTES;
                    geogram::Attribute<int> ff_found(interface.Triangulation()->facets.attributes(), ATTRIBUTE_INTERFACE_FACET_FOUND_AFTER_MERGE);
                    ff_found[f] = 1;
                    std::string fmt = std::format("found (({},{}), ({},{}), ({},{})",
                        -f0.x, -f0.y,
                        -f1.x, -f1.y,
                        -f2.x, -f2.y
                    );
                    f_found = true;
                    //std::cerr << fmt << std::endl;
                    found_facets++;
                    goto exit_;
                }
            }
            /*if (moist::predicates::facet_matches_cell(c, f, merged, *triangulation))
            {
                std::string fmt = std::format("found (({},{}), ({},{}), ({},{})",
                    -f0.x, -f0.y,
                    -f1.x, -f1.y,
                    -f2.x, -f2.y
                );
                f_found = true;
                std::cerr << fmt << std::endl;
                found_facets++;
            }*/
        }
exit_:
        if (!f_found)
        {
            missing_facets++;
        }
    }
        // RecordProperty("MissingEdge", std::format("missing [{}, {}, {}] -> [{}, {}, {}]", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z));
        //std::cerr << std::format("missing [{}, {}, {}] -> [{}, {}, {}]", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z) << std::endl;
    //});

    EXPECT_EQ(missing_facets, 0);
    EXPECT_EQ(found_facets, triangulation->facets.nb());

    moist::utils::geo::save("interface_testexport.geogram", *interface.Triangulation());
}
