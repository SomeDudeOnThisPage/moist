#include <iostream>
#include <limits>

#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/environment.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

struct Arguments
{
    std::filesystem::path a;
    std::filesystem::path b;
    std::filesystem::path interface;
};

// Load both Slice descriptors.
// Check the LocalInterface to be merged.
// Determine tets to be decimated by checking both InterfaceCellQuality lists.

// For each Slice:
// Load into memory
// Decimate required tets

static void decimate(geogram::Mesh& slice, const std::vector<g_index>& facets, moist::Interface interface, const size_t n)
{
    // find cells to be decimated...
    for (g_index i = 0; i < n; i++)
    {
        for (const g_index cell : slice.cells)
        {
            if (!moist::predicates::facet_matches_cell(cell, facets[i], slice, *interface.Triangulation()))
            {
                continue;
            }

            // shortest edge collapse (of the three interface edges!)...
            l_index shortest_edge = -1;
            double shortest_edge_length = std::numeric_limits<double>::max();
            bool already_decimated = false;
            for (l_index le = 0; le < slice.cells.nb_edges(cell); le++)
            {
                const vec3 p0 = slice.vertices.point(slice.cells.edge_vertex(cell, le, 0));
                const vec3 p1 = slice.vertices.point(slice.cells.edge_vertex(cell, le, 1));

                // TODO: Disallow collapsing edges that are boundary adjacent (at least one of the two vertices is a boundary vertex of one of the meshes)
                //       Make this an attribute in the interface mae
                if (moist::predicates::point_on_plane(p0, *interface.Plane()) && moist::predicates::point_on_plane(p1, *interface.Plane()))
                {
                    const double distance = geogram::distance(p0, p1);
                    if (distance == 0.0)
                    {
                        already_decimated = true;
                        break;
                    }
                    if (distance < shortest_edge_length)
                    {
                        shortest_edge_length = distance;
                        shortest_edge = le;
                    }
                }
            }

            // decimate le#p1 onto le#p0
            if (shortest_edge == -1 || already_decimated)
            {
                OOC_DEBUG("oi noi");
            }

            const vec3 p0 = slice.vertices.point(slice.cells.edge_vertex(cell, shortest_edge, 0));
            slice.vertices.point(slice.cells.edge_vertex(cell, shortest_edge, 1)).x = p0.x;
            slice.vertices.point(slice.cells.edge_vertex(cell, shortest_edge, 1)).y = p0.y;
            slice.vertices.point(slice.cells.edge_vertex(cell, shortest_edge, 1)).z = p0.z;

            if (p0.z != -1.0)
            {
                OOC_DEBUG("oi noi");
            }
            //break; // edge is decimated, we don't need to check other cells here, they will all be deleted in the next step at once
        }
    }

    geogram::vector<geogram::index_t> deleted(slice.cells.nb());
    for (const auto tet : slice.cells)
    {
        if (moist::geometry::has_duplicate_vertex(tet, slice))
        {
            deleted[tet] = 1;
        }
    }

    slice.cells.delete_elements(deleted, false);
}

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.a, "Mesh (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-i, --interface", arguments.interface, "Interface mesh (.geogram) file")
        ->required()
        ->check(cli::ExistingFile);

    // determine n to decimate...
    // mesh A has X vertices inserted into interface, mesh B Y
    // so nb_vertices = X + Y, target = (X + Y - same(X, Y)) / 2
    // or in other ways, n = (unique(X) + unique(Y)) / 2

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::CmdLine::import_arg_group("sys"); // needs to be called in order to be able to export .geogram meshes...
    //geogram::CmdLine::set_arg("sys:compression_level", "0");
    geogram::Logger::instance()->set_quiet(true);

    auto interface = moist::Interface(arguments.interface, moist::AxisAlignedInterfacePlane {
        moist::Axis::Z,
        -1.0,
        0.0
    });
    geogram::MeshIOFlags flags;
    flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS |
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS
    ));

    GEO::MeshIOFlags export_flags;
    export_flags.set_attribute(geogram::MESH_ALL_ATTRIBUTES);
    export_flags.set_dimension(3);
    export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    export_flags.set_verbose(true);

    auto triangulation = interface.Triangulation();
    geogram::Attribute<int> m_target_vertices(interface.Triangulation()->cells.attributes(), ATTRIBUTE_INTERFACE_TARGET_VERTICES);
    const size_t n = m_target_vertices[0];

    // Copy to new vector
    std::vector<g_index> facets;
    for (const g_index facet : triangulation->facets)
    {
        facets.push_back(facet);
    }

    std::sort(facets.begin(), facets.end(), [triangulation](const g_index& a, const g_index& b)
    {
        geogram::Attribute<double> f_quality(triangulation->facets.attributes(), ATTRIBUTE_INTERFACE_TETMERGE_QUALITY);
        return f_quality[a] > f_quality[b]; // smaller is better in this case actually
    });

    geogram::Mesh slice;
    if (!geogram::mesh_load(arguments.a, slice))
    {
        OOC_ERROR("Failed to load mesh A: " << arguments.a);
        return 1;
    }

    decimate(slice, facets, interface, n);
    geogram::mesh_save(slice, arguments.a.replace_extension(".decimated.msh"), export_flags);
}
