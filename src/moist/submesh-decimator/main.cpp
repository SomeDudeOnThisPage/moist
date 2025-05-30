#include <iostream>
#include <limits>
#include <vector>
#include <cmath>

#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/environment.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include "moist/core/defines.hpp"
#include "moist/core/utils.hpp"
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

static vec3 middle(const vec3 v0, const vec3 v1)
{
    return vec3(
        (v0.x + v1.x) / 2.0,
        (v0.y + v1.y) / 2.0,
        (v0.z + v1.z) / 2.0
    );
}

    struct _Edge
    {
        vec3 p0;
        vec3 p1;

        bool operator==(const _Edge& other) const
        {
            return (p0 == other.p0 && p1 == other.p1) ||
                   (p0 == other.p1 && p1 == other.p0);
        }
    };

static void decimate(geogram::Mesh& slice, const std::vector<g_index>& facets, moist::Interface interface, const size_t n)
{
    auto tri = interface.Triangulation();
    std::vector<_Edge> to_collapse;
    size_t to_decimate_n = 0;
    // for each interface facet to be decimated add the shortest edge of the facet...
    for (int i = 0; i < n; i++)
    {
        const g_index f = facets[i];
        _Edge edge {};
        double shortest_edge_length = std::numeric_limits<double>::max();
        for (l_index lv = 0; lv < 3; lv++)
        {
            const vec3 p0 = tri->vertices.point(tri->facets.vertex(f, lv));
            const vec3 p1 = tri->vertices.point(tri->facets.vertex(f, (lv + 1) % 3));
            double distance = std::fabs(geogram::distance(p0, p1));
            if (distance < shortest_edge_length)
            {
                shortest_edge_length = distance;
                edge.p0 = p0;
                edge.p1 = p1;
            }
        }
        to_collapse.push_back(edge);
        to_decimate_n++;
    }

    to_collapse.erase(std::unique(to_collapse.begin(), to_collapse.end()), to_collapse.end());

    size_t decimated = 0;
    for (const g_index cell : slice.cells)
    {
        // check if slice contains one of the edges to decimate
        for (l_index le = 0; le < slice.cells.nb_edges(cell); le++)
        {
            const _Edge edge
            {
                slice.vertices.point(slice.cells.edge_vertex(cell, le, 0)),
                slice.vertices.point(slice.cells.edge_vertex(cell, le, 1))
            };

            if (std::find(to_collapse.begin(), to_collapse.end(), edge) != to_collapse.end())
            {
                const vec3 mid = middle(edge.p0, edge.p1);
                for (g_index v : slice.vertices)
                {
                    vec3 p = slice.vertices.point(v);
                    if (p == edge.p0 || p == edge.p1)
                    {
                        slice.vertices.point(v).x = mid.x;
                        slice.vertices.point(v).y = mid.y;
                        slice.vertices.point(v).z = mid.z;

                        slice.vertices.point(v).x = mid.x;
                        slice.vertices.point(v).y = mid.y;
                        slice.vertices.point(v).z = mid.z;
                        decimated++;
                    }
                }
            }
        }
    }

    OOC_DEBUG("to_decimate_n " << to_decimate_n << " decimated " << decimated);

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

    moist::utils::geo::initialize();

    auto interface = moist::Interface(arguments.interface);

    auto triangulation = interface.Triangulation();
    geogram::Attribute<int> m_target_vertices(interface.Triangulation()->cells.attributes(), ATTRIBUTE_INTERFACE_TARGET_VERTICES);
    const size_t n = m_target_vertices[0];

    // Copy to new vector
    std::vector<g_index> facets;
    for (const g_index facet : triangulation->facets)
    {
        facets.push_back(facet);
    }

    std::sort(facets.begin(), facets.end(), [&triangulation](const g_index& a, const g_index& b)
    {
        geogram::Attribute<double> f_quality(triangulation->facets.attributes(), ATTRIBUTE_INTERFACE_TETMERGE_QUALITY);
        return f_quality[a] > f_quality[b]; // smaller is better in this case actually
    });

    geogram::Mesh slice;
    moist::utils::geo::load(arguments.a, slice);
    decimate(slice, facets, interface, n);
    moist::utils::geo::save(arguments.a.replace_extension(".decimated.msh"), slice);
}
