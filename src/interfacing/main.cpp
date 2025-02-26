#include <filesystem>

#include <CLI/CLI.hpp>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <iostream>

#include "../core.hpp"
#include "../utils.hpp"

#include "interface.hpp"

// Decide if this should be a wmtk "component" or just "use" wmtk. I'm leaning to option #2, as this project does a lot of stuff that isn't "core"
// wmtk related, like sizing field generation, surface extraction or mesh merging...
// That way, we also wouldn't need to build all of wmtk, like the applications or components we do not need, reducing project complexity.
int main(int argc, char* argv[])
{
    ooc::InterfaceExtractionOptions options;
    CLI::App app{argv[0]};

    app.add_option("-a, --input-a", options.path_mesh_a, "Mesh 'A' (.mesh) file")
        ->required(true)
        ->check(CLI::ExistingFile);
    app.add_option("-b, --input-b", options.path_mesh_b, "Mesh 'B' (.mesh) file")
        ->required(true)
        ->check(CLI::ExistingFile);
    app.add_option("-o, --output", options.path_mesh_out, "Output Triangulation (.obj) file")
        ->required(true);

    // TODO: This should realistically be some form of "Plane-Mesh", i.e. a mesh consisting of four edges.
    //       That way we could do some BSP fuckery and have multiple "intersecting" interfaces.
    app.add_option("-p, --plane", options.plane, "Interface plane defined by one vertex and 2 spanning vectors (vx, vy, vz) (x0, y0, z0), (x1, y1, z1)")
        ->required(true); // this option does nothing yet...

    // TODO: Define a unified "slice"-definition structure, so that all "sub-programs" have the same definition of what a data-slice is.
    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    // TODO: From here, all of this should be part of some library so it's also usable via call from the WildmeshingToolkit.
    //       This all should be moved into a wildmeshing-toolkit-component called "InterfaceExtraction" or something.
    GEO::initialize(GEO::GEOGRAM_INSTALL_NONE); // why can I not disable the geogram logger??

    GEO::MeshIOFlags flags;
    flags.set_elements(GEO::MeshElementsFlags(
        GEO::MeshElementsFlags::MESH_VERTICES |
        GEO::MeshElementsFlags::MESH_FACETS
    ));

    GEO::Mesh mesh_a, mesh_b;
    if (!GEO::mesh_load(options.path_mesh_a, mesh_a))
    {
        OOC_ERROR("Failed to load mesh A: " << options.path_mesh_a);
        return 1;
    }

    if (!GEO::mesh_load(options.path_mesh_b, mesh_b))
    {
        OOC_ERROR("Failed to load mesh B: " << options.path_mesh_b);
        return 1;
    }

    // TODO: Parse Plane from cmd, const for testing.
    //       Also, non AA-Planes w/ epsilon-tolerance along normal...
    //       This could also be a part of the "slice"-definition
    const auto plane = ooc::InterfacePlane {
        ooc::utils::parse_vector("(0.0, 0.0, -1.0)"),
        ooc::utils::parse_vector("(1.0, 0.0, 0.0)"),
        ooc::utils::parse_vector("(0.0, 1.0, 0.0)")
    };

    auto interface = ooc::Interface(plane);
    interface.AddConstraints(mesh_a);
    interface.AddConstraints(mesh_b);
    interface.Triangulate();

#ifndef NDEBUG
    ooc::export_delaunay("delaunay.obj", *interface.Triangulation());
#endif

    // TODO: insert triangulation (do this in separate process maybe?)... use wmtk for this...
    // TODO: decimate non_triangulation + non_boundary edges... use wmtk for this... here we kinda need both meshes again...
    //       so maybe just start 2 threads here and assume we have enough memory to have 2 slices in it at once.
    return 0;
}
