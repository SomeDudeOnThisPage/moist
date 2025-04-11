#include <filesystem>
#include <algorithm>
#include <regex>
#include <map>

#include <CLI/CLI.hpp>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>

#include <iostream>

#include "../core.hpp"
#include "../core-interface.hpp"
#include "../utils.hpp"

#include "ooc_mesh.hpp"
#include "attributes.inl"
#include "predicates.inl"

int main(int argc, char* argv[])
{
    std::map<std::string, incremental_meshing::Axis> axis_option_map{
        {"X", incremental_meshing::Axis::X},
        {"Y", incremental_meshing::Axis::Y},
        {"Z", incremental_meshing::Axis::Z}
    };

    incremental_meshing::InterfaceExtractionOptions options;
    CLI::App app{argv[0]};

    app.add_option("-m, --mesh", options.path_mesh_a, "Mesh (.mesh) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-i, --interface", options.path_mesh_b, "Interface mesh (.mesh) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o, --output", options.path_mesh_out, "Output Triangulation (.obj) file")
        ->required();
    app.add_option("-x, --axis", options.axis, "Interface axis (X|Y|Z)")
        ->required()
        ->transform(CLI::CheckedTransformer(axis_option_map, CLI::ignore_case));
    app.add_option("-e, --envelope", options.envelope_size, "Envelope size around the interface plane.")
        ->required();
    app.add_option("-p, --plane", options.plane, "Interface plane along the axis (-x)")
        ->required();

    // TODO: Define a unified "slice"-definition structure, so that all "sub-programs" have the same definition of what a data-slice is.
    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::Logger::instance()->set_quiet(true);

    // register custom attributes
    geogram::geo_register_attribute_type<incremental_meshing::InterfaceVertexStrategy>(incremental_meshing::INTERFACE_VERTEX_STRATEGY_ATTRIBUTE);
    geogram::geo_register_attribute_type<double>(incremental_meshing::INTERFACE);

    incremental_meshing::attributes::initialize();

    geogram::MeshIOFlags flags;
    flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS |
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS
    ));

    incremental_meshing::SubMesh mesh_a("a");
    if (!geogram::mesh_load(options.path_mesh_a, mesh_a))
    {
        OOC_ERROR("Failed to load mesh A: " << options.path_mesh_a);
        return 1;
    }

    auto interface = incremental_meshing::Interface(options.path_mesh_b, incremental_meshing::AxisAlignedInterfacePlane {
        options.axis,
        options.plane,
        options.envelope_size
    });

    mesh_a.InsertInterface(interface);
    mesh_a.assert_is_valid();

    return 0;
}
