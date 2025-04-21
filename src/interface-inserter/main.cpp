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

#include "core.hpp"
#include "core-interface.hpp"
#include "utils.hpp"

#include "ooc_mesh.hpp"
#include "attributes.inl"
#include "predicates.inl"

struct Arguments
{
    std::filesystem::path input;
    std::filesystem::path interface;
    std::filesystem::path output;
    double plane;
    double envelope_size;
    incremental_meshing::Axis axis;
};

std::map<std::string, incremental_meshing::Axis> AXIS_OPTION_ARGUMENT_MAP {
    {"X", incremental_meshing::Axis::X},
    {"Y", incremental_meshing::Axis::Y},
    {"Z", incremental_meshing::Axis::Z}
};

int main(int argc, char* argv[])
{
    Arguments arguments{};
    CLI::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.input, "Mesh (.mesh|.msh|.geogram) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-i, --interface", arguments.interface, "Interface mesh (.mesh|.msh|.obj|.off|.geogram) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o, --output", arguments.output, "Output Triangulation (.obj) file");
    app.add_option("-x, --axis", arguments.axis, "Interface axis (X|Y|Z)")
        ->required()
        ->transform(CLI::CheckedTransformer(AXIS_OPTION_ARGUMENT_MAP, CLI::ignore_case));
    app.add_option("-e, --envelope", arguments.envelope_size, "Envelope size around the interface plane.")
        ->required();
    app.add_option("-p, --plane", arguments.plane, "Interface plane along the axis (-x)")
        ->required();

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::Logger::instance()->set_quiet(true);
    incremental_meshing::attributes::initialize();

    geogram::MeshIOFlags flags;
    flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS |
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS
    ));

    incremental_meshing::MeshSlice slice;
    if (!geogram::mesh_load(arguments.input, slice))
    {
        OOC_ERROR("Failed to load mesh A: " << arguments.input);
        return 1;
    }

    auto interface = incremental_meshing::Interface(arguments.interface, incremental_meshing::AxisAlignedInterfacePlane {
        arguments.axis,
        arguments.plane,
        arguments.envelope_size
    });

    slice.InsertInterface(interface);
    slice.assert_is_valid();

    return 0;
}
