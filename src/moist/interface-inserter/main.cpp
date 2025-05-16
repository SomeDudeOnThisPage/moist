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
#include <geogram/basic/environment.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"

#include "slice.hpp"

struct Arguments
{
    std::filesystem::path input;
    std::filesystem::path interface;
    std::filesystem::path output;
    double plane;
    double envelope_size;
    moist::Axis axis;
};

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.input, "Mesh (.mesh|.msh|.geogram) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-i, --interface", arguments.interface, "Interface mesh (.mesh|.msh|.obj|.off|.geogram) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-o, --output", arguments.output, "Output Triangulation (.obj) file");
    app.add_option("-x, --axis", arguments.axis, "Interface axis (X|Y|Z)")
        ->required()
        ->transform(cli::CheckedTransformer(moist::AXIS_OPTION_ARGUMENT_MAP, cli::ignore_case));
    app.add_option("-e, --envelope", arguments.envelope_size, "Envelope size around the interface plane.")
        ->required();
    app.add_option("-p, --plane", arguments.plane, "Interface plane along the axis (-x)")
        ->required();

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::CmdLine::import_arg_group("sys"); // needs to be called in order to be able to export .geogram meshes...
    //geogram::CmdLine::set_arg("sys:compression_level", "0");
    geogram::Logger::instance()->set_quiet(true);
    moist::attributes::initialize();

    geogram::MeshIOFlags flags;
    flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS |
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS
    ));

    moist::MeshSlice slice;
    if (!geogram::mesh_load(arguments.input, slice))
    {
        OOC_ERROR("Failed to load mesh A: " << arguments.input);
        return 1;
    }

    auto interface = moist::Interface(arguments.interface, moist::AxisAlignedInterfacePlane {
        arguments.axis,
        arguments.plane,
        arguments.envelope_size
    });

    slice.InsertInterface(interface);
    slice.assert_is_valid();

    GEO::MeshIOFlags export_flags;
    export_flags.set_attribute(geogram::MESH_NO_ATTRIBUTES);
    export_flags.set_dimension(3);
    export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    export_flags.set_verbose(true);
    geogram::mesh_save(slice, "_test_export.geogram", export_flags);

    return 0;
}
