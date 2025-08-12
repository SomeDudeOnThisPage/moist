#include <iostream>
#include <fstream>

#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/command_line.h>

#include "moist/core/defines.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/descriptor.hpp"
#include "moist/core/core_interface.hpp"

#include "submesh_merger.hpp"

struct Arguments
{
    std::filesystem::path mesh;
    std::filesystem::path interface;
    std::filesystem::path slice;
};

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.mesh, "Mesh (.msh) file.")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-s, --slice", arguments.slice, "Slice (.msh) file.")
        ->required()
        ->check(cli::ExistingFile);

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geogram::initialize();

    geo::CmdLine::import_arg_group("algo");
    geo::CmdLine::import_arg_group("output");
    geo::CmdLine::set_arg("output:precision", "16");
    geo::CmdLine::set_arg("sys:use_doubles", true);

    // TODO: Insert Steiner Points of Slices into each other...
    //geogram::Mesh slice;
    //moist::utils::geo::load(arguments.slice, slice);

    moist::Merger merger(arguments.mesh, arguments.slice);
    merger.Merge();
    merger.CopyToOriginal();
}
