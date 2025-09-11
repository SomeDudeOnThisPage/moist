#include <filesystem>
#include <algorithm>
#include <regex>
#include <unordered_set>
#include <map>

#include <CLI/CLI.hpp>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/command_line.h>

#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/metrics.hpp"
#include "moist/core/mesh_quality.inl"

#include "interface_inserter.hpp"
#include "merger.hpp"

struct Arguments
{
    std::filesystem::path input_a;
    std::filesystem::path input_b;

    double epsilon;
    double extent;
    std::filesystem::path output_a;
    std::filesystem::path output_b;

    std::filesystem::path metrics;
    std::string tc_name;

    // Eval
    bool use_static_grid_size;
    float grid_size_factor;

    // Eval
    double hmin_factor;
    double hmax_factor;
};

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-a, --mesh-a", arguments.input_a, "Mesh A (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-b, --mesh-b", arguments.input_b, "Mesh A (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-e, --epsilon", arguments.epsilon, "Interface epsilon")
        ->required();
    app.add_option("-x, --extent", arguments.extent, "Interface extent along the z axis")
        ->required();
    app.add_option("--output-a", arguments.output_a, "Output Triangulation A (.msh) file")
        ->required();
    app.add_option("--output-b", arguments.output_b, "Output Triangulation B (.msh) file")
        ->required();
    app.add_option("--grid-factor", arguments.grid_size_factor, "Override grid size factor")
        ->default_val(1.0);
    const auto* use_metrics = app.add_option("--output-metrics", arguments.metrics, "Output Metrics (.csv) file");
    app.add_option("--output-metrics-tc-name", arguments.tc_name, "Test-Case Name written to output metrics file")
        ->default_str("you should probably add a tc name when outputting metrics");

    app.add_option("--hmin", arguments.hmin_factor, "Factor to apply to hmin while remeshing")
        ->default_val(10.0);
    app.add_option("--hmax", arguments.hmax_factor, "Factor to apply to hmin while remeshing")
        ->default_val(0.1);

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    auto metrics = moist::metrics::Metrics(arguments.tc_name);
    moist::Timer timer("InterfaceInserter::Main", metrics);

    moist::utils::geogram::initialize();
    geo::CmdLine::import_arg_group("algo");
    geo::CmdLine::import_arg_group("output");
    geo::CmdLine::set_arg("output:precision", "16");
    geo::CmdLine::set_arg("sys:use_doubles", true);

    // what do I need?
    // Overrefinement
        // Initial MeshA/B
        // After Reconnecting
        // After Coarsening
        // After Remeshing
        // After Remeshing w/o coarsening
    // Quality (Aspect Ratio, Dihedral Angles)
        // Initial MeshA/B
        // After Reconnecting
        // After Coarsening
        // After Remeshing
        // After Remeshing w/o coarsening
    // Runtimes
        // Time spent in Lookup based on different factors... start with 1x1 big cell (baseline), keep going exponentially (factor = 2, 4, 16 times sqrt cells)?

    geo::Mesh slice_a, slice_b;
    moist::utils::geogram::load(arguments.input_a, slice_a, 3, false);
    moist::utils::geogram::load(arguments.input_b, slice_b, 3, false);
    const auto plane = moist::AxisAlignedPlane { moist::Axis::Z, arguments.extent, arguments.epsilon };
    // moist::create_interface_mesh(slice_a, slice_b, plane, moist::RemeshingParameters {arguments.hmin_factor, arguments.hmax_factor}, metrics);
    moist::Merger merger(slice_a, slice_b, plane, arguments.grid_size_factor, moist::RemeshingParameters {arguments.hmin_factor, arguments.hmax_factor}, metrics);

#ifndef NDEBUG
    moist::ScopeTimer::Print();
#endif // NDEBUG

    moist::ScopeTimer::WriteMetrics(metrics);
    if (use_metrics->count() > 0)
    {
        metrics->AppendCSV(arguments.metrics);
    }
    return 0;
}
