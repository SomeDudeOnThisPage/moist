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

struct Arguments
{
    std::filesystem::path input_a;
    std::filesystem::path input_b;

    double epsilon;
    double extent;
    std::filesystem::path output_a;
    std::filesystem::path output_b;
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

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    auto metrics = moist::metrics::Metrics("InterfaceInserter");
    moist::Timer timer("InterfaceInserter::Main", metrics);

    moist::utils::geogram::initialize();
    geo::CmdLine::import_arg_group("algo");
    geo::CmdLine::import_arg_group("output");
    geo::CmdLine::set_arg("output:precision", "16");
    geo::CmdLine::set_arg("sys:use_doubles", true);

    // log general run information
    //*metrics << moist::metrics::Metric {"interface::extent", interface.Plane()->extent}
    //         << moist::metrics::Metric {"interface::epsilon", interface.Plane()->epsilon}
    //         << moist::metrics::Metric {"mesh::A", arguments.input_a}
    //         << moist::metrics::Metric {"mesh::B", arguments.input_b}
    //         ;

    moist::metrics::MeshQuality quality_before_a{"A::before"};
    moist::metrics::MeshQuality quality_after_a{"A::after"};
    moist::metrics::MeshQuality quality_interface_before_a{"A::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_a{"A::interface::after"};
    moist::metrics::MeshQuality quality_before_b{"B::before"};
    moist::metrics::MeshQuality quality_after_b{"B::after"};
    moist::metrics::MeshQuality quality_interface_before_b{"B::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_b{"B::interface::after"};

    // moist::MeshSlice slice_a, slice_b;
    geo::Mesh slice_a, slice_b;
    moist::utils::geogram::load(arguments.input_a, slice_a, 3, false);
    moist::utils::geogram::load(arguments.input_b, slice_b, 3, false);
    const auto plane = moist::AxisAlignedPlane { moist::Axis::Z, arguments.extent, arguments.epsilon };
    moist::create_interface_mesh(slice_a, slice_b, plane);

#ifndef NDEBUG
    moist::ScopeTimer::Print();
#endif // NDEBUG
    //metrics->AppendCSV("test.csv");
    return 0;
}
