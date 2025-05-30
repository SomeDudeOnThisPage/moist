#include <filesystem>
#include <algorithm>
#include <regex>
#include <map>

#include <CLI/CLI.hpp>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>

#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/predicates.inl"
#include "moist/core/metrics.hpp"
#include "moist/core/mesh_quality.inl"

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
    auto metrics = moist::metrics::TimeMetrics("InterfaceInserter");
    moist::Timer timer("InterfaceInserter::Main", metrics);

    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.input, "Mesh (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-i, --interface", arguments.interface, "Interface mesh (.geogram) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-o, --output", arguments.output, "Output Triangulation (.msh) file")
        ->required();

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geo::initialize();

    auto interface = moist::Interface(arguments.interface);

    moist::metrics::MeshQuality metrics_before{0};
    moist::metrics::MeshQuality metrics_after{0};

    moist::MeshSlice slice;
    moist::utils::geo::load(arguments.input, slice, 3, false);
    moist::mesh_quality::compute(metrics_before, slice);
    {
        moist::Timer _scope_timer("MeshSlice::InsertInterface", metrics);
        slice.InsertInterface(interface, metrics);
        geogram::mesh_repair(slice, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);
    }
    moist::mesh_quality::compute(metrics_after, slice);

    OOC_DEBUG("[METRICS] before insertion: " << metrics_before);
    OOC_DEBUG("[METRICS] after insertion: " << metrics_after);

    moist::utils::geo::save(arguments.output.replace_extension(".msh"), slice);
    moist::utils::geo::save(arguments.interface.replace_extension(".geogram"), *interface.Triangulation());

    return 0;
}
