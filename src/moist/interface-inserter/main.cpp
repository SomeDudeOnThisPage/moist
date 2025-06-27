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
    std::filesystem::path input_a;
    std::filesystem::path input_b;

    std::filesystem::path interface;
    std::filesystem::path output_a;
    std::filesystem::path output_b;

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

    app.add_option("-a, --mesh-a", arguments.input_a, "Mesh A (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-b, --mesh-b", arguments.input_b, "Mesh A (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-i, --interface", arguments.interface, "Interface mesh (.geogram) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("--output-a", arguments.output_a, "Output Triangulation A (.msh) file")
        ->required();
    app.add_option("--output-b", arguments.output_b, "Output Triangulation B (.msh) file")
        ->required();

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geo::initialize();

    auto interface = moist::Interface(arguments.interface);

    moist::metrics::MeshQuality metrics_before{0};
    moist::metrics::MeshQuality metrics_after{0};
    moist::SteinerPoints steiner_points;

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geo::load(arguments.input_a, slice_a, 3, false);
    moist::utils::geo::load(arguments.input_b, slice_b, 3, false);

    {
        moist::Timer _scope_timer("A::MeshSlice::InsertInterface", metrics);
        auto sps = slice_a.InsertInterface(interface, metrics);
        steiner_points.insert(sps.begin(), sps.end());
    }

    {
        moist::Timer _scope_timer("B::MeshSlice::InsertInterface", metrics);
        auto sps = slice_b.InsertInterface(interface, metrics);
        steiner_points.insert(sps.begin(), sps.end());
    }

    geogram::mesh_repair(slice_a);
    geogram::mesh_repair(slice_b);

    if (!steiner_points.empty())
    {
        slice_a.FlushTetrahedra(true);
        slice_b.FlushTetrahedra(true);
        OOC_DEBUG("inserting " << steiner_points.size() << " steiner points into both meshes...");
        for (const vec3& steiner_point : steiner_points)
        {
            slice_a.InsertVertex(steiner_point, *interface.Plane());
            slice_b.InsertVertex(steiner_point, *interface.Plane());
        }

    }

    slice_a.FlushTetrahedra(false);
    slice_b.FlushTetrahedra(false);

    OOC_DEBUG("[METRICS] before insertion: " << metrics_before);
    OOC_DEBUG("[METRICS] after insertion: " << metrics_after);

    moist::utils::geo::save(arguments.output_a.replace_extension(".mesh"), slice_a);
    moist::utils::geo::save(arguments.output_b.replace_extension(".mesh"), slice_b);



    moist::utils::geo::save(arguments.interface.replace_extension(".geogram"), *interface.Triangulation());

    return 0;
}
