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

    auto metrics = moist::metrics::Metrics("InterfaceInserter");
    moist::Timer timer("InterfaceInserter::Main", metrics);

    moist::utils::geo::initialize();

    auto interface = moist::Interface(arguments.interface);

    // log general run information
    *metrics << moist::metrics::Metric {"interface::extent", interface.Plane()->extent}
             << moist::metrics::Metric {"interface::epsilon", interface.Plane()->epsilon}
             << moist::metrics::Metric {"mesh::A", arguments.input_a}
             << moist::metrics::Metric {"mesh::B", arguments.input_b}
             ;

    moist::metrics::MeshQuality quality_before_a{"A::before"};
    moist::metrics::MeshQuality quality_after_a{"A::after"};
    moist::metrics::MeshQuality quality_interface_before_a{"A::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_a{"A::interface::after"};
    moist::metrics::MeshQuality quality_before_b{"B::before"};
    moist::metrics::MeshQuality quality_after_b{"B::after"};
    moist::metrics::MeshQuality quality_interface_before_b{"B::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_b{"B::interface::after"};

    moist::SteinerPoints steiner_points;

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geo::load(arguments.input_a, slice_a, 3, false);
    moist::utils::geo::load(arguments.input_b, slice_b, 3, false);

    moist::mesh_quality::compute(quality_before_a, slice_a);
    moist::mesh_quality::compute(quality_before_b, slice_b);
    moist::mesh_quality::compute(quality_interface_before_a, slice_a, interface);
    moist::mesh_quality::compute(quality_interface_before_b, slice_b, interface);

    {
        // log mesh a metrics...
        *metrics << moist::metrics::Metric {"mesh::A::before::nb_vertices", slice_a.vertices.nb()};
        moist::Timer _scope_timer("A::MeshSlice::InsertInterface", metrics);
        auto sps = slice_a.InsertInterface(interface, metrics);
        steiner_points.insert(sps.begin(), sps.end());
    }

    {
        *metrics << moist::metrics::Metric {"mesh::B::before::nb_vertices", slice_a.vertices.nb()};
        moist::Timer _scope_timer("B::MeshSlice::InsertInterface", metrics);
        auto sps = slice_b.InsertInterface(interface, metrics);
        steiner_points.insert(sps.begin(), sps.end());
    }

    geogram::mesh_repair(slice_a);
    geogram::mesh_repair(slice_b);

    *metrics << moist::metrics::Metric {"mesh::nb_steiner_vertices", steiner_points.size()}
             << moist::metrics::Metric {"A::after::nb_vertices", slice_a.vertices.nb()}
             << moist::metrics::Metric {"B::before::nb_vertices", slice_b.vertices.nb()};

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

    moist::mesh_quality::compute(quality_after_a, slice_a);
    moist::mesh_quality::compute(quality_after_b, slice_b);
    moist::mesh_quality::compute(quality_interface_after_a, slice_a, interface);
    moist::mesh_quality::compute(quality_interface_after_b, slice_b, interface);

    *metrics << moist::metrics::Metric {"A::steiner::nb_vertices", slice_a.vertices.nb()}
             << moist::metrics::Metric {"B::steiner::nb_vertices", slice_b.vertices.nb()};

    *metrics << quality_before_a
             << quality_interface_before_a
             << quality_after_a
             << quality_interface_after_a
             << quality_before_b
             << quality_interface_before_b
             << quality_after_b
             << quality_interface_after_b
             ;

    moist::utils::geo::save(arguments.output_a.replace_extension(".mesh"), slice_a);
    moist::utils::geo::save(arguments.output_b.replace_extension(".mesh"), slice_b);

    moist::utils::geo::save(arguments.interface.replace_extension(".geogram"), *interface.Triangulation());

    metrics->AppendCSV("test.csv");
    return 0;
}
