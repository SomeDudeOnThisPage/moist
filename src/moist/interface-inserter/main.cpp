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

static geo::index_t find_vertex_by_point(const vec3& p, const geo::Mesh& mesh)
{
    for (const geo::index_t v : mesh.vertices)
    {
        if (p == mesh.vertices.point(v))
        {
            return v;
        }
    }

    return geo::NO_VERTEX;
}

static void merge_edge_meshes(const geo::Mesh& a, const geo::Mesh& b, geo::Mesh& constraints)
{
    std::unordered_map<geo::index_t, geo::index_t> v_b_to_constraints;
    for (const geo::index_t e : a.edges)
    {
        const geo::index_t v0 = a.edges.vertex(e, 0);
        const geo::index_t v1 = a.edges.vertex(e, 1);
        const geo::vec3 p0 = a.vertices.point(v0);
        const geo::vec3 p1 = a.vertices.point(v1);

        const geo::index_t nv0 = constraints.vertices.create_vertex(p0);
        const geo::index_t nv1 = constraints.vertices.create_vertex(p1);

        const geo::index_t b_v0 = find_vertex_by_point(p0, b);
        const geo::index_t b_v1 = find_vertex_by_point(p1, b);

        if (b_v0 != geo::NO_VERTEX)
        {
            v_b_to_constraints.insert({b_v0, nv0});
        }

        if (b_v1 != geo::NO_VERTEX)
        {
            v_b_to_constraints.insert({b_v1, nv1});
        }

        constraints.edges.create_edge(nv0, nv1);
    }

    std::unordered_set<geo::index_t> to_delete;
    for (const geo::index_t e : b.edges)
    {
        const geo::index_t v0 = b.edges.vertex(e, 0);
        const geo::index_t v1 = b.edges.vertex(e, 1);
        const geo::vec3 p0 = b.vertices.point(v0);
        const geo::vec3 p1 = b.vertices.point(v1);

        const geo::index_t nv0 = (v_b_to_constraints.contains(v0)) ? v_b_to_constraints.at(v0) : constraints.vertices.create_vertex(p0);
        const geo::index_t nv1 = (v_b_to_constraints.contains(v1)) ? v_b_to_constraints.at(v1) : constraints.vertices.create_vertex(p1);

        // Check for intersections...
        bool intersecting = false;
        for (const geo::index_t e_other : constraints.edges)
        {
            if (to_delete.contains(e_other))
            {
                continue;
            }

            const geo::index_t v0_other = constraints.edges.vertex(e_other, 0);
            const geo::index_t v1_other = constraints.edges.vertex(e_other, 0);

            // Edges cannot intersect if they share one or more vertices (they can only be equal or intersecting in an endpoint, which is fine)...
            if (nv0 == v0_other || nv0 == v1_other || nv1 == v0_other || nv1 == v1_other)
            {
                continue;
            }

            const geo::vec3 p0_other = constraints.vertices.point(v0_other);
            const geo::vec3 p1_other = constraints.vertices.point(v1_other);

            const auto intersection = moist::predicates::xy::get_line_intersection(
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1),
                    reinterpret_cast<const vec2&>(p0_other),
                    reinterpret_cast<const vec2&>(p1_other)
            );


            if (!intersection.has_value())
            {
                continue;
            }
            else
            {
                intersecting = true;
                to_delete.insert(e_other);

                const geo::index_t v_split = constraints.vertices.create_vertex(intersection.value());
                constraints.edges.create_edge(v0_other, v_split);
                constraints.edges.create_edge(v1_other, v_split);
                constraints.edges.create_edge(nv0, v_split);
                constraints.edges.create_edge(nv1, v_split);
            }
        }

        if (!intersecting)
        {
            constraints.edges.create_edge(nv0, nv1);
        }
    }

    if (!to_delete.empty())
    {
        geo::vector<geo::index_t> del(constraints.edges.nb());
        for (const geo::index_t to_d : to_delete)
        {
            del[to_d] = true;
        }

        constraints.edges.delete_elements(del);
    }
}

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

    moist::utils::geogram::initialize();

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

    moist::MeshSlice slice_a, slice_b;
    moist::utils::geogram::load(arguments.input_a, slice_a, 3, false);
    moist::utils::geogram::load(arguments.input_b, slice_b, 3, false);

    moist::insert_constraints(slice_a, slice_b, interface, metrics);

    // First insertion...
    /*{
        // log mesh a metrics...
        *metrics << moist::metrics::Metric {"mesh::A::before::nb_vertices", slice_a.vertices.nb()};
        moist::Timer _scope_timer("A::MeshSlice::InsertInterface", metrics);
        //auto sps = slice_a.InsertInterface(interface, metrics);
        //steiner_points.insert(sps.begin(), sps.end());
        slice_a.InsertEdges(*interface.Triangulation(), *interface.Plane());
    }

    {
        *metrics << moist::metrics::Metric {"mesh::B::before::nb_vertices", slice_a.vertices.nb()};
        moist::Timer _scope_timer("B::MeshSlice::InsertInterface", metrics);
        //auto sps = slice_b.InsertInterface(interface, metrics);
        //steiner_points.insert(sps.begin(), sps.end());
        slice_b.InsertEdges(*interface.Triangulation(), *interface.Plane());
    }

    geo::mesh_repair(slice_a);
    geo::mesh_repair(slice_b);

    geo::Mesh steiner_a;
    geo::Mesh steiner_b;
    slice_a.GetFixedGeometry(steiner_a);
    slice_b.GetFixedGeometry(steiner_b);

    geo::Mesh edge_constraints;
    merge_edge_meshes(steiner_a, steiner_b, edge_constraints);
    moist::utils::geogram::save("merged_steiner_edges.obj", edge_constraints);

    slice_a.InsertEdges(edge_constraints, *interface.Plane());
    slice_b.InsertEdges(edge_constraints, *interface.Plane());

    *metrics << moist::metrics::Metric {"A::after::nb_vertices", slice_a.vertices.nb()}
             << moist::metrics::Metric {"B::before::nb_vertices", slice_b.vertices.nb()};

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
*/

    moist::utils::geogram::save(arguments.output_a.replace_extension(".mesh"), slice_a);
    moist::utils::geogram::save(arguments.output_b.replace_extension(".mesh"), slice_b);

    moist::utils::geogram::save(arguments.interface.replace_extension(".geogram"), *interface.Triangulation());

    metrics->AppendCSV("test.csv");
    return 0;
}
