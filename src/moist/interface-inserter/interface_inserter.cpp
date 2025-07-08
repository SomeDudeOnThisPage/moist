#include "interface_inserter.hpp"

#include <unordered_map>
#include <unordered_set>

#include <geogram/mesh/mesh_repair.h>

#include "moist/core/mesh_quality.inl"
#include "moist/core/timer.hpp"
#include "moist/core/utils.hpp"

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

void moist::insert_constraints(moist::MeshSlice& a, moist::MeshSlice& b, moist::Interface& interface, moist::metrics::Metrics_ptr metrics)
{
    moist::Timer timer("InterfaceInserter::InsertInterfaceConstraints", metrics);

    // Quality metrics for meshes at different points in time during the algorithm.
    moist::metrics::MeshQuality quality_before_a{"A::before"};
    moist::metrics::MeshQuality quality_after_a{"A::after"};
    moist::metrics::MeshQuality quality_interface_before_a{"A::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_a{"A::interface::after"};

    moist::metrics::MeshQuality quality_before_b{"B::before"};
    moist::metrics::MeshQuality quality_after_b{"B::after"};
    moist::metrics::MeshQuality quality_interface_before_b{"B::interface::before"};
    moist::metrics::MeshQuality quality_interface_after_b{"B::interface::after"};


    moist::mesh_quality::compute(quality_before_a, a);
    moist::mesh_quality::compute(quality_before_b, b);
    moist::mesh_quality::compute(quality_interface_before_a, a, interface);
    moist::mesh_quality::compute(quality_interface_before_b, b, interface);

    geo::Mesh nde_a, nde_b, nde_ab;
    a.InsertEdges(*interface.Triangulation(), *interface.Plane());
    b.InsertEdges(*interface.Triangulation(), *interface.Plane());

    geo::mesh_repair(a);
    geo::mesh_repair(b);

    a.GetFixedGeometry(nde_a);
    b.GetFixedGeometry(nde_b);

    size_t iterations = 0;
    if (nde_a.vertices.nb() > 0 || nde_b.vertices.nb() > 0)
    {
        do
        {
            iterations++;
            OOC_DEBUG("inserting " << (nde_a.edges.nb() + nde_b.edges.nb()) << " non-decimatable edges into both meshes (iteration #" << iterations << ")");

            merge_edge_meshes(nde_a, nde_b, nde_ab);
        #ifndef NDEBUG
            moist::utils::geogram::save("merged_steiner_edges.obj", nde_ab);
        #endif // NDEBUG

            a.InsertEdges(nde_ab, *interface.Plane());
            b.InsertEdges(nde_ab, *interface.Plane());

            geo::mesh_repair(a);
            geo::mesh_repair(b);

            nde_a.clear(false, false);
            nde_b.clear(false, false);
            nde_ab.clear(false, false);

            a.GetFixedGeometry(nde_a);
            b.GetFixedGeometry(nde_b);
        }
        while (nde_a.vertices.nb() > 0 || nde_b.vertices.nb() > 0);
    }

    moist::mesh_quality::compute(quality_after_a, a);
    moist::mesh_quality::compute(quality_after_b, b);
    moist::mesh_quality::compute(quality_interface_after_a, a, interface);
    moist::mesh_quality::compute(quality_interface_after_b, b, interface);

    *metrics << moist::metrics::Metric {"nde_insertion_iterations", iterations};
    *metrics << quality_before_a
             << quality_interface_before_a
             << quality_after_a
             << quality_interface_after_a
             << quality_before_b
             << quality_interface_before_b
             << quality_after_b
             << quality_interface_after_b
             ;
}
