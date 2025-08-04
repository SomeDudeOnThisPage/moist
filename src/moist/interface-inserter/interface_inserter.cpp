#include "interface_inserter.hpp"

#include <unordered_map>
#include <unordered_set>

#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/mesh_quality.inl"
#include "moist/core/timer.hpp"
#include "moist/core/utils.hpp"

#include "exact_mesh.hpp"
#include "geometry_exact.inl"

#ifdef MOIST_OPTION_EXACT_PREDICATES

static std::size_t find_vertex_by_point(const moist::exact::Point& p, const moist::ExactMesh& mesh)
{
    for (std::size_t v = 0; v < mesh.NbPoints(); v++)
    {
        if (mesh.Point(v) == p)
        {
            return v;
        }
    }

    return moist::exact::NO_VERTEX;
}

static void merge_edge_meshes(const moist::ExactMesh& a, const moist::ExactMesh& b, moist::ExactMesh& constraints)
{
    // Ensure proper indexing and no duplication of points.
    std::unordered_map<std::size_t, std::size_t> map_b_to_constraints;

    // Add all vertices and edges from mesh a to the constraints, and create a map from point indices shared between meshes to the constraint mesh point indices.
    for (const auto edge : a.Edges())
    {
        for (const auto v : edge._points)
        {
            const auto point = a.Point(v);
            const auto cv = constraints.Add(point);
            if (const auto bv = find_vertex_by_point(point, b) != moist::exact::NO_VERTEX)
            {
                map_b_to_constraints.insert({bv, cv});
            }
        }

        constraints.Add(moist::exact::Edge(edge));
    }

    for (const auto edge : b.Edges())
    {
        const auto v0 = map_b_to_constraints.contains(edge[0]) ? map_b_to_constraints.at(edge[0]) : constraints.Add(b.Point(edge[0]));
        const auto v1 = map_b_to_constraints.contains(edge[1]) ? map_b_to_constraints.at(edge[1]) : constraints.Add(b.Point(edge[1]));

        // Check if our edge intersects with another edge from mesh a.
        // If so, split both edges at the intersection point, creating four new edges.
        for (auto other_edge : constraints.Edges())
        {
            if (other_edge._deleted)
            {
                continue;
            }

            // Two edges cannot intersect themselves in a single point, if they share at least one point.
            const auto ov0 = other_edge[0];
            const auto ov1 = other_edge[1];
            if (v0 == ov0 || v0 == ov1 || v1 == ov0 || v1 == ov1)
            {
                continue;
            }

            const auto intersection = moist::geometry::exact::intersection(
                moist::exact::EdgePoints(constraints.Point(v0), constraints.Point(v1)),
                moist::exact::EdgePoints(constraints.Point(ov0), constraints.Point(ov1))
            );

            if (!intersection.has_value())
            {
                constraints.Add(moist::exact::Edge(v0, v1));
                continue;
            }
            else
            {
                other_edge._deleted = true;
                const auto v_split = constraints.Add(intersection.value());
                constraints.Add(moist::exact::Edge(v0, v_split));
                constraints.Add(moist::exact::Edge(v1, v_split));
                constraints.Add(moist::exact::Edge(ov0, v_split));
                constraints.Add(moist::exact::Edge(ov0, v_split));
            }
        }
    }
}
#else // ifndef MOIST_OPTION_EXACT_PREDICATES

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
#endif // MOIST_OPTION_EXACT_PREDICATES

enum class MeshQualityStage
{
    BEFORE, AFTER, INTERFACE_BEFORE, INTERFACE_AFTER
};

#ifdef MOIST_OPTION_EXACT_PREDICATES
static std::size_t empty(moist::ExactMesh& mesh) { return mesh.NbPoints() == 0; }
static void reset(moist::ExactMesh& mesh) { mesh.ResetMesh(); }

static void to_exact_edge_mesh(geo::Mesh& mesh, moist::ExactMesh& exact)
{
    std::unordered_map<geo::index_t, std::size_t> vertices;
    for (const auto v : mesh.vertices)
    {
        vertices[v] = exact.Add(moist::exact::Point(mesh.vertices.point(v)));
    }

    const auto edges = moist::geometry::collect_edges(mesh);
    for (const auto edge : edges)
    {
        exact.Add(moist::exact::Edge(vertices.at(edge.v0), vertices.at(edge.v1)));
    }
}


#else // ifndeef MOIST_OPTION_EXACT_PREDICATES
static std::size_t empty(geo::Mesh& mesh) { return mesh.vertices.nb() == 0; }
static void reset(moist::ExactMesh& mesh) { mesh.clear(false, false); }
#endif // MOIST_OPTION_EXACT_PREDICATES


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

#ifdef MOIST_OPTION_EXACT_PREDICATES
    moist::ExactMesh exact_interface;
    to_exact_edge_mesh(*interface.Triangulation(), exact_interface);
    exact_interface.DebugMesh("exact_interface.msh");
    a.InsertEdges(*interface.Triangulation(), *interface.Plane());
    b.InsertEdges(*interface.Triangulation(), *interface.Plane());
#else // ifndef MOIST_OPTION_EXACT_PREDICATES
    a.InsertEdges(*interface.Triangulation(), *interface.Plane());
    b.InsertEdges(*interface.Triangulation(), *interface.Plane());
#endif // MOIST_OPTION_EXACT_PREDICATES

    std::size_t iterations = 0;
#ifdef MOIST_OPTION_EXACT_PREDICATES
    moist::ExactMesh fixed_a, fixed_b, fixed_ab;
    a.GetFixedGeometry(fixed_a);
    b.GetFixedGeometry(fixed_b);
    if (fixed_a.NbPoints() > 0 || fixed_b.NbPoints() > 0)
    {
#else // ifndef MOIST_OPTION_EXACT_PREDICATES
    geo::Mesh fixed_a, fixed_b, fixed_ab;
    a.GetFixedGeometry(fixed_a);
    b.GetFixedGeometry(fixed_b);
    if (nde_a.vertices.nb() > 0 || nde_b.vertices.nb() > 0)
    {
#endif // MOIST_OPTION_EXACT_PREDICATES
        do
        {
            iterations++;
            merge_edge_meshes(fixed_a, fixed_b, fixed_ab);

    #ifndef NDEBUG
        #ifdef MOIST_OPTION_EXACT_PREDICATES
            fixed_ab.DebugMesh("merged_steiner_edges_" + to_string(iterations) + ".obj");
        #else
            moist::utils::geogram::save("merged_steiner_edges_" + to_string(iterations) + ".obj", fixed_ab);
        #endif //  MOIST_OPTION_EXACT_PREDICATES
    #endif // NDEBUG

            //a.InsertEdges(fixed_ab, *interface.Plane());
            //b.InsertEdges(fixed_ab, *interface.Plane());

            reset(fixed_a);
            reset(fixed_b);
            reset(fixed_ab);

            a.GetFixedGeometry(fixed_a);
            b.GetFixedGeometry(fixed_b);
        }
        while (!empty(fixed_a) || !empty(fixed_b));
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
