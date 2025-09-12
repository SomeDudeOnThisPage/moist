#include "merger.hpp"

#include <cmath>

#include <tetgen.h>

#include <geogram/mesh/mesh_repair.h>

#include "moist/core/predicates.inl"
#include "moist/core/utils.hpp"
#include "moist/core/timer.hpp"

#include "arrangement.hpp"
#include "exact_types.hpp"
#include "local_operations.hpp"
#include "new_predicates.inl"
#include "geometry_exact.inl"

static geo::Mesh dbg_triangles(3);

#define GRID_FACTOR 3.0
//#define GRID_BY_EXTENT
#define GRID_BY_NUM_CELLS
#define MIST_EVAL_GRID
// #define MIST_EVAL_QUALITY

// multiplied by shortest edge in the mesh BEFORE remeshing...
// gradually relax bounds
constexpr std::array<double, 4> EVAL_HMIN =
{
    1.0, 0.8, 0.6, 0.4
};

// multiplied by longest edge in the mesh BEFORE remeshing...
constexpr std::array<double, 4> EVAL_HMAX =
{
    1.0, 1.2, 1.4, 1.6
};

// hausd is unchanged we just use the suggested value from mmg docs

void moist::Merger::EvaluateRemeshing(moist::metrics::Metrics_ptr metrics)
{
    const geo::Box3d aabb = create_mesh_bbox3d_exact(_crown);
    const double hd_factor = std::cbrt((aabb.xyz_max[0] - aabb.xyz_min[0]) * (aabb.xyz_max[1] - aabb.xyz_min[1]) * (aabb.xyz_max[2] - aabb.xyz_min[2])) * 0.01;
    *metrics << moist::metrics::Metric("remeshing::hd_factor", hd_factor);

    for (std::size_t i = 0; i < 4; i++)
    {
        // for (std::size_t j = 0; j < 4; j++) // too many useless data points, as relaxing hmax without relaxing hmin usually does not result in different bounds on elements
        // {
            auto remeshing_metrics = moist::metrics::Metrics(metrics->name + "_remeshing_" + std::to_string((int) i) + "-" + std::to_string((int) i));
            const double hmin = EVAL_HMIN[i];
            const double hmax = EVAL_HMAX[i];
            moist::metrics::MeshQuality metrics_before("before");
            moist::metrics::MeshQuality metrics_after("after");

            MOIST_INFO("performing evaluation remeshing pass w/ hmin=" << hmin << ", hmax=" << hmax);
            *remeshing_metrics << moist::metrics::Metric("pass", i);
            *remeshing_metrics << moist::metrics::Metric("hmin", hmin);
            *remeshing_metrics << moist::metrics::Metric("hmax", hmax);
            *remeshing_metrics << moist::metrics::Metric("hmin_abs", _min_edge_length * hmin);
            *remeshing_metrics << moist::metrics::Metric("hmax_abs", _max_edge_length * hmax);
            _crown.ComputeMetrics(metrics_before);

            MMG5_pMesh mmg_mesh = NULL;
            MMG5_pSol mmg_solution = NULL;
            MMG3D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &mmg_mesh,
                    MMG5_ARG_ppMet, &mmg_solution,
                MMG5_ARG_end);

            moist::mmg3d::transform(_crown, mmg_mesh, mmg_solution);

            moist::ExactMesh remeshing_output;
            moist::mmg3d::set_solution(mmg_mesh, mmg_solution, _min_edge_length * hmin, _max_edge_length * hmax, hd_factor);
            moist::mmg3d::remesh(mmg_mesh, mmg_solution);
            moist::mmg3d::transform(mmg_mesh, mmg_solution, remeshing_output);
            remeshing_output.ComputeMetrics(metrics_after);
            *remeshing_metrics << metrics_before;
            *remeshing_metrics << metrics_after;
            remeshing_metrics->AppendCSV(metrics->name + "_metrics.csv");
        // }
    }
}


static void snap_histograms(moist::ExactMesh& mesh, moist::metrics::Metrics_ptr metrics, const std::string& suffix)
{
#ifdef MIST_EVAL_QUALITY
    // save histogram to metrics for qualities
    moist::metrics::Histogram hg_aspect_ratio(10);
    moist::metrics::Histogram hg_mean_ratio(10);
    moist::metrics::Histogram hg_mmg3d_quality(10);

    mesh.ComputeHistogram(hg_aspect_ratio, moist::metrics::QualityMetricType::ASPECT_RATIO);
    mesh.ComputeHistogram(hg_mean_ratio, moist::metrics::QualityMetricType::MEAN_RATIO);
    mesh.ComputeHistogram(hg_mmg3d_quality, moist::metrics::QualityMetricType::MMG3D_QUALITY);

    *metrics << moist::metrics::Metric("histogram::aspect_ratio::" + suffix, hg_aspect_ratio);
    *metrics << moist::metrics::Metric("histogram::mean_ratio" + suffix, hg_mean_ratio);
    *metrics << moist::metrics::Metric("histogram::mmg3d_quality" + suffix, hg_mmg3d_quality);
#endif
}

static void snap_histograms_2(moist::ExactMesh& a, moist::ExactMesh& b, moist::metrics::Metrics_ptr metrics, const std::string& suffix)
{
#ifdef MIST_EVAL_QUALITY
    moist::metrics::Histogram hg_aspect_ratio(10);
    moist::metrics::Histogram hg_mean_ratio(10);
    moist::metrics::Histogram hg_mmg3d_quality(10);
    a.ComputeHistogram(hg_aspect_ratio, moist::metrics::QualityMetricType::ASPECT_RATIO);
    a.ComputeHistogram(hg_mean_ratio, moist::metrics::QualityMetricType::MEAN_RATIO);
    a.ComputeHistogram(hg_mmg3d_quality, moist::metrics::QualityMetricType::MMG3D_QUALITY);
    b.ComputeHistogram(hg_aspect_ratio, moist::metrics::QualityMetricType::ASPECT_RATIO);
    b.ComputeHistogram(hg_mean_ratio, moist::metrics::QualityMetricType::MEAN_RATIO);
    b.ComputeHistogram(hg_mmg3d_quality, moist::metrics::QualityMetricType::MMG3D_QUALITY);
    *metrics << moist::metrics::Metric("histogram::aspect_ratio::" + suffix, hg_aspect_ratio);
    *metrics << moist::metrics::Metric("histogram::mean_ratio" + suffix, hg_mean_ratio);
    *metrics << moist::metrics::Metric("histogram::mmg3d_quality" + suffix, hg_mmg3d_quality);
#endif
}

moist::Merger::Merger(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane, float grid_factor, const moist::RemeshingParameters remeshing, moist::metrics::Metrics_ptr metrics) : _metrics(metrics)
{
    _min_edge_length = std::numeric_limits<double>::max();
    _max_edge_length = -std::numeric_limits<double>::max();
    _grid_factor = grid_factor;

    this->ConstructExactMesh(a, _mesh_a, plane);
    this->ConstructExactMesh(b, _mesh_b, plane);

    snap_histograms_2(_mesh_a, _mesh_b, metrics, "ab");

    *metrics << moist::metrics::Metric("grid_factor", _grid_factor);
    *metrics << moist::metrics::Metric("grid_resolution::a", _mesh_a.Grid().Resolution());
    *metrics << moist::metrics::Metric("grid_resolution::b", _mesh_b.Grid().Resolution());
#ifdef MIST_EVAL_QUALITY
    *metrics << moist::metrics::Metric("remeshing::hmin_factor", remeshing.hmin);
    *metrics << moist::metrics::Metric("remeshing::hmax_factor", remeshing.hmax);
    *metrics << moist::metrics::Metric("remeshing::min_edge_length::before", _min_edge_length);
    *metrics << moist::metrics::Metric("remeshing::max_edge_length::before", _max_edge_length);
#endif

    _crown.ResetGrid((_mesh_a.Grid().Resolution() + _mesh_b.Grid().Resolution()) / 2.0);
    moist::utils::geogram::save("crownless_a.mesh", a);
    moist::utils::geogram::save("crownless_b.mesh", b);

#ifdef MIST_EVAL_QUALITY
    moist::metrics::MeshQuality a_before("a_before");
    moist::metrics::MeshQuality b_before("b_before");
    _mesh_a.ComputeMetrics(a_before);
    _mesh_b.ComputeMetrics(b_before);
    moist::metrics::MeshQuality ab_before = a_before + b_before;
    ab_before.name = "ab_before"; // operator doesnt set name
    *metrics << ab_before;
#endif

    this->InsertAndMapPoints(_mesh_a, _mesh_b);
    this->InsertAndMapPoints(_mesh_b, _mesh_a);

#ifndef NDEBUG
    _mesh_a.DebugMesh("exact_a.mesh");
    _mesh_b.DebugMesh("exact_b.mesh");
#endif // NDEBUG

    this->InsertBToA();
#ifdef MIST_EVAL_QUALITY
    moist::metrics::MeshQuality crown_before_remeshing("crown_before_remeshing");
    _crown.ComputeMetrics(crown_before_remeshing);
    *metrics << crown_before_remeshing;
#endif

#ifndef NDEBUG
    _crown.DebugMesh("crown_before_remeshing.mesh");
    moist::utils::geogram::save("triangles.off", dbg_triangles);
#endif // NDEBUG
    snap_histograms(_crown, metrics, "crown::before_remeshing");

    MOIST_INFO("coarsening...");
    tetgenio coarsen_mesh;
    tetgenio output_coarsen_mesh;
    output_coarsen_mesh.initialize();
    moist::tetgen::transform(_crown, coarsen_mesh);
    moist::tetgen::coarsen(coarsen_mesh, output_coarsen_mesh);
    moist::ExactMesh coarsened_crown;
    moist::tetgen::transform(output_coarsen_mesh, coarsened_crown, plane);

#ifdef MIST_EVAL_QUALITY
    moist::metrics::MeshQuality crown_after_coarsening("crown_after_coarsening");
    coarsened_crown.ComputeMetrics(crown_after_coarsening);
    *metrics << crown_after_coarsening;
    coarsened_crown.DebugMesh("crown_after_coarsening.mesh");
    snap_histograms(coarsened_crown, metrics, "crown_after_coarsening");
#endif

#ifdef MIST_EVAL_QUALITY
    this->EvaluateRemeshing(metrics);
#endif

#ifndef MIST_EVAL_GRID
    const geo::Box3d aabb = create_mesh_bbox3d_exact(_crown);
    const double hd_factor = std::cbrt((aabb.xyz_max[0] - aabb.xyz_min[0]) * (aabb.xyz_max[1] - aabb.xyz_min[1]) * (aabb.xyz_max[2] - aabb.xyz_min[2])) * 0.01;

    MMG5_pMesh mmg_mesh = NULL;
    MMG5_pSol mmg_solution = NULL;
    MMG3D_Init_mesh(MMG5_ARG_start,
            MMG5_ARG_ppMesh, &mmg_mesh,
            MMG5_ARG_ppMet, &mmg_solution,
        MMG5_ARG_end);

    moist::mmg3d::transform(_crown, mmg_mesh, mmg_solution);

    moist::ExactMesh remeshing_output;
    moist::mmg3d::set_solution(mmg_mesh, mmg_solution, _min_edge_length, 2.0 * _max_edge_length, hd_factor);
    moist::mmg3d::remesh(mmg_mesh, mmg_solution);
    moist::mmg3d::transform(mmg_mesh, mmg_solution, remeshing_output);
    remeshing_output.DebugMesh("crown_after_remeshing.msh");

/*#ifndef NDEBUG
    if (MMG3D_saveMesh(mmg_mesh, "output_mmg3d_before_remesning.mesh") != 1)
    {
        OOC_ERROR("cannot debug mmg_mesh");
    }
#endif // NDEBUG


    moist::metrics::MeshQuality crown_after_remeshing("crown_after_remeshing");
    _crown.ComputeMetrics(crown_after_remeshing);
    *metrics << crown_after_remeshing;
    _crown.DebugMesh("crown_after_remeshing.mesh");

    snap_histograms(_crown, metrics, "crown::after_remeshing");*/
#endif MIST_EVAL_GRID
}

/**
 * @brief Preprocess rounding to avoid extremely small triangles in overlaps during processing (i.e. edge lengths of e-20 or less)
 *
 * @param value
 * @return double
 */
static double round16(double value)
{
    double rounded = std::round(value * 1e8) / 1e8;
    if (rounded == 0.0)
    {
        rounded = 0.0; // don't ask
    }
    return rounded;
}

void moist::Merger::ConstructExactMesh(geo::Mesh& mesh, moist::ExactMesh& exact_mesh, const moist::AxisAlignedPlane& plane)
{
    moist::ScopeTimer TIMER("moist::Merger::ConstructExactMesh");

    std::unordered_map<geo::index_t, std::size_t> added_vertices;
    geo::vector<geo::index_t> to_delete(mesh.cells.nb());
    for (const geo::index_t c : mesh.cells)
    {
        if (!moist::predicates::cell_on_plane(c, mesh, plane))
        {
            continue;
        }

        std::array<std::size_t, 4> vertices;
        for (geo::index_t lv = 0; lv < 4; lv++)
        {
            const geo::index_t v = mesh.cells.vertex(c, lv);
            if (added_vertices.contains(v))
            {
                vertices[lv] = added_vertices[v];
                continue;
            }
            else
            {
                // here we must keep _other of the point always on NO_VERTEX...
                auto& point = mesh.vertices.point(v);
                bool interface = moist::predicates::point_on_plane(point, plane);
                if (interface)
                {
                    point[2] = plane.extent;
                }

                vec3 point_rounded = vec3(round16(point[0]), round16(point[1]), round16(point[2]));
                vertices[lv] = exact_mesh.Add(moist::exact::Point(point_rounded, interface));
                added_vertices[v] = vertices[lv];
            }
        }

        if (geo::mesh_cell_volume(mesh, c) != 0.0)
        {
            exact_mesh.Add(moist::exact::Cell(
                vertices[0], vertices[1], vertices[2], vertices[3],
                moist::geometry::exact::get_cell_type(vertices, exact_mesh)
            ), true);
        }
        else
        {
            OOC_DEBUG("skipped cell " << c << " during interface plane mesh construction due to degeneracy in eps...");
        }

        to_delete[c] = true; // delete all these cells from the original mesh...

        // Update edge lengths...
        const auto lengths = moist::geometry::longest_shortest_edge(c, mesh);
        _min_edge_length = (lengths[0] < _min_edge_length) ? lengths[0] : _min_edge_length;
        _max_edge_length = (lengths[1] > _max_edge_length) ? lengths[0] : _max_edge_length;
    }
    mesh.cells.delete_elements(to_delete);

    // edge length does not make the most sense, since it must be based relative to the overall length of our bbox diag.
    // dumb hard coded value... use nb of cells on the interface maybe... times some factor controlled so we can evaluate performance based on it...
#ifdef GRID_BY_EXTENT
    const auto bbox = create_mesh_bbox3d_exact(exact_mesh);
    double x = bbox.xyz_max[0] - bbox.xyz_min[0];
    double y = bbox.xyz_max[1] - bbox.xyz_min[1];
    exact_mesh.ResetGrid(_grid_factor * (x + y) / 2.0);
    OOC_DEBUG("grid resolution: " << _grid_factor * (x + y) / 2.0);
#else
    exact_mesh.ResetGrid(_grid_factor * std::sqrt(exact_mesh.Cells().size()));
    OOC_DEBUG("grid resolution: " << _grid_factor * std::sqrt(exact_mesh.Cells().size()));
#endif
    OOC_DEBUG("edge lengths -> max: " << _max_edge_length << ", min: " << _min_edge_length);
}

void moist::Merger::InsertAndMapPoints(const moist::ExactMesh& from, moist::ExactMesh& to)
{
    for (std::size_t v = 0; v < from.Points().size(); v++)
    {
        const auto& point = from.Point(v);
        if (point._other != moist::exact::NO_VERTEX)
        {
            continue;
        }

        this->InsertPoint(v, from, to);
    }
}

static bool point_in_tet_on_edge(const moist::predicates::PointInTet& pit)
{
    return pit == moist::predicates::PointInTet::EDGE01
           || pit == moist::predicates::PointInTet::EDGE02
           || pit == moist::predicates::PointInTet::EDGE03
           || pit == moist::predicates::PointInTet::EDGE12
           || pit == moist::predicates::PointInTet::EDGE13
           || pit == moist::predicates::PointInTet::EDGE23;
}

// This is done in the "old" meshes, not in the one being currently constructed, because tets with an edge on the interface may be cut
// multiple times.
bool moist::Merger::InsertPointIntoEdges(const moist::exact::Point& point, moist::ExactMesh& mesh)
{
    moist::ScopeTimer TIMER("moist::Merger::InsertPointIntoEdges");

    const auto& grid_cells = mesh.Grid().GetCells(moist::create_point_box2d_exact(point));
    std::size_t v_1to2 = moist::exact::NO_VERTEX;

    bool found = false;
    std::unordered_set<std::size_t> inserted;
    for (const auto grid_cell : grid_cells)
    {
        const auto cells = mesh.Grid().GetMeshCells(grid_cell);
        for (const auto c : cells)
        {
            const auto& cell = mesh.Cell(c);
            if (cell._deleted || inserted.contains(c))
            {
                continue;
            }

            std::size_t nb_interface = 0;
            for (std::size_t lv = 0; lv < 4; lv++)
            {
                if (mesh.Point(cell._points[lv])._interface)
                {
                    nb_interface++;
                }
            }

            if (nb_interface != 2)
            {
                continue;
            }

            const auto location = moist::new_predicates::point_in_tet_exact(mesh, c, point, true);
            if (point_in_tet_on_edge(location))
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = mesh.Add(moist::exact::Point(point));
                    mesh.Point(v_1to2)._interface = true;
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, mesh);
                found = true;
                inserted.insert(c);
            }
            else if (location == moist::predicates::PointInTet::FACET)
            {
                // const std::size_t v = mesh.Add(point);
                // moist::operation::exact::InsertVertexOnCellBoundaryFacet(c, v, mesh);
            }
        }
    }

    return found;
}

void moist::Merger::InsertPoint(const std::size_t v, const moist::ExactMesh& from, moist::ExactMesh& to)
{
    // check where to insert the point, and update to accordingly. also map the point!
    const auto point = from.Point(v);
    const auto& grid_cells = to.Grid().GetCells(moist::create_point_box2d_exact(point));
    std::size_t v_1to2 = moist::exact::NO_VERTEX;

    for (const auto grid_cell : grid_cells)
    {
        const auto cells = to.Grid().GetMeshCells(grid_cell);
        for (const auto c : cells)
        {
            const auto& cell = to.Cell(c);
            if (cell._deleted)
            {
                continue;
            }

            const auto location = moist::new_predicates::point_in_tet_exact(to, c, point, true);
            if (location == moist::predicates::PointInTet::FACET)
            {
                const std::size_t v_new = to.Add(moist::exact::Point(point));
                moist::operation::exact::InsertVertexOnCellBoundaryFacet(c, v_new, to);
                to.Point(v_new)._other = v;
                to.Point(v_new)._interface = true;
            }
            else if (point_in_tet_on_edge(location))
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = to.Add(moist::exact::Point(point));
                    to.Point(v_1to2)._other = v;
                    to.Point(v_1to2)._interface = true;
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, to);
            }
        }
    }
}

static std::unordered_set<std::size_t> touched_a;
static geo::Mesh dbg_full(3);
void moist::Merger::InsertBToA()
{
    moist::ScopeTimer TIMER("moist::Merger::InsertBToA");

    const auto size_target = _mesh_b.Cells().size();
    for (std::size_t c = 0; c < _mesh_b.Cells().size(); c++)
    {
        if (c % 50 == 0)
        {
            OOC_DEBUG("merged " << c << " out of " << size_target << " cells...");
        }

        const auto cell = _mesh_b.Cell(c);
        if (cell._deleted)
        {
            continue;
        }

        if (!moist::geometry::exact::cell_has_interface_facet(cell, _mesh_b))
        {
            continue;
        }

        this->InsertCellBToA(c);
    }

    for (std::size_t c = 0; c < _mesh_a.Cells().size(); c++)
    {
        const auto& cell = _mesh_a.Cell(c);
        if (cell._deleted || touched_a.contains(c))
        {
            continue;
        }

        const auto p0 = _mesh_a.Point(cell._points[0]);
        const auto p1 = _mesh_a.Point(cell._points[1]);
        const auto p2 = _mesh_a.Point(cell._points[2]);
        const auto p3 = _mesh_a.Point(cell._points[3]);

        _crown.Add(moist::exact::Cell(
            _crown.Add(moist::exact::Point(p0)),
            _crown.Add(moist::exact::Point(p1)),
            _crown.Add(moist::exact::Point(p2)),
            _crown.Add(moist::exact::Point(p3)),
            cell._type
        ), true);
    }

    for (std::size_t c = 0; c < _mesh_b.Cells().size(); c++)
    {
        const auto& cell = _mesh_b.Cell(c);
        if (cell._deleted)
        {
            continue;
        }

        const auto p0 = _mesh_b.Point(cell._points[0]);
        const auto p1 = _mesh_b.Point(cell._points[1]);
        const auto p2 = _mesh_b.Point(cell._points[2]);
        const auto p3 = _mesh_b.Point(cell._points[3]);

        _crown.Add(moist::exact::Cell(
            _crown.Add(moist::exact::Point(p0)),
            _crown.Add(moist::exact::Point(p1)),
            _crown.Add(moist::exact::Point(p2)),
            _crown.Add(moist::exact::Point(p3)),
            cell._type
        ), true);
    }

    // resolve remaining intersections between type II cells...
    /*_crown.ResetGrid(GRID_FACTOR * std::sqrt(_crown.NbCells()));
    std::unordered_set<std::size_t> touched;
    MOIST_INFO("size of crown " << _crown.Cells().size());
    const std::size_t maxc = _crown.Cells().size();
    for (std::size_t c = 0; c < maxc; c++)
    {
        const auto& cell = _crown.Cell(c);
        if (cell._deleted) continue;

        const auto edges = moist::geometry::exact::get_interface_edges(c, _crown);
        if (edges.size() == 1)
        {
            const moist::exact::EdgePoints e0 = {_crown.Point(cell[edges[0][0]]), _crown.Point(cell[edges[0][1]])};
            const auto bbox = create_edge_box2d_exact(e0);
            const auto grid_cells = _crown.Grid().GetCells(bbox);

            for (const auto grid_cell : grid_cells)
            {
                for (const auto c_other : _crown.Grid().GetMeshCells(grid_cell))
                {
                    if (c_other == c) continue;
                    auto cell_other = _crown.Cell(c_other);
                    if (cell_other._deleted || touched.contains(c_other))
                    {
                        continue;
                    }

                    const auto edges_other = moist::geometry::exact::get_interface_edges(c_other, _crown);
                    for (const auto edge : edges_other)
                    {
                        const moist::exact::EdgePoints e1 = {_crown.Point(cell_other[edge[0]]), _crown.Point(cell_other[edge[1]])};
                        if (!e1.p0._interface || !e1.p1._interface) continue;
                        const auto intersection = moist::geometry::exact::intersection(e0, e1);
                        if (!intersection.has_value()) continue;

                        const auto point = intersection.value();
                        const auto p = _crown.Add(moist::exact::Point(point));
                        _crown.Point(p)._interface = true;
                        //this->InsertPoint(p, _crown, _crown);
                        cell_other._deleted = true;
                        touched.insert(c_other);
                        break;
                    }
                }
            }
        }
    }
    MOIST_INFO("size of crown " << _crown.Cells().size());*/
}

bool moist::Merger::InsertCellBToA(const std::size_t cb)
{
    moist::ScopeTimer TIMER("moist::Merger::InsertCellBToA");

    // since grid cells can overlap in which cells they contain, we need to track into which cells we already inserted, so as to not
    // insert a cell twice!
    std::unordered_set<std::size_t> inserted;

    auto cell = _mesh_b.Cell(cb);
    const auto grid_cells = _mesh_a.Grid().GetCells(moist::create_interface_cell_bbox2d_exact(cb, _mesh_b));

    for (const auto grid_cell : grid_cells)
    {
        for (const auto ca : _mesh_a.Grid().GetMeshCells(grid_cell))
        {
            const auto cell = _mesh_a.Cell(ca);
            if (cell._deleted || inserted.contains(ca))
            {
                continue;
            }

            if (!moist::geometry::exact::cell_has_interface_facet(cell, _mesh_a))
            {
                continue;
            }

            if (this->InsertCellIntoCell(ca, cb))
            {
                inserted.emplace(ca);
            }
        }
    }

    if (!inserted.empty())
    {
        _mesh_b.Cell(cb)._deleted = true;
    }
    return !inserted.empty();
}

static moist::exact::Triangle get_interface_triangle(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    std::vector<moist::exact::Point> points;
    for (const auto v : cell._points)
    {
        const auto point = mesh.Point(v);
        if (point._interface)
        {
            points.push_back(point);
        }
    }

    if (points.size() != 3)
    {
        OOC_ERROR("invalid interface cell triangle configuration");
    }

    return moist::exact::Triangle(points.at(0), points.at(1), points.at(2));
}

static std::size_t get_opposite_vertex(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    for (const auto v : cell._points)
    {
        const auto point = mesh.Point(v);
        if (!point._interface)
        {
            return v;
        }
    }

    OOC_ERROR("invalid opposite vertex");
    return -1U;
}

static moist::exact::VertexCorrespondence get_correspondence(const moist::exact::Triangle& triangle)
{
    auto correspondence = moist::exact::VertexCorrespondence::AB;

    for (std::size_t i = 0; i < 3; i++)
    {
        if (triangle._points.at(i)._correspondence != moist::exact::VertexCorrespondence::AB)
        {
            if (correspondence != moist::exact::VertexCorrespondence::AB && triangle._points.at(i)._correspondence != correspondence)
            {
                OOC_ERROR("failed to determine vertex correspondence");
            }

            correspondence = triangle._points.at(i)._correspondence; // technically we can return here already
        }
    }

    return correspondence;
}

static bool triangle_has_point(const moist::exact::Triangle& triangle, const moist::exact::Kernel::Point_3& point)
{
    for (std::size_t lv = 0; lv < 3; lv++)
    {
        if (triangle._t.vertex(lv) == point)
        {
            return true;
        }
    }

    return false;
};

// remove all cells of crown that intersect the new cell...
static void remove_intersecting_III_cells(const std::size_t c, moist::ExactMesh& mesh)
{
    const auto cell = mesh.Cell(c);
    const auto t0 = get_interface_triangle(cell, mesh);

    const auto grid_cells = mesh.Grid().GetCells(moist::create_cell_bbox2d_exact(c, mesh));
    for (const auto grid_cell : grid_cells)
    {
        for (auto lc : mesh.Grid().GetMeshCells(grid_cell))
        {
            auto& cell_other = mesh.Cell(lc);
            if (lc == c || cell_other._deleted || cell_other._type != moist::exact::CellType::III)
            {
                continue;
            }

            const auto t1 = get_interface_triangle(mesh.Cell(lc), mesh);
            if (moist::exact::arrangeable(t0, t1))
            {
                mesh.Cell(lc)._deleted = true;
            }
        }
    }
}

bool moist::Merger::InsertCellIntoCell(const std::size_t ca, const std::size_t cb)
{
    moist::ScopeTimer TIMER("moist::Merger::InsertCellIntoCell");

    auto& cell_a = _mesh_a.Cell(ca);
    auto& cell_b = _mesh_b.Cell(cb);

    const auto t0 = get_interface_triangle(cell_a, _mesh_a);
    const auto t1 = get_interface_triangle(cell_b, _mesh_b);
    const std::size_t v_opposite_a = get_opposite_vertex(cell_a, _mesh_a);
    const std::size_t v_opposite_b = get_opposite_vertex(cell_b, _mesh_b);

    if (!moist::exact::arrangeable(t0, t1))
    {
        return false;
    }

    const auto triangulation = moist::exact::arrange(t0, t1);

    if (triangulation.empty())
    {
        return false;
    }

    for (const auto triangle : triangulation)
    {
        const auto p_opposite_a = _mesh_a.Point(v_opposite_a);
        const auto p_opposite_b = _mesh_b.Point(v_opposite_b);

        // check which mesh(es) the triangle belongs to:
        // 1. if all points are AB, connect in both
        // 2. if some or no points are AB, connect to whatever the other points are
        // 3. if we have A and B in one point, fail, but this should not happen unless the CGAL routine messes up
        const auto correspondence = get_correspondence(triangle);
        if (correspondence == moist::exact::VertexCorrespondence::AB)
        {
            for (std::size_t lv = 0; lv < 3; lv++)
            {
                bool in_a = true;
                bool in_b = true;
                if (!triangle_has_point(t0, triangle._t.vertex(lv)))
                {
                    in_a = this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_a);
                }

                if (!triangle_has_point(t1, triangle._t.vertex(lv)))
                {
                    in_b = this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_b);
                }
            }

            auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
            auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
            auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
            points_triangle_0._interface = true;
            points_triangle_1._interface = true;
            points_triangle_2._interface = true;
            std::size_t v0 = _crown.Add(points_triangle_0);
            std::size_t v1 = _crown.Add(points_triangle_1);
            std::size_t v2 = _crown.Add(points_triangle_2);

            _crown.Add(moist::exact::Cell(
                v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_a)),
                moist::exact::CellType::III
            ), true);

            _crown.Add(moist::exact::Cell(
                v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_b)),
                moist::exact::CellType::III
            ), true);
        }
        else if (correspondence == moist::exact::VertexCorrespondence::A)
        {
            // the cell is only present in mesh b in this arragement -> check if it intersects with any cell of mesh a, if not, add it...
            const auto aabb = moist::create_triangle_bbox2d_exact(triangle);
            const auto grid_cells = _mesh_b.Grid().GetCells(aabb);
            bool do_add = true;
            for (const auto grid_cell : grid_cells)
            {
                for (const auto cb_check : _mesh_b.Grid().GetMeshCells(grid_cell))
                {
                    if (_mesh_b.Cell(cb_check)._type != moist::exact::CellType::III)
                    {
                        continue;
                    }

                    const auto t0 = get_interface_triangle(_mesh_b.Cell(cb_check), _mesh_b);
                    if (moist::exact::arrangeable(t0, triangle))
                    {
                        // if we have at least one triangular overlap, we must not add it to the crown.
                        do_add = false;
                        goto add_cell_b_into_crown;
                    }
                }
            }
add_cell_b_into_crown:
            if (do_add)
            {
                auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
                auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
                auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
                points_triangle_0._interface = true;
                points_triangle_1._interface = true;
                points_triangle_2._interface = true;
                std::size_t v0 = _crown.Add(points_triangle_0);
                std::size_t v1 = _crown.Add(points_triangle_1);
                std::size_t v2 = _crown.Add(points_triangle_2);
                std::size_t new_cell_crown = _crown.Add(moist::exact::Cell(
                    v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_a)),
                    moist::exact::CellType::III
                ), true);
            }
        }
else if (correspondence == moist::exact::VertexCorrespondence::B)
        {
            // the cell is only present in mesh b in this arragement -> check if it intersects with any cell of mesh a, if not, add it...
            const auto aabb = moist::create_triangle_bbox2d_exact(triangle);
            const auto grid_cells = _mesh_a.Grid().GetCells(aabb);
            bool do_add = true;
            for (const auto grid_cell : grid_cells)
            {
                for (const auto cb_check : _mesh_a.Grid().GetMeshCells(grid_cell))
                {
                    if (_mesh_a.Cell(cb_check)._type != moist::exact::CellType::III)
                    {
                        continue;
                    }

                    const auto t0 = get_interface_triangle(_mesh_a.Cell(cb_check), _mesh_a);
                    if (moist::exact::arrangeable(t0, triangle))
                    {
                        // if we have at least one triangular overlap, we must not add it to the crown.
                        do_add = false;
                        goto add_cell_a_into_crown;
                    }
                }
            }
add_cell_a_into_crown:
            if (do_add)
            {
                auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
                auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
                auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
                points_triangle_0._interface = true;
                points_triangle_1._interface = true;
                points_triangle_2._interface = true;
                std::size_t v0 = _crown.Add(points_triangle_0);
                std::size_t v1 = _crown.Add(points_triangle_1);
                std::size_t v2 = _crown.Add(points_triangle_2);
                std::size_t new_cell_crown = _crown.Add(moist::exact::Cell(
                    v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_b)),
                    moist::exact::CellType::III
                ), true);
            }
        }
    }
    touched_a.emplace(ca);
    return true;
}
