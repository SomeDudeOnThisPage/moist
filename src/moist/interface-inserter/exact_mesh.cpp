#include "exact_mesh.hpp"

#include <tetrahedron.hpp>
#include <triangle.hpp>
#include <cmath>

#include "geometry_exact.inl"
#include "remeshing.hpp"

#ifndef NDEBUG
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_repair.h>
#include "moist/core/utils.hpp"
#endif // NDEBUG

extern "C"
{
    #include "mmg/mmg3d/libmmg3d.h"
}

static void get_data(double* data, const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    #pragma unroll 4
    for (l_index lv = 0; lv < 4; lv++)
    {
        const auto p = mesh.Point(cell[lv]);
        data[3 * lv] = p.x();
        data[3 * lv + 1] = p.y();
        data[3 * lv + 2] = p.z();
    }
}

static double cell_aspect_ratio(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    double data[12];
    get_data(data, cell, mesh);

    // Sphere insphere;
    // Sphere circumsphere;
    //
    // tetrahedron_insphere(data, insphere.radius, insphere.center.data());
    // tetrahedron_circumsphere(data, circumsphere.radius, circumsphere.center.data());
    //
    // return circumsphere.radius / (3.0 * insphere.radius);

    return tetrahedron::tetrahedron_quality1(data);
}

static double cell_mean_ratio(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    double data[12];
    get_data(data, cell, mesh);
    return tetrahedron::tetrahedron_quality3(data);
}

static double cell_quality_mmg3d(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
}

// returns average, lower bound, higher bound
static std::array<double, 3> aspect_ratio(const moist::ExactMesh& mesh)
{
    double aspect_ratio = 0.0f;
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    std::size_t nb_cells = 0;
    for (const auto cell : mesh.Cells())
    {
        if (cell._deleted)
        {
            continue;
        }

        const double l_aspect_ratio = cell_aspect_ratio(cell, mesh);
        if (l_aspect_ratio < min)
        {
            min = l_aspect_ratio;
        }

        if (l_aspect_ratio > max)
        {
            max = l_aspect_ratio;
        }

        aspect_ratio += l_aspect_ratio;
        nb_cells++;
    }

    return { aspect_ratio / static_cast<double>(nb_cells), min, max };
}

// returns average, lower bound, higher bound
// tetrahedron actually uses the inverse mean ratio (lolol)
static std::array<double, 3> mean_ratio(const moist::ExactMesh& mesh)
{
    double mean_ratio = 0.0f;
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    std::size_t nb_cells = 0;
    for (const auto cell : mesh.Cells())
    {
        if (cell._deleted)
        {
            continue;
        }

        const double l_mean_ratio = cell_mean_ratio(cell, mesh);
        if (l_mean_ratio < min)
        {
            min = l_mean_ratio;
        }

        if (l_mean_ratio > max)
        {
            max = l_mean_ratio;
        }

        mean_ratio += l_mean_ratio;
        nb_cells++;
    }

    return { mean_ratio / static_cast<double>(nb_cells), min, max };
}

static std::array<double, 3> mmg_quality(moist::ExactMesh& mesh)
{
    MMG5_pMesh mmg_mesh = NULL;
    MMG5_pSol mmg_solution = NULL;
    MMG3D_Init_mesh(MMG5_ARG_start,
            MMG5_ARG_ppMesh, &mmg_mesh,
            MMG5_ARG_ppMet, &mmg_solution,
        MMG5_ARG_end);

    moist::mmg3d::transform(mesh, mmg_mesh, mmg_solution);
    int nb_points;
    int nb_cells;
    MMG3D_Get_meshSize(mmg_mesh, &nb_points, &nb_cells, NULL, NULL, NULL, NULL);

    double mmg_quality = 0.0f;
    double min = 100000.0;
    double max = -100000.0;
    for (std::size_t c = 1; c < nb_cells; c++)
    {
        const double mmg_quality_l = MMG3D_Get_tetrahedronQuality(mmg_mesh, nullptr, c);
        min = mmg_quality_l < min ? mmg_quality_l : min;
        max = mmg_quality_l > max ? mmg_quality_l : max;
        mmg_quality += mmg_quality_l;
    }
    return { mmg_quality / static_cast<double>(nb_cells), min, max };
}

moist::ExactMesh::ExactMesh() : _grid(moist::LookupGridExact()), _point_grid(moist::LookupPointGrid(10.0))
{
}

void moist::ExactMesh::ComputeMetrics(moist::metrics::MeshQuality& metrics)
{
    const auto aspect_ratios = aspect_ratio(*this);
    const auto mean_ratios = mean_ratio(*this);
    const auto mmg = mmg_quality(*this);

    std::size_t nv = 0;
    std::size_t nc = 0;
    for (const auto point : _points)
    {
        nv += (point._deleted) ? 0 : 1;
    }

    for (const auto cell : _cells)
    {
        nc += (cell._deleted) ? 0 : 1;
    }

    metrics.nb_vertices = nv;
    metrics.nb_cells = nc;
    metrics.aspect_ratio = aspect_ratios[0];
    metrics.aspect_ratio_lb = aspect_ratios[1];
    metrics.aspect_ratio_ub = aspect_ratios[2];
    metrics.mean_ratio = mean_ratios[0];
    metrics.mean_ratio_lb = mean_ratios[1];
    metrics.mean_ratio_ub = mean_ratios[2];
    metrics.mmg = mmg[0];
    metrics.mmg_lb = mmg[1];
    metrics.mmg_ub = mmg[2];
}

void moist::ExactMesh::ComputeHistogram(moist::metrics::Histogram& histogram, const moist::metrics::QualityMetricType& type, const std::size_t nb_bins)
{
    if (nb_bins == -1U && histogram.size() == 0)
    {
        throw std::out_of_range("invalid histogram size");
    }

    const int n = static_cast<int>((nb_bins == -1U) ? histogram.size() : nb_bins);
    // histogram.clear();
    if (n != histogram.size())
    {
        histogram.resize(n);
    }

    auto get_bin = [&n](const double q) -> std::size_t
    {
        return std::min(int(q * n), n - 1);
    };

    if (type == moist::metrics::QualityMetricType::MMG3D_QUALITY)
    {
        MMG5_pMesh mmg_mesh = NULL;
        MMG5_pSol mmg_solution = NULL;
        MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh, &mmg_mesh,
                MMG5_ARG_ppMet, &mmg_solution,
            MMG5_ARG_end);

        moist::mmg3d::transform(*this, mmg_mesh, mmg_solution);
        int nb_points;
        int nb_cells;
        MMG3D_Get_meshSize(mmg_mesh, &nb_points, &nb_cells, NULL, NULL, NULL, NULL);

        for (std::size_t c = 1; c < nb_cells; c++)
        {
            const std::size_t bin = get_bin(MMG3D_Get_tetrahedronQuality(mmg_mesh, nullptr, c));
            histogram[bin]++;
        }
    }
    else
    {
        for (const auto cell : this->Cells())
        {
            if (cell._deleted)
            {
                continue;
            }

            double metric = 0.0;
            switch (type)
            {
                case moist::metrics::QualityMetricType::ASPECT_RATIO: metric = cell_aspect_ratio(cell, *this); break;
                case moist::metrics::QualityMetricType::MEAN_RATIO: metric = cell_mean_ratio(cell, *this); break;
            }
            if (std::isnan(metric))
            {
                continue;
            }
            const std::size_t bin = get_bin(metric);
            histogram[bin]++;
        }
    }

}

void moist::ExactMesh::ResetGrid(const double resolution)
{
    _grid._grid.clear();
    _grid._resolution = resolution; // TODO: Calculate!
    _grid._min_bounds = vec2(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    _grid._max_bounds = vec2(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    for (const auto& point : this->Points())
    {
        _grid._min_bounds.x = std::min(_grid._min_bounds.x, point.x());
        _grid._min_bounds.y = std::min(_grid._min_bounds.y, point.y());
        _grid._max_bounds.x = std::max(_grid._max_bounds.x, point.x());
        _grid._max_bounds.y = std::max(_grid._max_bounds.y, point.y());
    }

    _grid._cell_size = (_grid._max_bounds - _grid._min_bounds) / _grid._resolution;

    for (std::size_t c = 0; c < this->Cells().size(); c++)
    {
        _grid.InsertCell(c, *this);
    }
}

void moist::ExactMesh::ResetMesh()
{
    _points.clear();
    _edges.clear();
    _cells.clear();
    _pid = 0;
    _tid = 0;
    _eid = 0;
    _grid._grid.clear();
}

std::size_t moist::ExactMesh::Add(const moist::exact::Point p, const bool initialization)
{
    if (!initialization)
    {
        // check if the point already exists...
        const auto existing = _point_grid.Get(p);
        if (existing != std::nullopt)
        {
            return existing.value();
        }
    }

    std::size_t index = _pid;
    _points.push_back(p);
    _point_grid.Add(p, _pid);
    _pid++;
    return index;
}

std::size_t moist::ExactMesh::Add(const moist::exact::Cell c, const bool initialization)
{
    std::size_t index = _tid;
    _cells.push_back(c);
    _tid++;
    if (!initialization)
    {
        _grid.InsertCell(index, *this);
    }
    return index;
}

std::size_t moist::ExactMesh::Add(const moist::exact::Edge e)
{
    std::size_t index = _eid;
    if (!_edges.contains(e))
    {
        _edges.insert(e);
        _eid++;
    }
    return index; // kinda meaningless lol
}

std::size_t moist::ExactMesh::Add(const moist::exact::Facet f)
{
    _facets.push_back(f);
    return _facets.size();
}

void moist::ExactMesh::DeletePoint(const std::size_t v)
{
    _points.at(v)._deleted = true;
}

void moist::ExactMesh::DeleteCell(const std::size_t c)
{
    _cells.at(c)._deleted = true;
}

void moist::ExactMesh::FlushDeletedElements()
{
    std::size_t del = 0;
    for (const auto cell : _cells)
    {
        if (cell._deleted)
        {
            del++;
        }
    }
    OOC_DEBUG(del);

    const std::size_t size_before = _cells.size();
    _cells.erase(
        std::remove_if(_cells.begin(), _cells.end(), [](const auto& cell) { return cell._deleted; }), _cells.end());
    _cells.shrink_to_fit();
    OOC_DEBUG("shrinking cells from " << size_before << " to " << _cells.size());
    // this->ResetGrid(this->_grid->_resolution);
}

moist::exact::Point &moist::ExactMesh::Point(const std::size_t& index)
{
    if (_points.size() < index)
    {
        throw std::runtime_error("Point index does not exist: " + index);
    }

    return _points.at(index);
}

const moist::exact::Point & moist::ExactMesh::Point(const std::size_t & index) const
{
    if (_points.size() < index)
    {
        throw std::runtime_error("Point index does not exist: " + index);
    }

    return _points.at(index);
}

moist::exact::Cell& moist::ExactMesh::Cell(const std::size_t& index)
{
    if (_cells.size() < index)
    {
        throw std::runtime_error("Cell index does not exist: " + index);
    }

    return _cells.at(index);
}

const moist::exact::Cell& moist::ExactMesh::Cell(const std::size_t& index) const
{
    if (_cells.size() < index)
    {
        throw std::runtime_error("Cell index does not exist: " + index);
    }

    return _cells.at(index);
}

#ifndef NDEBUG
void moist::ExactMesh::DebugMesh(const std::filesystem::path& file)
{
#ifndef NDEBUG
    std::set<std::array<double, 3>> points;
#endif // NDEBUG

    geo::Mesh mesh(3);
    for (const auto& p : this->_points)
    {
        const auto point = geo::vec3(CGAL::to_double(p._p.x()), CGAL::to_double(p._p.y()), CGAL::to_double(p._p.z()));
    #ifndef NDEBUG
        if (points.contains({ point.x, point.y, point.z }))
        {
            OOC_WARNING("degenerate vertices in R");
        }
        points.insert({ point.x, point.y, point.z });
    #endif // NDEBUG
        mesh.vertices.create_vertex(point);
    }

    geo::vector<geo::index_t> to_delete(this->NbCells());
    for (const auto c : this->_cells)
    {
        if (c._deleted || moist::geometry::exact::is_degenerate(c))
        {
            continue;
        }

        const auto t = mesh.cells.create_tet(
            geo::index_t(c._points[0]),
            geo::index_t(c._points[1]),
            geo::index_t(c._points[2]),
            geo::index_t(c._points[3])
        );

        if (std::fabs(geo::mesh_cell_volume(mesh, t)) == 0.0)
        {
            OOC_WARNING("degenerate tet " << t);
            to_delete[t] = true;
        }
    }
    mesh.cells.delete_elements(to_delete);

    for (const auto e : this->_edges)
    {
        if (e._deleted)
        {
            continue;
        }

        mesh.edges.create_edge(
            geo::index_t(e._points[0]),
            geo::index_t(e._points[1])
        );
    }

    for (const auto f : this->_facets)
    {
        mesh.facets.create_triangle(
            geo::index_t(f[0]),
            geo::index_t(f[1]),
            geo::index_t(f[2])
        );
    }

    // geo::mesh_repair(mesh, geo::MeshRepairMode::MESH_REPAIR_COLOCATE);
    to_delete.clear();
    to_delete.resize(mesh.cells.nb());
    double smallest_vol = 10000.0;
    for (const geo::index_t c : mesh.cells)
    {
        to_delete[c] = std::fabs(geo::mesh_cell_volume(mesh, c)) == 0.0 || moist::geometry::has_duplicate_vertex(c, mesh);
        const auto vol = std::fabs(geo::mesh_cell_volume(mesh, c));
        if (vol != 0.0 && vol < smallest_vol)
        {
            smallest_vol = vol;
        }
    }
    mesh.cells.delete_elements(to_delete);
    OOC_DEBUG("smallest volume: " << smallest_vol);

    moist::utils::geogram::save(file, mesh);
}
#endif // NDEBUG
