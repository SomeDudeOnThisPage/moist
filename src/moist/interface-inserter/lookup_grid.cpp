#include "lookup_grid.hpp"

#include <geogram/mesh/mesh_AABB.h>

#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"
#include "moist/core/timer.hpp"

#include "exact_mesh.hpp"

void moist::LookupGridExact::InsertCell(const std::size_t c, const moist::ExactMesh& mesh)
{
    auto covered_cells = this->GetCells(moist::create_cell_bbox2d_exact(c, mesh));

    for (const auto& cell : covered_cells)
    {
        _grid[cell].insert(c);
    }
}

std::vector<moist::LookupGridExact::GridCell> moist::LookupGridExact::GetCells(const geo::Box2d& aabb) const
{
    std::vector<GridCell> result;

    const GridCell min_cell
    {
        std::clamp(static_cast<int>((aabb.xy_min[0] - 0.1 - _min_bounds.x) / _cell_size.x), 0, static_cast<int>(_resolution) - 1),
        std::clamp(static_cast<int>((aabb.xy_min[1] - 0.1 - _min_bounds.y) / _cell_size.y), 0, static_cast<int>(_resolution) - 1)
    };

    const GridCell max_cell
    {
        std::clamp(static_cast<int>((aabb.xy_max[0] + 0.1 - _min_bounds.x) / _cell_size.x), 0, static_cast<int>(_resolution) - 1),
        std::clamp(static_cast<int>((aabb.xy_max[1] + 0.1 - _min_bounds.y) / _cell_size.y), 0, static_cast<int>(_resolution) - 1)
    };

    for (int i = min_cell.first; i <= max_cell.first; i++)
    {
        for (int j = min_cell.second; j <= max_cell.second; j++)
        {
            result.emplace_back(i, j);
        }
    }

    return result;
}

moist::LookupGridExact::MeshCells& moist::LookupGridExact::GetMeshCells(const moist::LookupGridExact::GridCell& grid_cell)
{
    return _grid.at(grid_cell);
}

std::optional<moist::exact::index_t> moist::LookupPointGrid::Get(moist::exact::Point point)
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::LookupPointGrid::Get");
#endif // NDEBUG

    const moist::LookupPointGrid::Coordinates key = this->GetCoordinates(point);
    auto it = _grid.find(key);
    if (it == _grid.end())
    {
        return std::nullopt;
    }

    for (const auto& existing : it->second)
    {
        if (existing == point)
        {
            return existing.Index(); // this point must contain a reference to its' own index in the mesh in question!
            // On second thougt it wouldve been better to also reference the mesh in here and just have this grid as part of exactmesh for implicit indexing...
        }
    }

    return std::nullopt;
}

void moist::LookupPointGrid::Add(moist::exact::Point point, const moist::exact::index_t index)
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::LookupPointGrid::Add");
#endif // NDEBUG

    const moist::LookupPointGrid::Coordinates key = this->GetCoordinates(point);
    auto& cell = _grid[key];

    for (const auto& existing : cell)
    {
        if (existing == point)
        {
            return;
        }
    }

    cell.push_back(moist::exact::IndexedPoint(point, index));
}

moist::LookupPointGrid::Coordinates moist::LookupPointGrid::GetCoordinates(const moist::exact::Point& point)
{
    double gx_ft = std::floor(point.x() / _resolution);
    double gy_ft = std::floor(point.y() / _resolution);

    long long gx = static_cast<long long>(gx_ft);
    long long gy = static_cast<long long>(gy_ft);

    return { gx, gy };
}
