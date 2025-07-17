#include "lookup_grid.hpp"

#include <geogram/mesh/mesh_AABB.h>

#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

moist::LookupGrid::LookupGrid(const geo::Mesh& mesh) : _mesh(mesh)
{

}

void moist::LookupGrid::Initialize(const double resolution)
{
    _grid.clear();
    _resolution = resolution;
    _min_bounds = vec2(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    _max_bounds = vec2(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    for (geo::index_t v = 0; v < _mesh.vertices.nb(); v++)
    {
        const double* p = _mesh.vertices.point_ptr(v);
        _min_bounds.x = std::min(_min_bounds.x, p[0]);
        _min_bounds.y = std::min(_min_bounds.y, p[1]);
        _max_bounds.x = std::max(_max_bounds.x, p[0]);
        _max_bounds.y = std::max(_max_bounds.y, p[1]);
    }

    _cell_size = (_max_bounds - _min_bounds) / _resolution;

    auto v_interface = geo::Attribute<bool>(_mesh.vertices.attributes(), "v_interface");
    for (geo::index_t c = 0; c < _mesh.cells.nb(); c++)
    {
        if (!moist::predicates::is_interface_cell(c, _mesh))
        {
            continue;
        }

        // Construct AABB for each cell on the interface...
        this->InsertCell(c);
    }
}

std::vector<moist::LookupGrid::GridCell> moist::LookupGrid::GetCells(const geo::Box2d& aabb) const
{
    std::vector<GridCell> result;

    auto to_cell = [&](const vec2& p)
    {
        int ix = std::clamp(int((p.x - _min_bounds.x) / _cell_size.x), 0, int(_resolution) - 1);
        int iy = std::clamp(int((p.y - _min_bounds.y) / _cell_size.y), 0, int(_resolution) - 1);
        return std::make_pair(ix, iy);
    };

    const GridCell min_cell {
        std::clamp(int((aabb.xy_min[0] - _min_bounds.x) / _cell_size.x), 0, int(_resolution) - 1),
        std::clamp(int((aabb.xy_min[1] - _min_bounds.y) / _cell_size.y), 0, int(_resolution) - 1)
    };
    const GridCell max_cell {
        std::clamp(int((aabb.xy_max[0] - _min_bounds.x) / _cell_size.x), 0, int(_resolution) - 1),
        std::clamp(int((aabb.xy_max[1] - _min_bounds.y) / _cell_size.y), 0, int(_resolution) - 1)
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

void moist::LookupGrid::InsertCell(const geo::index_t c)
{
    auto aabb = create_interface_cell_bbox2d(c, _mesh, geo::Attribute<bool>(_mesh.vertices.attributes(), "v_interface"));
    auto covered_cells = this->GetCells(aabb);

    for (const auto& cell : covered_cells)
    {
        _grid[cell].insert(c);
    }
}
