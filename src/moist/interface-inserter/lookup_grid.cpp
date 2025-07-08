#include "lookup_grid.hpp"

moist::LookupGrid::LookupGrid(const vec2 grid_size, const vec2 cell_size) : _grid_size(grid_size), _cell_size(cell_size)
{
}

static void get_interface_aabb()
{

}

void moist::LookupGrid::GetGridCells(const GEO::Box2d &aabb, std::vector<g_index> &cells)
{
    vec2 min = {aabb.xy_min[0], aabb.xy_min[1]};
    vec2 max = {aabb.xy_max[0], aabb.xy_max[1]};

    uint32_t start_x = std::max(static_cast<uint32_t>(0), static_cast<uint32_t>(std::floor(min.x / _cell_size.x)));
    uint32_t end_x = std::min(static_cast<uint32_t>(_grid_size.x) - 1, static_cast<uint32_t>(std::floor(max.x / _cell_size.x)));

    uint32_t start_y = std::max(static_cast<uint32_t>(0), static_cast<uint32_t>(std::floor(min.y / _cell_size.y)));
    uint32_t end_y = std::min(static_cast<uint32_t>(_grid_size.y) - 1, static_cast<uint32_t>(std::floor(max.y / _cell_size.y)));

    // Collect overlapping cells...
    for (size_t y = start_y; y <= end_y; y++)
    {
        for (size_t x = start_x; x <= end_x; x++)
        {
            g_index index = y * static_cast<uint32_t>(_grid_size.x) + x;
            cells.push_back(index);
        }
    }
}

void moist::LookupGrid::Insert(const g_index c, const GEO::Mesh &mesh)
{
    mesh.cells.points(c);
}

void moist::LookupGrid::Erase(const g_index c, const GEO::Mesh &mesh)
{
}

std::vector<g_index> &moist::LookupGrid::GetGridCell(const g_index &index)
{
    //return _cells.at(index);
}
