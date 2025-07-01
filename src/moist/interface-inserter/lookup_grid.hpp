#ifndef MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP
#define MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP

#include <vector>
#include <ranges>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

#include "moist/core/defines.hpp"

namespace moist
{

    class LookupGrid
    {
    public:
        LookupGrid(const vec2 grid_size, const vec2 cell_size);

        /**
         * @brief Get the Grid Cells overlapping the given aabb.
         *
         * @param aabb
         * @param cells
         */
        void GetGridCells(const geogram::Box2d &aabb, std::vector<g_index>& cells);

        void Insert(const g_index c, const geogram::Mesh& mesh);
        void Erase(const g_index c, const geogram::Mesh& mesh);

        std::vector<g_index>& GetGridCell(const g_index& index);
    private:
        vec2 _grid_size;
        vec2 _cell_size;

        std::vector<std::vector<g_index>> _cells;
    };
}

#endif
