#ifndef MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP
#define MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP

#include <vector>
#include <ranges>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

namespace moist
{
    class LookupGrid
    {
    public:
        LookupGrid(const vec2 grid_size, const vec2 cell_size);

        std::ranges::join_view<std::vector<std::span<g_index>>>> Test();
    private:
        vec2 _grid_size;
        vec2 _cell_size;

        std::vector<std::vector<g_index>> _cells;
    };
}

#endif
