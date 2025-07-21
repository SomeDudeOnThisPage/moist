#ifndef MOIST_CORE_GEOMETRY_EXACT_INL_
#define MOIST_CORE_GEOMETRY_EXACT_INL_

#include <array>

#include "moist/core/defines.hpp"

#include "exact_mesh.hpp"

namespace moist::geometry::exact
{
    PURE INLINE std::array<std::size_t, 3> other(const std::size_t& c, const std::size_t& v_opposite, const moist::ExactMesh& mesh)
    {
        const auto& cell = mesh.Cell(c);
        for (std::size_t lv = 0; lv < 4; lv++)
        {
            const auto& v = cell._points[lv];
            if (v_opposite == v)
            {
                return std::array
                {
                    cell._points[(lv + 1) % 4],
                    cell._points[(lv + 2) % 4],
                    cell._points[(lv + 3) % 4]
                };
            }
        }

        return std::array
        {
            std::size_t(-1U), std::size_t(-1U), std::size_t(-1U)
        };
    }
}

#endif // MOIST_CORE_GEOMETRY_EXACT_INL_
