#ifndef MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP
#define MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP

#include <ranges>
#include <unordered_set>
#include <vector>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/geometry.inl"

namespace moist
{
    inline geo::Box2d create_interface_cell_bbox2d(const geo::index_t& c, const geo::Mesh& mesh, const geo::Attribute<bool>& v_interface)
    {
        geo::Box2d aabb;
        aabb.xy_min[0] = std::numeric_limits<double>::max();
        aabb.xy_min[1] = std::numeric_limits<double>::max();
        aabb.xy_max[0] = std::numeric_limits<double>::lowest();
        aabb.xy_max[1] = std::numeric_limits<double>::lowest();

        for (const geo::index_t v : moist::geometry::cell_vertices(c, mesh))
        {
            if (!v_interface[v])
            {
                continue;
            }

            // Again this only works with z-growing...
            const double* p = mesh.vertices.point_ptr(v);
            aabb.xy_min[0] = std::min(aabb.xy_min[0], p[0]);
            aabb.xy_min[1] = std::min(aabb.xy_min[1], p[1]);
            aabb.xy_max[0] = std::max(aabb.xy_max[0], p[0]);
            aabb.xy_max[1] = std::max(aabb.xy_max[1], p[1]);
        }

        return aabb;
    }

    inline geo::Box2d create_edge_box2d(const geo::vec2& p1, const geo::vec2& p2)
    {
        geo::Box2d box;

        box.xy_min[0] = std::min(p1.x, p2.x);
        box.xy_min[1] = std::min(p1.y, p2.y);

        box.xy_max[0] = std::max(p1.x, p2.x);
        box.xy_max[1] = std::max(p1.y, p2.y);

        return box;
    }

    inline geo::Box2d create_point_box2d(const geo::vec2& p)
    {
        geo::Box2d box;

        box.xy_min[0] = p.x;
        box.xy_min[1] = p.y;

        box.xy_max[0] = p.x;
        box.xy_max[1] = p.y;

        return box;
    }

    class LookupGrid
    {
    public:
        using MeshCells = std::unordered_set<geo::index_t>;
        using GridCell = std::pair<uint32_t, uint32_t>;

        // You can probably make this nice in methods but im no c++ wizard and have no time
        struct GridCellHash
        {
            std::size_t operator()(const GridCell& p) const
            {
                return std::hash<uint32_t>()(p.first) ^ (std::hash<uint32_t>()(p.second) << 1);
            }
        };

        LookupGrid(const geo::Mesh& mesh);

        void Initialize(const double resolution);
        void InsertCell(const geo::index_t c);
        std::vector<moist::LookupGrid::GridCell> GetCells(const geo::Box2d& aabb) const;

        // You can probably make this nice in methods but im no c++ wizard and have no time
        std::unordered_map<GridCell, MeshCells, GridCellHash> _grid;
    private:

        const geo::Mesh& _mesh;
        double _resolution;
        geo::vec2 _min_bounds;
        geo::vec2 _max_bounds;
        geo::vec2 _cell_size;
    };
}

#endif
