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

#include "exact_mesh.hpp"

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

    inline geo::Box2d create_interface_cell_bbox2d_exact(const std::size_t& c, const moist::ExactMesh& mesh)
    {
        geo::Box2d aabb;
        aabb.xy_min[0] = std::numeric_limits<double>::max();
        aabb.xy_min[1] = std::numeric_limits<double>::max();
        aabb.xy_max[0] = std::numeric_limits<double>::lowest();
        aabb.xy_max[1] = std::numeric_limits<double>::lowest();

        const auto& cell = mesh.Cell(c);
        for (const std::size_t v : cell._points)
        {
            const auto& point = mesh.Point(v);
            if (point._v != geo::NO_VERTEX)
            {
                continue;
            }

            aabb.xy_min[0] = std::min(aabb.xy_min[0], point.x());
            aabb.xy_min[1] = std::min(aabb.xy_min[1], point.y());
            aabb.xy_max[0] = std::max(aabb.xy_max[0], point.x());
            aabb.xy_max[1] = std::max(aabb.xy_max[1], point.y());
        }

        return aabb;
    }

    inline geo::Box2d create_edge_box2d_exact(const std::size_t& v0, const std::size_t& v1, const moist::ExactMesh& mesh)
    {
        geo::Box2d box;

        const auto& p0 = mesh.Point(v0);
        const auto& p1 = mesh.Point(v1);

        box.xy_min[0] = std::min(p0.x(), p1.x());
        box.xy_min[1] = std::min(p0.y(), p1.y());

        box.xy_max[0] = std::max(p0.x(), p1.x());
        box.xy_max[1] = std::max(p0.y(), p1.y());

        return box;
    }

    inline geo::Box2d create_point_box2d_exact(const moist::ExactMesh::ExactPoint& p)
    {
        geo::Box2d box;

        box.xy_min[0] = p.x();
        box.xy_min[1] = p.y();

        box.xy_max[0] = p.x();
        box.xy_max[1] = p.y();

        return box;
    }

    class LookupGridExact
    {
    public:
        using MeshCells = std::unordered_set<std::size_t>;
        using GridCell = std::pair<std::size_t, std::size_t>;

        struct GridCellHash
        {
            std::size_t operator()(const GridCell& p) const
            {
                return std::hash<uint32_t>()(p.first) ^ (std::hash<uint32_t>()(p.second) << 1);
            }
        };

        LookupGridExact() = default;
        void Initialize(const moist::ExactMesh& mesh, const double resolution);
        void InsertCell(const std::size_t c);
        std::vector<moist::LookupGridExact::GridCell> GetCells(const geo::Box2d& aabb) const;

        std::unordered_map<GridCell, MeshCells, GridCellHash> _grid;
    private:
        const moist::ExactMesh* _mesh = nullptr; // ugly raw pointer never freed!!!
        double _resolution;
        geo::vec2 _min_bounds;
        geo::vec2 _max_bounds;
        geo::vec2 _cell_size;
    };
}

#endif
