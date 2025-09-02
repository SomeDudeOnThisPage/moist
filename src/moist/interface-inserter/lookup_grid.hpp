#ifndef MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP
#define MOIST_INTERFACE_INSERTER_LOOKUP_GRID_HPP

#include <ranges>
#include <unordered_set>
#include <vector>
#include <optional>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

#include "moist/core/defines.hpp"
#include "moist/core/core_interface.hpp"
#include "moist/core/geometry.inl"

#include "exact_types.hpp"

namespace moist
{
    class ExactMesh;
    class LookupGridExact
    {
    public:
        friend class ExactMesh;

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
        //void Initialize(const moist::ExactMesh& mesh, const double resolution);
        void InsertCell(const std::size_t c, const moist::ExactMesh& mesh);
        std::vector<moist::LookupGridExact::GridCell> GetCells(const geo::Box2d& aabb, const bool initialization = false) const;
        moist::LookupGridExact::MeshCells& GetMeshCells(const moist::LookupGridExact::GridCell& grid_cell);

        double Resolution() { return _resolution; }

        std::unordered_map<GridCell, MeshCells, GridCellHash> _grid;
    private:
        double _resolution;
        geo::vec2 _min_bounds;
        geo::vec2 _max_bounds;
        geo::vec2 _cell_size;
    };

    /**
     * @brief Lookup point binning for better performance to keep indexing during mesh creation
     * It stores 3D points but only on a 2d grid, since our meshes are very "flat" we don't really care.
     * Point equality is still handled in 3D though.
     */
    class LookupPointGrid
    {
    public:
        LookupPointGrid(const double resolution) : _resolution(resolution) {};
        std::optional<moist::exact::index_t> Get(moist::exact::Point point);
        void Add(moist::exact::Point point, const moist::exact::index_t index);
    private:
        struct Coordinates
        {
            long long gx, gy; // honestly our meshes won't be THAT large that we need exact arithmetic for resolving coordinates... same with LookupGridExact
            bool operator==(const Coordinates& other) const
            {
                return gx == other.gx && gy == other.gy;
            }
        };

        struct CoordinatesHash
        {
            std::size_t operator()(const Coordinates& c) const
            {
                std::size_t h1 = std::hash<long long>()(c.gx);
                std::size_t h2 = std::hash<long long>()(c.gy);
                return h1 ^ (h2 << 1);
            }
        };

        double _resolution;
        std::unordered_map<Coordinates, std::vector<moist::exact::IndexedPoint>, CoordinatesHash> _grid;
        Coordinates GetCoordinates(const moist::exact::Point& point);
    };
}

#endif
