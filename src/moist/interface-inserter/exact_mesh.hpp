#ifndef MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP
#define MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

#include <geogram/mesh/mesh.h>
#include "moist/core/defines.hpp"

namespace moist
{
    constexpr std::size_t NO_VERTEX = std::size_t(-1U);

    typedef CGAL::Exact_predicates_exact_constructions_kernel ExactKernel;
    typedef ExactKernel::Segment_3 Segment;
    typedef ExactKernel::Line_3 Line;

    class ExactMesh
    {
    public:
        struct ExactPoint
        {
            moist::ExactKernel::Point_3 _p;
            geo::index_t _v; // Only non-interface vertices must have this set, to retranslate the new geometry back into the main mesh later... All others must be geo::NO_VERTEX
            bool _deleted;
            ExactPoint(const geo::vec3 p) : _p(ExactKernel::Point_3(p.x, p.y, p.z)), _v(geo::NO_VERTEX), _deleted(false) {}
            ExactPoint(const geo::vec3 p, geo::index_t v) : _p(ExactKernel::Point_3(p.x, p.y, p.z)), _v(v), _deleted(false) {}

            const double x() const { return CGAL::to_double(_p.x()); }
            const double y() const { return CGAL::to_double(_p.y()); }
            const double z() const { return CGAL::to_double(_p.z()); }
        };

        struct ExactCell
        {
            std::array<std::size_t, 4> _points;
            geo::index_t _c; // maybe useful for metrics
            bool _deleted;
            ExactCell(const std::size_t v0, const std::size_t v1, const std::size_t v2, const std::size_t v3) : _points({v0, v1, v2, v3}), _c(geo::NO_CELL), _deleted(false) {};
        };

        ExactMesh() = default;
        ~ExactMesh() = default;

        std::size_t Add(const ExactPoint& p);
        std::size_t Add(const ExactCell& cell);

        void DeletePoint(const std::size_t v);
        void DeleteCell(const std::size_t c);

        const ExactPoint& Point(const std::size_t& index) const;
        const ExactCell& Cell(const std::size_t& index) const;

        const std::vector<ExactPoint>& Points() const { return _points; }
        const std::vector<ExactCell>& Cells() const { return _cells; }

    #ifndef NDEBUG
        void DebugMesh(const std::filesystem::path& file);
    #endif

    private:
        std::vector<ExactPoint> _points;
        std::vector<ExactCell> _cells;
        // std::unordered_map<std::size_t, ExactPoint> _points;
        // std::unordered_map<std::size_t, ExactCell> _cells;

        std::size_t _pid = 0;
        std::size_t _tid = 0;
    };
}

#endif // MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP
