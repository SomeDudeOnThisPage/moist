#ifndef MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP
#define MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

#include <tetgen.h>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

#include "exact_types.hpp"
#include "lookup_grid.hpp"

namespace moist
{
    class ExactMesh
    {
    public:
        ExactMesh();
        ~ExactMesh() = default;

        void ResetGrid(const double resolution);
        void ResetMesh();

        std::size_t Add(const moist::exact::Point p);
        std::size_t Add(const moist::exact::Cell cell, const bool initialization = false);
        std::size_t Add(const moist::exact::Edge e);
        std::size_t Add(const moist::exact::Facet f);

        void DeletePoint(const std::size_t v);
        void DeleteCell(const std::size_t c);
        void FlushDeletedElements();

        moist::exact::Cell& Cell(const std::size_t& index);
        moist::exact::Point& Point(const std::size_t& index);
        const moist::exact::Point& Point(const std::size_t& index) const;
        const moist::exact::Cell& Cell(const std::size_t& index) const;

        const std::vector<moist::exact::Point>& Points() const { return _points; }
        const std::vector<moist::exact::Cell>& Cells() const { return _cells; }
        const std::vector<moist::exact::Facet>& Facets() const { return _facets; }
        const std::unordered_set<moist::exact::Edge>& Edges() const { return _edges; }

        const std::size_t NbPoints() const { return _points.size(); }
        const std::size_t NbCells() const { return _cells.size(); }
        const std::size_t NbFacets() const { return _facets.size(); }
        const std::size_t NbEdges() const { return _edges.size(); }

        std::shared_ptr<moist::LookupGridExact> Grid() { return _grid; }

    #ifndef NDEBUG
        void DebugMesh(const std::filesystem::path& file);
    #endif

    private:
        std::vector<moist::exact::Point> _points;
        std::vector<moist::exact::Cell> _cells;
        std::vector<moist::exact::Facet> _facets;
        std::unordered_set<moist::exact::Edge> _edges;

        std::size_t _pid = 0;
        std::size_t _tid = 0;
        std::size_t _eid = 0;

        std::shared_ptr<moist::LookupGridExact> _grid;
        moist::LookupPointGrid _point_grid;
    };

    /**
     * @brief Creates a two dimensional bounding box of a cell, projected onto the xy-plane.
     *
     * @param c
     * @param mesh
     * @return geo::Box2d
     */
    inline geo::Box2d create_cell_bbox2d_exact(const std::size_t& c, const moist::ExactMesh& mesh)
    {
        geo::Box2d aabb;
        aabb.xy_min[0] = std::numeric_limits<double>::max();
        aabb.xy_min[1] = std::numeric_limits<double>::max();
        aabb.xy_max[0] = std::numeric_limits<double>::lowest();
        aabb.xy_max[1] = std::numeric_limits<double>::lowest();

        const auto cell = mesh.Cell(c);
        for (const std::size_t v : cell._points)
        {
            const auto point = mesh.Point(v);
            aabb.xy_min[0] = std::min(aabb.xy_min[0], point.x());
            aabb.xy_min[1] = std::min(aabb.xy_min[1], point.y());
            aabb.xy_max[0] = std::max(aabb.xy_max[0], point.x());
            aabb.xy_max[1] = std::max(aabb.xy_max[1], point.y());
        }

        return aabb;
    }

    inline geo::Box2d create_triangle_bbox2d_exact(const moist::exact::Triangle& triangle)
    {
        geo::Box2d aabb;
        aabb.xy_min[0] = std::numeric_limits<double>::max();
        aabb.xy_min[1] = std::numeric_limits<double>::max();
        aabb.xy_max[0] = std::numeric_limits<double>::lowest();
        aabb.xy_max[1] = std::numeric_limits<double>::lowest();

        for (const auto point : triangle._points)
        {
            aabb.xy_min[0] = std::min(aabb.xy_min[0], point.x());
            aabb.xy_min[1] = std::min(aabb.xy_min[1], point.y());
            aabb.xy_max[0] = std::max(aabb.xy_max[0], point.x());
            aabb.xy_max[1] = std::max(aabb.xy_max[1], point.y());
        }

        return aabb;
    }

    inline geo::Box2d create_interface_cell_bbox2d_exact(const std::size_t& c, const moist::ExactMesh& mesh)
    {
        geo::Box2d aabb;
        aabb.xy_min[0] = std::numeric_limits<double>::max();
        aabb.xy_min[1] = std::numeric_limits<double>::max();
        aabb.xy_max[0] = std::numeric_limits<double>::lowest();
        aabb.xy_max[1] = std::numeric_limits<double>::lowest();

        const auto cell = mesh.Cell(c);
        for (const std::size_t v : cell._points)
        {
            const auto point = mesh.Point(v);
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

    inline geo::Box2d create_edge_box2d_exact(const moist::exact::EdgePoints& edge)
    {
        geo::Box2d box;

        box.xy_min[0] = std::min(edge.p0.x(), edge.p1.x());
        box.xy_min[1] = std::min(edge.p0.y(), edge.p1.y());

        box.xy_max[0] = std::max(edge.p0.x(), edge.p1.x());
        box.xy_max[1] = std::max(edge.p0.y(), edge.p1.y());

        return box;
    }

    inline geo::Box2d create_point_box2d_exact(const moist::exact::Point& p)
    {
        geo::Box2d box;

        box.xy_min[0] = p.x() + 1e-12;
        box.xy_min[1] = p.y() + 1e-12;

        box.xy_max[0] = p.x() + 1e-12;
        box.xy_max[1] = p.y() + 1e-12;

        return box;
    }
}

#endif // MOIST_INTERFACE_INSERTER_EXACT_MESH_HPP
