#ifndef MOIST_INTERRACE_INSERTER_EXACT_TYPES_HPP_
#define MOIST_INTERRACE_INSERTER_EXACT_TYPES_HPP_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

#include <geogram/mesh/mesh.h>
#include "moist/core/defines.hpp"

namespace moist::exact
{
    using index_t = std::size_t;
    constexpr std::size_t NO_VERTEX = std::size_t(-1U);

    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::Segment_3 Segment;
    typedef Kernel::Line_3 Line;

    enum class VertexCorrespondence
    {
        A, B, AB, UNINIT
    };

    struct Point
    {
        moist::exact::Kernel::Point_3 _p;
        GEO::index_t _v; // Only non-interface vertices must have this set, to retranslate the new geometry back into the main mesh later... All others must be geo::NO_VERTEX
        bool _deleted;
        bool _fixed;
        bool _interface;
        moist::exact::VertexCorrespondence _correspondence;

        std::size_t _other;
        Point(const GEO::vec3 p) : _p(moist::exact::Kernel::Point_3(p.x, p.y, p.z)), _v(GEO::NO_VERTEX), _deleted(false), _fixed(false), _other(moist::exact::NO_VERTEX) {}
        Point(const double x, const double y, const double z) : _p(moist::exact::Kernel::Point_3(x, y, z)), _v(GEO::NO_VERTEX), _deleted(false), _fixed(false), _other(moist::exact::NO_VERTEX) {}
        Point(const GEO::vec3 p, GEO::index_t v) : _p(moist::exact::Kernel::Point_3(p.x, p.y, p.z)), _v(v), _deleted(false), _fixed(false), _other(moist::exact::NO_VERTEX) {}
        Point(const moist::exact::Kernel::Point_3& p) : _p(p), _v(GEO::NO_VERTEX), _deleted(false), _fixed(false), _other(moist::exact::NO_VERTEX) {}
        Point(const moist::exact::Kernel::FT x, const moist::exact::Kernel::FT y, const moist::exact::Kernel::FT z, const moist::exact::VertexCorrespondence correspondence) :
                            _p(moist::exact::Kernel::Point_3(x, y, z)),
                            _v(GEO::NO_VERTEX), _deleted(false),
                            _fixed(false),
                            _other(moist::exact::NO_VERTEX),
                            _correspondence(correspondence) {}
        // Point(const moist::exact::Point& p) : _p(p._p), _interface(p._interface) {};
        Point(const GEO::vec3& p, const bool interface) :
            _p(moist::exact::Kernel::Point_3(p.x, p.y, p.z)),
            _v(GEO::NO_VERTEX),
            _deleted(false),
            _fixed(false),
            _other(moist::exact::NO_VERTEX),
            _correspondence(moist::exact::VertexCorrespondence::A),
            _interface(interface) {};

        bool operator==(const moist::exact::Point& other) const { return _p == other._p; }

        const double x() const { return CGAL::to_double(_p.x()); }
        const double y() const { return CGAL::to_double(_p.y()); }
        const double z() const { return CGAL::to_double(_p.z()); }

        void set(const moist::exact::Point p)
        {
            _p = moist::exact::Kernel::Point_3(p._p);
            _v = p._v;
            _deleted = p._deleted;
            _fixed = p._fixed;
        }
    };

    struct Point2
    {
        moist::exact::Kernel::Point_2 _p;
        moist::exact::VertexCorrespondence _correspondence;

        Point2(const moist::exact::Kernel::Point_2& p, const moist::exact::VertexCorrespondence correspondence) : _p(moist::exact::Kernel::Point_2(p.x(), p.y())), _correspondence(correspondence) {};
        Point2(const moist::exact::Kernel::Point_3& p, const moist::exact::VertexCorrespondence correspondence) : _p(moist::exact::Kernel::Point_2(p.x(), p.y())), _correspondence(correspondence) {};

        bool operator==(const moist::exact::Point2& other) const { return _p == other._p; }

        const double x() const { return CGAL::to_double(_p.x()); }
        const double y() const { return CGAL::to_double(_p.y()); }
    };

    class IndexedPoint
    {
        friend class moist::exact::Point;
    public:
        IndexedPoint(moist::exact::Point& point, const moist::exact::index_t index) : _p(point._p), _index(index) {}
        IndexedPoint(moist::exact::Kernel::Point_3/*&*/ point, const moist::exact::index_t index) : _p(point), _index(index) {}

        const moist::exact::index_t Index() const { return _index; }
        const double x() const { return CGAL::to_double(_p.x()); }
        const double y() const { return CGAL::to_double(_p.y()); }
        const double z() const { return CGAL::to_double(_p.z()); }

        bool operator==(const moist::exact::IndexedPoint& other) const { return _p == other._p; }
        bool operator==(const moist::exact::Point& other) const { return _p == other._p; }
    private:
        moist::exact::Kernel::Point_3 _p;
        moist::exact::index_t _index;
    };

    struct Triangle
    {
        moist::exact::Kernel::Triangle_3 _t;
        std::array<moist::exact::Point, 3> _points;
        Triangle(const moist::exact::Point& p0, const moist::exact::Point& p1, const moist::exact::Point& p2) : _points({ p0, p1, p2 })
        {
            _t = {p0._p, p1._p, p2._p};
        }
    };

    struct Facet
    {
        std::array<std::size_t, 3> _points;
        Facet(const std::size_t v0, const std::size_t v1, const std::size_t v2) : _points({v0, v1, v2}) {}

        std::size_t& operator[](std::size_t index) { moist::assert::in_range(index, 0, 3); return _points.at(index); }
        const std::size_t& operator[](std::size_t index) const { moist::assert::in_range(index, 0, 3); return _points.at(index); }
    };

    enum class CellType
    {
        I, II, III
    };

    struct Cell
    {
        std::array<std::size_t, 4> _points;
        bool _deleted;
        CellType _type;

        Cell(const std::size_t v0, const std::size_t v1, const std::size_t v2, const std::size_t v3, const moist::exact::CellType type) :
            _points({v0, v1, v2, v3}),
            _type(type),
            _deleted(false) {}

        std::size_t& operator[](std::size_t index) { moist::assert::in_range(index, 0, 4); return _points.at(index); }
        const std::size_t& operator[](std::size_t index) const { moist::assert::in_range(index, 0, 4); return _points.at(index); }
    };

    struct Edge
    {
        std::array<std::size_t, 2> _points;
        bool _deleted;
        Edge(const std::size_t v0, const std::size_t v1) : _points({v0, v1}), _deleted(false) {};

        bool operator==(const moist::exact::Edge& other) const { return (_points[0] == other._points[0] && _points[1] == other._points[1]) || (_points[0] ==  other._points[1] && _points[1] ==  other._points[0]); }
        std::size_t& operator[](std::size_t index) { moist::assert::in_range(index, 0, 2); return _points.at(index); }
        const std::size_t& operator[](std::size_t index) const { moist::assert::in_range(index, 0, 2); return _points.at(index); }
    };

    struct EdgePoints
    {
        moist::exact::Point p0;
        moist::exact::Point p1;

        const bool IsInterface() const { return p0._v == GEO::NO_VERTEX && p1._v == GEO::NO_VERTEX; }
    };

}


namespace std
{
    template<>
    struct hash<moist::exact::Edge>
    {
        std::size_t operator()(const moist::exact::Edge& e) const
        {
            std::size_t a = std::min(e._points[0], e._points[1]);
            std::size_t b = std::max(e._points[0], e._points[1]);
            return std::hash<std::size_t>{}(a) ^ (std::hash<std::size_t>{}(b) << 1);
        }
    };
}

#endif // MOIST_INTERRACE_INSERTER_EXACT_TYPES_HPP_
