#ifndef MOIST_INTERRACE_INSERTER_EXACT_TYPES_HPP_
#define MOIST_INTERRACE_INSERTER_EXACT_TYPES_HPP_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

namespace moist::exact
{
    constexpr std::size_t NO_VERTEX = std::size_t(-1U);

    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::Segment_3 Segment;
    typedef Kernel::Line_3 Line;

    struct Point
    {
        moist::exact::Kernel::Point_3 _p;
        geo::index_t _v; // Only non-interface vertices must have this set, to retranslate the new geometry back into the main mesh later... All others must be geo::NO_VERTEX
        bool _deleted;
        bool _fixed;
        Point(const geo::vec3 p) : _p(moist::exact::Kernel::Point_3(p.x, p.y, p.z)), _v(geo::NO_VERTEX), _deleted(false), _fixed(false) {}
        Point(const geo::vec3 p, geo::index_t v) : _p(moist::exact::Kernel::Point_3(p.x, p.y, p.z)), _v(v), _deleted(false), _fixed(false) {}
        Point(const moist::exact::Kernel::Point_3& p) : _p(p), _v(geo::NO_VERTEX), _deleted(false), _fixed(false) {}

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

    struct Cell
    {
        std::array<std::size_t, 4> _points;
        bool _deleted;
        Cell(const std::size_t v0, const std::size_t v1, const std::size_t v2, const std::size_t v3) : _points({v0, v1, v2, v3}), _deleted(false) {};

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

        const bool IsInterface() const { return p0._v == geo::NO_VERTEX && p1._v == geo::NO_VERTEX; }
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
