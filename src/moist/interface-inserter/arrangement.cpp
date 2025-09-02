#include "arrangement.hpp"

#include <vector>

#include <CGAL/Arr_dcel_base.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>
#include <CGAL/intersections.h>

using Traits_2 = CGAL::Arr_segment_traits_2<moist::exact::Kernel>;
using Dcel = CGAL::Arr_extended_dcel<Traits_2, moist::exact::VertexCorrespondence, std::size_t, std::size_t>;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2, Dcel>;
typedef CGAL::Polygon_2<moist::exact::Kernel> Polygon_2;
typedef moist::exact::Kernel::Point_2 Point_2;

//namespace
//{
static bool triangles_share_vertex(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1)
{
    for (int i = 0; i < 3; i++)
    {
        const moist::exact::Kernel::Point_3& p1 = t0._t.vertex(i);
        for (int j = 0; j < 3; j++)
        {
            const moist::exact::Kernel::Point_3& p2 = t1._t.vertex(j);
            if (p1 == p2)
            {
                return true;
            }
        }
    }
    return false;
}

bool moist::exact::arrangeable(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1)
{
    auto result = CGAL::intersection(t0._t, t1._t);

    if (!result) return false; // no intersection at all

    if (const moist::exact::Kernel::Point_3* p = std::get_if<moist::exact::Kernel::Point_3>(&*result))
    {
        return false;
    }
    else if (const moist::exact::Kernel::Segment_3* s = std::get_if<moist::exact::Kernel::Segment_3>(&*result))
    {
        return false;
    }
    else if (const moist::exact::Kernel::Triangle_3* tri = std::get_if<moist::exact::Kernel::Triangle_3>(&*result))
    {
        return true;
    }
    else
    {
        return true; // polyhon
    }

    // if the triangles only intersect on their boundary or one vertex, we don't care.
    // return CGAL::do_intersect(t0._t, t1._t);
}

class VertexObserver : public CGAL::Arr_observer<Arrangement_2> {
    public:
        VertexObserver(Arrangement_2& arrangement, const std::array<moist::exact::Point2, 3> points_a, const std::array<moist::exact::Point2, 3> points_b)
            : CGAL::Arr_observer<Arrangement_2>(arrangement), _correspondence(moist::exact::VertexCorrespondence::A), _points_a(points_a), _points_b(points_b)
        {
            CGAL_precondition(arrangement.is_empty());
        }

        virtual void after_create_vertex(Arrangement_2::Vertex_handle v)
        {
            const auto& point = v->point();
            moist::exact::VertexCorrespondence correspondence;

            bool is_a = std::find_if(_points_a.begin(), _points_a.end(), [&](const auto p) -> bool { return p._p == v->point(); }) != _points_a.end();
            bool is_b = std::find_if(_points_b.begin(), _points_b.end(), [&](const auto p) -> bool { return p._p == v->point(); }) != _points_b.end();

            if (is_a && is_b || !(is_a || is_b)) // if the point is shared or new (edge split point)
            {
                correspondence = moist::exact::VertexCorrespondence::AB;
            }
            else
            {
                correspondence = is_a ? moist::exact::VertexCorrespondence::A : moist::exact::VertexCorrespondence::B;
            }

            // if the vertex lies inside the triangle defined by the other set of points, it must become AB
            if (correspondence == moist::exact::VertexCorrespondence::A && this->InTriangle(point, _points_b))
            {
                correspondence = moist::exact::VertexCorrespondence::AB;
            }
            else if (correspondence == moist::exact::VertexCorrespondence::B && this->InTriangle(point, _points_a))
            {
                correspondence = moist::exact::VertexCorrespondence::AB;
            }

            v->set_data(correspondence);
        }
    private:
        const std::array<moist::exact::Point2, 3> _points_a;
        const std::array<moist::exact::Point2, 3> _points_b;
        moist::exact::VertexCorrespondence _correspondence;

        bool InTriangle(const moist::exact::Kernel::Point_2& point, const std::array<moist::exact::Point2, 3>& triangle) const
        {
            const std::array<moist::exact::Kernel::Point_2, 3> it_polygon = {
                triangle[0]._p,
                triangle[1]._p,
                triangle[2]._p
            };
            const auto side = CGAL::bounded_side_2(it_polygon.begin(), it_polygon.end(), point);
            return side == CGAL::ON_BOUNDARY || side == CGAL::ON_BOUNDED_SIDE;
        }
};

static moist::exact::VertexCorrespondence get_correspondence(const moist::exact::Kernel::Point_2& p, const std::vector<moist::exact::Point2>& points)
{
    for (const auto point : points)
    {
        if (point._p == p)
        {
            return point._correspondence;
        }
    }
    OOC_ERROR("invalid point without correspondence in triangulation");
}

#include "moist/core/utils.hpp"

static geo::Mesh dbg_arrangement(3);
static geo::Mesh dbg_arrangement_triangulation(3);

moist::exact::Triangulation moist::exact::arrange(const moist::exact::Triangle& t0, const moist::exact::Triangle& t1)
{
    //dbg_arrangement.clear();
    //dbg_arrangement_triangulation.clear();
    moist::exact::Triangulation triangulation;

    Arrangement_2 arrangement;
    const std::array<moist::exact::Point2, 3> points_a =
    {{
        moist::exact::Point2(t0._t.vertex(0), moist::exact::VertexCorrespondence::A),
        moist::exact::Point2(t0._t.vertex(1), moist::exact::VertexCorrespondence::A),
        moist::exact::Point2(t0._t.vertex(2), moist::exact::VertexCorrespondence::A)
    }};

    const std::array<moist::exact::Point2, 3> points_b =
    {{
        moist::exact::Point2(t1._t.vertex(0), moist::exact::VertexCorrespondence::B),
        moist::exact::Point2(t1._t.vertex(1), moist::exact::VertexCorrespondence::B),
        moist::exact::Point2(t1._t.vertex(2), moist::exact::VertexCorrespondence::B)
    }};

    /*dbg_arrangement.facets.create_triangle(
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_a[0]._p.x()), CGAL::to_double(points_a[0]._p.y()), 0)),
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_a[1]._p.x()), CGAL::to_double(points_a[1]._p.y()), 0)),
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_a[2]._p.x()), CGAL::to_double(points_a[2]._p.y()), 0))
    );

    dbg_arrangement.facets.create_triangle(
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_b[0]._p.x()), CGAL::to_double(points_b[0]._p.y()), 0)),
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_b[1]._p.x()), CGAL::to_double(points_b[1]._p.y()), 0)),
        dbg_arrangement.vertices.create_vertex(geo::vec3(CGAL::to_double(points_b[2]._p.x()), CGAL::to_double(points_b[2]._p.y()), 0))
    );*/
    // moist::utils::geogram::save("dbg_arrangement.mesh", dbg_arrangement);

    VertexObserver observer(arrangement, points_a, points_b);

    moist::exact::Kernel::Segment_2 segments[] = {
        moist::exact::Kernel::Segment_2(points_a.at(0)._p, points_a.at(1)._p),
        moist::exact::Kernel::Segment_2(points_a.at(1)._p, points_a.at(2)._p),
        moist::exact::Kernel::Segment_2(points_a.at(2)._p, points_a.at(0)._p),
        moist::exact::Kernel::Segment_2(points_b.at(0)._p, points_b.at(1)._p),
        moist::exact::Kernel::Segment_2(points_b.at(1)._p, points_b.at(2)._p),
        moist::exact::Kernel::Segment_2(points_b.at(2)._p, points_b.at(0)._p)
    };

    CGAL::insert(arrangement, segments, segments + 6);

    // The triangulation we use does not support data on points like the arrangements... so keep a copy of all points with their corresponding mesh(es) here to reindex later...
    std::vector<moist::exact::Point2> points;

    for (auto vit = arrangement.vertices_begin(); vit != arrangement.vertices_end(); vit++)
    {
        points.push_back(moist::exact::Point2(vit->point(), vit->data()));
    }

    for (auto fit = arrangement.faces_begin(); fit != arrangement.faces_end(); fit++)
    {
        if (fit->is_unbounded())
        {
            continue;
        }

        // Disjoint
        if (fit->number_of_outer_ccbs() != 1)
        {
            continue;
        }

        Polygon_2 polygon;
        auto circ = fit->outer_ccb();
        auto start = circ;

        do
        {
            polygon.push_back(circ->source()->point());
            ++circ;
        }
        while (circ != start);

        if (!polygon.is_simple())
        {
            continue;
        }

        if (polygon.orientation() == CGAL::CLOCKWISE)
        {
            polygon.reverse_orientation();
        }

        // TODO: we don't always need to do a full triangulation, since most overlaps are by their nature already triangulated.
        //       some notable exception is when two overlapping triangles share zero
        std::vector<Polygon_2> triangles;
        CGAL::Polygon_triangulation_decomposition_2<moist::exact::Kernel> decomposition; // why is every single CGAL API structured COMPLETELY differently from all others??
        decomposition(polygon, std::back_inserter(triangles));

        const auto z = t0._t.vertex(0).z();
        for (const auto& triangle : triangles)
        {
            const auto c0 = get_correspondence(triangle.vertex(0), points);
            const auto c1 = get_correspondence(triangle.vertex(1), points);
            const auto c2 = get_correspondence(triangle.vertex(2), points);

            triangulation.push_back(moist::exact::Triangle(
                moist::exact::Point(triangle.vertex(0).x(), triangle.vertex(0).y(), z, c0),
                moist::exact::Point(triangle.vertex(1).x(), triangle.vertex(1).y(), z, c1),
                moist::exact::Point(triangle.vertex(2).x(), triangle.vertex(2).y(), z, c2)
            ));

            /*dbg_arrangement_triangulation.facets.create_triangle(
                dbg_arrangement_triangulation.vertices.create_vertex(geo::vec3(CGAL::to_double(triangle.vertex(0).x()), CGAL::to_double(triangle.vertex(0).y()), 0)),
                dbg_arrangement_triangulation.vertices.create_vertex(geo::vec3(CGAL::to_double(triangle.vertex(1).x()), CGAL::to_double(triangle.vertex(1).y()), 0)),
                dbg_arrangement_triangulation.vertices.create_vertex(geo::vec3(CGAL::to_double(triangle.vertex(2).x()), CGAL::to_double(triangle.vertex(2).y()), 0))
            );*/
        }
    }
    // moist::utils::geogram::save("dbg_arrangement_triangulation.mesh", dbg_arrangement_triangulation);
    return triangulation;
}
