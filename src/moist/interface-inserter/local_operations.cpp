#include "local_operations.hpp"

#include "moist/core/geometry.inl"
#include "moist/core/attributes.inl"
#include "moist/core/utils.hpp"

#include "new_predicates.inl"
#include "geometry_exact.inl"

#ifndef NDEBUG
#define VECE(a, b, c, d) ((a) == (b) || (a) == (c) || (a) == (d) || \
                          (b) == (c) || (b) == (d) || \
                          (c) == (d))
#endif // NDEBUG

static constexpr std::array<std::pair<std::size_t, std::size_t>, 6> TET_EDGE_DESCRIPTOR =
{{
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
}};

std::vector<moist::CrossedEdgeExact> moist::operation::exact::FindIntersectedEdges(const moist::exact::EdgePoints& edge, const std::size_t& c, const moist::ExactMesh& mesh)
{
    const auto& cell = mesh.Cell(c);
    auto edges = std::vector<moist::CrossedEdgeExact>();

    for (const auto& [i, j] : TET_EDGE_DESCRIPTOR)
    {
        const auto& v0 = cell._points[i];
        const auto& v1 = cell._points[j];
        const auto& p0 = mesh.Point(v0);
        const auto& p1 = mesh.Point(v1);

        if (p0._v != moist::exact::NO_VERTEX || p1._v != moist::exact::NO_VERTEX)
        {
            continue;
        }

        const auto o_intersection = moist::geometry::exact::intersection(edge, {p0, p1});
        if (!o_intersection.has_value())
        {
            continue;
        }

        bool cont = false;
        for (const std::size_t v : cell._points)
        {
            if (moist::geometry::exact::points_are_equal(mesh.Point(v)._p, o_intersection.value()._p))
            {
                cont = true;
                break;
            }
        }

        if (cont)
        {
            continue;
        }

        if (std::find_if(edges.begin(), edges.end(), [&](const moist::CrossedEdgeExact& edge) { return edge.p == o_intersection.value(); }) == edges.end())
        {
            edges.push_back(moist::CrossedEdgeExact { v0, v1, o_intersection.value(), moist::exact::NO_VERTEX });
        }
    }

    return edges;
}

void moist::operation::exact::InsertVertexOnCellBoundaryFacet(const std::size_t &c, const std::size_t &v, moist::ExactMesh &mesh)
{
    std::size_t v_opposite;
    for (const std::size_t& cv : mesh.Cell(c)._points)
    {
        const auto& cp = mesh.Point(cv);
        if (cp._v != geo::NO_VERTEX)
        {
            v_opposite = cv;
            break;
        }
    }

    const auto [v0, v1, v2] = moist::geometry::exact::other(c, v_opposite, mesh);
    mesh.Add(moist::exact::Cell(v_opposite, v0, v1, v));
    mesh.Add(moist::exact::Cell(v_opposite, v1, v2, v));
    mesh.Add(moist::exact::Cell(v_opposite, v2, v0, v));
    mesh.DeleteCell(c);
}

bool moist::operation::exact::InsertVertexOnCellBoundaryEdgeOld(const std::size_t &c, const std::size_t &v, moist::ExactMesh &mesh, const double eps)
{
    const auto& point = mesh.Point(v);
    const auto& cell = mesh.Cell(c);
    bool ret = true;

    for (const auto& [i, j] : TET_EDGE_DESCRIPTOR)
    {
        const std::size_t v0 = cell._points[i];
        const std::size_t v1 = cell._points[j];
        const auto p0 = mesh.Point(v0);
        const auto p1 = mesh.Point(v1);
        if (p0._v != moist::exact::NO_VERTEX || p1._v == moist::exact::NO_VERTEX)
        {
            continue;
        }

        bool collinear = (eps == 0.0) ? CGAL::collinear(p0._p, p1._p, point._p) : moist::new_predicates::approx_collinear(p0._p, p1._p, point._p, eps);
        if (collinear)
        {
            //if (CGAL::collinear_are_ordered_along_line(p0._p, point._p, p1._p))
            //{
                std::array<int, 2> v_other;
                int index = 0;
                for (int k = 0; k < 4; k++)
                {
                    if (k != i && k != j)
                    {
                        v_other[index++] = k;
                    }
                }

                const std::size_t v2 = cell._points[v_other[0]];
                const std::size_t v3 = cell._points[v_other[1]];

                mesh.Add(moist::exact::Cell(v, v0, v2, v3));
                mesh.Add(moist::exact::Cell(v, v1, v2, v3));

                if (mesh.Point(v0).z() == 10.0 && mesh.Point(v2).z() == 10.0 && mesh.Point(v3).z() == 10.0)
                {
                    OOC_DEBUG("wat??");
                }
                if (mesh.Point(v1).z() == 10.0 && mesh.Point(v2).z() == 10.0 && mesh.Point(v3).z() == 10.0)
                {
                    OOC_DEBUG("wat??");
                }

            #ifndef NDEBUG
                double vol0 = moist::new_predicates::approx_volume(mesh.Point(v), mesh.Point(v0), mesh.Point(v2), mesh.Point(v3));
                double vol1 = moist::new_predicates::approx_volume(mesh.Point(v), mesh.Point(v1), mesh.Point(v2), mesh.Point(v3));
                if (vol0 == 0.0)
                {
                    double check = moist::new_predicates::approx_volume(mesh.Point(v), mesh.Point(v0), mesh.Point(v2), mesh.Point(v3));
                    OOC_DEBUG("created zero   tet!");
                    OOC_DEBUG("point v: " << mesh.Point(v)._p);
                    OOC_DEBUG("point v0: " << mesh.Point(v0)._p);
                    OOC_DEBUG("point v1: " << mesh.Point(v1)._p);
                    OOC_DEBUG("point v2: " << mesh.Point(v2)._p);
                    OOC_DEBUG("point v3: " << mesh.Point(v3)._p);
                    ret = false;
                }

                if (vol1 == 0.0)
                {
                    double check = moist::new_predicates::approx_volume(mesh.Point(v), mesh.Point(v1), mesh.Point(v2), mesh.Point(v3));
                    OOC_DEBUG("created zero   tet!");
                    OOC_DEBUG("point v: " << mesh.Point(v)._p);
                    OOC_DEBUG("point v0: " << mesh.Point(v0)._p);
                    OOC_DEBUG("point v1: " << mesh.Point(v1)._p);
                    OOC_DEBUG("point v2: " << mesh.Point(v2)._p);
                    OOC_DEBUG("point v3: " << mesh.Point(v3)._p);

                    ret = false;
                }
            #endif // NDEBUG
                mesh.DeleteCell(c);
                break;
            //}
        }
    }

    return ret;
}

static inline std::array<std::size_t, 2> get_edge_vertices(const moist::predicates::PointInTet& location)
{
    switch (location)
    {
        case moist::predicates::PointInTet::EDGE01: return {0, 1};
        case moist::predicates::PointInTet::EDGE02: return {0, 2};
        case moist::predicates::PointInTet::EDGE03: return {0, 3};
        case moist::predicates::PointInTet::EDGE12: return {1, 2};
        case moist::predicates::PointInTet::EDGE13: return {1, 3};
        case moist::predicates::PointInTet::EDGE23: return {2, 3};
        default:
            throw std::invalid_argument("Invalid or non-edge enum value");
    }
}

void moist::operation::exact::InsertVertexOnCellBoundaryEdge(const std::size_t & c, const std::size_t & v, const moist::predicates::PointInTet& location, moist::ExactMesh & mesh)
{
    const auto& point = mesh.Point(v);
    const auto& cell = mesh.Cell(c);

    const auto v_edge = get_edge_vertices(location);
    std::array<std::size_t, 2> v_other;
    int index = 0;
    for (int k = 0; k < 4; k++)
    {
        if (k != v_edge[0] && k != v_edge[1])
        {
            v_other[index++] = k;
        }
    }

    const std::size_t v0 = cell._points[v_edge[0]];
    const std::size_t v1 = cell._points[v_edge[1]];
    const std::size_t v2 = cell._points[v_other[0]];
    const std::size_t v3 = cell._points[v_other[1]];

    const auto z0 = mesh.Point(v0).z();
    const auto z1 = mesh.Point(v1).z();

    if (mesh.Point(v0)._v != moist::exact::NO_VERTEX || mesh.Point(v1)._v != moist::exact::NO_VERTEX)
    {
        return;
    }

    mesh.Add(moist::exact::Cell(v, v0, v2, v3));
    mesh.Add(moist::exact::Cell(v, v1, v2, v3));
    mesh.DeleteCell(c);
}

void moist::operation::exact::SplitEdge1_2(const std::size_t &c, const moist::CrossedEdgeExact &edge, moist::ExactMesh &mesh, std::vector<moist::exact::Cell>& created_cells)
{
    std::vector<std::size_t> others; // should prolly be static but how does this impact possible parallelization
    for (const auto& v : mesh.Cell(c)._points)
    {
        if (v != edge.v0 && v != edge.v1)
        {
            others.push_back(v);
        }
    }

    created_cells.push_back(moist::exact::Cell(edge.vp, edge.v0, others[0], others[1]));
    created_cells.push_back(moist::exact::Cell(edge.vp, edge.v1, others[0], others[1]));
    mesh.DeleteCell(c);
}

void moist::operation::exact::SplitEdge1_3(const std::size_t &c, const moist::CrossedEdgeExact &edge0, const moist::CrossedEdgeExact &edge1, moist::ExactMesh &mesh, std::vector<moist::exact::Cell>& created_cells)
{
    std::vector<std::size_t> others; // should prolly be static but how does this impact possible parallelization
    for (const auto& v : mesh.Cell(c)._points)
    {
        if (v != edge0.v0 && v != edge0.v1 && v != edge1.v0 && v != edge1.v1)
        {
            others.push_back(v);
        }
    }

#ifndef NDEBUG
    if (others.size() != 1)
    {
        const auto& pp = mesh.Cell(c)._points;
        OOC_ERROR("invalid 1->3 split -- found " << others.size() << " uninvolved vertices");
    }
#endif // NDEBUG

    const auto v_shared = (edge0.v0 == edge1.v0 || edge0.v0 == edge1.v1) ? edge0.v0 : edge0.v1;
    const auto v_other_e1 = (edge1.v0 == v_shared) ? edge1.v1 : edge1.v0;
    const auto v_other_e0 = (edge0.v0 == v_shared) ? edge0.v1 : edge0.v0;
    const auto v_other = others.at(0);

    created_cells.push_back(moist::exact::Cell {v_other, edge0.vp, edge1.vp, v_shared});
    created_cells.push_back(moist::exact::Cell {v_other, edge0.vp, edge1.vp, v_other_e0});
    created_cells.push_back(moist::exact::Cell {v_other, edge1.vp, v_other_e0, v_other_e1});
    mesh.DeleteCell(c);
}
