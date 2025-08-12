#include "slice.hpp"

#include <CGAL/squared_distance_3.h>

#include "moist/core/geometry.inl"

#include "exact_types.hpp"

#include "local_operations.hpp"
#include "new_predicates.inl"
#include "geometry_exact.inl"

#ifdef MOIST_OPTION_EXACT_PREDICATES

inline bool is_equal_test(double a, double b, double epsilon = 1e-4) {
    return std::fabs(a - b) <= epsilon;
}

void moist::MeshSlice::InsertTetrahedra(moist::MeshSlice& other, const moist::AxisAlignedPlane& plane)
{
    const auto start = other.ReorderInterfaceCells(plane);
    other.ConstructExactOverlayMesh(static_cast<std::size_t>(start));

    // begin inserting all other tets into this exactmesh.
    for (const auto cell : other._overlay.Cells())
    {
        if (cell._deleted)
        {
            continue;
        }

        // new idea... for each pair of intersecting tets, find the points they intersect in, and insert them in some order...
        // remember that order, and insert into both, but make a copy, dont use the entire mesh.
        // then, finally, delete the one existig tet, and insert all created tets into the mesh.
        // that way we keep perfect element locality... and... enable parallelization 👀?

        // in the first step, only handle tetrahedra which have a full facet on the interface
        // we need to find the corresponding points in this mesh, which are from the cell of the other mesh
        // this is in order to avoid problems with floating point inaccuracy -- after insertion, we only use these...
        const auto vertices = moist::geometry::exact::interface_vertices(cell, other._overlay);
        if (vertices.size() != 3)
        {
            continue;
        }



    }
}

/* private */ void moist::MeshSlice::ConstructExactOverlayMesh(const std::size_t start)
{
    auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
    std::unordered_map<geo::index_t, std::size_t> added_vertices;

    for (geo::index_t c = start; c < this->cells.nb(); c++)
    {
        size_t points[4];
        for (geo::index_t lv = 0; lv < 4; lv++)
        {
            const geo::index_t v = this->cells.vertex(c, lv);
            const auto po = this->vertices.point(v);

            if (added_vertices.contains(v))
            {
                points[lv] = added_vertices[v];
                continue;
            }
            else
            {
                points[lv] = _overlay.Add(moist::exact::Point(this->vertices.point(v), v_interface[v] ? geo::NO_VERTEX : v));
                added_vertices[v] = points[lv];
            }
        }

        this->_overlay.Add(moist::exact::Cell(points[0], points[1], points[2], points[3]), true);
    }

    this->_overlay.ResetGrid(15.0);
#ifndef NDEBUG
    this->_overlay.DebugMesh("exact_mesh.mesh");
#endif // NDEBUG
}

/* public */ void moist::MeshSlice::InsertEdges(const geo::Mesh& edge_mesh, const moist::AxisAlignedPlane& plane)
{
    auto v_interface = geo::Attribute<bool>(this->vertices.attributes(), "v_interface");
    for (const geo::index_t v : this->vertices)
    {
        if (!moist::predicates::point_on_plane(this->vertices.point(v), plane))
        {
            v_interface[v] = false;
            continue;
        }

        v_interface[v] = true;
        switch (plane.axis)
        {
            case moist::Axis::X:
                this->vertices.point(v).x = plane.extent;
                break;
            case moist::Axis::Y:
                this->vertices.point(v).y = plane.extent;
                break;
            case moist::Axis::Z:
                this->vertices.point(v).z = plane.extent;
                break;
            default:
                throw std::runtime_error("unknown axis definition");
        }
    }

    this->_start_interface_cell = this->ReorderInterfaceCells(plane);

    this->ConstructExactOverlayMesh(this->_start_interface_cell);

    long long timer_iv, timer_ie;
    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceVertices");
        this->InsertVertices(edge_mesh);
        timer_iv = timer.Elapsed();
    }

    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceEdges");
        this->InsertEdgesPrivate(edge_mesh);
        timer_ie = timer.Elapsed();
    }

    OOC_DEBUG("Inserting Vertices took " << timer_iv << "ms");
    OOC_DEBUG("Inserting Edges took " << timer_ie << "ms");
}

/* private */ void moist::MeshSlice::InsertVertices(const geo::Mesh& edge_mesh)
{
    //OOC_DEBUG("inserting " << edge_mesh.NbPoints() << " vertices");
    //for (const auto point : edge_mesh.Points())
    for (const auto v : edge_mesh.vertices)
    {
        this->InsertVertex(moist::exact::Point(edge_mesh.vertices.point(v)));
        _created_exact_cell_ids.clear();
    }

#ifndef NDEBUG
    this->_overlay.DebugMesh("exact_mesh_after_insertion.mesh");
#endif // NDEBUG
}

/* private */ void moist::MeshSlice::InsertVertex(const moist::exact::Point& point)
{
    const auto& grid_cells = _overlay.Grid()->GetCells(moist::create_point_box2d_exact(point));
    std::size_t v_1to2 = moist::exact::NO_VERTEX;

    bool dbg_cell = false;
    // 62.243751 10.2449047 10
    if (is_equal_test(point.x(), 62.243751) && is_equal_test(point.y(), 10.2449047))
    {
        dbg_cell = true;
    }

    if (is_equal_test(point.x(), 62.9993) && is_equal_test(point.y(), 11))
    {
        dbg_cell = true;
    }

    // 62 10.0013
    if (is_equal_test(point.x(), 62) && is_equal_test(point.y(), 10.0013))
    {
        dbg_cell = true;
    }

    for (const auto grid_cell : grid_cells)
    {
        if (!_overlay.Grid()->_grid.contains(grid_cell)) // cringe
        {
            continue;
        }

        if (dbg_cell)
        {
            moist::ExactMesh dbg_cells;
            const auto cells = _overlay.Grid()->GetMeshCells(grid_cell);
            for (const auto c : cells)
            {
                const auto cell = this->_overlay.Cell(c);
                if (cell._deleted) continue;
                dbg_cells.Add(moist::exact::Cell(
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
                ), true);
            }

            dbg_cells.Add(moist::exact::Cell(
                dbg_cells.Add(moist::exact::Point(point)),
                dbg_cells.Add(moist::exact::Point(point)),
                dbg_cells.Add(moist::exact::Point(point)),
                dbg_cells.Add(moist::exact::Point(point))
            ), true);
            dbg_cells.Add(moist::exact::Point(point));
            dbg_cells.DebugMesh("dbg_cell_" + std::to_string(grid_cell.first) + "_" + std::to_string(grid_cell.second) + ".msh");
        }

        const auto nb_cells = _overlay.Grid()->GetMeshCells(grid_cell).size();
        for (const std::size_t& c : _overlay.Grid()->GetMeshCells(grid_cell))
        {
            const auto& cell = this->_overlay.Cell(c);
            if (dbg_cell)
            {
                moist::ExactMesh dbg_cells;
                dbg_cells.Add(moist::exact::Cell(
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
                ), true);
                dbg_cells.Add(moist::exact::Cell(
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point))
                ), true);
                //dbg_cells.DebugMesh("dbg_cell_" + std::to_string(c) + ".msh");
            }

            if (c == 317)
            {
                int x = 43242;
            }

            const auto location = moist::new_predicates::point_in_tet_exact(this->_overlay, c, point, true);
            if (cell._deleted)
            {
                continue;
            }

            if (location == moist::predicates::PointInTet::FACET)
            {
                const std::size_t v = this->_overlay.Add(point);
                moist::operation::exact::InsertVertexOnCellBoundaryFacet(c, v, this->_overlay);
            }
            else if (location == moist::predicates::PointInTet::EDGE)
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = this->_overlay.Add(point);
                }

                if (!moist::operation::exact::InsertVertexOnCellBoundaryEdgeOld(c, v_1to2, this->_overlay))
                {
                                moist::ExactMesh dbg_cells;
            const auto cells = _overlay.Grid()->GetMeshCells(grid_cell);
            for (const auto c : cells)
            {
                const auto cell = this->_overlay.Cell(c);
                if (cell._deleted) continue;
                dbg_cells.Add(moist::exact::Cell(
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                    dbg_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
                    ), true);
                }

                dbg_cells.Add(moist::exact::Cell(
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point)),
                    dbg_cells.Add(moist::exact::Point(point))
                ), true);
                dbg_cells.Add(moist::exact::Point(point));
                dbg_cells.DebugMesh("dbg_cell_" + std::to_string(grid_cell.first) + "_" + std::to_string(grid_cell.second) + ".msh");
                OOC_DEBUG("point is " << point._p);
                }
            }
            /*else if (location == moist::predicates::PointInTet::EDGE_APPROXIMATED)
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = this->_overlay.Add(point);
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, this->_overlay, 1e-14);
            }*/
        }
    }
}



/* private */ void moist::MeshSlice::InsertEdgesPrivate(const geo::Mesh& edge_mesh)
{
    const auto edges = moist::geometry::collect_edges(edge_mesh);
    std::size_t i = 0;
    for (const auto edge : edges)
    {
        this->_created_exact_cell_ids.clear();
        this->InsertEdge(moist::exact::EdgePoints { edge_mesh.vertices.point(edge.v0), edge_mesh.vertices.point(edge.v1) });

        if (/*i % 25 == 0 ||*/ (i == 182))
        {
            this->_overlay.DebugMesh("dbg_edge_insertion_" + std::to_string(i) + ".msh");
        }
        i++;
    }

    this->_overlay.DebugMesh("exact_mesh_after_edge_insertion.mesh");
}

/* private */ void moist::MeshSlice::InsertEdge(const moist::exact::EdgePoints& edge)
{
    std::unordered_set<std::size_t> edge_vertices;
    std::vector<moist::exact::Cell> created_cells;

    moist::ExactMesh dbg_crossed_cells;
    moist::ExactMesh dbg_created_cells;

    const auto grid_cells = _overlay.Grid()->GetCells(moist::create_edge_box2d_exact(edge));
    int i = 0;
    for (const auto grid_cell : grid_cells)
    {
        if (!_overlay.Grid()->_grid.contains(grid_cell))
        {
            continue;
        }

        for (const std::size_t c : _overlay.Grid()->GetMeshCells(grid_cell))
        {
            const auto cell = _overlay.Cell(c);
            if (cell._deleted || moist::geometry::exact::is_degenerate(cell))
            {
                continue;
            }

            auto intersected_edges = moist::operation::exact::FindIntersectedEdges(edge, c, _overlay);
            for (auto& intersected_edge : intersected_edges)
            {
                auto& point = intersected_edge.p;
                const auto& existing = std::find_if(edge_vertices.begin(), edge_vertices.end(), [&](const std::size_t v)
                {
                    return _overlay.Point(v) == point;
                });

                if (existing == edge_vertices.end())
                {
                    const auto& v = _overlay.Add(point);
                    edge_vertices.insert(v);
                    intersected_edge.vp = v;
                }
                else
                {
                    intersected_edge.vp = *existing;
                }
            }

            if (intersected_edges.size() > 0)
            {
                dbg_crossed_cells.Add(moist::exact::Cell(
                    dbg_crossed_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                    dbg_crossed_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                    dbg_crossed_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                    dbg_crossed_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
                ), true);
            }

            switch (intersected_edges.size())
            {
                case 0:
                    break;
                case 1:
                    moist::operation::exact::SplitEdge1_2(c, intersected_edges.at(0), _overlay, created_cells);
                    break;
                case 2:
                    moist::operation::exact::SplitEdge1_3(c, intersected_edges.at(0), intersected_edges.at(1), _overlay, created_cells);
                    break;
                default:
                    OOC_ERROR("invalid amount of edge intersections: " << intersected_edges.size());
                    break;
            }

            if (!created_cells.empty() && !intersected_edges.empty())
            {
                moist::ExactMesh dbg_during_creation;
                for (const auto& cell : created_cells)
                {
                    if (cell._deleted) continue;
                    dbg_during_creation.Add(moist::exact::Cell(
                        dbg_during_creation.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                        dbg_during_creation.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                        dbg_during_creation.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                        dbg_during_creation.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
                    ), true);
                }

                // dbg_during_creation.DebugMesh("dbg_during_creation.msh");
            }
        }
    }

    for (const auto cell : created_cells)
    {
        if (cell._deleted || moist::geometry::exact::is_degenerate(cell)) continue;
        dbg_created_cells.Add(moist::exact::Cell(
            dbg_created_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
           dbg_created_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
            dbg_created_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
            dbg_created_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
        ), true);

        _created_exact_cell_ids.push_back(this->_overlay.Add(cell));
    }

    if (!created_cells.empty())
    {
        dbg_created_cells.Add(moist::exact::Point(edge.p0));
        dbg_created_cells.Add(moist::exact::Point(edge.p1));
        dbg_crossed_cells.Add(moist::exact::Point(edge.p0));
        dbg_crossed_cells.Add(moist::exact::Point(edge.p1));
    }

    auto e0x = edge.p0.x();
    auto e0y = edge.p0.y();
    auto e1x = edge.p1.x();
    auto e1y = edge.p1.y();

    double ppx, ppy;
    if (!edge_vertices.empty())
    {
        auto first = *edge_vertices.begin();
        ppx = _overlay.Point(first).x();
        ppy = _overlay.Point(first).y();
    }
    //dbg_crossed_cells.DebugMesh("crossed.msh");
    //dbg_created_cells.DebugMesh("created.msh");
    this->DecimateEdges(edge, edge_vertices);

    moist::ExactMesh dbg_decimated_cells;
    if (!created_cells.empty())
    {
        for (const auto& c : _created_exact_cell_ids)
        {
            const auto cell = _overlay.Cell(c);
            if (cell._deleted) continue;
            dbg_decimated_cells.Add(moist::exact::Cell(
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
            ), true);
        }

        //dbg_decimated_cells.DebugMesh("dbg_after_decimation.msh");
    }
}

static std::array<std::size_t, 2> find_edge_vertices(const std::vector<std::size_t> cells, const moist::exact::EdgePoints& edge, moist::ExactMesh& mesh)
{
    std::array<std::size_t, 2> points
    {
        moist::exact::NO_VERTEX,
        moist::exact::NO_VERTEX
    };

    double nearest_e0 = std::numeric_limits<double>::max();
    double nearest_e1 = std::numeric_limits<double>::max();;

    for (const auto& c : cells)
    {
        if (nearest_e0 == -1.0 && nearest_e1 == -1.0)
        {
            break;
        }

        const auto cell = mesh.Cell(c);
        if (cell._deleted || moist::geometry::exact::is_degenerate(cell)) { continue; }

        for (const auto v : cell._points)
        {
            const auto point = mesh.Point(v);
            const bool eq_p0 = moist::geometry::exact::points_are_equal(point._p, edge.p0._p);
            const bool eq_p1 = moist::geometry::exact::points_are_equal(point._p, edge.p1._p);
            if (eq_p0)
            {
                nearest_e0 = -1.0;
                points[0] = v;
            }

            if (eq_p1)
            {
                nearest_e1 = -1.0;
                points[1] = v;
            }

            const double d_e0 = std::sqrt(CGAL::to_double(CGAL::squared_distance(point._p, edge.p0._p)));
            const double d_e1 = std::sqrt(CGAL::to_double(CGAL::squared_distance(point._p, edge.p1._p)));

            if (d_e0 < nearest_e0)
            {
                nearest_e0 = d_e0;
                points[0] = v;
            }

            if (d_e1 < nearest_e1)
            {
                nearest_e1 = d_e1;
                points[1] = v;
            }
        }
    }

    if (points[0] == moist::exact::NO_VERTEX || points[1] == moist::exact::NO_VERTEX)
    // if (std::find(points.begin(), points.end(), moist::exact::NO_VERTEX))
    {
        OOC_ERROR("edge endpoints are missing in crossed cells - p0 = " << edge.p0._p << ", p1 = " << edge.p1._p);
    }

    if (nearest_e0 != -1.0)
    {
        OOC_DEBUG("edge endpoint was not found in mesh - filling onto " << points[0]);
        mesh.Point(points[0])._fixed = true;
    }

    if (nearest_e1 != -1.0)
    {
        OOC_DEBUG("edge endpoint was not found in mesh - filling onto " << points[1]);
        mesh.Point(points[1])._fixed = true;
    }

    return points;
}

static std::unordered_set<std::size_t> find_vertex_adjacent_neighbours(const std::size_t v, const std::vector<std::size_t> cells, const std::unordered_set<std::size_t>& vertices, const moist::ExactMesh& mesh)
{
    std::unordered_set<std::size_t> neighbours;
    for (const auto c : cells)
    {
        const auto cell = mesh.Cell(c);
        // Exclude any cells that do not contain v, we only want direct neighbours (i.e. vertices that are in vertices set, and are incident to v in at least one cell)!
        if (!std::find_if(cell._points.begin(), cell._points.end(), [&](const std::size_t cv) { return cv == v; }))
        {
            continue;
        }

        for (const auto cv : cell._points)
        {
            if (vertices.contains(cv) && cv != v)
            {
                neighbours.insert(cv);
            }
        }
    }

    return neighbours;
}

/* private */ void moist::MeshSlice::DecimateEdges(const moist::exact::EdgePoints& edge, const std::unordered_set<std::size_t>& edge_vertices)
{
    if (edge_vertices.empty())
    {
        return;
    }

    std::unordered_set<std::size_t> vertices;
    for (const auto& v : edge_vertices)
    {
        vertices.insert(v);
    }

    // Find edge vertices in cells, and add them to the cluster...
    const auto endpoints = find_edge_vertices(_created_exact_cell_ids, edge, _overlay);
    if (endpoints[0] == moist::exact::NO_VERTEX && endpoints[1] == moist::exact::NO_VERTEX)
    {
        return;
    }
    else if (endpoints[0] == moist::exact::NO_VERTEX || endpoints[1] == moist::exact::NO_VERTEX)
    {
        return;
    }

    for (const auto endpoint : endpoints)
    {
        vertices.insert(endpoint);
    }

    int i = 0;
    // For each vertex, attempt to move it (and any vertices clustered onto it) onto a neighbour it is attached to by an edge shared by one or multiple cells.
    for (const auto v : edge_vertices)
    {
        bool decimatable = false;
        const auto neighbours = find_vertex_adjacent_neighbours(v, _created_exact_cell_ids, vertices, _overlay);
        const auto vp = _overlay.Point(v);

        for (const auto endpoint : endpoints)
        {
            if (neighbours.contains(endpoint) && this->CanMoveVertex(v, endpoint))
            {
                for (const std::size_t c : _created_exact_cell_ids)
                {
                    auto& cell = _overlay.Cell(c);
                    for (std::size_t lv = 0; lv < 4; lv++)
                    {
                        if (cell._points[lv] == v)
                        {
                            cell._points[lv] = endpoint;
                        }
                    }
                }

                decimatable = true;
                goto dbg;
            }
        }

        for (const std::size_t neighbour : neighbours)
        {
            if (this->CanMoveVertex(v, neighbour))
            {
                for (const auto c : _created_exact_cell_ids)
                {
                    auto& cell = _overlay.Cell(c);
                    for (std::size_t lv = 0; lv < 4; lv++)
                    {
                        if (cell._points[lv] == v)
                        {
                            cell._points[lv] = neighbour;
                        }
                    }
                }

                decimatable = true;
                goto dbg;
            }
        }

        if (!decimatable)
        {
            auto& p = _overlay.Point(v);
            p._fixed = true;

            OOC_DEBUG(std::setprecision(18) << "inserting non decimatable point " << v << " at [" << p.x() << ", " << p.y() << ", " << p.z() << "]");
            bool same_e0 = p._p == edge.p0._p;
            bool same_e1 = p._p == edge.p1._p;
            OOC_DEBUG(std::setprecision(50) << p._p);
            OOC_DEBUG(std::setprecision(50) << edge.p0._p);
            OOC_DEBUG(std::setprecision(50) << edge.p1._p);
        }

dbg:
        moist::ExactMesh dbg_decimated_cells;
        for (const auto& c : _created_exact_cell_ids)
        {
            const auto cell = _overlay.Cell(c);
            if (cell._deleted) continue;
            dbg_decimated_cells.Add(moist::exact::Cell(
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[0]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[1]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[2]))),
                dbg_decimated_cells.Add(moist::exact::Point(_overlay.Point(cell._points[3])))
            ), true);
        }
        //dbg_decimated_cells.DebugMesh("dbg_during_decimation" + std::to_string(i) + ".msh");
        i++;
    }

    for (const auto c : _created_exact_cell_ids)
    {
        std::unordered_set<std::size_t> elements;
        auto cell = _overlay.Cell(c);
        if (cell._deleted) { continue; }

        for (const auto v : cell._points)
        {
            if (elements.contains(v))
            {
                cell._deleted = true;
                break;
            }
            elements.insert(v);
        }

        const auto sign = CGAL::sign(CGAL::volume(
            _overlay.Point(cell._points[0])._p,
            _overlay.Point(cell._points[1])._p,
            _overlay.Point(cell._points[2])._p,
            _overlay.Point(cell._points[3])._p
        ));

        if (sign == CGAL::Sign::ZERO)
        {
            cell._deleted = true;
        }
    }
}

/* private */ bool moist::MeshSlice::CanMoveVertex(const std::size_t& v_from, const std::size_t& v_to)
{
    for (std::size_t c = 0; c < this->_created_exact_cell_ids.size(); c++)
    {
        const auto cell = _overlay.Cell(_created_exact_cell_ids[c]);
        if (cell._deleted || moist::geometry::exact::is_degenerate(cell))
        {
            continue;
        }

        const auto before = CGAL::sign(CGAL::volume(
            _overlay.Point(cell._points[0])._p,
            _overlay.Point(cell._points[1])._p,
            _overlay.Point(cell._points[2])._p,
            _overlay.Point(cell._points[3])._p
        ));

        const auto after = CGAL::sign(CGAL::volume(
            _overlay.Point(cell._points[0] == v_from ? v_to : cell._points[0])._p,
            _overlay.Point(cell._points[1] == v_from ? v_to : cell._points[1])._p,
            _overlay.Point(cell._points[2] == v_from ? v_to : cell._points[2])._p,
            _overlay.Point(cell._points[3] == v_from ? v_to : cell._points[3])._p
        ));

        if (after != CGAL::Sign::ZERO && before != CGAL::Sign::ZERO && before != after)
        {
            return false;
        }
    }

    return true;
}

/* public */ void moist::MeshSlice::GetFixedGeometry(moist::ExactMesh& mesh, const moist::AxisAlignedPlane& plane)
{
    std::unordered_map<std::size_t, std::size_t> vertices;
    for (const auto cell : _overlay.Cells())
    {
        // find the fixed vertices in each cell...
        if (std::none_of(cell._points.begin(), cell._points.end(), [&](const g_index v) { return _overlay.Point(v)._fixed; }))
        {
            continue;
        }

        for (const auto& [i, j] : moist::geometry::exact::TET_EDGE_DESCRIPTOR)
        {
            const auto v0 = cell._points[i];
            const auto v1 = cell._points[j];

            // if the edge lies on the interface plane, we add it to the edge mesh...
            const moist::exact::EdgePoints edge_points { _overlay.Point(v0), _overlay.Point(v1) };
            if (!edge_points.IsInterface())
            {
                continue;
            }

            if (edge_points.p0._fixed || edge_points.p1._fixed)
            {
                if (!vertices.contains(v0))
                {
                    vertices.emplace(v0, mesh.Add(edge_points.p0));
                }

                if (!vertices.contains(v1))
                {
                    vertices.emplace(v1, mesh.Add(edge_points.p1));
                }

                mesh.Add(moist::exact::Edge(vertices[v0], vertices[v1]));
            }
        }
    }
}

#endif // MOIST_OPTION_EXACT_PREDICATES
