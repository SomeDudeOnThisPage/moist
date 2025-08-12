#include "merger.hpp"

#include <geogram/mesh/mesh_repair.h>

#include "moist/core/predicates.inl"
#include "moist/core/utils.hpp"

#include "arrangement.hpp"
#include "exact_types.hpp"
#include "local_operations.hpp"
#include "new_predicates.inl"
#include "geometry_exact.inl"

#include <tetgen.h>

static geo::Mesh dbg_triangles(3);
static void geo_to_tetgen(const geo::Mesh& geo_mesh, tetgenio& in)
{
    in.initialize();
    in.numberofpoints = (int) geo_mesh.vertices.nb();
    in.pointlist = new double[in.numberofpoints * 3];

    for (geo::index_t v = 0; v < geo_mesh.vertices.nb(); v++)
    {
        const double* p = geo_mesh.vertices.point_ptr(v);
        in.pointlist[3 * v] = p[0];
        in.pointlist[3 * v + 1] = p[1];
        in.pointlist[3 * v + 2] = p[2];
    }

    in.numberoftetrahedra = (int) geo_mesh.cells.nb();
    in.tetrahedronlist = new int[in.numberoftetrahedra * 4];

    for (geo::index_t c = 0; c < geo_mesh.cells.nb(); c++)
    {
        for (int lv = 0; lv < 4; lv++)
        {
            in.tetrahedronlist[4 * c + lv] = (int) geo_mesh.cells.vertex(c, lv);
        }
    }

    if (geo_mesh.vertices.attributes().is_defined("v_boundary_marker"))
    {
        //in.numberofpointattributes = 1;
        in.pointmarkerlist = new int[in.numberofpoints];
        auto v_boundary_marker = geo::Attribute<bool>(geo_mesh.vertices.attributes(), "v_boundary_marker");
        for (geo::index_t v = 0; v < geo_mesh.vertices.nb(); v++)
        {
            if (v_boundary_marker[v])
            {
                in.pointmarkerlist[v] = -1;
            }
        }
    }
}

static void tetgen_to_geo(tetgenio& tetgen_io, geo::Mesh& geo_mesh)
{
    geo_mesh.clear(false, false);
    for (int v = 0; v < tetgen_io.numberofpoints; v++)
    {
        geo_mesh.vertices.create_vertex(geo::vec3(
            tetgen_io.pointlist[3 * v],
            tetgen_io.pointlist[3 * v + 1],
            tetgen_io.pointlist[3 * v + 2]
        ));
    }

    for (int c = 0; c < tetgen_io.numberoftetrahedra; c++)
    {
        geo_mesh.cells.create_tet(
            tetgen_io.tetrahedronlist[4 * c],
            tetgen_io.tetrahedronlist[4 * c + 1],
            tetgen_io.tetrahedronlist[4 * c + 2],
            tetgen_io.tetrahedronlist[4 * c + 3]
        );
    }
}

moist::Merger::Merger(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane)
{
    this->ConstructExactMesh(a, _mesh_a, plane);
    this->ConstructExactMesh(b, _mesh_b, plane);
#ifndef NDEBUG
    _mesh_a.DebugMesh("exact_mesh_a.mesh");
    _mesh_b.DebugMesh("exact_mesh_b.mesh");

    moist::utils::geogram::save("crownless_a.msh", a);
    moist::utils::geogram::save("crownless_b.msh", b);
    moist::utils::geogram::save("crownless_a.mesh", a);
    moist::utils::geogram::save("crownless_b.mesh", b);
#endif // NDEBUG

    this->InsertBToA();
    this->HollowCrown();
    _crown.DebugMesh("crown.mesh");
    _crown_surface.DebugMesh("crown_surface.mesh");

    // check how the interface looks for debugging after all ops
    /*auto v_boundary_marker = geo::Attribute<bool>(_crown.vertices.attributes(), "v_boundary_marker");
    for (const geo::index_t v : _crown.vertices)
    {
        v_boundary_marker[v] = _crown.vertices.point(v).z == plane.extent;
    }

    tetgenio in;
    tetgenio out;
    tetgenbehavior behavior;
    //behavior.plc = 1;
    behavior.refine = 1;
    //behavior.quality = 1;
    behavior.coarsen = 1;
    //behavior.opt_scheme = 7;
    behavior.coarsen_percent = 0.1;

    geo_to_tetgen(_crown, in);
    tetrahedralize(&behavior, &in, &out);
    tetgen_to_geo(out, _crown);

    moist::utils::geogram::save("crown.msh", _crown);*/
}

void moist::Merger::HollowCrown()
{
    _crown.ResetGrid(10.0);
    // this is SUPER bad but for now it must be enough
    std::vector<moist::exact::Point> inserted_points;
    for (std::size_t c = 0; c < _crown.NbCells(); c++)
    {
        std::array<std::size_t, 4> facet_adjacencies {0, 0, 0, 0};
        const auto aabb = moist::create_cell_bbox2d_exact(c, _crown);
        const auto grid_cells = _crown.Grid()->GetCells(aabb);
        const auto cell = _crown.Cell(c);

        const std::array<moist::exact::Point, 4> points = {
            _crown.Point(cell[0]),
            _crown.Point(cell[1]),
            _crown.Point(cell[2]),
            _crown.Point(cell[3])
        };

        if (points[0] == points[1] || points[0] == points[2] || points[0] == points[3] || points[1] == points[2] || points[1] == points[3])
        {
            OOC_WARNING("Degenerate cell " << c << "!");
            _crown.DeleteCell(c);
            continue;
        }

        // check all cells in the vicinity for a shared facet. If it exists, increment the corresponding facet_adjacenties
        // since we don't have any internal mesh connectivity, just brute force this via our grid optimization
        std::unordered_set<std::size_t> checked{c};
        for (std::size_t lf = 0; lf < 4; lf++) // local facet of the original cell
        {
            bool found_another = false;
            const auto& lf_indices = moist::geometry::exact::TET_FACET_DESCRIPTOR.at(lf);
            const std::array<moist::exact::Point, 3> facet
            {
                points[lf_indices[0]],
                points[lf_indices[1]],
                points[lf_indices[2]]
            };

            for (const auto grid_cell : grid_cells)
            {
                if (!_crown.Grid()->_grid.contains(grid_cell))
                {
                    continue;
                }

                for (const std::size_t c_other : _crown.Grid()->GetMeshCells(grid_cell))
                {
                    if (checked.contains(c_other))
                    {
                        continue;
                    }

                    const auto& cell_other = _crown.Cell(c_other);
                    // if c_other contains all points of the facet we are currently checking, they match
                    bool found[4] { false, false, false, false };
                    for (const std::size_t v_other : cell_other._points)
                    {
                        const auto p_other = _crown.Point(v_other);
                        if (p_other == points[0])
                        {
                            found[0] = true;
                        }
                        else if (p_other == points[1])
                        {
                            found[1] = true;
                        }
                        else if (p_other == points[2])
                        {
                            found[2] = true;
                        }
                        else if (p_other == points[3])
                        {
                            found[3] = true;
                        }
                    }
                    /*for (std::size_t lv = 0; lv < 4; lv++)
                    {
                        const auto& point_other = _crown.Point(cell_other[lv]);
                        for (std::size_t lf_v = 0; lf_v < 3; lf_v++)
                        {
                            if (facet[lf_v] == point_other)
                            {
                                found[lf_v] = true;
                            }
                        }
                    }*/

                    checked.insert(c_other);

                    std::size_t nb_found = 0;
                    for (std::size_t i = 0; i < 4; i++)
                    {
                        if (found[i]) nb_found++;
                    }

                    if (nb_found >= 3)
                    {
                        found_another = true;
                        break;
                    }
                }

                if (found_another)
                {
                    break;
                }
            }

            // create LF in _crown_surface
            if (!found_another)
            {
                const auto& lf_indices = moist::geometry::exact::TET_FACET_DESCRIPTOR.at(lf);
                std::size_t f_points[3];
                for (std::size_t lv = 0; lv < 3; lv++)
                {
                    auto it = std::find(inserted_points.begin(), inserted_points.end(), points[lf_indices[lv]]);
                    if (it == inserted_points.end())
                    {
                        f_points[lv] = _crown_surface.Add(moist::exact::Point(points[lf_indices[lv]]));
                        inserted_points.push_back(_crown_surface.Point(f_points[lv]));
                    }
                    else
                    {
                        f_points[lv] = std::distance(inserted_points.begin(), it);
                    }
                }

                if (f_points[0] == f_points[1] || f_points[0] == f_points[2] || f_points[1] == f_points[2])
                {
                    OOC_DEBUG(points[0]._p);
                    OOC_DEBUG(points[1]._p);
                    OOC_DEBUG(points[2]._p);
                    OOC_DEBUG(points[3]._p);

                    OOC_ERROR("invalid triangle!");
                }
                _crown_surface.Add(moist::exact::Facet(
                    f_points[0],
                    f_points[1],
                    f_points[2]
                ));
            }
        }
    }
}

std::vector<std::size_t> shared_facet(const std::array<std::size_t, 4>& a, const std::array<std::size_t, 4>& b)
{
    std::vector<size_t> common;
    std::array<size_t, 4> b_copy = b; // so we can mark used elements

    for (size_t va : a) {
        auto it = std::find(b_copy.begin(), b_copy.end(), va);
        if (it != b_copy.end()) {
            common.push_back(va);
            *it = static_cast<size_t>(-1); // mark as used
        }
    }

    if (common.size() == 3) {
        return common; // exactly 3 shared points
    }
    return {}; // empty if not exactly 3
}

std::array<std::size_t, 2> longest_edge_on_interface(const moist::exact::Cell& cell, moist::ExactMesh& mesh)
{
    moist::exact::Kernel::FT distance = std::numeric_limits<double>::lowest();
    std::size_t min0, min1;
    for (const auto& [i, j] : moist::geometry::exact::TET_EDGE_DESCRIPTOR)
    {
        const auto& p0 = mesh.Point(cell._points[i]);
        const auto& p1 = mesh.Point(cell._points[j]);

        if (p0._v != moist::exact::NO_VERTEX || p1._v != moist::exact::NO_VERTEX)
        {
            continue;
        }

        const auto edge_length = CGAL::squared_distance(p0._p, p1._p);
        if (edge_length > distance)
        {
            distance = edge_length;
            min0 = cell._points[i];
            min1 = cell._points[j];
        }
    }

    return {min0, min1};
}

void moist::Merger::Prune(moist::ExactMesh& mesh, const double min_volume)
{
    // for each tet, with a face on the interface:
    // then, for each tet, check how many faces are free.
    // if it has at least one more free non-interface face, and bad volume, we can prune it.
    // this ensures we don't prune non-boundary facets we need to keep for connectivity.

    std::size_t pruned = 0;
    for (std::size_t c = 0; c < mesh.NbCells(); c++)
    {
        auto& cell = mesh.Cell(c);
        auto longest_interface_edge = longest_edge_on_interface(cell, mesh);
        if (cell._deleted || !moist::geometry::exact::cell_has_interface_facet(cell, mesh))
        {
            // we cannot prune cells without a face on the interface, since they may be required for connectivity inside the mesh
            continue;
        }

        std::unordered_set<std::size_t> neighbours;
        bool has_longest_edge_neighbour = false;
        const auto bbox = moist::create_interface_cell_bbox2d_exact(c, mesh);
        const auto grid_cells = mesh.Grid()->GetCells(bbox);
        for (const auto grid_cell : grid_cells)
        {
            if (!mesh.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
            {
                continue;
            }

            for (const auto c_other : mesh.Grid()->GetMeshCells(grid_cell))
            {
                const auto& cell_other = mesh.Cell(c_other);
                if (c_other == c || cell_other._deleted || neighbours.contains(c_other))
                {
                    continue;
                }

                const auto face = shared_facet(cell._points, cell_other._points);
                if (!face.empty())
                {
                    if (std::find(face.begin(), face.end(), longest_interface_edge[0]) != face.end() ||  std::find(face.begin(), face.end(), longest_interface_edge[1]) != face.end())
                    {
                        has_longest_edge_neighbour = true;
                    }
                    neighbours.insert(c_other);
                }
            }
        }

        if (neighbours.size() == 1 || (neighbours.size() == 2 && has_longest_edge_neighbour)) // we can check 2 here, since we made sure it's not wharing its longest facet!
        {
            const auto volume = std::abs(CGAL::to_double(CGAL::volume(
                mesh.Point(cell._points[0])._p,
                mesh.Point(cell._points[1])._p,
                mesh.Point(cell._points[2])._p,
                mesh.Point(cell._points[3])._p
            )));

            if (volume < min_volume)
            {
                cell._deleted = true;
                pruned++;
            }
        }
full_continue:
    }

    OOC_DEBUG("pruned " << pruned << " cells");
}

/**
 * @brief Preprocess rounding to avoid extremely small triangles in overlaps during processing (i.e. edge lengths of e-20 or less)
 *
 * @param value
 * @return double
 */
static double round16(double value)
{
    double rounded = std::round(value * 1e8) / 1e8;
    if (rounded == 0.0)
    {
        rounded = 0.0; // don't ask
    }
    return rounded;
}

void moist::Merger::ConstructExactMesh(geo::Mesh& mesh, moist::ExactMesh& exact_mesh, const moist::AxisAlignedPlane& plane)
{
    std::unordered_map<geo::index_t, std::size_t> added_vertices;
    geo::vector<geo::index_t> to_delete(mesh.cells.nb());
    for (const geo::index_t c : mesh.cells)
    {
        if (!moist::predicates::cell_on_plane(c, mesh, plane))
        {
            continue;
        }

        std::size_t points[4];
        for (geo::index_t lv = 0; lv < 4; lv++)
        {
            const geo::index_t v = mesh.cells.vertex(c, lv);
            if (added_vertices.contains(v))
            {
                points[lv] = added_vertices[v];
                continue;
            }
            else
            {
                // here we must keep _other of the point always on NO_VERTEX...
                auto point = mesh.vertices.point(v);
                bool interface = moist::predicates::point_on_plane(point, plane);
                if (interface)
                {
                    point[2] = plane.extent;
                }

                vec3 point_rounded = vec3(round16(point[0]), round16(point[1]), round16(point[2]));
                points[lv] = exact_mesh.Add(moist::exact::Point(point_rounded, interface ? geo::NO_VERTEX : v));
                added_vertices[v] = points[lv];
            }
        }

        exact_mesh.Add(moist::exact::Cell(points[0], points[1], points[2], points[3]), true);
        to_delete[c] = true;
    }
    mesh.cells.delete_elements(to_delete);

    // edge length does not make the most sense, since it must be based relative to the overall length of our bbox diag.
    // dumb hard coded value... use nb of cells on the interface maybe... times some factor controlled so we can evaluate performance based on it...
    constexpr double FACTOR = 0.5;
    exact_mesh.ResetGrid(/*FACTOR * std::sqrt(exact_mesh.NbCells())*/ 10.0);
}

void moist::Merger::InsertAndMapPoints(const moist::ExactMesh& from, moist::ExactMesh& to)
{
    for (std::size_t v = 0; v < from.NbPoints(); v++)
    {
        const auto& point = from.Point(v);
        if (point._other != moist::exact::NO_VERTEX)
        {
            continue;
        }

        this->InsertPoint(v, from, to);
    }
}

static bool point_in_tet_on_edge(const moist::predicates::PointInTet& pit)
{
    return pit == moist::predicates::PointInTet::EDGE01
           || pit == moist::predicates::PointInTet::EDGE02
           || pit == moist::predicates::PointInTet::EDGE03
           || pit == moist::predicates::PointInTet::EDGE12
           || pit == moist::predicates::PointInTet::EDGE13
           || pit == moist::predicates::PointInTet::EDGE23;
}

// This is done in the "old" meshes, not in the one being currently constructed, because tets with an edge on the interface may be cut
// multiple times.
void moist::Merger::InsertPointIntoEdges(const moist::exact::Point& point, moist::ExactMesh& mesh)
{
    const auto& grid_cells = mesh.Grid()->GetCells(moist::create_point_box2d_exact(point));
    std::size_t v_1to2 = moist::exact::NO_VERTEX;

    bool found = false;
    std::unordered_set<std::size_t> inserted;
    for (const auto grid_cell : grid_cells)
    {
        if (!mesh.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
        {
            continue;
        }

        const auto cells = mesh.Grid()->GetMeshCells(grid_cell);
        for (const auto c : cells)
        {
            const auto& cell = mesh.Cell(c);
            if (cell._deleted || inserted.contains(c))
            {
                continue;
            }

            std::size_t nb_interface = 0;
            for (std::size_t lv = 0; lv < 4; lv++)
            {
                if (mesh.Point(cell._points[lv])._v != moist::exact::NO_VERTEX)
                {
                    nb_interface++;
                }
            }

            if (nb_interface != 2)
            {
                continue;
            }

            const auto location = moist::new_predicates::point_in_tet_exact(mesh, c, point, true);
            if (point_in_tet_on_edge(location))
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = mesh.Add(moist::exact::Point(point));
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, mesh);
                found = true;
                inserted.insert(c);
            }
        }
    }

    if (!found)
    {
        //OOC_DEBUG("Oh noes! " << point._p);
    }
    else
    {
        //OOC_DEBUG("Oh yese! " << point._p);
    }
}

void moist::Merger::InsertPoint(const std::size_t v, const moist::ExactMesh& from, moist::ExactMesh& to)
{
    // check where to insert the point, and update to accordingly. also map the point!
    const auto point = from.Point(v);
    const auto& grid_cells = to.Grid()->GetCells(moist::create_point_box2d_exact(point));
    std::size_t v_1to2 = moist::exact::NO_VERTEX;

    for (const auto grid_cell : grid_cells)
    {
        if (!to.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
        {
            continue;
        }

        const auto cells = to.Grid()->GetMeshCells(grid_cell);
        for (const auto c : cells)
        {
            const auto& cell = to.Cell(c);
            if (cell._deleted)
            {
                continue;
            }

            const auto location = moist::new_predicates::point_in_tet_exact(to, c, point, true);
            if (location == moist::predicates::PointInTet::FACET)
            {
                const std::size_t v_new = to.Add(moist::exact::Point(point));
                moist::operation::exact::InsertVertexOnCellBoundaryFacet(c, v_new, to);
                to.Point(v_new)._other = v;
            }
            else if (point_in_tet_on_edge(location))
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = to.Add(moist::exact::Point(point));
                    to.Point(v_1to2)._other = v;
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, to);
            }
        }
    }
}

static geo::Mesh dbg_full(3);
void moist::Merger::InsertBToA()
{
    for (std::size_t c = 0; c < _mesh_b.NbCells(); c++)
    {
        const auto cell = _mesh_b.Cell(c);
        if (cell._deleted)
        {
            continue;
        }

        if (!moist::geometry::exact::cell_has_interface_facet(cell, _mesh_b))
        {
            continue;
        }

        this->InsertCellBToA(c);

        if (c % 50 == 0)
        {
            //geo::mesh_repair(dbg_full, geo::MeshRepairMode::MESH_REPAIR_DEFAULT, 1e-14);
            moist::utils::geogram::save("dbg_full_" + std::to_string(c) + ".mesh", dbg_full);
            moist::utils::geogram::save("dbg_triangles_" + std::to_string(c) + ".mesh", dbg_triangles);
        }
    }

    for (std::size_t c = 0; c < _mesh_a.NbCells(); c++)
    {
        const auto& cell = _mesh_a.Cell(c);
        if (cell._deleted)
        {
            continue;
        }

        std::size_t nb_non_interface = 0;
        for (std::size_t lv = 0; lv < 4; lv++)
        {
            if (_mesh_a.Point(cell._points[lv])._v != moist::exact::NO_VERTEX)
            {
                nb_non_interface++;
            }
        }

        if (nb_non_interface == 1)
        {
            continue;
        }

        const auto p0 = _mesh_a.Point(cell._points[0]);
        const auto p1 = _mesh_a.Point(cell._points[1]);
        const auto p2 = _mesh_a.Point(cell._points[2]);
        const auto p3 = _mesh_a.Point(cell._points[3]);

        dbg_full.cells.create_tet(
            dbg_full.vertices.create_vertex(vec3(p0.x(), p0.y(), p0.z())),
            dbg_full.vertices.create_vertex(vec3(p1.x(), p1.y(), p1.z())),
            dbg_full.vertices.create_vertex(vec3(p2.x(), p2.y(), p2.z())),
            dbg_full.vertices.create_vertex(vec3(p3.x(), p3.y(), p3.z()))
        );

        _crown.Add(moist::exact::Cell(
            _crown.Add(moist::exact::Point(p0)),
            _crown.Add(moist::exact::Point(p1)),
            _crown.Add(moist::exact::Point(p2)),
            _crown.Add(moist::exact::Point(p3))
        ), true);
    }

    for (std::size_t c = 0; c < _mesh_b.NbCells(); c++)
    {
        const auto& cell = _mesh_b.Cell(c);
        if (cell._deleted)
        {
            continue;
        }

        std::size_t nb_non_interface = 0;
        for (std::size_t lv = 0; lv < 4; lv++)
        {
            if (_mesh_b.Point(cell._points[lv])._v != moist::exact::NO_VERTEX)
            {
                nb_non_interface++;
            }
        }

        if (nb_non_interface == 1)
        {
            continue;
        }

        const auto p0 = _mesh_b.Point(cell._points[0]);
        const auto p1 = _mesh_b.Point(cell._points[1]);
        const auto p2 = _mesh_b.Point(cell._points[2]);
        const auto p3 = _mesh_b.Point(cell._points[3]);

        _crown.Add(moist::exact::Cell(
            _crown.Add(moist::exact::Point(p0)),
            _crown.Add(moist::exact::Point(p1)),
            _crown.Add(moist::exact::Point(p2)),
            _crown.Add(moist::exact::Point(p3))
        ), true);

        dbg_full.cells.create_tet(
            dbg_full.vertices.create_vertex(vec3(p0.x(), p0.y(), p0.z())),
            dbg_full.vertices.create_vertex(vec3(p1.x(), p1.y(), p1.z())),
            dbg_full.vertices.create_vertex(vec3(p2.x(), p2.y(), p2.z())),
            dbg_full.vertices.create_vertex(vec3(p3.x(), p3.y(), p3.z()))
        );
    }

    geo::mesh_repair(dbg_full, geo::MeshRepairMode::MESH_REPAIR_DEFAULT, 1e-14);
    moist::utils::geogram::save("dbg_crown.msh", dbg_full);
    moist::utils::geogram::save("dbg_crown.mesh", dbg_full);
}

void moist::Merger::InsertCellBToA(const std::size_t cb)
{
    // since grid cells can overlap in which cells they contain, we need to track into which cells we already inserted, so as to not
    // insert a cell twice!
    std::unordered_set<std::size_t> inserted;

    auto cell = _mesh_b.Cell(cb);
    const auto grid_cells = _mesh_a.Grid()->GetCells(moist::create_interface_cell_bbox2d_exact(cb, _mesh_b));

    for (const auto grid_cell : grid_cells)
    {
        if (!_mesh_a.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
        {
            continue;
        }

        for (const auto ca : _mesh_a.Grid()->GetMeshCells(grid_cell))
        {
            if (_mesh_a.Cell(ca)._deleted || inserted.contains(ca))
            {
                continue;
            }

            if (!moist::geometry::exact::cell_has_interface_facet(_mesh_a.Cell(ca), _mesh_a))
            {
                continue;
            }

            this->InsertCellIntoCell(ca, cb);
            inserted.emplace(ca);
        }
    }
}

static moist::exact::Triangle get_interface_triangle(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    std::vector<moist::exact::Point> points;
    for (const auto v : cell._points)
    {
        const auto point = mesh.Point(v);
        if (point._v == moist::exact::NO_VERTEX)
        {
            points.push_back(point);
        }
    }

    if (points.size() != 3)
    {
        OOC_ERROR("invalid interface cell triangle configuration");
    }

    return moist::exact::Triangle(points.at(0), points.at(1), points.at(2));
}

static std::size_t get_opposite_vertex(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    for (const auto v : cell._points)
    {
        const auto point = mesh.Point(v);
        if (point._v != moist::exact::NO_VERTEX)
        {
            return v;
        }
    }

    OOC_ERROR("invalid opposite vertex");
    return -1U;
}

static moist::exact::VertexCorrespondence get_correspondence(const moist::exact::Triangle& triangle)
{
    auto correspondence = moist::exact::VertexCorrespondence::AB;

    for (std::size_t i = 0; i < 3; i++)
    {
        if (triangle._points.at(i)._correspondence != moist::exact::VertexCorrespondence::AB)
        {
            if (correspondence != moist::exact::VertexCorrespondence::AB && triangle._points.at(i)._correspondence != correspondence)
            {
                OOC_ERROR("failed to determine vertex correspondence");
            }

            correspondence = triangle._points.at(i)._correspondence; // technically we can return here already
        }
    }

    return correspondence;
}

static bool triangle_has_point(const moist::exact::Triangle& triangle, const moist::exact::Kernel::Point_3& point)
{
    for (std::size_t lv = 0; lv < 3; lv++)
    {
        if (triangle._t.vertex(lv) == point)
        {
            return true;
        }
    }

    return false;
};

static bool marked_for_debug = false;
void moist::Merger::InsertCellIntoCell(const std::size_t ca, const std::size_t cb)
{
    auto& cell_a = _mesh_a.Cell(ca);
    auto& cell_b = _mesh_b.Cell(cb);

    const auto t0 = get_interface_triangle(cell_a, _mesh_a);
    const auto t1 = get_interface_triangle(cell_b, _mesh_b);
    const std::size_t v_opposite_a = get_opposite_vertex(cell_a, _mesh_a);
    const std::size_t v_opposite_b = get_opposite_vertex(cell_b, _mesh_b);

    if (!moist::exact::arrangeable(t0, t1))
    {
        //return;
    }

    const auto triangulation = moist::exact::arrange(t0, t1);

    geo::Mesh dbg(3);
    for (const auto triangle : triangulation)
    {
        const auto p_opposite_a = _mesh_a.Point(v_opposite_a);
        const auto p_opposite_b = _mesh_b.Point(v_opposite_b);

        // check which mesh(es) the triangle belongs to:
        // 1. if all points are AB, connect in both
        // 2. if some or no points are AB, connect to whatever the other points are
        // 3. if we have A and B in one point, fail, but this should not happen unless the CGAL routine messes up
        const auto correspondence = get_correspondence(triangle);
        if (correspondence == moist::exact::VertexCorrespondence::AB)
        {
            for (std::size_t lv = 0; lv < 3; lv++)
            {
                if (!triangle_has_point(t0, triangle._t.vertex(lv)))
                {
                    this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_a);
                }

                if (!triangle_has_point(t1, triangle._t.vertex(lv)))
                {
                    this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_b);
                }
            }

        #ifndef NDEBUG
            const vec3 vt0(CGAL::to_double(triangle._t.vertex(0).x()), CGAL::to_double(triangle._t.vertex(0).y()), CGAL::to_double(triangle._t.vertex(0).z()));
            const vec3 vt1(CGAL::to_double(triangle._t.vertex(1).x()), CGAL::to_double(triangle._t.vertex(1).y()), CGAL::to_double(triangle._t.vertex(1).z()));
            const vec3 vt2(CGAL::to_double(triangle._t.vertex(2).x()), CGAL::to_double(triangle._t.vertex(2).y()), CGAL::to_double(triangle._t.vertex(2).z()));
            dbg_triangles.facets.create_triangle(
                dbg_triangles.vertices.create_vertex(vt0),
                dbg_triangles.vertices.create_vertex(vt1),
                dbg_triangles.vertices.create_vertex(vt2)
            );
        #endif // NDEBUG

            _crown.Add(moist::exact::Cell(
                _crown.Add(moist::exact::Point(triangle._t.vertex(0))),
                _crown.Add(moist::exact::Point(triangle._t.vertex(1))),
                _crown.Add(moist::exact::Point(triangle._t.vertex(2))),
                _crown.Add(moist::exact::Point(p_opposite_a))
            ), true);

            _crown.Add(moist::exact::Cell(
                _crown.Add(moist::exact::Point(triangle._t.vertex(0))),
                _crown.Add(moist::exact::Point(triangle._t.vertex(1))),
                _crown.Add(moist::exact::Point(triangle._t.vertex(2))),
                _crown.Add(moist::exact::Point(p_opposite_b))
            ), true);
        }
        else
        {
            OOC_WARNING("invalid triangle correspondence!");
        }
    }
}
