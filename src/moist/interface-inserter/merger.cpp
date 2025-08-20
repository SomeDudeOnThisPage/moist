#include "merger.hpp"

#include <tetgen.h>

#include <geogram/mesh/mesh_repair.h>

#include "moist/core/predicates.inl"
#include "moist/core/utils.hpp"
#include "moist/core/timer.hpp"

#include "arrangement.hpp"
#include "exact_types.hpp"
#include "local_operations.hpp"
#include "new_predicates.inl"
#include "geometry_exact.inl"
#include "tetgen_utils.hpp"

static geo::Mesh dbg_triangles(3);

moist::Merger::Merger(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane)
{
    // this->Prune(_crown, plane.extent);

    this->ConstructExactMesh(a, _mesh_a, plane);
    this->ConstructExactMesh(b, _mesh_b, plane);
    _crown.ResetGrid(10.0);
    moist::utils::geogram::save("crownless_a.mesh", a);
    moist::utils::geogram::save("crownless_b.mesh", b);

    this->InsertAndMapPoints(_mesh_a, _mesh_b);
    this->InsertAndMapPoints(_mesh_b, _mesh_a);

#ifndef NDEBUG
    _mesh_a.DebugMesh("exact_a.mesh");
    _mesh_b.DebugMesh("exact_b.mesh");
#endif // NDEBUG

    this->InsertBToA();

    #ifndef NDEBUG
    _crown.DebugMesh("crown.mesh");
    moist::utils::geogram::save("triangles.off", dbg_triangles);
#endif // NDEBUG

    this->Prune(_crown, plane.extent);
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

        if (!p0._interface || !p1._interface)
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

static std::pair<moist::predicates::PointInTet, double> get_shortest_edge(moist::ExactMesh& mesh, const std::size_t c)
{
    const auto cell = mesh.Cell(c);
    double distance = std::numeric_limits<double>::max();
    moist::predicates::PointInTet shortest_edge = moist::predicates::PointInTet::NONE;
    for (const auto& [i, j] : moist::geometry::exact::TET_EDGE_DESCRIPTOR)
    {
        const auto p0 = mesh.Point(cell[i]);
        const auto p1 = mesh.Point(cell[j]);

        if (!p0._interface || !p1._interface)
        {
            continue;
        }

        const auto _distance = CGAL::to_double(CGAL::squared_distance(p0._p, p1._p));
        if (_distance < distance)
        {
            distance = _distance;
            if (i == 0 && j == 1)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE01;
            }
            else if (i == 0 && j == 2)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE02;
            }
            else if (i == 0 && j == 3)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE03;
            }
            else if (i == 1 && j == 2)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE12;
            }
            else if (i == 1 && j == 3)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE13;
            }
            else if (i == 2 && j == 3)
            {
                shortest_edge = moist::predicates::PointInTet::EDGE23;
            }
        }
    }

    if (shortest_edge != moist::predicates::PointInTet::NONE)
    {
        return std::make_pair(shortest_edge, distance);
    }
    return std::make_pair(moist::predicates::PointInTet::NONE, 0.0);
    // else
    // {
    //     OOC_DEBUG("cell without an interface edge...");
    // }
}

extern "C"
{
  #include "mmg/mmg3d/libmmg3d.h"
}

static void mark_required_triangles(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, const double extent)
{
    constexpr int faces[4][3] = {
        {1,2,3},
        {0,2,3},
        {0,1,3},
        {0,1,2}
    };

    int ne, nv, na, nt, np;
    MMG3D_Get_meshSize(mmgMesh, &np, &ne, &na, &nt, NULL, NULL);

    int num_elements = ne;
    int num_vertices = np;
    int num_triangles = 0;
    std::vector<std::array<int, 3>> triangles;
    std::vector<std::array<int, 4>> cells;
    std::vector<std::array<double, 3>> points;

    for (int v = 1; v <= np; v++)
    {
        points.push_back({mmgMesh->point[v].c[0], mmgMesh->point[v].c[1], mmgMesh->point[v].c[2]});
    }

    // we must count boundary faces ourselves :(
    std::map<std::array<int, 3>, int> face_count;
    for (int c = 1; c <= ne; c++)
    {
        const int *v = mmgMesh->tetra[c].v;
        cells.push_back({mmgMesh->tetra[c].v[0], mmgMesh->tetra[c].v[1], mmgMesh->tetra[c].v[2], mmgMesh->tetra[c].v[3]});
        for (int lf=0; lf < 4; lf++)
        {
            std::array<int, 3> key = { v[faces[lf][0]], v[faces[lf][1]], v[faces[lf][2]] };
            std::sort(key.begin(), key.end());
            face_count[key] += 1;
            if (face_count[key] > 2)
            {
                OOC_WARNING("non manifold face");
            }
        }
    }

    auto is_z_extent = [&](int vidx) -> bool
    {
        return mmgMesh->point[vidx].c[2] == extent;
    };

    int triCount = 0;
    // Iterate over tets
    for (int k = 1; k <= ne; k++)
    {
        const int* v = mmgMesh->tetra[k].v;
        int z_extent_count = 0;
        for (int j = 0; j < 4; j++)
        {
            if (is_z_extent(v[j]))
            {
                z_extent_count++;
            }
        }

        if (z_extent_count == 1)
        {
            // its👏a👏boundary👏tet👏
            int lv_extent = -1;
            for (int lv = 0; lv < 4; lv++)
            {
                if (is_z_extent(v[lv]))
                {
                    lv_extent = lv;
                }
            }

            num_triangles++;
            const auto v0 = v[(lv_extent + 1) % 4];
            const auto v1 = v[(lv_extent + 2) % 4];
            const auto v2 = v[(lv_extent + 3) % 4];
            triangles.push_back({v[(lv_extent + 1) % 4], v[(lv_extent + 2) % 4], v[(lv_extent + 3) % 4]});
        }
    }

    MMG3D_Set_meshSize(mmgMesh, num_vertices, num_elements, 0, num_triangles, 0, 0);
    const auto points_size = points.size();
    const auto cells_size = cells.size();
    for (std::size_t v = 1; v <= num_vertices; v++)
    {
        MMG3D_Set_vertex(mmgMesh, points[v - 1][0], points[v - 1][1], points[v - 1][2], 0, v);
    }
    for (std::size_t f = 1; f <= num_triangles; f++)
    {
        MMG3D_Set_triangle(mmgMesh, triangles[f - 1][0], triangles[f - 1][1], triangles[f - 1][2], 0, f);
        MMG3D_Set_requiredTriangle(mmgMesh, f);
    }
    for (std::size_t c = 1; c <= num_elements; c++)
    {
        if (c == 494 || c == 495 || c == 496)
        {
            const std::array<int, 4> cc = {{
                cells[c - 1][0], cells[c - 1][1], cells[c - 1][2], cells[c - 1][3]
            }};

            std::vector<std::array<double, 3>> pp;
            for (int i = 0; i < 4; i++)
            {
                const auto cp = points[cells[c - 1][i]];
                pp.push_back({cp.at(0), cp.at(1), cp.at(2)});
            }
            int ggg = 0;
        }
        MMG3D_Set_tetrahedron(mmgMesh, cells[c - 1][0], cells[c - 1][1], cells[c - 1][2], cells[c - 1][3], 0, c);
    }

    double min_quality = 100000000.0;
    std::size_t worst_tet = 0;
    for (std::size_t c = 1; c <= num_elements; c++)
    {
        if (c == 494 || c == 495 || c == 496)
        {
            double q = MMG3D_Get_tetrahedronQuality(mmgMesh, mmgSol, c);
            int ggwgw = 1;
        }
        double q = MMG3D_Get_tetrahedronQuality(mmgMesh, mmgSol, c);
        if (q < min_quality)
        {
            min_quality = q;
            worst_tet = c;
        }
        if (q < 1.e-30)
        {
            OOC_WARNING("NULKAL");
        }
        if (q < 1.e-15)
        {
            OOC_WARNING("EPSOK " << c);
        }
    }

    OOC_DEBUG("min_quality: " << min_quality);
}

void moist::Merger::Prune(moist::ExactMesh& mesh, const double extent)
{
    MOIST_INFO("remeshing...");

    tetgenio coarsen_mesh;
    tetgenio output_coarsen_mesh;
    output_coarsen_mesh.initialize();
    moist::tetgen::transform(_crown, coarsen_mesh);
    tetgenbehavior coarsen_pass;
    coarsen_pass.nobisect = true;
    coarsen_pass.supsteiner_level = 0;
    coarsen_pass.refine = true;
    coarsen_pass.coarsen = true;
    tetrahedralize(&coarsen_pass, &coarsen_mesh, &output_coarsen_mesh);
    geo::Mesh geo_coarsen_mesh(3);
    moist::tetgen::transform(coarsen_mesh, geo_coarsen_mesh);
    moist::utils::geogram::save("dbg_after_coarsening.mesh", geo_coarsen_mesh);

    MMG5_pMesh      mmgMesh;
    MMG5_pSol       mmgSol;

    mmgMesh = NULL;
    mmgSol  = NULL;
    MMG3D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh,&mmgMesh,
                    MMG5_ARG_ppMet,&mmgSol,
                    MMG5_ARG_end);

    MMG3D_loadMesh(mmgMesh, "crown.mesh");
    mark_required_triangles(mmgMesh, mmgSol, extent);

    //if (MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1)
    //{
    //    OOC_WARNING("bad mesh data");
    //}

    //MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nosplit, 1);
    //MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_angle, 1);
    // MMG3D_SET_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_)
    int nv;
    MMG3D_Get_meshSize(mmgMesh, &nv, NULL, NULL, NULL, NULL, NULL);
    MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, nv, MMG5_Scalar);
    //MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_opnbdy, 1);
    //MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nofem, 1);
    //MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_angle, 1);
    for (int i = 1; i <= nv; i++)
    {
        MMG3D_Set_scalarSol(mmgSol, 1.0, i);
    }

    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, 0.1);
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, 0.05);

    if (MMG3D_saveMesh(mmgMesh, "output_mmg3d_before_remesning.mesh") != 1)
    {
        std::cerr << "Error: cannot save output.mesh\n";
    }

    //if (MMG3D_mmg3dls(mmgMesh, mmgSol, NULL) == MMG5_STRONGFAILURE )
    //{
    //    std::cerr << "Error: MMG3D remeshing failed\n";
    //}
    if (MMG3D_mmg3dlib(mmgMesh, mmgSol) == MMG5_STRONGFAILURE )
    {
        std::cerr << "Error: MMG3D remeshing failed\n";
    }

    if (MMG3D_saveMesh(mmgMesh, "output_mmg3d.mesh") != 1)
    {
        std::cerr << "Error: cannot save output.mesh\n";
    }

    MOIST_INFO("finished remeshing...");
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
    double lowest = 10000.0;
    double highest = -100000.0;

#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::Merger::ConstructExactMesh");
#endif // NDEBUG

    std::unordered_map<geo::index_t, std::size_t> added_vertices;
    geo::vector<geo::index_t> to_delete(mesh.cells.nb());
    for (const geo::index_t c : mesh.cells)
    {
        if (!moist::predicates::cell_on_plane(c, mesh, plane))
        {
            continue;
        }

        std::array<std::size_t, 4> vertices;
        for (geo::index_t lv = 0; lv < 4; lv++)
        {
            const geo::index_t v = mesh.cells.vertex(c, lv);
            if (added_vertices.contains(v))
            {
                vertices[lv] = added_vertices[v];
                continue;
            }
            else
            {
                // here we must keep _other of the point always on NO_VERTEX...
                auto& point = mesh.vertices.point(v);
                if (point.z < lowest)
                {
                    lowest = point.z;
                }
                if (point.z > highest)
                {
                    highest = point.z;
                }
                bool interface = moist::predicates::point_on_plane(point, plane);
                if (interface)
                {
                    point[2] = plane.extent;
                }

                vec3 point_rounded = vec3(round16(point[0]), round16(point[1]), round16(point[2]));
                vertices[lv] = exact_mesh.Add(moist::exact::Point(point_rounded, interface));
                added_vertices[v] = vertices[lv];
            }
        }

        if (geo::mesh_cell_volume(mesh, c) != 0.0)
        {
            exact_mesh.Add(moist::exact::Cell(
                vertices[0], vertices[1], vertices[2], vertices[3],
                moist::geometry::exact::get_cell_type(vertices, exact_mesh)
            ), true);
        }
        else
        {
            OOC_DEBUG("skipped cell " << c << " during interface plane mesh construction due to degeneracy in eps...");
        }

        to_delete[c] = true;
    }
    mesh.cells.delete_elements(to_delete);

    // edge length does not make the most sense, since it must be based relative to the overall length of our bbox diag.
    // dumb hard coded value... use nb of cells on the interface maybe... times some factor controlled so we can evaluate performance based on it...
    constexpr double FACTOR = 0.5;
    exact_mesh.ResetGrid(FACTOR * std::sqrt(exact_mesh.NbCells()));
    OOC_DEBUG("grid resolution: " << FACTOR * std::sqrt(exact_mesh.NbCells()));
    OOC_DEBUG("lowest point: " << lowest);
    OOC_DEBUG("highest point: " << highest);
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
bool moist::Merger::InsertPointIntoEdges(const moist::exact::Point& point, moist::ExactMesh& mesh)
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::Merger::InsertPointIntoEdges");
#endif // NDEBUG
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
                if (mesh.Point(cell._points[lv])._interface)
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
                    mesh.Point(v_1to2)._interface = true;
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, mesh);
                found = true;
                inserted.insert(c);
            }
            else if (location == moist::predicates::PointInTet::FACET)
            {
                // const std::size_t v = mesh.Add(point);
                // moist::operation::exact::InsertVertexOnCellBoundaryFacet(c, v, mesh);
            }
        }
    }

    return found;
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
                to.Point(v_new)._interface = true;
            }
            else if (point_in_tet_on_edge(location))
            {
                if (v_1to2 == moist::exact::NO_VERTEX)
                {
                    v_1to2 = to.Add(moist::exact::Point(point));
                    to.Point(v_1to2)._other = v;
                    to.Point(v_1to2)._interface = true;
                }

                moist::operation::exact::InsertVertexOnCellBoundaryEdge(c, v_1to2, location, to);
            }
        }
    }
}

static std::unordered_set<std::size_t> touched_a;
static geo::Mesh dbg_full(3);
void moist::Merger::InsertBToA()
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::Merger::InsertBToA");
#endif // NDEBUG

    for (std::size_t c = 0; c < _mesh_b.NbCells(); c++)
    {
        if (c % 50 == 0)
        {
            OOC_DEBUG("merged " << c << " out of " << _mesh_b.NbCells() << " cells...");
        }

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
        _mesh_b.Cell(c)._deleted = true;
    }

    for (std::size_t c = 0; c < _mesh_a.NbCells(); c++)
    {
        const auto& cell = _mesh_a.Cell(c);
        if (cell._deleted || touched_a.contains(c))
        {
            continue;
        }

        const auto p0 = _mesh_a.Point(cell._points[0]);
        const auto p1 = _mesh_a.Point(cell._points[1]);
        const auto p2 = _mesh_a.Point(cell._points[2]);
        const auto p3 = _mesh_a.Point(cell._points[3]);

        _crown.Add(moist::exact::Cell(
            _crown.Add(moist::exact::Point(p0)),
            _crown.Add(moist::exact::Point(p1)),
            _crown.Add(moist::exact::Point(p2)),
            _crown.Add(moist::exact::Point(p3)),
            moist::exact::CellType::I
        ), true);
    }

    for (std::size_t c = 0; c < _mesh_b.NbCells(); c++)
    {
        const auto& cell = _mesh_b.Cell(c);
        if (cell._deleted)
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
            _crown.Add(moist::exact::Point(p3)),
            moist::exact::CellType::I
        ), true);
    }
}

bool moist::Merger::InsertCellBToA(const std::size_t cb)
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::Merger::InsertCellBToA");
#endif // NDEBUG

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
            const auto cell = _mesh_a.Cell(ca);
            if (cell._deleted || inserted.contains(ca))
            {
                continue;
            }

            if (!moist::geometry::exact::cell_has_interface_facet(cell, _mesh_a))
            {
                continue;
            }

            if (this->InsertCellIntoCell(ca, cb))
            {
                inserted.emplace(ca);
            }
        }
    }
    return !inserted.empty();
}

static moist::exact::Triangle get_interface_triangle(const moist::exact::Cell& cell, const moist::ExactMesh& mesh)
{
    std::vector<moist::exact::Point> points;
    for (const auto v : cell._points)
    {
        const auto point = mesh.Point(v);
        if (point._interface)
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
        if (!point._interface)
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

// remove all cells of crown that intersect the new cell...
static void remove_intersecting_III_cells(const std::size_t c, moist::ExactMesh& mesh)
{
    const auto cell = mesh.Cell(c);
    const auto t0 = get_interface_triangle(cell, mesh);

    const auto grid_cells = mesh.Grid()->GetCells(moist::create_cell_bbox2d_exact(c, mesh));
    for (const auto grid_cell : grid_cells)
    {
        if (!mesh.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
        {
            continue;
        }

        for (auto lc : mesh.Grid()->GetMeshCells(grid_cell))
        {
            auto& cell_other = mesh.Cell(lc);
            if (lc == c || cell_other._deleted || cell_other._type != moist::exact::CellType::III)
            {
                continue;
            }

            const auto t1 = get_interface_triangle(mesh.Cell(lc), mesh);
            if (moist::exact::arrangeable(t0, t1))
            {
                mesh.Cell(lc)._deleted = true;
            }
        }
    }
}

static bool marked_for_debug = false;
bool moist::Merger::InsertCellIntoCell(const std::size_t ca, const std::size_t cb)
{
#ifndef NDEBUG
    moist::ScopeTimer TIMER("moist::Merger::InsertCellIntoCell");
#endif // NDEBUG

    auto& cell_a = _mesh_a.Cell(ca);
    auto& cell_b = _mesh_b.Cell(cb);

    const auto t0 = get_interface_triangle(cell_a, _mesh_a);
    const auto t1 = get_interface_triangle(cell_b, _mesh_b);
    const std::size_t v_opposite_a = get_opposite_vertex(cell_a, _mesh_a);
    const std::size_t v_opposite_b = get_opposite_vertex(cell_b, _mesh_b);

    if (!moist::exact::arrangeable(t0, t1))
    {
        return false;
    }

    const auto triangulation = moist::exact::arrange(t0, t1);

    if (triangulation.empty())
    {
        return false;
    }

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
                bool in_a = true;
                bool in_b = true;
                if (!triangle_has_point(t0, triangle._t.vertex(lv)))
                {
                    in_a = this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_a);
                }

                if (!triangle_has_point(t1, triangle._t.vertex(lv)))
                {
                    in_b = this->InsertPointIntoEdges(moist::exact::Point(triangle._t.vertex(lv)), _mesh_b);
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

            auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
            auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
            auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
            points_triangle_0._interface = true;
            points_triangle_1._interface = true;
            points_triangle_2._interface = true;
            std::size_t v0 = _crown.Add(points_triangle_0);
            std::size_t v1 = _crown.Add(points_triangle_1);
            std::size_t v2 = _crown.Add(points_triangle_2);

            _crown.Add(moist::exact::Cell(
                v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_a)),
                moist::exact::CellType::III
            ), true);

            _crown.Add(moist::exact::Cell(
                v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_b)),
                moist::exact::CellType::III
            ), true);
        }
        else if (correspondence == moist::exact::VertexCorrespondence::A)
        {
            // the cell is only present in mesh b in this arragement -> check if it intersects with any cell of mesh a, if not, add it...
            const auto aabb = moist::create_triangle_bbox2d_exact(triangle);
            const auto grid_cells = _mesh_b.Grid()->GetCells(aabb);
            bool do_add = true;
            for (const auto grid_cell : grid_cells)
            {
                if (!_mesh_b.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
                {
                    continue;
                }

                for (const auto cb_check : _mesh_b.Grid()->GetMeshCells(grid_cell))
                {
                    if (_mesh_b.Cell(cb_check)._type != moist::exact::CellType::III)
                    {
                        continue;
                    }

                    const auto t0 = get_interface_triangle(_mesh_b.Cell(cb_check), _mesh_b);
                    if (moist::exact::arrangeable(t0, triangle))
                    {
                        // if we have at least one triangular overlap, we must not add it to the crown.
                        do_add = false;
                        goto add_cell_b_into_crown;
                    }
                }
            }
add_cell_b_into_crown:
            if (do_add)
            {
                auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
                auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
                auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
                points_triangle_0._interface = true;
                points_triangle_1._interface = true;
                points_triangle_2._interface = true;
                std::size_t v0 = _crown.Add(points_triangle_0);
                std::size_t v1 = _crown.Add(points_triangle_1);
                std::size_t v2 = _crown.Add(points_triangle_2);
                std::size_t new_cell_crown = _crown.Add(moist::exact::Cell(
                    v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_a)),
                    moist::exact::CellType::III
                ), true);
                // we also must delete any intersecting cell in a...
                //remove_intersecting_III_cells(new_cell_crown, _crown);
            }
        }
else if (correspondence == moist::exact::VertexCorrespondence::B)
        {
            // the cell is only present in mesh b in this arragement -> check if it intersects with any cell of mesh a, if not, add it...
            const auto aabb = moist::create_triangle_bbox2d_exact(triangle);
            const auto grid_cells = _mesh_a.Grid()->GetCells(aabb);
            bool do_add = true;
            for (const auto grid_cell : grid_cells)
            {
                if (!_mesh_a.Grid()->_grid.contains(grid_cell)) // cringe bug from old code...
                {
                    continue;
                }

                for (const auto cb_check : _mesh_a.Grid()->GetMeshCells(grid_cell))
                {
                    if (_mesh_a.Cell(cb_check)._type != moist::exact::CellType::III)
                    {
                        continue;
                    }

                    const auto t0 = get_interface_triangle(_mesh_a.Cell(cb_check), _mesh_a);
                    if (moist::exact::arrangeable(t0, triangle))
                    {
                        // if we have at least one triangular overlap, we must not add it to the crown.
                        do_add = false;
                        goto add_cell_a_into_crown;
                    }
                }
            }
add_cell_a_into_crown:
            if (do_add)
            {
                auto points_triangle_0 = moist::exact::Point(triangle._t.vertex(0));
                auto points_triangle_1 = moist::exact::Point(triangle._t.vertex(1));
                auto points_triangle_2 = moist::exact::Point(triangle._t.vertex(2));
                points_triangle_0._interface = true;
                points_triangle_1._interface = true;
                points_triangle_2._interface = true;
                std::size_t v0 = _crown.Add(points_triangle_0);
                std::size_t v1 = _crown.Add(points_triangle_1);
                std::size_t v2 = _crown.Add(points_triangle_2);
                std::size_t new_cell_crown = _crown.Add(moist::exact::Cell(
                    v0, v1, v2, _crown.Add(moist::exact::Point(p_opposite_b)),
                    moist::exact::CellType::III
                ), true);
                // we also must delete any intersecting cell in a...
                //remove_intersecting_III_cells(new_cell_crown, _crown);
            }
        }
    }
    touched_a.emplace(ca);
    return true;
}
