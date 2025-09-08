#include "remeshing.hpp"

#include <geogram/mesh/mesh_repair.h>

#include <moist/core/utils.hpp>
#include <moist/core/geometry.inl>

#include "geometry_exact.inl"

static bool is_degenerate_in_double(const std::array<moist::exact::Point, 4>& points)
{
    vec3 discretized[4];
    for (std::size_t i = 0; i < 4; i++)
    {
        discretized[i] = vec3(points[i].x(), points[i].y(), points[i].z());
    }

    return std::fabs(geo::Geom::tetra_signed_volume(discretized[0], discretized[1], discretized[2], discretized[3])) == 0.0;
}

void moist::tetgen::transform(moist::ExactMesh& mesh, tetgenio& tetgen_io)
{
    tetgen_io.initialize();
    tetgen_io.numberofpoints = mesh.NbPoints();
    tetgen_io.pointlist = new double[static_cast<std::size_t>(tetgen_io.numberofpoints) * 3]();
    tetgen_io.pointmarkerlist = new int[tetgen_io.numberofpoints];

    std::size_t v = 0;
    for (const auto point : mesh.Points())
    {
        if (point._deleted)
        {
            continue;
        }

        tetgen_io.pointlist[3 * v] = point.x();
        tetgen_io.pointlist[3 * v + 1] = point.y();
        tetgen_io.pointlist[3 * v + 2] = point.z();
        tetgen_io.pointmarkerlist[v] = (point._interface ? -1 : 1); // -1 -> can be coarsened
        v++;
    }

    tetgen_io.numberoftetrahedra = mesh.NbCells();
    tetgen_io.tetrahedronlist = new int[static_cast<std::size_t>(tetgen_io.numberoftetrahedra) * 4]();

    std::size_t c = 0;
    for (const auto cell : mesh.Cells())
    {
        if (cell._deleted)
        {
            continue;
        }

        #pragma unroll 4
        for (int lv = 0; lv < 4; lv++)
        {
            tetgen_io.tetrahedronlist[4 * c + lv] = cell[lv];
            if (cell._type == moist::exact::CellType::I)
            {
                const auto point = mesh.Point(cell[lv]);
                if (!point._interface)
                {
                    tetgen_io.pointmarkerlist[cell[lv]] = 22;
                }
            }
        }
        c++;
    }
}

void moist::tetgen::transform(tetgenio& tetgen_io, moist::ExactMesh& to, const moist::AxisAlignedPlane& plane)
{
    to.ResetMesh();
    for (int v = 0; v < tetgen_io.numberofpoints; v++)
    {
        auto point = moist::exact::Point(tetgen_io.pointlist[3 * v], tetgen_io.pointlist[3 * v + 1], tetgen_io.pointlist[3 * v + 2]);
        // we use the interface bool now in the way that we set it to false only for points with marker, everything that is not a constraint for remeshing is set to interface!
        point._interface = tetgen_io.pointmarkerlist[v] == 22; // 2 are marked points that are from the "old" type I cells, since coarsening removes our implicit check for cell types!
        to.Add(point, true);
    }

    for (int c = 0; c < tetgen_io.numberoftetrahedra; c++)
    {
        constexpr std::array<moist::exact::CellType, 5> CELL_TYPE_TABLE = {
            /*0 interf is now I too */ moist::exact::CellType::I, moist::exact::CellType::I, moist::exact::CellType::II, moist::exact::CellType::III, moist::exact::CellType::III
        };

        std::size_t nb_marked_two = 0;
        for (std::size_t lv = 0; lv < 4; lv++)
        {
            const auto point = to.Point(tetgen_io.tetrahedronlist[4 * c + lv]);
            nb_marked_two += (point._interface) ? 1 : 0;
        }

        to.Add(moist::exact::Cell(
            tetgen_io.tetrahedronlist[4 * c],
            tetgen_io.tetrahedronlist[4 * c + 1],
            tetgen_io.tetrahedronlist[4 * c + 2],
            tetgen_io.tetrahedronlist[4 * c + 3],
            nb_marked_two == 3 ? moist::exact::CellType::I : moist::exact::CellType::III
        ), true);
    }
}

void moist::tetgen::coarsen(tetgenio& from, tetgenio& to)
{
    tetgenbehavior coarsen_pass;
    coarsen_pass.nobisect = true;
    coarsen_pass.supsteiner_level = 1;
    coarsen_pass.refine = true;
    coarsen_pass.coarsen = true;

    tetrahedralize(&coarsen_pass, &from, &to);
}

void moist::mmg3d::transform(moist::ExactMesh& from, MMG5_pMesh to, MMG5_pSol solution)
{
    std::vector<std::array<int, 3>> triangles; // constraints
    std::vector<std::array<int, 4>> cells;
    std::vector<std::array<double, 3>> points;

    std::size_t nb_points = 0;
    std::size_t nb_cells = 0;
    std::size_t nb_triangles = 0;

    for (const auto point : from.Points())
    {
        if (!point._deleted)
        {
            nb_points++;
            points.push_back({ point.x(), point.y(), point.z() });
        }
    }

    std::map<std::array<int, 3>, int> face_count;
    for (const auto cell : from.Cells())
    {
        if (cell._deleted)
        {
            continue;
        }

        nb_cells++;
        cells.push_back({ static_cast<int>(cell[0] + 1), static_cast<int>(cell[1] + 1), static_cast<int>(cell[2] + 1), static_cast<int>(cell[3] + 1) });

        if (cell._type != moist::exact::CellType::I)
        {
            continue;
        }

        for (std::size_t lf = 0; lf < 4; lf++)
        {
            std::array<int, 3> facet =
            {
                cell[moist::geometry::TET_FACET_DESCRIPTOR[lf][0]],
                cell[moist::geometry::TET_FACET_DESCRIPTOR[lf][1]],
                cell[moist::geometry::TET_FACET_DESCRIPTOR[lf][2]]
            };

            if (!from.Point(facet[0])._interface && !from.Point(facet[1])._interface && !from.Point(facet[2])._interface)
            {
                // its👏a👏boundary👏facet👏
                nb_triangles++;
                triangles.push_back({facet[0] + 1, facet[1] + 1, facet[2] + 1});
                break;
            }
        }
    }

    MMG3D_Set_meshSize(to, nb_points, nb_cells, 0, nb_triangles, 0, 0);

    // mmg data structures are 1 indexed!!!
    for (std::size_t v = 1; v <= nb_points; v++)
    {
        MMG3D_Set_vertex(to, points[v - 1][0], points[v - 1][1], points[v - 1][2], 0, v);
    }

    for (std::size_t f = 1; f <= nb_triangles; f++)
    {
        MMG3D_Set_triangle(to, triangles[f - 1][0], triangles[f - 1][1], triangles[f - 1][2], 0, f);
        MMG3D_Set_requiredTriangle(to, f);
    }

    for (std::size_t c = 1; c <= nb_cells; c++)
    {
        MMG3D_Set_tetrahedron(to, cells[c - 1][0], cells[c - 1][1], cells[c - 1][2], cells[c - 1][3], 0, c);
    }

    // keep initial size of solution
    MMG3D_Set_solSize(to, solution, MMG5_Vertex, nb_points, MMG5_Scalar);
    for (int i = 1; i <= nb_points; i++)
    {
        MMG3D_Set_scalarSol(solution, 1.0, i);
    }
}

void moist::mmg3d::transform(MMG5_pMesh from, MMG5_pSol solution, moist::ExactMesh& to)
{
    to.ResetMesh();

    int nb_points;
    int nb_cells;
    MMG3D_Get_meshSize(from, &nb_points, &nb_cells, NULL, NULL, NULL, NULL);

    for (int v = 1; v <= nb_points; v++)
    {
        to.Add(moist::exact::Point(from->point[v].c[0], from->point[v].c[1], from->point[v].c[2]), true);
    }

    for (int c = 1; c <= nb_cells; c++)
    {
        const int *v = from->tetra[c].v;
        to.Add(moist::exact::Cell(
            from->tetra[c].v[0] - 1,
            from->tetra[c].v[1] - 1,
            from->tetra[c].v[2] - 1,
            from->tetra[c].v[3] - 1,
            moist::exact::CellType::I), // dont care about type here anymore, its only used fro metrics and saving to fs after this...
        true);
    }
}

void moist::mmg3d::set_solution(MMG5_pMesh mesh, MMG5_pSol solution, double hmin, double hmax, double hausd)
{
    MMG3D_Set_dparameter(mesh, solution, MMG3D_DPARAM_hmin, hmin);
    MMG3D_Set_dparameter(mesh, solution, MMG3D_DPARAM_hmax, hmax);
    MMG3D_Set_dparameter(mesh, solution, MMG3D_DPARAM_hausd, hausd);
}

void moist::mmg3d::remesh(MMG5_pMesh mesh, MMG5_pSol solution)
{
    if (MMG3D_mmg3dlib(mesh, solution) == MMG5_STRONGFAILURE)
    {
        MOIST_INFO("Error: MMG3D remeshing failed");
    }
}
