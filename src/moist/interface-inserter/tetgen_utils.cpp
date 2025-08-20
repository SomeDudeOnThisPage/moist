#include "tetgen_utils.hpp"

#include <geogram/mesh/mesh_repair.h>

#include <moist/core/utils.hpp>

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
    // pre-transform to a geogram mesh to utilize its' mesh repair functions
    geo::Mesh geo_mesh(3);
    auto v_interface = geo::Attribute<bool>(geo_mesh.vertices.attributes(), "v_interface");
    for (const auto& p : mesh.Points())
    {
        geo::index_t v = geo_mesh.vertices.create_vertex(geo::vec3(CGAL::to_double(p._p.x()), CGAL::to_double(p._p.y()), CGAL::to_double(p._p.z())));
        v_interface[v] = p._interface;
    }

    geo::vector<geo::index_t> to_delete(mesh.NbCells());
    for (const auto c : mesh.Cells())
    {
        if (c._deleted || moist::geometry::exact::is_degenerate(c))
        {
            continue;
        }

        const auto t = geo_mesh.cells.create_tet(
            geo::index_t(c._points[0]),
            geo::index_t(c._points[1]),
            geo::index_t(c._points[2]),
            geo::index_t(c._points[3])
        );

        if (std::fabs(geo::mesh_cell_volume(geo_mesh, t)) <= 1e-14)
        {
            OOC_WARNING("degenerate tet " << t);
            to_delete[t] = true;
        }
    }
    geo_mesh.cells.delete_elements(to_delete);

    geo::mesh_repair(geo_mesh);

    OOC_DEBUG("number of points: " << geo_mesh.vertices.nb());
    OOC_DEBUG("number of cells: " << geo_mesh.cells.nb());

    tetgen_io.initialize();
    tetgen_io.numberofpoints = geo_mesh.vertices.nb();
    tetgen_io.pointlist = new double[static_cast<std::size_t>(tetgen_io.numberofpoints) * 3]();
    tetgen_io.pointmarkerlist = new int[tetgen_io.numberofpoints];

    for (std::size_t v = 0; v < tetgen_io.numberofpoints; v++)
    {
        const auto& point = geo_mesh.vertices.point(v);
        tetgen_io.pointlist[3 * v] = point.x;
        tetgen_io.pointlist[3 * v + 1] = point.y;
        tetgen_io.pointlist[3 * v + 2] = point.z;
        tetgen_io.pointmarkerlist[v] = (v_interface[v] ? -1 : 1); // -1 -> can be coarsened
    }

    tetgen_io.numberoftetrahedra = geo_mesh.cells.nb();
    tetgen_io.tetrahedronlist = new int[static_cast<std::size_t>(tetgen_io.numberoftetrahedra) * 4]();

    for (std::size_t c = 0; c < tetgen_io.numberoftetrahedra; c++)
    {
        #pragma unroll 4
        for (int lv = 0; lv < 4; lv++)
        {
            tetgen_io.tetrahedronlist[4 * c + lv] = geo_mesh.cells.vertex(c, lv);
        }
    }
}

void moist::tetgen::transform(const tetgenio& tetgen_io, geo::Mesh& mesh)
{
    mesh.clear(false, false);
    for (int v = 0; v < tetgen_io.numberofpoints; v++)
    {
        mesh.vertices.create_vertex(geo::vec3(
            tetgen_io.pointlist[3 * v],
            tetgen_io.pointlist[3 * v + 1],
            tetgen_io.pointlist[3 * v + 2]
        ));
    }

    for (int c = 0; c < tetgen_io.numberoftetrahedra; c++)
    {
        mesh.cells.create_tet(
            tetgen_io.tetrahedronlist[4 * c],
            tetgen_io.tetrahedronlist[4 * c + 1],
            tetgen_io.tetrahedronlist[4 * c + 2],
            tetgen_io.tetrahedronlist[4 * c + 3]
        );
    }
}
