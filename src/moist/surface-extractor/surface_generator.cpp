#include "surface_generator.hpp"

#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/utils.hpp"

moist::SurfaceGenerator::SurfaceGenerator(const Tiff &tiff, const uint32_t offset, const Axis axis, const bool center, const bool invert)
{
    _width = tiff.width();
    _height = tiff.height();
    _depth = tiff.depth();
    _offset = offset;
    _axis = axis;
    _invert = invert;
    _center = center;

    // create a "skirt" of 0-values around the entire mesh to close it off
    // _grid.set_grid_dimensions(tiff.width() + 2, tiff.height() + 2, tiff.depth() + 2);

    const size_t sx = tiff.width();
    const size_t sy = tiff.height();
    const size_t sz = tiff.depth();

    _grid.set_grid_dimensions(sx/* + 2*/, sy/* + 2*/, sz/* + 2*/);
    _mc.set_grid3d(_grid);

    for (size_t x = 0; x < sx/* + 2*/; x++)
    {
        for (size_t y = 0; y < sy/* + 2*/; y++)
        {
            for (size_t z = 0; z < sz/* + 2*/; z++)
            {
                _grid.set_grid_value(x, y, z, static_cast<double>(0.0));
            }
        }
    }

    for (size_t x = 0; x < sx; x++)
    {
        for (size_t y = 0; y < sy; y++)
        {
            for (size_t z = 0; z < sz; z++)
            {
                const double datapoint = tiff.data[(z) * tiff.width() * tiff.height() + (y) * tiff.width() + x];
                _grid.set_grid_value(x, y, z, static_cast<double>((invert) ? 1.0 - datapoint : datapoint));
            }
        }
    }

    // for (size_t x = 1; x < sx + 1; x++)
    // {
    //     for (size_t y = 1; y < sy + 1; y++)
    //     {
    //         for (size_t z = 1; z < sz + 1; z++)
    //         {
    //             const double datapoint = tiff.data[(z - 1) * tiff.width() * tiff.height() + (y - 1) * tiff.width() + x - 1];
    //             _grid.set_grid_value(x, y, z, static_cast<double>((invert) ? 1.0 - datapoint : datapoint));
    //         }
    //     }
    // }
}

#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/NL/nl.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/progress.h>
#include <geogram/bibliography/bibliography.h>

static void remesh(geo::Mesh& M_in, geo::Mesh& M_out,
        geo::index_t nb_points,
        geo::coord_index_t dim = 0,
        geo::index_t nb_Lloyd_iter = 1,
        geo::index_t nb_Newton_iter = 1,
        geo::index_t Newton_m = 7,
        bool adjust = true,
        double adjust_max_edge_distance = 0.5)
{
    if(dim == 0) {
        dim = geo::coord_index_t(M_in.vertices.dimension());
    }

    geo::geo_argused(dim);

    geo::Stopwatch W("Remesh");

    geo::CentroidalVoronoiTesselation CVT(&M_in);

    if(nb_points == 0) {
        nb_points = M_in.vertices.nb();
    }
    //CVT.compute_initial_sampling(nb_points, true); // true: for verbose
    CVT.set_points(M_in.vertices.nb(), M_in.vertices.point_ptr(0));
    double lowest_z = 1000.0;
    double highest_z = -1000.0;
    std::size_t locked = 0;
    for (geo::index_t v = 0; v < M_in.vertices.nb(); v++)
    {
        const auto point = M_in.vertices.point(v);
        if (point.z < lowest_z)
        {
            lowest_z = point.z;
        }
        if (point.z > highest_z)
        {
            highest_z = point.z;
        }
        if (point.z == 0.9999 || point.z == 1.0 || point.z == 50.0 || point.z == 50.0001)
        {
            CVT.lock_point(v);
            locked++;
        }
    }
    OOC_DEBUG("lowest point " << lowest_z);
    OOC_DEBUG("highest point " << highest_z);
    OOC_DEBUG("locked " << locked << " points");


    try {
        geo::ProgressTask progress("Lloyd", 100);
        CVT.set_progress_logger(&progress);
        CVT.Lloyd_iterations(nb_Lloyd_iter);
    }
    catch(const geo::TaskCanceled&) {
        // TODO_CANCEL
    }

    if(nb_Newton_iter != 0) {
        try {
            geo::ProgressTask progress("Newton", 100);
            CVT.set_progress_logger(&progress);
            CVT.Newton_iterations(nb_Newton_iter, Newton_m);
        }
        catch(const geo::TaskCanceled&) {
            // TODO_CANCEL
        }
    }

    if(M_in.vertices.dimension() == 6 &&
        geo::CmdLine::get_arg_bool("dbg:save_6d")
        ) {
        geo::Logger::out("Remesh")
            << "Saving source mesh into mesh6.obj6" << std::endl;
        mesh_save(M_in, "mesh6.obj6");
        geo::Logger::out("Remesh")
            << "Saving sampling into points6.txt" << std::endl;
        std::ofstream out("points6.txt");
        out << CVT.delaunay()->nb_vertices() << std::endl;
        for(geo::index_t i = 0; i < CVT.delaunay()->nb_vertices(); i++) {
            for(geo::coord_index_t c = 0; c < 6; c++) {
                out << CVT.delaunay()->vertex_ptr(i)[c] << " ";
            }
            out << std::endl;
        }
    }

    // Delete auxiliary storage used for each threads (it uses a lot of RAM,
    //   we need this RAM to create the surface now...)
    CVT.RVD()->delete_threads();

    CVT.set_use_RVC_centroids(
        geo::CmdLine::get_arg_bool("remesh:RVC_centroids")
    );
    bool multi_nerve = geo::CmdLine::get_arg_bool("remesh:multi_nerve");

    geo::Logger::out("Remesh") << "Computing RVD..." << std::endl;

    CVT.compute_surface(&M_out, multi_nerve);
    if(geo::CmdLine::get_arg_bool("dbg:save_ANN_histo")) {
        geo::Logger::out("ANN")
            << "Saving histogram to ANN_histo.dat" << std::endl;
        std::ofstream out("ANN_histo.dat");
        CVT.delaunay()->save_histogram(out);
    }

    if(adjust) {
        mesh_adjust_surface(M_out, M_in, adjust_max_edge_distance);
    }
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <unordered_set>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/property_map.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index Vertex_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

static void geo_to_cgal(const geo::Mesh& geo_mesh, Mesh& cgal_mesh)
{
    std::vector<Vertex_index> index_map(geo_mesh.vertices.nb());

    for (geo::index_t v = 0; v < geo_mesh.vertices.nb(); v++)
    {
        const double* p = geo_mesh.vertices.point_ptr(v);
        Point pt(p[0], p[1], p[2]);
        index_map[v] = cgal_mesh.add_vertex(pt);
    }

    for (geo::index_t f = 0; f < geo_mesh.facets.nb(); f++)
    {
        geo::index_t nb = geo_mesh.facets.nb_vertices(f);
        std::vector<Vertex_index> face_vertices;
        for (geo::index_t lv = 0; lv < nb; lv++)
        {
            face_vertices.push_back(index_map[geo_mesh.facets.vertex(f, lv)]);
        }
        cgal_mesh.add_face(face_vertices);
    }
}


void cgal_to_geo(const Mesh& cgal_mesh, geo::Mesh& geo_mesh)
{
    using vertex_descriptor = typename CGAL::Surface_mesh<Point>::Vertex_index;
    using face_descriptor = typename CGAL::Surface_mesh<Point>::Face_index;

    geo_mesh.clear();
    geo_mesh.vertices.clear();
    geo_mesh.facets.clear();

    /*std::unordered_map<vertex_descriptor, geo::index_t> vertex_map;

    for (vertex_descriptor vd : cgal_mesh.vertices())
    {
        const Point& p = cgal_mesh.point(vd);
        geo::index_t vid = geo_mesh.vertices.create_vertex(geo::vec3(p.x(), p.y(), p.z()));
        vertex_map[vd] = vid;
    }

    for (face_descriptor fd : cgal_mesh.faces())
    {
        std::vector<geo::index_t> facet;
        for (auto hd : CGAL::halfedges_around_face(cgal_mesh.halfedge(fd), cgal_mesh))
        {
            vertex_descriptor vd = cgal_mesh.target(hd);
            facet.push_back(vertex_map.at(vd));
        }
        geo_mesh.facets.create_triangle(facet[0], facet[1], facet[2]);
    }*/
    //geo_mesh.facets.connect();

    std::map<vertex_descriptor, geo::index_t> vertex_to_index;
    geo::vector<double> vertices(cgal_mesh.number_of_vertices() * 3);
    int vv = 0;
    for (auto v : cgal_mesh.vertices())
    {
        const Point& p = cgal_mesh.point(v);
        vertices[3 * vv] = p.x();
        vertices[3 * vv + 1] = p.y();
        vertices[3 * vv + 2] = p.z();
        vertex_to_index[v] = vv;
        vv++;
    }

    geo::vector<GEO::index_t> triangles;
    geo::index_t t = 0;
    for (face_descriptor fd : cgal_mesh.faces())
    {
        std::vector<geo::index_t> facet;
        for (auto hd : CGAL::halfedges_around_face(cgal_mesh.halfedge(fd), cgal_mesh))
        {
            facet.push_back(vertex_to_index[cgal_mesh.target(hd)]);
        }
        triangles.resize(triangles.size() + 3);
        triangles[3 * t] = facet[0];
        triangles[3 * t + 1] = facet[1];
        triangles[3 * t + 2] = facet[2];
        t++;
    }

    geo_mesh.facets.assign_triangle_mesh((GEO::coord_index_t) 3, vertices, triangles, true);
    moist::utils::geogram::save("mesh_converted.off", geo_mesh);
}

#include <boost/property_map/property_map.hpp>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
static void remesh_CGAL(geo::Mesh& mesh, geo::Mesh& output)
{
    Mesh cgal_mesh;
    geo_to_cgal(mesh, cgal_mesh);

    std::set<Mesh::Vertex_index> constrained_vertices;
    for(auto v : cgal_mesh.vertices())
    {
        const auto point = cgal_mesh.point(v);
        if (point.z() == 0.999 || point.z() == 1.0 || point.z() == 10.001 || point.z() == 10.0)
        {
            constrained_vertices.insert(v);
        }
    }

    OOC_DEBUG("Constraining: " << constrained_vertices.size() << " vertices");
    CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);
    //auto vertex_constrained_map = CGAL::make_function_property_map<vertex_descriptor, bool>(
    //    [&locked_vertices](vertex_descriptor v)
    //    {
    //        return locked_vertices.count(v) > 0;
    //    }
    //);
    //CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);
    //CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
    //CGAL::merge_border_vertices(mesh, 1e-12);
    //CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(cgal_mesh);
    bool valid = CGAL::is_valid_polygon_mesh(cgal_mesh);

    OOC_DEBUG("begin surface smoothing...");
    CGAL::Polygon_mesh_processing::smooth_shape(cgal_mesh, 0.01,
    CGAL::parameters::
        number_of_iterations(10)
        .vertex_is_constrained_map(vcmap)
        //.prevent_surface_shrinkage(true)
    );
    OOC_DEBUG("finalized surface smoothing...");

    cgal_to_geo(cgal_mesh, output);
}

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

static bool is_hole_on_plane(halfedge_descriptor h, Mesh& mesh, double z_plane, double tolerance = 1e-6)
{
    int num_hole_edges = 0;
    for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh))
    {
        const Point& p = mesh.point(target(hc, mesh));
        if (std::abs(p.z() - z_plane) < tolerance)
        {
            return false;
        }
    }

  return true;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

typedef bool FaceInfo;

template <typename K_>
class My_face_base :
  public CGAL::Triangulation_face_base_with_info_2<
      bool, K_,
      CGAL::Constrained_triangulation_face_base_2<K_,
          CGAL::Triangulation_face_base_2<K_>>> {
    typedef CGAL::Triangulation_face_base_with_info_2<
        bool, K_,
        CGAL::Constrained_triangulation_face_base_2<K_,
            CGAL::Triangulation_face_base_2<K_>>> Base;

public:
    My_face_base() : Base() {}

    bool is_in_domain() const { return this->info(); }
    void set_in_domain(bool b) { this->info() = b; }
};

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef My_face_base<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS> CDT;
typedef CDT::Point Point2;

static void fill_holes_on_z(Mesh& mesh, double z)
{
    unsigned int nb_holes = 0;
    std::vector<halfedge_descriptor> border_cycles;

    // collect one halfedge per boundary cycle
    OOC_DEBUG("collection border cycles on z = " << z << "...");
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    OOC_DEBUG("collected " << border_cycles.size() << " border cycles on z = " << z << "...");

    for(halfedge_descriptor h : border_cycles)
    {
        if (is_hole_on_plane(h, mesh, z))
        {
            std::vector<face_descriptor>  patch_facets;
            std::vector<vertex_descriptor> patch_vertices;

            std::vector<Point2> polygon_2d;
            for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh))
            {
                const auto& p = mesh.point(target(hc, mesh));
                polygon_2d.emplace_back(p.x(), p.y());
            }

            CDT cdt;
            for (std::size_t i = 0; i < polygon_2d.size(); i++)
            {
                cdt.insert_constraint(polygon_2d[i], polygon_2d[(i + 1) % polygon_2d.size()]);
            }
            CGAL::mark_domain_in_triangulation(cdt);

            /*for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++)
            {
                OOC_DEBUG(fit->info());
                if (fit->info() == 2)
                {
                    Point2 p0 = fit->vertex(0)->point();
                    Point2 p1 = fit->vertex(1)->point();
                    Point2 p2 = fit->vertex(2)->point();

                    Point P0(p0.x(), p0.y(), z);
                    Point P1(p1.x(), p1.y(), z);
                    Point P2(p2.x(), p2.y(), z);

                    auto v0 = mesh.add_vertex(P0);
                    auto v1 = mesh.add_vertex(P1);
                    auto v2 = mesh.add_vertex(P2);
                    mesh.add_face(v0, v1, v2);
                }
            }*/
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++)
            {
                if (fit->info() == 1) {
                    Point2 p0 = fit->vertex(0)->point();
                    Point2 p1 = fit->vertex(1)->point();
                    Point2 p2 = fit->vertex(2)->point();

                    Point P0(p0.x(), p0.y(), z);
                    Point P1(p1.x(), p1.y(), z);
                    Point P2(p2.x(), p2.y(), z);

                    auto v0 = mesh.add_vertex(P0);
                    auto v1 = mesh.add_vertex(P1);
                    auto v2 = mesh.add_vertex(P2);
                    mesh.add_face(v0, v1, v2);
                }
            }
                        ++nb_holes;
            //bool success = std::get<0>(
            //PMP::triangulate_refine_and_fair_hole(mesh, h,
            //CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
            //.vertex_output_iterator(std::back_inserter(patch_vertices))));
            if (nb_holes % 25 == 0)
            {
                OOC_DEBUG("filled " << nb_holes << " on z = " << z << "...");
            }
        }
    }

}

static void fill_holes(geo::Mesh& mesh, geo::Mesh& output, double z0, double z1)
{
    Mesh cgal_mesh;
    geo_to_cgal(mesh, cgal_mesh);
    fill_holes_on_z(cgal_mesh, static_cast<double>(z0));
    fill_holes_on_z(cgal_mesh, static_cast<double>(z1));
    cgal_to_geo(cgal_mesh, output);
}

void moist::SurfaceGenerator::Generate(GEO::Mesh &mesh, const float isovalue)
{
    if (mesh.vertices.dimension() != 3)
    {
        OOC_ERROR("invalid mesh dimension, expected 3");
    }

    if (isovalue > 1.0f || isovalue < 0.0f)
    {
        OOC_ERROR("invalid isosurface bounds");
    }

    _mc.calculate_isosurface(_surface, static_cast<double>(isovalue));

    GEO::vector<double> vertices(_surface.get_num_vertices() * 3);
    const float c_offset_x = (_center) ? static_cast<float>(_width)  / 2.0f : 0.0f;
    const float c_offset_y = (_center) ? static_cast<float>(_height) / 2.0f : 0.0f;
    const float c_offset_z = (_center) ? static_cast<float>(_depth)  / 2.0f : 0.0f;

    for(g_index v = 0; v < _surface.get_num_vertices(); v++)
    {
        const auto vertex = _surface.getVertex(v);
        // after all this time, why shouldn't I just move everything here? :)
        vertices[3 * v] = vertex[0] - c_offset_x + ((_axis == moist::Axis::X) ? _offset : 0);
        vertices[3 * v + 1] = vertex[1] - c_offset_y + ((_axis == moist::Axis::Y) ? _offset : 0);
        vertices[3 * v + 2] = vertex[2] - c_offset_z + ((_axis == moist::Axis::Z) ? _offset : 0);
    }

    GEO::vector<GEO::index_t> triangles(_surface.get_num_triangles() * 3);
    for(g_index t = 0; t < _surface.get_num_triangles(); t++)
    {
        const auto triangle = _surface.getTriangle(t);
        triangles[3 * t] = triangle[0];
        triangles[3 * t + 1] = triangle[1];
        triangles[3 * t + 2] = triangle[2];
    }

    geo::Mesh before(3);
    before.facets.assign_triangle_mesh((GEO::coord_index_t) 3, vertices, triangles, true);
    moist::utils::geogram::save("before_.off", before);
    fill_holes(before, mesh, _offset, _offset + _depth);

    //remesh_CGAL(before, mesh);
}
