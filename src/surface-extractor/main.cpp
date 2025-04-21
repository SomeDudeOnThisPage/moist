#include <CLI/CLI.hpp>

#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/ImageIO.h>
#include <CGAL/Image_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

#include "core.hpp"
#include "tiff_data.hpp"

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_Image_3;
typedef CGAL::Implicit_surface_3<GT, Gray_level_Image_3> Surface_3;

/**
 * Parameters:
 *  TODO
 *
 * Output:
 *  TODO
 *
 * Description:
 *  TODO
 */
int main(int argc, char* argv[])
{
    struct Arguments
    {
        std::filesystem::path directory;
        std::string pattern;
        uint32_t start;
        uint32_t num_files;
        double scale_factor;
        bool generate_sizing_field;
    };

    Arguments arguments{};
    CLI::App app{argv[0]};

    app.add_option("-d, --directory", arguments.directory, "Path of .tiff files")
        //->required()
        ->check(CLI::ExistingDirectory);
    app.add_option("-p, --pattern", arguments.pattern, "Pattern for imported files, with '*'-characters denoting numerals")
        //->required()
        ;
    app.add_option("-s, --start", arguments.start, "Start file (number resolved in pattern)")
        //->required()
        ;
    app.add_option("-n, --num-files", arguments.num_files, "Total number of files to mesh")
        //->required()
        ;
    //app.add_option("-s, --scale-factor", arguments.scale_factor, "Scale-Factor for the whole mesh (defaults to 1)")
        //->default_val((uint32_t) 1)
    //    ;
    app.add_option("-f, --generate-sizing-field", arguments.generate_sizing_field, "Generates a sizing field - TODO: More configs for this")
        //->default_val(false)
        ;

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));


    Tr tr;
    C2t3 c2t3(tr);

    incremental_meshing::TiffData tiff_data = incremental_meshing::TiffData(std::filesystem::path("0.tiff"), 1);
    _image pi{};
    pi.xdim = 10;
    pi.ydim = 10;
    pi.zdim = 1;
    pi.vdim = 1;
    pi.vx = 1.0;
    pi.vy = 1.0;
    pi.vz = 1.0;
    pi.tx = 0;
    pi.ty = 0;
    pi.tz = 0;
    pi.rx = 0;
    pi.ry = 0;
    pi.rz = 0;
    pi.cx = 5;
    pi.cy = 5;
    pi.cz = 0;
    pi.data = &tiff_data._data;
    pi.wdim = 2;
    pi.wordKind = WORD_KIND::WK_FIXED;
    pi.sign = SIGN::SGN_UNSIGNED;

    Gray_level_Image_3 image(CGAL::Image_3(&pi), 1.f);

    GT::Point_3 bounding_sphere_center(5., 5., 0.5);
    GT::FT bounding_sphere_squared_radius = 10.*10.*2.;
    GT::Sphere_3 bounding_sphere(bounding_sphere_center, bounding_sphere_squared_radius);

    Surface_3 surface(image, bounding_sphere, 1e-5);

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 5., 5.);

    // meshing surface, with the "manifold without boundary" algorithm
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    std::ofstream out("out.off");
    CGAL::output_surface_facets_to_off(out, c2t3);

    return 0;
}
