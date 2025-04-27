#include <CLI/CLI.hpp>

// TODO: For now, since I had to fix some small things to make it compile with cpp 20, include it directly as files here...
// Before finishing, fork MC33 and fix it / compile a static library with cmake, then fetch content in extern folder
#include "mc33/include/MC33.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <geogram/mesh/mesh_geometry.h>

#include "core.hpp"
#include "utils.hpp"
#include "tiff_data.hpp"

static CLI::Validator PatternPlaceholderCountValidator(const std::size_t placeholders)
{
    return CLI::Validator([placeholders](const std::string &input) -> std::string {
        const size_t amount = incremental_meshing::utils::count_placeholders(input);
        return (amount == placeholders)
            ? ""
            : std::format("invalid pattern {}: expected {} placeholders, got {}", input, placeholders, amount);
    }, "PATTERN");
}

int main(int argc, char* argv[])
{
    struct Arguments
    {
        std::string input;
        std::string output;

        uint32_t first;
        uint32_t amount;
        uint32_t sample_bits;

        incremental_meshing::Axis axis;

        int32_t dir_offset;
        float isovalue;
        double scale_factor;
        bool f_center;
        bool f_generate_sizing_field;
        bool f_invert;
    };

    Arguments arguments{};
    CLI::App app{argv[0]};

    app.add_option("-i, --input", arguments.input, "CPP pattern describing input files. Arguments: [number:uint32_t] (e.g. data/{:03}.tif).")
        ->required()
        ->check(PatternPlaceholderCountValidator(1));
    app.add_option("-o, --output", arguments.output, "CPP pattern describing the output file. Arguments: [from:uint32|to:uint32t] (e.g. slice{}-{}.obj).")
        ->required()
        ->check(PatternPlaceholderCountValidator(2));
    app.add_option("-f, --first", arguments.first, "Initial file number.")
        ->check(CLI::PositiveNumber)
        ->required();
    app.add_option("-n, --amount", arguments.amount, "Total number of files to mesh into this slice.")
        ->check(CLI::PositiveNumber)
        ->required();
    app.add_option("-v, --isovalue", arguments.isovalue, "Isovalue.")
        ->check(CLI::PositiveNumber)
        ->required();
    app.add_option("-b, --sample-bits", arguments.sample_bits, "Bits per sample [8, 16, 32, 64].")
        ->check(CLI::IsMember({8, 16, 32, 64}))
        ->default_val(16);
    app.add_option("-x, --axis", arguments.axis, "Interface axis (X|Y|Z)")
        ->transform(CLI::CheckedTransformer(incremental_meshing::AXIS_OPTION_ARGUMENT_MAP, CLI::ignore_case))
        ->default_val(incremental_meshing::Axis::Z); // TODO: Maybe a "growth-direction" (positive/negative)?
    app.add_option("--directional-offset", arguments.dir_offset, "Directional offset along the growth axis (Z), e.g. if your files start with File500.tif, set this to -500. If you have 100 files and want them centered around 0 in the growth direction, set this to -50.")
        ->default_val(0);
    app.add_flag("--center", arguments.f_center, "Center around 0,0,0")
        ->default_val(false);
    app.add_flag("--generate-sizing-field", arguments.f_generate_sizing_field, "Generates a sizing field - TODO: More configs for this")
        ->default_val(false);
    app.add_flag("--invert", arguments.f_invert, "Invert output (mesh negative space).")
        ->default_val(true);

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    incremental_meshing::TiffData tiff_data(arguments.input, arguments.first, arguments.amount);

    MC33 mc;
    surface surface;
    grid3d grid;
    // create a "skirt" of 0-values around the entire mesh to close it off
    grid.set_grid_dimensions(tiff_data.width() + 2, tiff_data.height() + 2, tiff_data.depth() + 2);
    mc.set_grid3d(grid);

    // since we must transform the data anyway from the scanline-datatype, set it directly in this loop instead of creating another array.
    const float samples_max = static_cast<float>(2 ^ static_cast<uint64_t>(arguments.sample_bits));
    for (uint32_t x = 1; x < tiff_data.width() + 1; x++)
    {
        for (uint32_t y = 1; y < tiff_data.height() + 1; y++)
        {
            for (uint32_t z = 1; z < tiff_data.depth() + 1; z++)
            {
                // TODO: ugly
                const auto datapoint = (((float) tiff_data.data[(z - 1) * tiff_data.width() * tiff_data.height() + (y - 1) * tiff_data.width() + x - 1]) / samples_max);
                grid.set_grid_value(x, y, z, (arguments.f_invert) ? 1.0 - datapoint : datapoint);
            }
        }
    }

    mc.calculate_isosurface(surface, arguments.isovalue);

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::CmdLine::import_arg_group("sys");
    geogram::Logger::instance()->set_quiet(true);

    geogram::vector<double> vertices(surface.get_num_vertices() * 3);
    const float c_offset_x = (arguments.f_center) ? static_cast<float>(tiff_data.width()) / 2.0f : 0.0f;
    const float c_offset_y = (arguments.f_center) ? static_cast<float>(tiff_data.height()) / 2.0f : 0.0f;
    const float c_offset_z = (arguments.f_center) ? static_cast<float>(tiff_data.depth()) / 2.0f : 0.0f;

    for(g_index v = 0; v < surface.get_num_vertices(); v++)
    {
        const auto vertex = surface.getVertex(v);
        // after all this time, why shouldn't I just move everything here? :)
        vertices[3 * v] = vertex[0] - c_offset_x + ((arguments.axis == incremental_meshing::Axis::X) ? arguments.first + arguments.dir_offset : 0);
        vertices[3 * v + 1] = vertex[1] - c_offset_y + ((arguments.axis == incremental_meshing::Axis::Y) ? arguments.first + arguments.dir_offset : 0);
        vertices[3 * v + 2] = vertex[2] - c_offset_z + ((arguments.axis == incremental_meshing::Axis::Z) ? arguments.first + arguments.dir_offset : 0);
    }

    geogram::vector<geogram::index_t> triangles(surface.get_num_triangles() * 3);
    for(g_index t = 0; t < surface.get_num_triangles(); t++)
    {
        const auto triangle = surface.getTriangle(t);
        triangles[3 * t] = triangle[0];
        triangles[3 * t + 1] = triangle[1];
        triangles[3 * t + 2] = triangle[2];
    }

    geogram::Mesh mesh(3);
    mesh.facets.assign_triangle_mesh((geogram::coord_index_t) 3, vertices, triangles, true);

    GEO::MeshIOFlags export_flags;
    export_flags.set_dimension(3);
    export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    export_flags.set_verbose(true);
    geogram::mesh_save(mesh, std::vformat(arguments.output, std::make_format_args(arguments.first, arguments.amount + arguments.first)), export_flags);

    if (!arguments.f_generate_sizing_field)
    {
        return 0;
    }

    /*geogram::compute_sizing_field(mesh);
    geogram::mesh_save(
        mesh,
        std::filesystem::path(std::vformat(arguments.output, std::make_format_args(arguments.first, arguments.amount + arguments.first))).replace_extension(".sf.geogram"),
        export_flags
    );*/

    return 0;
}
