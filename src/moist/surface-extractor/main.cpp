#include <format>

// works only on linux
#include <sys/resource.h>

#include <CLI/CLI.hpp>

#include <MC33.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/metrics.hpp"

#include "tiff_data.hpp"
#include "surface_generator.hpp"
// #include "sizing-field/scalar_field.hpp"

struct Arguments
{
    std::string input;
    std::string output;

    uint32_t first;
    uint32_t amount;
    uint32_t sample_bits;

    moist::Axis axis;

    int32_t dir_offset;
    float isovalue;
    double scale_factor;
    bool f_center;
    bool f_invert;
    bool f_generate_sizing_field;

    float sizing_field_scale_denominator;

    std::filesystem::path csv;
    std::string csv_name;
};

static cli::Validator PatternPlaceholderCountValidator(const std::size_t placeholders)
{
    return cli::Validator([placeholders](const std::string &input) -> std::string {
        const size_t amount = moist::utils::count_placeholders(input);
        return (amount == placeholders)
            ? ""
            : std::format("invalid pattern {}: expected {} placeholders, got {}", input, placeholders, amount);
    }, "PATTERN");
}

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-i, --input", arguments.input, "CPP pattern describing input files. Arguments: [number:uint32_t] (e.g. data/{:03}.tif).")
        ->required()
        ->check(PatternPlaceholderCountValidator(1));
    app.add_option("-o, --output", arguments.output, "CPP pattern describing the output file. Arguments: [from:uint32|to:uint32t] (e.g. slice{}-{}.obj).")
        ->required()
        ->check(PatternPlaceholderCountValidator(2));
    app.add_option("-f, --first", arguments.first, "Initial file number.")
        ->check(cli::NonNegativeNumber)
        ->required();
    app.add_option("-n, --amount", arguments.amount, "Total number of files to mesh into this slice.")
        ->check(cli::PositiveNumber)
        ->required();
    app.add_option("-v, --isovalue", arguments.isovalue, "Isovalue.")
        //->check(cli::PositiveNumber)
        //->check(cli::Bound(0.0f, 1.0f))
        ->required();
    app.add_option("-b, --sample-bits", arguments.sample_bits, "Bits per sample [8, 16, 32, 64].")
        ->check(cli::IsMember({8, 16, 32, 64}))
        ->default_val(16);
    app.add_option("-x, --axis", arguments.axis, "Interface axis (X|Y|Z)")
        ->transform(cli::CheckedTransformer(moist::AXIS_OPTION_ARGUMENT_MAP, cli::ignore_case))
        ->default_val(moist::Axis::Z); // TODO: Maybe a "growth-direction" (positive/negative)?
    app.add_option("--directional-offset", arguments.dir_offset, "Directional offset along the growth axis (Z), e.g. if your files start with File500.tif, set this to -500.")
        ->default_val(0);
    app.add_flag("--center", arguments.f_center, "Center around 0,0,0. Applies to the entire mesh, not just this slice.")
        ->default_val(false);
    app.add_flag("--invert", arguments.f_invert, "Invert output (mesh negative space).")
        ->default_val(true);
    app.add_flag("--generate-sizing-field", arguments.f_generate_sizing_field, "Generates a sizing field - TODO: More configs for this")
        ->default_val(false);

    // sizing-field specific, how do I make this conditional if f_generate_sizing_field is set?
    app.add_option("--sizing-field-scale-denominator", arguments.sizing_field_scale_denominator, "Inverse minimum sizing field octree node size, in relation to the width of the data. E.g., with a width of 500, setting this to 50 will set the minimum size of an octree cell to 500/50 = 10x10x10px")
        ->default_val(100.0f);

    app.add_option("--csv", arguments.csv, "CSV file for metrics")
        ->required();
    app.add_option("--run-name", arguments.csv_name, "Run name for metrics")
        ->required();

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geogram::initialize();
    // geo::CmdLine::import_arg_group("algo");

    auto metrics = moist::metrics::Metrics(arguments.csv_name);
    moist::Timer runtime_timer("Total", metrics);

    moist::Timer tiff_timer("Tiff::Tiff", metrics);
    const moist::Tiff tiff_data(arguments.input, arguments.first, arguments.amount);
    tiff_timer.End();
    *metrics << moist::metrics::Metric("volume_x", tiff_data.width());
    *metrics << moist::metrics::Metric("volume_y", tiff_data.height());
    *metrics << moist::metrics::Metric("volume_z", tiff_data.depth());
    *metrics << moist::metrics::Metric("volume_xyz", tiff_data.width() * tiff_data.height() * tiff_data.depth());

    moist::SurfaceGenerator generator(tiff_data, int(arguments.first) + arguments.dir_offset, arguments.axis, arguments.f_center, arguments.f_invert);

    geo::Mesh mesh(3);
    {
        moist::Timer timer("SurfaceGenerator::Generate", metrics);
        generator.Generate(mesh, arguments.isovalue);
    }

    {
        moist::Timer timer("SurfaceGenerator::Save", metrics);
        const std::string path = std::vformat(arguments.output, std::make_format_args(arguments.first, arguments.amount + arguments.first - 1));
        moist::utils::geogram::save(std::filesystem::path(path), mesh);
        const std::uintmax_t size = std::filesystem::file_size(std::filesystem::path(path));
        double size_mb = static_cast<double>(size) / (1024 * 1024);
        MOIST_INFO("saved " << path << ", filesize = " << size_mb << "mb");
        *metrics << moist::metrics::Metric("filesize_mb", size_mb);
        *metrics << moist::metrics::Metric("num_triangles", mesh.facets.nb());
        *metrics << moist::metrics::Metric("num_vertices", mesh.vertices.nb());
    }

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    *metrics << moist::metrics::Metric("peak_ram_usage_mb", static_cast<double>(usage.ru_maxrss) / 1024);

    runtime_timer.End();
    metrics->AppendCSV(arguments.csv);

    if (!arguments.f_generate_sizing_field)
    {
        return 0;
    }

    // moist::generate_scalar_field(tiff_data, moist::IsoValue(arguments.isovalue), metrics);

    // TODO: Generate Sizing Field Mesh...
    return 0;
}
