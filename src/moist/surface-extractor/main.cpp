#include <format>

#include <CLI/CLI.hpp>

#include <MC33.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/utils.hpp"

#include "tiff_data.hpp"
#include "surface_generator.hpp"
#include "sizing-field/scalar_field.hpp"

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
    auto metrics = moist::metrics::TimeMetrics("SurfaceExtractor");
    moist::Timer timer("SurfaceExtractor::Main", metrics);

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
        ->check(cli::PositiveNumber)
        ->check(cli::Bound(0.0f, 1.0f))
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

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geo::initialize();

    const moist::Tiff tiff_data(arguments.input, arguments.first, arguments.amount);
    moist::SurfaceGenerator generator(tiff_data, arguments.first + arguments.dir_offset, arguments.axis, arguments.f_center, arguments.f_invert);

    geogram::Mesh mesh(3);
    {
        moist::Timer _scope_timer("SurfaceGenerator::Generate", metrics);
        generator.Generate(mesh, arguments.isovalue);
    }

    geogram::mesh_repair(mesh);

    const std::string path = std::vformat(arguments.output, std::make_format_args(arguments.first, arguments.amount + arguments.first - 1));
    moist::utils::geo::save(std::filesystem::path(path), mesh);

    if (!arguments.f_generate_sizing_field)
    {
        return 0;
    }

    moist::generate_scalar_field(tiff_data, moist::IsoValue(arguments.isovalue), metrics);

    // TODO: Generate Sizing Field Mesh...
    return 0;
}
