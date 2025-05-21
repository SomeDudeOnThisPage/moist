#include <limits>

#include <CLI/CLI.hpp>

#include <geogram/basic/environment.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include "moist/core/defines.hpp"
#include "moist/core/attributes.inl"
#include "interface_generator.hpp"

static bool load_mesh(const std::filesystem::path& path, moist::InterfaceGenerator& generator)
{
    geogram::MeshIOFlags import_flags;
    import_flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS ||
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS
    ));

    geogram::Mesh mesh;
    if (!geogram::mesh_load(path, mesh))
    {
        OOC_ERROR("Failed to load mesh: " << path);
        return false;
    }
    generator.AddConstraints(mesh);
    return true;
}

int main(int argc, char* argv[])
{
    struct Arguments
    {
        std::filesystem::path path_mesh_a;
        std::filesystem::path path_mesh_b;
        std::filesystem::path path_mesh_out;
        double plane;
        double envelope_size;
        moist::Axis axis;
    };

    const std::map<std::string, moist::Axis> axis_option_map{
        {"X", moist::Axis::X},
        {"Y", moist::Axis::Y},
        {"Z", moist::Axis::Z}
    };

    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-a, --input-a", arguments.path_mesh_a, "Mesh 'A' (.mesh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-b, --input-b", arguments.path_mesh_b, "Mesh 'B' (.mesh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-o, --output", arguments.path_mesh_out, "Output Triangulation (.obj) file")
        ->required();
    app.add_option("-x, --axis", arguments.axis, "Interface axis (X|Y|Z)")
        ->required()
        ->default_str("Y")
        ->transform(cli::CheckedTransformer(axis_option_map, cli::ignore_case));
    app.add_option("-p, --plane", arguments.plane, "Interface plane along the axis (-x)")
        ->required();
    app.add_option("-e, --envelope", arguments.envelope_size, "Envelope size around the interface plane.")
        ->default_val(0.00001);

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize();
    geogram::CmdLine::import_arg_group("sys"); // needs to be called in order to be able to export .geogram meshes...
    geogram::CmdLine::set_arg("sys:compression_level", "0");
    geogram::Logger::instance()->set_quiet(false);

    auto generator = moist::InterfaceGenerator(moist::AxisAlignedInterfacePlane {
        arguments.axis,
        arguments.plane,
        arguments.envelope_size
    });

    if (!load_mesh(arguments.path_mesh_a, generator))
    {
        return -1;
    }

    if (!load_mesh(arguments.path_mesh_b, generator))
    {
        return -1;
    }

    generator.Triangulate();

    // inject "quality" attribute on facets for later use...
    geogram::Attribute<double> f_quality(generator.Triangulation()->facets.attributes(), ATTRIBUTE_INTERFACE_TETMERGE_QUALITY);
    for (const g_index f : generator.Triangulation()->facets)
    {
        f_quality[f] = -std::numeric_limits<double>::max();
    }

    GEO::MeshIOFlags export_flags;
    export_flags.set_dimension(3);
    export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    export_flags.set_verbose(true);
    geogram::mesh_save(*generator.Triangulation(), arguments.path_mesh_out.replace_extension(".geogram"), export_flags);

#ifndef NDEBUG
    geogram::mesh_save(*generator.Triangulation(), arguments.path_mesh_out.replace_extension(".debug.obj"), export_flags);
#endif

    return 0;
}
