#include <limits>

#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

#include "moist/core/defines.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "interface_generator.hpp"

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

    moist::utils::geo::initialize();

    auto generator = moist::InterfaceGenerator(moist::AxisAlignedInterfacePlane {
        arguments.axis,
        arguments.plane,
        arguments.envelope_size
    });

    geogram::Mesh mesh, mesh2;
    moist::utils::geo::load(arguments.path_mesh_a, mesh);
    geogram::mesh_repair(mesh, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);
    generator.AddConstraints(mesh);
    mesh.clear(false, false);
    moist::utils::geo::load(arguments.path_mesh_b, mesh);
    geogram::mesh_repair(mesh, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);
    generator.AddConstraints(mesh);
    generator.Triangulate();
    //generator.Decimate();

    // inject "quality" attribute on facets for later use...
    // geogram::Attribute<double> f_quality(generator.Triangulation()->facets.attributes(), ATTRIBUTE_INTERFACE_TETMERGE_QUALITY);
    // for (const g_index f : generator.Triangulation()->facets)
    // {
    //     f_quality[f] = -std::numeric_limits<double>::max();
    // }

    moist::utils::geo::save(arguments.path_mesh_out.replace_extension(".geogram"), *generator.Triangulation());
    moist::utils::geo::save(arguments.path_mesh_out.replace_extension(".debug.msh"), *generator.Triangulation());

    return 0;
}
