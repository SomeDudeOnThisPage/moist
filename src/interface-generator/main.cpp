#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include "core.hpp"
#include "interface-generator.hpp"

static bool load_mesh(const std::filesystem::path& path, incremental_meshing::InterfaceGenerator& generator)
{
    geogram::MeshIOFlags import_flags;
    import_flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_ALL_ELEMENTS ||
        geogram::MeshElementsFlags::MESH_ALL_SUBELEMENTS
    ));

    geogram::Mesh mesh;
    if (!geogram::mesh_load(path, mesh))
    {
        OOC_ERROR("Failed to load mesh A: " << path);
        return false;
    }
    generator.AddConstraints(mesh);
    return true;
}

/**
 * Parameters:
 *  a {path}: Tetrahedal Mesh A (.msh/.mesh format)
 *  b {path}: Tetrahedal Mesh B (.msh/.mesh format)
 *  o {path}: Output-Triangle-Mesh containing the interface (geogram mesh format with metadata)
 *  axis {x|y|z}: Axis the interface plane lies on
 *  plane {double}: Extent of the interface plane
 *  eps {double, optional}: Epsilon to snap vertices from/to
 *
 * Output:
 *  Interface-Mesh in geogram format, with additional metadata.
 *
 * Description:
 *  Loads meshes sequentially, gathering all vertices on the interface plane, outputs a (constrained) triangulation of these vertices.
 *  Memory requirements are estimated as (max(mem(a), mem(b)) + max(mem(a), mem(b)) / slice_amount), where slice_amount is the amount of
 *  the max amount of slices of CT-data used to generate mesh a and b.
 */
int main(int argc, char* argv[])
{
    incremental_meshing::InterfaceExtractionOptions options;
    const std::map<std::string, incremental_meshing::Axis> axis_option_map{
        {"X", incremental_meshing::Axis::X},
        {"Y", incremental_meshing::Axis::Y},
        {"Z", incremental_meshing::Axis::Z}
    };

    CLI::App app{argv[0]};
    app.add_option("-a, --input-a", options.path_mesh_a, "Mesh 'A' (.mesh) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-b, --input-b", options.path_mesh_b, "Mesh 'B' (.mesh) file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o, --output", options.path_mesh_out, "Output Triangulation (.obj) file")
        ->required();
    app.add_option("-x, --axis", options.axis, "Interface axis (X|Y|Z)")
        ->required()
        ->transform(CLI::CheckedTransformer(axis_option_map, CLI::ignore_case));
    app.add_option("-p, --plane", options.plane, "Interface plane along the axis (-x)")
        ->required();
    app.add_option("-e, --envelope", options.envelope_size, "Envelope size around the interface plane.");

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
    geogram::Logger::instance()->set_quiet(false);

    auto generator = incremental_meshing::InterfaceGenerator(incremental_meshing::AxisAlignedInterfacePlane {
        options.axis,
        options.plane,
        options.envelope_size
    });

    if (!load_mesh(options.path_mesh_a, generator))
    {
        return -1;
    }

    if (!load_mesh(options.path_mesh_b, generator))
    {
        return -1;
    }

    generator.Triangulate();

    GEO::MeshIOFlags export_flags;
    export_flags.set_dimension(3);
    export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
    export_flags.set_verbose(true);
    geogram::mesh_save(*generator.Triangulation(), options.path_mesh_out.string(), export_flags);

    return 0;
}
