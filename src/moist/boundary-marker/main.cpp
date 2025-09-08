#include <limits>

#include <CLI/CLI.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

#include <moist/core/defines.hpp>
#include <moist/core/utils.hpp>

enum BoundaryMarker : int
{
    SOURCE = 1,
    SINK = 2,
    BOUNDARY = 3
};

struct Arguments
{
    std::filesystem::path mesh;
    std::filesystem::path mesh_out;
    double z_in;
    double z_out;
    double epsr;
};

// this is super scuffed but works...
// really we'd like to make the phys names configurable, etc.
void save_with_physgroups(const geo::Mesh& mesh, const std::filesystem::path& path)
{
    std::ofstream os(path);
    if(!os) return;
    os << std::fixed << std::setprecision(16);

    os << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    os << "$PhysicalNames\n1\n";
    os << "2 1 \"source\"\n";
    os << "2 2 \"sink\"\n";
    os << "2 3 \"boundary\"\n";
    os << "3 4 \"fluid\"\n";
    os << "$EndPhysicalNames\n";

    geo::index_t nb_vertices = mesh.vertices.nb();
    os << "$Nodes\n" << nb_vertices << "\n";
    for(geo::index_t v = 0; v < nb_vertices; v++)
    {
        const double* p = mesh.vertices.point_ptr(v);
        os << (v + 1) << " " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    os << "$EndNodes\n";

    geo::index_t elem_count = mesh.facets.nb() + mesh.cells.nb();

    os << "$Elements\n" << elem_count << "\n";
    geo::index_t id = 1;
    geo::Attribute<int> boundary_marker(mesh.facets.attributes(), "boundary_marker");
    for(geo::index_t f = 0; f < mesh.facets.nb(); f++)
    {
        if(mesh.facets.nb_vertices(f) != 3) continue;

        int phys_tag = boundary_marker[f];
        geo::index_t v0 = mesh.facets.vertex(f, 0) + 1;
        geo::index_t v1 = mesh.facets.vertex(f, 1) + 1;
        geo::index_t v2 = mesh.facets.vertex(f, 2) + 1;

        os << id++ << " 2 2 " << phys_tag << " " << phys_tag << " " << v0 << " " << v1 << " " << v2 << "\n";
    }

    for(geo::index_t c = 0; c < mesh.cells.nb(); c++)
    {
        int phys = 4; // when we have at least one physgroup, every element needs one, even cells...

        geo::index_t v0 = mesh.cells.vertex(c,0) + 1;
        geo::index_t v1 = mesh.cells.vertex(c,1) + 1;
        geo::index_t v2 = mesh.cells.vertex(c,2) + 1;
        geo::index_t v3 = mesh.cells.vertex(c,3) + 1;

        os << id << " 4 2 " << phys << " " << phys
           << " " << v0 << " " << v1 << " " << v2 << " " << v3 << "\n";
        id++;
    }
    os << "$EndElements\n";

    os.close();
    MOIST_INFO("saved mesh to " << path);
}

void collect_boundary_facets(const geo::Mesh& mesh, std::vector<geo::index_t>& boundary_facets)
{
    boundary_facets.clear();
    for(geo::index_t c = 0; c < mesh.cells.nb(); c++)
    {
        geo::index_t nlf = mesh.cells.nb_facets(c); // should be 4 for a tet
        for(geo::index_t lf = 0; lf < nlf; ++lf) {
            if (mesh.cells.adjacent(c, lf) == geo::NO_CELL)
            {
                geo::index_t f = mesh.cells.facet(c, lf);
                boundary_facets.push_back(f);
            }
        }
    }
}

inline bool in_eps(double x, double ref, double eps)
{
    return (x >= ref - eps) && (x <= ref + eps);
}

int main(int argc, char* argv[])
{
    Arguments arguments{};
    cli::App app{argv[0]};

    app.add_option("-m, --mesh", arguments.mesh, "Mesh (.msh) file")
        ->required()
        ->check(cli::ExistingFile);
    app.add_option("-x, --mesh-out", arguments.mesh_out, "Mesh output (.msh) file")
        ->required();
    app.add_option("-i, --input-z", arguments.z_in)
        ->required();
    app.add_option("-o, --output-z", arguments.z_out)
        ->required();
    app.add_option("-e, --epsr", arguments.epsr, "Epsilon, in case the mesh was generated via envelope meshing")
        ->default_val(0.0);

    CLI11_PARSE(app, argc, app.ensure_utf8(argv));

    moist::utils::geogram::initialize();
    geo::Mesh mesh(3);
    moist::utils::geogram::load(arguments.mesh, mesh, 3, true);
    mesh.facets.clear(false, false);
    mesh.cells.connect();
    mesh.cells.compute_borders();

    geo::Attribute<int> boundary_marker(mesh.facets.attributes(), "boundary_marker");
    for (const geo::index_t f : mesh.facets)
    {
        const geo::vec3 p0 = mesh.vertices.point(mesh.facets.vertex(f, 0));
        const geo::vec3 p1 = mesh.vertices.point(mesh.facets.vertex(f, 1));
        const geo::vec3 p2 = mesh.vertices.point(mesh.facets.vertex(f, 2));

        if (in_eps(p0.z, arguments.z_in, arguments.epsr) && in_eps(p1.z, arguments.z_in, arguments.epsr) && in_eps(p2.z, arguments.z_in, arguments.epsr))
        {
            boundary_marker[f] = BoundaryMarker::SOURCE;
        }
        else if (in_eps(p0.z, arguments.z_out, arguments.epsr) && in_eps(p1.z, arguments.z_out, arguments.epsr) && in_eps(p2.z, arguments.z_out, arguments.epsr))
        {
            boundary_marker[f] = BoundaryMarker::SINK;
        }
        else
        {
            boundary_marker[f] = BoundaryMarker::BOUNDARY;
        }
    }

    save_with_physgroups(mesh, arguments.mesh_out);

    return EXIT_SUCCESS;
}
