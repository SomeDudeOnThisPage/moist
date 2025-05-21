#include "slice_io.hpp"

#include <geogram/basic/line_stream.h>
#include <geogram/mesh/mesh_io.h>

struct MshMoistSliceLocalMetadata
{
    size_t v;
    size_t e;
    std::vector<double/*acts as spatial index "extent"*/> interfaces;
};

struct MshMoistSliceMetadata
{
    size_t v;
    size_t e;
    size_t v_b_start;
    size_t v_e_start;
    std::vector<double/*acts as spatial index "extent"*/> interfaces;
};

struct MshMoistMetadata
{
    size_t v_total;
    size_t e_total;
    std::vector<MshMoistSliceMetadata> slices;
};

// modified from geogram meshio
void moist::slice_io::msh::save(const std::filesystem::path file, const geogram::Mesh& slice, const geogram::MeshIOFlags ioflags)
{
    const geogram::index_t id_offset_msh = 1;
    const geogram::index_t msh2geo_hex[8] = {1, 3, 7, 5, 0, 2, 6, 4 };
    const geogram::index_t msh2geo_def[8] = {0, 1, 2, 3, 4, 5, 6, 7 };
    const geogram::index_t celltype_geo2msh[5] = {4, 5, 6, 7};

    geogram::Mesh M(slice.vertices.dimension());
    M.copy(slice, true);

    M.vertices.remove_isolated();

    geogram::Attribute<int> region;
    geogram::Attribute<int> bdr_region;
    if (M.cells.attributes().is_defined("region"))
    {
        region.bind(M.cells.attributes(), "region");
    }
    if (M.facets.attributes().is_defined("bdr_region"))
    {
        bdr_region.bind(M.facets.attributes(), "bdr_region");
    }

    std::ofstream out(file.c_str()) ;

    if(!out)
    {
        OOC_ERROR("Fail to open \"" << file << "\" for writing");
        return;
    }

    out.precision(16);

    /* Header */
    out << "$MeshFormat\n";
    out << "2.2 0 " << sizeof(double) << std::endl;
    out << "$EndMeshFormat\n";

    out << "$MoistMetadata\n";
    out << "$EndMoistMetadata\n";

    /* Vertices */
    out << "$Nodes" << std::endl ;
    out << M.vertices.nb() << std::endl ;
    for(geogram::index_t v = 0; v < M.vertices.nb(); v++)
    {
        out << v + id_offset_msh << " "
            << M.vertices.point_ptr(v)[0] << " "
            << M.vertices.point_ptr(v)[1] << " "
            << M.vertices.point_ptr(v)[2] << '\n' ;
    }
    out << "$EndNodes" << std::endl ;

    /* Elements */
    geogram::index_t nb_tet = 0;
    geogram::index_t nb_hex = 0;
    geogram::index_t nb_pyr = 0;
    geogram::index_t nb_pri = 0;
    for(geogram::index_t c = 0; c != M.cells.nb(); ++c)
    {
        if(M.cells.type(c) == GEO::MESH_TET)
        {
            ++nb_tet;
        }
        else if(M.cells.type(c) == GEO::MESH_HEX)
        {
            ++nb_hex;
        }
        else if(M.cells.type(c) == GEO::MESH_PYRAMID)
        {
            ++nb_pyr;
        }
        else if(M.cells.type(c) == GEO::MESH_PRISM)
        {
            ++nb_pri;
        }
    }
    geogram::index_t nb_elt = nb_tet + nb_hex + nb_pyr + nb_pri;
    if (ioflags.has_element(geogram::MeshElementsFlags::MESH_FACETS)) nb_elt += M.facets.nb();

    out << "$Elements" << std::endl ;
    out << nb_elt << std::endl ;
    geogram::index_t elt_id = 0; /* starts at 1, common for faces and cells */
    if (ioflags.has_element(geogram::MeshElementsFlags::MESH_FACETS))
    {
        for (geogram::index_t f = 0; f < M.facets.nb(); ++f)
        {
            int attr_value = 0;
            if (bdr_region.is_bound())
            {
                attr_value = bdr_region[f];
            }
            int type = -1;
            if (M.facets.nb_vertices(f) == 3)
            {
                type = 2;
            } else if (M.facets.nb_vertices(f) == 4)
            {
                type = 3;
            } else
            {
                geo_assert_not_reached
            }

            elt_id += 1;
            out << elt_id << " " << type << " " << "2" << " "
                << attr_value << " " << attr_value << " ";
            for (geogram::index_t li = 0; li < M.facets.nb_vertices(f); ++li)
            {
                out << M.facets.vertex(f, li) + id_offset_msh << " ";
            }
            out << std::endl;
        }
    }
    for (geogram::index_t c = 0; c < M.cells.nb(); ++c)
    {
        if (M.cells.type(c) == GEO::MESH_CONNECTOR)
        {
            continue;
        }
        int attr_value = 0;
        if (region.is_bound()) attr_value = region[c];
        const geogram::index_t* msh2geo =
            (M.cells.type(c) == GEO::MESH_HEX) ?
            msh2geo_hex : msh2geo_def;
        elt_id += 1;

        /* Write to file, format is:
            *   elm-number elm-type number-of-tags < tag > ...
            *   node-number-list
            */
        out << elt_id << " " << celltype_geo2msh[M.cells.type(c)]
            << " " << "1" << " " << attr_value << " ";
        for (geogram::index_t li = 0; li < M.cells.nb_vertices(c); ++li)
        {
            out << M.cells.vertex(c, msh2geo[li]) + id_offset_msh
                << " ";
        }
        out << std::endl;
    }
    out << "$EndElements" << std::endl;

    out.close();
}

/*void moist::slice_io::msh::read_metadata(const std::filesystem::path file, const geogram::Mesh& slice)
{

}*/
