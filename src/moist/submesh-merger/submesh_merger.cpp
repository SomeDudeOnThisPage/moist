#include "submesh_merger.hpp"

#include <format>

#include "moist/core/utils.hpp"

moist::Merger::Merger(const std::filesystem::path& mesh, const std::filesystem::path& slice, const std::string& tmp_pattern) : _mesh_path(mesh)
{
    _tmp_folder_path = std::filesystem::path(std::vformat(tmp_pattern, std::make_format_args(mesh.filename().string())));

    _v_swap_path = std::filesystem::path(_tmp_folder_path / "vertices");
    _e_swap_path = std::filesystem::path(_tmp_folder_path / "elements");

    moist::utils::geogram::load(slice, _slice, 3);

    if (std::filesystem::exists(_tmp_folder_path))
    {
        OOC_ERROR("cannot merge multiple slices into the same mesh simultaneously");
    }

    std::filesystem::create_directory(_tmp_folder_path);
}

moist::Merger::~Merger()
{
    std::filesystem::remove_all(_tmp_folder_path);
}

void moist::Merger::Merge()
{
    geo::LineInput input(_mesh_path);
    geo::index_t nb_nodes = this->WriteNodes(input);
    geo::index_t nb_elements = this->WriteElements(input);

    const auto tmp_msh_path = _tmp_folder_path / "merge";
    std::ofstream msh_tmp(tmp_msh_path);
    msh_tmp.precision(16);
    msh_tmp << "$MeshFormat" << std::endl;
    msh_tmp << "2.2 0 8" << std::endl;
    msh_tmp << "$EndMeshFormat" << std::endl;
    msh_tmp << "$Nodes" << std::endl;
    msh_tmp << nb_nodes << std::endl;

    geo::LineInput v_in(_v_swap_path);
    geo::LineInput e_in(_e_swap_path);

    geo::index_t i = 1;
    while (v_in.get_line())
    {
        msh_tmp << i++ << " " << v_in.current_line();
    }

    if (i - 1 != nb_nodes)
    {
        OOC_DEBUG("wrong number of nodes: expected " << nb_nodes << ", got " << i);
    }

    msh_tmp << "$EndNodes" << std::endl;
    msh_tmp << "$Elements" << std::endl;
    msh_tmp << nb_elements << std::endl;

    i = 1;
    while (e_in.get_line())
    {
        msh_tmp << i++ << " " << e_in.current_line();
    }

    if (i - 1 != nb_elements)
    {
        OOC_DEBUG("wrong number of elements: expected " << nb_elements << ", got " << i);
    }


    msh_tmp << "$EndElements" << std::endl;
    msh_tmp.close();
}

void moist::Merger::CopyToOriginal(const std::filesystem::path &destination_override)
{
    std::filesystem::path destination = destination_override.string().empty()
        ? _mesh_path
        : destination_override;

    destination.replace_extension(".msh");
    std::filesystem::copy_file(_tmp_folder_path / "merge", destination, std::filesystem::copy_options::overwrite_existing);
#ifndef NDEBUG
    //load/save as .mesh to quickly validate with tetgen if it is a valid mesh
    geo::Mesh dbg(3);
    moist::utils::geogram::load(destination, dbg);
    moist::utils::geogram::save("custom_merge.mesh", dbg);
#endif // NDEBUG
}

geo::index_t moist::Merger::WriteNodes(geo::LineInput& input)
{
    geo::index_t nb_vertices;
    std::ofstream v_swapfile(_v_swap_path);
    v_swapfile.precision(16);

    while (input.get_line())
    {
        input.get_fields();
        if (input.field_matches(0, "$Nodes"))
        {
            input.get_line();
            input.get_fields();
            nb_vertices = input.field_as_uint(0);
            _nb_vertices_mesh = nb_vertices;

            for (g_index v = 1; v < nb_vertices + 1; v++)
            {
                input.get_line();
                input.get_fields();
                const vec3 point(input.field_as_double(1), input.field_as_double(2), input.field_as_double(3));

                //geogram::parallel_for(0, slice.vertices.nb(), [v, &point, &slice, &v_s2m](const g_index _v)
                //{
                for (const auto _v : _slice.vertices)
                {
                    if (point == _slice.vertices.point(_v))
                    {
                        if (_v_translated_indices.contains(_v))
                        {
                            OOC_DEBUG("invalid point mapping: multiple definitions of v = [" << point.x << ", " << point.y << ", " << point.z << "] - mesh has colocated vertices within DOUBLE_EPSILON");
                        }

                        _v_translated_indices[_v] = v;
                    }
                }
                //});

                v_swapfile << point.x << " " << point.y << " " << point.z << std::endl;
            }
        }
        else if (input.field_matches(0, "$EndNodes"))
        {
            size_t offset = 0;
            for (const g_index v : _slice.vertices)
            {
                if (_v_translated_indices.contains(v))
                {
                    // vertex is translated, vertex will not be written to final file, increment inverse prefix offset
                    _v_inverse_prefix_offset[v] = ++offset;
                }
                else
                {
                    // vertex is not translated, vertex will be written to final file
                    _v_inverse_prefix_offset[v] = offset;

                    const vec3 point = _slice.vertices.point(v);
                    v_swapfile << point.x << " " << point.y << " " << point.z << std::endl;
                    nb_vertices++;
                }
            }

            break;
        }
    }

    v_swapfile.close();

    return nb_vertices;
}

geo::index_t moist::Merger::WriteElements(geo::LineInput& input)
{
    geo::index_t nb_elements = 0;
    std::ofstream e_swapfile(_e_swap_path);
    e_swapfile.precision(16);

    while (input.get_line())
    {
        input.get_fields();
        if (input.field_matches(0, "$Elements"))
        {
            input.get_line();
            input.get_fields();
            nb_elements = input.field_as_uint(0);

            for (g_index e = 0; e < nb_elements; e++)
            {
                input.get_line();
                input.get_fields();

                for (geo::index_t i = 1; i < input.nb_fields(); i++)
                {
                    e_swapfile << input.field(i) << " ";
                }
                e_swapfile << std::endl;
            }
        }
        else if (input.field_matches(0, "$EndElements"))
        {
            const geo::index_t offset = nb_elements;
            for (const g_index c : _slice.cells)
            {
                g_index vertices[4];
                for (l_index lv = 0; lv < 4; lv++)
                {
                    g_index v = _slice.cells.vertex(c, lv);

                    if (_v_translated_indices.contains(v))
                    {
                        vertices[lv] = _v_translated_indices[v];
                    }
                    else
                    {
                        v -= _v_inverse_prefix_offset[v];
                        vertices[lv] = v + _nb_vertices_mesh + /* msh is 1-indexed, vs geogram 0-based index */ 1;
                    }
                }

                nb_elements++;
                e_swapfile << "4 1 0 "; // only write no-attribute cells for now
                e_swapfile << vertices[0] << " " << vertices[1] << " " << vertices[2] << " " << vertices[3] << std::endl;
            }
        }
    }

    e_swapfile.close();

    return nb_elements;
}
