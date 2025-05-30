#ifndef MOIST_SUBMESH_MERGER_SUBMESH_MERGER_HPP_
#define MOIST_SUBMESH_MERGER_SUBMESH_MERGER_HPP_

#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <geogram/basic/line_stream.h>
#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

namespace moist
{
    class Merger
    {
    public:
        Merger(const std::filesystem::path& mesh, const std::filesystem::path& slice, const std::string& tmp_pattern = "./.tmp{}");
        ~Merger();

        void Merge();
        void CopyToOriginal(const std::filesystem::path& destination = "");
    private:
        std::filesystem::path _tmp_folder_path;
        std::filesystem::path _v_swap_path;
        std::filesystem::path _e_swap_path;

        std::filesystem::path _mesh_path;

        geogram::index_t _nb_vertices_mesh;
        std::unordered_map<g_index, g_index> _v_translated_indices;
        std::unordered_map<g_index, g_index> _v_inverse_prefix_offset;

        geogram::Mesh _slice;

        geogram::index_t WriteNodes(geogram::LineInput& input);
        geogram::index_t WriteElements(geogram::LineInput& input);
    };
}

#endif // MOIST_SUBMESH_MERGER_SUBMESH_MERGER_HPP_
