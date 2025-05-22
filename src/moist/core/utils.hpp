#ifndef MOIST_CORE_UTILS_HPP_
#define MOIST_CORE_UTILS_HPP_

#include <sstream>
#include <string>
#include <unordered_set>

#ifdef GEOGRAM_API
#include <filesystem>

#include <geogram/basic/environment.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#endif // GEOGRAM_API

namespace moist::utils {
#ifdef GEOGRAM_API
    namespace geo
    {
        inline geogram::vec3 parse_vector(std::string str)
        {
            double x, y, z;
            char none;

            std::istringstream stream(str);
            stream >> none >> x >> none >> y >> none >> z >> none;

            return geogram::vec3(x, y, z);
        }

        inline void initialize()
        {
            geogram::initialize(geogram::GEOGRAM_INSTALL_NONE);
            geogram::CmdLine::import_arg_group("sys"); // needs to be called in order to be able to export .geogram meshes...
            geogram::CmdLine::set_arg("sys:compression_level", "0");
            geogram::Logger::instance()->set_quiet(true);
        }

        inline void load(const std::filesystem::path& file, geogram::Mesh& mesh, const geogram::index_t dimension = 0, const bool attributes = true)
        {
            geogram::MeshIOFlags flags;
            flags.set_dimension(dimension != 0 ? dimension : mesh.vertices.dimension());
            flags.set_attribute(attributes ? geogram::MESH_ALL_ATTRIBUTES : geogram::MESH_NO_ATTRIBUTES);
            flags.set_elements(geogram::MeshElementsFlags(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS));
            flags.set_verbose(false);

            if (!geogram::mesh_load(file.string(), mesh, flags))
            {
                // TODO: Exception!
            }
        }

        inline void save(const std::filesystem::path& file, const geogram::Mesh& mesh, bool attributes = true)
        {
            GEO::MeshIOFlags flags;
            flags.set_dimension(mesh.vertices.dimension());
            flags.set_attribute(attributes ? geogram::MESH_ALL_ATTRIBUTES : geogram::MESH_NO_ATTRIBUTES);
            flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
            flags.set_verbose(false);
            geogram::mesh_save(mesh, file.string(), flags);
        }
    }
#endif // GEOGRAM_API

    template <typename T, typename H>
    inline std::unordered_set<T, H> symmetric_difference(const std::unordered_set<T, H>& a, const std::unordered_set<T, H>& b)
    {
        std::unordered_set<T, H> result;

        for (const auto& item : a)
        {
            if (b.find(item) == b.end())
            {
                result.insert(item);
            }
        }
        for (const auto& item : b)
        {
            if (a.find(item) == a.end())
            {
                result.insert(item);
            }
        }

        return result;
    }

    inline std::size_t count_placeholders(const std::string& pattern)
    {
        std::size_t count = 0;
        for (std::size_t i = 0; i < pattern.size(); ++i)
        {
            if (pattern[i] == '{')
            {
                if (i + 1 < pattern.size() && pattern[i + 1] == '{')
                {
                    ++i;
                }
                else
                {
                    ++count;
                }
            }
            else if (pattern[i] == '}')
            {
                if (i + 1 < pattern.size() && pattern[i + 1] == '}')
                {
                    ++i;
                }
            }
        }
        return count;
    }
}

#endif // MOIST_CORE_UTILS_HPP_
