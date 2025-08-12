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

#include "moist/core/defines.hpp"

namespace moist::utils {
#ifdef GEOGRAM_API
    namespace geogram
    {
        inline geo::vec3 parse_vector(std::string str)
        {
            double x, y, z;
            char none;

            std::istringstream stream(str);
            stream >> none >> x >> none >> y >> none >> z >> none;

            return geo::vec3(x, y, z);
        }

        inline void initialize()
        {
            geo::initialize(geo::GEOGRAM_INSTALL_NONE);
            geo::CmdLine::import_arg_group("sys"); // needs to be called in order to be able to export .geogram meshes...
            geo::CmdLine::set_arg("sys:compression_level", "0");
            // geogram::Logger::instance()->set_quiet(false);
        }

        inline void load(const std::filesystem::path& file, geo::Mesh& mesh, const geo::index_t dimension = 0, const bool attributes = true)
        {
            geo::MeshIOFlags flags;
            flags.set_dimension(dimension != 0 ? dimension : mesh.vertices.dimension());
            flags.set_attribute(attributes ? geo::MESH_ALL_ATTRIBUTES : geo::MESH_NO_ATTRIBUTES);
            flags.set_elements(geo::MeshElementsFlags(geo::MeshElementsFlags::MESH_ALL_ELEMENTS));
            flags.set_verbose(true);

            if (!geo::mesh_load(file.string(), mesh, flags))
            {
                OOC_ERROR("Failed to load mesh '" << file << "'");
            }
        }

        inline void save(const std::filesystem::path& file, const geo::Mesh& mesh, bool attributes = true, geo::MeshElementsFlags elements = geo::MeshElementsFlags::MESH_ALL_ELEMENTS)
        {
            geo::MeshIOFlags flags;
            flags.set_dimension(mesh.vertices.dimension());
            flags.set_attribute(attributes ? geo::MESH_ALL_ATTRIBUTES : geo::MESH_NO_ATTRIBUTES);
            flags.set_elements(elements);
            flags.set_verbose(true);
            geo::mesh_save(mesh, file.string(), flags);
        }

        /**
         * @brief Calculate epsilon based on mesh bbox diag and eps value (similar to tetwild), and add another boundary eps value (e-14)
         *
         * @param mesh
         * @param eps
         * @return double
         */
        inline double calulate_epsilon(const geo::Mesh& mesh, const double eps)
        {
            //return geo::bbox_diagonal(mesh) * eps + 2.0 * std::numeric_limits<double>::epsilon();
            return 0.0;
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
