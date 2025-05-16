#ifndef MOIST_CORE_UTILS_HPP_
#define MOIST_CORE_UTILS_HPP_

#include <sstream>
#include <string>
#include <unordered_set>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

namespace moist::utils {
    /**
     * Transforms string representations "(x,y,z)" into geogram vector.
     */
    inline geogram::vec3 parse_vector(std::string str) {
        double x, y, z;
        char none;

        std::istringstream stream(str);
        stream >> none >> x >> none >> y >> none >> z >> none;

        return geogram::vec3(x, y, z);
    }

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

#ifndef NDEBUG
    inline void dump_mesh(geogram::Mesh& mesh, std::string file)
    {
        GEO::MeshIOFlags export_flags;
        export_flags.set_attribute(geogram::MESH_NO_ATTRIBUTES);
        export_flags.set_dimension(3);
        export_flags.set_elements(geogram::MeshElementsFlags::MESH_ALL_ELEMENTS);
        export_flags.set_verbose(true);
        geogram::mesh_save(mesh, file, export_flags);
    }
#endif
}

#endif // MOIST_CORE_UTILS_HPP_
