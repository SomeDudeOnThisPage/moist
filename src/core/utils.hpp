#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <sstream>
#include <string>

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

#endif // __UTILS_HPP
