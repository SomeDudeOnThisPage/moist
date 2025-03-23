#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <sstream>
#include <string>

namespace incremental_meshing::utils {
    /**
     * Transforms string representations "(x,y,z)" into geogram vector.
     */
    geogram::vec3 parse_vector(std::string str) {
        double x, y, z;
        char none;

        std::istringstream stream(str);
        stream >> none >> x >> none >> y >> none >> z >> none;

        return geogram::vec3(x, y, z);
    }
}

#endif // __UTILS_HPP
