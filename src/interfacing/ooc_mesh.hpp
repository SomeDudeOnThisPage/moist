#ifndef __OOC_MESH_HPP
#define __OOC_MESH_HPP

#include <unordered_set>

#include <geogram/mesh/mesh.h>

#include "../core.hpp"
#include "interface.hpp"

namespace ooc
{
    class OOCMesh : public geogram::Mesh
    {
    public:
        OOCMesh(std::string identifier, geogram::index_t dimension = 3, bool single_precision = false);
        void InsertInterface(ooc::Interface& interface);

    private:
        std::string _identifier;
        std::shared_ptr<ooc::Interface> _interface;
        std::unordered_set<geogram::index_t> _deleted_tets;

        void InsertVertex(geogram::vec3& point);
    };
}

#endif // __OOC_MESH_HPP
