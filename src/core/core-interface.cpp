#include <geogram/mesh/mesh_io.h>

#include "core-interface.hpp"

incremental_meshing::Interface::Interface(const std::filesystem::path mesh_path, const AxisAlignedInterfacePlane plane) : _plane(std::make_shared<AxisAlignedInterfacePlane>(plane))
{
    geogram::MeshIOFlags flags;
    flags.set_elements(geogram::MeshElementsFlags(
        geogram::MeshElementsFlags::MESH_VERTICES |
        geogram::MeshElementsFlags::MESH_EDGES
    ));
    flags.set_dimension(2);
    flags.set_verbose(false);

    this->_triangulation = std::make_shared<geogram::Mesh>();
    if (!geogram::mesh_load(mesh_path, *this->_triangulation))
    {
        OOC_ERROR("Failed to load mesh: " << mesh_path);
        return;
    }
}

std::shared_ptr<geogram::Mesh> incremental_meshing::Interface::Triangulation()
{
    return this->_triangulation;
}

std::shared_ptr<incremental_meshing::AxisAlignedInterfacePlane> incremental_meshing::Interface::Plane()
{
    return this->_plane;
}
