#include "core_interface.hpp"

#include <geogram/mesh/mesh_io.h>

#include "moist/core/attributes.inl"

moist::Interface::Interface(const std::filesystem::path mesh_path, const AxisAlignedInterfacePlane plane) : _plane(std::make_shared<AxisAlignedInterfacePlane>(plane))
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

    // project
    for (const g_index v : this->_triangulation->vertices)
    {
        switch (plane.axis)
        {
            case Axis::X:
                break;
            case Axis::Y:
                break;
            case Axis::Z:
                this->_triangulation->vertices.point(v).z = plane.extent;
                break;
        }
    }

    if (this->_triangulation->cells.nb() == 0)
    {
        // create a "dummy" cell for storing metadata attributes (the plane extent, i.e. the "Interface Index")
        // this is kind of a hack so I don't have to write an entire mesh saver / loader and can still retain metadata about this interface inside this file...
        this->_triangulation->cells.create_tet(0, 0, 0, 0);
        geogram::Attribute<double> m_plane(this->_triangulation->cells.attributes(), ATTRIBUTE_INTERFACE_INDEX_A);
        m_plane[0] = plane.extent;
    }
}

std::shared_ptr<geogram::Mesh> moist::Interface::Triangulation()
{
    return this->_triangulation;
}

std::shared_ptr<moist::AxisAlignedInterfacePlane> moist::Interface::Plane()
{
    return this->_plane;
}
