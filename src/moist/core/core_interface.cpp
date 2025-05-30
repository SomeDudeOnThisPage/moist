#include "core_interface.hpp"

#include <geogram/mesh/mesh_io.h>

#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"

moist::Interface::Interface(const std::filesystem::path file, const AxisAlignedInterfacePlane plane) : _plane(std::make_shared<AxisAlignedInterfacePlane>(plane))
{
    _triangulation = std::make_shared<geogram::Mesh>(3);
    moist::utils::geo::load(file, *_triangulation, 3); // load with dim 2, but mesh has dim 3 for projection

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

    this->WriteMetadata();
}

moist::Interface::Interface(const std::filesystem::path file)
{
    _triangulation = std::make_shared<geogram::Mesh>(3);
    moist::utils::geo::load(file, *_triangulation, 3); // load with dim 2, but mesh has dim 3 for projection

    _plane = std::make_shared<AxisAlignedInterfacePlane>(this->ReadMetadata());

    // project
    for (const g_index v : this->_triangulation->vertices)
    {
        switch (_plane->axis)
        {
            case Axis::X:
                break;
            case Axis::Y:
                break;
            case Axis::Z:
                this->_triangulation->vertices.point(v).z = _plane->extent;
                break;
        }
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

moist::AxisAlignedInterfacePlane moist::Interface::ReadMetadata() const
{
    if (this->_triangulation->cells.nb() == 1)
    {
        geogram::Attribute<double> m_extent(this->_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_EXTENT);
        geogram::Attribute<double> m_envelope(this->_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_ENVELOPE);
        geogram::Attribute<char> m_axis(this->_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_AXIS);

        return
        {
            moist::AXIS_OPTION_ARGUMENT_MAP.at(std::string(1, m_axis[0])),
            m_extent[0],
            m_envelope[0]
        };
    }
    else
    {
        OOC_ERROR("failed to deserialize metadata - expected 1 cell, got " << this->_triangulation->cells.nb());
        return AxisAlignedInterfacePlane {}; // TODO: exception
    }
}

void moist::Interface::WriteMetadata()
{
    if (_triangulation->cells.nb() == 0)
    {
        _triangulation->cells.create_tet(0, 0, 0, 0);
    }

    geogram::Attribute<double> m_extent(_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_EXTENT);
    geogram::Attribute<double> m_envelope(_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_ENVELOPE);
    geogram::Attribute<char> m_axis(_triangulation->cells.attributes(), ATTRIBUTE_META_PLANE_AXIS);

    m_extent[0] = _plane->extent;
    m_envelope[0] = _plane->epsilon;
    m_axis[0] = 'Z'; // TODO?
}
