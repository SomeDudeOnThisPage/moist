#include "interface_inserter.hpp"

#include "merger.hpp"

void moist::create_interface_mesh(geo::Mesh& a, geo::Mesh& b, const moist::AxisAlignedPlane& plane)
{
    moist::Merger merger(a, b, plane);
}
