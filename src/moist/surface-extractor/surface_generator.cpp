#include "surface_generator.hpp"

moist::SurfaceGenerator::SurfaceGenerator(const Tiff &tiff, const uint32_t offset, const Axis axis, const bool center, const bool invert)
{
    _width = tiff.width();
    _height = tiff.height();
    _depth = tiff.depth();
    _offset = offset;
    _axis = axis;
    _invert = invert;
    _center = center;

    // create a "skirt" of 0-values around the entire mesh to close it off
    // _grid.set_grid_dimensions(tiff.width() + 2, tiff.height() + 2, tiff.depth() + 2);

    //const size_t sx = tiff.width();
    //const size_t sy = tiff.height();
    const size_t sz = tiff.depth();

    const size_t sx = 100;
    const size_t sy = 100;

    _grid.set_grid_dimensions(sx + 2, sy + 2, sz + 2);
    _mc.set_grid3d(_grid);

    for (size_t x = 1; x < sx + 1; x++)
    {
        for (size_t y = 1; y < sy + 1; y++)
        {
            for (size_t z = 1; z < sz + 1; z++)
            {
                const float datapoint = tiff.data[(z - 1) * tiff.width() * tiff.height() + (y - 1) * tiff.width() + x - 1];
                _grid.set_grid_value(x, y, z, (invert) ? 1.0 - datapoint : datapoint);
            }
        }
    }
}

void moist::SurfaceGenerator::Generate(GEO::Mesh &mesh, const float isovalue)
{
    if (mesh.vertices.dimension() != 3)
    {
        OOC_ERROR("invalid mesh dimension, expected 3");
    }

    if (isovalue > 1.0f || isovalue < 0.0f)
    {
        OOC_ERROR("invalid isosurface bounds");
    }

    _mc.calculate_isosurface(_surface, isovalue);

    GEO::vector<double> vertices(_surface.get_num_vertices() * 3);
    const float c_offset_x = (_center) ? static_cast<float>(_width)  / 2.0f : 0.0f;
    const float c_offset_y = (_center) ? static_cast<float>(_height) / 2.0f : 0.0f;
    const float c_offset_z = (_center) ? static_cast<float>(_depth)  / 2.0f : 0.0f;

    for(g_index v = 0; v < _surface.get_num_vertices(); v++)
    {
        const auto vertex = _surface.getVertex(v);
        // after all this time, why shouldn't I just move everything here? :)
        vertices[3 * v] = vertex[0] - c_offset_x + ((_axis == moist::Axis::X) ? _offset : 0);
        vertices[3 * v + 1] = vertex[1] - c_offset_y + ((_axis == moist::Axis::Y) ? _offset : 0);
        vertices[3 * v + 2] = vertex[2] - c_offset_z + ((_axis == moist::Axis::Z) ? _offset : 0);
    }

    GEO::vector<GEO::index_t> triangles(_surface.get_num_triangles() * 3);
    for(g_index t = 0; t < _surface.get_num_triangles(); t++)
    {
        const auto triangle = _surface.getTriangle(t);
        triangles[3 * t] = triangle[0];
        triangles[3 * t + 1] = triangle[1];
        triangles[3 * t + 2] = triangle[2];
    }

    mesh.facets.assign_triangle_mesh((GEO::coord_index_t) 3, vertices, triangles, true);
}
