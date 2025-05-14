#ifndef __SURFACE_GENERATOR_HPP
#define __SURFACE_GENERATOR_HPP

#include <MC33.h>
#include <geogram/mesh/mesh.h>

#include <core.hpp>

#include "tiff_data.hpp"

namespace moist
{
    class SurfaceGenerator
    {
    public:
        SurfaceGenerator(const TiffData& tiff, const uint32_t offset, const Axis axis = Axis::Z, const bool center = true, const bool invert = false);
        ~SurfaceGenerator() = default;

        void generate(geogram::Mesh& mesh, const float isovalue);
    private:
        MC33 _mc;
        surface _surface;
        grid3d _grid;

        uint32_t _width;
        uint32_t _height;
        uint32_t _depth;
        uint32_t _offset;
        Axis _axis;
        bool _invert;
        bool _center;
    };
}

#endif // __SURFACE_GENERATOR_HPP
