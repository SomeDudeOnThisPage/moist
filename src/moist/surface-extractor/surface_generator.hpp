#ifndef MOIST_SURFACE_EXTRACTOR_SURFACE_GENERATOR_HPP_
#define MOIST_SURFACE_EXTRACTOR_SURFACE_GENERATOR_HPP_

#include <MC33.h>
#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

#include "tiff_data.hpp"

namespace moist
{
    class SurfaceGenerator
    {
    public:
        SurfaceGenerator(const Tiff& tiff, const uint32_t offset, const Axis axis = Axis::Z, const bool center = true, const bool invert = false);
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

#endif // MOIST_SURFACE_EXTRACTOR_SURFACE_GENERATOR_HPP_
