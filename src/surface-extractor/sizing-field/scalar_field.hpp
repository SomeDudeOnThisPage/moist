#ifndef __SCALAR_FIELD_HPP
#define __SCALAR_FIELD_HPP

#include <memory>

#include <cuda_runtime.h>

#include "../tiff_data.hpp"

namespace incremental_meshing
{
    struct IsoValue // TODO: Make this core and use in more places...
    {
        float isovalue;

        __host__ __device__ static float check(const float v)
        {
            if (v < 0.0f)
            {
                return 0.0f;
            }
            if (v > 1.0f)
            {
                return 1.0f;
            }
            return v;
        }

        __host__ __device__ IsoValue(const float v) : isovalue(check(v)){}
        __host__ __device__ float get() const { return isovalue; }
        __host__ __device__ operator float() const { return isovalue; }
    };

    struct PitchedMatrix3d
    {
        size_t width;
        size_t height;
        size_t depth;
        size_t pitch;
        uint32_t* data;
    };

    void generate_scalar_field(const TiffData& tiff, const IsoValue isovalue);
}

#endif // __SCALAR_FIELD_HPP
