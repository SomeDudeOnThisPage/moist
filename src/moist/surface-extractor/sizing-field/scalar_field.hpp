#ifndef __SCALAR_FIELD_HPP
#define __SCALAR_FIELD_HPP

#include <memory>

#include <cuda_runtime.h>

#include "moist/core/timer.hpp"

#include "../tiff_data.hpp"

namespace moist
{
    struct IsoValue
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
        __host__ __device__ operator float() const { return isovalue; } // TODO: Does a static_cast use implicit conversion operator?
    };

    struct PitchedMatrix3d
    {
        size_t width;
        size_t height;
        size_t depth;
        size_t pitch;
        uint32_t* data;
    };

    void generate_scalar_field(const Tiff& tiff, const IsoValue isovalue, moist::metrics::Metrics_ptr metrics = nullptr);
}

#endif // __SCALAR_FIELD_HPP
