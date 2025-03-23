#ifndef __SCALAR_FIELD_HPP
#define __SCALAR_FIELD_HPP

#include <memory>

#include <cuda_runtime.h>

#include "tiff_data.hpp"

namespace incremental_meshing
{
    struct PitchedMatrix3d
    {
        size_t width;
        size_t height;
        size_t depth;
        size_t pitch;
        uint32_t* data;
    };

    void generate_scalar_field(std::shared_ptr<TiffData> data);
}

#endif // __SCALAR_FIELD_HPP
