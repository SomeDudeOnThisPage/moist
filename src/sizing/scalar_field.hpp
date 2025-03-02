#ifndef __SCALAR_FIELD_HPP
#define __SCALAR_FIELD_HPP

#include <memory>

#include <cuda_runtime.h>

#include "tiff_data.hpp"

namespace ooc
{
    struct PitchedMatrix3d
    {
        int width;
        int height;
        int depth;
        cudaPitchedPtr ptr;
    };

    void generate_scalar_field(std::shared_ptr<TiffData> data);
}

#endif // __SCALAR_FIELD_HPP
