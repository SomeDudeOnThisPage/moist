#ifndef __SCALAR_FIELD_CUH
#define __SCALAR_FIELD_CUH

#include <cstdint>
#include <iostream>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "scalar_field.hpp"
#include "utils.cuh"

namespace ooc
{
    // https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
    #define cuda_error(ans) { cuda_assert((ans), __FILE__, __LINE__); }
    static inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
    {
       if (code != cudaSuccess)
       {
          std::cout << (stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line) << std::endl;
          if (abort) exit(code);
       }
    }

    /**
     * Unoptimized Kernel counting sign changes (0, 1) in a 3d-Array. Each thread counts sign changes in a cube around /radius/ of one pixel.
     * TODO: Set a min-size of some cuboid, so as to not have to store each pixel, but only an average of areas to save memory.
     * TODO: Shared memory 👀. Currently, there's loads of global memory accesses, as each thread currently accesses radius^3 values!
     */
    __global__ static void averageScalarField(const PitchedMatrix3d input, PitchedMatrix3d output, int size_x, int size_y, int size_z, int radius)
    {
        const int local_x = cuda_utils::tidx();
        const int local_y = cuda_utils::tidy();
        const int local_z = cuda_utils::tidz();

        const bool is_debug = local_x == 1 && local_y == 1 && local_z == 0;

        if (is_debug)
        {
            printf("thread %d %d %d\n", local_x, local_y, local_z);
        }

        int sign_changes = 0;
        int visited = 0;
        // TODO: min with max. width/height/depth
        uint16_t last = cuda_utils::pitched::get<uint16_t>(input, max(-radius + local_x, 0), max(-radius + local_y, 0), max(-radius + local_z, 0));

        for (int dz = -radius + local_z; dz <= radius + local_x; dz++)
        {
            for (int dy = -radius + local_y; dy <= radius + local_y; dy++)
            {
                for (int dx = -radius + local_x; dx <= radius + local_z; dx++)
                {
                    if (dx < 0 || dy < 0 || dz < 0)
                    {
                        continue;
                    }

                    if (dx >= size_x || dy >= size_y || dz >= size_z)
                    {
                        continue;
                    }

                    const auto local_data = cuda_utils::pitched::get<uint16_t>(input, dx, dy, dz);
                    if (is_debug)
                    {
                        printf("[1, 1, 0] processing %d %d %d - last sign: %d, local: %d\n", dx, dy, dz, last, local_data);
                    }

                    if (last != local_data)
                    {
                        sign_changes++;
                        last = local_data;
                    }

                    visited++;
                }
            }
        }

        printf("thread %d %d %d counted %d sign changes, visited %d\n", cuda_utils::tidx(), cuda_utils::tidy(), cuda_utils::tidz(), sign_changes, visited);
    }
}

void ooc::generate_scalar_field(std::shared_ptr<TiffData> data)
{
    printf("running...\n");
    PitchedMatrix3d matrix;

    cudaExtent extent = make_cudaExtent(data->width() * sizeof(uint16_t), data->height(), 1);
    cuda_error(cudaMalloc3D(&matrix.ptr, extent));

    uint16_t* host_data = new uint16_t[data->width() * data->height() * 1];
    for (int z = 0; z < 1; z++)
    {
        for (int y = 0; y < data->height(); y++)
        {
            for (int x = 0; x < data->width(); x++)
            {
                int index = z * data->height() * data->width() + y * data->width() + x;
                host_data[index] = data->_data[z][y][x];
            }
        }
    }

    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr(host_data, data->width() * sizeof(uint16_t), data->width(), data->height());
    copyParams.dstPtr = matrix.ptr;
    copyParams.extent = extent;
    copyParams.kind = cudaMemcpyHostToDevice;

    cuda_error(cudaMemcpy3D(&copyParams));

    dim3 dim_block(2, 2, 1);
    dim3 dim_grid(1, 1, 1);
    printf("call...\n");
    ooc::averageScalarField<<<dim_grid, dim_block>>>(matrix, /* temp until output*/ matrix, data->width(), data->height(), 1, 2);
    cuda_error(cudaPeekAtLastError());
    printf("done...\n");
    cuda_error(cudaDeviceSynchronize());
}
#endif // __SCALAR_FIELD_CUH
