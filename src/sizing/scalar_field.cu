#ifndef __SCALAR_FIELD_CUH
#define __SCALAR_FIELD_CUH

#include <cstdint>
#include <iostream>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "scalar_field.hpp"

namespace ooc
{
    #define cuda_error(ans) { cuda_assert((ans), __FILE__, __LINE__); }
    static inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
    {
       if (code != cudaSuccess)
       {
          std::cout << (stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line) << std::endl;
          if (abort) exit(code);
       }
    }

    struct PitchedMatrix3d
    {
        int width;
        int height;
        int depth;
        cudaPitchedPtr ptr;
    };

    static __device__ unsigned int tidx()
    {
        return blockIdx.x * blockDim.x + threadIdx.x;
    }

    static __device__ __inline__ unsigned int tidy()
    {
        return blockIdx.y * blockDim.y + threadIdx.y;
    }

    static __device__ __inline__ unsigned int tidz()
    {
        return blockIdx.y * blockDim.y + threadIdx.y;
    }

    namespace pitched
    {
        template <typename T>
        static __device__ T getValue(const PitchedMatrix3d& matrix, unsigned int x, unsigned int y, unsigned int z)
        {
            if (x >= matrix.width || y >= matrix.height || z >= matrix.depth)
            {
                return 0.0f;
            }

            char *ptr = (char*) matrix.ptr.ptr + z * matrix.ptr.ysize + y * matrix.ptr.pitch;
            T *row = (T*) ptr;
            return row[x];
        }

        /**
         * See https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1g32bd7a39135594788a542ae72217775c
         * for information on what this actually does.
         *
         * Also this helpful stackoverflow answer: https://stackoverflow.com/a/16119944/13728996 :)
         */
        /*template <typename T>
        __device__ T get(const T* elements, unsigned int pitch, unsigned int row, unsigned int column, unsigned int depth)
        {
            return ((T*) ((uint8_t*) elements + row * column * pitch))[depth];
        }

        template <typename T>
        __device__ void put(const T* elements, unsigned int pitch, unsigned int row, unsigned int column, unsigned int depth, T data)
        {
            ((T*) ((uint8_t*) elements + row * column * pitch))[depth] = data;
        }*/
    }

    /**
     * Unoptimized Kernel counting sign changes (0, 1) in a 3d-Array. Each thread counts sign changes in a cube around /radius/ of one pixel.
     * TODO: Set a min-size of some cuboid, so as to not have to store each pixel, but only an average of areas to save memory.
     * TODO: Shared memory 👀. Currently, there's loads of global memory accesses, as each thread currently accesses radius^3 values!
     */
    __global__ static void averageScalarField(/*const PitchedMatrix3d input, PitchedMatrix3d output,*/ /*int size_x, int size_y, int size_z, int radius*/)
    {
        printf("Hello from the Kernel!");
        /*const int local_x = blockIdx.x * blockDim.x + threadIdx.x;
        const int local_y = blockIdx.y * blockDim.y + threadIdx.y;
        const int local_z = blockIdx.z * blockDim.z + threadIdx.z;

        // Clamp.
        if (local_x >= size_x || local_y >= size_y || local_z >= size_z)
        {
            return;
        }

        int sign_changes = 0;
        int visited = 0;
        bool last = 0;

        for (int dz = -radius; dz <= radius; dz++)
        {
            for (int dy = -radius; dy <= radius; dy++)
            {
                for (int dx = -radius; dx <= radius; dx++)
                {
                    if (dx < 0 || dy < 0 || dz < 0)
                    {
                        continue;
                    }

                    if (dx >= size_x || dy >= size_y || dz >= size_z)
                    {
                        continue;
                    }

                    // const auto local_data = cuda_utils::pitched::get<int>(input.elements, input.pitch, cuda_utils::tidx() + dx, cuda_utils::tidy() + dy, cuda_utils::tidz() + dz);

                    const auto local_data = pitched::getValue<int>(input, dx, dy, dz);
                    if (last != local_data)
                    {
                        sign_changes++;
                        last = local_data;
                    }

                    visited++;
                }
            }
        }

        printf("thread %d %d %d counted %d sign changed\n", tidx(), tidy(), tidz(), sign_changes);*/
    }
}

void ooc::do_it()
{
    printf("running...\n");
    PitchedMatrix3d matrix;
    PitchedMatrix3d output;
    //cuda_error(cudaMalloc3D(&matrix.ptr, make_cudaExtent(10 /*width*/ * sizeof(int), 10 /*height*/, 10 /*depth*/)));
    //cuda_error(cudaMalloc3D(&output.ptr, make_cudaExtent(10 /*width*/ * sizeof(int), 10 /*height*/, 10 /*depth*/)));

    dim3 dim_block(32, 32);
    dim3 dim_grid(1, 1);
    printf("call...\n");
    ooc::averageScalarField<<<dim_grid, dim_block>>>();
    cuda_error(cudaPeekAtLastError());
    printf("done...\n");
    cuda_error(cudaDeviceSynchronize());
}
#endif // __SCALAR_FIELD_CUH
