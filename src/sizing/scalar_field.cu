#ifndef __SCALAR_FIELD_CUH
#define __SCALAR_FIELD_CUH

#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <tiffio.h>

#include "scalar_field.hpp"
#include "utils.cuh"

#define DIV_UP(x, y) (x + y - 1) / y

// #define KERNEL_DEBUG(msg) printf("[%d,%d,%d]: %s\n", cuda_utils::tidx(), cuda_utils::tidy(), cuda_utils::tidz(), msg)

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
    __global__ static void averageScalarField(uint16_t* i, size_t i_pitch, uint16_t* o, size_t o_pitch, int size_x, int size_y, int size_z, int radius)
    {
        const int local_x = cuda_utils::tidx();
        const int local_y = cuda_utils::tidy();
        const int local_z = cuda_utils::tidz();

        const bool is_debug = local_x == 100 && local_y == 100 && local_z == 0;
        if (local_x >= size_x || local_y >= size_y || local_z >= size_z)
        {
            return;
        }

        if (is_debug)
        {
            printf("thread %d %d %d\n", local_x, local_y, local_z);
        }

        int sign_changes = 0;
        int visited = 0;
        // TODO: min with max. width/height/depth
        uint16_t last = i[max(-radius + local_y, 0) * i_pitch * sizeof(uint16_t) + max(-radius + local_x, 0)];  // cuda_utils::pitched::get<uint16_t>(&input, max(-radius + local_x, 0), max(-radius + local_y, 0), max(-radius + local_z, 0));

        for (int dz = max(0, -radius + local_z); dz < min(size_z, radius + local_z); dz++)
        {
            // TODO: better kernel error macro
            if (is_debug)
            {
                printf("[15, 15, 0] dx_from=%d dx_to=%d dy_from=%d dy_to=%d dz_from=%d dz_to=%d\n",
                    -radius + local_x, radius + local_x, -radius + local_y, radius + local_y, max(0, -radius + local_z), min(size_z, radius + local_z));
            }
            for (int dy = max(0, -radius + local_y); dy < min(size_y, radius + local_y); dy++)
            {
                for (int dx = max(0, -radius + local_x); dx < min(size_x, radius + local_x); dx++)
                {
                    const auto local_data = i[dy * i_pitch * sizeof(uint16_t) + dx]; // cuda_utils::pitched::get<uint16_t>(&input, dx, dy, dz);
                    if (is_debug)
                    {
                        printf("[15, 15, 0] processing %d %d %d - last sign: %d, local: %d\n", dx, dy, dz, last, local_data);
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

        if (is_debug)
        {
            printf("thread %d %d %d counted %d sign changes, visited %d\n", cuda_utils::tidx(), cuda_utils::tidy(), cuda_utils::tidz(), sign_changes, visited);
        }

        o[local_y * i_pitch * sizeof(uint16_t) + local_x] = sign_changes * 10000 / radius; //
        // cuda_utils::pitched::set<uint16_t>(&o, o_pitch, local_x, local_y, local_z, sign_changes * 100);
    }
}

void ooc::generate_scalar_field(std::shared_ptr<TiffData> data)
{
    cudaExtent extent = make_cudaExtent(data->width() * sizeof(uint16_t), data->height(), data->depth());

    uint16_t* host_data = new uint16_t[data->width() * data->height() * 1];
    for (int z = 0; z < data->depth(); z++)
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

    cudaPitchedPtr d_input;
    cudaPitchedPtr d_output;
    cuda_error(cudaMalloc3D(&d_input, extent));
    cuda_error(cudaMalloc3D(&d_output, extent));

    cudaMemcpy3DParms copy_host_to_device = {0};
    copy_host_to_device.srcPtr = make_cudaPitchedPtr(host_data, data->width() * sizeof(uint16_t), data->width(), data->height());
    copy_host_to_device.dstPtr = d_input;
    copy_host_to_device.extent = extent;
    copy_host_to_device.kind = cudaMemcpyHostToDevice;
    cuda_error(cudaMemcpy3D(&copy_host_to_device));

    dim3 dim_block(std::min((uint32_t) 16, data->width()), std::min((uint32_t) 16, data->height()), std::min((uint32_t) 16, data->depth()));
    dim3 dim_grid(DIV_UP(data->width(), 16), DIV_UP(data->height(), 16), DIV_UP(data->depth(), 16));

    OOC_DEBUG("calling w/ grid dimensions " << DIV_UP(data->width(), 16) << " " << DIV_UP(data->height(), 16) << " " << DIV_UP(data->depth(), 16));
    ooc::averageScalarField<<<dim_grid, dim_block>>>(
        (uint16_t*) d_input.ptr,
        d_input.pitch,
        (uint16_t*) d_output.ptr,
        d_output.pitch,
        data->width(),
        data->height(),
        data->depth(),
        5
    );
    cuda_error(cudaPeekAtLastError());
    cuda_error(cudaDeviceSynchronize());

    cudaMemcpy3DParms copy_device_to_host = {0};
    copy_device_to_host.srcPtr = d_output;
    copy_device_to_host.dstPtr = make_cudaPitchedPtr(host_data, data->width() * sizeof(uint16_t), data->width(), data->height());
    copy_device_to_host.extent = extent;
    copy_device_to_host.kind = cudaMemcpyDeviceToHost;
    cuda_error(cudaMemcpy3D(&copy_device_to_host));

    TIFF* tif = TIFFOpen("test.tif", "w");
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, data->width());
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, data->height());
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, data->width() * sizeof(uint16_t)));

    std::cout << copy_device_to_host.dstPtr.pitch << " " << data->width() << std::endl;
    for (uint32_t depth = 0; depth < data->depth(); depth++)
    {
        for (uint32_t row = 0; row < data->height(); row++)
        {
            if (TIFFWriteScanline(tif,  &host_data[row * data->width()], row, 0) < 0)
            {
                std::cerr << "Error writing row " << row << "\n";
                TIFFClose(tif);
            }
        }
    }

    TIFFClose(tif);
    free(host_data);
    cuda_error(cudaFree(d_input.ptr));
    cuda_error(cudaFree(d_output.ptr));
}
#endif // __SCALAR_FIELD_CUH
