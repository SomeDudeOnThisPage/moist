#ifndef __SCALAR_FIELD_CUH
#define __SCALAR_FIELD_CUH

#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cmath>

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <tiffio.h>

#include "scalar_field.hpp"
#include "utils.cuh"

#define DIV_UP(x, y) (x + y - 1) / y

// #define KERNEL_DEBUG(msg) printf("[%d,%d,%d]: %s\n", cuda_utils::tidx(), cuda_utils::tidy(), cuda_utils::tidz(), msg)

namespace moist
{
    // https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
    #define cuda_error(ans) { cuda_assert((ans), __FILE__, __LINE__); }
    static inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
    {
       if (code != cudaSuccess)
       {
          OOC_ERROR("cuda_error " << cudaGetErrorString(code));
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

struct CU_TiffParameters
{
    size_t size;
    size_t width;
    size_t height;
    size_t depth;
};

// TODO: Should I make this configurable or only use it for internal testing, settling on the (tested) best values?
enum CU_ScalarFieldOverflowWrappingMode : uint8_t
{
    OVERFLOW_MODE_ZERO = 0,             /** Set over/underflow to zero. */
    OVERFLOW_MODE_WRAP_AROUND = 1,      /** Wrap around the entire given data in the over/underflow direction. */
    OVERFLOW_MODE_WRAP_INTO_SELF = 2    /** Wrap around the current cell being processed - cells must be sized at least `2 * CU_ScalarFieldParameters#kernel_radius`. */
};

struct CU_ScalarFieldParameters
{
    /**
     * @brief Minimal size of an octree node in the scalar field in pixels (one thread will walk one 3d-block of extent `min_size * min_size * min_size`).
     *
     * This also defines the size of the output array.
     */
    size_t octree_node_size;
    size_t kernel_radius;
    CU_ScalarFieldOverflowWrappingMode overflow_mode;
    moist::IsoValue isovalue;
};

__global__ static void scalar_field_kernel(const float* __restrict__ input, uint16_t* __restrict__ output, const CU_TiffParameters tiff_parameters, const CU_ScalarFieldParameters scalar_field_parameters)
{
    constexpr uint16_t UINT16_MAX_VALUE = cuda_utils::math::c_pow(2, 16);

    // gates
    const int tidx = cuda_utils::tidx();
    const int tidy = cuda_utils::tidy();
    const int tidz = cuda_utils::tidz();

    if (tidx >= tiff_parameters.width || tidy >= tiff_parameters.height || tidz >= tiff_parameters.depth) return;

#ifndef NDEBUG
    const bool is_debug_thread = tidx == 100 && tidy == 100 && tidz == 10;
    if (is_debug_thread) { printf("DEBUG_CUDA: debug info from thread (%llu,%llu,%llu)\r\n", (unsigned long long) tidx, (unsigned long long) tidy, (unsigned long long) tidz); }
    if (is_debug_thread) { printf("DEBUG_CUDA: radius: %d, tiff-size: (%llu,%llu,%llu)\r\n", scalar_field_parameters.kernel_radius, (unsigned long long) tiff_parameters.width, (unsigned long long) tiff_parameters.height, (unsigned long long) tiff_parameters.depth); }
#endif

    const float pixel = input[cuda_utils::index(tidx, tidy, tidz, tiff_parameters.width, tiff_parameters.height)];

#define OPTION_GENERATE_VISUAL_SIZING_FIELD
#ifdef OPTION_GENERATE_VISUAL_SIZING_FIELD // ignore z changes so the sizing-field-visualization of a slice is unimpacted
    const int radius = scalar_field_parameters.kernel_radius;
    const int radius_z = 0;
#else
    const int radius = scalar_field_parameters.kernel_radius;
#endif //OPTION_GENERATE_VISUAL_SIZING_FIELD
    unsigned int count = 0;

    if (is_debug_thread) { printf("DEBUG_CUDA: tidz=%d, local_delta_z=%d, dz=%d\n", (int)tidz, -radius, ((int)tidz)+(-radius)); }
#ifdef OPTION_GENERATE_VISUAL_SIZING_FIELD
    for (int local_delta_z = 0; local_delta_z <= 0; local_delta_z++)
    {
        const int dz = tidz;
    #else
    for (int local_delta_z = -radius; local_delta_z <= radius; local_delta_z++)
    {
        const int dz = tidz + local_delta_z;
    #endif
        if (is_debug_thread) { printf("DEBUG_CUDA: %d\n", dz); }

        if (dz < 0 || dz >= tiff_parameters.depth) continue;

        for (int local_delta_y = -radius; local_delta_y <= radius; local_delta_y++) {
            const int dy = tidy + local_delta_y;
            if (dy < 0 || dy >= tiff_parameters.height) continue;

            for (int local_delta_x = -radius; local_delta_x <= radius; local_delta_x++)
            {
                const int dx = tidx + local_delta_x;
                if (dx < 0 || dx >= tiff_parameters.width) continue;

                const float pixel_other = input[cuda_utils::index(dx, dy, dz, tiff_parameters.width, tiff_parameters.height)];
            #ifndef NDEBUG
                if (is_debug_thread) { printf("DEBUG_CUDA: center (%d,%d,%d): %.2f, other (%d,%d,%d): %.2f\r\n", tidx, tidy, tidz, pixel, dx, dy, dz, pixel_other); }
            #endif
                if (pixel != pixel_other)
                {
                    count++;
                }
            }
        }
    }

    output[cuda_utils::index(tidx, tidy, tidz, tiff_parameters.width, tiff_parameters.height)] = count;
}

void moist::generate_scalar_field(const Tiff& tiff, const IsoValue isovalue, moist::metrics::Metrics_ptr metrics)
{
    moist::Timer timer("generate_scalar_field", metrics);
    const size_t size = tiff.width() * tiff.height() * tiff.depth();

    const CU_TiffParameters tiff_parameters
    {
        size,
        tiff.width(),
        tiff.height(),
        tiff.depth(),
    };

    const CU_ScalarFieldParameters scalar_field_parameters
    {
        5,
        25,
        CU_ScalarFieldOverflowWrappingMode::OVERFLOW_MODE_ZERO,
        isovalue
    };

    const dim3 dim_block(
        scalar_field_parameters.octree_node_size,
        scalar_field_parameters.octree_node_size,
        scalar_field_parameters.octree_node_size
    );

    const dim3 dim_grid(
        DIV_UP(tiff.width(), dim_block.x),
        DIV_UP(tiff.height(), dim_block.y),
        DIV_UP(tiff.depth(), dim_block.z)
    );

    float* d_input;
    uint16_t* d_output;
    std::vector<uint16_t> h_output(size);

    cuda_error(cudaMalloc(&d_input, size * sizeof(float)));
    cuda_error(cudaMalloc(&d_output, size * sizeof(uint16_t)));

    cuda_error(cudaMemcpy(d_input, tiff.data.data(), size * sizeof(float), cudaMemcpyHostToDevice));

    OOC_DEBUG("calling w/ block dimensions " << dim_block.x << ", " << dim_block.y << ", " << dim_block.z);
    OOC_DEBUG("calling w/ grid dimensions " << dim_grid.x << ", " << dim_grid.y << ", " << dim_grid.z);

    cudaEvent_t start, stop;
    cuda_error(cudaEventCreate(&start));
    cuda_error(cudaEventCreate(&stop));

    cuda_error(cudaEventRecord(start));
    scalar_field_kernel<<<dim_grid, dim_block>>>(d_input, d_output, tiff_parameters, scalar_field_parameters);
    cuda_error(cudaEventRecord(stop));
    cuda_error(cudaEventSynchronize(stop));

    cuda_error(cudaPeekAtLastError());

    float ms = 0;
    cuda_error(cudaEventElapsedTime(&ms, start, stop));
    OOC_DEBUG("scalar_field_kernel runtime: " << ms << "ms");

    cuda_error(cudaMemcpy(h_output.data(), d_output, size * sizeof(uint16_t), cudaMemcpyDeviceToHost));
    cuda_error(cudaFree(d_input));
    cuda_error(cudaFree(d_output));

    // normalize h_output in-place to span [0, 65535 (or rather max() of data type)] for easy use in the actual sizing mesh generation, and visual export, without needing a separate float vector.
    // TODO: Make this work for other data types...
    // adapted from: https://t4tutorials.com/min-max-normalization-of-data-in-data-mining/
    // v' = (v-min)/(max-min)(new_max - new_min)+new_min // new_min is 0 so it's irrelevant... => v' = (v-min)/(max-min)(new_max)
    auto [min_it, max_it] = std::minmax_element(h_output.begin(), h_output.end());
    const uint16_t min = *min_it;
    const uint16_t max = *max_it;

    if (max != min)
    {
        constexpr uint16_t uint16_max = std::numeric_limits<uint16_t>::max(); // TODO: Get max from tiff-data.
        std::transform(h_output.begin(), h_output.end(), h_output.begin(), [max, min, uint16_max](const uint16_t v) -> uint16_t
        {
            // cast internally to uint32_t to avoid overflow shenanigans
            return static_cast<uint16_t>(((static_cast<uint32_t>(v) - min) * uint16_max) / (max - min));
        });
    }

#ifndef NDEBUG
    TIFF* tif = TIFFOpen("test.tif", "w");
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, tiff.width());
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, tiff.height());
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, tiff.width() * sizeof(uint16_t)));

    for (size_t depth = 0; depth < 1; depth++)
    {
        for (size_t row = 0; row < tiff.height(); row++)
        {
            if (TIFFWriteScanline(tif,  &h_output[row * tiff.width()], row, 0) < 0)
            {
                std::cerr << "Error writing row " << row << "\n";
                TIFFClose(tif);
            }
        }
    }

    TIFFClose(tif);
#endif

    cuda_error(cudaEventDestroy(start));
    cuda_error(cudaEventDestroy(stop));

    /*incremental_meshing::averageScalarField<<<dim_grid, dim_block>>>(
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
    cuda_error(cudaDeviceSynchronize());*/

    /*cudaMemcpy3DParms copy_device_to_host = {0};
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

    TIFFClose(tif);*/

    //free(host_data);
    //cuda_error(cudaFree(d_input.ptr));
    //cuda_error(cudaFree(d_output.ptr));
}
#endif // __SCALAR_FIELD_CUH
