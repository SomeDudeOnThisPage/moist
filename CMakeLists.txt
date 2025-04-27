cmake_minimum_required(VERSION 3.20..3.31.5)

project(OutOfCore_Meshing
    DESCRIPTION "out of core surface extraction, tetrahedral meshing and sizing field generation tools"
)

# Options
option(OOC_WITH_TETWILD "Bundle fTetWild" ON)

set(TIFF_DATA_FORMAT "u16" CACHE STRING "Define data format: u8 u16 u32 f32 f64")
set_property(CACHE TIFF_DATA_FORMAT PROPERTY STRINGS u8 u16 u32 f32 f64)

add_subdirectory(extern)
add_subdirectory(src/core)

# generate as own "subprojects", that may be started to work in parallel on different "slices"
# memory-requirement per "subpoject" is the memory-footprint of ONE slice (the one being worked on), plus some constant
#     e.g., if a slice-mesh takes N memory, and consists of M tiff-slices, its a factor of N + (N/M * x), with x being some factor
# add_subdirectory(src/interfacing) # todo rename to "interface-inserter" or something
add_subdirectory(src/surface-extractor)
add_subdirectory(src/interface-generator)
add_subdirectory(src/interface-inserter)

if(OOC_WITH_TETWILD)
    #add_subdirectory(src/tetwild-wrapper)
endif()

# Make work again after core lib refactor
#set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
#add_subdirectory(src/sizing)
