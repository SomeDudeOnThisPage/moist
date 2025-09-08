cmake_minimum_required(VERSION 3.20..3.31.5)

project(OutOfCore_Meshing
    DESCRIPTION "out of core surface extraction, tetrahedral meshing and sizing field generation tools"
)

# Options
# option(OOC_WITH_TETWILD "Bundle fTetWild" ON)
option(MOIST_TESTS                        "Build Tests" ON)
option(CMAKE_OPTION_USE_LOOKUP_GRID       "Use LookupGrid optimizations" ON)

#set(CUDA_DIR "/usr/local/cuda-12.8/bin/")
#set(CUDA_INCLUDE_DIR "/usr/local/cuda-12.8/include")

include_directories(/usr/include/x86_64-linux-gnu)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
list(INSERT CMAKE_PREFIX_PATH 0 "$ENV{HOME}/local/cgal")

add_subdirectory(extern)
add_subdirectory(src/moist/core)

# generate as own "subprojects", that may be started to work in parallel on different "slices"
# memory-requirement per "subpoject" is the memory-footprint of ONE slice (the one being worked on), plus some constant
#     e.g., if a slice-mesh takes N memory, and consists of M tiff-slices, its a factor of N + (N/M * x), with x being some factor
# add_subdirectory(src/interfacing) # todo rename to "interface-inserter" or something
add_subdirectory(src/moist/surface-extractor)
add_subdirectory(src/moist/interface-generator)
add_subdirectory(src/moist/interface-inserter)
add_subdirectory(src/moist/submesh-merger)
add_subdirectory(src/moist/boundary-marker)

if(MOIST_TESTS)
    add_subdirectory(src/tests)
endif()
