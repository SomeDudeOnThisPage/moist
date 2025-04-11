cmake_minimum_required(VERSION 3.20..3.31.5)

project(OutOfCore_Meshing
    DESCRIPTION "out of core surface extraction, tetrahedral meshing and sizing field generation tools"
)

add_subdirectory(extern)

# shared stuff as static library
add_subdirectory(src/core)

# generate as own "subprojects", that may be started to work in parallel on different "slices"
# add_subdirectory(src/interfacing) # todo rename to "interface-inserter" or something
add_subdirectory(src/interface-generator)

# Make work again after core lib refactor
#set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
#add_subdirectory(src/sizing)
