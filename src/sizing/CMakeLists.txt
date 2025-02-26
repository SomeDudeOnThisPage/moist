cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# CUDA standard?

project(OutOfCore_Meshing_SizingFieldGenerator
    DESCRIPTION "coplanar tetrahedral facet merger - sizing mesh generator."
    LANGUAGES CXX CUDA
)

include_directories("/usr/local/cuda/include")

file(GLOB_RECURSE SOURCES "*.cpp" "*.cu")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl" "*.cuh")

set_target_properties(OutOfCore_Meshing_SizingFieldGenerator PROPERTIES CUDA_ARCHITECTURES "75")
