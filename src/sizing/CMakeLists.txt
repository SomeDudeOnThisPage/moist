cmake_minimum_required(VERSION 3.20..3.31.5)

#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(OutOfCore_Meshing_SizingFieldGeneratorProject
    DESCRIPTION "coplanar tetrahedral facet merger - sizing mesh generator."
    LANGUAGES CXX CUDA
)

include_directories("/usr/local/cuda-12.8/include")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE CUDA_SOURCES "*.cuh" "*.cu")
file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)
set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)

add_executable(OutOfCore_Meshing_SizingFieldGenerator ${SOURCES} ${HEADERS} ${CUDA_SOURCES})
target_link_libraries(OutOfCore_Meshing_SizingFieldGenerator PRIVATE CLI11::CLI11)

# still support 50 but offline comp. will be removed soon...
set_target_properties(OutOfCore_Meshing_SizingFieldGenerator PROPERTIES CUDA_ARCHITECTURES "75;80")
set_target_properties(OutOfCore_Meshing_SizingFieldGenerator PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
