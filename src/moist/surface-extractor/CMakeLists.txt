cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CUDA_STANDARD 20)
#set(CMAKE_CUDA_STANDARD_REQUIRED ON)
#set(CMAKE_CUDA_COMPILER ${CUDA_DIR}/nvcc)
#include_directories(${CUDA_INCLUDE_DIR})

project(
    OutOfCore_Meshing_SurfaceExtractor
    DESCRIPTION "coplanar tetrahedral facet merger - surface extractor and sizing mesh generator."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE CUDA_SOURCES "*.cuh" "*.cu")
file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.h" "*.hpp" "*.inl")
find_package(TIFF REQUIRED)

#set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CXX)
#set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)

add_executable(OutOfCore_Meshing_SurfaceExtractor ${SOURCES} ${HEADERS}) #${CUDA_SOURCES}
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE CLI11::CLI11)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE TIFF::TIFF)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE MC33::MC33)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE OutOfCore_Meshing_Core)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE CGAL::CGAL)
#target_include_directories(OutOfCore_Meshing_SurfaceExtractor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#set_target_properties(OutOfCore_Meshing_SurfaceExtractor PROPERTIES CUDA_ARCHITECTURES "75;80")
#set_target_properties(OutOfCore_Meshing_SurfaceExtractor PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(OutOfCore_Meshing_SurfaceExtractor PROPERTIES OUTPUT_NAME "moist-extract")
