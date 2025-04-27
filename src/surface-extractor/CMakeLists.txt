cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_SurfaceExtractor
    DESCRIPTION "coplanar tetrahedral facet merger - surface extractor and sizing mesh generator."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.h" "*.hpp" "*.inl")

find_package(CGAL REQUIRED)
find_package(TIFF REQUIRED)

add_executable(OutOfCore_Meshing_SurfaceExtractor ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE CLI11::CLI11)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE CGAL::CGAL)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE TIFF::TIFF)
#target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE MC33)
target_link_libraries(OutOfCore_Meshing_SurfaceExtractor PRIVATE OutOfCore_Meshing_Core)
#target_include_directories(OutOfCore_Meshing_SurfaceExtractor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
