cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_MeshMerger
    DESCRIPTION "coplanar tetrahedral facet merger - interface mesh generator."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

add_executable(OutOfCore_Meshing_MeshMerger ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE CLI11::CLI11)

target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE OutOfCore_Meshing_Core)
target_include_directories(OutOfCore_Meshing_MeshMerger PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(OutOfCore_Meshing_MeshMerger PROPERTIES OUTPUT_NAME "moist-merge")
