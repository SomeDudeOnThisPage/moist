cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_MeshMerger
    DESCRIPTION "coplanar tetrahedral facet merger - interface mesh generator."
    LANGUAGES CXX
)

if(OOC_WITH_TESTS)
    add_library(MOIST_Merger_Units OBJECT
        submesh_merger.hpp
        submesh_merger.cpp
    )
    target_link_libraries(MOIST_Merger_Units PRIVATE geogram)
    target_link_libraries(MOIST_Merger_Units PRIVATE OutOfCore_Meshing_Core)
    set(SOURCES
        main.cpp
    )
    set(HEADERS)
else()
    file(GLOB_RECURSE SOURCES "*.cpp")
    file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")
endif()



list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (OOC_WITH_TESTS)
    add_executable(OutOfCore_Meshing_MeshMerger ${SOURCES} ${HEADERS} $<TARGET_OBJECTS:MOIST_Merger_Units>)
else()
    add_executable(OutOfCore_Meshing_MeshMerger ${SOURCES} ${HEADERS})
endif()

target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE CLI11::CLI11)

target_link_libraries(OutOfCore_Meshing_MeshMerger PRIVATE OutOfCore_Meshing_Core)
target_include_directories(OutOfCore_Meshing_MeshMerger PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(OutOfCore_Meshing_MeshMerger PROPERTIES OUTPUT_NAME "moist-merge")
