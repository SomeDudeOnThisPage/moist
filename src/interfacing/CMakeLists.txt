cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(OutOfCore_Meshing_Interface
    DESCRIPTION "coplanar tetrahedral facet merger - interface plane extraction and retetrahedralization."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

add_executable(OutOfCore_Meshing_Interface ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_Interface PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_Interface PRIVATE CLI11::CLI11)
