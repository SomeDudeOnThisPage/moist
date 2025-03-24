cmake_minimum_required(VERSION 3.20..3.31.5)

option(CMAKE_OPTION_UNROLL_LOOPS "Unroll simple loops." ON)
option(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS "Process local operations in parallel if applicable." ON)
option(CMAKE_OPTION_INSERT_TEST_INTERFACE "Inserts a test-interface (test/interface.obj) into both SubMeshes." OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_Interface
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

if(CMAKE_OPTION_INSERT_TEST_INTERFACE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(OutOfCore_Meshing_Interface PRIVATE OPTION_DEBUG_TEST_INTERFACE)
endif()

if(CMAKE_OPTION_UNROLL_LOOPS)
    message("Option: UNROLL_LOOPS")
    target_compile_definitions(OutOfCore_Meshing_Interface PRIVATE OPTION_UNROLL_LOOPS)
endif()

if(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS)
    message("Option: PARALLEL_LOCAL_OPERATIONS")
    target_compile_definitions(OutOfCore_Meshing_Interface PRIVATE OPTION_PARALLEL_LOCAL_OPERATIONS)
endif()
