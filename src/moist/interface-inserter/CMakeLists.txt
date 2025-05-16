cmake_minimum_required(VERSION 3.20..3.31.5)

option(CMAKE_OPTION_UNROLL_LOOPS "Unroll simple loops." ON)
option(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS "Process local operations in parallel if applicable." ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_InterfaceGenerator
    DESCRIPTION "coplanar tetrahedral facet merger - interface inserter."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

add_executable(OutOfCore_Meshing_InterfaceInserter ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE CLI11::CLI11)

target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE OutOfCore_Meshing_Core)
target_include_directories(OutOfCore_Meshing_InterfaceInserter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_OPTION_UNROLL_LOOPS)
    message("Option: UNROLL_LOOPS")
    target_compile_definitions(OutOfCore_Meshing_InterfaceInserter PRIVATE OPTION_UNROLL_LOOPS)
endif()

if(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS)
    message("Option: PARALLEL_LOCAL_OPERATIONS")
    target_compile_definitions(OutOfCore_Meshing_InterfaceInserter PRIVATE OPTION_PARALLEL_LOCAL_OPERATIONS)
endif()

set_target_properties(OutOfCore_Meshing_InterfaceInserter PROPERTIES OUTPUT_NAME "moist-insert")
