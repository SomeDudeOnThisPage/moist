cmake_minimum_required(VERSION 3.20..3.31.5)

option(CMAKE_OPTION_UNROLL_LOOPS "Unroll simple loops." ON)
option(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS "Process local operations in parallel if applicable." ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
    OutOfCore_Meshing_InterfaceInserter
    DESCRIPTION "coplanar tetrahedral facet merger - interface inserter."
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(OOC_WITH_TESTS)
    add_library(MOIST_InterfaceInserter_Units OBJECT
        slice.hpp
        slice.cpp
        local_operations.hpp
        local_operations.cpp
    )
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE geogram)
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE OutOfCore_Meshing_Core)
    set(SOURCES
        main.cpp
    )
    set(HEADERS
        interface_inserter.hpp
    )
else()
    file(GLOB_RECURSE SOURCES "*.cpp")
    file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")
endif()

if (OOC_WITH_TESTS)
    add_executable(OutOfCore_Meshing_InterfaceInserter ${SOURCES} ${HEADERS} $<TARGET_OBJECTS:MOIST_InterfaceInserter_Units>)
else()
    add_executable(OutOfCore_Meshing_InterfaceInserter ${SOURCES} ${HEADERS})
endif()

target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE CLI11::CLI11)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE nlohmann_json::nlohmann_json)
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
