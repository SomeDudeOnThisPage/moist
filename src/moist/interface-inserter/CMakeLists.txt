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

find_library(GMPXX_LIB gmpxx REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_library(MPFR_LIB mpfr REQUIRED)
# honestly im grasping at straws here WHY IS THIS SH NOT WORKING
message(STATUS "GMPXX_LIB = ${GMPXX_LIB}")
message(STATUS "GMP_LIB = ${GMP_LIB}")
message(STATUS "MPFR_LIB = ${MPFR_LIB}")

if(MOIST_TESTS)
    add_library(MOIST_InterfaceInserter_Units OBJECT
        arrangement.hpp
        arrangement.cpp
        #slice.hpp
        #slice.cpp
        merger.hpp
        merger.cpp
        tetgen_utils.hpp
        tetgen_utils.cpp
        #slice_exact.cpp
        local_operations.hpp
        local_operations.cpp
        lookup_grid.hpp
        lookup_grid.cpp
        exact_mesh.hpp
        exact_mesh.cpp
        interface_inserter.hpp
        interface_inserter.cpp
        new_predicates.inl
        geometry_exact.inl
        exact_types.hpp
    )
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE geogram)
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE CLI11::CLI11)
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE OutOfCore_Meshing_Core)
    #target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE ${GMP_LIB} ${GMPXX_LIB} ${MPFR_LIB})
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE CGAL::CGAL)
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE tetgen::tetgen)
    target_link_libraries(MOIST_InterfaceInserter_Units PRIVATE mmg3d)
    target_include_directories(MOIST_InterfaceInserter_Units PRIVATE ${MMG3D_INCLUDE_DIR})

    set(SOURCES
        main.cpp
    )
    set(HEADERS
    )
else()
    file(GLOB_RECURSE SOURCES "*.cpp")
    file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")
endif()

if (MOIST_TESTS)
    add_executable(OutOfCore_Meshing_InterfaceInserter ${SOURCES} ${HEADERS} $<TARGET_OBJECTS:MOIST_InterfaceInserter_Units>)
else()
    add_executable(OutOfCore_Meshing_InterfaceInserter ${SOURCES} ${HEADERS})
endif()


target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE CLI11::CLI11)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE OutOfCore_Meshing_Core)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE ${GMPXX_LIB} ${GMP_LIB} ${MPFR_LIB})
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE CGAL::CGAL)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE tetgen::tetgen)
target_link_libraries(OutOfCore_Meshing_InterfaceInserter PRIVATE mmg3d)

target_include_directories(OutOfCore_Meshing_InterfaceInserter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(OutOfCore_Meshing_InterfaceInserter PRIVATE ${MMG3D_INCLUDE_DIR})

if(CMAKE_OPTION_UNROLL_LOOPS)
    message("Option: UNROLL_LOOPS")
    target_compile_definitions(OutOfCore_Meshing_InterfaceInserter PRIVATE OPTION_UNROLL_LOOPS)
endif()

if(CMAKE_OPTION_PARALLEL_LOCAL_OPERATIONS)
    message("Option: PARALLEL_LOCAL_OPERATIONS")
    target_compile_definitions(OutOfCore_Meshing_InterfaceInserter PRIVATE OPTION_PARALLEL_LOCAL_OPERATIONS)
endif()

if(CMAKE_OPTION_USE_LOOKUP_GRID)
    message("Option: USE_LOOKUP_GRID")
    target_compile_definitions(OutOfCore_Meshing_InterfaceInserter PRIVATE OPTION_LOOKUP_GRID)
endif()

set_target_properties(OutOfCore_Meshing_InterfaceInserter PROPERTIES OUTPUT_NAME "moist-insert")
