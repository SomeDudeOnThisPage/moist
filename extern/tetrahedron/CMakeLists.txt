cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(tetrahedron
    DESCRIPTION "tetrahedron, a C++ code which computes properties, including the centroid, circumsphere, dihedral angles, edge lengths, face angles, face areas, insphere, quality, solid angles, and volume, of a tetrahedron in 3D."
    LANGUAGES CXX
)

add_library(tetrahedron SHARED
    source/tetrahedron.cpp
)

add_library(tetrahedron::tetrahedron ALIAS tetrahedron)

target_compile_options(tetrahedron PRIVATE -Wall -I)
target_compile_features(tetrahedron PUBLIC cxx_std_23)
target_include_directories(tetrahedron PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)
