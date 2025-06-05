cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(triangle
    DESCRIPTION " triangle, a C++ code which can compute properties of a triangle, including angles, area, centroid, circumcircle, edge lengths, incircle, orientation, orthocenter, and quality."
    LANGUAGES CXX
)

add_library(triangle SHARED
    source/triangle.cpp
)

add_library(triangle::triangle ALIAS triangle)

#target_compile_options(triangle PRIVATE -Wall -I)
target_compile_features(triangle PRIVATE cxx_std_23)
target_include_directories(triangle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)
