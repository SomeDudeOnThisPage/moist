cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TETGEN_CXX_FLAGS "-O3")
set(PREDICATES_CXX_FLAGS "-O0")

project(tetgen
    DESCRIPTION ""
    LANGUAGES CXX
)

add_library(tetgen SHARED
    source/tetgen.cxx
    source/predicates.cxx
    source/tetgen.h
)

add_library(tetgen::tetgen ALIAS tetgen)

target_compile_features(tetgen PRIVATE cxx_std_23)
target_include_directories(tetgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)

#target_compile_options(tetgen PRIVATE $<$<COMPILE_FILE:source/tetgen.cxx>:${TETGEN_CXX_FLAGS}>)
#target_compile_options(tetgen PRIVATE $<$<COMPILE_FILE:source/predicates.cxx>:${PREDICATES_CXX_FLAGS}>)

#target_compile_definitions(tetgen PRIVATE TETLIBRARY)

find_library(M_LIB m) # 100% honest here, I don't really know why this is required on linux systems. I scratched my head for 5 hours and just asked chatgpt to check the tetgen source code and see what im missing. now it works.
if(M_LIB)
    target_link_libraries(tetgen PRIVATE ${M_LIB})
endif()
