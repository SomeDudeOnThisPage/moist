cmake_minimum_required(VERSION 3.20..3.31.5)

include(FetchContent)

set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/extern)

#######################################
### GEOGRAM
#######################################
set(GEOGRAM_LIB_ONLY ON CACHE BOOL "" FORCE)

set(GEOGRAM_WITH_GRAPHICS ON CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_LEGACY_NUMERICS OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_HLBFGS OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_TETGEN OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_LUA OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_FPG OFF CACHE BOOL "" FORCE)
set(GEOGRAM_USE_SYSTEM_GLFW3 ON CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_GARGANTUA OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_TBB OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_TRIANGLE ON CACHE BOOL "" FORCE) # uses triangle for CDT

FetchContent_Declare(
    geogram
    GIT_REPOSITORY https://github.com/BrunoLevy/geogram.git
    GIT_TAG main
)

FetchContent_MakeAvailable(geogram)
if(NOT geogram_POPULATED)
    FetchContent_Populate(geogram)
    add_subdirectory(
      ${geogram_SOURCE_DIR}
      ${geogram_BINARY_DIR}
      EXCLUDE_FROM_ALL
    )
endif()

#######################################
### CLI11
#######################################
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG main
)

FetchContent_MakeAvailable(cli11)
if(NOT cli11_POPULATED)
    FetchContent_Populate(cli11)
    add_subdirectory(
      ${cli11_SOURCE_DIR}
      ${cli11_BINARY_DIR}
      EXCLUDE_FROM_ALL
    )
endif()

#######################################
### fTetWild
#######################################
#if(OOC_WITH_TETWILD)
#    FetchContent_Declare(
#        fTetWild
#        GIT_REPOSITORY https://github.com/wildmeshing/fTetWild.git
#        GIT_TAG        master
#    )

#    FetchContent_MakeAvailable(fTetWild) # goofy aaah since this also downloads and builds it's own geogram...
#endif()

# TODO: Make this better with a fork of MC33 where I fix the problems it has with modern cpp and cmake...
#######################################
### MC33
#######################################
#FetchContent_Declare(
#    mc33
#    GIT_REPOSITORY https://github.com/dvega68/MC33_cpp_library.git
#    GIT_TAG main
#)
#
#FetchContent_MakeAvailable(mc33)
#
#add_library(MC33 STATIC
#    ${mc33_SOURCE_DIR}/include/MC33.h
#    ${mc33_SOURCE_DIR}/source/MC33.cpp
#    ${mc33_SOURCE_DIR}/source/libMC33++.cpp
#    ${mc33_SOURCE_DIR}/source/MC33_LookUpTable.h
#    ${mc33_SOURCE_DIR}/source/surface.cpp
#    ${mc33_SOURCE_DIR}/source/grid3d.cpp
#)
#
## TODO: Probably best to compile a different static version for each of these types, and then choose based on parameter...
#if(TIFF_DATA_FORMAT STREQUAL "u8")
#    message(STATUS "compiling mc33 with u8 data type")
#    target_compile_definitions(MC33 PRIVATE integer_GRD)
#    target_compile_definitions(MC33 PRIVATE GRD_type_size=1)
#elseif(TIFF_DATA_FORMAT STREQUAL "u16")
#    message(STATUS "compiling mc33 with u16 data type")
#    target_compile_definitions(MC33 PRIVATE integer_GRD)
#    target_compile_definitions(MC33 PRIVATE GRD_type_size=2)
#elseif(TIFF_DATA_FORMAT STREQUAL "u32")
#    message(STATUS "compiling mc33 with u32 data type")
#    target_compile_definitions(MC33 PRIVATE integer_GRD)
#    target_compile_definitions(MC33 PRIVATE GRD_type_size=4)
#elseif(TIFF_DATA_FORMAT STREQUAL "f32")
#    message(STATUS "compiling mc33 with f32 data type")
#    target_compile_definitions(MC33 PRIVATE GRD_type_size=4)
#elseif(TIFF_DATA_FORMAT STREQUAL "f64")
#    message(STATUS "compiling mc33 with f64 data type")
#    target_compile_definitions(MC33 PRIVATE GRD_type_size=8)
#else()
#    message(FATAL_ERROR "Invalid tiff data type: ${TIFF_DATA_FORMAT}")
#endif()
#
## removed Wextra from original compiler flags in mc33 makefile
#target_compile_options(MC33 PRIVATE -std=c++11 -Ofast -m64 -Wall -funroll-loops)
#target_include_directories(MC33 PUBLIC ${mc33_SOURCE_DIR}/include)

# this doesnt fckn work. just make cgal a requirement to be installed...

#######################################
### CGAL
#######################################
FetchContent_Declare(
  CGAL
  GIT_REPOSITORY https://github.com/CGAL/cgal.git
  GIT_TAG        v5.6.1  # cgal 6 not on github currently, but we use deprecated 5.x features anyway...
)

FetchContent_MakeAvailable(CGAL)

#if(NOT CGAL_POPULATED)
#    FetchContent_Populate(CGAL)
#    add_subdirectory(
#      ${CGAL_SOURCE_DIR}
#      ${CGAL_BINARY_DIR}
#      EXCLUDE_FROM_ALL
#    )
#endif()

#######################################
### TinyTIFF
#######################################
#FetchContent_Declare(
#    TinyTIFF
#    GIT_REPOSITORY https://github.com/jkriege2/TinyTIFF
#    GIT_TAG master
#)

#FetchContent_MakeAvailable(TinyTIFF)
#if(NOT TinyTIFF_POPULATED)
#    FetchContent_Populate(TinyTIFF)
#    add_subdirectory(
#      ${TinyTIFF_SOURCE_DIR}
#      ${TinyTIFF_BINARY_DIR}
#      EXCLUDE_FROM_ALL
#    )
#endif()


#######################################
### EIGEN (why do I build this myself??, part of wmtk)
#######################################
# FetchContent_Declare(
#     eigen
#     GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#     GIT_TAG master
# )

# FetchContent_MakeAvailable(eigen)
# if(NOT eigen_POPULATED)
#     FetchContent_Populate(eigen)
#     add_subdirectory(
#       ${eigen_SOURCE_DIR}
#       ${eigen_BINARY_DIR}
#       EXCLUDE_FROM_ALL
#     )
# endif()
