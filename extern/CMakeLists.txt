cmake_minimum_required(VERSION 3.20..3.31.5)

include(FetchContent)

set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/extern)

#######################################
### GEOGRAM
#######################################
set(GEOGRAM_LIB_ONLY ON CACHE BOOL "" FORCE)

set(GEOGRAM_WITH_GRAPHICS OFF CACHE BOOL "" FORCE)
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
    GIT_TAG v1.9.5
)

FetchContent_MakeAvailable(geogram)
#if(NOT geogram_POPULATED)
#    FetchContent_Populate(geogram)
#    add_subdirectory(
#      ${geogram_SOURCE_DIR}
#      ${geogram_BINARY_DIR}
#      EXCLUDE_FROM_ALL
#    )
#endif()

#######################################
### CLI11
#######################################
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)

FetchContent_MakeAvailable(cli11)
#if(NOT cli11_POPULATED)
#    FetchContent_Populate(cli11)
#    add_subdirectory(
#      ${cli11_SOURCE_DIR}
#      ${cli11_BINARY_DIR}
#      EXCLUDE_FROM_ALL
#    )
#endif()

#######################################
### fTetWild
#######################################
if(OOC_WITH_TETWILD)
#    FetchContent_Declare(
#        fTetWild
#        GIT_REPOSITORY https://github.com/wildmeshing/fTetWild.git
#        GIT_TAG        master
#    )
#    message("using ${CMAKE_SOURCE_DIR}/extern as third party folder for tetwild")
#    set(LIBIGL_USE_STATIC_LIBRARY ON CACHE BOOL "Force libigl static library" FORCE)
#    set(INPUT_THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/extern)
#    FetchContent_MakeAvailable(fTetWild) # goofy aaah since this also downloads and builds it's own geogram...
endif()

#######################################
### MC33
#######################################
FetchContent_Declare(
    mc33
    GIT_REPOSITORY https://github.com/SomeDudeOnThisPage/mc33-cmake.git
    GIT_TAG v5.2.0
)

set(MC33_ORTHOGONAL_GRID OFF)
set(MC33_NEGATE_NORMALS OFF)
set(MC33_DATA_FORMAT f32)

FetchContent_MakeAvailable(mc33)

if (MOIST_TESTS)
  FetchContent_Declare(
      gtest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG  v1.17.0
  )

  FetchContent_MakeAvailable(gtest)
endif()

FetchContent_Declare(
    libtiff
    GIT_REPOSITORY https://gitlab.com/libtiff/libtiff.git
    GIT_TAG v4.7.0
)

#FetchContent_MakeAvailable(libtiff)

#FetchContent_Declare(
#  json
#  GIT_REPOSITORY https://github.com/nlohmann/json.git
#  GIT_TAG v3.12.0
#)

#FetchContent_MakeAvailable(json)

set(CDT_USE_AS_COMPILED_LIBRARY ON CACHE BOOL "") # what a mess...
FetchContent_Declare(
  cdt
  GIT_REPOSITORY https://github.com/artem-ogre/CDT.git
  GIT_TAG 1.4.4
  SOURCE_SUBDIR CDT
)

FetchContent_MakeAvailable(cdt)

# Idea: Use cgals exact predicates??
#FetchContent_Declare(
#  boost
#  GIT_REPOSITORY https://github.com/boostorg/boost.git
#  GIT_TAG        boost-1.84.0
#)
#
#set(BOOST_ENABLE_CMAKE ON CACHE BOOL "")
#FetchContent_MakeAvailable(boost)

find_package(Boost REQUIRED)
#find_package(GMP REQUIRED)
#find_package(MPFR REQUIRED)
find_package(CGAL REQUIRED)
#FetchContent_Declare(
#  CGAL
#  GIT_REPOSITORY https://github.com/CGAL/cgal.git
#  GIT_TAG        v5.6.1
#)

set(CGAL_HEADER_ONLY OFF)
set(CGAL_USE_GMP ON)
set(CGAL_USE_MPFR ON)

#FetchContent_MakeAvailable(CGAL)

# TODO: Maybe create a "forked" github repo, port tetrahedron to cmake and c++23, and include via fetchcontent...
add_subdirectory(tetrahedron)
add_subdirectory(triangle)
add_subdirectory(tetgen)
