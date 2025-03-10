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
