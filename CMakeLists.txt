cmake_minimum_required(VERSION 3.20..3.31.5)

project(OutOfCore_Meshing
    DESCRIPTION "out of core surface extraction, tetrahedral meshing and sizing field generation tools"
)

add_subdirectory(extern)

# generate as own "subprojects", that may be started to work in parallel on different "slices"
add_subdirectory(src/interfacing)
# add_subdirectory(src/sizing)
