cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()
add_executable(MOIST_Tests
    merger_test.gtest.cpp
     $<TARGET_OBJECTS:MOIST_Merger_Units>
)

target_link_libraries(MOIST_Tests PRIVATE gtest_main)
target_link_libraries(MOIST_Tests PRIVATE geogram)
target_link_libraries(MOIST_Tests PRIVATE OutOfCore_Meshing_Core)

include(GoogleTest)
gtest_discover_tests(MOIST_Tests)
