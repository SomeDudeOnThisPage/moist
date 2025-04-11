cmake_minimum_required(VERSION 3.20..3.31.5)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

add_library(OutOfCore_Meshing_Core STATIC ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_Core PRIVATE geogram)

target_include_directories(OutOfCore_Meshing_Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
