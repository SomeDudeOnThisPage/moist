cmake_minimum_required(VERSION 3.20..3.31.5)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.hpp" "*.inl")

add_library(OutOfCore_Meshing_Core STATIC ${SOURCES} ${HEADERS})
target_link_libraries(OutOfCore_Meshing_Core PRIVATE geogram)
target_link_libraries(OutOfCore_Meshing_Core PRIVATE nlohmann_json::nlohmann_json)

target_include_directories(OutOfCore_Meshing_Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../)
set_target_properties(OutOfCore_Meshing_Core PROPERTIES OUTPUT_NAME "moist-core")
