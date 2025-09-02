#ifndef MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_
#define MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_

#include <array>
#include <string>
#include <vector>

namespace moist::testing::fixture
{
    using InterfaceInserterFiles = std::array<std::string, 3>;

    std::vector<InterfaceInserterFiles> CORE_GROUP_INTERFACE_INSERTION = {
        {"Cube/Cylinder","simple/cube_cylinder/0.mesh","simple/cube_cylinder/1.mesh"},
        {"Cube/Cylinder_2","cube_cylinder/interface.geogram","cube_cylinder/raw/cube.mesh"},
        {"Tight Carbonate","tight_carbonate/interface.geogram","tight_carbonate/raw/0-4.msh"},
        {"Packed Bed (res5)","target/interfaces/res5/425.geogram","target/volumes/res5/416-425.msh"},
        {"Packed Bed (res23)","res23-416_465-465_514/interface.geogram","res23-416_465-465_514/raw/416-465.msh"}
    };
}

#endif // MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_
