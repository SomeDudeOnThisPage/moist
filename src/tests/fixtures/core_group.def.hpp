#ifndef MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_
#define MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_

#include <array>
#include <string>
#include <vector>

namespace moist::testing::fixture
{
    using InterfaceInserterFiles = std::array<std::string, 4>;

    std::vector<InterfaceInserterFiles> CORE_GROUP_INTERFACE_INSERTION = {
        {"Cube/Cylinder","cube_cylinder/interface.geogram","cube_cylinder/raw/cube.mesh","cube_cylinder/raw/cylinder.mesh"},
        {"Tight Carbonate","tight_carbonate/interface.geogram","tight_carbonate/raw/0-4.msh","tight_carbonate/raw/4-8.msh"},
        {"Packed Bed (res5)","res5-406_455-455-504/interface.geogram","res5-406_455-455-504/raw/406-455.msh","res5-406_455-455-504/raw/455-504.msh"},
        {"Packed Bed (res23)","res23-416_465-465_514/interface.geogram","res23-416_465-465_514/raw/416-465.msh","res23-416_465-465_514/raw/465-514.msh"}
    };
}

#endif // MOIST_TESTS_FIXTURES_CORE_GROUP_DEF_HPP_
