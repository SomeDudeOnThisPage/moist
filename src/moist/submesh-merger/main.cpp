#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

struct Arguments
{
    std::filesystem::path mesh;
    std::filesystem::path slice;
};

int main(int argc, char* argv[])
{
    std::cout << "Hi!" << std::endl;
}
