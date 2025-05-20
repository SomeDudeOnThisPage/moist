#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

struct Arguments
{
    std::filesystem::path mesh;
    std::filesystem::path slice;
};

// Load both Slice descriptors.
// Check the LocalInterface to be merged.
// Determine tets to be decimated by checking both InterfaceCellQuality lists.

// For each Slice:
// Load into memory
// Decimate required tets

int main(int argc, char* argv[])
{
    std::cout << "Hi!" << std::endl;
}
