#include <iostream>

#include "moist/core/defines.hpp"
#include "moist/core/descriptor.hpp"

struct Arguments
{
    std::filesystem::path mesh;
    std::filesystem::path slice;
};

// Load interface tets of BigMesh (could be more than 1 interface!!!), load all tets of Slice
// "merge", i.e. create one continuous mesh
// decimate worst tet combinations along the interface plane
// save back into part of the file...
int main(int argc, char* argv[])
{
    // keep file for vertices and file for elements... like tetgen format...

    // slice...
    // find interface(s)...
    // slice to files

    // merge vertices and elements files into one .msh file
}
