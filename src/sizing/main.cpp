#include <iostream>

#include "scalar_field.hpp"

int main(int argc, char* argv[])
{
    std::cout << "TODO: Sizing Field Generation..." << std::endl;

    ooc::do_it();
    // 1. Load a "Slice" + Buffer data into memory
    // 2. Perform Scalar Field generation for this slice (e.g. counting sign changes)
    // 3. Generate Regular Grid (easy in beginning), later Octree
    // TWO OPTIONS: OPTION ONE:
    // 3.1. Place a Vertex in the center of each Cell with attribute "value" (TODO: was it "values"?)
    // 3.2. "value" = Average of cell contents with some factor based on total input pixel size to get realistic sizing values... lots of testing required...
    // 3.3. Delaunay-Tetrahedralize Vertices?
    // OPTION TWO:
    // 3.1. Create a "Cell" for each part of the octree, with 8 Vertices.
    // 3.2. Weigh "value" in vertices.
    // => This would need to Tetrahedralization of the Background Mesh => TODO: Can TetWild handle non-Manifold Sizing Meshes? It should, as each positon in the
    //    Mesh is locally checked in which "Cell" of the Sizing Mesh it lies, and linearly interpolated based on Distance.
    // 4. Save to .mesh file.

    // TODOS:
    //      1. How to make the sizing "configurable" to an extent?
    //      2. What about efficiency when being used with TetWild? Which "Option" leads to better runtime and, crucially, memory usage.
    //      3. How does the sizing field affect mesh quality? Especially with the interface merging.
    //      4. How to "chunk" the sizing mesh, we need only a "small" part of the overall mesh at one time, so we should also generate it in "slices".
    //         => e.g. some kind of "skirt" up and down, e.g. at least one cell above/below the boundaries

    // First implement Delaunay-Variant, as it's easier, later try to see if "Regular-Octree-Grid"-Approach leads to better results. Essentially Voronoi
    // vs. "Regular" sizing mesh value interpolation.
}
