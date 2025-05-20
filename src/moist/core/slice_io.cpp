#include "slice_io.hpp"

void moist::slice_io::extract(const std::filesystem::path file, const moist::descriptor::MeshSlice descriptor, geogram::Mesh &slice)
{
}

// probably unecessarily complicated, why not just do this "linearly" at the end, and just iterate all slices that should exist then?
// that way I don't have to worry about "out of order" insertions of slices...
// and just iteratively flush slice after slice into the main file...
void moist::slice_io::insert(const std::filesystem::path file, const moist::descriptor::MeshSlice descriptor, geogram::Mesh &slice)
{
    // create a map vec3 -> [g_index (local slice), g_index (global mesh)], and populate the local field (can also be done as attribute)
    // check each "block" around the current slice, i.e. iterate the bigmesh descriptor, find all already inserted slices S, check their interfaces.
    // if an interface matches the given slices' interfaces, go through all the vertices of S, and insert the global index into the map if positions match
    // write all remaining vertices without a global index, begin at current global.num_vertices read from the header or metadata
    // when writing cells, begin at global.num_tets read from the header or metadata
    //      when an index is written that has a global
}
