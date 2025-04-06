// void mark_vertices() // mark all edges to be decimated // TODO: maybe mark faces?`that way mark cells?
// void collapse_edges() // collapse all edges (or cells along an edge?) // TODO: Make sure to collapse ONTO an interface vertex if possible
// how to effectively handle maps? Use google's stuff??
// TODO: we have no grid to cluster onto, but edge connectivity to interface vertices...

// for each non-interface vertex on p:
// compute interface facet, on which point lies
// merge point with closest interface facet vertex
// delete all 0 volume tets
// this is pretty easily parallelizable

#ifndef __OOC_VERTEX_CLUSTERING
#define __OOC_VERTEX_CLUSTERING

#include "../core.hpp"

namespace incremental_meshing
{

}

#endif // __OOC_VERTEX_CLUSTERING
