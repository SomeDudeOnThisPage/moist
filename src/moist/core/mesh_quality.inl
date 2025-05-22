#ifndef MOIST_CORE_MESH_QUALITY_INL_
#define MOIST_CORE_MESH_QUALITY_INL_

#include <iostream>
#include <string>
#include <format>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include "moist/core/defines.hpp"
#include "moist/core/metrics.hpp"

namespace moist::mesh_quality
{
    /**
     * @brief Computes the distance from point p to a plane spanned by p0, a and b.
     *
     * @see https://stackoverflow.com/questions/3860206/signed-distance-between-plane-and-point
     *
     * @param p
     * @param p0
     * @param a
     * @param b
     * @return float
     */
    static float height(const vec3 p, const vec3 p0, const vec3 a, const vec3 b)
    {
        const vec3 normal = geogram::cross(p0 - a, p0 - b);
        return std::abs(geogram::dot(normal, p0 - p));
    }

    PURE INLINE float tetrahedron_aspect_ratio(const g_index cell, const geogram::Mesh& mesh)
    {
        // assert cell.nb_vertices() == 4;
        // assert cell.nb_edges() == 6;

        const geogram::index_t nb_vertices = mesh.cells.nb_vertices(cell);
        const geogram::index_t nb_edges = mesh.cells.nb_edges(cell);

        float heights[nb_vertices];
        for (l_index lv = 0; lv < nb_vertices; lv++)
        {
            heights[lv] = height(
                mesh.vertices.point(mesh.cells.vertex(cell, lv % nb_vertices)),
                mesh.vertices.point(mesh.cells.vertex(cell, (lv + 1) % nb_vertices)),
                mesh.vertices.point(mesh.cells.vertex(cell, (lv + 2) % nb_vertices)),
                mesh.vertices.point(mesh.cells.vertex(cell, (lv + 3) % nb_vertices))
            );
        }

        float edge_length[nb_edges];
        for (l_index le = 0; le < nb_edges; le++)
        {
            edge_length[le] = geogram::distance(mesh.vertices.point(mesh.cells.edge_vertex(cell, le, 0)), mesh.vertices.point(mesh.cells.edge_vertex(cell, le, 1)));
        }

        float max_edge_length = *std::max_element(edge_length, edge_length + nb_edges);
        float min_height = *std::min_element(heights, heights + nb_vertices);

        return std::sqrt(max_edge_length / min_height);
    }

    /**
     * @brief Computes the average aspect ratio of all tetrahedral elements of a given mesh.
     *
     * The aspect ratio defines the ratio of the longest edge to the shortest height (from a vertex to the opposite face) in a tetrahedron.
     * Ideal aspect ratio: 1.0.
     *
     * @param mesh Mesh.
     * @return float Average aspect-ratio.
     */
    PURE INLINE float aspect_ratio(const geogram::Mesh& mesh)
    {
        float aspect_ratio = 0.0f;
        for (const g_index cell : mesh.cells)
        {
            aspect_ratio += tetrahedron_aspect_ratio(cell, mesh);
        }

        return aspect_ratio / static_cast<float>(mesh.cells.nb());
    }

    /**
     * @brief Computes the average skewness of all tetrahedral elements of a given mesh.
     *
     * Skewness defines the deviation from an ideal equilateral tetrahedron.
     * Ideal skewness: 0.0, meaning that all angles are around 70.53 degrees.
     *
     * @param mesh
     * @return float
     */
    PURE INLINE float skewness(const geogram::Mesh& mesh)
    {
        return 0.0f;
    }

    /**
     * @brief TODO (defines "slivers")
     *
     * @param mesh TODO
     * @return float TODO
     */
    PURE INLINE float volume_ratio(const geogram::Mesh& mesh)
    {
        return 0.0f;
    }

    inline void compute(moist::metrics::MeshQuality& metrics, const geogram::Mesh& mesh)
    {
        metrics.aspect_ratio = 0; // aspect_ratio(mesh);
        metrics.skewness = 0;// skewness(mesh);
    }
}

#endif // MOIST_CORE_MESH_QUALITY_INL_
