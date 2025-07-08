#ifndef MOIST_CORE_MESH_QUALITY_INL_
#define MOIST_CORE_MESH_QUALITY_INL_

#include <iostream>
#include <string>
#include <format>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <tetrahedron.hpp>
#include <triangle.hpp>

#include "moist/core/defines.hpp"
#include "moist/core/metrics.hpp"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

namespace moist::mesh_quality
{

    namespace tet
    {
        // struct Sphere
        // {
        //     double radius;
        //     vec3 center;
        // };

        static void get_data(double* data, const g_index cell, const geo::Mesh& mesh)
        {
            #pragma unroll 4
            for (l_index lv = 0; lv < 4; lv++)
            {
                const vec3 p = mesh.vertices.point(mesh.cells.vertex(cell, lv));
                data[3 * lv] = p.x;
                data[3 * lv + 1] = p.y;
                data[3 * lv + 2] = p.z;
            }
        }

        /**
         * @brief Computes the aspect ratio (circumsphere-radius / (3 * inner-sphere-radius)) of a tetrahedron.
         *
         * @param cell
         * @param mesh
         * @return float
         */
        PURE inline float aspect_ratio(const g_index cell, const geo::Mesh& mesh)
        {
            const geo::index_t nb_vertices = mesh.cells.nb_vertices(cell);
            if (nb_vertices != 4) throw std::runtime_error("cell#nb_vertices != 4");

            double data[12];
            get_data(data, cell, mesh);

            // Sphere insphere;
            // Sphere circumsphere;
            //
            // tetrahedron_insphere(data, insphere.radius, insphere.center.data());
            // tetrahedron_circumsphere(data, circumsphere.radius, circumsphere.center.data());
            //
            // return circumsphere.radius / (3.0 * insphere.radius);

            return tetrahedron::tetrahedron_quality1(data);
        }

        /**
         * @brief Computes the mean ratio of a tetrahedron.
         *
         * The mean ratio defines, how close a tetrahedron is to an "ideal" or "regular" shape.
         *
         * mean_ratio = 1:   perfect tetrahedron
         * mean_ratio -> 0:  degenerate tetrahedron
         *
         * @param cell
         * @param mesh
         * @return float
         */
        PURE inline float mean_ratio(const g_index cell, const geo::Mesh& mesh)
        {
            const geo::index_t nb_vertices = mesh.cells.nb_vertices(cell);
            if (nb_vertices != 4) throw std::runtime_error("cell#nb_vertices != 4");

            double data[12];
            get_data(data, cell, mesh);

            return tetrahedron::tetrahedron_quality3(data);
        }
    }

    namespace tri
    {
        PURE inline float aspect_ratio(const g_index facet, const geo::Mesh& mesh)
        {
            const geo::index_t nb_vertices = mesh.facets.nb_vertices(facet);
            if (nb_vertices != 3) throw std::runtime_error("facet#nb_vertices != 3");

            double data[6];
            #pragma unroll 3
            for (l_index lv = 0; lv < 3; lv++)
            {
                const vec3 p = mesh.vertices.point(mesh.facets.vertex(facet, lv));
                data[2 * lv] = p.x;
                data[2 * lv + 1] = p.y;
                // TODO [Axis-Support]: Use a library / function that works in 3D.
                // crap, triangle only supports 2d triangles... somehow project? or just adapt and use own code...
                // since, for now, we only grow in the -z-direction, we can drop the z coordinate ("project" onto xy) and still keep the quality metric valid...
                // data[lv + 2] = p.z;
            }

            return triangle::triangle_quality(data);
        }
    }

    PURE INLINE float aspect_ratio(const geo::Mesh& mesh, /* ugly but fine */ moist::Interface* interface = nullptr)
    {
        float aspect_ratio = 0.0f;
        float max = -std::numeric_limits<double>::max();
        float min = std::numeric_limits<double>::max();

        size_t nb_cells = 0;
        for (const g_index cell : mesh.cells)
        {
            if (interface != nullptr && !moist::predicates::cell_on_plane(cell, mesh, *interface->Plane()))
            {
                continue;
            }

            const double l_aspect_ratio = tet::aspect_ratio(cell, mesh);
            if (l_aspect_ratio < min)
            {
                min = l_aspect_ratio;
            }

            if (l_aspect_ratio > max)
            {
                max = l_aspect_ratio;
            }

            aspect_ratio += l_aspect_ratio;
            nb_cells++;
        }

        return aspect_ratio / static_cast<float>(nb_cells);
    }

    PURE INLINE float mean_ratio(const geo::Mesh& mesh, moist::Interface* interface = nullptr)
    {
        float mean_ratio = 0.0f;
        size_t nb_cells = 0;
        for (const g_index cell : mesh.cells)
        {
            if (interface != nullptr && !moist::predicates::cell_on_plane(cell, mesh, *interface->Plane()))
            {
                continue;
            }

            mean_ratio += tet::mean_ratio(cell, mesh);
            nb_cells++;
        }

        return mean_ratio / static_cast<float>(nb_cells);
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
    PURE INLINE float skewness(const geo::Mesh& mesh)
    {
        return 0.0f;
    }

    inline void compute(moist::metrics::MeshQuality& metrics, const geo::Mesh& mesh)
    {
        metrics.aspect_ratio = aspect_ratio(mesh);
        metrics.mean_ratio = mean_ratio(mesh);
        metrics.skewness = skewness(mesh);

        metrics.nb_vertices = mesh.vertices.nb();
        metrics.nb_edges = 0; // not needed for eval
        metrics.nb_cells = mesh.cells.nb();
    }

    inline void compute(moist::metrics::MeshQuality& metrics, const geo::Mesh& mesh, moist::Interface& interface)
    {
        metrics.aspect_ratio = aspect_ratio(mesh, &interface);
        metrics.mean_ratio = mean_ratio(mesh, &interface);
        metrics.skewness = skewness(mesh);

        metrics.nb_vertices = 0;
        for (const g_index v : mesh.vertices)
        {
            if (moist::predicates::point_on_plane(mesh.vertices.point(v), *interface.Plane()))
            {
                metrics.nb_vertices++;
            }
        }

        metrics.nb_edges = moist::geometry::collect_interface_edges(mesh, *interface.Plane()).size();

        metrics.nb_cells = 0;
        for (const g_index c : mesh.cells)
        {
            if (moist::predicates::cell_on_plane(c, mesh, *interface.Plane()))
            {
                metrics.nb_cells++;
            }
        }
    }
}

#endif // MOIST_CORE_MESH_QUALITY_INL_
