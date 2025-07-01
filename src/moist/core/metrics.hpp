#ifndef MOIST_CORE_METRICS_HPP_
#define MOIST_CORE_METRICS_HPP_

#include <iostream>
#include <string>
#include <format>
#include <variant>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include "moist/core/defines.hpp"

namespace moist::metrics
{
    namespace internal
    {
        using timestamp_ms_t = long long;
        using metric_value_t = std::variant<timestamp_ms_t, double, unsigned int, int, size_t, std::string>;

        struct metric_t
        {
            std::string name;
            metric_value_t metric;

            metric_t(const std::string& name, metric_value_t metric) : name(name), metric(metric) {};
        };

        struct metrics_t
        {
            /** @brief Name of the metric object to be dumped into a filename, one line per dump. Name is the first column. */
            std::string name;
            /** @brief List of TimeMetrics taken. Need to be equal to the columns in the CSV, if it exists. If not, the headers will be written. */
            std::vector<metric_t> metrics;

            metrics_t(const std::string& name);
            void AppendCSV(const std::filesystem::path& csv);

            metrics_t& operator<<(const metric_t& m)
            {
                metrics.push_back(m);
                return *this;
            }
        };
    }

    using Metrics_ptr = std::shared_ptr<internal::metrics_t>;
    Metrics_ptr Metrics(const std::string& name);
    using Metric = internal::metric_t;

    struct MeshQuality
    {
        std::string name;
        size_t nb_vertices;
        size_t nb_cells;
        size_t nb_edges;

        float aspect_ratio;
        float mean_ratio;
        float skewness;
        MeshQuality(const std::string& name) : name(name) {};
    };

    inline std::ostream& operator<<(std::ostream& os, const MeshQuality& metrics)
    {
        // os << std::format("aspect ratio: {:.5f}; skewness: {:.5f}", metrics.aspect_ratio, metrics.skewness);
        return os;
    }

    inline internal::metrics_t& operator<<(internal::metrics_t& metrics, const MeshQuality& mesh_quality_metrics)
    {
        metrics << internal::metric_t {mesh_quality_metrics.name + "::nb_vertices", mesh_quality_metrics.nb_vertices}
                << internal::metric_t {mesh_quality_metrics.name + "::nb_edges", mesh_quality_metrics.nb_edges}
                << internal::metric_t {mesh_quality_metrics.name + "::nb_cells", mesh_quality_metrics.nb_cells}
                << internal::metric_t {mesh_quality_metrics.name + "::aspect_ratio", mesh_quality_metrics.aspect_ratio}
                << internal::metric_t {mesh_quality_metrics.name + "::mean_ratio", mesh_quality_metrics.mean_ratio};
        return metrics;
    }
}

#endif // MOIST_CORE_METRICS_HPP_
