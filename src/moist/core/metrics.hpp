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
    using Histogram = std::vector<std::size_t>;

    namespace internal
    {
        using timestamp_ms_t = long long;
        using metric_value_t = std::variant<timestamp_ms_t, double, unsigned int, int, size_t, std::string, Histogram>;

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

    // very simple json plopper
    inline std::ostream& operator<<(std::ostream& os, const Histogram& h)
    {
        os << "{ \"bins\": [";
        for (std::size_t i = 0; i < h.size(); i++)
        {
            os << h[i];
            if (i + 1 < h.size()) os << ",";
        }
        os << "] }";
        return os;
    }

    enum class QualityMetricType
    {
        ASPECT_RATIO,
        MEAN_RATIO,
        MMG3D_QUALITY
    };

    struct MeshQuality
    {
        std::string name;
        size_t nb_vertices;
        size_t nb_cells;

        double aspect_ratio;
        double mean_ratio;
        double aspect_ratio_ub;
        double mean_ratio_ub;
        double aspect_ratio_lb;
        double mean_ratio_lb;
        double mmg;
        double mmg_lb;
        double mmg_ub;
        double skewness;
        MeshQuality(const std::string& name) : name(name) {};
    };

    inline std::ostream& operator<<(std::ostream& os, const MeshQuality& metrics)
    {
        // os << std::format("aspect ratio: {:.5f}; skewness: {:.5f}", metrics.aspect_ratio, metrics.skewness);
        return os;
    }

    inline MeshQuality operator+(const MeshQuality& a, const MeshQuality& b)
    {
        MeshQuality q(a.name + b.name);
        q.nb_vertices = a.nb_vertices + b.nb_vertices;
        q.nb_cells = a.nb_cells + b.nb_cells;
        q.aspect_ratio = (a.aspect_ratio + b.aspect_ratio) / 2.0;
        q.aspect_ratio_lb = std::min(a.aspect_ratio_lb, b.aspect_ratio_lb);
        q.aspect_ratio_ub = std::max(a.aspect_ratio_ub, b.aspect_ratio_ub);
        q.mean_ratio = (a.mean_ratio + b.mean_ratio) / 2.0;
        q.mean_ratio_lb = std::min(a.mean_ratio_lb, b.mean_ratio_lb);
        q.mean_ratio_ub = std::max(a.mean_ratio_ub, b.mean_ratio_ub);
        q.mmg = (a.mmg + b.mmg) / 2.0;
        q.mmg_lb = std::min(a.mmg_lb, b.mmg_lb);
        q.mmg_ub = std::max(a.mmg_ub, b.mmg_ub);
        return q;
    }

    inline internal::metrics_t& operator<<(internal::metrics_t& metrics, const MeshQuality& mesh_quality_metrics)
    {
        metrics << internal::metric_t {mesh_quality_metrics.name + "::nb_vertices", mesh_quality_metrics.nb_vertices}
                << internal::metric_t {mesh_quality_metrics.name + "::nb_cells", mesh_quality_metrics.nb_cells}
                << internal::metric_t {mesh_quality_metrics.name + "::aspect_ratio", mesh_quality_metrics.aspect_ratio}
                << internal::metric_t {mesh_quality_metrics.name + "::aspect_ratio::upper_bound", mesh_quality_metrics.aspect_ratio_ub}
                << internal::metric_t {mesh_quality_metrics.name + "::aspect_ratio::lower_bound", mesh_quality_metrics.aspect_ratio_lb}
                << internal::metric_t {mesh_quality_metrics.name + "::mean_ratio", mesh_quality_metrics.mean_ratio}
                << internal::metric_t {mesh_quality_metrics.name + "::mean_ratio::upper_bound", mesh_quality_metrics.mean_ratio_ub}
                << internal::metric_t {mesh_quality_metrics.name + "::mean_ratio::lower_bound", mesh_quality_metrics.mean_ratio_lb}
                << internal::metric_t {mesh_quality_metrics.name + "::mmg", mesh_quality_metrics.mmg}
                << internal::metric_t {mesh_quality_metrics.name + "::mmg::upper_bound", mesh_quality_metrics.mmg_ub}
                << internal::metric_t {mesh_quality_metrics.name + "::mmg::lower_bound", mesh_quality_metrics.mmg_lb};

        return metrics;
    }
}

#endif // MOIST_CORE_METRICS_HPP_
