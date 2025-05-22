#ifndef MOIST_CORE_METRICS_HPP_
#define MOIST_CORE_METRICS_HPP_

#include <iostream>
#include <string>
#include <format>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include "moist/core/defines.hpp"

namespace moist::metrics
{
    namespace internal
    {
        struct time_metric_t
        {
            /** @brief Identifier for this statistic. */
            std::string name;
            /** @brief Duration. */
            long long ms;
        };

        struct time_metrics_t
        {
            /** @brief Name of the metric object to be dumped into a filename. */
            std::string name;
            /** @brief List of TimeMetrics taken. */
            std::vector<time_metric_t> metrics;

            time_metrics_t(const std::string& name);

            void AppendCSV(const std::filesystem::path& csv);
        };
    }

    using TimeMetrics_ptr = std::shared_ptr<internal::time_metrics_t>;
    TimeMetrics_ptr TimeMetrics(const std::string& name);

    struct MeshQuality
    {
        float aspect_ratio;
        float skewness;
    };

    inline std::ostream& operator<<(std::ostream& os, const MeshQuality& metrics)
    {
        // os << std::format("aspect ratio: {:.5f}; skewness: {:.5f}", metrics.aspect_ratio, metrics.skewness);
        return os;
    }
}

#endif // MOIST_CORE_METRICS_HPP_
