#include "metrics.hpp"

moist::metrics::internal::time_metrics_t::time_metrics_t(const std::string &name) : name(name), metrics(std::vector<moist::metrics::internal::time_metric_t>())
{

}

void moist::metrics::internal::time_metrics_t::AppendCSV(const std::filesystem::path &csv)
{
}

moist::metrics::TimeMetrics_ptr moist::metrics::TimeMetrics(const std::string &name)
{
    return std::make_shared<internal::time_metrics_t>(name);
}
