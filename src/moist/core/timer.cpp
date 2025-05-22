#include "timer.hpp"

moist::Timer::Timer(const std::string& name, moist::metrics::TimeMetrics_ptr metrics) : _name(name), _start(std::chrono::system_clock::now()), _metrics(metrics)
{
}

moist::Timer::~Timer()
{
    if (_metrics != nullptr)
    {
        const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _start).count();
        _metrics->metrics.push_back(moist::metrics::internal::time_metric_t {_name, ms});
    #ifdef PRINT_TIMES
        OOC_DEBUG(_name << ": " << ms << "ms");
    #endif // PRINT_TIMES
    }
}

long long moist::Timer::Elapsed()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _start).count();
}
