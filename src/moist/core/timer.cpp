#include "timer.hpp"

moist::Timer::Timer(const std::string& name, moist::metrics::Metrics_ptr metrics) : _name(name), _start(std::chrono::system_clock::now()), _ended(false), _metrics(metrics)
{
}

moist::Timer::~Timer()
{
    if (_metrics != nullptr && !_ended)
    {
        const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _start).count();
        _metrics->metrics.push_back(moist::metrics::internal::metric_t {"timer::" + _name, ms});
    #ifdef PRINT_TIMES
        OOC_DEBUG(_name << ": " << ms << "ms");
    #endif // PRINT_TIMES
    }
}

void moist::Timer::End(bool write_to_metrics)
{
    if (!_ended && _metrics != nullptr && write_to_metrics)
    {
        const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _start).count();
        _metrics->metrics.push_back(moist::metrics::internal::metric_t {"timer::" + _name, ms});
    }
    _ended = true;
}

long long moist::Timer::Elapsed()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _start).count();
}

moist::ScopeTimer::ScopeTimer(const std::string& name) : Timer(name, nullptr)
{
}

moist::ScopeTimer::~ScopeTimer()
{
    const auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - _start).count();

    auto& data = moist::ScopeTimer::Timings()[_name];
    data.total += ms;
    data.count++;
}
