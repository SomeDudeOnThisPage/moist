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

// dirty hack but no time...
std::ofstream moist::ScopeTimer::_csv_stream;
std::chrono::time_point<std::chrono::high_resolution_clock> moist::ScopeTimer::_timing_start;

moist::ScopeTimer::ScopeTimer(const std::string& name, const bool write_timing_direct) : Timer(name, nullptr), _write(write_timing_direct)
{
    if (write_timing_direct && !_csv_stream.is_open())
    {
        const std::filesystem::path csv("./" + name + "_timings.csv");
        bool exists = std::filesystem::exists(csv);
        _csv_stream = std::ofstream(csv, std::ios::app); // i never actually close these D:
        if (!_csv_stream.is_open())
        {
            throw std::runtime_error("Failed to open CSV file " + csv.string() + " for dumping timings...");
        }
        _csv_stream.precision(10); // let's keep it reasonable...
        _csv_stream << "name;timing;timing_runtime";
        _csv_stream << std::endl;
        _timing_start = std::chrono::high_resolution_clock::now();
    }
}

moist::ScopeTimer::~ScopeTimer()
{
    const auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - _start).count();
    double ms_runtime = std::chrono::duration<double, std::milli>(end - _timing_start).count();

    auto& data = moist::ScopeTimer::Timings()[_name];
    data.total += ms;
    data.count++;

    if (_write && _csv_stream.is_open())
    {
        auto end = std::chrono::steady_clock::now();
        _csv_stream << _name << ";" << ms << ";" << ms_runtime;
        _csv_stream << std::endl;
    }
}
