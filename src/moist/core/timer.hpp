#ifndef MOIST_CORE_TIMER_HPP_
#define MOIST_CORE_TIMER_HPP_

#include <string>
#include <vector>
#include <memory>

#include "moist/core/defines.hpp"
#include "moist/core/metrics.hpp"

#define PRINT_TIMES

namespace moist
{
    /**
     * @brief Simple scoped timer writing durations in ms to a stats-object when it goes out of scope or is destructed.
     */
    class Timer
    {
    public:
        Timer(const std::string& name, metrics::Metrics_ptr metrics = nullptr);
        ~Timer();

        long long Elapsed();
        void End(bool write_to_metrics = true);
    protected:
        std::string _name;
        std::chrono::time_point<std::chrono::system_clock> _start;
    private:
        metrics::Metrics_ptr _metrics;
        bool _ended;
    };

#ifndef NDEBUG
    class ScopeTimer : Timer
    {
    public:
        ScopeTimer(const std::string& name);
        ~ScopeTimer();

        static void Print()
        {
            OOC_DEBUG("=========== RUNTIMES ===========");
            for (auto& [name, data] : ScopeTimer::Timings())
            {
                OOC_DEBUG(name << ": " << (data.total / data.count) << "ms x " << data.count << " => " << data.total << "ms");
            }
        }
    private:
        struct ScopeTime
        {
            double total = 0.0;
            std::size_t count = 0;
        };

        static std::unordered_map<std::string, ScopeTime>& Timings()
        {
            static std::unordered_map<std::string, ScopeTime> data;
            return data;
        }
    };
#endif // NDEBUG
}

#endif // MOIST_CORE_TIMER_HPP_
