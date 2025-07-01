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
    private:
        std::string _name;
        std::chrono::time_point<std::chrono::system_clock> _start;
        metrics::Metrics_ptr _metrics;
    };
}

#endif // MOIST_CORE_TIMER_HPP_
