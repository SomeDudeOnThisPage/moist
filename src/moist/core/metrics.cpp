#include "metrics.hpp"

#include <iostream>
#include <fstream>
#include <variant>

moist::metrics::internal::metrics_t::metrics_t(const std::string& name): name(name), metrics(std::vector<moist::metrics::internal::metric_t>())
{

}

void moist::metrics::internal::metrics_t::AppendCSV(const std::filesystem::path &csv)
{
    bool exists = std::filesystem::exists(csv);
    std::ofstream out(csv, std::ios::app);
    if (!out.is_open())
    {
        throw std::runtime_error("Failed to open CSV file " + csv.string() + " for dumping metrics...");
    }

    out.precision(10); // let's keep it reasonable...

    if (!exists)
    {
        out << "metric_name";
        for (const auto& m : metrics)
        {
            out << "," << m.name;
        }
        out << std::endl;
    }

    out << this->name; // name identifying the run, more for debugging purposes than actually useful

    // Just assume the order is always the same (as it should be when running the program without modifications back to back).
    for (const auto& m : metrics)
    {
        std::visit([&out](auto&& metric)
        {
            out << "," << metric; // all used types in variant defined in metrics.hpp are safe to just dump into ostream...
        }, m.metric);
    }
    out << std::endl;
}

moist::metrics::Metrics_ptr moist::metrics::Metrics(const std::string &name)
{
    return std::make_shared<internal::metrics_t>(name);
}
