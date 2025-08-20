#ifndef MOIST_CORE_CORE_HPP_
#define MOIST_CORE_CORE_HPP_

#include <format>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <source_location>

#ifndef NDEBUG
    #define OPTION_LOOKUP_GRID
    #define MOIST_OPTION_EXACT_PREDICATES
// #define OPTION_PARALLEL_LOCAL_OPERATIONS // dev
// #define OPTION_DEBUG_TEST_INTERFACE
#endif // NDEBUG

#ifndef NDEBUG
  //#define debug(...) do { __VA_ARGS__; } while (0)
#else
  //#define debug(...) ((void)0)
#endif

namespace moist
{
    /**
     * @brief Axis along which a mesh can be incrementally constructed.
     */
    enum class Axis
    {
        X = 1 << 0,
        Y = 2 << 1,
        Z = 3 << 2
    };

//! @internal
    const double __DOUBLE_EPSILON = std::numeric_limits<double>::epsilon() * 3.0; // https://stackoverflow.com/questions/48133572/what-can-stdnumeric-limitsdoubleepsilon-be-used-for

    const std::map<std::string, moist::Axis> AXIS_OPTION_ARGUMENT_MAP
    {
        {"X", moist::Axis::X},
        {"Y", moist::Axis::Y},
        {"Z", moist::Axis::Z}
    };
//! @endinternal

    class Context
    {

    };
}

namespace moist::assert
{
    inline void in_range(std::size_t index, std::size_t lower, std::size_t upper, const std::source_location& loc = std::source_location::current())
    {
        if (!(index >= lower && index < upper))
        {
            throw std::out_of_range
            (
                "index " + std::to_string(index) + " out of range [" + std::to_string(lower) + ", " + std::to_string(upper) + ") at "
                + loc.file_name() + ":" + std::to_string(loc.line()) + " in function " + loc.function_name()
            );
        }
    }
}

//! @internal
#ifdef GEOGRAM_API
    #define OOC_ERROR(msg) \
        geo::Logger::err("MOIST") << "[" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl; std::exit(EXIT_FAILURE)

    #define OOC_WARNING(msg) \
        geo::Logger::warn("MOIST") << "[" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl
    #define MOIST_INFO(msg) \
        geo::Logger::out("MOIST") << "[" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl

    #ifndef NDEBUG
    #define OOC_DEBUG(msg) \
        geo::Logger::out("MOIST") << "[" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl
    #else
    #define OOC_DEBUG(msg)
    #endif
#else
    #define OOC_ERROR(msg) \
        std::cout << "ERROR: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl; std::exit(EXIT_FAILURE)

    #define OOC_WARNING(msg) \
        std::cout << "WARNING: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl

    #ifndef NDEBUG
    #define OOC_DEBUG(msg) \
        std::cout << "DEBUG: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl
    #else
    #define OOC_DEBUG(msg)
    #endif
#endif

#ifdef CLI11_INLINE
    namespace cli = CLI;
#endif

#ifdef GEOGRAM_API
#include <geogram/api/defs.h>
#include <geogram/basic/common.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/stopwatch.h>
// not parallelizable, either go back to raw geogram stopwatch or make sure to only time on the main thread!
#define TIMER_START(x) geogram::Stopwatch __SCOPE_TIMER(x, false); std::string __SCOPE_TIMER_TASK_NAME = x;
#define TIMER_END OOC_DEBUG(__SCOPE_TIMER_TASK_NAME << " took " << __SCOPE_TIMER.elapsed_time() << "s");

namespace geo = GEO;

// custom index aliases to make differentiating between local and global vertices easier
using g_index = geo::index_t;
using l_index = geo::index_t;

// in case we switch to real predicates, similar to, for instance, marco attenes stuff...
using vec3 = geo::vec3;
using vec2 = geo::vec2;

// additional operators
/*inline bool operator==(const vec3& a, const vec3& b)
{
    return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
           std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON &&
           std::fabs(a.z - b.z) < moist::__DOUBLE_EPSILON;
}*/

inline bool operator==(const geo::vec3& a, const geo::vec3& b)
{
    return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
           std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON &&
           std::fabs(a.z - b.z) < moist::__DOUBLE_EPSILON;
    // return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool operator==(const geo::vec2& a, const geo::vec2& b)
{
    return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
           std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON;
}

#endif // GEOGRAM_API

// compiler-specific macros
#ifdef _MSC_VER
    #define INLINE __forceinline
    #define PURE
#elif defined(__GNUC__) || defined(__clang__)
    #define INLINE __attribute__((always_inline)) inline
    #define PURE __attribute__((pure))
#else
    #define INLINE inline
    #define PURE
#endif

#ifdef __CUDACC__
    #define restricted __restrict__
#endif

#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    #define PARALLEL_CONTINUE return
    #define PARALLEL_BREAK return
#else
//    #warning "compiling with non-parallel local operations"
    #define PARALLEL_CONTINUE continue
    #define PARALLEL_BREAK break
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

// Since float and double are aliased in cmake options, alias them in code aswell.
using float32_t = float;
using float64_t = double;
//! @endinternal

#endif // MOIST_CORE_CORE_HPP_
