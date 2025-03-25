#ifndef __CORE_HPP
#define __CORE_HPP

#include <string>
#include <filesystem>
#include <iostream>

#define OOC_ERROR(msg) \
    std::cerr << "ERROR: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl; std::exit(EXIT_FAILURE)

#define OOC_WARNING(msg) \
    std::cerr << "WARNING: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl

#ifndef NDEBUG
#define OOC_DEBUG(msg) \
    std::cout << "DEBUG: [" << __FILE__ << ":" << __LINE__ << "]: " << msg  << std::endl
#else
#define OOC_DEBUG(msg)
#endif

#ifdef GEOGRAM_API
#include <geogram/api/defs.h>
#include <geogram/basic/stopwatch.h>
// not parallelizable, either go back to raw geogram stopwatch or make sure to only time on the main thread!
#define TIMER_START(x) geogram::Stopwatch __SCOPE_TIMER(x, false); std::string __SCOPE_TIMER_TASK_NAME = x;
#define TIMER_END OOC_DEBUG(__SCOPE_TIMER_TASK_NAME << " took " << __SCOPE_TIMER.elapsed_time() << "s");

namespace geogram = GEO;

using g_index = geogram::index_t;
using l_index = geogram::index_t;

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

// parallel
#define OPTION_PARALLEL_LOCAL_OPERATIONS
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    #define PARALLEL_CONTINUE return
    #define PARALLEL_BREAK return
#else
    #warning "compiling with non-parallel local operations"
    #define PARALLEL_CONTINUE continue
    #define PARALLEL_BREAK break
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

namespace incremental_meshing
{
    // something like this?
    // save sizing field in here too?
    // must be serializable / passable between steps
    // typedef struct
    // {
    //     /* bbox */ bounding_box;
    //     /* plane[] */ interface_planes;
    // } Slice;

    enum class Axis
    {
        X,
        Y,
        Z
    };

    typedef struct
    {
        std::filesystem::path path_mesh_a;
        std::filesystem::path path_mesh_b;
        std::filesystem::path path_mesh_out;
        double plane;
        double envelope_size;
        Axis axis;
    } InterfaceExtractionOptions;

}

#endif // __CORE_HPP
