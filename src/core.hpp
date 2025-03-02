#ifndef __CORE_HPP
#define __CORE_HPP

#include <string>
#include <filesystem>
#include <iostream>

#define OOC_ERROR(msg) \
    std::cerr << "ERROR: " << msg << " (" << __FILE__ << ":" << __LINE__ << ")" << std::endl

#ifndef NDEBUG
#define OOC_DEBUG(msg) \
    std::cerr << msg << " (" << __FILE__ << ":" << __LINE__ << ")" << std::endl
#else
#define OOC_DEBUG(msg)
#endif

#ifdef GEOGRAM_API
#include <geogram/api/defs.h>
namespace geogram = GEO;
#endif // GEOGRAM_API

namespace ooc
{
    // something like this?
    // save sizing field in here too?
    // must be serializable / passable between steps
    // typedef struct
    // {
    //     /* bbox */ bounding_box;
    //     /* plane[] */ interface_planes;
    // } Slice;

    typedef struct
    {
        std::filesystem::path path_mesh_a;
        std::filesystem::path path_mesh_b;
        std::filesystem::path path_mesh_out;
        std::string plane;
    } InterfaceExtractionOptions;

}

#endif // __CORE_HPP
