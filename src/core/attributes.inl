#ifndef __OOC_ATTRIBUTES_HPP
#define __OOC_ATTRIBUTES_HPP

#include <mutex>

#include <geogram/mesh/mesh.h>

#include "core.hpp"

#define ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS "VertexDescriptorFlags"
#define ATTRIBUTE_DISCARD "Discard"
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    #define LOCK_ATTRIBUTES std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR)
#else
    #define LOCK_ATTRIBUTES
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

namespace incremental_meshing::attributes
{
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    inline std::mutex _MUTEX_VERTEX_DESCRIPTOR;
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    /// TODO: Doxygen...
    enum class VertexDescriptorFlags : uint8_t
    {
        INTERFACE = 1 << 0, // marks vertex as interface vertex for the currently processed plane
        DISCARD   = 1 << 1  // marks vertex as discardable in the decimation stage
        // DELETED = 1 << 3 // TODO: use this instead of a "Deleted"-Map, so it's 1 bit instead of 32 bit to mark a deleted tet?
    };

    // needs to be defined for the compiler not to complain
    inline std::istream& operator>>(std::istream& is, VertexDescriptorFlags& flags){ return is; }
    inline std::ostream& operator<<(std::ostream& os, const VertexDescriptorFlags& flags) { return os; }

    // must be an enum class, geogram doesn't seem to like "normal" pre-c11 enums
    // https://stackoverflow.com/questions/32578638/how-to-use-c11-enum-class-for-flags
    inline VertexDescriptorFlags operator|(const VertexDescriptorFlags lhs, const VertexDescriptorFlags rhs)
    {
        return static_cast<VertexDescriptorFlags>(
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(lhs) |
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(rhs)
        );
    }

    inline VertexDescriptorFlags operator&(VertexDescriptorFlags lhs, VertexDescriptorFlags rhs)
    {
        return static_cast<VertexDescriptorFlags>(
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(lhs) &
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(rhs)
        );
    }

    inline bool operator==(const VertexDescriptorFlags lhs, const VertexDescriptorFlags rhs)
    {
        return (1 !=
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(lhs) &
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(rhs)
        );
    }

    inline VertexDescriptorFlags& operator|=(VertexDescriptorFlags& lhs, VertexDescriptorFlags rhs)
    {
        lhs = static_cast<VertexDescriptorFlags>(
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(lhs) |
            static_cast<std::underlying_type<VertexDescriptorFlags>::type>(rhs)
        );
        return lhs;
    }

    inline void initialize()
    {
        geogram::geo_register_attribute_type<VertexDescriptorFlags>("VertexDescriptorFlags");
        geogram::geo_register_attribute_type<int>("Discard");
    }
}
#endif // __OOC_ATTRIBUTES_HPP
