#ifndef __OOC_ATTRIBUTES_HPP
#define __OOC_ATTRIBUTES_HPP

#include <mutex>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"

#define ATTRIBUTE_VERTEX_DESCRIPTOR_FLAGS "VertexDescriptorFlags"
#define ATTRIBUTE_DISCARD "Discard"
#define ATTRIBUTE_INTERFACE "Interface" // TODO: Make all this flags after testing, maybe keep "int" attributes for visualization in debug mode vorpaview
#define ATTRIBUTE_CLUSTER_ONTO "ClusterOnto" // TODO: Make all this flags after testing, maybe keep "int" attributes for visualization in debug mode vorpaview
#define ATTRIBUTE_DEBUG_V "Debug"

#define ATTRIBUTE_CONSTRAINT_EDGE "ConstraintEdge"
#define ATTRIBUTE_INTERFACE_INDEX_A "InterfaceIndexA"
#define ATTRIBUTE_INTERFACE_INDEX_B "InterfaceIndexB"
#define ATTRIBUTE_INTERFACE_TARGET_VERTICES "InterfaceTargetVertices"
#define ATTRIBUTE_INTERFACE_TETMERGE_QUALITY "TetrahedralMergeQuality"

// Interface Meta Attributes
#define ATTRIBUTE_META_PLANE_EXTENT "meta.extent"
#define ATTRIBUTE_META_PLANE_ENVELOPE "meta.envelope"
#define ATTRIBUTE_META_PLANE_AXIS "meta.axis"

#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    #define LOCK_ATTRIBUTES std::lock_guard<std::mutex> lock(incremental_meshing::attributes::_MUTEX_VERTEX_DESCRIPTOR)
#else
    #define LOCK_ATTRIBUTES
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS
#include <geogram/basic/numeric.h>
namespace moist::attributes
{
#ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
    inline std::mutex _MUTEX_VERTEX_DESCRIPTOR;
#endif // OPTION_PARALLEL_LOCAL_OPERATIONS

    /// TODO: Doxygen...
    enum class VertexDescriptorFlags : uint8_t
    {
        DELETED = 1 << 1, // marks vertex as deleted
        DISCARD = 1 << 2  // marks vertex as discardable
    };

    // needs to be defined for the compiler not to complain
    inline std::istream& operator>>(std::istream& is, VertexDescriptorFlags& flags){
        return is;
    }

    inline std::ostream& operator<<(std::ostream& os, const VertexDescriptorFlags& flags) {
        return os;
    }

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
        //geogram::geo_register_attribute_type<geogram::Numeric::uint8>("VertexDescriptorFlags");
        //geogram::geo_register_attribute_type<int>("Discard");
    }
}
#endif // __OOC_ATTRIBUTES_HPP
