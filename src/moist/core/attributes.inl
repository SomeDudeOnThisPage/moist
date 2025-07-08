#ifndef MOIST_CORE__ATTRIBUTES_HPP_
#define MOIST_CORE__ATTRIBUTES_HPP_

#include <mutex>

#include <geogram/basic/attributes.h>
#include <geogram/basic/numeric.h>

#include "moist/core/defines.hpp"

#define VERTEX_FLAGS "v_flags"
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
#define ATTRIBUTE_INTERFACE_FACET_FOUND_AFTER_MERGE "FoundAfterMerge"
#define ATTRIBUTE_VERTEX_STEINER_POINT "VertexIsSteinerPoint"

// Interface Meta Attributes
#define ATTRIBUTE_META_PLANE_EXTENT "meta.extent"
#define ATTRIBUTE_META_PLANE_ENVELOPE "meta.envelope"
#define ATTRIBUTE_META_PLANE_AXIS "meta.axis"

namespace moist::attributes
{
    constexpr std::string VertexFlags = "v_flags";

    inline std::mutex _MUTEX_VERTEX_DESCRIPTOR;
    struct LockAttributes
    {
        LockAttributes() : lock(_MUTEX_VERTEX_DESCRIPTOR) {}
    private:
        std::lock_guard<std::mutex> lock;
    };

    #define LOCK_ATTRIBUTES std::lock_guard<std::mutex> ___lock(moist::attributes::_MUTEX_VERTEX_DESCRIPTOR)
    constexpr std::string V_INTERFACE = "v_interface";

    enum class VertexDescriptorFlags : uint8_t
    {
        /** @brief Marks a vertex as deleted. */
        DELETED     = 1 << 1,
        /** @brief Marks a vertex as discardable for decimation. */
        DISCARD     = 1 << 2,
        /** @brief Marks a vertex as part of the interface plane. */
        INTERFACE   = 1 << 3
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
        geo::geo_register_attribute_type<geo::Numeric::uint8>("VertexDescriptorFlags");
    }
}
#endif // MOIST_CORE__ATTRIBUTES_HPP_
