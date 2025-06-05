#include "test_utils.hpp"

#include <string>

#include <geogram/basic/process.h>

#include "moist/core/predicates.inl"
#include "moist/core/attributes.inl"

bool moist::test::contains_constraints(geogram::Mesh& mesh, geogram::Mesh& constraints)
{
    constexpr std::string_view attribute = "__contains_constraints";
    geogram::parallel_for(0, mesh.cells.nb(), [&mesh, &constraints, &attribute] (const g_index c)
    {
        for (const g_index f : constraints.facets)
        {
            if (!moist::predicates::facet_matches_cell(c, f, mesh, constraints))
            {
                continue;
            }

            LOCK_ATTRIBUTES;
            geogram::Attribute<double> f_found(constraints.facets.attributes(), std::string(attribute));
            if (!f_found[f])
            {
                f_found[f] = 1;
            }
            break;
        }
    });

    geogram::Attribute<double> f_found(constraints.facets.attributes(), std::string(attribute));
    for (const g_index f : constraints.facets)
    {
        if (!f_found[f])
        {
            return false;
        }
    }

    return true;
}
