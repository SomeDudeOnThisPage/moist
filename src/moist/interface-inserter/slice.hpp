#ifndef MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_
#define MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_

#include <unordered_set>
#include <initializer_list>

#include <geogram/mesh/mesh.h>

#include "moist/core/defines.hpp"
#include "moist/core/timer.hpp"
#include "moist/core/descriptor.hpp"
#include "moist/core/core_interface.hpp"

#include "interface_inserter.hpp"

namespace moist
{
    class MeshSlice : public geogram::Mesh
    {
    public:
        /**
         * @brief Constructs a MeshSlice object.
         *
         * Initializes a new MeshSlice with the given dimension and precision settings.
         *
         * @param dimension The dimension of the mesh (default `3`).
         * @param single_precision If true, use `float`; otherwise, use `double` (default `false`).
         */
        MeshSlice(geogram::index_t dimension = 3, bool single_precision = false);

        /**
         * @brief Inserts a created Interface into this MeshSlice.
         *
         * @param interface The (initialized) interface reference.
         */
        void InsertInterface(moist::Interface& interface, moist::metrics::TimeMetrics_ptr metrics = nullptr);

        void CreateTetrahedra(const CreatedTetrahedon tet) { this->CreateTetrahedra({tet}); }
        void CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra);
        void DeleteTetrahedra(const g_index tet) { this->DeleteTetrahedra({tet}); };
        void DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra);
    private:

        // deprecated
        std::string _identifier;

        std::vector<g_index> _created_tets_idx;
        std::vector<CreatedTetrahedon> _created_tets;

        std::unordered_set<geogram::index_t> _deleted_tets;

        std::vector<moist::descriptor::LocalInterface> _interfaces;

        // global operations
        void InsertInterfaceVertices(moist::Interface& interface);
        void InsertInterfaceEdges(moist::Interface& interface);
        void InsertTetQuality(moist::Interface& interface);

        void InsertVertex(const geogram::vec3& point, const moist::AxisAlignedInterfacePlane& plane);
        void CreateTetrahedra();
        void FlushTetrahedra();

        void Validate(moist::Interface& interface);
    };
}

#endif // MOIST_INTERFACE_INSERTER_MESH_SLICE_HPP_
