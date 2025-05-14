#ifndef __OOC_MESH_HPP
#define __OOC_MESH_HPP

#include <unordered_set>
#include <initializer_list>

#include <geogram/mesh/mesh.h>

#include "core.hpp"
#include "core-interface.hpp"

namespace moist
{
    typedef struct CROSSED_EDGE_STRUCT
    {
        g_index e_v0;
        g_index e_v1;
        g_index p;

        g_index e_interface; // (custom) interface edge global index, used in decimation to find a point to decimate onto

        bool operator==(const CROSSED_EDGE_STRUCT& other) const
        {
            return p == other.p && ((e_v0 == other.e_v0 && e_v1 == other.e_v1) || (e_v0 == other.e_v1 && e_v1 == other.e_v0));
        }
    } CrossedEdge;

    typedef struct
    {
        g_index v0; // alias v_line
        g_index v1;
        g_index v2;
        g_index v3;
    } CreatedTetrahedon;

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
        void InsertInterface(moist::Interface& interface);

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

        // global operations
        void InsertInterfaceVertices(moist::Interface& interface);
        void InsertInterfaceEdges(moist::Interface& interface);

        // utils
        void InsertVertex(const geogram::vec3& point, const moist::AxisAlignedInterfacePlane& plane);
        void CreateTetrahedra();
        void FlushTetrahedra();
    };
}

#endif // __OOC_MESH_HPP
