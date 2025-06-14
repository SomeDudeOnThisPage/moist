#include "slice.hpp"

#include <unordered_set>
#include <optional>
#include <functional>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/process.h>
#include <limits>
#include <format>
#include <mutex>

#include "moist/core/metrics.hpp"
#include "moist/core/utils.hpp"
#include "moist/core/attributes.inl"
#include "moist/core/mesh_quality.inl"
#include "moist/core/predicates.inl"
#include "moist/core/geometry.inl"

#include "local_operations.hpp"

struct Vec3HashOperator
{
    std::size_t operator()(const vec3& v) const
    {
        std::hash<double> hasher;
        std::size_t hx = hasher(v.x);
        std::size_t hy = hasher(v.y);
        std::size_t hz = hasher(v.z);

        return hx ^ (hy << 1) ^ (hz << 2);
    }
};

struct Vec3EqualOperator
{
    bool operator()(const GEO::vecng<3, double>& a, const GEO::vecng<3, double>& b) const
    {
        return std::fabs(a.x - b.x) < moist::__DOUBLE_EPSILON &&
            std::fabs(a.y - b.y) < moist::__DOUBLE_EPSILON &&
            std::fabs(a.z - b.z) < moist::__DOUBLE_EPSILON;
    }
};

moist::MeshSlice::MeshSlice(geogram::index_t dimension, bool single_precision) : geogram::Mesh(dimension, single_precision)
{
    _deleted_tets = std::unordered_set<geogram::index_t>();
}

void moist::MeshSlice::CreateTetrahedra(const std::initializer_list<CreatedTetrahedon> tetrahedra)
{
    for (const CreatedTetrahedon tet : tetrahedra)
    {
        _created_tets.push_back(tet);
    }
}

void moist::MeshSlice::DeleteTetrahedra(const std::initializer_list<g_index> tetrahedra)
{
    for (const g_index tet : tetrahedra)
    {
        _deleted_tets.insert(tet);
    }
}

void moist::MeshSlice::InsertInterface(moist::Interface& interface, moist::metrics::TimeMetrics_ptr metrics)
{
    moist::descriptor::LocalInterface descriptor;
    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceVertices", metrics);
        this->InsertInterfaceVertices(interface);
    }

    {
        auto timer = moist::Timer("MeshSlice::InsertInterfaceEdges", metrics);
        this->InsertInterfaceEdges(interface);
    }

    {
        auto timer = moist::Timer("MeshSlice::InsertTetQuality", metrics);
        this->InsertTetQuality(interface);
    }

#ifndef NDEBUG
    // this->Validate(interface);
#endif
}

void moist::MeshSlice::InsertTetQuality(moist::Interface& interface)
{
    for (const g_index cell : this->cells)
    {
        if (!moist::predicates::cell_on_plane(cell, *this, *interface.Plane()))
        {
            continue;
        }

        // find corresponding interface facet to attach quality to
        // TODO: make a reducing list like "unmatched_facets" so we don't need to iterate all for each cell?
        // TODO: or simply parallelize this...
        for (const g_index facet : interface.Triangulation()->facets)
        {
            if (!moist::predicates::facet_matches_cell(cell, facet, *this, *interface.Triangulation()))
            {
                continue;
            }
        }
    }
}


void moist::MeshSlice::InsertInterfaceVertices(moist::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    for (const g_index v : triangulation->vertices)
    {
        const vec3 p = triangulation->vertices.point(v);
        std::mutex m_deleted_tets;
        size_t edge_insert_cells = 0;
    #ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        geogram::parallel_for(0, this->cells.nb(), [this, point, plane](const g_index cell)
        {
    #else
        for (const g_index cell : this->cells)
        {
    #endif // OPTION_PARALLEL_LOCAL_OPERATIONS

            {
                std::lock_guard<std::mutex> lock(m_deleted_tets);
                if (_deleted_tets.contains(cell))
                {
                    PARALLEL_CONTINUE;
                }
            }

            // 0, -0.583800, -1
            moist::predicates::PointInTet pit = moist::predicates::point_in_tet(*this, cell, p, true);
            if (pit == moist::predicates::PointInTet::FACET)
            {
                moist::operation::vertex_insert_1to3(*this, cell, p, *interface.Plane());
                PARALLEL_CONTINUE;
            }

            if (pit == moist::predicates::PointInTet::EDGE)
            {
                moist::operation::vertex_insert_1to2(*this, cell, p, *interface.Plane());
                edge_insert_cells++;
                PARALLEL_CONTINUE;
            }
        }
    #ifdef OPTION_PARALLEL_LOCAL_OPERATIONS
        );
    #endif // OPTION_PARALLEL_LOCAL_OPERATIONS

        if (edge_insert_cells > 0 && edge_insert_cells == 1)
        {
            OOC_WARNING("inserted vertex on edge - only one cell was split... possible boundary cell, otherwise invalid geometry...");
        }
        this->CreateTetrahedra();
    }

    this->FlushTetrahedra();
    OOC_DEBUG("mesh after vertex insertion: " << this->vertices.nb() << " vertices, " << this->cells.nb() << " cells");

    geogram::mesh_repair(*this, geogram::MeshRepairMode::MESH_REPAIR_COLOCATE);
    for (const g_index cell : this->cells)
    {
        if (geogram::mesh_cell_volume(*this, cell) == 0.0 || moist::geometry::has_duplicate_vertex(cell, *this))
        {
            OOC_DEBUG("0 volume cell after insert...");
        }
    }
    moist::utils::geo::save("after_insert.msh", *this);
}

bool moist::MeshSlice::CanMoveVertex(const g_index v, const vec3& p)
{
    std::vector<geogram::Sign> signs(this->_created_tets_idx.size());
    for (size_t c = 0; c < this->_created_tets_idx.size(); c++)
    {
        signs.push_back(geogram::geo_sgn(geogram::mesh_cell_volume(*this, _created_tets_idx[c])));
    }

    const auto original_pos = this->vertices.point(v);

    this->vertices.point(v) = p;
    for (size_t c = 0; c < this->_created_tets_idx.size(); c++)
    {
        if (_deleted_tets.contains(_created_tets_idx[c])) continue;

        const auto sign = geogram::geo_sgn(geogram::mesh_cell_volume(*this, _created_tets_idx[c]));
        if (sign != geogram::Sign::ZERO && sign != signs[c])
        {
            this->vertices.point(v) = original_pos;
            return false;
        }
    }

    this->vertices.point(v) = original_pos;
    return true;
}

void moist::MeshSlice::InsertInterfaceEdges(moist::Interface& interface)
{
    const auto triangulation = interface.Triangulation();
    const auto plane = interface.Plane();
    const auto edges = moist::geometry::collect_edges(*triangulation);

    uint32_t i = 0;
    uint32_t highest = 0;
    for (const auto edge : edges)
    {
        i++;
        this->_created_tets_idx.clear();

        vec3 p0 = triangulation->vertices.point(edge.v0);
        vec3 p1 = triangulation->vertices.point(edge.v1);
        p0.z = plane->extent;
        p1.z = plane->extent;

        // find all tetrahedra that lie on the line between the two points
        const auto cells = this->cells.nb();
        std::vector<CreatedTetrahedon> crossed_cells;

        size_t nb_created_cells = 0;
        size_t nb_crossed_cells = 0;
        geogram::Mesh dbg_crossed_cells(3);
        geogram::Mesh dbg_created_cells(3);
        bool p_eql = false;

        for (g_index cell = 0; cell < cells; cell++)
        {
            nb_created_cells = 0;
            if (_deleted_tets.contains(cell))
            {
                continue;
            }

            if (geogram::mesh_cell_volume(*this, cell) == 0.0)
            {
                p_eql = true;
                _deleted_tets.insert(cell);
                continue;
            }

            // insert vertices where the line crosses edges of the tetrahedra
            std::vector<CrossedEdge> crossed_edges;
            std::unordered_map<vec3, g_index, Vec3HashOperator, Vec3EqualOperator> created_vertices;
            for (g_index e_id = 0; e_id < this->cells.nb_edges(cell); e_id++)
            {
                const g_index v0 = this->cells.edge_vertex(cell, e_id, 0);
                const g_index v1 = this->cells.edge_vertex(cell, e_id, 1);
                const vec3 cp0 = this->vertices.point(v0);
                const vec3 cp1 = this->vertices.point(v1);

                if (!moist::predicates::edge_on_plane(cp0, cp1, *plane))
                {
                    continue;
                }

                // internally, vec3 and vec2 are just represented by double[{2|3}], so we can efficiently reinterpret vec2 from vec3
                if (!moist::predicates::xy::check_lines_aabb(
                    reinterpret_cast<const vec2&>(cp0),
                    reinterpret_cast<const vec2&>(cp1),
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1)
                ))
                {
                    continue;
                }

                const auto intersection_opt = moist::predicates::xy::get_line_intersection(
                    reinterpret_cast<const vec2&>(cp0),
                    reinterpret_cast<const vec2&>(cp1),
                    reinterpret_cast<const vec2&>(p0),
                    reinterpret_cast<const vec2&>(p1)
                );

                if (!intersection_opt.has_value())
                {
                    continue;
                }

                const vec3 new_vertex = vec3(intersection_opt.value().x, intersection_opt.value().y, plane->extent);
                const g_index p = created_vertices.find(new_vertex) != created_vertices.end()
                    ? created_vertices.at(new_vertex)
                    : this->vertices.create_vertex(new_vertex.data());

                created_vertices.emplace(new_vertex, p);
                crossed_edges.push_back({ v0, v1, p });
            #ifndef NDEBUG
                dbg_crossed_cells.cells.create_tet(
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 0))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 1))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 2))),
                    dbg_crossed_cells.vertices.create_vertex(this->vertices.point(this->cells.vertex(cell, 3)))
                );
            #endif // NDEBUG
            }

            if (crossed_edges.size() > 0)
            {
                nb_crossed_cells++;
            }

            switch (crossed_edges.size())
            {
                case 1:
                    moist::operation::edge_split_1to2(*this, cell, crossed_edges[0], *plane);
                    nb_created_cells += 2;
                    break;
                case 2:
                    if (crossed_edges[0].p == crossed_edges[1].p)
                    {
                        vec3 pp0 = this->vertices.point(this->cells.vertex(cell, 0));
                        vec3 pp1 = this->vertices.point(this->cells.vertex(cell, 1));
                        vec3 pp2 = this->vertices.point(this->cells.vertex(cell, 2));
                        vec3 pp3 = this->vertices.point(this->cells.vertex(cell, 3));
                        p_eql = true;
                        //moist::operation::edge_split_1to3(*this, cell, crossed_edges[0], crossed_edges[1], *plane);
                    }
                    else
                    {
                        moist::operation::edge_split_1to3(*this, cell, crossed_edges[0], crossed_edges[1], *plane);
                        nb_created_cells += 3;
                    }
                    break;
            default:
                    break;
            }

            this->CreateTetrahedra();

            if (nb_created_cells > highest)
            {
                highest = nb_created_cells;
                OOC_DEBUG("new highest: " << i);
            }
        }

        if (i == 97) moist::utils::geo::save("dbg-crossed-cells.msh", dbg_crossed_cells);
        if (i == 97) this->DebugMesh("dbg-created-cells.msh", this->_created_tets_idx);

        std::vector<g_index> moved_vertices;
        std::unordered_map<vec3, std::unordered_set<g_index>, Vec3HashOperator, Vec3EqualOperator> cluster;
        geogram::Attribute<int> v_discard(this->vertices.attributes(), ATTRIBUTE_DISCARD);
        for (const g_index cc : this->_created_tets_idx)
        {
            for (l_index lv = 0; lv < 4; lv++)
            {
                const g_index v = this->cells.vertex(cc, lv);
                const vec3 p = this->vertices.point(v);

                if (!moist::predicates::point_on_plane(p, *plane)) continue;
                if (!cluster.contains(p))
                {
                    cluster.emplace(p, std::unordered_set<g_index>());
                }
                cluster[p].insert(v);
            }
        }
        // move a vertex of every cell in a way that leads to no inversions (zeroeing tets is allowed...)
        for (const g_index cc : this->_created_tets_idx)
        {
            // find candidate vertices
            std::unordered_set<g_index> candidates;
            g_index candidate_edge_v = geogram::NO_VERTEX;

            for (l_index lv = 0; lv < 4; lv++)
            {
                const g_index v = this->cells.vertex(cc, lv);
                const vec3 p = this->vertices.point(v);
                if (v_discard[v])
                {
                    candidates.insert(v);
                }

                if (p == p0 || p == p1)
                {
                    candidate_edge_v = v;
                }
            }

            // """"simple"""" case: Edge ends in one of the cells' vertices
            if (candidate_edge_v != geogram::NO_VERTEX)
            {
                const auto to = this->vertices.point(candidate_edge_v);
                for (const g_index v : candidates)
                {
                    if (v == candidate_edge_v)
                    {
                        continue;
                    }

                    const auto from = this->vertices.point(v);
                    for (const g_index vm : cluster[from])
                    {
                        if (vm >= this->vertices.nb()) continue;
                        this->vertices.point(vm) = to;
                        cluster[to].insert(vm);
                        v_discard[vm] = false;
                    }
                }
            }
            // difficult case: tet needs to be decimated to not cause inversion to any other tet
            else if (!candidates.empty())
            {
                int ffffff = 0;

                if (candidates.size() > 1)
                {
                    g_index first = geogram::NO_VERTEX;
                    for (const g_index candidate : candidates)
                    {
                        if (v_discard[candidate])
                        {
                            first = candidate;
                            break;
                        }
                    }

                    if (first == geogram::NO_VERTEX)
                    {
                        OOC_WARNING("invalid geometry...");
                        continue;
                    }

                    for (const g_index candidate : candidates)
                    {
                        const vec3 to = this->vertices.point(first);
                        for (const g_index vm : cluster[this->vertices.point(candidate)])
                        {
                            const vec3 from = this->vertices.point(vm);
                            if (vm >= this->vertices.nb()) continue;

                            if (this->CanMoveVertex(vm, to))
                            {
                                this->vertices.point(vm) = this->vertices.point(first);
                                cluster[this->vertices.point(first)].insert(vm);
                                v_discard[vm] = true;//v_discard[first];
                            }
                        }
                    }
                }
            }
            else
            {
                OOC_WARNING("did not find suitable decimation direction...");
            }

            if (geogram::mesh_cell_volume(*this, cc) == 0.0 || moist::geometry::has_duplicate_vertex(cc, *this))
            {
                //_deleted_tets.insert(cc);
            }
        }

        for (const g_index c : this->_created_tets_idx)
        {
            if (geogram::mesh_cell_volume(*this, c) == 0.0 || moist::geometry::has_duplicate_vertex(c, *this))
            {
                _deleted_tets.insert(c);
            }
        }

        if (i == 67 || i == 68)
        {
            this->FlushTetrahedra();
            moist::utils::geo::save(std::format("dbg-{}.msh", i), *this);
        }
        // 64 is teh bad one
        if (i == 97) this->DebugMesh("dbg-after-decimation.msh", this->_created_tets_idx);

    #ifndef NDEBUG
        bool _found = false;
        for (g_index _cell = 0; _cell < this->cells.nb(); _cell++)
        {
            for (g_index e_id = 0; e_id < this->cells.nb_edges(_cell); e_id++)
            {
                const g_index v0 = this->cells.edge_vertex(_cell, e_id, 0);
                const g_index v1 = this->cells.edge_vertex(_cell, e_id, 1);
                const vec3 cp0 = this->vertices.point(v0);
                const vec3 cp1 = this->vertices.point(v1);

                if (cp0 == p0 && cp1 == p1 || cp0 == p1 && cp1 == p0)
                {
                    _found = true;
                    break;
                }
            }

            if (_found)
            {
                break;
            }
        }
    #endif // NDEBUG

    }

    this->FlushTetrahedra();

    moist::utils::geo::save("before_other_decimation.msh", *this);

    for (const g_index cell : this->cells)
    {
        if (moist::geometry::has_duplicate_vertex(cell, *this))
        {
            _deleted_tets.insert(cell);
        }
    }

    this->FlushTetrahedra();
}

void moist::MeshSlice::CreateTetrahedra()
{
    for (const auto tet : this->_created_tets)
    {
        const auto t = this->cells.create_tet(tet.v0, tet.v1, tet.v2, tet.v3);
        _created_tets_idx.push_back(t);

        // TODO [Bugs]: Figure out why this happens. Also, figure out why this leads to a validatably correct output...
        const auto volume = geogram::mesh_cell_volume(*this, t);
        if (volume == 0.0)
        {
        #ifndef NDEBUG
            const auto p0 = this->vertices.point(tet.v0);
            const auto p1 = this->vertices.point(tet.v1);
            const auto p2 = this->vertices.point(tet.v2);
            const auto p3 = this->vertices.point(tet.v3);
            OOC_WARNING("cell t0 " << t << " has zero volume: " << volume);
            geogram::Mesh dbg(3);
            dbg.cells.create_tet(
                dbg.vertices.create_vertex(p0.data()),
                dbg.vertices.create_vertex(p1.data()),
                dbg.vertices.create_vertex(p2.data()),
                dbg.vertices.create_vertex(p3.data())
            );
            moist::utils::geo::save(std::format("debug/dbg-zerotet-{}.msh", t), dbg);
        #endif // NDEBUG
            _deleted_tets.insert(t);
        }
    }

    this->_created_tets.clear();
}

void moist::MeshSlice::FlushTetrahedra()
{
    geogram::vector<geogram::index_t> tets_to_delete(this->cells.nb());
    for (geogram::index_t deleted_tet : _deleted_tets)
    {
        tets_to_delete[deleted_tet] = 1;
    }
    this->cells.delete_elements(tets_to_delete, false);
    _deleted_tets.clear();
}

// TODO [Testing]: Move this code into a gtest module...
void moist::MeshSlice::Validate(moist::Interface& interface)
{
    // validate vertices mesh -> interface
    for (const g_index v : this->vertices)
    {
        const vec3 p = this->vertices.point(v);
        if (moist::predicates::point_on_plane(p, *interface.Plane()))
        {
            bool found = false;
            for (const g_index v_i : interface.Triangulation()->vertices)
            {
                const vec3 p_i = interface.Triangulation()->vertices.point(v_i);
                if (p == p_i)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                OOC_DEBUG("validation: invalid point " << p << " does not exist in the interface");
            }
        }
    }

    // validate vertices interface -> mesh
    for (const g_index v_i : interface.Triangulation()->vertices)
    {
        const vec3 p_i = interface.Triangulation()->vertices.point(v_i);

        bool found = false;
        for (const g_index v : this->vertices)
        {
            const vec3 p = this->vertices.point(v);
            if (p == p_i)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            //OOC_DEBUG("validation: missing point " << p_i << " does not exist in the mesh");
        }
    }

    const auto interface_edges = moist::geometry::collect_edges(*interface.Triangulation());
    const auto mesh_interface_edges = moist::geometry::collect_edges(*interface.Triangulation(), *interface.Plane());

    for (const auto interface_edge : interface_edges)
    {
        const auto p0 = interface.Triangulation()->vertices.point(interface_edge.v0);
        const auto p1 = interface.Triangulation()->vertices.point(interface_edge.v1);

        bool found = false;
        for (const auto mesh_interface_edge : mesh_interface_edges)
        {
            const auto ep0 = this->vertices.point(mesh_interface_edge.v0);
            const auto ep1 = this->vertices.point(mesh_interface_edge.v1);

            if (ep0 == p0 && ep1 == p1 || ep0 == p1 && ep1 == p0)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_DEBUG("validation: missing edge (" << p0 << " -> " << p1 << ") in slice");
        }
    }

    for (const auto mesh_interface_edge : mesh_interface_edges)
    {
        const auto p0 = this->vertices.point(mesh_interface_edge.v0);
        const auto p1 = this->vertices.point(mesh_interface_edge.v1);

        bool found = false;
        for (const auto interface_edge : interface_edges)
        {


            const auto ep0 = interface.Triangulation()->vertices.point(interface_edge.v0);
            const auto ep1 = interface.Triangulation()->vertices.point(interface_edge.v1);

            if (ep0 == p0 && ep1 == p1 || ep0 == p1 && ep1 == p0)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            OOC_DEBUG("validation: additional edge (" << p0 << " -> " << p1 << ") in slice");
        }
    }
}

#ifndef NDEBUG
void moist::MeshSlice::DebugMesh(std::string file, std::vector<g_index>& tetrahedra)
{
    geogram::Mesh dbg(3);
    for (const g_index c : tetrahedra)
    {
        if (!_deleted_tets.contains(c) && c < this->cells.nb())
        {
            dbg.cells.create_tet(
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(c, 0))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(c, 1))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(c, 2))),
                dbg.vertices.create_vertex(this->vertices.point(this->cells.vertex(c, 3)))
            );
        }
    }

    moist::utils::geo::save(file, dbg);
}
#endif
