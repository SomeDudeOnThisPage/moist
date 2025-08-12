#include "exact_mesh.hpp"

#include "geometry_exact.inl"

#ifndef NDEBUG
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_repair.h>
#include "moist/core/utils.hpp"
#endif // NDEBUG

moist::ExactMesh::ExactMesh() : _grid(std::make_shared<moist::LookupGridExact>()), _point_grid(moist::LookupPointGrid(100.0))
{
}

void moist::ExactMesh::ResetGrid(const double resolution)
{
    _grid->_grid.clear();
    _grid->_resolution = resolution; // TODO: Calculate!
    _grid->_min_bounds = vec2(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    _grid->_max_bounds = vec2(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    for (const auto& point : this->Points())
    {
        _grid->_min_bounds.x = std::min(_grid->_min_bounds.x, point.x());
        _grid->_min_bounds.y = std::min(_grid->_min_bounds.y, point.y());
        _grid->_max_bounds.x = std::max(_grid->_max_bounds.x, point.x());
        _grid->_max_bounds.y = std::max(_grid->_max_bounds.y, point.y());
    }

    _grid->_cell_size = (_grid->_max_bounds - _grid->_min_bounds) / _grid->_resolution;

    for (std::size_t c = 0; c < this->Cells().size(); c++)
    {
        _grid->InsertCell(c, *this);
    }
}

void moist::ExactMesh::ResetMesh()
{
    _points.clear();
    _edges.clear();
    _cells.clear();
    _pid = 0;
    _tid = 0;
    _eid = 0;
    _grid->_grid.clear();
}

std::size_t moist::ExactMesh::Add(const moist::exact::Point p)
{
    // check if the point already exists...
    const auto existing = _point_grid.Get(p);
    if (existing != std::nullopt)
    {
        return existing.value();
    }

    std::size_t index = _pid;
    _points.push_back(p);
    _point_grid.Add(p, _pid);
    _pid++;
    return index;
}

std::size_t moist::ExactMesh::Add(const moist::exact::Cell c, const bool initialization)
{
    std::size_t index = _tid;
    _cells.push_back(c);
    _tid++;
    if (!initialization)
    {
        _grid->InsertCell(index, *this);
    }
    return index;
}

std::size_t moist::ExactMesh::Add(const moist::exact::Edge e)
{
    std::size_t index = _eid;
    if (!_edges.contains(e))
    {
        _edges.insert(e);
        _eid++;
    }
    return index; // kinda meaningless lol
}

std::size_t moist::ExactMesh::Add(const moist::exact::Facet f)
{
    _facets.push_back(f);
    return _facets.size();
}

void moist::ExactMesh::DeletePoint(const std::size_t v)
{
    _points.at(v)._deleted = true;
}

void moist::ExactMesh::DeleteCell(const std::size_t c)
{
    _cells.at(c)._deleted = true;
}

moist::exact::Point &moist::ExactMesh::Point(const std::size_t& index)
{
    if (_points.size() < index)
    {
        throw std::runtime_error("Point index does not exist: " + index);
    }

    return _points.at(index);
}

const moist::exact::Point & moist::ExactMesh::Point(const std::size_t & index) const
{
    if (_points.size() < index)
    {
        throw std::runtime_error("Point index does not exist: " + index);
    }

    return _points.at(index);
}

moist::exact::Cell& moist::ExactMesh::Cell(const std::size_t& index)
{
    if (_cells.size() < index)
    {
        throw std::runtime_error("Cell index does not exist: " + index);
    }

    return _cells.at(index);
}

const moist::exact::Cell& moist::ExactMesh::Cell(const std::size_t& index) const
{
    if (_cells.size() < index)
    {
        throw std::runtime_error("Cell index does not exist: " + index);
    }

    return _cells.at(index);
}

#ifndef NDEBUG
void moist::ExactMesh::DebugMesh(const std::filesystem::path& file)
{
    geo::Mesh mesh(3);
    for (const auto& p : this->_points)
    {
        mesh.vertices.create_vertex(geo::vec3(CGAL::to_double(p._p.x()), CGAL::to_double(p._p.y()), CGAL::to_double(p._p.z())));
    }

    geo::vector<geo::index_t> to_delete(this->NbCells());
    for (const auto c : this->_cells)
    {
        if (c._deleted || moist::geometry::exact::is_degenerate(c))
        {
            continue;
        }

        const auto t = mesh.cells.create_tet(
            geo::index_t(c._points[0]),
            geo::index_t(c._points[1]),
            geo::index_t(c._points[2]),
            geo::index_t(c._points[3])
        );

        if (geo::mesh_cell_volume(mesh, t) == 0.0)
        {
            OOC_WARNING("degenerate tet " << t);
            to_delete[t] = true;
        }
    }
    mesh.cells.delete_elements(to_delete);

    for (const auto e : this->_edges)
    {
        if (e._deleted)
        {
            continue;
        }

        mesh.edges.create_edge(
            geo::index_t(e._points[0]),
            geo::index_t(e._points[1])
        );
    }

    for (const auto f : this->_facets)
    {
        mesh.facets.create_triangle(
            geo::index_t(f[0]),
            geo::index_t(f[1]),
            geo::index_t(f[2])
        );
    }

    //geo::mesh_repair(mesh, geo::MESH_REPAIR_COLOCATE);
    moist::utils::geogram::save(file, mesh);
}
#endif // NDEBUG
