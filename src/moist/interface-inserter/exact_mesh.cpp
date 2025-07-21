#include "exact_mesh.hpp"

#ifndef NDEBUG
#include <geogram/mesh/mesh.h>
#include "moist/core/utils.hpp"
#endif // NDEBUG

std::size_t moist::ExactMesh::Add(const ExactPoint& p)
{
    size_t index = _pid;
    // _points[index] = p;
    // _points.emplace(index, p);
    _points.push_back(p);
    _pid++;
    return index;
}

std::size_t moist::ExactMesh::Add(const ExactCell& c)
{
    size_t index = _tid;
    // _cells[index] = c;
    // _cells.emplace(index, c);
    _cells.push_back(c);
    _tid++;
    return index;
}

void moist::ExactMesh::DeletePoint(const std::size_t v)
{
    _points.at(v)._deleted = true;
}

void moist::ExactMesh::DeleteCell(const std::size_t c)
{
    _cells.at(c)._deleted = true;
}

const moist::ExactMesh::ExactPoint &moist::ExactMesh::Point(const std::size_t &index) const
{
    if (_points.size() < index)
    {
        throw std::runtime_error("Point index does not exist");
    }

    return _points.at(index);
}

const moist::ExactMesh::ExactCell& moist::ExactMesh::Cell(const std::size_t& index) const
{
    if (_cells.size() < index)
    {
        throw std::runtime_error("Point index does not exist");
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

    for (const auto& c : this->_cells)
    {
        if (c._deleted)
        {
            continue;
        }

        mesh.cells.create_tet(
            geo::index_t(c._points[0]),
            geo::index_t(c._points[1]),
            geo::index_t(c._points[2]),
            geo::index_t(c._points[3])
        );
    }

    moist::utils::geogram::save(file, mesh);
}
#endif // NDEBUG
