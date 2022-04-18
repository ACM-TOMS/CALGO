#include "Vertex.h"
#include <iostream>
#include <cstddef>

using namespace Application;

MPI_Datatype Vertex::datatype = -1;

Vertex::Vertex() : id(-1) { }

Vertex::Vertex(std::array<double, 2> point, bool boundary, long id)
        : coordinates(point), boundary(boundary), id(id) {}

void Vertex::setCoordinate(int index, double coordinate) {
    coordinates[index] = coordinate;
}

void Vertex::setBoundary(bool b) {
    boundary = b;
}

long Vertex::getID() const {
    return id;
}

double Vertex::getCoordinate(int index) const {
    return coordinates[index];
}

bool Vertex::getBoundary() const {
    return boundary;
}

point_t Vertex::getPoint() const {
    return coordinates;
}

bool Vertex::orientation2D(const Vertex* a, const Vertex& b) const {
    return (((a->projected[0] - projected[0]) * (b.projected[1] - projected[1])) -
            ((a->projected[1] - projected[1]) * (b.projected[0] - projected[0]))) < -Application::precision;
}

void Vertex::calculateProjected(Vertex* base, bool axis) {
    projected[0] = coordinates[!axis] - base->coordinates[!axis];
    projected[1] = pow((coordinates[0]-base->coordinates[0]), 2) + pow((coordinates[1]-base->coordinates[1]), 2);
}

double Vertex::getProjected(int index) const {
    return projected[index];
}

bool Vertex::useLowerConvexHull() {
    if(lower_convex_hull) {
        lower_convex_hull = false;
        return true;
    }

    return false;
}

void Vertex::setLowerConvexHull(bool exist) {
    lower_convex_hull = exist;
}

void Vertex::print() const {
    std::cout << "Vertex " << id << ": (" << coordinates[0] << ", " << coordinates[1] << ") is ";
    if(not boundary)
        std::cout << "not ";
    std::cout << "on boundary";
}

MPI_Datatype Vertex::getMPIDataType() {
    return Vertex::datatype;
}

bool vertexCompareX(const Vertex& a, const Vertex& b) {
    return (a.getCoordinate(0) < (b.getCoordinate(0) - precision)) or
           (areEqual(a.getCoordinate(0), b.getCoordinate(0)) and (a.getCoordinate(1) < b.getCoordinate(1)));
}

bool vertexCompareY(const Vertex& a, const Vertex& b) {
    return (a.getCoordinate(1) < (b.getCoordinate(1)) - precision) or
           (areEqual(a.getCoordinate(1), b.getCoordinate(1)) and (a.getCoordinate(0) < b.getCoordinate(0)));
}

void Vertex::createMPIDataType() {
    const int items = 3;
    int block_lengths[items] = {2, 1, 1};
    MPI_Datatype types[items] = {MPI_DOUBLE, MPI_CXX_BOOL, MPI_LONG};
    MPI_Aint offsets[items] = {offsetof(Vertex, coordinates), offsetof(Vertex, boundary), offsetof(Vertex, id)};
    MPI_Datatype mpi_datatype;
    MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_datatype);
    MPI_Type_commit(&mpi_datatype);
    Vertex::datatype = mpi_datatype;
}
