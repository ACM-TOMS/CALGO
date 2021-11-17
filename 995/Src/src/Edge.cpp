//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include <iostream>
#include <cstddef>
#include "Edge.h"

using namespace std;

Edge::Edge() : id(-1), vertices(array<long, 2>{{-1, -1}}), type(-1) {}

Edge::Edge(int edge_id, const Vertex& vertex_a, const Vertex& vertex_b, int edge_type)
        : id(edge_id), vertices(array<long, 2>{{vertex_a.getID(), vertex_b.getID()}}), type(edge_type) {}

Edge::Edge(int edge_id, long vertex_a, long vertex_b, int edge_type)
        : id(edge_id), vertices(array<long, 2> {vertex_a, vertex_b}), type(edge_type) {}

void Edge::setVertices(Vertex& vertex_a, Vertex& vertex_b) {
    vertices = {vertex_a.getID(), vertex_b.getID()};
}

void Edge::setType(int edge_type) {
    type = edge_type;
}

int Edge::getId() const {
    return id;
}

long Edge::getVertexId(int index) const {
    return vertices[index];
}

int Edge::getType() const {
    return type;
}

std::array<long, 2> Edge::getVertices() const {
    return vertices;
}

void Edge::print() const {
    cout << "Edge " << id << ": Comprised of vertices (" << vertices[0] << " & " << vertices[1] << ") and is ";

    if(type == 0) cout << "constrained";
    else if(type == 1) cout << "part of the mesh geometry";
    else if(type == 2) cout << "part of the outer boundary";
    else if(type == 3) cout << "part of the enclosing boundary layer";
    else throw logic_error("Invalid edge type = " + to_string(type));

    cout << endl;
}

MPI_Datatype Edge::createMPIDatatype() {
    const int items = 3;
    int block_lengths[items] = {1, 2, 1};
    MPI_Datatype types[items] = {MPI_INT, MPI_LONG, MPI_INT};
    MPI_Aint offsets[items] = {offsetof(Edge, id), offsetof(Edge, vertices), offsetof(Edge, type)};
    MPI_Datatype mpi_datatype;
    MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_datatype);
    MPI_Type_commit(&mpi_datatype);
    return mpi_datatype;
}

bool Edge::shouldBeInFinalMesh() const {
    return type == 1 or type == 2;
}
