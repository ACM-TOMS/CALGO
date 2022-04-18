//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_EDGE_H
#define FEA_MESHER2D_EDGE_H
#include "Application.h"
#include "Vertex.h"
#include <mpi.h>

//The Edge class is responsible for storing information about mesh edges.

class Edge {
public:
    //Default constructor, sets everything to -1
    Edge();

    //Alternate constructor
    //Input:
    //edge_id - unique identifier
    //vertex_a - reference to starting endpoint
    //vertex_b - reference to ending endpoint
    //edge_type - classification of the type of edge
    Edge(int edge_id, const Vertex& vertex_a, const Vertex& vertex_b, int edge_type);

    //Alternate constructor
    //Input
    //edge_id - unique identifier
    //vertex_a - id of starting endpoint
    //vertex_b - id of ending endpoint
    //edge_type - classification of the type of edge
    Edge(int edge_id, long vertex_a, long vertex_b, int edge_type);

    ~Edge() = default;

    //Sets the endpoint ids
    //Input:
    //vertex_a - reference to starting endpoint
    //vertex_b - reference to ending endpoint
    void setVertices(Vertex& vertex_a, Vertex& vertex_b);

    void setType(int edge_type);
    int getId() const;
    long getVertexId(int index) const;
    int getType() const;
    std::array<long, 2> getVertices() const;

    //Returns true if the edge is part of the geometry or farfield boundary
    //Used to apply boundary conditions
    bool shouldBeInFinalMesh() const;

    void print() const;

    //Defines the memory layout used for MPI communications
    static MPI_Datatype createMPIDatatype();

private:
    //Unique identifier
    int id;

    //Identifiers of its endpoints
    std::array<long, 2> vertices;

    //Classification
    //0 for constrained
    //1 for geometry
    //2 for farfield boundary
    //3 for boundary layer outer border
    int type;
}; //class Edge


#endif //FEA_MESHER2D_EDGE_H
