//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_VERTEX_H
#define FEA_MESHER2D_VERTEX_H

#include <array>
#include "Application.h"
#include <mpi.h>

//The Vertex class is responsible for storing information about mesh vertices.

class Vertex {
public:
    //Default constructor, sets id to -1
    Vertex();

    //Alternate constructor
    //Input:
    //point - the x and y coordinates
    //boundary - true if the vertex is on a boundary edge
    //id - unique identifier
    Vertex(std::array<double, 2> point, bool boundary, long id);

    //Destructor
    ~Vertex() = default;

    //Sets the x-coordinate to coordinate if index is 0
    //Sets the y-coordinate to coordinate if index is 1
    void setCoordinate(int index, double coordinate);

    //x-coordinate is index 0
    //y-coordinate is index 1
    double getCoordinate(int index) const;

    void setBoundary(bool b);
    long getID() const;
    bool getBoundary() const;
    std::array<double, 2> getPoint() const;

    //Used to determine if a vertex should be part of the lower convex hull when decomposing the boundary layer
    //Answers the question: which side of the directed line a->b does this vertex lie towards
    //Input:
    //a - pointer to the starting vertex of the directed line
    //b - reference to the ending vertex of the directed line
    bool orientation2D(const Vertex* a, const Vertex& b) const;

    //Calculates the project coordinates of this vertex on the paraboloid
    //Input:
    //base - the base of the paraboloid
    //axis - the non-z axis that is used in the vertical plane
    void calculateProjected(Vertex* base, bool axis);

    //Projected x-coordinate is index 0
    //Projected y-coordinate is index 1
    double getProjected(int index) const;

    //Returns true and sets to false if this vertex is on the lower convex hull when decomposing the boundary layer
    bool useLowerConvexHull();

    void setLowerConvexHull(bool exist);

    //Prints to std::cout
    void print() const;

    //Sets datatype
    //Only needs to be called once by each process
    static void createMPIDataType();

    static MPI_Datatype getMPIDataType();

private:
    //The x and y coordinates of the point
    std::array<double, 2> coordinates;
    //The coordinates of the projected point on the vertical plane for decomposing the boundary layer
    std::array<double, 2> projected;
    //Flag for if the vertex is on a boundary edge
    bool boundary;
    //Flag for if the vertex is on the lower convex hull when decomposing the boundary layer
    bool lower_convex_hull;
    //Unique identifier
    long id;
    //Defines the memory layout used for MPI communications
    static MPI_Datatype datatype;
}; //class Vertex

//Lexicographical sorting functions
bool vertexCompareX(const Vertex& a, const Vertex& b);
bool vertexCompareY(const Vertex& a, const Vertex& b);
static std::array<bool (*) (const Vertex& a, const Vertex& b), 2> vertex_compare{vertexCompareX, vertexCompareY};


#endif //FEA_MESHER2D_VERTEX_H
