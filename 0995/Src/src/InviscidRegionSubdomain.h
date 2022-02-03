//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_INVISCIDREGIONSUBDOMAIN_H
#define FEA_MESHER2D_INVISCIDREGIONSUBDOMAIN_H

#include <unordered_set>
#include "Application.h"
#include <mpi.h>

struct triangulateio;
class InviscidRegionMesh;

//The InviscidRegionSubdomain class is responsible for storing the information to refine an isotropic subdomain
//Vertices are wound counter-clockwise and form a simple, closed polygon

class InviscidRegionSubdomain {
public:
    //Default constructor, sets nearbody to false
    InviscidRegionSubdomain();

    //Alternate constructor
    InviscidRegionSubdomain(bool nearbody_subdomain);

    //Destructor
    ~InviscidRegionSubdomain() = default;

    //Computes the bounding box, used to compute the estimated number of triangles in the subdomain and to determine the
    //center point
    std::array<point_t, 2> getExtentBox() const;

    //Gets the center point of this subdomain's bounding box
    point_t getCenterPoint() const;

    //Creates the counter-clockwise edges
    void createBorder();

    //Refines this subdomain and passes the resulting mesh back to owning_mesh
    //This subdomain can be destroyed after this function returns
    void mesh(InviscidRegionMesh& owning_mesh);

    void sendSubdomain(int destination);
    void recvSubdomain(int source);

    //Defines the memory layout used for MPI communications
    MPI_Datatype createMPIDataType();

    //Wound counter-clockwise
    std::vector<Vertex> vertices;
    std::vector<std::array<long, 2>> edges;

    //True if the subdomain contains a portion of the input geometry
    bool nearbody;

    //Estimated number of triangles that will be in this subdomain
    int cost = -1;
}; //InviscidRegionSubdomain class

//Comparison function used for sorting the subdomains in a priority queue
inline bool operator<(const inviscid_subdomain_t& subdomain_a, const inviscid_subdomain_t& subdomain_b) {
    return subdomain_a->cost < subdomain_b->cost;
}


#endif //FEA_MESHER2D_INVISCIDREGIONSUBDOMAIN_H
