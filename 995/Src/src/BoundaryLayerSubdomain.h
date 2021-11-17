//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_BOUNDARYLAYERSUBDOMAIN_H
#define FEA_MESHER2D_BOUNDARYLAYERSUBDOMAIN_H

#include <unordered_set>
#include "Application.h"
#include <mpi.h>

struct triangulateio;

class BoundaryLayerMesh;

//The BoundaryLayerSubdomain class is responsible for performing the paraboloid and lower convex hull decomposition
//steps in order to split a BoundaryLayerSubdomain into two new subdomains
//The BoundaryLayerSubdomain class is also responsible for calling Triangle to triangulate its vertices

class BoundaryLayerSubdomain {
public:
    //Default constructor
    BoundaryLayerSubdomain();

    //Destructor
    ~BoundaryLayerSubdomain() = default;

    //Decomposes this subdomain into two new subdomains
    Decomposition decompose();

    //Triangulates this subdomain and returns the output to owning_mesh
    void mesh(BoundaryLayerMesh& owning_mesh);

    void sendSubdomain(int destination);
    void recvSubdomain(int source);

    //The vertices sorted lexicographically by their x-coordinates
    std::vector<Vertex> x_vertices;

    //The vertices sorted lexicographically by their y-coordinates
    std::vector<Vertex> y_vertices;

    std::vector<std::array<long, 2>> edges;

    //Represents how many times this subdomain has been decomposed
    int decomposition_level;

    //The coordinate-axis orthogonal to the cut axis
    bool axis;

    //The memory addresses of x_vertices and y_vertices
    std::array<std::vector<Vertex>*, 2> vertices;

    //The vertices that lie on the lower convex hull of the flattened paraboloid in the vertical plane
    std::vector<Vertex> lower_convex_hull;

    //The other subdomain that is formed by a decomposition step
    std::shared_ptr<BoundaryLayerSubdomain> sub_subdomain;

    //The median vertex along the coordinate-axis specified by axis
    Vertex* median_vertex;

    //The maximum number of times a subdomain can be decomposed
    static int decomposition_threshold;

    //The maximum id for all of the boundary layer vertices
    static long max_vertex_id;

private:
    void initialize();
    void determineCut();
    bool containsInternalPoints();
    void computeMedian();
    void projectPointsOnParaboloid();
    void lowerConvexHull();
    void finalizeLowerConvexHull(std::vector<Vertex*>& lch);

    Decomposition maintainSorted();
    void maintainPrimarySorted(std::vector<Vertex>& lch, std::vector<Vertex>::iterator lch_median);
    void maintainCutSorted();
    void maintainOriginalCutSorted();
    void createSubCutSorted();

    void cleanSubdomain();
    std::vector<bool> removeOldEdges();
    std::vector<long> removePinchPoints();

    MPI_Datatype createMPIDataType();
}; //BoundaryLayerSubdomain class


#endif //FEA_MESHER2D_BOUNDARYLAYERSUBDOMAIN_H
