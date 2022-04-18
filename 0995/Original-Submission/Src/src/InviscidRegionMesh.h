//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_INVISCIDREGIONMESH_H
#define FEA_MESHER2D_INVISCIDREGIONMESH_H

#include "Vertex.h"
#include "GeoPrimitives.h"
#include "Application.h"
#include "Edge.h"
#include <queue>

class MeshGenerator;
struct triangulateio;

//The InviscidRegionMesh class is responsible for creating the isotropic nearbody and inviscid subdomains
//It also contains tunable parameters to control the function that determines the desired triangle size for a point in
//space
//These parameters are fast_growth which controls the distance where the size of triangles will grow at a more rapid
//rate and uniform which uses the size of the triangles in the nearbody region for the entirety of the inviscid region
//The last parameter is decoupling_work_threshold. If you want smaller subdomains during the initial decoupling
//procedure, then decrease this number. The smaller the subdomains, the more subdomains there are, and the more
//concurrency that can be exploited.

class InviscidRegionMesh {
    friend class MeshGenerator;
    friend class InviscidRegionSubdomain;
public:
    //No default constructor because some functions require communicating with MeshGenerator
    InviscidRegionMesh() = delete;

    //Alternate constructor
    InviscidRegionMesh(MeshGenerator& owning_mesher);

    //Destructor
    ~InviscidRegionMesh() = default;

    //Creates the initial nearbody subdomain and four inviscid subdomains that extent all the way to the farfield
    void createInitialSubdomains();

    void setNearBodyBoundingBox(std::vector<std::array<point_t, 2>>& extents);

    //A square is generated that encloses the domain and each side is chord_lengths in each direction (+x, -x, +y, -y)
    //from the center of the mesh. One chord length is the length of the input geometry in the x-direction
    void setFarFieldDistance(double chord_lengths);

    //Called by all processes so everyone has the same values of the parameters used for the triangle sizing functions
    void synchronizeInviscidParameters();

    //Decouples all subdomains larger than decoupling_work_threshold
    //Processes will send and receive subdomains from each other until everyone has some subdomains
    void decoupleInitialSubdomains();

    //Called by all other processes except the root
    void receiveInitialSubdomains();

    void useUniformTriangles();

private:
    //The MeshGenerator that this inviscid region mesh belongs to
    MeshGenerator& mesher;

    //The average area of triangles in the nearbody region
    //Used as the base size to grow from when creating a graded inviscid region, where triangles are larger the further
    //they are away from the boundary layer
    double isotropic_area;

    //True to generate uniform triangles with size isotropic_area in the entirety of the inviscid region
    bool uniform = false;

    //Number of chord lengths where triangles past this distance will grow at a faster rate
    double fast_growth = 10.0;

    //The distance in chord lengths for the domain away from the geometry that will be meshed
    double farfield;

    long next_vertex_id;
    int next_edge_id;

    //The bounding box of the nearbody region
    std::array<point_t, 4> nearbody_box;

    //The center points of each input geometry element's bounding box
    std::vector<point_t> aabb_centers;

    //Used to call Triangle, contains the same data as aabb_centers
    std::vector<double> triangle_aabb_centers;

    //The edges that make up the outer border of the domain
    //Used to apply boundary conditions for the flow solver
    std::vector<Edge> farfield_edges;

    //Stores the subdomain with the largest estimated number of triangles on top
    //Yes std::less should be used for this because std::greater would put the smallest estimated subdomain on top
    std::priority_queue<inviscid_subdomain_t, std::vector<inviscid_subdomain_t>, std::less<inviscid_subdomain_t>>
            initial_subdomains;

    //The largest an inviscid subdomain can be in terms of estimated number of triangles
    const int decoupling_work_threshold = 25000;

    //-----------------------Private Functions-------------------------
    void initializeInviscidParamters();
    void createNearBodySubdomain(const std::array<std::vector<Vertex>, 4>& sides);
    std::vector<Vertex> createAxisAlignedDecouplingPath(Vertex first, Axis axis, double last);
    std::vector<Vertex> createDecouplingPath(Vertex& first, point_t last);
    void createFarfieldSubdomains(const std::array<std::vector<Vertex>, 4>& sides);
    std::tuple<bool, std::array<inviscid_subdomain_t, 4>> splitInviscidSubdomain(inviscid_subdomain_t& subdomain);
    std::tuple<bool, std::array<std::vector<Vertex>::iterator, 4>> determineSplitVertices(
            inviscid_subdomain_t& subdomain, const point_t& joint);
    std::array<inviscid_subdomain_t, 4> partitionParentVertices(
            const inviscid_subdomain_t& parent, const std::array<std::vector<Vertex>::iterator, 4>& spliterators);
    std::array<std::vector<Vertex>, 4> createSplitBorders(
            const std::array<std::vector<Vertex>::iterator, 4>& splits, const point_t& sub_center);

    double getEdgeConstraintAtPoint(const point_t& point);
    double getAreaConstraintAtPoint(const point_t& point);

    void enlargeNearBodyBoundingBox();
    void translateNearBodyBoundingBox();
    std::array<std::vector<Vertex>, 4> decoupleNearBodyBoundingBox();
    long addFarfieldEdges(long first_id, std::vector<Vertex>::iterator next, std::vector<Vertex>::iterator stop);
    const std::vector<double>& getBoundaryLayerHoles() const;

    void globalizeTriangleOutput(struct triangulateio& out, std::vector<long>& global_ids);

    void computeSubdomainCost(inviscid_subdomain_t& s);
    int approximateTriangleCount(const point_t& a, const point_t& b, double area);

    void distributeInitialUniformSubdomains();
    void distributeInitialGradedSubdomains();
    void overDecoupleLocalSubdomains();
}; //InviscidRegionMesh class


#endif //FEA_MESHER2D_INVISCIDREGIONMESH_H
