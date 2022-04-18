//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_BOUNDARYLAYERMESH_H
#define FEA_MESHER2D_BOUNDARYLAYERMESH_H

#include <iostream>
#include <vector>
#include "Application.h"
#include "Vertex.h"
#include "Edge.h"
#include "Ray.h"
#include "GeoPrimitives.h"
#include "BoundaryLayerSubdomain.h"

class MeshGenerator;
class InviscidRegionMesh;

//The BoundaryLayerMesh class is responsible for generating the high-fidelity anisotropic boundary layer points, edges,
//and initial subdomains that will be triangulated. This class checks for intersections in the boundary layer and
//resolves them, smooths poor quality ray regions, adds points to ensure a smooth transition to the inviscid region, and
//removes points that would result in poor quality triangles

class BoundaryLayerMesh {
    friend class MeshGenerator;
    friend class Ray;
    friend class InviscidRegionMesh;
    friend class BoundaryLayerSubdomain;
public:
    //No default constructor because a MeshGenerator is needed
    BoundaryLayerMesh() = delete;

    //Alternate constructor
    BoundaryLayerMesh(MeshGenerator& owning_mesher);

    //Destructor
    ~BoundaryLayerMesh() = default;

    //Reads the input model referenced by the path filename
    void initializeModelSurface(std::string filename);

    //Sets the growth function used to generate the anisotropic triangles
    //Input:
    //first_layer_thickness - the distance from the geometry that the points in the first layer will be placed
    //layer_growth_rate - the geometric rate at which each layer's thickness will grow from the previous layer
    //                    A value of 1 causes all layers to have the same thickness as first_layer_thickness
    //                    A value greater than 1 causes subsequent layers to grow in thickness
    //initial_layers - the number of anisotropic layers that will be originally grown from the geometry
    //                 Some triangles in some layers will be removed due to intersections or a poor quality shape
    void setGrowthFunction(double first_layer_thickness, double layer_growth_rate, int initial_layers);

    void insertBoundaryLayerPoints();

    //Creates the initial boundary layer subdomains
    void createBoundaryLayerSubdomains();

    //Decomposes the initial subdomains of the boundary layer until each process has a subdomain or a subdomain
    //cannot be further decomposed
    void decomposeInitialSubdomains();

    void receiveInitialSubdomains();

    void outputVTK() const;

private:
    //The owning MeshGenerator
    MeshGenerator& mesher;

    //The neighboring inviscid region
    InviscidRegionMesh* inviscid_region;

    int num_model_vertices;
    int num_model_edges;
    int num_elements;

    //The distance along the x-axis of the input model
    double chord_length;

    point_t mesh_center;

    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<double> holes;

    //The id of the first vertex of each element
    std::vector<int> surface_element_start;

    //The vertices that are on the enclosing border of the boundary layer
    std::vector<std::vector<Vertex*>> transition_vertices;
    std::vector<Ray> rays;

    //The index of the first ray of each element
    std::vector<int> ray_element_start;
    int num_enclosing_edges;

    //If the angle between two rays is greater than this value, then those rays will be added to large_angles
    double ray_angle_tolerance = 4.0;

    //Ray angles larger than this value, but less than trailing_edge_angle_tolerance will be marked as a cusp
    double cusp_angle_tolerance = 30.0;

    //Ray angles larger than this value will be marked as being a sharp trailing edge
    double trailing_edge_angle_tolerance = 60.0;

    std::vector<int> cusps;
    std::vector<int> sharp_trailing_edges;
    std::vector<std::tuple<int, int, double>> large_angles;

    //A process' subset of the boundary layer vertices that it created
    std::vector<Vertex> local_vertices;

    long next_vertex_id;
    int next_edge_id;
    long max_boundary_layer_vertex_id;

    //The id of the first vertex of each element that is not on the model's surface
    std::vector<long> boundary_layer_element_start;

    //The thickness of the final layer of the boundary layer
    double last_thickness;

    double growth_rate;
    int num_layers;

    //The maximum number of layers that can exist in the boundary layer
    //Used for gradation control
    //Set to 1.25 * num_layers
    int max_layers;

    //The thicknesses of each layer
    std::vector<double> layer_offsets;

    //The extent box that contains the entire boundary layer
    AABB max_extent_box;

    //The largest extent box of each element
    std::vector<AABB> max_element_extent_boxes;

    std::vector<bl_subdomain_t> initial_subdomains;

    //Used by the root process to distribute initial subdomains
    int next_recv_process = 0;

    int total_initial_subdomains;

    //-----------------------Private Functions-------------------------
    void checkElementsWindingOrder();

    const Vertex& getPreviousSurfaceVertex(int index) const;
    const Vertex& getNextSurfaceVertex(int index) const;
    int getElement(long vertex_id) const;
    int getPreviousRayId(int index) const;
    int getNextRayId(int index) const;
    Ray& getPreviousRay(int index);
    const Ray& getPreviousRay(int index) const;
    Ray& getNextRay(int index);
    const Ray& getNextRay(int index) const;

    std::vector<Ray> initializeLocalRays();

    void refineBoundaryLayer();
    void detectCuspsAndLargeAngles();
    void determineVertexInsertID(const std::vector<Ray>& local_rays);
    void determineLocalBoundaryLayerSize(const std::vector<Ray>& local_rays);
    void createBoundaryLayerVertices(std::vector<Ray>& local_rays);
    void gatherLocalVerticesOnRoot();
    void closeBoundaryLayerExterior(std::vector<Ray>& local_rays);
    void createEnclosingEdges();
    void createEdgesAlongRay(const Ray& ray);

    double getAverageEnclosingTriangleArea() const;
    const std::vector<double>& getModelHoles() const;

    std::vector<bl_subdomain_t> createSubdomainsForElement(int element);
    std::vector<int> getElementRaySplits(int element);
    std::vector<bool> populateSubdomainVertices(bl_subdomain_t& subdomain, int ray_index, int last_ray);
    void determineSubdomainEdges(bl_subdomain_t& subdomain, const std::vector<bool>& exist);
    std::array<point_t, 2> finalizeSubdomains(std::vector<bl_subdomain_t>& subdomains);

    void addToMesherOutbox(struct triangulateio& out, std::vector<long> &global_ids);

    void gradationReductionControl();
    void removeRaySpikes();
    void gradationGrowthControl();

    void handleSelfIntersections();

    std::vector<point_t> computeBoundaryPoints(int element) const;
    std::vector<point_t> computeExtendedBoundaryPoints(int element) const;
    std::vector<Segment> createRaySegments(const std::vector<point_t>& end_points, int element) const;
    std::vector<Segment> createSurfaceSegments(int element) const;
    std::vector<Segment> createOuterBorderSegments(const std::vector<Segment>& ray_segments);
    std::vector<Segment> createOuterBorderSegments(const std::vector<point_t>& boundary_points);
    void setMaximumBoundaryLayerDomain(int element);
    void resolveSurfaceIntersections(int element, std::vector<Segment>& test_segments);
    void resolveOuterBorderIntersections(int element, std::vector<Segment>& test_segments);

    std::vector<std::array<int, 2>> findExtentIntersections(
            const std::vector<Segment>& test_segments, const std::vector<Segment>& target_segments, int test_element,
            int target_element) const;
    std::vector<std::array<int, 2>> findExtentIntersections(const std::vector<Segment>& test_segments,
                                                            int element) const;
    std::vector<std::tuple<int, int, point_t>> findIntersections(
            const std::vector<Segment>& test_segments, const std::vector<Segment>& target_segments,
            const std::vector<std::array<int, 2>>& candidates) const;
    void clipSurfaceIntersections(const std::vector<std::tuple<int, int, point_t>>& intersections,
                                  std::vector<Segment>& test_segments);
    void clipOuterBorderIntersections(const std::vector<std::tuple<int, int, point_t>>& intersections,
                                      std::vector<Segment>& test_segments, std::vector<Segment>& border, int element);
    void prioritizeSelfIntersections(std::vector<std::tuple<int, int, point_t>>& intersections, int element);
    void smoothPoorQualityRays(int element, std::vector<Segment>& test_segments);

    std::vector<std::array<int, 2>> determineSmoothingRegions(
            const std::vector<std::tuple<int, int, point_t>>& intersections, int element);
    void setNeighborhoodsForContiguousIntersections(std::vector<bool>& crosses, int element);
    void setNeighborhoodsForCuspsAndLargeAngles(std::vector<bool>& crosses, int element);
    std::vector<std::array<int, 2>> determineRayRangesToSmooth(const std::vector<bool>& crosses, int element);
    std::vector<int> getFixedRays(int start, int stop);
    void raySmoothing(int start, int stop);

    void carveOutElementCavities();
    std::vector<std::array<int, 2>> detectElementCavities(int element);
    void smoothCavities(const std::vector<std::array<int, 2>>& cavities);

    void handleMultiElementIntersections();
    std::vector<int> pruneRayAABBIntersections(const std::vector<Segment>& test_segments, const AABB& target) const;
    std::vector<std::array<int, 2>> findExtentIntersections(
            const std::vector<Segment>& test_segments, const std::vector<int>& candidates,
            const std::vector<Segment>& target_segments) const;
    void clipMultiElementIntersections(const std::vector<std::tuple<int, int, point_t>>& intersections,
                                       std::vector<Segment>& test_segments, int test_element,
                                       std::vector<Segment>& target_border, int target_element);

    void addFansAtSharpTrailingEdges();
    std::vector<Ray> createHalfFan(const Ray* a, const Ray* b, int size, const Vertex& origin, bool reverse);
    void incrementElementRayStarts(std::vector<std::array<int, 3>>& ray_increments);
}; //BoundaryLayerMesh class


#endif //FEA_MESHER2D_BOUNDARYLAYERMESH_H
