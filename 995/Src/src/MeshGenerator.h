//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_MESHGENERATOR_H
#define FEA_MESHER2D_MESHGENERATOR_H

#include <iostream>
#include <array>
#include <vector>
#include <unordered_map>
#include "BoundaryLayerMesh.h"
#include "InviscidRegionMesh.h"
#include "TriangleData.h"
#include "MeshingManager.h"
#include <mpi.h>
#include <pthread.h>
#include <queue>

struct triangulateio;

//The MeshGenerator class is what the end user will interact with in order to generate a mesh
//The MeshGenerator reads in an input geometry file, has its meshing parameters initialized: boundary layer growth
//funcdtion and inviscid region farfield distance

class MeshGenerator {
    friend class BoundaryLayerMesh;
    friend class InviscidRegionMesh;
    friend class MeshingManager;
public:

    //Default constructor
    MeshGenerator();

    //Destructor
    ~MeshGenerator() = default;

    //Reads the input model referenced by the path filename
    void readInputModel(std::string filename);

    //Sets the growth function used to generate the anisotropic boundary layer
    //Input:
    //first_thickness - the distance from the geometry that the points in the first layer will be placed
    //growth_rate - the geometric rate at which each layer's thickness will grow from the previous layer
    //              A value of 1 causes all layers to have the same thickness as first_thickness
    //              A value greater than 1 causes subsequent layers to grow in thickness, which is the intended usage
    //initial_layers - the number of anisotropic layers that will be originally grown from the geometry
    //                 Some triangles in some layers will be removed due to intersections or a poor quality shape
    void setBoundaryLayerGrowthFunction(double first_thickness, double growth_rate, int initial_layers);

    //Sets the distance for the domain away from the geometry that will be meshed
    //Input:
    //chord_lengths - the distance in each direction (+x, -x, +y, -y) from the center of the mesh
    //                One chord length is the length of the input geometry in the x-direction
    void setFarFieldDistance(double chord_lengths);

    //Uses the same sized triangles to mesh the entirety of the inviscid region
    //By default, the inviscid region uses smaller triangles near the geometry and grades to larger triangles farther
    //away from the geometry
    void useUniformTrianglesForInviscidRegion();

    //This function creates all the mesh vertices and mesh triangles using the number of processes launched by MPI
    void meshDomain();

    //Collects all of the mesh entities on the root process so that the mesh can be output
    void collectFinalMesh();

    //Outputs the final mesh as a Tecplot file
    void outputTecplot() const;

    //Outputs the final mesh files for Shewchuk's "Show Me" application
    void outputShowMe() const;

    //Outputs the final mesh as a .msh file used for the flow solver FUN3D
    void outputFUN3D() const;

    //Outputs the final mesh in the .vtk format, suitable for ParaView
    void outputVTK() const;

private:
    //Unique id for this process
    const int rank;

    //Number of distributed MeshGenerators working
    const int processes;

    BoundaryLayerMesh boundary_layer;
    InviscidRegionMesh inviscid_region;
    MeshingManager manager;

    //The name of the input file
    //Used to name the output files
    std::string input_name;

    //The x and y coordinates of the resulting mesh vertices
    std::vector<std::array<double, 2>> final_vertices;

    //Number of edges on the input geometry plus the edges on the farfield
    //Used to apply boundary conditions
    int final_edges;

    //The endpoint vertex ids of the resulting mesh triangles
    std::vector<std::array<long, 3>> final_triangles;

    //Maps the non-sequential global id of a mesh vertex to its id in the resulting mesh
    std::unordered_map<long, long> global_to_final;

    //Unique sequential identifier to assign to the next final vertex that will be registered
    long next_final_id = 0;

    //Number of triangles in each resulting subdomain mesh
    std::vector<long> subdomain_num_triangles;

    std::vector<bl_subdomain_t> boundary_layer_subdomains;

    //Used to synchronize access to inviscid_subdomains between the manager and worker thread
    pthread_mutex_t subdomain_mutex;

    //Holds all of the inviscid subdomains for a process
    //std::less makes it so the subdomain at the top of the queue is the most expensive
    std::priority_queue<inviscid_subdomain_t, std::vector<inviscid_subdomain_t>, std::less<inviscid_subdomain_t>>
            inviscid_subdomains;

    //Holds the triangulated or refined subdomain meshes
    std::queue<TriangleData> outbox;

    //Flag that denotes if the worker thread is currently triangulating or refining a subdomain
    bool meshing = false;

    //-----------------------Private Functions-------------------------
    bool isIdle();

    void addMeshToOutbox(struct triangulateio& out, std::vector<long>& global_ids);
    void processOutbox();
    void addToFinalMesh(double* vertices, int num_triangles, int* triangles, std::vector<long>& global_ids);
    void sendMeshToRoot(double* vertices, int num_triangles, int* triangles, std::vector<long>& global_ids);
    void recvMesh(int source);

    std::tuple<inviscid_subdomain_t, bool> getSubdomainFromQueue();
    void addSubdomainToQueue(inviscid_subdomain_t& subdomain);
    void meshSubdomains();
}; //MeshGenerator class


#endif //FEA_MESHER2D_MESHGENERATOR_H
