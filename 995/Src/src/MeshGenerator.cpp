//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "MeshGenerator.h"
#include "MPICommunications.h"
#include "BoundaryLayerSubdomain.h"
#include "InviscidRegionSubdomain.h"
#include <fstream>
#include <unistd.h>
#include <numeric>

extern "C" {
#include "triangle/triangle.h"
}

using namespace std;
using namespace Application;
using namespace MPICommunications;

MeshGenerator::MeshGenerator()
        : rank(myRank()), processes(numberOfProcesses()), boundary_layer(*this), inviscid_region(*this), manager(*this) {
    Vertex::createMPIDataType();
    pthread_mutex_init(&subdomain_mutex, NULL);
}

void MeshGenerator::readInputModel(std::string filename) {
    input_name = filename;
    string::size_type index = input_name.find(".poly");
    if(index != string::npos)
        input_name.erase(index, string(".poly").length());
    index = input_name.find("/");
    if(index != string::npos)
        input_name.erase(0, index+1);

    boundary_layer.initializeModelSurface(filename);
}

void MeshGenerator::setBoundaryLayerGrowthFunction(double first_thickness, double growth_rate, int initial_layers) {
    boundary_layer.setGrowthFunction(first_thickness, growth_rate, initial_layers);
}

void MeshGenerator::setFarFieldDistance(double chord_lengths) {
    inviscid_region.setFarFieldDistance(chord_lengths);
}

void MeshGenerator::useUniformTrianglesForInviscidRegion() {
    inviscid_region.useUniformTriangles();
}

void MeshGenerator::meshDomain() {
    boundary_layer.insertBoundaryLayerPoints();

    if(rank == 0) {
        boundary_layer.createBoundaryLayerSubdomains();
        inviscid_region.createInitialSubdomains();
    } else {
        boundary_layer.receiveInitialSubdomains();
        inviscid_region.receiveInitialSubdomains();
    }

    boundary_layer.decomposeInitialSubdomains();
    inviscid_region.synchronizeInviscidParameters();
    inviscid_region.decoupleInitialSubdomains();

    meshing = true;

    pthread_t manager_thread;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int ls = pthread_create(&manager_thread, &attr, [](void* mesher)->void* {
        ((MeshGenerator*) mesher)->manager.manageMeshingProgress();
        return 0;
    }, &this[0]);
    if(ls != 0) cerr << "Launch status error " << ls << " for communicator thread " << endl;

    meshSubdomains();

    pthread_join(manager_thread, NULL);
}

void MeshGenerator::meshSubdomains() {
    inviscid_subdomain_t subdomain;

    while(not boundary_layer_subdomains.empty()) {
        boundary_layer_subdomains.back()->mesh(boundary_layer);
        boundary_layer_subdomains.pop_back();
    }

    meshing = false;

    while(not manager.finished) {
        tie(subdomain, meshing) = getSubdomainFromQueue();

        if(meshing) {
            subdomain->mesh(inviscid_region);
            subdomain.reset();
            meshing = false;
        } else {

        }
    }
}

bool MeshGenerator::isIdle() {
    pthread_mutex_lock(&subdomain_mutex);
    bool idle = boundary_layer_subdomains.empty() and inviscid_subdomains.empty() and not meshing;
    pthread_mutex_unlock(&subdomain_mutex);
    return idle;
}

void MeshGenerator::addMeshToOutbox(struct triangulateio& out, std::vector<long>& global_ids) {
    outbox.emplace(out.pointlist, out.numberoftriangles, out.trianglelist, global_ids);
}

std::tuple<inviscid_subdomain_t, bool> MeshGenerator::getSubdomainFromQueue() {
    inviscid_subdomain_t subdomain;
    bool has_subdomain = false;

    pthread_mutex_lock(&subdomain_mutex);
    if(not inviscid_subdomains.empty()) {
        subdomain = move(inviscid_subdomains.top());
        inviscid_subdomains.pop();
        manager.work_units -= subdomain->cost;
        has_subdomain = true;
    }
    pthread_mutex_unlock(&subdomain_mutex);

    return make_tuple(subdomain, has_subdomain);
}

void MeshGenerator::addSubdomainToQueue(inviscid_subdomain_t& subdomain) {
    manager.work_units += subdomain->cost;
    pthread_mutex_lock(&subdomain_mutex);
    inviscid_subdomains.push(move(subdomain));
    pthread_mutex_unlock(&subdomain_mutex);
}

void MeshGenerator::collectFinalMesh() {
    vector<int> outbox_sizes;
    Gather(((int) outbox.size()), outbox_sizes, 0);

    processOutbox();

    if(rank == 0)
        for(int p=1; p<processes; ++p)
            for(int i=0; i<outbox_sizes[p]; ++i)
                recvMesh(p);
}

void MeshGenerator::processOutbox() {
    while(not outbox.empty()) {
        TriangleData data(move(outbox.front()));
        outbox.pop();
        if(rank == 0)
            addToFinalMesh(data.vertices, data.num_triangles, data.triangles, data.global_ids);
        else sendMeshToRoot(data.vertices, data.num_triangles, data.triangles, data.global_ids);
    }
}

void MeshGenerator::addToFinalMesh(double* vertices, int num_triangles, int* triangles, std::vector<long>& global_ids) {
    for(int i=0; i<global_ids.size(); ++i) {
        if(not global_to_final.count(global_ids[i])) {
            global_to_final[global_ids[i]] = next_final_id++;
            final_vertices.push_back({{vertices[2*i], vertices[(2*i)+1]}});
        }
    }

    bool valid;
    for(int t=0; t<num_triangles; ++t) {
        valid = true;
        for(int i=0; i<3; ++i)
            if(not global_to_final.count(global_ids[triangles[i+(t*3)]])) {
                cout << "ERROR: global id not registered for local node " << global_ids[triangles[i+(t*3)]] << endl;
                valid = false;
            }

        if(valid)
            final_triangles.push_back({{global_to_final[global_ids[triangles[t*3]]],
                                               global_to_final[global_ids[triangles[(t*3)+1]]],
                                               global_to_final[global_ids[triangles[(t*3)+2]]]}});
        else final_triangles.push_back({{0, 0, 0}});
    }

    subdomain_num_triangles.push_back(num_triangles);
}

void MeshGenerator::sendMeshToRoot(double* vertices, int num_triangles, int* triangles, std::vector<long>& global_ids) {
    array<int, 2> sizes {(int)global_ids.size(), num_triangles};
    MPI_Send(sizes.data(), 2, MPI_INT, 0, MessageTag::TRANSFER_MESH, MPI_COMM_WORLD);
    MPI_Send(vertices, (int)(2*global_ids.size()), MPI_DOUBLE, 0, MessageTag::FINAL_VERTICES, MPI_COMM_WORLD);
    MPI_Send(triangles, 3*num_triangles, MPI_INT, 0, MessageTag::FINAL_TRIANGLES, MPI_COMM_WORLD);
    MPI_Send(global_ids.data(), (int)global_ids.size(), MPI_LONG, 0, MessageTag::GLOBAL_IDS, MPI_COMM_WORLD);
}

void MeshGenerator::recvMesh(int source) {
    array<int, 2> sizes = array<int, 2> ();
    MPI_Recv(sizes.data(), 2, MPI_INT, source, MessageTag::TRANSFER_MESH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double* vertices = new double[2*sizes[0]];
    int* triangles = new int[3*sizes[1]];
    vector<long> global_ids((unsigned long)sizes[0], 0);

    MPI_Recv(vertices, 2*sizes[0], MPI_DOUBLE, source, MessageTag::FINAL_VERTICES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(triangles, 3*sizes[1], MPI_INT, source, MessageTag::FINAL_TRIANGLES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(global_ids.data(), sizes[0], MPI_LONG, source, MessageTag::GLOBAL_IDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    addToFinalMesh(vertices, sizes[1], triangles, global_ids);
}

void MeshGenerator::outputTecplot() const {
    ofstream out("output/" + input_name + ".dat");

    out << "VARIABLES = \"X\" \"Y\"" << endl;
    out << "ZONE NODES=" << final_vertices.size() << ", ELEMENTS=" << final_triangles.size()
    << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE" << endl;

    out.precision(16);

    for(long v=0; v<final_vertices.size(); v++)
        out << final_vertices[v][0] << " " << final_vertices[v][1] << endl;

    for(long t=0; t<final_triangles.size(); t++)
        out << final_triangles[t][0]+1 << " " << final_triangles[t][1]+1 << " " << final_triangles[t][2]+1 << endl;

    out.close();
}

void MeshGenerator::outputShowMe() const {
    string filename = "output/" + input_name;
    ofstream node((filename + ".node"));
    ofstream ele((filename + ".ele"));
    ofstream part((filename + ".part"));

    node.precision(16);

    node << final_vertices.size() << " 2 0 1" << endl;
    for(int v=0; v<final_vertices.size(); ++v)
        node << v << " " << final_vertices[v][0] << " " << final_vertices[v][1] << " 1" << endl;
    node.close();

    ele << final_triangles.size() << " 3 0" << endl;
    for(int t=0; t<final_triangles.size(); ++t)
        ele << t << " " << final_triangles[t][0] << " " << final_triangles[t][1] << " " << final_triangles[t][2] << endl;
    ele.close();

    if(not subdomain_num_triangles.empty()) {
        part << final_triangles.size() << " " << subdomain_num_triangles.size() << endl;
        long id = 0;
        for(long i=0; i<subdomain_num_triangles.size(); ++i)
            for(long j=0; j<subdomain_num_triangles[i]; ++j)
                part << id++ << " " << i << endl;
    }
    part.close();
}

void MeshGenerator::outputFUN3D() const {
    ofstream out("output/" + input_name + ".msh");

    out << "MeshVersionFormatted 0" << endl << endl;
    out << "Dimension" << endl << "2" << endl << endl;

    out.precision(16);

    out << "Vertices" << endl << final_vertices.size() << endl;
    for(long v=0; v<final_vertices.size(); v++)
        out << final_vertices[v][0] << " " << final_vertices[v][1] << " 1" << endl;

    out << endl;

    out << "Edges" << endl << final_edges << endl;
    for(const Edge& edge : boundary_layer.edges)
        if(edge.shouldBeInFinalMesh())
            out << global_to_final.at(edge.getVertexId(0)) + 1 << " " << global_to_final.at(edge.getVertexId(1)) + 1
                << 3 << endl;
    for(const Edge& edge : inviscid_region.farfield_edges)
        if(edge.shouldBeInFinalMesh())
            out << global_to_final.at(edge.getVertexId(0)) + 1 << " " << global_to_final.at(edge.getVertexId(1)) + 1
                << 4 << endl;
    out << endl;

    out << "Triangles" << endl << final_triangles.size() << endl;
    for(long t=0; t<final_triangles.size(); t++)
        out << final_triangles[t][0]+1 << " " << final_triangles[t][1]+1 << " " << final_triangles[t][2]+1 << " 1" << endl;

    out << endl << "End" << endl;

    out.close();
}

void MeshGenerator::outputVTK() const {
    ofstream out("output/" + input_name + ".mesh.vtk");
    out << "# vtk DataFile Version 3.0" << endl;
    out << input_name << " Mesh" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;

    out.precision(16);

    out << "POINTS " << final_vertices.size() << " double" << endl;
    for(const auto& vertex : final_vertices)
        out << vertex[0] << " " << vertex[1] << " 0" << endl;

    out << endl << "CELLS " << final_edges + final_triangles.size() << " " << 3*final_edges + 4*final_triangles.size()
    << endl;

    for(const Edge& edge : boundary_layer.edges)
        if(edge.shouldBeInFinalMesh())
            out << "2 " << global_to_final.at(edge.getVertexId(0)) << " " << global_to_final.at(edge.getVertexId(1)) << endl;
    for(const Edge& edge : inviscid_region.farfield_edges)
        if(edge.shouldBeInFinalMesh())
            out << "2 " << global_to_final.at(edge.getVertexId(0)) << " " << global_to_final.at(edge.getVertexId(1)) << endl;

    for(const auto& triangle : final_triangles)
        out << "3 " << triangle[0] << " " << triangle[1] << " " << triangle[2] << endl;

    out << endl << "CELL_TYPES " << final_edges + final_triangles.size() << endl;
    for(long i=0; i<final_edges; ++i)
        out << "3" << endl;
    for(long i=0; i<final_triangles.size(); ++i)
        out << "5" << endl;

    out << endl;
    out.close();
}
