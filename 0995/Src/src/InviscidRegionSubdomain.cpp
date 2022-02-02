//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "InviscidRegionSubdomain.h"
#include "InviscidRegionMesh.h"
#include "Vertex.h"
#include "MPICommunications.h"
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <cstddef>

extern "C" {
#include "triangle/triangle.h"
}

using namespace std;
using namespace Application;
using namespace MPICommunications;

InviscidRegionSubdomain::InviscidRegionSubdomain() : nearbody(false) {}

InviscidRegionSubdomain::InviscidRegionSubdomain(bool nearbody_subdomain) : nearbody(nearbody_subdomain) {}

std::array<point_t, 2> InviscidRegionSubdomain::getExtentBox() const {
    std::array<point_t, 2> extent_box {{{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}}};

    for(const Vertex& vertex : vertices) {
        extent_box[0][0] = min(vertex.getCoordinate(0), extent_box[0][0]);
        extent_box[0][1] = min(vertex.getCoordinate(1), extent_box[0][1]);
        extent_box[1][0] = max(vertex.getCoordinate(0), extent_box[1][0]);
        extent_box[1][1] = max(vertex.getCoordinate(1), extent_box[1][1]);
    }

    return extent_box;
}

point_t InviscidRegionSubdomain::getCenterPoint() const {
    auto extent(getExtentBox());
    return midPoint(extent[0], extent[1]);
}

//Since the vertices are wound to form a simple polygon, we can offput creating the edge segments until we are ready
//to refine this subdomain
//This is mainly an optimization to reduce communication costs, but it also reduces the amount of memory used because
//each process will only be working on one subdomain at a time, and the edges are created and then destroyed once the
//subdomain has been refined. This is better than having all of the subdomains' edges be created and freed one-by-one
void InviscidRegionSubdomain::createBorder() {
    for(auto it=vertices.begin(); it<vertices.end()-1; ++it)
        edges.push_back({{(*it).getID(), (*(it+1)).getID()}});
    edges.push_back({{vertices.back().getID(), vertices.front().getID()}});
}

void InviscidRegionSubdomain::mesh(InviscidRegionMesh& owning_mesh) {
    if(not nearbody)
        createBorder();
    sort(vertices.begin(), vertices.end(), vertexCompareX);

    unordered_map<long, int> id_translator;
    vector<long> global_ids;
    global_ids.reserve(vertices.size());

    struct triangulateio in;

    in.numberofpoints = (int)vertices.size();
    in.numberofpointattributes = 0;
    in.pointlist = (double*) malloc(in.numberofpoints * 2 * sizeof(double));
    in.pointmarkerlist = (int*) malloc(in.numberofpoints * sizeof(int));

    for(int i=0; i<vertices.size(); ++i) {
        id_translator[vertices[i].getID()] = i;
        global_ids.push_back(vertices[i].getID());

        in.pointlist[i*2] = vertices[i].getCoordinate(0);
        in.pointlist[(i*2)+1] = vertices[i].getCoordinate(1);
        in.pointmarkerlist[i] = vertices[i].getBoundary();
    }

    in.numberofsegments = (int)edges.size();
    in.segmentlist = (int*) malloc(in.numberofsegments * 2 * sizeof(int));
    in.segmentmarkerlist = (int*) malloc(in.numberofsegments * sizeof(int));

    for(int i=0; i<edges.size(); ++i) {
        in.segmentlist[i*2] = id_translator[edges[i][0]];
        in.segmentlist[(i*2)+1] = id_translator[edges[i][1]];
        in.segmentmarkerlist[i] = 1;
    }

    if(nearbody) {
        const vector<double>& holes(owning_mesh.getBoundaryLayerHoles());
        in.numberofholes = (int)holes.size()/2;
        in.holelist = (double*) malloc(in.numberofholes * 2 * sizeof(double));
        for(int i=0; i<holes.size(); ++i)
            in.holelist[i] = holes[i];
    } else {
        in.numberofholes = 0;
        in.holelist = (double*) NULL;
    }

    in.numberofregions = 0;

    struct triangulateio out;
    out.pointlist = (double*) NULL;
    out.pointmarkerlist = (int*) NULL;
    out.trianglelist = (int*) NULL;
    out.segmentlist = (int*) NULL;
    out.segmentmarkerlist = (int*) NULL;

    triangulate((char*)"q33YYlpQjzu", &in, &out, (struct triangulateio*) NULL, owning_mesh.isotropic_area,
                owning_mesh.uniform, (int)owning_mesh.aabb_centers.size(), owning_mesh.triangle_aabb_centers.data());
    owning_mesh.globalizeTriangleOutput(out, global_ids);

    if(in.pointlist != NULL) free(in.pointlist);
    if(in.pointmarkerlist != NULL) free(in.pointmarkerlist);
    if(in.segmentlist != NULL) free(in.segmentlist);
    if(in.segmentmarkerlist != NULL) free(in.segmentmarkerlist);
    if(in.holelist != NULL) free(in.holelist);
    if(out.pointlist != NULL) free(out.pointlist);
    if(out.pointmarkerlist != NULL) free(out.pointmarkerlist);
    if(out.trianglelist != NULL) free(out.trianglelist);
    if(out.segmentlist != NULL) free(out.segmentlist);
    if(out.segmentmarkerlist != NULL) free(out.segmentmarkerlist);
}

void InviscidRegionSubdomain::sendSubdomain(int destination) {
    vector<int> sizes {(int)vertices.size(), (int)edges.size()};
    Send(sizes, 2, destination);
    Send(createMPIDataType(), *this, destination);
}

void InviscidRegionSubdomain::recvSubdomain(int source) {
    vector<int> sizes;
    Recv(sizes, 2, source);
    vertices.resize(sizes[0]);
    edges.resize(sizes[1]);
    Recv(createMPIDataType(), *this, source);
}

MPI_Datatype InviscidRegionSubdomain::createMPIDataType() {
    const int items = 4;
    int block_lengths[items] = {1, 1, 1, 1};

    MPI_Datatype vertex_t = Vertex::getMPIDataType();
    MPI_Datatype vec_t;
    MPI_Type_contiguous((int)vertices.size(), vertex_t, &vec_t);

    MPI_Datatype edge_t;
    MPI_Type_contiguous(2, MPI_LONG, &edge_t);
    MPI_Datatype edge_vec_t;
    MPI_Type_contiguous((int)edges.size(), edge_t, &edge_vec_t);

    MPI_Datatype types[items] = {vec_t, edge_vec_t, MPI_CXX_BOOL, MPI_INT};

    MPI_Aint offsets[items];
    MPI_Aint front_address;
    MPI_Get_address(&vertices, &front_address);
    MPI_Get_address(vertices.data(), &offsets[0]);
    offsets[0] -= front_address;
    MPI_Get_address(edges.data(), &offsets[1]);
    offsets[1] -= front_address;
    offsets[2] = offsetof(InviscidRegionSubdomain, nearbody);
    offsets[3] = offsetof(InviscidRegionSubdomain, cost);

    MPI_Datatype mpi_datatype;
    MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_datatype);
    MPI_Type_commit(&mpi_datatype);
    return mpi_datatype;
}
