//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "BoundaryLayerSubdomain.h"
#include "GeoPrimitives.h"
#include "Application.h"
#include "BoundaryLayerMesh.h"
#include "Vertex.h"
#include "MPICommunications.h"
#include "MeshingManager.h"
#include <unordered_map>
#include <map>
#include <cstddef>
#include <fstream>

extern "C" {
#include "triangle/triangle.h"
}

using namespace std;
using namespace Application;
using namespace MPICommunications;

int BoundaryLayerSubdomain::decomposition_threshold = 7;
long BoundaryLayerSubdomain::max_vertex_id = 0;

BoundaryLayerSubdomain::BoundaryLayerSubdomain()
        : x_vertices(), y_vertices(), vertices({&x_vertices, &y_vertices}), decomposition_level(0) {}

void BoundaryLayerSubdomain::initialize() {
    sub_subdomain = make_shared<BoundaryLayerSubdomain>();
    sub_subdomain->decomposition_level = ++decomposition_level;
    determineCut();
}

Decomposition BoundaryLayerSubdomain::decompose() {
    if(not containsInternalPoints())
        return Decomposition::ALREADY;

    initialize();
    computeMedian();
    projectPointsOnParaboloid();
    lowerConvexHull();
    auto decomposition = maintainSorted();

    cleanSubdomain();
    sub_subdomain->cleanSubdomain();

    if(sub_subdomain->x_vertices.empty() and sub_subdomain->y_vertices.empty()) {
        return Decomposition::ALREADY;
    } else if(x_vertices.empty() and y_vertices.empty()) {
        swap(x_vertices, sub_subdomain->x_vertices);
        swap(edges, sub_subdomain->edges);
        return Decomposition::ALREADY;
    } else {
        return decomposition;
    }
}

void BoundaryLayerSubdomain::determineCut() {
    if(((x_vertices.back().getCoordinate(0) - x_vertices.front().getCoordinate(0)) <
        (y_vertices.back().getCoordinate(1) - y_vertices.front().getCoordinate(1))))
        axis = Y_AXIS;
    else axis = X_AXIS;
}

bool BoundaryLayerSubdomain::containsInternalPoints() {
    for(const Vertex& vertex : x_vertices)
        if(not vertex.getBoundary())
            return true;
    return false;
}

void BoundaryLayerSubdomain::computeMedian() {
    median_vertex = &((*vertices[axis])[(static_cast<long>(floor(vertices[axis]->size()/2)))]);
}

void BoundaryLayerSubdomain::projectPointsOnParaboloid() {
    for(Vertex& vertex : *(vertices[!axis]))
        vertex.calculateProjected(median_vertex, axis);
}

void BoundaryLayerSubdomain::lowerConvexHull() {
    auto it = vertices[!axis]->begin();
    auto end = vertices[!axis]->end();

    vector<Vertex*> lch;

    while(it < end) {
        while((lch.size() >= 2) and lch[lch.size()-2]->orientation2D(lch[lch.size()-1], *it))
            lch.pop_back();

        lch.push_back(&(*it));
        ++it;
    }

    finalizeLowerConvexHull(lch);
}

void BoundaryLayerSubdomain::finalizeLowerConvexHull(std::vector<Vertex*>& lch) {
    if(lch.size() > 1)
        for(int i=0; i<lch.size()-1; ++i)
            if(not (lch[i]->getBoundary() and lch[i+1]->getBoundary()))
                edges.push_back(array<long, 2> {lch[i]->getID(), lch[i+1]->getID()});

    for(auto l_it=lch.begin(); l_it<lch.end(); ++l_it) {
        (*l_it)->setLowerConvexHull(true);
        (*l_it)->setBoundary(true);
    }

    lower_convex_hull.clear();
    lower_convex_hull.reserve(lch.size());
    transform(lch.begin(), lch.end(), back_inserter(lower_convex_hull), [](Vertex* v) {return *v;});
}

Decomposition BoundaryLayerSubdomain::maintainSorted() {
    copy(edges.begin(), edges.end(), back_inserter(sub_subdomain->edges));
    sub_subdomain->lower_convex_hull = lower_convex_hull;
    vector<Vertex> lch_primary_sorted(lower_convex_hull);
    sort(lch_primary_sorted.begin(), lch_primary_sorted.end(), vertex_compare[axis]);

    auto it = vertices[axis]->begin();
    for(const Vertex& vertex : lch_primary_sorted) {
        while(it->getID() != vertex.getID())
            ++it;
        it->setBoundary(true);
    }

    auto lch_median = lower_bound(lch_primary_sorted.begin(), lch_primary_sorted.end(), *median_vertex, 
                                  vertex_compare[axis]);
    auto vertex_median = vertices[axis]->begin() + (vertices[axis]->size()/2);

    long left_lch = distance(lch_median, lch_primary_sorted.end());
    double half_vertices = x_vertices.size()/2.0;
    long right_size = (long)(ceil(half_vertices) + lch_primary_sorted.size() - left_lch);

    bool left_decomposed = ((floor(half_vertices) + left_lch < 250) or (decomposition_level > decomposition_threshold));
    bool right_decomposed = ((right_size < 250) or (sub_subdomain->decomposition_level > decomposition_threshold));

    if(left_decomposed and right_decomposed) {
        if(not axis) {
            maintainPrimarySorted(lch_primary_sorted, lch_median);
        } else {
            sub_subdomain->x_vertices.reserve(right_size);
            maintainCutSorted();
        }

        y_vertices.clear();
        return Decomposition::BOTH;
    }
    else if(left_decomposed) {
        if(not axis) {
            maintainPrimarySorted(lch_primary_sorted, lch_median);
            sub_subdomain->y_vertices.reserve(right_size);
            createSubCutSorted();
        } else {
            move(lch_primary_sorted.begin(), lch_median, back_inserter(sub_subdomain->y_vertices));
            move(vertex_median, y_vertices.end(), back_inserter(sub_subdomain->y_vertices));
            sub_subdomain->x_vertices.reserve(right_size);
            maintainCutSorted();
        }

        y_vertices.clear();
        return Decomposition::LEFT;
    }
    else if(right_decomposed) {
        if(not axis) {
            maintainPrimarySorted(lch_primary_sorted, lch_median);
            maintainOriginalCutSorted();
        } else {
            move(lch_median, lch_primary_sorted.end(), vertex_median);
            sub_subdomain->x_vertices.reserve(right_size);
            maintainCutSorted();
        }

        sub_subdomain->y_vertices.clear();
        return Decomposition::RIGHT;
    }
    else {
        maintainPrimarySorted(lch_primary_sorted, lch_median);
        sub_subdomain->vertices[!axis]->reserve(right_size);
        maintainCutSorted();
        return Decomposition::NONE;
    }
}

void BoundaryLayerSubdomain::maintainPrimarySorted(std::vector<Vertex> &lch, std::vector<Vertex>::iterator lch_median) {
    auto vertex_median = vertices[axis]->begin() + (vertices[axis]->size()/2);
    move(lch.begin(), lch_median, back_inserter(*(sub_subdomain->vertices[axis])));
    move(vertex_median, vertices[axis]->end(), back_inserter(*(sub_subdomain->vertices[axis])));
    auto end = move(lch_median, lch.end(), vertex_median);
    vertices[axis]->resize(distance(vertices[axis]->begin(), end));
}

void BoundaryLayerSubdomain::maintainCutSorted() {
    double median_line = median_vertex->getCoordinate(axis);

    auto it = vertices[!axis]->begin();
    auto end = vertices[!axis]->end();
    auto s_it = vertices[!axis]->begin();
    vector<Vertex>& sub_vertices = *(sub_subdomain->vertices[!axis]);

    while(it < end) {
        while((it < end) and ((*it).useLowerConvexHull())) {
            sub_vertices.push_back(*it);
            (*s_it++) = move(*it++);
        }

        if(it < end) {
            if((*it).getCoordinate(axis) < median_line)
                (*s_it++) = move(*it++);
            else sub_vertices.push_back(*it++);
        }
    }

    vertices[!axis]->resize(distance(vertices[!axis]->begin(), s_it));
}

void BoundaryLayerSubdomain::maintainOriginalCutSorted() {
    double median_line = median_vertex->getCoordinate(axis);

    auto it = vertices[!axis]->begin();
    auto end = vertices[!axis]->end();
    auto s_it = vertices[!axis]->begin();

    while(it < end) {
        while((it < end) and ((*it).useLowerConvexHull()))
            (*s_it++) = move(*it++);

        if(it < end) {
            if((*it).getCoordinate(axis) < median_line)
                (*s_it++) = move(*it);

            ++it;
        }
    }

    vertices[!axis]->resize(distance(vertices[!axis]->begin(), s_it));
}

void BoundaryLayerSubdomain::createSubCutSorted() {
    double median_line = median_vertex->getCoordinate(axis);

    auto it = vertices[!axis]->begin();
    auto end = vertices[!axis]->end();
    vector<Vertex>& sub_vertices = *(sub_subdomain->vertices[!axis]);

    while(it < end) {
        while((it < end) and ((*it).useLowerConvexHull()))
            sub_vertices.push_back(*it++);

        if(it < end) {
            if(not ((*it).getCoordinate(axis) < median_line))
                sub_vertices.push_back(*it);

            ++it;
        }
    }
}

void BoundaryLayerSubdomain::cleanSubdomain() {
    auto existing_vertices(removeOldEdges());
    auto pinch_points(removePinchPoints());
    for_each(pinch_points.begin(), pinch_points.end(), [&existing_vertices](long id) {existing_vertices[id] = false;});
    edges.erase(remove_if(edges.begin(), edges.end(), [&existing_vertices](array<long, 2>& edge) {
        return not (existing_vertices[edge[0]] and existing_vertices[edge[1]]);
    }), edges.end());
}

std::vector<bool> BoundaryLayerSubdomain::removeOldEdges() {
    vector<bool> exist(max_vertex_id, false);
    for_each(x_vertices.begin(), x_vertices.end(), [&exist](const Vertex& v) {exist[v.getID()] = true;});
    edges.erase(remove_if(edges.begin(), edges.end(), [&exist](array<long, 2>& e) {
        return not (exist[e[0]] and exist[e[1]]);
    }), edges.end());
    return exist;
}

std::vector<long> BoundaryLayerSubdomain::removePinchPoints() {
    unordered_map<long, vector<long>> vertex_degrees;
    vector<long> pinch_points;

    for(const auto& edge : edges) {
        for(int i=0; i<2; ++i) {
            auto it = vertex_degrees.find(edge[i]);
            if(it == vertex_degrees.end())
                vertex_degrees.insert({edge[i], {edge[!i]}});
            else it->second.push_back(edge[!i]);
        }
    }

    auto it = vertex_degrees.begin();
    for(auto current=vertex_degrees.begin(); current!=vertex_degrees.end(); ++current) {
        it = current;
        while(current->second.size() < 2) {
            pinch_points.push_back(current->first);
            if(current->second.empty())
                break;
            auto neighbor = vertex_degrees.find(current->second.front());
            current->second.clear();
            neighbor->second.erase(remove(neighbor->second.begin(), neighbor->second.end(), current->first),
                                   neighbor->second.end());
            current = neighbor;
        }
        current = it;
    }

    for(const Vertex& vertex : lower_convex_hull)
        if(vertex_degrees.find(vertex.getID()) == vertex_degrees.end())
            pinch_points.push_back(vertex.getID());

    sort(pinch_points.begin(), pinch_points.end());

    x_vertices.erase(remove_if(x_vertices.begin(), x_vertices.end(), [&pinch_points](const Vertex& v) {
        return binary_search(pinch_points.begin(), pinch_points.end(), v.getID());
    }), x_vertices.end());
    y_vertices.erase(remove_if(y_vertices.begin(), y_vertices.end(), [&pinch_points](const Vertex& v) {
        return binary_search(pinch_points.begin(), pinch_points.end(), v.getID());
    }), y_vertices.end());

    return pinch_points;
}

void BoundaryLayerSubdomain::mesh(BoundaryLayerMesh& owning_mesh) {
    unordered_map<long, int> id_translator;
    vector<long> global_ids;
    global_ids.reserve(x_vertices.size());

    struct triangulateio in;

    in.numberofpoints = (int)x_vertices.size();
    in.numberofpointattributes = 0;
    in.pointlist = (double*) malloc(in.numberofpoints * 2 * sizeof(double));
    in.pointmarkerlist = (int*) malloc(in.numberofpoints * sizeof(int));

    for(int i=0; i<x_vertices.size(); ++i) {
        id_translator[x_vertices[i].getID()] = i;
        global_ids.push_back(x_vertices[i].getID());

        in.pointlist[i*2] = x_vertices[i].getCoordinate(0);
        in.pointlist[(i*2)+1] = x_vertices[i].getCoordinate(1);
        in.pointmarkerlist[i] = x_vertices[i].getBoundary();
    }

    in.numberofsegments = (int)edges.size();
    in.segmentlist = (int*) malloc(in.numberofsegments * 2 * sizeof(int));
    in.segmentmarkerlist = (int*) malloc(in.numberofsegments * sizeof(int));

    for(int i=0; i<edges.size(); ++i) {
        in.segmentlist[i*2] = id_translator[edges[i][0]];
        in.segmentlist[(i*2)+1] = id_translator[edges[i][1]];
        in.segmentmarkerlist[i] = 1;
    }

    const vector<double>& holes(owning_mesh.getModelHoles());
    in.numberofholes = (int)holes.size()/2;
    in.holelist = (double*) malloc(in.numberofholes * 2 * sizeof(double));
    for(int i=0; i<holes.size(); ++i)
        in.holelist[i] = holes[i];

    in.numberofregions = 0;

    struct triangulateio out;
    out.pointlist = (double*) NULL;
    out.pointmarkerlist = (int*) NULL;
    out.trianglelist = (int*) NULL;
    out.segmentlist = (int*) NULL;
    out.segmentmarkerlist = (int*) NULL;

    triangulate((char*)"lpQz", &in, &out, (struct triangulateio*)NULL, 0, false, 0, (double*)NULL);
    owning_mesh.addToMesherOutbox(out, global_ids);

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

void BoundaryLayerSubdomain::sendSubdomain(int destination) {
    array<int, 3> sizes {(int)x_vertices.size(), (int)y_vertices.size(), (int)edges.size()};
    MPI_Send(sizes.data(), 3, MPI_INT, destination, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
    Send(createMPIDataType(), *this, destination);
}

void BoundaryLayerSubdomain::recvSubdomain(int source) {
    array<int, 3> sizes;
    MPI_Recv(sizes.data(), 3, MPI_INT, source, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    x_vertices.resize(sizes[0]);
    y_vertices.resize(sizes[1]);
    edges.resize(sizes[2]);
    Recv(createMPIDataType(), *this, source);
}

MPI_Datatype BoundaryLayerSubdomain::createMPIDataType() {
    const int items = 4;
    int block_lengths[items] = {1, 1, 1, 1};

    MPI_Datatype vertex_t = Vertex::getMPIDataType();
    MPI_Datatype x_vec_t;
    MPI_Type_contiguous((int)x_vertices.size(), vertex_t, &x_vec_t);
    MPI_Datatype y_vec_t;
    MPI_Type_contiguous((int)y_vertices.size(), vertex_t, &y_vec_t);

    MPI_Datatype edge_t;
    MPI_Type_contiguous(2, MPI_LONG, &edge_t);
    MPI_Datatype edge_vec_t;
    MPI_Type_contiguous((int)edges.size(), edge_t, &edge_vec_t);

    MPI_Datatype types[items] = {x_vec_t, y_vec_t, edge_vec_t, MPI_INT};

    MPI_Aint offsets[items];
    MPI_Aint front_address;
    MPI_Get_address(&x_vertices, &front_address);
    MPI_Get_address(x_vertices.data(), &offsets[0]);
    offsets[0] -= front_address;
    MPI_Get_address(y_vertices.data(), &offsets[1]);
    offsets[1] -= front_address;
    MPI_Get_address(edges.data(), &offsets[2]);
    offsets[2] -= front_address;
    offsets[3] = offsetof(BoundaryLayerSubdomain, decomposition_level);

    MPI_Datatype mpi_datatype;
    MPI_Type_create_struct(items, block_lengths, offsets, types, &mpi_datatype);
    MPI_Type_commit(&mpi_datatype);
    return mpi_datatype;
}
