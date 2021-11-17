//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "BoundaryLayerMesh.h"
#include "Ray.h"
#include "BoundaryLayerSubdomain.h"
#include "MeshGenerator.h"
#include "MPICommunications.h"
#include "ADT2DExtent.h"
#include <fstream>
#include <float.h>
#include <cmath>
#include <assert.h>
#include <map>
#include <unistd.h>

extern "C" {
#include "triangle/triangle.h"
}

using namespace std;
using namespace Application;
using namespace MPICommunications;

BoundaryLayerMesh::BoundaryLayerMesh(MeshGenerator& owning_mesher)
        : mesher(owning_mesher), inviscid_region(&owning_mesher.inviscid_region) {}

void BoundaryLayerMesh::initializeModelSurface(std::string filename) {
    double x_coord;
    double y_coord;
    int vertex0;
    int vertex1;
    int id;
    int element = 0;
    int num_holes;
    double min_x = DBL_MAX;
    double max_x = -DBL_MAX;
    double min_y = DBL_MAX;
    double max_y = -DBL_MAX;

    ifstream input(filename);

    input >> num_model_vertices >> num_elements;
    vertices.reserve(num_model_vertices);

    for(int i=0; i<num_model_vertices; ++i) {
        input >> id >> x_coord >> y_coord >> element;
        assert(id == i);

        if(element == surface_element_start.size())
            surface_element_start.push_back(i);

        min_x = min(min_x, x_coord);
        max_x = max(max_x, x_coord);
        min_y = min(min_y, y_coord);
        max_y = max(max_y, y_coord);

        vertices.emplace_back(point_t{x_coord, y_coord}, true, id);
    }

    chord_length = max_x - min_x;
    max_extent_box.setDomain({{{{min_x, min_y}}, {{max_x, max_y}}}});
    mesh_center = midPoint(max_extent_box.getLowPoint(), max_extent_box.getHighPoint());

    transition_vertices.resize(num_elements);
    surface_element_start.push_back(num_model_vertices);
    next_vertex_id = num_model_vertices;

    input >> num_model_edges;
    mesher.final_edges = num_model_edges;
    edges.reserve(num_model_edges);

    for(int i=0; i< num_model_edges; ++i) {
        input >> id >> vertex0 >> vertex1;
        assert(id == i);

        if((vertex0 < num_model_vertices) and (vertex1 < num_model_vertices))
            edges.emplace_back(id, vertices[vertex0].getID(), vertices[vertex1].getID(), true);
        else throw logic_error("Vertex does not exist");
    }

    input >> num_holes;
    holes.reserve(2*num_holes);

    if(num_holes < 1)
        throw logic_error("No hole specified");

    for(int i=0; i<num_holes; ++i) {
        input >> id >> x_coord >> y_coord;
        holes.push_back(x_coord);
        holes.push_back(y_coord);
    }

    input.close();

    checkElementsWindingOrder();
}

void BoundaryLayerMesh::checkElementsWindingOrder() {
    for(int i=0; i<num_elements; ++i)
        if(isClockwise(vertices.begin()+surface_element_start[i], vertices.begin()+surface_element_start[i+1]-1))
            throw logic_error("Points should be ordered counter-clockwise");
}


const Vertex& BoundaryLayerMesh::getPreviousSurfaceVertex(int index) const {
    return vertices[index == 0 ? surface_element_start[1]-1 :
                    (getElement(vertices[index].getID()) == getElement(vertices[index-1].getID()) ?
                     index-1 : surface_element_start[getElement(vertices[index].getID())+1]-1)];
}

const Vertex& BoundaryLayerMesh::getNextSurfaceVertex(int index) const {
    return vertices[index == num_model_vertices-1 ? surface_element_start[getElement(vertices[index].getID())] :
                    (getElement(vertices[index].getID()) == getElement(vertices[index+1].getID()) ?
                     index+1 : surface_element_start[getElement(vertices[index].getID())])];
}

int BoundaryLayerMesh::getElement(long vertex_id) const {
    if(vertex_id < num_model_vertices) {
        for(int e=0; e<num_elements; ++e)
            if(vertex_id < surface_element_start[e+1])
                return e;
    } else {
        for(int e=0; e<num_elements; ++e)
            if(vertex_id < boundary_layer_element_start[e+1])
                return e;
    }
    return -1;
}

int BoundaryLayerMesh::getPreviousRayId(int index) const {
    return index == 0 ? ray_element_start[1]-1 :
           (rays[index].element == rays[index-1].element ? index-1 : ray_element_start[rays[index].element+1]-1);
}

int BoundaryLayerMesh::getNextRayId(int index) const {
    return index == rays.size()-1 ? ray_element_start[rays[index].element] :
           (rays[index].element == rays[index+1].element ? index+1 : ray_element_start[rays[index].element]);
}

Ray& BoundaryLayerMesh::getPreviousRay(int index) {
    return rays[getPreviousRayId(index)];
}

const Ray& BoundaryLayerMesh::getPreviousRay(int index) const {
    return rays[getPreviousRayId(index)];
}

const Ray& BoundaryLayerMesh::getNextRay(int index) const {
    return rays[getNextRayId(index)];
}

Ray& BoundaryLayerMesh::getNextRay(int index) {
    return rays[getNextRayId(index)];
}

void BoundaryLayerMesh::setGrowthFunction(double first_layer_thickness, double layer_growth_rate, int initial_layers) {
    last_thickness = first_layer_thickness;
    growth_rate = layer_growth_rate;
    num_layers = initial_layers;
    max_layers = static_cast<int>(ceil(1.25*num_layers));

    double offset = last_thickness;
    layer_offsets.resize(max_layers+1);
    for(int i=1; i<max_layers; ++i) {
        layer_offsets[i] = offset;
        last_thickness *= growth_rate;
        offset += last_thickness;
    }
    layer_offsets[max_layers] = offset;

    max_extent_box.inflateDomain(offset);
}

void BoundaryLayerMesh::insertBoundaryLayerPoints() {
    vector<Ray> local_rays(initializeLocalRays());
    Gatherv(local_rays, rays, 0);

    if(mesher.rank == 0)
        refineBoundaryLayer();

    Scatterv(rays, local_rays, 0);

    determineVertexInsertID(local_rays);
    createBoundaryLayerVertices(local_rays);
    gatherLocalVerticesOnRoot();
    closeBoundaryLayerExterior(local_rays);
}

std::vector<Ray> BoundaryLayerMesh::initializeLocalRays() {
    ray_element_start.insert(ray_element_start.end(), surface_element_start.begin(), surface_element_start.end());
    vector<Ray> local_rays;

    int start = num_model_vertices * mesher.rank / mesher.processes;
    int stop = (mesher.rank == mesher.processes-1) ? num_model_vertices :
               num_model_vertices * (mesher.rank+1)/mesher.processes;

    local_rays.reserve(stop-start);
    for(int i=start; i<stop; ++i)
        local_rays.emplace_back(vertices[i], num_layers, *this);
    return local_rays;
}

void BoundaryLayerMesh::determineVertexInsertID(const std::vector<Ray>& local_rays) {
    vector<int> sizes(mesher.processes);
    int local_size = 0;
    for(const Ray& ray : local_rays)
        local_size += ray.last_layer;

    MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    next_vertex_id = num_model_vertices;
    for(int p=0; p<mesher.rank; ++p)
        next_vertex_id += sizes[p];

    max_boundary_layer_vertex_id = next_vertex_id;
    for(int p=mesher.rank; p<mesher.processes; ++p)
        max_boundary_layer_vertex_id += sizes[p];

    BoundaryLayerSubdomain::max_vertex_id = max_boundary_layer_vertex_id;
}

void BoundaryLayerMesh::determineLocalBoundaryLayerSize(const std::vector<Ray>& local_rays) {
    unsigned long local_bl_size = 0;
    for(const Ray& ray : local_rays)
        local_bl_size += ray.last_layer;
    local_vertices.reserve(local_bl_size);
}

void BoundaryLayerMesh::createBoundaryLayerVertices(std::vector<Ray>& local_rays) {
    determineLocalBoundaryLayerSize(local_rays);
    long id;

    for(Ray& ray : local_rays) {
        id = next_vertex_id;
        next_vertex_id += ray.last_layer;

        for(int a=1; a<=ray.last_layer; ++a)
            local_vertices.emplace_back(ray.pointAtDistance(layer_offsets[a]), a==ray.last_layer, id++);

        if(ray.last_layer == 0)
            ray.last_vertex_id = ray.endpoint_id;
        else ray.last_vertex_id = id-1;
    }
}

void BoundaryLayerMesh::gatherLocalVerticesOnRoot() {
    vector<Vertex> boundary_layer;
    Gatherv(Vertex::getMPIDataType(), local_vertices, boundary_layer, 0);

    if(mesher.rank == 0)
        vertices.insert(vertices.end(), boundary_layer.begin(), boundary_layer.end());
}

void BoundaryLayerMesh::closeBoundaryLayerExterior(std::vector<Ray>& local_rays) {
    vector<long> last_vertices;
    vector<long> recv_vertices;
    last_vertices.reserve(local_rays.size());
    transform(local_rays.begin(), local_rays.end(), back_inserter(last_vertices),
              [](const Ray &r) {return r.last_vertex_id;});
    Gatherv(last_vertices, recv_vertices, 0);

    if(mesher.rank == 0) {
        for(int i=0; i<rays.size(); i++)
            rays[i].last_vertex_id = recv_vertices[i];
        createEnclosingEdges();
    }

    Broadcast(edges, 0);
}

void BoundaryLayerMesh::createEnclosingEdges() {
    Ray* current;
    Ray* neighbor;
    int difference;
    int element;
    next_edge_id = num_model_edges;
    num_enclosing_edges = 0;

    for(int i=0; i<rays.size(); ++i) {
        current = &rays[i];
        element = current->element;
        neighbor = &getNextRay(i);
        difference = current->last_layer - neighbor->last_layer;
        transition_vertices[element].push_back(&vertices[current->last_vertex_id]);

        if(fabs(difference) <= 1) {
            edges.emplace_back(next_edge_id++, current->last_vertex_id, neighbor->last_vertex_id, 3);
            ++num_enclosing_edges;
        } else {
            if(difference < 0) {
                swap(current, neighbor);
                difference *= -1;
            }

            if(neighbor->last_layer == 0) {
                edges.emplace_back(next_edge_id++, current->last_vertex_id - difference + 1, current->endpoint_id, 3);
                ++num_enclosing_edges;
            } else {
                edges.emplace_back(next_edge_id++, current->last_vertex_id - difference + 1,
                                   neighbor->last_vertex_id, 3);
                ++num_enclosing_edges;
            }

            for(int e=0; e<difference-1; e++) {
                edges.emplace_back(next_edge_id++, current->last_vertex_id-e, current->last_vertex_id-e-1, 3);
                ++num_enclosing_edges;
                vertices[current->last_vertex_id-e-1].setBoundary(true);
                transition_vertices[element].push_back(&vertices[current->last_vertex_id-e-1]);
            }
        }
    }

    for(int cusp : cusps)
        createEdgesAlongRay(rays[cusp]);
    for(int te : sharp_trailing_edges)
        createEdgesAlongRay(rays[te]);
}

double BoundaryLayerMesh::getAverageEnclosingTriangleArea() const {
    double area = 0;
    auto calculate_area = [this](const Ray& a, const Ray& b) -> double {
        int min_layer = min(a.last_layer, b.last_layer);
        long id = a.layerVertexID(min_layer);
        return triangleArea(vertices[id].getPoint(), vertices[id-1].getPoint(),
                            vertices[b.layerVertexID(min_layer)].getPoint());
    };

    for(int e=0; e<num_elements; ++e) {
        area += calculate_area(*(rays.begin() + ray_element_start[e+1] - 1), *(rays.begin() + ray_element_start[e]));

        const auto stop = rays.cbegin()+ray_element_start[e+1]-1;
        for(auto it=rays.cbegin()+ray_element_start[e]; it<stop; ++it)
            area += calculate_area(*it, *(it + 1));
    }

    return area/rays.size();
}

void BoundaryLayerMesh::createBoundaryLayerSubdomains() {
    total_initial_subdomains = 0;
    vector<array<point_t, 2>> element_extents;

    for(int e=0; e<num_elements; ++e) {
        vector<bl_subdomain_t> subdomains(createSubdomainsForElement(e));
        element_extents.push_back(finalizeSubdomains(subdomains));
    }

    inviscid_region->setNearBodyBoundingBox(element_extents);
}

void BoundaryLayerMesh::createEdgesAlongRay(const Ray& ray) {
    long last_id = ray.last_vertex_id;
    long first_id = last_id - ray.last_layer + 1;

    for(long i=first_id; i<last_id; ++i)
        vertices[i].setBoundary(true);

    edges.emplace_back(next_edge_id++, ray.endpoint_id, first_id, 0);
    for(long i=first_id; i<last_id; ++i)
        edges.emplace_back(next_edge_id++, i, i+1, 0);
}

std::vector<bl_subdomain_t> BoundaryLayerMesh::createSubdomainsForElement(int element) {
    vector<int> splits(getElementRaySplits(element));
    vector<bl_subdomain_t> subdomains;

    for(int i=0; i<splits.size(); ++i) {
        subdomains.push_back(make_shared<BoundaryLayerSubdomain>());
        subdomains.back()->decomposition_level = 1;
        vector<bool> exist(populateSubdomainVertices(subdomains.back(), splits[i], splits[(i+1) % splits.size()]));
        determineSubdomainEdges(subdomains.back(), exist);
    }

    return subdomains;
}

std::vector<int> BoundaryLayerMesh::getElementRaySplits(int element) {
    vector<int> element_splits;
    for(int cusp : cusps)
        if(rays[cusp].element == element)
            element_splits.push_back(cusp);
    for(int te : sharp_trailing_edges)
        if(rays[te].element == element)
            element_splits.push_back(te);
    sort(element_splits.begin(), element_splits.end());

    if(element_splits.empty()) {
        element_splits.push_back(ray_element_start[element]);
        element_splits.push_back((ray_element_start[element+1] + ray_element_start[element]) / 2);
        createEdgesAlongRay(rays[element_splits.front()]);
        createEdgesAlongRay(rays[element_splits.back()]);
    } else if(element_splits.size() == 1) {
        int half_num_rays = (ray_element_start[element+1] - ray_element_start[element]) / 2;
        if(element_splits.front() < (ray_element_start[element+1] + ray_element_start[element]) / 2)
            element_splits.push_back(element_splits.front() + half_num_rays);
        else element_splits.push_back(element_splits.front() - half_num_rays);
        createEdgesAlongRay(rays[element_splits.back()]);
    }

    return element_splits;
}

std::vector<bool> BoundaryLayerMesh::populateSubdomainVertices(bl_subdomain_t& subdomain, int ray_index, int last_ray) {
    vector<bool> exist(max_boundary_layer_vertex_id, false);
    Ray* ray = &rays[ray_index];

    while(true) {
        if(not exist[ray->endpoint_id]) {
            subdomain->x_vertices.push_back(vertices[ray->endpoint_id]);
            subdomain->y_vertices.push_back(vertices[ray->endpoint_id]);
            exist[ray->endpoint_id] = true;
        }

        const long last_id = ray->last_vertex_id;
        for(long id=last_id-ray->last_layer+1; id<=last_id; ++id) {
            subdomain->x_vertices.push_back(vertices[id]);
            subdomain->y_vertices.push_back(vertices[id]);
            exist[id] = true;
        }

        if(ray_index == last_ray)
            break;
        ray_index = getNextRayId(ray_index);
        ray = &rays[ray_index];
    }

    return exist;
}

void BoundaryLayerMesh::determineSubdomainEdges(bl_subdomain_t& subdomain, const std::vector<bool>& exist) {
    array<long, 2> vertex_ids;
    for(const Edge& edge : edges) {
        vertex_ids = edge.getVertices();

        if(exist[vertex_ids[0]] and exist[vertex_ids[1]])
            subdomain->edges.push_back(vertex_ids);
    }
}

std::array<point_t, 2> BoundaryLayerMesh::finalizeSubdomains(std::vector<bl_subdomain_t>& subdomains) {
    array<point_t, 2> extent {{{{DBL_MAX, DBL_MAX}}, {{-DBL_MAX, -DBL_MAX}}}};
    vector<double> recv_extent;
    vector<int> recv_processes;
    vector<bl_subdomain_t> my_subdomains;

    for(bl_subdomain_t& subdomain : subdomains) {
        ++total_initial_subdomains;
        if(next_recv_process != 0) {
            subdomain->sendSubdomain(next_recv_process);
            recv_processes.push_back(next_recv_process);
        } else {
            my_subdomains.push_back(move(subdomain));
        }

        next_recv_process = (next_recv_process+1) % mesher.processes;
    }

    for(bl_subdomain_t& subdomain : my_subdomains) {
        sort(subdomain->x_vertices.begin(), subdomain->x_vertices.end(), vertexCompareX);
        sort(subdomain->y_vertices.begin(), subdomain->y_vertices.end(), vertexCompareY);
        extent[0][0] = min(subdomain->x_vertices.front().getCoordinate(0), extent[0][0]);
        extent[1][0] = max(subdomain->x_vertices.back().getCoordinate(0), extent[1][0]);
        extent[0][1] = min(subdomain->y_vertices.front().getCoordinate(1), extent[0][1]);
        extent[1][1] = max(subdomain->y_vertices.back().getCoordinate(1), extent[1][1]);
    }

    initial_subdomains.insert(initial_subdomains.end(), my_subdomains.begin(), my_subdomains.end());

    MPI_Request request;
    MPI_Ibarrier(MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    for(int source : recv_processes) {
        Recv(recv_extent, 4, source);
        extent[0][0] = min(recv_extent[0], extent[0][0]);
        extent[1][0] = max(recv_extent[1], extent[1][0]);
        extent[0][1] = min(recv_extent[2], extent[0][1]);
        extent[1][1] = max(recv_extent[3], extent[1][1]);
    }

    return extent;
}

void BoundaryLayerMesh::decomposeInitialSubdomains() {
    Broadcast(total_initial_subdomains, 0);

    if(total_initial_subdomains >= mesher.processes) {
        mesher.boundary_layer_subdomains.insert(mesher.boundary_layer_subdomains.end(),
                                                make_move_iterator(initial_subdomains.begin()),
                                                make_move_iterator(initial_subdomains.end()));
        initial_subdomains.clear();
        return;
    }

    int iterations = static_cast<int>(ceil(log(mesher.processes)/log(total_initial_subdomains)));
    int increment;
    int partner;

    for(int i=1; i<iterations; ++i) {
        increment = static_cast<int>(pow(total_initial_subdomains, i));

        if(mesher.rank < increment) {
            partner = mesher.rank + increment;

            if(partner >= mesher.processes) {
                mesher.boundary_layer_subdomains.insert(mesher.boundary_layer_subdomains.end(),
                                                        make_move_iterator(initial_subdomains.begin()),
                                                        make_move_iterator(initial_subdomains.end()));
                initial_subdomains.clear();
            } else if(initial_subdomains.empty()) {
                int has_subdomain = 0;
                MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
            } else {
                auto subdomain(initial_subdomains.back());
                initial_subdomains.pop_back();

                switch(subdomain->decompose()) {
                    case Decomposition::BOTH: {
                        int has_subdomain = 1;
                        MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
                        subdomain->sub_subdomain->sendSubdomain(partner);
                        mesher.boundary_layer_subdomains.push_back(subdomain);
                        break;
                    } case Decomposition::LEFT: {
                        initial_subdomains.push_back(subdomain->sub_subdomain);
                        int has_subdomain = 1;
                        MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
                        subdomain->sendSubdomain(partner);
                        break;
                    } case Decomposition::RIGHT: {
                        int has_subdomain = 1;
                        MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
                        subdomain->sub_subdomain->sendSubdomain(partner);
                        initial_subdomains.push_back(subdomain);
                        break;
                    } case Decomposition::ALREADY: {
                        int has_subdomain = 0;
                        MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
                        mesher.boundary_layer_subdomains.push_back(subdomain);
                        break;
                    } case Decomposition::NONE: {
                        int has_subdomain = 2;
                        MPI_Send(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD);
                        subdomain->sub_subdomain->sendSubdomain(partner);
                        initial_subdomains.push_back(subdomain);
                        break;
                    }
                }
            }
        } else if(mesher.rank < 2*increment) {
            partner = mesher.rank - increment;
            int has_subdomain;
            MPI_Recv(&has_subdomain, 1, MPI_INT, partner, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(has_subdomain != 0) {
                bl_subdomain_t subdomain(make_shared<BoundaryLayerSubdomain>());
                subdomain->recvSubdomain(partner);
                if(has_subdomain == 1)
                    mesher.boundary_layer_subdomains.push_back(subdomain);
                else initial_subdomains.push_back(subdomain);
            }
        }
    }

    mesher.boundary_layer_subdomains.insert(mesher.boundary_layer_subdomains.end(),
                                            make_move_iterator(initial_subdomains.begin()),
                                            make_move_iterator(initial_subdomains.end()));
    initial_subdomains.clear();
}

const std::vector<double>& BoundaryLayerMesh::getModelHoles() const {
    return holes;
}

void BoundaryLayerMesh::addToMesherOutbox(struct triangulateio& out, std::vector<long>& global_ids) {
    mesher.addMeshToOutbox(out, global_ids);
}

void BoundaryLayerMesh::refineBoundaryLayer() {
    detectCuspsAndLargeAngles();
    handleSelfIntersections();
    gradationReductionControl();
    handleMultiElementIntersections();
    carveOutElementCavities();
    gradationGrowthControl();
    addFansAtSharpTrailingEdges();
}

void BoundaryLayerMesh::detectCuspsAndLargeAngles() {
    vector<double> thetas;
    thetas.reserve(num_model_vertices);
    double prev_theta;
    int next;

    for(auto e=ray_element_start.begin(); e<ray_element_start.end()-1; ++e)
        for(auto it=rays.begin()+*e; it<rays.begin()+*(e+1); ++it)
            thetas.push_back(angleBetweenVectors(it->normal_vector, getNextRay(it->endpoint_id).normal_vector));

    for(int i=0; i<num_model_vertices; ++i) {
        prev_theta = thetas[getPreviousRayId(i)];
        if(prev_theta > trailing_edge_angle_tolerance and thetas[i] > trailing_edge_angle_tolerance) {
            sharp_trailing_edges.push_back(i);
            rays[i].can_grade = false;
            getPreviousRay(i).can_grade = false;
            getNextRay(i).can_grade = false;
        } else if((prev_theta > cusp_angle_tolerance and thetas[i] > cusp_angle_tolerance) or
                  (prev_theta < -cusp_angle_tolerance and thetas[i] < -cusp_angle_tolerance)) {
            cusps.push_back(i);
        }
    }

    for(int i=0; i<num_model_vertices; ++i) {
        if(fabs(thetas[i]) > ray_angle_tolerance) {
            next = getNextRayId(i);
            if(not binary_search(cusps.begin(), cusps.end(), i) and not binary_search(cusps.begin(), cusps.end(), next)
                    and not binary_search(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), i)
               and not binary_search(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), next))
                large_angles.emplace_back(i, next, thetas[i]);
        }
    }
}

void BoundaryLayerMesh::gradationReductionControl() {
    const double tolerance = 0.55;
    point_t final;
    double distance;
    double layer_height;

    auto should_reduce = [&](Ray& ray, int index) -> bool {
        return calculateDistance(final, getPreviousRay(index).pointAtDistance(distance))/layer_height < tolerance and
               calculateDistance(final, getNextRay(index).pointAtDistance(distance))/layer_height < tolerance and
               ray.last_layer > 2;
    };

    for(int i=0; i<rays.size(); ++i) {
        Ray& ray = rays[i];
        distance = layer_offsets[ray.last_layer];
        layer_height = distance - layer_offsets[ray.last_layer - 1];
        final = ray.pointAtDistance(distance);

        if(should_reduce(ray, i)) {
            ray.can_grade = false;
            do {
                ray.last_layer--;
                distance = layer_offsets[ray.last_layer];
                layer_height = distance - layer_offsets[ray.last_layer - 1];
                final = ray.pointAtDistance(distance);
            } while(should_reduce(ray, i));
        }
    }

    removeRaySpikes();
}

void BoundaryLayerMesh::removeRaySpikes() {
    for(int i=0; i<rays.size(); ++i)
        rays[i].last_layer = min(rays[i].last_layer, max(getPreviousRay(i).last_layer, getNextRay(i).last_layer)+1);
}

void BoundaryLayerMesh::handleSelfIntersections() {
    for(int e=0; e<num_elements; ++e) {
        vector<Segment> test_segments(createRaySegments(computeBoundaryPoints(e), e));
        resolveSurfaceIntersections(e, test_segments);
        smoothPoorQualityRays(e, test_segments);
        resolveOuterBorderIntersections(e, test_segments);
    }
}

void BoundaryLayerMesh::resolveSurfaceIntersections(int element, std::vector<Segment>& test_segments) {
    setMaximumBoundaryLayerDomain(element);
    vector<Segment> surface(createSurfaceSegments(element));
    auto extent_intersections(findExtentIntersections(test_segments, surface, element, element));
    auto intersections(findIntersections(test_segments, surface, extent_intersections));
    clipSurfaceIntersections(intersections, test_segments);
}

std::vector<point_t> BoundaryLayerMesh::computeBoundaryPoints(int element) const {
    vector<point_t> points;
    points.reserve(ray_element_start[element+1] - ray_element_start[element]);
    transform(rays.begin()+ray_element_start[element], rays.begin()+ray_element_start[element+1], back_inserter(points),
              [this](const Ray& r) {return r.pointAtDistance(layer_offsets[r.last_layer]);});
    return points;
}

std::vector<point_t> BoundaryLayerMesh::computeExtendedBoundaryPoints(int element) const {
    vector<point_t> points;
    points.reserve(ray_element_start[element+1] - ray_element_start[element]);
    transform(rays.begin()+ray_element_start[element], rays.begin()+ray_element_start[element+1], back_inserter(points),
              [this](const Ray& r) {
                  return r.pointAtDistance(layer_offsets[r.last_layer] +
                                                   0.75*(layer_offsets[r.last_layer] - layer_offsets[r.last_layer-1]));
              });
    return points;
}

std::vector<Segment> BoundaryLayerMesh::createRaySegments(const std::vector<point_t>& end_points, int element) const {
    vector<Segment> segments;
    segments.reserve(surface_element_start[element+1] - surface_element_start[element]);
    transform(vertices.begin()+surface_element_start[element], vertices.begin()+surface_element_start[element+1],
              end_points.begin(), back_inserter(segments),
              [](const Vertex& v, const point_t& p) {return Segment(v.getPoint(), p);});
    return segments;
}

std::vector<Segment> BoundaryLayerMesh::createSurfaceSegments(int element) const {
    vector<Segment> segments;
    segments.reserve(surface_element_start[element+1] - surface_element_start[element]);
    transform(vertices.begin()+surface_element_start[element], vertices.begin()+surface_element_start[element+1]-1,
              vertices.begin()+surface_element_start[element]+1, back_inserter(segments),
              [](const Vertex& a, const Vertex& b) {return Segment(a.getPoint(), b.getPoint());});
    segments.push_back(Segment(vertices[surface_element_start[element+1]-1].getPoint(),
                               vertices[surface_element_start[element]].getPoint()));
    return segments;
}

std::vector<std::tuple<int, int, point_t>> BoundaryLayerMesh::findIntersections(
        const std::vector<Segment>& test_segments, const std::vector<Segment>& target_segments,
        const std::vector<std::array<int, 2>>& candidates) const {
    vector<tuple<int, int, point_t>> intersections;

    for(auto candidate : candidates) {
        auto points = test_segments[candidate[0]].intersectsAt(target_segments[candidate[1]]);

        for(auto point : points)
            intersections.push_back(make_tuple(candidate[0], candidate[1], point));
    }

    return intersections;
}

void BoundaryLayerMesh::clipSurfaceIntersections(const std::vector<std::tuple<int, int, point_t>>& intersections,
                                                 std::vector<Segment>& test_segments) {
    Ray *ray;
    for(const auto& it : intersections) {
        ray = &rays[get<0>(it)];
        int layers = (unsigned)distance(layer_offsets.begin(),
                                        lower_bound(layer_offsets.begin(), layer_offsets.end(),
                                                    calculateDistance(ray->point, get<2>(it)))) - 1;

        ray->decreaseLayers(layers);
        ray->can_grade = false;
        test_segments[get<0>(it)].setB(ray->pointAtDistance(layer_offsets[ray->last_layer]));
    }
}

void BoundaryLayerMesh::resolveOuterBorderIntersections(int element, std::vector<Segment>& test_segments) {
    vector<Segment> border(createOuterBorderSegments(test_segments));
    auto extent_intersections(findExtentIntersections(test_segments, border, element, element));
    auto intersections(findIntersections(test_segments, border, extent_intersections));
    prioritizeSelfIntersections(intersections, element);
    clipOuterBorderIntersections(intersections, test_segments, border, element);
}

std::vector<Segment> BoundaryLayerMesh::createOuterBorderSegments(const std::vector<Segment>& ray_segments) {
    vector<Segment> segments;
    segments.reserve(ray_segments.size());
    transform(ray_segments.begin(), ray_segments.end()-1, ray_segments.begin()+1, back_inserter(segments),
              [](const Segment& a, const Segment& b) {return Segment(a.getPointA(), b.getPointB());});
    segments.push_back(Segment(ray_segments.back().getPointA(), ray_segments.front().getPointB()));
    return segments;
}

std::vector<Segment> BoundaryLayerMesh::createOuterBorderSegments(const std::vector<point_t>& boundary_points) {
    vector<Segment> segments;
    segments.reserve(boundary_points.size());
    transform(boundary_points.begin(), boundary_points.end()-1, boundary_points.begin()+1, back_inserter(segments),
              [](const point_t& a, const point_t& b) {return Segment(a, b);});
    segments.push_back(Segment(boundary_points.back(), boundary_points.front()));
    return segments;
}

void BoundaryLayerMesh::setMaximumBoundaryLayerDomain(int element) {
    auto bounding_box = computeBoundingBox(vertices.begin()+surface_element_start[element],
                                           vertices.begin()+surface_element_start[element+1]);
    const double bl_thickness = layer_offsets.back();
    bounding_box[0][0] -= bl_thickness;
    bounding_box[0][1] -= bl_thickness;
    bounding_box[1][0] += bl_thickness;
    bounding_box[1][1] += bl_thickness;
    max_element_extent_boxes.emplace_back(bounding_box);
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::findExtentIntersections(
        const std::vector<Segment>& test_segments, const std::vector<Segment>& target_segments, int test_element,
        int target_element) const {
    const AABB* domain;
    if(test_element == target_element)
        domain = &max_element_extent_boxes[test_element];
    else domain = &max_extent_box;

    ADT2DExtent adt(domain->getExtent());
    vector<array<int, 2>> extent_intersections;

    for(int i=0; i<target_segments.size(); ++i)
        adt.store(i, target_segments[i].getExtent());

    for(int i=0; i<test_segments.size(); ++i) {
        auto extents = adt.retrieve(test_segments[i].getExtent());
        for(auto extent : extents)
            if(test_element != target_element or (i != extent and i != (extent + 1)%test_segments.size()))
                extent_intersections.push_back({{i, extent}});
    }

    return extent_intersections;
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::findExtentIntersections(const std::vector<Segment>& test_segments,
                                                                           int element) const {
    ADT2DExtent adt(max_element_extent_boxes[element].getExtent());
    vector<array<int, 2>> extent_intersections;

    for(int i=0; i<test_segments.size(); ++i)
        adt.store(i, test_segments[i].getExtent());

    for(int i=0; i<test_segments.size(); ++i) {
        auto extents = adt.retrieve(test_segments[i].getExtent());
        for(auto extent : extents)
            if(i < extent)
                extent_intersections.push_back({{i, extent}});
    }

    return extent_intersections;
}

void BoundaryLayerMesh::clipOuterBorderIntersections(
        const std::vector<std::tuple<int, int, point_t>>& intersections, std::vector<Segment>& test_segments,
        std::vector<Segment>& border, int element) {
    int ray_id;
    int boundary_id;
    Segment *ray_segment;
    Segment *boundary;
    Ray *ray;
    array<Ray*, 2> boundary_rays;
    const int offset_id = ray_element_start[element];

    for(const auto& it : intersections) {
        ray_id = get<0>(it);
        ray = &rays[offset_id + ray_id];
        ray_segment = &test_segments[ray_id];
        ray->can_grade = false;

        boundary_id = get<1>(it);
        boundary = &border[boundary_id];
        boundary_rays[0] = &rays[offset_id + boundary_id];
        boundary_rays[1] = &getNextRay(offset_id + boundary_id);
        boundary_rays[0]->can_grade = false;
        boundary_rays[1]->can_grade = false;

        while(ray_segment->doesIntersect(*boundary)) {
            if(ray->last_layer > max(boundary_rays[0]->last_layer, boundary_rays[1]->last_layer)) {
                --(ray->last_layer);
                ray_segment->setB(ray->pointAtDistance(layer_offsets[ray->last_layer]));
                border[ray_id].setA(ray_segment->getPointB());
                border[(ray_id+border.size()-1)%border.size()].setB(ray_segment->getPointB());
            } else {
                if(boundary_rays[0]->last_layer >= boundary_rays[1]->last_layer) {
                    --(boundary_rays[0]->last_layer);
                    boundary->setA(boundary_rays[0]->pointAtDistance(layer_offsets[boundary_rays[0]->last_layer]));
                    test_segments[boundary_id].setB(boundary->getPointA());
                    border[(boundary_id+border.size()-1)%border.size()].setB(boundary->getPointA());
                } else {
                    --(boundary_rays[1]->last_layer);
                    boundary->setB(boundary_rays[1]->pointAtDistance(layer_offsets[boundary_rays[1]->last_layer]));
                    test_segments[(boundary_id+1)%border.size()].setB(boundary->getPointB());
                    border[(boundary_id+1)%border.size()].setA(boundary->getPointB());
                }
            }
        }
    }
}

void BoundaryLayerMesh::prioritizeSelfIntersections(std::vector<std::tuple<int, int, point_t>>& intersections,
                                                    int element) {
    const int offset_id = ray_element_start[element];
    sort(intersections.begin(), intersections.end(), [&](const tuple<int, int, point_t>& a,
                                                         const tuple<int, int, point_t>& b) {
        point_t intersection_a = get<2>(a);
        point_t intersection_b = get<2>(b);
        double distance_a = calculateDistance(rays[offset_id + get<0>(a)].point, intersection_a) +
                calculateDistance(rays[offset_id + get<1>(a)].point, intersection_a);
        double distance_b = calculateDistance(rays[offset_id + get<0>(b)].point, intersection_b) +
                calculateDistance(rays[offset_id + get<1>(b)].point, intersection_b);
        return distance_a < distance_b;
    });
}

void BoundaryLayerMesh::smoothPoorQualityRays(int element, std::vector<Segment>& test_segments) {
    auto extent_intersections(findExtentIntersections(test_segments, element));
    auto intersections(findIntersections(test_segments, test_segments, extent_intersections));
    auto ranges(determineSmoothingRegions(intersections, element));

    for(const auto& range : ranges)
        raySmoothing(range[0], range[1]);
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::determineSmoothingRegions(
        const std::vector<std::tuple<int, int, point_t>>& intersections, int element) {
    vector<bool> crosses(ray_element_start[element+1] - ray_element_start[element], false);
    for(const auto& it : intersections) {
        crosses[get<0>(it)] = true;
        crosses[get<1>(it)] = true;
    }

    setNeighborhoodsForContiguousIntersections(crosses, element);
    setNeighborhoodsForCuspsAndLargeAngles(crosses, element);

    return determineRayRangesToSmooth(crosses, element);
}

void BoundaryLayerMesh::setNeighborhoodsForContiguousIntersections(std::vector<bool>& crosses, int element) {
    const int offset_id = ray_element_start[element];
    int start;
    int stop;
    int size;
    int neighborhood;
    int id;
    vector<array<int, 3>> ranges;

    if(crosses[0]) {
        start = static_cast<int>(crosses.size() - 1);
        while(crosses[start])
            --start;
        start = getNextRayId(start + offset_id) - offset_id;
    } else {
        start = -1;
    }

    for(int i=1; i<crosses.size(); ++i) {
        if(crosses[i]) {
            if(start == -1) {
                start = i;
            }
        } else if(start != -1) {
            stop = i - 1 + offset_id;
            start += offset_id;

            if(start <= stop)
                size = stop - start;
            else size = ray_element_start[element+1] - offset_id - (start - stop);

            if(size > 3) {
                neighborhood = static_cast<int>(ceil(2.0 * size));
                ranges.push_back({{start, stop, neighborhood}});
            }

            start = -1;
        }
    }

    for(const auto& range : ranges) {
        neighborhood = range[2];
        id = range[1];
        for(int r=0; r<neighborhood; ++r) {
            id = getNextRayId(id);
            crosses[id - offset_id] = true;
        }
        id = range[0];
        for(int r=0; r<neighborhood; ++r) {
            id = getPreviousRayId(id);
            crosses[id - offset_id] = true;
        }
    }
}

void BoundaryLayerMesh::setNeighborhoodsForCuspsAndLargeAngles(std::vector<bool>& crosses, int element) {
    const int offset_id = ray_element_start[element];
    int id;
    int neighborhood;

    for(const auto& it : large_angles) {
        if(rays[get<0>(it)].element != element)
            continue;

        if(fabs(get<2>(it)) < 15.0)
            neighborhood = static_cast<int>(ceil(fabs(get<2>(it))));
        else neighborhood = static_cast<int>(ceil(fabs(get<2>(it)) / 2.0));

        id = get<0>(it);
        crosses[id-offset_id] = true;
        for(int i=0; i<neighborhood; ++i) {
            id = getPreviousRayId(id);
            crosses[id-offset_id] = true;
        }
        id = get<1>(it);
        crosses[id-offset_id] = true;
        for(int i=0; i<neighborhood; ++i) {
            id = getNextRayId(id);
            crosses[id-offset_id] = true;
        }
    }

    neighborhood = 20;
    for(int cusp : cusps) {
        if(rays[cusp].element != element)
            continue;

        crosses[cusp-offset_id] = true;
        id = cusp;
        for(int i=0; i<neighborhood; ++i) {
            id = getNextRayId(id);
            crosses[id-offset_id] = true;
        }
        id = cusp;
        for(int i=0; i<neighborhood; ++i) {
            id = getPreviousRayId(id);
            crosses[id-offset_id] = true;
        }
    }
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::determineRayRangesToSmooth(const std::vector<bool>& crosses,
                                                                              int element) {
    const int offset_id = ray_element_start[element];
    vector<array<int, 2>> ranges;
    int start;
    int stop;
    int size;

    if(crosses[0]) {
        start = static_cast<int>(crosses.size() - 1);
        while(crosses[start])
            --start;
        start = getNextRayId(start + offset_id) - offset_id;
    } else {
        start = -1;
    }

    for(int i=1; i<crosses.size(); ++i) {
        if(crosses[i]) {
            if(start == -1)
                start = i;
        } else if(start != -1) {
            stop = i - 1 + offset_id;
            start += offset_id;

            if(start <= stop)
                size = stop - start;
            else size = ray_element_start[element+1] - offset_id - (start - stop);

            if(size > 3) {
                auto fixed = getFixedRays(start, stop);

                for(int r : fixed) {
                    if(start != r and start != getPreviousRayId(r))
                        ranges.push_back({{start, getPreviousRayId(r)}});
                    start = getNextRayId(r);
                }

                ranges.push_back({{start, stop}});
            }

            start = -1;
        }
    }

    return ranges;
}

std::vector<int> BoundaryLayerMesh::getFixedRays(int start, int stop) {
    vector<int> fixed;

    if(start < stop) {
        fixed.insert(fixed.end(), lower_bound(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), start),
                     upper_bound(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), stop));
    } else {
        fixed.insert(fixed.end(), sharp_trailing_edges.begin(),
                     upper_bound(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), stop));
        fixed.insert(fixed.end(), lower_bound(sharp_trailing_edges.begin(), sharp_trailing_edges.end(), start),
                     sharp_trailing_edges.end());
    }

    vector<int> te_neighbors;
    for(int te : fixed) {
        te_neighbors.push_back(getPreviousRayId(te));
        te_neighbors.push_back(getNextRayId(te));
    }
    fixed.insert(fixed.end(), te_neighbors.begin(), te_neighbors.end());

    if(start < stop) {
        fixed.insert(fixed.end(), lower_bound(cusps.begin(), cusps.end(), start),
                     upper_bound(cusps.begin(), cusps.end(), stop));
    } else {
        fixed.insert(fixed.end(), cusps.begin(), upper_bound(cusps.begin(), cusps.end(), stop));
        fixed.insert(fixed.end(), lower_bound(cusps.begin(), cusps.end(), start), cusps.end());
    }

    sort(fixed.begin(), fixed.end());

    if(start > stop) {
        const int element = rays[start].element;
        fixed.erase(remove_if(fixed.begin(), fixed.end(), [&](int id) {return rays[id].element != element;}),
                    fixed.end());
        rotate(fixed.begin(), lower_bound(fixed.begin(), fixed.end(), start), fixed.end());
    }

    return fixed;
}

void BoundaryLayerMesh::raySmoothing(int start, int stop) {
    int id;
    Ray* current;
    Ray* prev_ray;
    Ray* next_ray;
    vector<vector_t> normal_prime;
    int distance;
    if(stop < start) {
        distance = 1;
        int temp = start;
        while(temp != stop) {
            distance++;
            temp = getNextRayId(temp);
        }
    } else {
        distance = stop - start + 1;
    }
    normal_prime.resize(distance);

    double mag;
    double residual;

    int layers = max_layers;
    id = start;
    for(int i=0; i<distance; ++i) {
        layers = min(layers, rays[id].last_layer);
        id = getNextRayId(id);
    }

    double height = layer_offsets[layers];

    do {
        id = start;
        for(int i=0; i<distance; ++i) {
            current = &rays[id];
            prev_ray = &getPreviousRay(id);
            next_ray = &getNextRay(id);

            point_t mid = midPoint(prev_ray->pointAtDistance(height), next_ray->pointAtDistance(height));
            point_t temp = current->point;
            temp[0] -= mid[0];
            temp[1] -= mid[1];
            temp[0] /= height;
            temp[1] /= height;

            mag = current->getMagnitude()/calculateMagnitude(temp);
            temp[0] *= mag;
            temp[1] *= mag;

            normal_prime[i] = temp;
            id = getNextRayId(id);
        }

        residual = 0.0;
        id = start;
        for(int i=0; i<distance; i++) {
            residual += vectorDifference(normal_prime[i], rays[id].normal_vector);
            id = getNextRayId(id);
        }

        id = start;
        for(int i=0; i<distance; i++) {
            current = &rays[id];
            current->normal_vector = normal_prime[i];
            id = getNextRayId(id);
        }
    } while(residual > 1.0e-12);
}

void BoundaryLayerMesh::handleMultiElementIntersections() {
    if(num_elements == 1)
        return;

    vector<vector<point_t>> element_boundary_points;
    for(int i=0; i<num_elements; ++i)
        element_boundary_points.push_back(computeBoundaryPoints(i));

    vector<AABB> element_aabbs;
    for(int i=0; i<num_elements; ++i)
        element_aabbs.emplace_back(computeBoundingBox(element_boundary_points[i].begin(),
                                                      element_boundary_points[i].end()));
    for(int test=num_elements-1; test>=0; --test) {
        auto test_segments(createRaySegments(computeExtendedBoundaryPoints(test), test));

        for(int target=test-1; target>=0; --target) {
            if(max_element_extent_boxes[test].intersects(max_element_extent_boxes[target])) {
                auto aabb_intersections(pruneRayAABBIntersections(test_segments, element_aabbs[target]));
                auto target_segments(createOuterBorderSegments(element_boundary_points[target]));
                auto extent_intersections(findExtentIntersections(test_segments, aabb_intersections, target_segments));
                auto intersections(findIntersections(test_segments, target_segments, extent_intersections));
                clipMultiElementIntersections(intersections, test_segments, test, target_segments, target);

                for(int i=0; i<test_segments.size(); ++i)
                    element_boundary_points[test][i] = test_segments[i].getPointB();
                for(int i=0; i<target_segments.size(); ++i)
                    element_boundary_points[target][i] = target_segments[i].getPointA();

                //auto surface_segments = createSurfaceSegments(target);
                //surface_element_start[element]
            }
        }
    }
}

void BoundaryLayerMesh::gradationGrowthControl() {
    const double add_tolerance = 1.5;
    //const double angle_tolerance = 2.5;

    int offset_id;
    point_t final;
    Segment test;
    double distance;
    double layer_height;

    bool continue_grading;

    vector<vector<point_t>> boundary_points;
    vector<vector<Segment>> boundary;
    ADT2DExtent boundary_adt(max_extent_box.getExtent());

    for(int e=0; e<num_elements; ++e) {
        boundary_points.push_back(computeBoundaryPoints(e));
        boundary.push_back(createOuterBorderSegments(boundary_points[e]));
        offset_id = ray_element_start[e];
        for(int i=0; i<boundary.size(); ++i)
            boundary_adt.store(i + offset_id, boundary[e][i].getExtent());
    }

    for(int e=0; e<num_elements; ++e) {
        ADT2DExtent ray_adt(max_extent_box.getExtent());

        auto ray_segments(createRaySegments(boundary_points[e], e));
        offset_id = ray_element_start[e];
        for(int i=0; i<ray_segments.size(); ++i)
            ray_adt.store(i + offset_id, ray_segments[i].getExtent());

        do {
            continue_grading = false;
            for(int i=ray_element_start[e]; i<ray_element_start[e+1]; ++i) {
                Ray& ray = rays[i];
                if((not ray.can_grade) or ray.last_layer == max_layers)
                    continue;

                Ray& prev_ray = getPreviousRay(i);
                Ray& next_ray = getNextRay(i);

                if(ray.last_layer > prev_ray.last_layer or ray.last_layer > next_ray.last_layer)
                    continue;

                distance = layer_offsets[ray.last_layer];
                final = ray.pointAtDistance(distance);
                layer_height = calculateDistance(final, ray.pointAtDistance(layer_offsets[ray.last_layer-1]));

                if(add_tolerance > calculateDistance(final, prev_ray.pointAtDistance(distance))/layer_height and
                        add_tolerance > calculateDistance(final, next_ray.pointAtDistance(distance))/layer_height) {
                    ray.can_grade = false;
                    continue;
                }

                /*if(angle_tolerance < angleBetweenVectors(prev_ray.normal_vector, ray.normal_vector) and
                        angle_tolerance < angleBetweenVectors(ray.normal_vector, next_ray.normal_vector)) {
                    if(ray.last_layer > min(prev_ray.last_layer, next_ray.last_layer) - 1)
                        continue;
                } else {
                    if(ray.last_layer > min(prev_ray.last_layer, next_ray.last_layer) - 0)
                        continue;
                }*/

                test.setA(ray.point);
                test.setB(ray.pointAtDistance(2*layer_offsets[ray.last_layer+1] - layer_offsets[ray.last_layer]));

                auto extents = ray_adt.retrieve(test.getExtent());
                for(auto extent : extents)
                    if(i != extent)
                        if(test.doesIntersect(ray_segments[extent])) {
                            ray.can_grade = false;
                            continue;
                        }

                extents = boundary_adt.retrieve(test.getExtent());
                for(auto extent : extents)
                    if(i != extent and i != getNextRayId(extent))
                        if(test.doesIntersect(boundary[e][extent])) {
                            ray.can_grade = false;
                            continue;
                        }

                continue_grading = true;

                ++ray.last_layer;
                ray_adt.removeFirst(i, ray_segments[i-offset_id].getExtent());
                boundary_adt.removeFirst(i, boundary[e][i-offset_id].getExtent());
                boundary_adt.removeFirst(getPreviousRayId(i), boundary[e][getPreviousRayId(i)-offset_id].getExtent());

                final = ray.pointAtDistance(layer_offsets[ray.last_layer]);
                ray_segments[i-offset_id].setB(final);
                ray_adt.store(i, ray_segments[i-offset_id].getExtent());
                boundary[e][i-offset_id].setA(final);
                boundary[e][getPreviousRayId(i)-offset_id].setB(final);
                boundary_adt.store(i, boundary[e][i - offset_id].getExtent());
                boundary_adt.store(getPreviousRayId(i), boundary[e][getPreviousRayId(i) - offset_id].getExtent());
            }
        } while(continue_grading);
    }
}

void BoundaryLayerMesh::addFansAtSharpTrailingEdges() {
    int ray_id;
    Ray* te_ray;
    array<Ray*, 2> next_rays;
    array<Ray*, 2> prev_rays;
    double desired_width;
    double te_width;
    int num_rays;
    vector<array<int, 3>> increments;

    for(auto te=sharp_trailing_edges.rbegin(); te<sharp_trailing_edges.rend(); ++te) {
        te_ray = &rays[*te];
        ray_id = getNextRayId(*te);
        next_rays[0] = &rays[ray_id];
        next_rays[1] = &getNextRay(ray_id);
        ray_id = getPreviousRayId(*te);
        prev_rays[0] = &rays[ray_id];
        prev_rays[1] = &getPreviousRay(ray_id);

        te_ray->last_layer = min(next_rays[1]->last_layer, prev_rays[1]->last_layer);
        next_rays[0]->last_layer = te_ray->last_layer;
        next_rays[1]->last_layer = te_ray->last_layer;
        prev_rays[0]->last_layer = te_ray->last_layer;
        prev_rays[1]->last_layer = te_ray->last_layer;

        desired_width = layer_offsets[te_ray->last_layer-1] - layer_offsets[te_ray->last_layer-2];
        te_width = min(calculateDistance(next_rays[0]->pointAtDistance(layer_offsets[te_ray->last_layer]),
                                         te_ray->pointAtDistance(layer_offsets[te_ray->last_layer])),
                       calculateDistance(prev_rays[0]->pointAtDistance(layer_offsets[te_ray->last_layer]),
                                         te_ray->pointAtDistance(layer_offsets[te_ray->last_layer])));

        num_rays = static_cast<int>(floor(te_width/desired_width)) - 1;
        auto next_added_rays(createHalfFan(next_rays[0], te_ray, num_rays, vertices[te_ray->endpoint_id], true));
        auto prev_added_rays(createHalfFan(prev_rays[0], te_ray, num_rays, vertices[te_ray->endpoint_id], false));

        increments.push_back({{te_ray->element, num_rays, *te}});
        rays.insert(rays.begin() + *te + 1, next_added_rays.begin(), next_added_rays.end());
        rays.insert(rays.begin() + *te, prev_added_rays.begin(), prev_added_rays.end());
    }

    incrementElementRayStarts(increments);
}

std::vector<Ray> BoundaryLayerMesh::createHalfFan(const Ray* a, const Ray* b, int size, const Vertex& origin,
                                                  bool reverse) {
    vector<Ray> fan;
    fan.reserve(size);
    const int layers = a->last_layer;
    const int element = a->element;
    array<double, 2> normal = a->normal_vector;
    array<double, 2> normal_delta {(b->normal_vector[0]-a->normal_vector[0])/(size),
                                  (b->normal_vector[1]-a->normal_vector[1])/(size)};

    fan.emplace_back(origin, normal, layers, element);
    for(int i=0; i<size-1; ++i) {
        normal[0] += normal_delta[0];
        normal[1] += normal_delta[1];
        fan.emplace_back(origin, normal, layers, element);
    }

    if(reverse)
        std::reverse(fan.begin(), fan.end());

    return fan;
}

void BoundaryLayerMesh::incrementElementRayStarts(std::vector<std::array<int, 3>>& ray_increments) {
    int element;
    int half_added;
    int trailing_edge;

    for(auto& increment : ray_increments) {
        element = increment[0];
        half_added = increment[1];
        trailing_edge = increment[2];

        for(auto cusp=cusps.rbegin(); cusp<cusps.rend(); ++cusp)
            if(*cusp > trailing_edge)
                (*cusp) += (2*half_added);
            else if(*cusp == trailing_edge)
                (*cusp) += half_added;
            else break;
        for(auto te=sharp_trailing_edges.rbegin(); te<sharp_trailing_edges.rend(); ++te)
            if(*te > trailing_edge)
                (*te) += (2*half_added);
            else if(*te == trailing_edge)
                (*te) += half_added;
            else break;
        for(auto angle_tuple=large_angles.rbegin(); angle_tuple<large_angles.rend(); ++angle_tuple) {
            if(get<0>(*angle_tuple) > trailing_edge)
                get<0>(*angle_tuple) += (2*half_added);
            else if(get<0>(*angle_tuple) == trailing_edge)
                get<0>(*angle_tuple) += half_added;
            if(get<1>(*angle_tuple) > trailing_edge)
                get<1>(*angle_tuple) += (2*half_added);
            else if(get<1>(*angle_tuple) == trailing_edge)
                get<1>(*angle_tuple) += half_added;
        }

        for(int e=element+1; e<=num_elements; ++e)
            ray_element_start[e] += (2*half_added);
    }
}

void BoundaryLayerMesh::carveOutElementCavities() {
    for(int e=0; e<num_elements; ++e)
        smoothCavities(detectElementCavities(e));
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::detectElementCavities(int element) {
    const int offset_id = ray_element_start[element];
    const int stop_id = ray_element_start[element+1];
    const int neighbor_distance = 5;
    vector<array<int, 2>> cavities;
    int start;
    int size;
    int id;
    bool found;

    if(rays[offset_id].last_layer < num_layers) {
        start = offset_id;

        do {
            found = false;
            id = getPreviousRayId(start);
            for(int i=0; i<neighbor_distance; ++i) {
                if(rays[id].last_layer < num_layers) {
                    start = id;
                    found = true;
                    break;
                } else {
                    id = getPreviousRayId(id);
                }
            }
        } while(not found);
    } else {
        start = -1;
    }

    for(int i=offset_id+1; i<stop_id; ++i) {
        if(rays[i].last_layer < num_layers) {
            if(start == -1)
                start = i;
        } else if(start != -1) {
            found = false;
            if(i != stop_id-1) {
                id = getNextRayId(i);
                for(int j=0; j<neighbor_distance; ++j) {
                    if(id == stop_id-1) {
                        i = id;
                        break;
                    } else if(rays[id].last_layer < num_layers) {
                        i = id;
                        found = true;
                        break;
                    } else {
                        id = getNextRayId(id);
                    }
                }
            }

            if(not found) {
                if(start < i)
                    size = i - start;
                else size = stop_id - offset_id - start + i + 1;

                if(size > 10)
                    cavities.push_back({{start, i-1}});

                start = -1;
            }
        }
    }

    return cavities;
}

void BoundaryLayerMesh::smoothCavities(const std::vector<std::array<int, 2>>& cavities) {
    int min_layers;
    int min_layers_id;
    int neighbor_layers;

    for(const auto& cavity : cavities) {
        min_layers = rays[cavity[1]].last_layer;
        min_layers_id = cavity[1];

        for(int i=cavity[0]; i!=cavity[1]; i=getNextRayId(i)) {
            if(rays[i].last_layer < min_layers) {
                min_layers = rays[i].last_layer;
                min_layers_id = i;
            }
        }

        neighbor_layers = rays[cavity[0]].last_layer;
        for(int i=getNextRayId(cavity[0]); i!=min_layers_id; i=getNextRayId(i)) {
            if(rays[i].last_layer > neighbor_layers)
                rays[i].last_layer = neighbor_layers;
            else neighbor_layers = rays[i].last_layer;
        }

        neighbor_layers = rays[cavity[1]].last_layer;
        for(int i=getPreviousRayId(cavity[1]); i!=min_layers_id; i=getPreviousRayId(i)) {
            if(rays[i].last_layer > neighbor_layers)
                rays[i].last_layer = neighbor_layers;
            else neighbor_layers = rays[i].last_layer;
        }
    }
}

std::vector<int> BoundaryLayerMesh::pruneRayAABBIntersections(const std::vector<Segment>& test_segments,
                                                              const AABB& target) const {
    vector<int> aabb_intersections;
    for(int i=0; i<test_segments.size(); ++i)
        if(target.containsPortionOf(test_segments[i]))
            aabb_intersections.push_back(i);
    return aabb_intersections;
}

std::vector<std::array<int, 2>> BoundaryLayerMesh::findExtentIntersections(
        const std::vector<Segment>& test_segments, const std::vector<int>& candidates,
        const std::vector<Segment>& target_segments) const {
    vector<std::array<int, 2>> extent_intersections;
    ADT2DExtent adt(max_extent_box.getExtent());
    for(int i=0; i<target_segments.size(); ++i)
        adt.store(i, target_segments[i].getExtent());

    for(int candidate : candidates) {
        auto extents = adt.retrieve(test_segments[candidate].getExtent());
        for(int extent : extents)
            extent_intersections.push_back({candidate, extent});
    }

    return extent_intersections;
}

void BoundaryLayerMesh::clipMultiElementIntersections(const std::vector<std::tuple<int, int, point_t>>& intersections,
                                                      std::vector<Segment>& test_segments, int test_element,
                                                      std::vector<Segment>& target_border, int target_element) {
    int ray_id;
    int boundary_id;
    Segment *ray_segment;
    Segment *boundary;
    Ray *ray;
    array<Ray*, 2> boundary_rays;
    const int test_offset_id = ray_element_start[test_element];
    const int target_offset_id = ray_element_start[target_element];

    point_t intersection;
    point_t clip;

    for(const auto& it : intersections) {
        ray_id = get<0>(it);
        ray = &rays[test_offset_id + ray_id];
        ray_segment = &test_segments[ray_id];
        ray->can_grade = false;

        boundary_id = get<1>(it);
        boundary = &target_border[boundary_id];
        boundary_rays[0] = &rays[target_offset_id + boundary_id];
        boundary_rays[1] = &getNextRay(target_offset_id + boundary_id);
        boundary_rays[0]->can_grade = false;
        boundary_rays[1]->can_grade = false;

        intersection = get<2>(it);

        if(ray_segment->doesIntersect(*boundary)) {
            int layers = (unsigned)distance(layer_offsets.begin()+1,
                                            lower_bound(layer_offsets.begin(), layer_offsets.end(),
                                                        calculateDistance(ray_segment->getPointA(), intersection)) - 2);

            ray->decreaseLayers(layers==ray->last_layer?layers-1:layers);
            ray_segment->setB(ray->pointAtDistance(layer_offsets[ray->last_layer]));
        }
    }
}

void BoundaryLayerMesh::receiveInitialSubdomains() {
    MPI_Request request;
    MPI_Status status;
    int message_available;
    int received;

    vector<vector<double>> extents;

    for(int e=0; e<num_elements; ++e) {
        MPI_Ibarrier(MPI_COMM_WORLD, &request);

        while(true) {
            MPI_Iprobe(MPI_ANY_SOURCE, MessageTag::BL_SUBDOMAIN, MPI_COMM_WORLD, &message_available, &status);
            if(message_available) {
                bl_subdomain_t subdomain(make_shared<BoundaryLayerSubdomain>());
                subdomain->recvSubdomain(status.MPI_SOURCE);
                sort(subdomain->x_vertices.begin(), subdomain->x_vertices.end(), vertexCompareX);
                sort(subdomain->y_vertices.begin(), subdomain->y_vertices.end(), vertexCompareY);
                extents.push_back({subdomain->x_vertices.front().getCoordinate(0),
                                   subdomain->x_vertices.back().getCoordinate(0),
                                   subdomain->y_vertices.front().getCoordinate(1),
                                   subdomain->y_vertices.back().getCoordinate(1)});
                initial_subdomains.push_back(move(subdomain));
                message_available = 0;
            }

            MPI_Test(&request, &received, MPI_STATUS_IGNORE);
            if(received) {
                received = 0;
                break;
            }
        }

        for(auto& extent : extents)
            Send(extent, 4, 0);
        extents.clear();
    }
}

void BoundaryLayerMesh::outputVTK() const {
    ofstream out("output/" + mesher.input_name + ".boundary_layer.vtk");
    out << "# vtk DataFile Version 3.0" << endl;
    out << mesher.input_name << " Boundary Layer PSLG" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;

    out.precision(16);

    out << "POINTS " << vertices.size() << " double" << endl;
    for(const Vertex& vertex : vertices)
        out << vertex.getCoordinate(0) << " " << vertex.getCoordinate(1) << " 0" << endl;

    out << endl << "CELLS " << edges.size() << " " << 3*edges.size() << endl;
    for(const Edge& edge : edges)
        out << "2 " << edge.getVertexId(0) << " " << edge.getVertexId(1) << endl;

    out << endl << "CELL_TYPES " << edges.size() << endl;
    for(long i=0; i<edges.size(); ++i)
        out << "3" << endl;

    out << endl;
    out.close();
}
