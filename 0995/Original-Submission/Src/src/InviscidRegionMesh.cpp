//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "InviscidRegionMesh.h"
#include "MeshGenerator.h"
#include "Application.h"
#include "InviscidRegionSubdomain.h"
#include "MPICommunications.h"
#include <float.h>
#include <climits>
#include <fstream>
#include <iomanip>
#include <stack>
#include <functional>

extern "C" {
#include "triangle/triangle.h"
}

using namespace std;
using namespace Application;
using namespace MPICommunications;

InviscidRegionMesh::InviscidRegionMesh(MeshGenerator& owning_mesher) : mesher(owning_mesher) {}

void InviscidRegionMesh::createInitialSubdomains() {
    initializeInviscidParamters();
    array<vector<Vertex>, 4> nearbody_border(decoupleNearBodyBoundingBox());
    createNearBodySubdomain(nearbody_border);
    createFarfieldSubdomains(nearbody_border);
    mesher.final_edges += static_cast<int>(farfield_edges.size());
}

void InviscidRegionMesh::initializeInviscidParamters() {
    next_vertex_id = mesher.boundary_layer.max_boundary_layer_vertex_id;
    next_edge_id = mesher.boundary_layer.next_edge_id;
    fast_growth *= mesher.boundary_layer.chord_length;
    farfield *= mesher.boundary_layer.chord_length;
    isotropic_area = 4.0 * mesher.boundary_layer.getAverageEnclosingTriangleArea();
    enlargeNearBodyBoundingBox();
}

void InviscidRegionMesh::synchronizeInviscidParameters() {
    vector<double> params;
    if(mesher.rank == 0) {
        for(point_t &point : aabb_centers) {
            params.push_back(point[0]);
            params.push_back(point[1]);
        }
        params.push_back(isotropic_area);
    }

    const int num_elements = mesher.boundary_layer.num_elements;
    Broadcast(params, num_elements*2 + 1, 0);
    if(mesher.rank != 0) {
        isotropic_area = params.back();
        for(int e=0; e<num_elements; ++e)
            aabb_centers.emplace_back(point_t{{params[e*2], params[e*2+1]}});
    }

    triangle_aabb_centers.insert(triangle_aabb_centers.end(), params.begin(), params.end() - 1);

    array<long, 2> start_and_step = array<long, 2> ();
    if(mesher.rank == 0)
        start_and_step = {next_vertex_id, (LONG_MAX-next_vertex_id) / mesher.processes};
    MPI_Bcast(start_and_step.data(), 2, MPI_LONG, 0, MPI_COMM_WORLD);
    if(mesher.rank != 0)
        next_vertex_id = start_and_step[0] + mesher.rank * start_and_step[1];
}

void InviscidRegionMesh::createNearBodySubdomain(const std::array<std::vector<Vertex>, 4>& sides) {
    inviscid_subdomain_t subdomain = make_shared<InviscidRegionSubdomain>(InviscidRegionSubdomain(true));

    for(const auto& side : sides)
        subdomain->vertices.insert(subdomain->vertices.end(), side.begin(), side.end());

    subdomain->createBorder();

    for(const auto& bl_border : mesher.boundary_layer.transition_vertices)
        transform(bl_border.cbegin(), bl_border.cend(), back_inserter(subdomain->vertices),
                  [](const Vertex* vertex) {return *vertex;});

    transform(mesher.boundary_layer.edges.cbegin() + mesher.boundary_layer.num_model_edges,
              mesher.boundary_layer.edges.cbegin() + mesher.boundary_layer.num_model_edges +
                      mesher.boundary_layer.num_enclosing_edges,
              back_inserter(subdomain->edges), [](const Edge &e) {return e.getVertices();});

    computeSubdomainCost(subdomain);
    mesher.addSubdomainToQueue(subdomain);
}

std::vector<Vertex> InviscidRegionMesh::createAxisAlignedDecouplingPath(Vertex first, Axis axis, double last) {
    vector<Vertex> path(1, first);
    point_t point = first.getPoint();
    function<bool ()> progress;
    function<void (double)> increment;

    if(point[axis] < last) {
        progress = [&]() {return point[axis]<last;};
        increment = [&](double v) {point[axis]+=v;};
    } else {
        progress = [&]() {return point[axis]>last;};
        increment = [&](double v) {point[axis]-=v;};
    }

    while(progress()) {
        increment(getEdgeConstraintAtPoint(point));
        path.emplace_back(point, true, next_vertex_id++);
    }

    if(not path.empty()) {
        if(fabs(path[path.size() - 2].getCoordinate(axis) - last) < fabs(path.back().getCoordinate(axis) - last))
            path.pop_back();

        if(path.size() > 1)
            path.pop_back();
    }

    return path;
}

double InviscidRegionMesh::getEdgeConstraintAtPoint(const point_t &point) {
    return sqrt(getAreaConstraintAtPoint(point)/sqrt(2.0));
}

double InviscidRegionMesh::getAreaConstraintAtPoint(const point_t &point) {
    if(uniform)
        return isotropic_area;

    double distance = DBL_MAX;
    for(const point_t& center_point : aabb_centers)
        distance = min(distance, calculateDistance(point, center_point));

    if(distance < fast_growth)
        return isotropic_area*pow(1.0+distance, 4);
    else return pow(distance, 2)*isotropic_area*pow(1.0+fast_growth, 2);
}

void InviscidRegionMesh::setNearBodyBoundingBox(std::vector<std::array<point_t, 2>> &extents) {
    std::array<point_t, 2> nearbody_extent {{{{DBL_MAX, DBL_MAX}}, {{-DBL_MAX, -DBL_MAX}}}};

    for(const auto& extent : extents) {
        aabb_centers.push_back(midPoint(extent[0], extent[1]));
        nearbody_extent[0][0] = min(extent[0][0], nearbody_extent[0][0]);
        nearbody_extent[1][0] = max(extent[1][0], nearbody_extent[1][0]);
        nearbody_extent[0][1] = min(extent[0][1], nearbody_extent[0][1]);
        nearbody_extent[1][1] = max(extent[1][1], nearbody_extent[1][1]);
    }

    nearbody_box[0][0] = nearbody_extent[0][0];
    nearbody_box[0][1] = nearbody_extent[0][1];
    nearbody_box[1][0] = nearbody_extent[1][0];
    nearbody_box[1][1] = nearbody_extent[0][1];
    nearbody_box[2][0] = nearbody_extent[1][0];
    nearbody_box[2][1] = nearbody_extent[1][1];
    nearbody_box[3][0] = nearbody_extent[0][0];
    nearbody_box[3][1] = nearbody_extent[1][1];
}

void InviscidRegionMesh::enlargeNearBodyBoundingBox() {
    const double gap = 4.0*sqrt(2.0*isotropic_area);
    nearbody_box[0][0] -= gap;
    nearbody_box[0][1] -= gap;
    nearbody_box[1][0] += gap;
    nearbody_box[1][1] -= gap;
    nearbody_box[2][0] += gap;
    nearbody_box[2][1] += gap;
    nearbody_box[3][0] -= gap;
    nearbody_box[3][1] += gap;
}

void InviscidRegionMesh::translateNearBodyBoundingBox() {
    const point_t& center = mesher.boundary_layer.mesh_center;
    nearbody_box[0][0] += center[0];
    nearbody_box[0][1] += center[1];
    nearbody_box[1][0] += center[0];
    nearbody_box[1][1] += center[1];
    nearbody_box[2][0] += center[0];
    nearbody_box[2][1] += center[1];
    nearbody_box[3][0] += center[0];
    nearbody_box[3][1] += center[1];
}

void InviscidRegionMesh::setFarFieldDistance(double chord_lengths) {
    farfield = chord_lengths;
}

std::array<std::vector<Vertex>, 4> InviscidRegionMesh::decoupleNearBodyBoundingBox() {
    return array<vector<Vertex>, 4>
            {{
                     createAxisAlignedDecouplingPath(Vertex(nearbody_box[0], true, next_vertex_id++), X_AXIS,
                                                     nearbody_box[1][0]),
                     createAxisAlignedDecouplingPath(Vertex(nearbody_box[1], true, next_vertex_id++), Y_AXIS,
                                                     nearbody_box[2][1]),
                     createAxisAlignedDecouplingPath(Vertex(nearbody_box[2], true, next_vertex_id++), X_AXIS,
                                                     nearbody_box[3][0]),
                     createAxisAlignedDecouplingPath(Vertex(nearbody_box[3], true, next_vertex_id++), Y_AXIS,
                                                     nearbody_box[0][1])
             }};
}

void InviscidRegionMesh::createFarfieldSubdomains(const array<vector<Vertex>, 4>& interior) {
    point_t center = mesher.boundary_layer.mesh_center;

    array<Vertex, 4> shared_corners
            {{
                     Vertex({{interior[0].back().getCoordinate(0), -farfield + center[1]}}, true, next_vertex_id++),
                     Vertex({{farfield + center[0], interior[1].back().getCoordinate(1)}}, true, next_vertex_id++),
                     Vertex({{interior[2].back().getCoordinate(0), farfield + center[1]}}, true, next_vertex_id++),
                     Vertex({{-farfield + center[0], interior[3].back().getCoordinate(1)}}, true, next_vertex_id++)
             }};

    const array<array<Axis, 2>, 2> axis {{{{X_AXIS, Y_AXIS}}, {{Y_AXIS, X_AXIS}}}};
    const array<point_t, 4> chord {{
                                           {{-farfield + center[0], -farfield + center[1]}},
                                           {{farfield + center[0], -farfield + center[1]}},
                                           {{farfield + center[0], farfield + center[1]}},
                                           {{-farfield + center[0], farfield + center[1]}}
                                   }};
    vector<Vertex> path;
    vector<Vertex> prev_path;
    vector<Vertex> last_path;
    long last_id;
    long previous_id = -1;
    inviscid_subdomain_t subdomain;

    for(int i=0; i<4; i++) {
        subdomain = make_shared<InviscidRegionSubdomain>();
        subdomain->vertices.insert(subdomain->vertices.end(), interior[i].rbegin(), interior[i].rend()-1);

        if(i == 0) {
            path = createAxisAlignedDecouplingPath(interior[i].front(), axis[i%2][0], chord[i][axis[i%2][0]]);
            subdomain->vertices.insert(subdomain->vertices.end(), path.begin(), path.end());
            last_path = move(path);
        } else {
            subdomain->vertices.insert(subdomain->vertices.end(), prev_path.begin(), prev_path.end());
        }

        path = createAxisAlignedDecouplingPath(shared_corners[(i+4-1)%4], axis[i%2][1], chord[i][axis[i%2][1]]);
        path.emplace_back(point_t{{chord[i][0], chord[i][1]}}, true, next_vertex_id++);
        subdomain->vertices.insert(subdomain->vertices.end(), path.begin(), path.end());
        if(i == 0)
            last_id = path.front().getID();
        previous_id = addFarfieldEdges(previous_id, path.begin(), path.end()-1);

        path = createAxisAlignedDecouplingPath(shared_corners[i], axis[i%2][0], chord[i][axis[i%2][0]]);
        reverse(path.begin(), path.end());
        subdomain->vertices.insert(subdomain->vertices.end(), path.begin(), path.end());
        previous_id = addFarfieldEdges(previous_id, path.begin(), path.end()-1);

        if(i == 3) {
            subdomain->vertices.insert(subdomain->vertices.end(), last_path.rbegin(), last_path.rend());
        } else {
            path = createAxisAlignedDecouplingPath(interior[(i+1)%4].front(), axis[i%2][1], chord[i][axis[i%2][1]]);
            subdomain->vertices.insert(subdomain->vertices.end(), path.rbegin(), path.rend());
            prev_path = move(path);
        }

        computeSubdomainCost(subdomain);
        if(mesher.processes == 1 or i % mesher.processes == 0)
            initial_subdomains.push(move(subdomain));
        else subdomain->sendSubdomain(i % mesher.processes);
    }

    farfield_edges.emplace_back(next_edge_id++, previous_id, last_id, 2);
}

long InviscidRegionMesh::addFarfieldEdges(long first_id, vector<Vertex>::iterator next, vector<Vertex>::iterator stop) {
    if(first_id != -1)
        farfield_edges.emplace_back(next_edge_id++, first_id, next->getID(), 2);

    for(; next<stop; ++next)
        farfield_edges.emplace_back(next_edge_id++, next->getID(), (next+1)->getID(), 2);

    return stop->getID();
}

std::vector<Vertex> InviscidRegionMesh::createDecouplingPath(Vertex& first, point_t last) {
    vector<Vertex> border;
    array<double, 2> total_vec;
    array<double, 2> vec;
    point_t point;
    bool reversed;
    double distance;
    double length;
    double ratio;

    if(getEdgeConstraintAtPoint(last) < getEdgeConstraintAtPoint(first.getPoint())) {
        point = last;
        total_vec[0] = first.getCoordinate(0) - point[0];
        total_vec[1] = first.getCoordinate(1) - point[1];
        reversed = true;
    } else {
        point = first.getPoint();
        total_vec[0] = last[0] - point[0];
        total_vec[1] = last[1] - point[1];
        border.push_back(first);
        reversed = false;
    }

    distance = calculateMagnitude(total_vec);

    while(distance > 0) {
        length = getEdgeConstraintAtPoint(point);
        ratio = length/distance;
        vec[0] = ratio*total_vec[0];
        vec[1] = ratio*total_vec[1];
        point[0] += vec[0];
        point[1] += vec[1];
        border.emplace_back(point, true, next_vertex_id++);

        distance -= length;
        total_vec[0] -= vec[0];
        total_vec[1] -= vec[1];
    }

    if(reversed) {
        if(fabs(calculateDistance(border[border.size()-2].getPoint(), first.getPoint())) <
           fabs(calculateDistance(border.back().getPoint(), first.getPoint())))
            border.pop_back();

        border.pop_back();

        border.push_back(first);
        reverse(border.begin(), border.end());
    } else {
        if(fabs(calculateDistance(border[border.size()-2].getPoint(), last)) <
           fabs(calculateDistance(border.back().getPoint(), last)))
            border.pop_back();

        border.pop_back();
    }

    return border;
}

void InviscidRegionMesh::decoupleInitialSubdomains() {
    if(mesher.processes > 4) {
        if(uniform)
            distributeInitialUniformSubdomains();
        else distributeInitialGradedSubdomains();
    }

    overDecoupleLocalSubdomains();
}

const std::vector<double>& InviscidRegionMesh::getBoundaryLayerHoles() const {
    return mesher.boundary_layer.getModelHoles();
}

std::tuple<bool, std::array<inviscid_subdomain_t, 4>> InviscidRegionMesh::splitInviscidSubdomain(
        inviscid_subdomain_t& subdomain) {
    const point_t joint = subdomain->getCenterPoint();
    bool can_split;
    array<vector<Vertex>::iterator, 4> split_iterators;
    tie(can_split, split_iterators) = determineSplitVertices(subdomain, joint);
    if(not can_split)
        return make_tuple(false, array<inviscid_subdomain_t, 4>());

    array<inviscid_subdomain_t, 4> subs(partitionParentVertices(subdomain, split_iterators));

    auto borders = createSplitBorders(split_iterators, joint);
    for(int s=0; s<4; ++s) {
        subs[s]->vertices.insert(subs[s]->vertices.end(), borders[(s+1)%4].rbegin(), borders[(s+1)%4].rend()-1);
        subs[s]->vertices.insert(subs[s]->vertices.end(), borders[s].begin(), borders[s].end());
        computeSubdomainCost(subs[s]);
    }

    return make_tuple(true, subs);
}

std::tuple<bool, std::array<std::vector<Vertex>::iterator, 4>> InviscidRegionMesh::determineSplitVertices(
        inviscid_subdomain_t& subdomain, const point_t& joint) {
    array<vector<Vertex>::iterator, 4> split_iterators;
    int prev_quadrant = relativeQuadrantToControlPoint(subdomain->vertices.back().getPoint(), joint);
    int quadrant;
    int visited_quadrants = 0;

    for(auto it=subdomain->vertices.begin(); it<subdomain->vertices.end(); ++it) {
        quadrant = relativeQuadrantToControlPoint((*it).getPoint(), joint);

        if(quadrant == ((prev_quadrant+1)%4)) {
            if(angleBetweenEdgeAndAxis((*(it-1)).getPoint(), joint, quadrant%2) <
               angleBetweenEdgeAndAxis((*it).getPoint(), joint, quadrant%2)) {
                split_iterators[quadrant] = it-1;
            } else {
                split_iterators[quadrant] = it;
            }

            if(++visited_quadrants == 4)
                return make_tuple(true, split_iterators);
        }

        prev_quadrant = quadrant;
    }

    return make_tuple(false, split_iterators);
}

std::array<std::vector<Vertex>, 4> InviscidRegionMesh::createSplitBorders(
        const std::array<std::vector<Vertex>::iterator, 4>& splits, const point_t& sub_center) {
    array<vector<Vertex>, 4> borders;
    Vertex center_vertex(sub_center, true, next_vertex_id++);
    for(int b=0; b<4; b++)
        borders[b] = createDecouplingPath(center_vertex, splits[b]->getPoint());
    return borders;
}

void InviscidRegionMesh::globalizeTriangleOutput(struct triangulateio& out, std::vector<long>& global_ids) {
    for(int i=(int)global_ids.size(); i<out.numberofpoints; ++i)
        global_ids.push_back(next_vertex_id++);
    mesher.addMeshToOutbox(out, global_ids);
}

void InviscidRegionMesh::computeSubdomainCost(inviscid_subdomain_t& s) {
    auto extent = s->getExtentBox();
    double area = 0.25 * (extent[1][0] - extent[0][0]) * (extent[1][1] - extent[0][1]);
    if(uniform) {
        s->cost = (int)ceil(area/isotropic_area);
    } else {
        point_t center = midPoint(extent[0], extent[1]);

        s->cost = approximateTriangleCount(extent[0], center, area) +
                  approximateTriangleCount({{extent[1][0], extent[0][1]}}, center, area) +
                  approximateTriangleCount(extent[1], center, area) +
                  approximateTriangleCount({{extent[0][0], extent[1][1]}}, center, area);
    }
}

int InviscidRegionMesh::approximateTriangleCount(const point_t& a, const point_t& b, double area) {
    point_t center = midPoint(a, b);
    return (int)(((area/getAreaConstraintAtPoint(midPoint(a, center))) +
                  (area/getAreaConstraintAtPoint(midPoint({{b[0], a[1]}}, center))) +
                  (area/getAreaConstraintAtPoint(midPoint(b, center))) +
                  (area/getAreaConstraintAtPoint(midPoint({{a[0], b[1]}}, center))))/4);
}

std::array<inviscid_subdomain_t, 4> InviscidRegionMesh::partitionParentVertices(
        const inviscid_subdomain_t& parent, const std::array<std::vector<Vertex>::iterator, 4>& spliterators) {
    array<inviscid_subdomain_t, 4> subdomains{{make_shared<InviscidRegionSubdomain>(),
                                                      make_shared<InviscidRegionSubdomain>(),
                                                      make_shared<InviscidRegionSubdomain>(),
                                                      make_shared<InviscidRegionSubdomain>()}};
    for(int s=0; s<4; ++s) {
        if(spliterators[s] < spliterators[(s+1)%4]) {
            for(auto it=spliterators[s]; it<=spliterators[(s+1)%4]; ++it)
                subdomains[s]->vertices.push_back(*it);
        } else {
            for(auto it=spliterators[s]; it<parent->vertices.end(); ++it)
                subdomains[s]->vertices.push_back(*it);
            for(auto it=parent->vertices.begin(); it<=spliterators[(s+1)%4]; ++it)
                subdomains[s]->vertices.push_back(*it);
        }
    }

    return subdomains;
}

void InviscidRegionMesh::receiveInitialSubdomains() {
    for(int i=mesher.rank; i<4; i+=mesher.processes) {
        inviscid_subdomain_t subdomain(make_shared<InviscidRegionSubdomain>());
        subdomain->recvSubdomain(0);
        initial_subdomains.push(subdomain);
    }
}

void InviscidRegionMesh::useUniformTriangles() {
    uniform = true;
}

void InviscidRegionMesh::distributeInitialUniformSubdomains() {
    inviscid_subdomain_t subdomain;
    array<inviscid_subdomain_t, 4> split_subdomains;
    bool split_successful;
    const int iterations = static_cast<int>(ceil(log(mesher.processes)/log(4)));

    int increment;
    for(int i=1; i<iterations; i++) {
        increment = static_cast<int>(pow(4, i));

        if(mesher.rank < increment) {
            do {
                if(initial_subdomains.empty())
                    throw logic_error("Producer has no subdomains. There are too many processes for a small domain");
                subdomain = move(initial_subdomains.top());
                initial_subdomains.pop();
                tie(split_successful, split_subdomains) = splitInviscidSubdomain(subdomain);
                if(split_successful) {
                    initial_subdomains.push(split_subdomains[0]);

                    for(int p=1; p<4; ++p) {
                        if(mesher.rank + p*increment < mesher.processes)
                            split_subdomains[p]->sendSubdomain(mesher.rank + p*increment);
                        else if(mesher.rank + (p-1)*increment < mesher.processes)
                            split_subdomains[p]->sendSubdomain(mesher.rank + (p-1)*increment);
                        else initial_subdomains.push(split_subdomains[p]);
                    }
                } else {
                    mesher.addSubdomainToQueue(subdomain);
                }
            } while(not split_successful);
        } else if(mesher.rank < 4*increment) {
            int source = mesher.rank % increment;
            subdomain = make_shared<InviscidRegionSubdomain>();
            subdomain->recvSubdomain(source);
            initial_subdomains.push(subdomain);

            if(source + 2*increment >= mesher.processes) {
                subdomain = make_shared<InviscidRegionSubdomain>();
                subdomain->recvSubdomain(source);
                initial_subdomains.push(subdomain);
            }
        }
    }
}

void InviscidRegionMesh::distributeInitialGradedSubdomains() {
    inviscid_subdomain_t subdomain;
    array<inviscid_subdomain_t, 4> split_subdomains;
    bool split_successful;
    const int iterations = static_cast<int>(ceil(log(mesher.processes)/log(4)));

    if(mesher.rank < 4) {
        for(int i=1; i<=iterations; ++i) {
            do {
                if(initial_subdomains.empty())
                    throw logic_error("Producer has no subdomains. There are too many processes for a small domain");
                subdomain = move(initial_subdomains.top());
                initial_subdomains.pop();
                tie(split_successful, split_subdomains) = splitInviscidSubdomain(subdomain);
                if(split_successful) {
                    for(auto&& sub : split_subdomains)
                        initial_subdomains.push(move(sub));
                } else {
                    mesher.addSubdomainToQueue(subdomain);
                }
            } while(not split_successful);
        }

        for(int p=mesher.rank+4; p<mesher.processes; p+=4) {
            initial_subdomains.top()->sendSubdomain(p);
            initial_subdomains.pop();
        }
    } else {
        subdomain = make_shared<InviscidRegionSubdomain>();
        subdomain->recvSubdomain(mesher.rank % 4);
        initial_subdomains.push(subdomain);
    }
}

void InviscidRegionMesh::overDecoupleLocalSubdomains() {
    inviscid_subdomain_t subdomain;
    array<inviscid_subdomain_t, 4> split_subdomains;
    bool split_successful;

    while(not initial_subdomains.empty()) {
        subdomain = move(initial_subdomains.top());
        initial_subdomains.pop();

        tie(split_successful, split_subdomains) = splitInviscidSubdomain(subdomain);
        if(split_successful) {
            for(inviscid_subdomain_t& split : split_subdomains)
                if(split->cost < decoupling_work_threshold)
                    mesher.addSubdomainToQueue(split);
                else initial_subdomains.push(move(split));
        } else {
            mesher.addSubdomainToQueue(subdomain);
        }
    }
}
