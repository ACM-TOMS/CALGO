//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "Ray.h"
#include "BoundaryLayerMesh.h"

using namespace std;
using namespace Application;

Ray::Ray(const Vertex& end_point, vector_t normal, int layers, int element_id)
        : point(end_point.getPoint()), normal_vector(normal), endpoint_id(static_cast<int>(end_point.getID())),
          last_layer(layers), element(element_id), can_grade(false) {}

Ray::Ray(const Vertex& end_point, int layers, const BoundaryLayerMesh& owning_mesh)
        : point(end_point.getPoint()), endpoint_id((int)end_point.getID()), last_layer(layers),
          element(owning_mesh.getElement(end_point.getID())) {
    calculateNormalVector(owning_mesh.getPreviousSurfaceVertex(endpoint_id).getPoint(),
                          owning_mesh.getNextSurfaceVertex(endpoint_id).getPoint());
}

void Ray::calculateNormalVector(const point_t& prev_point, const point_t& next_point) {
    const double tolerance = 1.0;

    double theta = atan(-1*(next_point[0]-prev_point[0])/(next_point[1]-prev_point[1]));
    double vx = cos(theta);
    double vy = sin(theta);

    if(pointRightOfEdge(point, next_point, {point[0]-(vx*tolerance), point[1]-(vy*tolerance)}))
        normal_vector = {-1.0 * vx, -1.0 * vy};
    else normal_vector = {vx, vy};
}

point_t Ray::pointAtDistance(double distance) const {
    return point_t{point[0]-distance*normal_vector[0], point[1]-distance*normal_vector[1]};
}

void Ray::decreaseLayers(int desired) {
    if(desired < 0)
        last_layer = 0;
    else if(desired < last_layer)
        last_layer = desired;
}

double Ray::getMagnitude() const {
    return calculateMagnitude(normal_vector);
}

long Ray::layerVertexID(int layer) const {
    int difference = last_layer - layer;
    assert(difference >= 0);
    return last_vertex_id - difference;
}

void Ray::print() const {
    cout << "Ray " << endpoint_id << " layers (" << last_layer << ") at point (" << point[0] << ", " << point[1]
        << ") with normal <" << normal_vector[0] << ", " << normal_vector[1] << ">" << endl;
}
