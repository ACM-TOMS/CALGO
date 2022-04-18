//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_APPLICATION_H
#define FEA_MESHER2D_APPLICATION_H
#include <iostream>
#include <memory>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <assert.h>
#include <float.h>

//Enumerations
enum Axis {
    X_AXIS,
    Y_AXIS
};
enum class Decomposition {
    BOTH,
    LEFT,
    RIGHT,
    NONE,
    ALREADY
};

//Forward declarations
class Vertex;
class SubdomainVertex;
class Edge;
class Segment;
class BoundaryLayerSubdomain;
class InviscidRegionSubdomain;

//Typedefs
typedef std::shared_ptr<BoundaryLayerSubdomain> bl_subdomain_t;
typedef std::shared_ptr<InviscidRegionSubdomain> inviscid_subdomain_t;
typedef std::array<double, 2> vector_t;
typedef std::array<double, 2> point_t;

//The Application namespace is used as set of helper functions to perform floating-point comparisons and basic
//mathematical point and Euclidean vector operations
//std::array<double, 2> is used to represent point and vector types, point_t and vector_t, respectively

namespace Application {
    //Constant used as the floating-point precision for floating-point comparisons
    static const double precision = 1.0e-32;

    //Constant used to convert between radians and degrees
    static const double degree_radian_ratio = 3.1415926535897932/180.0;

    //Returns true if the difference is less than precision
    inline bool areEqual(const double& a, const double& b) {
        return fabs(a - b) < precision;
    }

    //Returns true if value is less than precision
    inline bool isZero(const double& value) {
        return fabs(value) < precision;
    }

    inline bool areEqual(const point_t& a, const point_t& b) {
        return areEqual(a[0], b[0]) and areEqual(a[1], b[1]);
    }

    inline double elapsedMsecs(const timeval& start, const timeval& stop) {
        return ((stop.tv_sec + stop.tv_usec/1000000.0) - (start.tv_sec + start.tv_usec/1000000.0)) * 1000.0;
    }

    inline double calculateDistance(const point_t& a, const point_t& b) {
        return sqrt(pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2));
    }

    inline bool pointRightOfEdge(const point_t& a, const point_t& b, const point_t& t) {
        return ((((a[0]-t[0]) * (b[1]-t[1])) - ((a[1]-t[1]) * (b[0]-t[0]))) > precision);
    }

    std::array<point_t, 2> computeBoundingBox(std::vector<Vertex>::iterator first, std::vector<Vertex>::iterator last);

    inline std::array<point_t, 2> computeBoundingBox(std::vector<point_t>::iterator first,
                                                     std::vector<point_t>::iterator last) {
        double min_x = (*first)[0];
        double max_x = (*first)[0];
        double min_y = (*first)[1];
        double max_y = (*first)[1];

        first++;

        while(first < last) {
            min_x = std::min(min_x, (*first)[0]);
            max_x = std::max(max_x, (*first)[0]);
            min_y = std::min(min_y, (*first)[1]);
            max_y = std::max(max_y, (*first)[1]);

            first++;
        }

        return {{{min_x, min_y}, {max_x, max_y}}};
    }

    bool isClockwise(std::vector<Vertex>::iterator first, std::vector<Vertex>::iterator last);

    inline double degreesToRadians(double angle) {
        return angle*degree_radian_ratio;
    }

    inline double radiansToDegrees(double radians) {
        return radians/degree_radian_ratio;
    }

    inline point_t midPoint(const point_t& a, const point_t& b) {
        return point_t{(a[0]+b[0])/2.0, (a[1]+b[1])/2.0};
    }

    //The control point is treated as the origin
    //The top-right quadrant is 0, the top-left quadrant is 1, the bottom-left quadrant is 2, and
    //the bottom-right quadrant is 3
    inline int relativeQuadrantToControlPoint(const point_t& p, const point_t& control) {
        return ((p[0] > control[0]) ? ((p[1] > control[1]) ? 0 : 3) : ((p[1] > control[1]) ? 1 : 2));
    }

    //Computes the angle between edge ab and the x-axis or y-axis
    inline double angleBetweenEdgeAndAxis(const point_t& a, const point_t& b, bool axis) {
        return (areEqual(a[axis], b[axis]) ? 1.570796326794897 : atan((b[!axis]-a[!axis]) / (b[axis]-a[axis])));
    }

    //The triangle defined by abc should be wound counter-clockwise for a positive area
    inline double triangleArea(const point_t& a, const point_t& b, const point_t& c) {
        return fabs(a[0]*b[1] + b[0]*c[1] + c[0]*a[1] - a[0]*c[1] - c[0]*b[1] - b[0]*a[1])/2.0;
    }

    inline double calculateMagnitude(const vector_t& v) {
        return sqrt(pow(v[0], 2) + pow(v[1], 2));
    }

    inline double dotProduct(const vector_t& a, const vector_t& b) {
        return a[0]*b[0]+a[1]*b[1];
    }

    inline double crossProduct(const vector_t& a, const vector_t& b) {
        return a[0]*b[1]-a[1]*b[0];
    }

    inline double vectorDifference(const vector_t& a, const vector_t& b) {
        return sqrt(pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2));
    }

    inline double angleBetweenVectors(const vector_t& u, const vector_t& v) {
        return radiansToDegrees(atan2(crossProduct(u, v), dotProduct(u,v)));
    }

    //void parallelMergeSortGroup(std::vector<SubdomainVertex>& vertices, Axis axis);
}; //Application namespace


#endif //FEA_MESHER2D_APPLICATION_H
