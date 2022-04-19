//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "GeoPrimitives.h"
#include "Application.h"

using namespace std;
using namespace Application;

Segment::Segment(point_t p1, point_t p2) : a(p1), b(p2) {}

Orientation Segment::orientation2D(const point_t& t) const {
    double area = (((a[0] - t[0]) * (b[1] - t[1])) - ((a[1] - t[1]) * (b[0] - t[0])));

    if(fabs(area) < Application::precision)
        return Orientation::ON;
    else if(area > 0)
        return Orientation::LEFT;
    else return Orientation::RIGHT;
}

AABB::AABB(point_t lo, point_t hi) : low(lo), high(hi) {}

AABB::AABB(std::array<point_t, 2> bounding_box) : low(bounding_box[0]), high(bounding_box[1]) {}

void AABB::setDomain(std::array<point_t, 2> bounding_box) {
    low = bounding_box[0];
    high = bounding_box[1];
}

bool AABB::intersects(const AABB& other) const {
    return (low[0] < other.high[0]) and (low[1] < other.high[1]) and
            (high[0] > other.low[0]) and (high[1] > other.low[1]);
}

int AABB::determineCohenSutherlandRegion(point_t p) const {
    int region = Region::INSIDE;

    if(p[0] < low[0])
        region += Region::LEFT;
    else if(p[0] > high[0])
        region += Region::RIGHT;

    if(p[1] < low[1])
        region += Region::BOTTOM;
    else if(p[1] > high[1])
        region += Region::TOP;

    return region;
}

bool AABB::containsPortionOf(Segment s) const {
    int region_a = determineCohenSutherlandRegion(s.a);
    int region_b = determineCohenSutherlandRegion(s.b);
    double x;
    double y;

    while(true) {
        if(!(region_a | region_b)) {
            return true;
        } else if(region_a & region_b) {
            return false;
        } else {
            int region = region_a ? region_a : region_b;

            if(region & Region::TOP) {
                x = s.a[0] + (s.b[0] - s.a[0]) * (high[1] - s.a[1]) / (s.b[1] - s.a[1]);
                y = high[1];
            } else if(region & Region::BOTTOM) {
                x = s.a[0] + (s.b[0] - s.a[0]) * (low[1] - s.a[1]) / (s.b[1] - s.a[1]);
                y = low[1];
            } else if(region & Region::RIGHT) {
                y = s.a[1] + (s.b[1] - s.a[1]) * (high[0] - s.a[0]) / (s.b[0] - s.a[0]);
                x = high[0];
            } else if(region & Region::LEFT) {
                y = s.a[1] + (s.b[1] - s.a[1]) * (low[0] - s.a[0]) / (s.b[0] - s.a[0]);
                x = low[0];
            }

            if(region == region_a) {
                s.a = {x, y};
                region_a = determineCohenSutherlandRegion(s.a);
            } else {
                s.b = {x, y};
                region_b = determineCohenSutherlandRegion(s.b);
            }
        }
    }
}

point_t AABB::getLowPoint() const {
    return low;
}

point_t AABB::getHighPoint() const {
    return high;
}

point_t Segment::getPointA() const {
    return a;
}

point_t Segment::getPointB() const {
    return b;
}

Extent Segment::getExtent() const {
    return ((a[0] < b[0]) ?
     ((a[1] < b[1]) ? Extent {a, b} : Extent {{a[0], b[1]}, {b[0], a[1]}}) :
     ((a[1] < b[1]) ? Extent {{b[0], a[1]}, {a[0], b[1]}} : Extent {b, a}));
}

vector<point_t> Segment::intersectsAt(const Segment& s) const {
    double dx1 = b[0] - a[0];
    double dx2 = s.b[0] - s.a[0];
    double dy1 = b[1] - a[1];
    double dy2 = s.b[1] - s.a[1];
    double denominator = (dx1 * dy2) - (dy1 * dx2);

    if(isZero(denominator))
        return vector<point_t>();

    double dxa = a[0] - s.a[0];
    double dya = a[1] - s.a[1];

    double numerator = (dya * dx2) - (dxa * dy2);
    double r = numerator/denominator;

    numerator = (dya * dx1) - (dxa * dy1);
    double r2 = numerator/denominator;

    if(((r < precision) or (r > (1.0-precision))) or ((r2 < precision) or (r2 > (1.0-precision))))
        return vector<point_t>();

    return vector<point_t> {{a[0]+(r*dx1), a[1]+(r*dy1)}};
}

bool Segment::doesIntersect(const Segment& s) const {
    array<Orientation, 4> orientations;

    orientations[0] = orientation2D(s.a);
    orientations[1] = orientation2D(s.b);

    if((orientations[0] == Orientation::ON) or (orientations[1] == Orientation::ON))
        return true;
    else if(orientations[0] == orientations[1])
        return false;

    orientations[2] = s.orientation2D(a);
    orientations[3] = s.orientation2D(b);

    if((orientations[2] == Orientation::ON) or (orientations[3] == Orientation::ON))
        return true;
    else if(orientations[2] == orientations[3])
        return false;

    return true;
}

void Segment::setA(point_t p) {
    a = p;
}

void Segment::setB(point_t p) {
    b = p;
}

Extent AABB::getExtent() const {
    return Extent{low, high};
}

void AABB::inflateDomain(double inflation) {
    low[0] -= inflation;
    low[1] -= inflation;
    high[0] += inflation;
    high[1] += inflation;
}

Extent::Extent(const point_t& low, const point_t& high) : lo(low), hi(high) {}

bool Extent::contains(const Extent& extent) const {
    for(int i=0; i<2; ++i) {
        if(extent.hi[i] < lo[i])
            return false;
        if(extent.lo[i] > hi[i])
            return false;
    }

    return true;
}
