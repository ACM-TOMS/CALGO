//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_GEOPRIMITIVES_H
#define FEA_MESHER2D_GEOPRIMITIVES_H

#include "Application.h"

enum class Orientation {
    LEFT,
    ON,
    RIGHT
};

enum Region {
    INSIDE = 0,
    LEFT = 1,
    RIGHT = 2,
    BOTTOM = 4,
    TOP = 8
};

//This set of classes are used as helper objects for performing the boundary layer intersections tests

class AABB;

//The Extent class is used with the ADT when performing the boundary layer intersection checks

class Extent {
public:
    //Default constructor
    Extent() = default;

    //Alternate constructor
    Extent(const point_t& low, const point_t& high);

    //Destructor
    ~Extent() = default;

    bool contains(const Extent& extent) const;

    point_t lo;
    point_t hi;
}; //Extent class


//The Segment class is used for the boundary layer intersection checks

class Segment {
friend class AABB;
public:
    //Default constructor
    Segment() = default;

    //Copy constructor
    Segment(const Segment& s) = default;

    //Alternate constructor
    Segment(point_t p1, point_t p2);

    //Destructor
    ~Segment() = default;

    //Determines if a test point lies on, left, or right of the directed line ab
    Orientation orientation2D(const point_t& t) const;
    bool doesIntersect(const Segment& s) const;

    std::vector<point_t> intersectsAt(const Segment& s) const;

    void setA(point_t p);
    void setB(point_t p);

    point_t getPointA() const;
    point_t getPointB() const;
    Extent getExtent() const;

private:
    point_t a;
    point_t b;
}; //Segment class


//The AABB class is used to prune the search space of candidate rays when checking for boundary layer intersections
//This class is also used for the sizing function for the inviscid region triangles

class AABB {
friend class Segment;
public:
    //Default constructor
    AABB() = default;

    //Alternate constructor
    AABB(point_t lo, point_t hi);

    //Alternate constructor
    AABB(std::array<point_t, 2> bounding_box);

    //Destructor
    ~AABB() = default;

    void setDomain(std::array<point_t, 2> bounding_box);

    //Expands the domain by inflation units in each direction (+x, -x, +y, -y)
    void inflateDomain(double inflation);

    bool intersects(const AABB& other) const;
    bool containsPortionOf(Segment s) const;
    point_t getLowPoint() const;
    point_t getHighPoint() const;
    Extent getExtent() const;

protected:
    int determineCohenSutherlandRegion(point_t p) const;

    point_t low;
    point_t high;
}; //AABB class


#endif //FEA_MESHER2D_GEOPRIMITIVES_H
