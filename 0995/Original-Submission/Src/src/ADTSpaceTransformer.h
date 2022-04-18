//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_ADTSPACETRANSFORMER_H
#define FEA_MESHER2D_ADTSPACETRANSFORMER_H

#include "GeoPrimitives.h"

//The ADTSpaceTransformer class is used as a helper class with the ADT2D class
//This class provides the coordinate transformations from unit space to real space and real space to unit space

class ADTSpaceTransformer {
public:
    //No default constructor because a real domain is needed
    ADTSpaceTransformer() = delete;

    //Alternate constructor
    //Input:
    //domain - the 2D bounding box that will represent the real domain
    ADTSpaceTransformer(const Extent& domain)
            : extent(domain), scale(std::max(extent.hi[0] - extent.lo[0], extent.hi[1] - extent.lo[1])),
              over_scale(1.0/scale) {}

    //Returns the domain in real space
    inline const Extent& getDomain() const {
        return extent;
    }

    inline point_t toUnitSpace(const point_t& real_point) const {
        point_t unit = real_point;
        unit[0] -= extent.lo[0];
        unit[1] -= extent.lo[1];
        unit[0] *= over_scale;
        unit[1] *= over_scale;
        return unit;
    }

    inline point_t toRealSpace(const point_t& unit_point) const {
        point_t real = unit_point;
        real[0] *= scale;
        real[1] *= scale;
        real[0] += extent.lo[0];
        real[1] += extent.lo[1];
        return real;
    }

private:
    //The real domain
    Extent extent;

    //Used to convert from unit space to real space
    //The length of the longest side of the domain
    double scale;

    //Used to convert from real space to unit space
    //Equal to 1/scale
    double over_scale;
}; //ADTSpaceTransformer class


#endif //FEA_MESHER2D_ADTSPACETRANSFORMER_H
