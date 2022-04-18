//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "ADT2DElement.h"
#include "Application.h"

using namespace Application;

ADT2DElement::ADT2DElement(int adt_level, const std::array<double, 4>& x_minimum,
                           const std::array<double, 4>& x_maximum, int element_id, const double* object_coordinates)
        : level(adt_level), id(element_id), x_min(x_minimum), x_max(x_maximum) {
    for(int i=0; i<4; ++i)
        object[i] = object_coordinates[i];
}

ADT2DElement::~ADT2DElement() {
    if(left_child != nullptr)
        delete left_child;
    if(right_child != nullptr)
        delete right_child;
}

bool ADT2DElement::containsHyperRectangle(const std::array<double, 4>& a, const std::array<double, 4>& b) const {
    for(int i=0; i<4; ++i)
        if(b[i] < x_min[i] - precision or a[i] > x_max[i] + precision)
            return false;

    return true;
}

bool ADT2DElement::hyperRectangleContainsObject(const std::array<double, 4>& a, const std::array<double, 4>& b) const {
    for(int i=0; i<4; ++i)
        if(b[i] < object[i] - precision or a[i] > object[i] + precision)
            return false;

    return true;
}
