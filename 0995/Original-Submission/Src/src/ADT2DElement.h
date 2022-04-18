//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_ADT2DELEMENT_H
#define FEA_MESHER2D_ADT2DELEMENT_H

#include <array>

//The ADT2DElement class is used to store individual nodes of the alternating digital tree

class ADT2DElement {
public:
    //No default constructor because all ADT2DElements should have a user-assigned id
    ADT2DElement() = delete;

    //Constructor
    ADT2DElement(int adt_level, const std::array<double, 4>& x_minimum, const std::array<double, 4>& x_maximum,
                 int element_id, const double* object_coordinates);

    //Destructor
    ~ADT2DElement();

    //Returns true if the hyper-rectangle is inside this element's domain
    //Input:
    //a - the lower point of the hyper-rectangle
    //b - the upper point of the hyper-rectangle
    bool containsHyperRectangle(const std::array<double, 4>& a, const std::array<double, 4>& b) const;

    //Returns true if the user-provided 4D point lies inside the hyper-rectangle
    //Input:
    //a - the lower point of the hyper-rectangle
    //b - the upper point of the hyper-rectangle
    bool hyperRectangleContainsObject(const std::array<double, 4>& a, const std::array<double, 4>& b) const;

    //Used to determine which child to pick
    int level;

    //Non-unique identifier
    int id;

    ADT2DElement* left_child = nullptr;
    ADT2DElement* right_child = nullptr;

    //The lower bound of this element's domain
    std::array<double, 4> x_min;

    //The upper bound of this element's domain
    std::array<double, 4> x_max;

    //The 4D point that the user provided
    std::array<double, 4> object;
}; //ADT2DElement class


#endif //FEA_MESHER2D_ADT2DELEMENT_H
