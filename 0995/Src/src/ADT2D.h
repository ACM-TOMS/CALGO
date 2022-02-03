//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_ADT2D_H
#define FEA_MESHER2D_ADT2D_H

#include "ADT2DElement.h"
#include <vector>
#include <tuple>

//The ADT2D class is used to efficiently determine if there are intersections while generating the boundary layer
//Each ray segment or border segment has their axis-aligned bounding box projected to a 4-dimensional halfspace
//The 4D points are stored in an alternating digital tree (ADT) and recursive searching is used to check for overlaps
//in the axis-aligned bounding boxes, or extent boxes

//More information:
//An Alternating Digial Tree (ADT) Algorithm for 3D Geometric Searching and Intersection Problems.
//Jaime Peraire, Javier Bonet
//International Journal for Numerical Methods In Eng, vol 31, 1-17, (1991)

class ADT2D {
public:
    //Default constructor
    //Sets root to nullptr
    ADT2D() = default;

    //No copy constructor since we are using raw pointers
    ADT2D(const ADT2D&) = delete;

    //Move constructor
    ADT2D(ADT2D&& other);

    //Destructor, cascades deletion of all ADT2DElements in the tree by starting with the root
    ~ADT2D();

    //Stores an extent box as a 4D point
    //Input:
    //id - non-unique identifier
    //x - starting memory address for the 4D point
    void store(int id, const double* x);

    //Retrieves the 4D points that overlap with a provided extent box
    //Input:
    //extent - starting memory address for the testing extent box
    //Output:
    //vector of ids of the 4D points that overlap
    std::vector<int> retrieve(const double* extent) const;

    //Removes the first occurrence of a 4D point in the tree that overlaps with a provided extent box
    //Input:
    //id - the identifier of the 4D point to remove
    //extent - starting memory address for the testing extent box
    void removeFirst(int id, const double* extent);

private:
    //Root element of the tree, default initialized to nullptr
    ADT2DElement* root = nullptr;

    //-----------------------Private Functions-------------------------

    void store(ADT2DElement* element, int id, const double* x);
    void retrieve(ADT2DElement* element, std::vector<int>& ids, const std::array<double, 4>& a,
                  const std::array<double, 4>& b) const;
    void removeFirst(ADT2DElement* element, ADT2DElement* parent, int id, const std::array<double, 4>& a,
                     const std::array<double, 4>& b);

    std::tuple<std::array<double, 4>, std::array<double, 4>> createHyperRectangle(const double* extent) const;
    int determineChild(ADT2DElement* element, double const* x);
    void replaceElementWithLeaf(ADT2DElement* element, ADT2DElement* parent);
    std::tuple<ADT2DElement*, ADT2DElement*> findFirstLeaf(ADT2DElement* element, ADT2DElement* parent);
    bool isLeaf(const ADT2DElement* element) const;
}; //class ADT2D


#endif //FEA_MESHER2D_ADT2D_H
