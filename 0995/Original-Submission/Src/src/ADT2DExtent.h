//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_ADT2DEXTENT_H
#define FEA_MESHER2D_ADT2DEXTENT_H

#include <vector>
#include "ADT2D.h"
#include "GeoPrimitives.h"
#include "ADTSpaceTransformer.h"

//The ADT2DExtent class is used so that the user can interact with the ADT without having to worry about the coordinate
//transformations

class ADT2DExtent {
public:
    //No default constructor since we need a domain to create the extent
    ADT2DExtent() = delete;

    //Alternate constructor
    //Input:
    //domain - the 2D bounding box that will be projected to a 4D point and used as the root of the ADT
    ADT2DExtent(const Extent& domain);

    ~ADT2DExtent() = default;

    //Stores a 4D point in unit space in the ADT
    //Input:
    //id - the non-unique identifier for the new element
    //extent - the 2D bounding box that will be projected to a 4D point
    void store(int id, const Extent& extent);

    //Searches the ADT for elements that contain the provided bounding box
    //Input:
    //domain - the 2D bounding box that will be projected to a 4D point and checked for overlaps with other elements
    //Output:
    //vector of ids of the 4D points that overlap
    std::vector<int> retrieve(const Extent& domain) const;

    //Removes the first occurrence of a 4D point in the ADT that overlaps with a provided bounding box
    //Input:
    //id - the identifier of the 4D point to remove
    //extent - the 2D bounding box that will be projected to a 4D point
    void removeFirst(int id, const Extent& extent);

private:
    //Performs the coordinate transformations to unit space or real space
    ADTSpaceTransformer space_transformer;

    //The underlying ADT used to store the elements
    ADT2D adt;
}; //ADT2DExtent class


#endif //FEA_MESHER2D_ADT2DEXTENT_H
