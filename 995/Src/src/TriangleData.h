//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_TRIANGLEDATA_H
#define FEA_MESHER2D_TRIANGLEDATA_H

#include <vector>

//The TriangleData struct is responsible for holding the output from calls to Triangle

struct TriangleData {
    TriangleData() = delete;

    //Alternate Constructor
    //Acquires ownership of vertices_in, triangles_in, and global_ids_in and sets the arguments to a null state
    TriangleData(double*& vertices_in, int num_triangles_in, int*& triangles_in, std::vector<long>& global_ids_in)
            : vertices(vertices_in), num_triangles(num_triangles_in), triangles(triangles_in),
              global_ids(std::move(global_ids_in)) {
        vertices_in = 0;
        triangles_in = 0;
    }

    //No copy constructor because we are using raw pointers
    TriangleData(const TriangleData& td) = delete;

    //Move Constructor
    TriangleData(TriangleData&& td)
            : vertices(td.vertices), num_triangles(td.num_triangles), triangles(td.triangles),
              global_ids(std::move(td.global_ids)) {
        td.vertices = 0;
        td.triangles = 0;
    }

    //Destructor
    ~TriangleData() {
        if(vertices != NULL)
            free(vertices);
        if(triangles != NULL)
            free(triangles);
    }

    //Starting memory address for the coordinates of the mesh vertices
    //Length is equal to twice the size of global_ids
    double* vertices;

    int num_triangles;

    //Starting memory address for the endpoints of the mesh triangles
    //Length is equal to three times num_triangles
    int* triangles;

    //Unique global identifiers for the vertices
    std::vector<long> global_ids;
}; //struct TriangleData


#endif //FEA_MESHER2D_TRIANGLEDATA_H
