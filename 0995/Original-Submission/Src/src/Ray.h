//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_RAY_H
#define FEA_MESHER2D_RAY_H

#include "Application.h"
#include "Edge.h"

class BoundaryLayerMesh;

//The Ray class is responsible for storing information about the normals emitting from the input geometry
//Rays are used to generate the anisotropic boundary layer
//Rays are smoothed, clipped, and grown to create valid and high-fidelity triangles

class Ray {
public:
    //Default constructor - only used to allocate memory for MPI communications
    Ray() = default;

    //Alternate constructor - used when creating the fans at trailing edges
    //Input:
    //end_point - the input geometry vertex that will be the base of the ray
    //normal - the unit vector that points will be inserted along
    //layers - the number of points to insert along normal_vector
    //element_id - the input geometry element that the ray is incident upon
    Ray(const Vertex& end_point, vector_t normal, int layers, int element_id);

    //Alternate constructor - used when creating the initial rays from the input geometry
    //Input:
    //end_point - the input geometry vertex that will be the base of the ray
    //layers - the number of points to insert along normal_vector
    //owning_mesh - the mesh that this ray is a part of, used to calculate normal_vector
    Ray(const Vertex& end_point, int layers, const BoundaryLayerMesh& owning_mesh);

    //Destructor
    ~Ray() = default;

    //Calculates the unique, topological normal that points outwards from the model
    //Input:
    //prev_point - the neighboring point before this ray's endpoint
    //next_point - the neighboring point after this ray's endpoint
    void calculateNormalVector(const point_t& prev_point, const point_t& next_point);

    //Calculates the location of a point inserted along normal_vector from point
    //Input:
    //distance - the distance along normal_vector from point
    point_t pointAtDistance(double distance) const;

    //Bounds checking to make sure the new value of last_layer is less than the previous value, but not less than zero
    //Input:
    //desired - the new value of last_layer
    void decreaseLayers(int desired);

    //Gets the length of normal_vector
    double getMagnitude() const;

    //Gets the vertex id at the requested layer
    long layerVertexID(int layer) const;

    void print() const;

    //The endpoint or base of the ray
    point_t point;

    //The vector where points will be inserted along
    vector_t normal_vector;

    //The mesh id for the surface vertex that this ray is emitted from
    int endpoint_id;

    //The index of the last layer, essentially, the number of layers
    int last_layer;

    //The input geometry element that this ray is incident upon
    int element;

    //The id of the vertex at last_layer of this ray
    long last_vertex_id;

    //Flag to denote if this ray can have more points past last_layer
    bool can_grade = true;
}; //Ray class


#endif //FEA_MESHER2D_RAY_H

