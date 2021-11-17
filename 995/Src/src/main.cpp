#include <iostream>
#include "MPICommunications.h"
#include "MeshGenerator.h"

using namespace std;
using namespace MPICommunications;

int naca0012_main();
int wright1903_main();
int airfoil30p30n_main();

int main(int argc, char* argv[]) {
    initialize();

    airfoil30p30n_main();

    finalize();
    return 0;
}

//Read the surface description in the file naca0012.poly in the "geometry" folder
//Set the growth function used to generate the boundary layer
//Set the distance to the outer edges (farfield) of the domain
//Specify to use the same sized triangles for the entirety of the inviscid region
//Mesh the domain, collect the resulting mesh, and save the file in the VTK format in the "output" folder
int naca0012_main() {
    MeshGenerator mesher;
    mesher.readInputModel("geometry/naca0012.poly");
    mesher.setBoundaryLayerGrowthFunction(0.0001, 1.15, 30);
    mesher.setFarFieldDistance(2.0);
    mesher.useUniformTrianglesForInviscidRegion();
    mesher.meshDomain();
    mesher.collectFinalMesh();

    if(myRank() == 0)
        mesher.outputVTK();

    return 0;
}

//Read the surface description in the file wright1903.poly in the "geometry" folder
//Set the growth function used to generate the boundary layer
//Set the distance to the outer edges (farfield) of the domain
//Mesh the domain, collect the resulting mesh, and save the file in the VTK format in the "output" folder
int wright1903_main() {
    MeshGenerator mesher;
    mesher.readInputModel("geometry/wright1903.poly");
    mesher.setBoundaryLayerGrowthFunction(0.00005, 1.10, 50);
    mesher.setFarFieldDistance(30.0);
    mesher.meshDomain();
    mesher.collectFinalMesh();

    if(myRank() == 0)
        mesher.outputVTK();

    return 0;
}

//Read the surface description in the file airfoil_30p30n.poly in the "geometry" folder
//Set the growth function used to generate the boundary layer
//Set the distance to the outer edges (farfield) of the domain
//Mesh the domain, collect the resulting mesh, and save the file in the VTK format in the "output" folder
int airfoil30p30n_main() {
    MeshGenerator mesher;
    mesher.readInputModel("geometry/airfoil30p30n.poly");
    mesher.setBoundaryLayerGrowthFunction(0.000075, 1.10, 30);
    mesher.setFarFieldDistance(10.0);
    mesher.meshDomain();
    mesher.collectFinalMesh();

    if(myRank() == 0)
        mesher.outputVTK();

    return 0;
}
