FEA_Mesher2D - By Juliette Pardue

This is a parallel, triangular mesh generator that is capable of generating boundary layer meshes and isotropic mesh regions suitable for finite element analysis.
FEA_Mesher2D generates a high-fidelity, anisotropic boundary layer mesh from a user-defined growth function, generates a globally Delaunay, graded, isotropic mesh region in parallel, resolves potential interpolation errors in the boundary layer caused by the local mesh density, resolves self intersections and multi-element intersections in the boundary layer, is a push-button mesh generator so no human interaction is required after startup, and is scalable and efficient.

The class MeshGenerator is the object that the end user will interact with. main.cpp shows an example of the calls.
The functions the end user should be concerned with are
-readInputModel(string filename)
-setBoundaryLayerGrowthFunction(double first_thickness, double growth_rate, int initial_layers)
-setFarFieldDistance(double chord_lengths)
-useUniformTrianglesForInviscidRegion() [optional]
-meshDomain()
-collectFinalMesh()

Provided output options
outputTecplot() - used for Tecplot
outputShowMe()  - used for Shewchuk's "Show Me" application
outputFUN3D()   - .msh file used for the flow solver FUN3D
outputVTK()     - .vtk format, can be used with ParaView

Included Input Files in "geometry" folder
naca0012.poly       - Standard benchmark airfoil
wright1903.poly     - Wright Brother's 1903 airfoil used for the first flight
airfoil_30p30n.poly - Complex airfoil with three elements

Suggested Input Parameters
Input file          - First Layer Thickness, Layer Thickness Growth Rate, Number of Initial Layers, Farfield Distance, Uniform Inviscid
naca0012.poly       - 0.0001, 1.15, 30, 2.0, Yes
wright1903.poly     - 0.00005, 1.10, 50, 30.0, No
airfoil_30p30n.poly - 0.000075, 1.10, 30, 50.0, No

The “output” folder is where mesh files will be written to, and the output meshes in the VTK format of the suggested input parameters of the three provided geometries are located in the “output” folder

Required packages
-C++11 compiler or later
-C11 compiler or later
-MPI Version 3 Implementation
-POSIX Threads Implementation
-CMake Version 3 or later

To compile the application, run the compile.sh script
The build script assumes that you have gcc and g++ set to the C and C++ compiler of your choice. The build script will find your system’s MPI libraries if they are correctly installed and loaded. The application will be built in the FEA_Mesher2D directory where the compilation script is. If you want to change the directory where the application is built, you can set the location in src/CMakeLists.txt as the CMAKE_RUNTIME_OUTPUT_DIRECTORY

Licensing information available in LICENSE.txt
