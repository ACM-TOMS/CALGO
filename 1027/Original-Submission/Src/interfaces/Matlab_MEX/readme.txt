The Matlab MEX interface allows to run NOMAD within the command line of Matlab.

Creating the Matlab MEX Interface to NOMAD requires to build source codes.
Building the interface requires compatibility of the versions of Matlab and
 the compiler that you want to use. Check the compatibility at:

https://www.mathworks.com/support/requirements/supported-compilers.html


************************
HOW TO BUILD AND INSTALL
************************
The interface build is managed by CMake that can be run *at NOMAD root*. 

The configuration command:
   cmake -DTEST_OPENMP=OFF -DBUILD_INTERFACE_MATLAB=ON -S . -B build/release.

Building the Matlab MEX interface is disabled when NOMAD uses OpenMP. 
Hence, the option -DTEST_OPENMP=OFF must be passed during configuration.

The command 
   cmake --build build/release 
or cmake --build build/release --config Release (for Windows) 
is used for building the selected configuration. 

The command
   cmake --install build/release 
must be run before using the Matlab nomadOpt function. 

Also, the Matlab command 
   addpath(strcat(getenv('NOMAD_HOME'),'/build/release/lib')) 
or addpath(strcat(getenv('NOMAD_HOME'),'/build/release/lib64')) 
must be executed to have access to the libraries and run the examples.

**********
HOW TO USE
**********
Some tests are proposed in the directory to check that everything 
is up and running.

All functionalities of NOMAD are available by using the nomadOpt 
function in Matlab command line. 

NOMAD parameters are provided in a Matlab structure with keywords and 
values using the same syntax as used in NOMAD parameter files. 
For example, 
   params = struct('initial_mesh_size','* 10','MAX_BB_EVAL','100');


