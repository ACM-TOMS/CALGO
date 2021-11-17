Readme for the algorithm implementation of the paper

"Optimal Expression Evaluation Through the Use of Expression Templates
When Evaluating Dense Linear Algebra Operators"
----------------------------------------------------------

Algorithm implementation is in the directory "Algorithm
Implementation".  

All required ThirdParty code is located in the "ThirdParty" directory.

Most documentation is contained in the source files themselves, with
additional high-level information contained in this file.


Building the Testing Material
-----------------------

To support cross-platform compilation, we use CMake to generate the
appropriate build files.  The following steps are required to build
the testing material.


1.  Third Party dependencies

The testing material relies on the Nektar++ finite simulation
package.  While the expression template library does not itself
require Nektar++, many of the tests use the linear algebra library
that is a part of Nektar++.

Nektar++ is located in the Nektar directory.  Full instructions for
how to build Nektar++ can be found in 
compile_linux.html, compile_osx.html, or compile_windows.html.  

- It may not be necessary to buid boost on many Linux distributions,
  as they often come with boost already.  However, boost 1.48.0 is
  required, and many distributions distribute older versions.  

After building Nektar++, all relevant files will have been installed
to Nektar/build/dist directory.  Since Nektar++ uses shared libraries,
the path for your system must be adjusted to point to
Nektar/build/dist/bin.  

The Testing and example drivers are in the Src directory.  Go to
Src/Build and execute cmake, either from the command line using 

Src/Build > ccmake ../

or in the gui with the Build directory set to Src/Builds, and the
source directory set to Src.

You will likely see an error relating to Nektar++, indicating that
Nektar++_DIR can't be found.  Set it to Nektar/builds/dist/bin.

If boost can't be found, the varaible BOOST_ROOT can be specified to
indicate its location.

The testing library consists of compile-time tests and runtime-tests.
If the library compiles without errors, then it has passed all of the
compile time tests relating to parse tree manipulations and temporary
reduction.  To run the runtime tests, simply execute

Src/Build/dist> ExpressionTemplateTesting

The performance tests provided in the paper can be executed with 

Src/Build/dist> ExpressionTemplatePerformance <TestName> <ProblemSize> <ProblemInputs> <Iterations>

TestName is one of:

VectorAddition
MatrixMultiplication
SetExpression

Corresponding to the tests in the paper.

ProblemSize is the size of the matrix/vectors to be tests.
ProblemInputs is the number of operands in the expression.
Iterations is the number of times to evaluate the expression.




