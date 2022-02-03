-----
About
-----
This folder contains the supplementary materials of the article:

  Agoston Roth, 2019. Algorithm 992: An OpenGL- and C++-based 
  function library for curve and surface modeling in a large 
  class of extended Chebyshev spaces, 
  ACM Transactions on Mathematical Software, 
  Vol. 45, No. 1, Article 13, https://doi.org/10.1145/3284979

All subfolders provide important readme-files.

In order to compile, build and install the proposed function 
library (•CAGD•), please follow the platform-dependent instructions
provided in the readme-file •01_Library/readme.txt•.

Your first task should be to install the free •OpenGL Extension
Wrangler Library• (GLEW) that is used by our library as a dependency 
in order to render any geometry by means of vertex buffer objects
and vertex/fragment shader programs.

Your second task should be to compile, build and install the 
proposed function library •CAGD•.

Finally, you can try to install the community edition of the
•Qt Software Development Kit• (Qt SDK) in order to compile,
build and run the test applications provided in the folder
•03_Examples•.

-------
Remarks
-------
We assume that the user has:
- a multi-core CPU; and
- a graphics adapter managed by a driver that supports at least 
  OpenGL 3.0.

For multi-threading we rely on a •C++ compiler• that supports at 
least •OpenMP 2.0• and which is compatible at least with the 
standard •C++ 11•. Please verify that your compiler is compatible
with •C++ 11•, otherwise the compilation of the source files
will fail.

Apart from •GLEW• no other external dependencies are used in
the proposed function library.

Although our library does not depend on the internal libraries 
of the •Qt SDK•, the five examples provided in the folder 
•03_Examples• rely on Qt-widgets in order to handle the graphical
user interfaces (GUIs) of the cross-platform test applications.

The main parts of the source codes of the 5 test applications 
(•03_Examples/Example_0x•, where x = 1,2,3,4,5) are also 
described in the •User Manual• without assuming Qt-based 
software development. The •User Manual• can be found in the 
folder •02_User_manual•. Therefore, if you do not want to 
install the •Qt SDK•, you can skip the testing step and, based
on the third chapter of the •User Manual•, you can write your 
own test applications.

For more details on the installation and testing processes,
consider the previously mentioned file •01_Library/readme.txt•.

