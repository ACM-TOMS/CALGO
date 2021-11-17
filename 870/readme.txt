
------------------------------------------------------
THE MADD LIBRARY

MADD is a  method for decomposing 2D  geometries.  The method targets, in  addition to the
load balancing and the small size of the separators, the creation of good angles, close to
$\pi/2$,  formed by  between the  separators and  the external  boundary.  The  program is
written in c99 standard C++ and the library is wrapped into the madd class.

The program makes use of two  libraries; Triangle, a triangulation library by J. Shewchuk,
and Metis, a graph partitioning library by  G. Karypis, both written in C.  The copyrights
for Triangle belong to  J.  Shewchuk and for Metis to the  University of Minnesota.  These
libraries and are not part of our code, but  they are in the public domain and can be used
under the conditions  described in their source code and  their web-pages.  Please consult
the web-pages

     http://www.cs.cmu.edu/~quake/triangle.html  for Triangle
     http://www-users.cs.umn.edu/~karypis/metis/  for Metis

for more information.

The madd folder includes the following files:
  
  madd.cc      The source madd library.
  madd.h       The header of the library
  maddi.cc     A user interface that reads, decomposes and writes 
               geometries. It is also a good reference if you want to use 
               the madd library directly.
  L_Timer.*    Routines for timing. 
  
The data folder includes a set of test  geometries and the runtests script that runs a set
of decomposition tests.  The file format of the geometries is the same as the ones used by
Triangle (see the web page above).  The maddi program receives a set of simple commands to
read, decompose and write geometries, and in this folder you will find some simple sets of
these commands.

------------------------------------------------------
INSTALATION

The file install  contains a script that  creates the madd library.  Before  you start the
instalation  make sure  that the  file make.settings  defines the  correct  compilers. The
install script will download and compile the Triangle and the Metis library, so you should
have an active  internet connection before you  run the script.  Then it  will compile and
link the madd library and run a set of tests.


------------------------------------------------------ 
I will be glad to hear from you if you  have any comments or bugs to report.  If you use a
decomposition produced by madd in a publication, please include an acknowledgment as well. 

Leonidas Linardakis
College of William & Mary
lxlina@wm.edu
