--------
Contents
--------
0. System requirements
1. Compilers and external dependencies (GLEW)
	1.a. Installing GLEW under Unix/Linux
	1.b. Installing GLEW under Macintosh
	1.c. Installing GLEW under Windows
2. Installation instructions
	2.1. Installation by using command prompts and make-files
		2.1.a. Under Unix/Linux
		2.1.b. Under Macintosh
		2.1.c. Under Windows
	2.2. Installation by using the community edition of the Qt Creator
		2.2.1. Installing the community edition of the Qt
		       Software Development Kit
			2.2.1.a. Under Unix/Linux
			2.2.1.b. Under Macintosh
			2.2.1.c. Under Windows
3. About the examples/test applications
	3.1. Short description of the examples

----------------------
0. System requirements
----------------------
We assume that the user has:
- a multi-core CPU; and
- a graphics adapter managed by a driver that supports at least 
  OpenGL 3.0.

---------------------------------------------
1. Compilers and external dependencies (GLEW)
---------------------------------------------
In order to render the geometry, we use vertex buffer and shader
objects through the •OpenGL Extension Wrangler Library• (GLEW), 
and for multi-threading we rely on a C++ compiler that supports 
at least •OpenMP 2.0• and which is compatible at least with the
standard •C++ 11•. Please verify that your compiler is compatible
with •C++ 11•, otherwise the compilation of the source files
will fail.

Apart from •GLEW• no other external dependencies are used.

The source (for all platforms) and the zipped Windows-binaries of 
the library •GLEW• can be downloaded from the link 
http://glew.sourceforge.net/

-------------------------------------
1.a. Installing GLEW under Unix/Linux
-------------------------------------
In case of Ubuntu, one can install the library •GLEW• either by
following the instructions of the link 

http://glew.sourceforge.net/ 

or simply by calling the next commands in a terminal:

  sudo apt update
  sudo apt install libglew-dev
  
In case of success, the release build of the library •GLEW• 
and its header file will be installed to the folders 

  /usr/lib 
  
and 

  /usr/include/GL, 
  
respectively.

------------------------------------
1.b. Installing GLEW under Macintosh
------------------------------------
On •Mac Os X 10.10+• one can install the library •GLEW• by running 
(as a super user) the following commands in a terminal:

  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
  brew install glew

The first of these two commands installs the package manager •Homebrew•,
while the second one installs the library •GLEW• by means of •Homebrew•.
If your system already has •Homebrew•, skip the first command.
  
In case of success, the release build of the library •GLEW• 
and its header file will be installed to the folders and 

  /usr/local/Cellar/glew/x.y.z/lib/

and

  /usr/local/Cellar/glew/x.y.z/include/GL,

respectively, where x.y.z denotes the version number of the 
installed library •GLEW•.

These paths will be important in the case of the Qt-based project 
file (*.pro) file of each test application provided in the 
directory •../03_Examples•. After opening such a project file in 
Qt Creator, one should edit within the Macintosh platform related 
scope mac{...} the version number x.y.z which appears in the lines

  INCLUDEPATH += "/usr/local/Cellar/glew/x.y.z/include/"
  LIBS        += -L"/usr/local/Cellar/glew/x.y.z/lib/" -lGLEW

of the just loaded .pro file. (When this readme-file was created, 
the version number of the library •GLEW• was 2.1.0, but this may 
change over the years.)

----------------------------------
1.c. Installing GLEW under Windows
----------------------------------
The link 

http://glew.sourceforge.net/ 

also provides the pre-compiled 32- and 64-bit release binaries 
glew32.{lib|dll} together with their common include file (glew.h)
in a zipped file that is named similarly to 

glew-x.y.z-win32.zip, 

where x.y.z denotes the current version number of the library 
•GLEW•.

We assume that the user has downloaded, extracted and copied 
the previously mentioned files as follows:

glew-x.y.z/bin/Release/Win32/glew32.dll --> ../00_Dependencies/Lib/GL/x86/glew32.dll
glew-x.y.z/bin/Release/x64/glew32.dll   --> ../00_Dependencies/Lib/GL/x64/glew32.dll

glew-x.y.z/lib/Release/Win32/glew32.lib --> ../00_Dependencies/Lib/GL/x86/glew32.lib
glew-x.y.z/lib/Release/x64/glew32.lib   --> ../00_Dependencies/Lib/GL/x64/glew32.lib

glew-x.y.z/include/GL/glew.h            --> ../00_Dependencies/Include/GL/glew.h

----------------------------
2. Installation instructions
----------------------------
In order to compile, build and install the proposed function 
library (CAGD), we provide two possibilities. The first of these 
is based on the usage of command prompts and traditional make-files,
while the other one relies on a (configurable) Qt-based project.

---------------------------------------------------------
2.1. Installation by using command prompts and make-files
---------------------------------------------------------

-----------------------
2.1.a. Under Unix/Linux
-----------------------
One should:
- open a terminal and navigate into the folder •01_Library•;
- call the commands:

  make -f Makefile.Linux.Debug
  make -f Makefile.Linux.Release

  sudo make -f Makefile.Linux.Debug install
  sudo make -f Makefile.Linux.Release install

In case of success, these commands will install the obtained
debug and release builds of the function library and the contents 
of their common include folder CAGD to the directories

  /usr/local/bin 

and 

  /usr/local/include/CAGD,

respectively. (Both in debug and release modes, all Qt-based test 
applications assume the latter library and include paths.)

The 32- and 64-bit debug builds will be named as libcagd32d.a and
libcagd64d.a, respectively.

The 32- and 64-bit release builds will be named as libcagd32.a 
and libcagd64.a, respectively.

If one modifies the variables PREFIX and INSTALL_DESTINATION_DIR
of the makefiles specified above, one can also install these 
libraries and include files into different folders, but in this 
case one has to manually modify the variables INCLUDEPATH and LIBS
of the .pro files of all examples as well. (The test applications 
can be found in the folder •../03_Examples•.)

If the debug/release libraries and their common include files are 
not required anymore, one can also uninstall them, by invoking the
commands:

  sudo make -f Makefile.Linux.Debug uninstall
  sudo make -f Makefile.Linux.Release uninstall

----------------------
2.1.b. Under Macintosh
----------------------
It may happen that you need to update the compiler •gcc•. This 
can be done by means of the package manager •Homebrew•. 

If you do not already have •Homebrew•, you can install it by 
calling (as a super user) the command

  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

in a terminal. After this, you can  update your compiler by 
running the command

  brew install gcc
  
as a super user.  

In order to avoid conflicts with already existing compilers, 
•Homebrew• usually installs the updated compiler to the path

  /usr/local/opt/gcc/bin/gcc-x
  
where x denotes the newly downloaded version number of •gcc• 
(e.g., when this readme-file was created, the number 9 appeared
instead of the letter x). If you have updated the compiler •gcc•,
then you will also need to edit the value of the variable •CXX• 
in both make-files •Makefile.Mac.Debug• and •Makefile.Mac.Release• 
appearing in the current folder. By default, the variable •CXX• 
is set to •gcc• (i.e., to your old compiler) in these files. You 
will have to change the line

  CXX = gcc
  
to

  CXX = /usr/local/opt/gcc/bin/gcc-x
  
in both of these files, where the letter x depends on the version 
number of the new compiler •gcc•.

After this, one should:
- open a terminal and navigate into the folder •01_Library•;
- call the commands:

  make -f Makefile.Mac.Debug
  make -f Makefile.Mac.Release

  sudo make -f Makefile.Mac.Debug install
  sudo make -f Makefile.Mac.Release install

In case of success, these commands will install the obtained
debug and release builds of the function library and the contents
of their common include folder CAGD to the directories

  /usr/local/bin 

and 

  /usr/local/include/CAGD,

respectively. (Both in debug and release modes, all Qt-based test
applications assume the latter library and include paths.)

The 32- and 64-bit debug builds will be named as libcagd32d.a and
libcagd64d.a, respectively.

The 32- and 64-bit release builds will be named as libcagd32.a 
and libcagd64.a, respectively.

If one modifies the variables PREFIX and INSTALL_DESTINATION_DIR 
of the makefiles specified above, one can also install these 
libraries and include files into different folders, but in this 
case one has to manually modify the variables INCLUDEPATH and LIBS 
of the .pro files of all examples as well. (The test applications 
can be found in the folder •../03_Examples•.)

If the debug/release libraries and their common include files are 
not required anymore, one can also uninstall them, by invoking 
the commands:

  sudo make -f Makefile.Mac.Debug uninstall
  sudo make -f Makefile.Mac.Release uninstall

--------------------
2.1.c. Under windows
--------------------
We assume that the user:
- has already installed the community edition of the •Microsoft 
  Visual Studio C++• (release year 2015 or above);
- has launched both the •x86• and •x64 Native Tools Command Prompts 
  for Visual Studio• and has navigated in both command prompts into
  the folder •01_Library• (if you do not know where these command
  prompts are on your system, the links
  
  https://docs.microsoft.com/en-us/cpp/build/building-on-the-command-line?view=vs-2019
  https://docs.microsoft.com/en-us/dotnet/framework/tools/developer-command-prompt-for-vs
  
  may be helpful).

In order to obtain both the debug and release builds of the 32-bit
variant of the proposed function library, one should call the 
commands

  nmake -f Makefile.Windows.Debug.x86
  nmake -f Makefile.Windows.Release.x86

in the •x86 Native Tools Command Prompt for Visual Studio•. 

In order to obtain both the debug and release builds of the 64-bit
variant of the proposed function library, one should call the
commands

  nmake -f Makefile.Windows.Debug.x64
  nmake -f Makefile.Windows.Release.x64

in the •x64 Native Tools Command Prompt for Visual Studio•.

In case of success:
- the 32- and 64-bit debug builds will be named as cagd32d.lib and 
  cagd64d.lib, respectively;
- the 32- and 64-bit release builds will be named as cagd32.lib and
  cagd64.lib, respectively.

These files will be copied to the folder •../00_Dependencies/Lib/CAGD•.

The include folder CAGD will be copied into the directory 
•../00_Dependencies/Include/•.

On Windows platforms, all examples in •../03_Examples• rely on the 
previously mentioned two directories.

------------------------------------------------------------------
2.2. Installation by using the community edition of the Qt Creator
------------------------------------------------------------------
The folder •cagd• contains the Qt-based project file •cagd.pro• 
that can be used to compile and build both the 32- and 64-bit 
variants of the proposed static function library.

The project file •cagd/cagd.pro• does not use Qt-specific (internal)
libraries, and the •Qt Creator• is only used as an integrated
development environment (IDE) that can be installed on •Unix/Linux•,
•Macintosh•, and •Windows•.

However, the test applications included in the folder 
•../03_Examples• are based on Qt-specific internal libraries 
in order to create a cross-platform graphical user interface (GUI).
Therefore, in order to build and run these test applications, 
one should install the free community edition of the •Qt Software
Development Kit• (Qt SDK), even if our library was built by 
using make-files.

----------------------------------------------------------
2.2.1. Installing the community edition of the Qt Software 
       Development Kit
----------------------------------------------------------
The latest community edition of the •Qt SDK• can be installed as
follows:

-------------------------
2.2.1.a. Under Unix/Linux
-------------------------
Either use the link

http://download.qt.io/official_releases/online_installers/qt-unified-linux-x64-online.run

or call the commands

  sudo apt update
  sudo apt install qtcreator qt5-default

where the parameter qt5-default should be changed based on the 
latest Qt version number.

------------------------
2.2.1.b. Under Macintosh
------------------------
Use the link

http://download.qt.io/official_releases/online_installers/qt-unified-mac-x64-online.dmg

----------------------
2.2.1.c. Under Windows
----------------------
We assume that the user has already installed that variant of 
the community edition of the •Microsoft Visual Studio C++• which 
is compatible with the newest variant of the community edition of 
the •Qt SDK•. (It may happen that the latest •Qt SDK• was compiled 
with a prior version of the •Microsoft Visual Studio C++•.)

By using one of the links

https://www.visualstudio.com/thank-you-downloading-visual-studio/?sku=Community&rel=15
https://visualstudio.microsoft.com/vs/older-downloads/

at first install the corresponding variant of the •Microsoft Visual
Studio C++•, then install the latest community edition of the 
•Qt SDK•, by accessing the link

http://download.qt.io/official_releases/online_installers/qt-unified-windows-x86-online.exe

---------------------------------------
3. About the examples/test applications
---------------------------------------
Although the proposed library does not depend on Qt-specific 
internal libraries, the five examples provided in the folder 
•../03_Examples• rely on Qt-widgets in order to handle the GUIs
of the cross-platform test applications.

The main parts of the source codes of the 5 test applications 
(•../03_Examples/Example_0x•, where x = 1,2,3,4,5) are also 
described in the •User Manual• without assuming Qt-based software
development. The •User Manual• can be found in the folder 
•../02_User_manual•.

--------------------------------------
3.1. Short description of the examples
--------------------------------------

------------
•Example_01•
------------
Illustrates the creation of different types of EC space objects,
differentiates and renders the unique normalized B-bases of the
constructed spaces.

------------
•Example_02•
------------
Generates different types of B-curves, and performs order elevation
division on them.

------------
•Example_03•
------------
Illustrates the B-representation (i.e., the control point
exact description) of several arcs of a logarithmic spiral given 
in traditional parametric form. Note that this exponential-
trigonometric curve can be approximated but cannot be exactly 
described by means of the standard (rational) Bézier and 
non-uniform B-spline curve modeling tools.

------------
•Example_04•
------------
Randomly generates a B-surface by using different types of 
normalized B-bases in directions u and v, then performs order 
elevation and subdivision on it.

------------
•Example_05•
------------
Illustrates the B-representation (i.e., the control point based
exact description) of several patches of an exponential-trigonometric
transcendental surface given in traditional parametric form that 
can be approximated but cannot be exactly represented by the 
standard (rational) Bézier and non-uniform B-spline surface 
modeling tools.


