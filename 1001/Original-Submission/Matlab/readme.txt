
IPscatt---a MATLAB Toolbox for the Inverse Medium Problem in Scattering

-----

# Overview

IPscatt is a free, open-source toolbox in MATLAB that solves the inverse medium problem in time-independent scattering (also known as time-harmonic scattering) in two- and three-dimensional setting.

Requirements:
Apart from MATLAB the software library IPscatt requires the "Signal Processing Toolbox" and the "Image Processing Toolbox". It was extensively tested in MATLAB R2016b under Linux, but short tests under Windows 7 and 10 were also successful.

Overview:
The folder "code" contains the source code, the folder "doc" the source code documentation and the folder "guide" the user guide.

Update:
The most recent version of IPscatt is provided on http://www.fbuergel.de/ipscatt/

-----

# Demonstration and Tests

Make sure to be in the folder code if you run one of the following routines, that are described in detail below:

demo              % Demonstrates the main program (run time: less than 1 min).

start;            % Representative demonstration of IPscatt (run time: 3 min)

seti = runtests() % Tests of the internal functions (recommended after code modification) (run time: 45 min).

For comparison purposes the expected outputs of these routines are available in the folder outputExpected in the subfolders demo, start and tests. In addition, the output from the MATLAB Command Window is available in commandWindow.txt in each case. These computations were carried out on a workstation with an Intel(R) Core(TM) i7-3770 CPU with 3.40 GHz and 32 GByte RAM.

For practical guides for using the library with code snippets as well as resulting figures we refer to the user guide of IPscatt.


# Demonstration

demo              % Demonstrates the main program (run time: less than 1 min).

For a first impression of IPscatt we use the routine demo in MATLAB with code as current folder. This executes a demonstration script in which all important parts of IPscatt are executed. However, to keep the run time low the number of transmitters/receivers, discretization points and reconstruction steps was set very low. This means that the resulting reconstruction is not representative for IPscatt. The generated figures are saved in a folder with the current date as prefix inside the directory output.


# Start: Representative demonstration of IPscatt

start;            % Representative demonstration of IPscatt (run time: 3 min)

In comparison to demo, the routine start results in a reconstruction that is representative for IPscatt. The generated figures are saved in a folder with the current date as prefix inside the directory output.


# Test of Internal Functions

The toolbox IPscatt was designed for further development. Therefore IPscatt supports the user by providing tests of internal functions. In particular, it tests the proper working of the forward operator and its derivative. These tests may take a long time and are started by 

seti = runtests() % Tests of the internal functions (recommended after code modification) (run time: 45 min).

Again, the current folder should be code. Note that some (but not all) generated figures are saved in a folder with the current date as prefix inside the directory output. For more information see the source code documentation of the routine runtests.

-----

# Authors

The toolbox IPscatt was developed and is maintained by

Florian Buergel <fbuergel@uni-bremen.de>
  University of Bremen
  Center for Industrial Mathematics
  Bibliothekstr. 5
  28359 Bremen
  GERMANY

Dr. Kamil S. Kazimierski <kazimier@uni-graz.at>
  Institute of Mathematics and Scientific Computing
  Karl-Franzens-University Graz
  Heinrichstr. 36/III
  8010 Graz
  AUSTRIA

(The Institute of Mathematics and Scientific Computing is a member of NAWI Graz (www.nawigraz.at) and BioTechMed Graz (www.biotechmed.at).)

Prof. Dr. Armin Lechleiter
  University of Bremen
  Center for Industrial Mathematics
  Bibliothekstr. 5
  28359 Bremen
  GERMANY

-----

# Feedback

If you have comments, questions or suggestions regarding IPscatt, do not hesitate to contact Florian Buergel <fbuergel@uni-bremen.de>. Alternatively, you might contact Kamil S. Kazimierski <kazimier@uni-graz.at>.

If IPscatt is useful for you, please report us. We are really interested in what application you are using IPscatt for.

-----

# Legal Information & Credits


# Third Party Code:

This toolbox contains third-party code in the folder 3rdparty.
The corresponding licenses are provided in this folder.


# License of IPscatt:

All remaining files contained in this package and its sub-directories are published with the implicit license:

IPscatt---a MATLAB Toolbox for the Inverse Medium Problem in Scattering
Copyright (C) 2017 Florian Buergel, Kamil S. Kazimierski, and Armin Lechleiter

This software was written by Florian Buergel, Kamil S. Kazimierski and Armin Lechleiter. It was developed at the Center for Industrial Mathematics, University of Bremen, Germany, and the Institute of Mathematics and Scientific Computing, University of Graz, Austria.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

-----

# Directory Structure

As mentioned the package contains the folders "code" with the source code, the folder "doc" with the source code documentation and the folder "guide" with this user guide. The directory structure of "code" is described in the following list.

* 3rdparty		(third-party source code)
* conv			(convenience functions)
* docCreate		(creates the documentation)
* guides		(code for the guides)
* incontrasts	(input of predefined contrasts)
* inseti		(input of structure array seti)
* output		(place to store files and figures)
* proc			(process)
  - auxi		(auxiliary routines)
  - expData		(working with real-world data)
  - expSetup	(experimental set-up)
  - intOps		(integral operators)
  - norms		(definitions of norms)
  - operators	(operators)
  - plots		(functions for plots)
  - plotsAux	(auxilary functions for plots)
  - recon		(reconstruction process)
  - reconAux	(associated auxiliary functions)
  - setData		(sets geometry and simulation)
  - setInput	(general input, e.g. make folders)
* tests			(test functions)
* auxi			(auxiliary functions)

-----

