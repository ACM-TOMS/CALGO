```
    _________   _____ ___    __
   / ____/   | / ___//   |  / /
  / __/ / /| | \__ \/ /| | / /
 / /___/ ___ |___/ / ___ |/ /___
/_____/_/  |_/____/_/  |_/_____/


```

Description
===========
This software is the open-source project EASAL: Efficient Atlasing, Analysis and Search of Molecular Assembly Landscapes. For a detailed explanation of the theory behind the concepts in this program, refer to the papers referenced in Publications.

This project leverages the convexity of the Cayley parameter space for atlasing and sampling molecular assembly configuration spaces. With respect to many performance measures the method is provably superior to the typical Monte Carlo approach. In particular it achieves comparable coverage of the assembly configuration space far more efficiently and permits flexible tradeoff between accuracy and efficiency. Hybridizing monte carlo with EASAL promises the best of both worlds.

This project should be viewed primarily as mathematical software. 

This README corresponds to the TOMS submission, namely, the backend of EASAL, without GUI and with text input and output. The source code for this can be found in the Source/backend_TOMS_Submission folder.
The experimental results in 4.1 of the paper, can be reproduced with this version using the sample input files given in the files directory. See section 5 of the included 'TOMSUserGuide.pdf' for detailed instructions on how to run the test driver.

An optional GUI has been included for intuitive visual verification of the
results.  This part of the code is not part of the TOMS submission. The source 
code for this can be found in the Source/optional_GUI folder.
Instructions on how to install, how to use and
major functionalities offered by the GUI have been detailed in the `Complete User
Guide' and the 'CompleteREADME.md'. 


Requirements
============
- Operating system: 
	- Ubuntu 12.04 or higher, 
	- Fedora 23 or higher, 
	- OSX Darwin or higher.
- C++ compiler
	- g++ version 4.8 or higher
	- clang++ version 3.3 or higher

Installation
============
- For the backend, all the required third party headers are included in the submission
	- Run ‘make‘ from the root/build directory.
	- To run EASAL run ‘bin/EASAL’ from the root/build directory.

Resources
=========
1. Input Files
	The 'files' directory in the root folder contains all the example input molecular data.
2. Data Directory 
	The data directory chosen by the user in the Input Window stores the atlas.
3. path_matrix.txt - The file in which the path matrix of a particular user input length for all pairs of 0D and 1D nodes is stored.
4. paths.txt - The file in which we output the shortest paths between 0D nodes that were clicked in the atlas view.

Usage
=====
- Command-line arguments
    - `-settings` Gives the settings file from which to load the input.

Example
=======
- Run EASAL from the command line by running the following command.
    `$bin/EASAL' - for the backend


Basic Code Overview
===================
Note that ^ means, explained further below.
```
.main
|
|__PointSet
|
|__^Atlas
|
|__^AtlasBuilder
|
|__.Settings
|
|__.SaveLoader
|
|__.readIn
```

```
.Atlas
|
|__.AtlasNode
   |
   |__.ActiveConstraintRegion
   |  |
   |  |__.CayleyPoint
   |     |
   |     |__.Orientation
   |
   |__.ActiveConstraintGraph
      |
      |__.ConvexChart
      |
      |__.CayleyParameterization
      |
      |__.CgMarker
```

```
.AtlasBuilder
|
|__.PointSet A/B
|
|__.PredefinedInteractions
|
|__.SaveLoader
|
|__.Atlas
|
|__.CgMarker
|
|__.ConstraintCheck
|
|__.CartesianRealizer
```

License
=========
EASAL is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
The GNU GPL license can be found at  <http://www.gnu.org/licenses/>.

Authors
=======
- James Pence
- Aysegul Ozkan
- Rahul Prabhu
- Troy Baker

Contact
=======
Meera Sitharam, CISE @ UF


Publications
============
 - Aysegul Ozkan and Meera Sitharam. 2011. ``EASAL: Efficient Atlasing, Analysis and Search of Molecular Assembly Landscapes''. In Proceedings of the ISCA 3rd International Conference on Bioinformatics and Computational Biology (BICoB-2011).
 - Aysegul Ozkan and Meera Sitharam. 2014. ``Best of Both Worlds: Uniform sampling in Cartesian and Cayley Molecular Assembly Configuration Space''. (2014). (on arxiv).
 - Aysegul Ozkan, Ruijin Wu, Jorg Peters, and Meera Sitharam. 2014b. ``Efficient Atlasing and Sampling of Assembly Free Energy Landscapes using EASAL: Stratification and Convexification via Customized Cayley Parametrization.'' (2014). (on arxiv).
