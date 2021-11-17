Thank you for downloading the 
MATLAB Tensor Classes for Fast Algorithm Prototyping. 

There files are distributed under the following agreement:

Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
Sandia National Laboratories, 2004.  Please address questions or
comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
distribute, display, and use this work, under the following
conditions. (1) You must give the original authors credit. (2) You may
not use or redistribute this work for commercial purposes. (3) You may
not alter, transform, or build upon this work. (4) For any reuse or
distribution, you must make clear to others the license terms of this
work. (5) Any of these conditions can be waived if you get permission
from the authors.

Documentation on using these classes can be found in:

Brett W. Bader and Tamara G. Kolda, MATLAB Tensor Classes for Fast
Algorithm Prototyping, SAND2004-5187, Sandia National Labs, Livermore,
CA, October 2004.

Please cite our ACM TOMS publication (should be published 2006) in any
resulting publications.

Installation instructions:

1. Unpack the tar/zip file into the directory of your choice.  Note
   the directory for addition to your MATLAB path in step 3.

2. Start MATLAB.

3. Add the Tensor directory to your path with the "addpath" or "path"
   commands.  For example:

   addpath ~/Tensor
      or
   path('~/Tensor', path);

   The directory provided must have access to the class directories
   (@tensor, @cp_tensor, @tucker_tensor, and @tensor_as_matrix) for
   MATLAB to recognize them.  The command "addpath" prepends the
   specified directory to the current matlabpath whereas "path" has
   more flexibility.


