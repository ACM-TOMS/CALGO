KLU and BTF, versions 1.1.1, April 2010.  Timothy A. Davis and Ekanathan
Palamadai Natarajan.  Copyright (c) 2009, University of Florida.  All Rights
Reserved.  Licensed under the GNU LGPL license (see */Doc/lesser.txt).
For non-GNU license, please contact T. Davis at davis@cise.ufl.edu.

To compile and install KLU and BTF for MATLAB (version 7.3 or later), type
klupackage_install in the MATLAB command window while in this directory.
This works on Windows, the Mac, Linux, and any other operating system that
can run MATLAB.

To compile and install the C-callable KLU and BTF libraries, first edit the
UFconfig/UFconfig.mk file as needed.  Then type "make" in this directory in the
Unix/Linux/Cygwin command shell.

KLU requires two prior Collected Algorithms of the ACM: AMD (Algorithm 837) and
COLAMD (Algorithm 836).  These two Algorithms are included in this
distribution.  KLU requires BTF.  BTF does not require any other package.

KLU optionally makes use of CAMD, CCOLAMD, CHOLMOD, and METIS to provide
additional fill-reducing ordering options.  CAMD, CCOLAMD, and CHOLMOD are
contained in Algorithm 887 of the Collected Algorithms of the AMD, and are part
of the SuiteSparse meta-package, at http://www.cise.ufl.edu/research/sparse .
CHOLMOD requires LAPACK and the BLAS.  METIS can be obtained at
http://www-users.cs.umn.edu/~karypis/metis .  The exhaustive test suites for
KLU and BTF do require all of these optional packages.

For more information, see the KLU/BTF user guide provide with this distrubtion,
in KLU/Doc/KLU_UserGuide.pdf

