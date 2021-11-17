-------------------------------------------------------------------
readme.txt for the Fortran 90 implementation of MINRES-QLP
(real and complex versions)

09 Sep 2013: Version 27
-------------------------------------------------------------------

The software for MINRES-QLP is provided by SOL, Stanford University
under the terms of the OSI Common Public License (CPL)
   http://www.opensource.org/licenses/cpl1.0.php
or the BSD License
   http://www.opensource.org/licenses/bsd-license.php

Please send comments to
   Sou-Cheng Choi <sctchoi@uchicago.edu>
   Computation Institute
   University of Chicago
   Chicago, IL 60637, USA

   Michael Saunders <saunders@stanford.edu>
   Systems Optimization Laboratory (SOL)
   Dept of Management Science and Engineering
   Stanford University
   Stanford, CA 94305-4121, USA

-------------------------------------------------------------------
SOURCE CODE
-------------------------------------------------------------------
The source code for MINRESQLP is the following files:

Real:                               Complex:
   minresqlpDataModule.f90             minresqlpDataModule.f90
   minresqlpBlasModule.f90             minresqlpBlasModule.f90
   minresqlpModule.f90                zminresqlpDataModule.f90
   mm_ioModule.f90                     mm_ioModule.f90
   minresqlpReadMtxModule.f90          minresqlpReadMtxModule.f90
   minresqlpTestModule.f90            zminresqlpTestModule.f90  
   minresqlpTestProgram.f90           zminresqlpTestProgram.f90
					 
-------------------------------------------------------------------
SYMMETRIC PROBLEMS
-------------------------------------------------------------------
To compile the real code and run the test program on Linux or Unix,
proceed as follows:

   make
   ./minresqlptest
   grep appears MINRESQLP.txt | cat -n

   "symortho  appears to be successful" should occur 10 times.
   "minresqlp appears to be successful" should occur 78 times.

Some of the tests use data in Matrix Market (MM) format
downloaded from Tim Davis's UFL Sparse Matrix Collection
   http://www.cise.ufl.edu/research/sparse/matrices/
and Leslie Foster's SJSU Singular Sparse Matrix Collection
   http://www.math.sjsu.edu/singular/matrices/

Three formats are tested using examples in the following directories:
   DataMtx/CDS/ (coordinate double  symmetric)
   DataMtx/CPS/ (coordinate pattern symmetric)
   DataMtx/CRS/ (coordinate real    symmetric)
Each matrix is used twice with variable "consis" true or false.
If consis=true, the rhs is computed as b = A*xtrue
and the computed solution should have small norm(r).
If consis=false, the rhs is altered in b(n-1:n).
If A is singular, norm(r) will not be small, but norm(A*r) will be.

-------------------------------------------------------------------
HERMITIAN PROBLEMS
-------------------------------------------------------------------
To compile the complex code and run the test program on Linux or Unix,
proceed as follows:

   make -f zMakefile
   ./zminresqlptest
   grep appears zMINRESQLP.txt | cat -n

   "zsymortho  appears to be successful" should occur 20 times.
   "zminresqlp appears to be successful" should occur 49 times.

Some Matrix Market examples are used from the following directory:
   DataMtx/CCH/ (coordinate complex Hermitian)
Again, each matrix is used twice with different rhs's.

-------------------------------------------------------------------
REFERENCES
-------------------------------------------------------------------
MINRESQLP is an implementation of the algorithm described in
the following references:

Sou-Cheng Choi,
Iterative Methods for Singular Linear Equations and Least-Squares
Problems, PhD thesis, ICME, Stanford University, 2006.

Sou-Cheng Choi, Christopher Paige, and Michael Saunders,
MINRES-QLP: A Krylov subspace method for indefinite or singular
symmetric systems, SIAM Journal of Scientific Computing 33:4,
1810-1836, 2011.

Sou-Cheng Choi and Michael Saunders,
ALGORITHM & DOCUMENTATION: MINRES-QLP for singular symmetric and
Hermitian linear equations and least-squares problems,
Technical Report, ANL/MCS-P3027-0812, Computation Institute,
University of Chicago/Argonne National Laboratory, 2012.

Sou-Cheng Choi and Michael Saunders,
ALGORITHM xxx: MINRES-QLP for singular Symmetric and Hermitian
linear equations and least-squares problems,
ACM Transactions on Mathematical Software, accepted Sep 2013.

The above documents and FORTRAN 90 and MATLAB implementations are
downloadable from
   http://www.stanford.edu/group/SOL/software.html
or http://home.uchicago.edu/sctchoi/
