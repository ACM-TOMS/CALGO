SuiteSparse:  A Suite of Sparse matrix packages
(SuiteSparseQR subset)

------------------
SuiteSparse/README
------------------

================================================================================
QUICK START FOR MATLAB USERS:  uncompress the SuiteSparseQR.zip or
SuiteSparseQR.tar.gz archive file (they contain the same thing), then in the
MATLAB Command Window, cd to the SuiteSparseQR directory and type
SuiteSparse_install.  All packages will be compiled, and several demos will be
run.
================================================================================

May 10, 2011.  SuiteSparse VERSION 3.6.1
(SuiteSparseQR subset)

    AMD         approximate minimum degree ordering

    CAMD        constrained approximate minimum degree ordering

    COLAMD      column approximate minimum degree ordering

    CCOLAMD     constrained column approximate minimum degree ordering

    CHOLMOD     sparse Cholesky factorization.  Requires AMD, COLAMD, CCOLAMD,
                the BLAS, and LAPACK.  Optionally uses METIS.

    UFconfig    configuration file for all the above packages.  The
                UFconfig/UFconfig.mk is included in the Makefile's of all
                packages.  CSparse and RBio do not use UFconfig.

    SuiteSparseQR   sparse QR factorization

Some codes optionally use METIS 4.0.1
(http://www-users.cs.umn.edu/~karypis/metis).  To use METIS, place a copy of
the metis-4.0 directory in the same directory containing this README file.
Be sure that you do not have a nested metis-4.0/metis-4.0 directory; SuiteSparse
won't find METIS if you do this, which can happen with a zip file of metis-4.0
on Windows.  The use of METIS will improve the ordering quality.

By default, the UFconfig/UFconfig.mk file assumes that METIS is not installed.

Refer to each package for license, copyright, and author information.  All
codes are authored or co-authored by Timothy A. Davis, CISE Dept., Univ. of
Florida.  email: my last name @ cise dot ufl dot edu.

================================================================================
If you use SuiteSparse_install in MATLAB, stop reading here.
================================================================================



----------------------------
To use "make" in Unix/Linux:
----------------------------

(1) Use the right BLAS and LAPACK libraries

    Edit your UFconfig/UFconfig.mk file to point to the right compilers,
    and to the correct BLAS and LAPACK libraries.  There are many examples
    of different computer architectures there.  Scroll through to find yours,
    and uncomment those lines.

    By default, the UFconfig.mk file assumes -lblas and -llapack.

(2) Install Intel's Threading Building Blocks (TBB)

    This is optionally used by SuiteSparseQR.  Refer to the User Guide in 
    SuiteSparse/SPQR/Doc/spqr_user_guide.pdf for details.

(3) Configure METIS (or don't use METIS)

    cd to metis-4.0 and edit the Makefile.in file.  I recommend making these
    changes to metis-4.0/Makefile.in:

        CC = gcc
        OPTFLAGS = -O3

    And, if you want to use METIS in MATLAB and compile with "make" instead
    of using SuiteSparse_install.m:

        COPTIONS = -fexceptions -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

    Next, cd to metis-4.0 and type "make".

    If you do not wish to use METIS, then edit the UFconfig/UFconfig.mk file,
    and change the lines

        CHOLMOD_CONFIG =
        SPQR_CONFIG =

    to

        CHOLMOD_CONFIG = -DNPARTITION
        SPQR_CONFIG = -DNPARTITION

    Also change the line

        METIS = ../../metis-4.0/libmetis.a

    to

        METIS =

(4) Make other changes to UFconfig/UFconfig.mk as needed

    Edit the UFconfig/UFconfig.mk file as needed.  Directions are in that file.
    If you have compiled SuiteSparse already (partially or completely), then
    whenever you edit the UFconfig/UFconfig.mk file, you should then type
    "make purge" (or "make realclean") in this directory.

(5) Type "make" in this directory.  All packages will be be compiled.  METIS
    will be compiled if you have it.  Several demos will be run.

    To compile just the libraries, without running any demos, use
    "make library".

    The libraries will appear in */Lib/*.a.  Include files, as needed by user
    programs that use CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD, BTF, KLU, UMFPACK,
    LDL, etc. are in */Include/*.h.

    The METIS library is in metis-4.0/libmetis.a.  METIS Include files (not
    needed by the end user of SuiteSparse) are in located in metis-4.0/Lib/*.h.

(6) To install, type "sudo make install".  This will place copies of all
    libraries in /usr/local/lib, and all include files in /usr/local/include.
    You can change the install location by editting UFconfig.mk.  These
    directories must already exist.

(7) To uninstall, type "sudo make uninstall", which reverses "make install"
    by removing the SuiteSparse libraries from /usr/local/lib, and the
    include files from /usr/local/include.
