This SuiteSparse_SPQR directory is a subset of SuiteSparse,
version 4.5.5, for consideration as a Collected Algorithm
of the ACM.  It contains SPQR version 2.0.8, and all dependent
packages.

To compile, just type:

    make

in this directory to compile SPQR and its dependent packages.
The Makefile in this directory has been modified so that
it doesn't attempt to compile packages not in this subset.

Next, to compile and test the GPU acceleration for SPQR, do:

    cd SPQR/Demo
    make gpu

Also included is a copy of metis-5.1.0, whose license permits
its redistribution.  The use of METIS is optional, and the codes
will compile without it if it is not present.  If users do not
wish to use it, the folder can be deleted.

To install SPQR for MATLAB, just type

    SuiteSparse_install

in the MATLAB command window, while you are in this directory.  Note
that the MATLAB interface for SPQR does not yet support the GPU.

Packages in this subset:

    For ordering and analysis phases:
    AMD                         approx min degree
    CAMD                        constrained approx min degree
    CCOLAMD                     constrained column approx min degree
    COLAMD                      column approx min degree
    CHOLMOD                     sparse Cholesky
    metis-5.1.0                 METIS version 5.1.0

    Support packages:
    SuiteSparse_config          SuiteSparse-wide configuration
    UFget                       MATLAB interface for the SuiteSparse Matrix
                                    Collection (formerly the UF Sparse Matrix
                                    Collection)

    SPQR itself:
    SPQR                        multifrontal sparse QR
    SuiteSparse_GPURuntime      runtime support for SuiteSparse on the GPU
    GPUQREngine                 GPU kernels for SPQR

    Files in this top-level directory:
    ChangeLog                   SuiteSparse v4.5.4 ChangeLog
    Contents.m                  Contents of this subset, for MATLAB
    Makefile                    SuiteSparse Makefile, some packages removed
    README.txt                  SuiteSparse v4.5.4 README file, unmodified
    README_for_ACM_TOMS.txt     this file
    SuiteSparse_demo.m          MATLAB demo
    SuiteSparse_install.m       Compile SuiteSparse_SPQR for MATLAB

Packages in SuiteSparse but in this subset (not required by SPQR):

    BTF, CSparse, CXSparse, KLU, LDL, MATLAB_Tools, RBio, UMFPACK

Tim Davis, Nuri Yeralan, Sanjay Ranka, Wissam Sid-Lakhdar, Apr 17, 2016
davis@tamu.edu
http://faculty.cse.tamu.edu/davis

