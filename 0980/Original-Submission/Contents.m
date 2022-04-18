% Welcome to SuiteSparse : a Suite of Sparse matrix packages, containing a
% collection of sparse matrix packages authored or co-authored by Tim Davis.
% Only the primary MATLAB functions are listed below.
%
% This software bundle includes only SuiteSparseQR (SPQR) and its dependent
% packages.
%
% Example:
%   SuiteSparse_install
% compiles and installs all of SuiteSparse, and runs several demos and tests.
%
%-------------------------------------------------------------------------------
% Ordering methods:
%-------------------------------------------------------------------------------
%
%   amd2         - approximate minimum degree ordering.
%   colamd2      - column approximate minimum degree ordering.
%   symamd2      - symmetrix approximate min degree ordering based on colamd.
%   camd         - constrained amd.
%   ccolamd      - constrained colamd.
%   csymamd      - constrained symamd.
%
%-------------------------------------------------------------------------------
% CHOLMOD: a sparse supernodal Cholesky update/downdate package:
%-------------------------------------------------------------------------------
%
%   cholmod2     - computes x=A\b when A is symmetric and positive definite.
%   chol2        - same as MATLAB chol(sparse(A)), just faster.
%   lchol        - computes an LL' factorization.
%   ldlchol      - computes an LDL' factorization.
%   ldlupdate    - updates an LDL' factorization.
%   resymbol     - recomputes symbolic LL or LDL' factorization.
%   ldlsolve     - solves Ax=b using an LDL' factorization.
%   ldlsplit     - splits LD into L and D.
%   metis        - interface to METIS node-nested-dissection.
%   nesdis       - interface to CHOLMOD's nested-dissection (based on METIS).
%   septree      - prune a separator tree.
%   bisect       - interface to METIS' node bisector.
%   analyze      - order and analyze using CHOLMOD.
%   etree2       - same as MATLAB "etree", just faster and more reliable.
%   sparse2      - same as MATLAB "sparse", just faster.
%   symbfact2    - same as MATLAB "symbfact", just faster and more reliable.
%   sdmult       - same as MATLAB S*F or S'*F (S sparse, F full), just faster.
%   ldl_normest  - compute error in LDL' factorization.
%   lu_normest   - compute error in LU factorization.
%   mread        - read a sparse matrix in Matrix Market format
%   mwrite       - write a sparse matrix in Matrix Market format
%   spsym        - determine the symmetry of a sparse matrix
%
%-------------------------------------------------------------------------------
% SuiteSparseQR: multifrontal rank-revealing sparse QR
%-------------------------------------------------------------------------------
%
%   spqr            - sparse QR
%   spqr_solve      - x=A\b using SuiteSparseQR
%   spqr_qmult      - y=Q*x, Q'*x, x*Q, or x*Q' using Q in Householder form
%
%-------------------------------------------------------------------------------
% Other packages:
%-------------------------------------------------------------------------------
%
%   UFget           MATLAB interface to the SuiteSparse Matrix Collection
%
%-------------------------------------------------------------------------------
%
% For help on compiling SuiteSparse or the demos, testing functions, etc.,
% please see the help for each individual package.
%
% Copyright 2014, Timothy A. Davis, http://www.suitesparse.com.

help SuiteSparse
