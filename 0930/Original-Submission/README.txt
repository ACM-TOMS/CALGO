The Factorize package has an optional dependency on another Collected Algorithm
of the ACM (SuiteSparseQR), which itself depends on other Collected Algorithms
of the ACM (AMD, COLAMD, CHOLMOD, BLAS, LAPACK).  The Factorize directory
contains just the Factorize package itself.  To install it, just add the
Factorize directory to your MATLAB path.  This is the primary code being
submitted to ACM TOMS for this paper.  The sparse COD (complete orthogonal
decomposition) will not be available, since it requires SPQR.  The Factorize
package gracefully skips the sparse COD if SPQR is not available.  By default,
the sparse COD is used only for rank-deficient sparse matrices.

Factorize/Doc/factorize_article.pdf is a copy of the ACM TOMS journal
submission.

The SuiteSparse directory is the entire SuiteSparseQR package, including the
Factorize package (see SuiteSparse/MATLAB_Tools/Factorize).  This is included
here as a convenience for the reviewer, since normally an end-user who wants
both the Factorize package and SPQR would be asked to download all of
SuiteSparse. To install SuiteSparse in MATLAB, cd to the SuiteSparse directory
and run the SuiteSparseQR_install command at the MATLAB command line. 

The Factorize/Demo will work in both versions.  The full test in
Factorize/Test/test_all.m requires SPQR, however, so it will only work in the
SuiteSparse/MATLAB_Tools/Factorize version (and only if SuiteSparse_install is
used first).

Tim Davis
Dec 13, 2012
