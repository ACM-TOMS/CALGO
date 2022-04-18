/* ========================================================================== */
/* === symamd mexFunction =================================================== */
/* ========================================================================== */

/*
    Usage:

	P = symamd (A) ;

	P = symamd (A, knobs) ;

	[ P, stats ] = symamd (A) ;

	[ P, stats ] = symamd (A, knobs) ;

    Returns a permutation vector P such that the Cholesky factorization of
    A (P,P) tends to be sparser than that of A.  This routine provides the same
    functionality as SYMMMD, but tends to be much faster and tends to return a
    better permutation vector.  Note that the SYMMMD m-file in
    MATLAB 5.2 also performs a symmetric elimination tree post-ordering.  This
    mexFunction does not do this post-ordering.  This mexFunction is a
    replacement for the p = sparsfun ('symmmd', A) operation.

    A must be square, and is assummed to have a symmetric nonzero pattern.
    Only the nonzero pattern of the lower triangular portion of A is accessed.
    This routine constructs a matrix M such that the nonzero pattern of M'M is
    equal to A (assuming A has symmetric pattern), and then performs a column
    ordering of M using colamd.

    The knobs and stats vectors are optional:

	knobs (1)	rows and columns with more than (knobs(1))*n entries
			are removed prior to ordering, and placed last in
			the output ordering.  If knobs is not present, then the
			default of 0.5 is used.

	knobs (2)	print level, similar to spparms ('spumoni')

	stats (1)	the number of dense (or empty) rows and columms.  These
			are ordered last, in their natural order.

	stats (2)	(same as stats (1))

	stats (3)	the number of garbage collections performed.

	stats (4)	return status:

			0:  matrix is a valid MATLAB matrix.

			1:  matrix has duplicate entries or unsorted columns.
			    This should not occur in a valid MATLAB matrix,
			    but the ordering proceeded anyway by ignoring the
			    duplicate row indices in each column.  See
			    stats (5:7) for more information.

	stats (5)	highest numbered column that is unsorted or has
			duplicate entries (zero if none)

	stats (6)	last seen duplicate or unsorted row index
			(zero if none)

	stats (7)	number of duplicate or unsorted row indices

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (davis@cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Date:

	September 8, 2003.  Version 2.3.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Notice:

	Copyright (c) 1998-2003 by the University of Florida.
	All Rights Reserved.

	See http://www.cise.ufl.edu/research/sparse/colamd (the colamd.c
	file) for the License.

    Availability:

	The colamd/symamd library is available at

	    http://www.cise.ufl.edu/research/sparse/colamd/

	This is the http://www.cise.ufl.edu/research/sparse/colamd/symamdmex.c
	file.  It requires the colamd.c and colamd.h files.

*/

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include "colamd.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

/* ========================================================================== */
/* === symamd mexFunction =================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
    /* === Local variables ================================================== */

    int *perm ;			/* column ordering of M and ordering of A */
    int *A ;			/* row indices of input matrix A */
    int *p ;			/* column pointers of input matrix A */
    int n_col ;			/* number of columns of A */
    int n_row ;			/* number of rows of A */
    int full ;			/* TRUE if input matrix full, FALSE if sparse */
    double knobs [COLAMD_KNOBS] ; /* colamd user-controllable parameters */
    double *out_perm ;		/* output permutation vector */
    double *out_stats ;		/* output stats vector */
    double *in_knobs ;		/* input knobs vector */
    int i ;			/* loop counter */
    mxArray *Ainput ;		/* input matrix handle */
    int spumoni ;		/* verbosity variable */
    int stats [COLAMD_STATS] ;	/* stats for symamd */

    /* === Check inputs ===================================================== */

    if (nrhs < 1 || nrhs > 2 || nlhs < 0 || nlhs > 2)
    {
	mexErrMsgTxt (
	"symamd: incorrect number of input and/or output arguments.") ;
    }

    /* === Get knobs ======================================================== */

    colamd_set_defaults (knobs) ;
    spumoni = 0 ;

    /* check for user-passed knobs */
    if (nrhs == 2)
    {
	in_knobs = mxGetPr (prhs [1]) ;
	i = mxGetNumberOfElements (prhs [1]) ;
	if (i > 0) knobs [COLAMD_DENSE_ROW] = in_knobs [COLAMD_DENSE_ROW] ;
	if (i > 1) spumoni = (int) in_knobs [1] ;
    }

    /* print knob settings if spumoni > 0 */
    if (spumoni > 0)
    {
	mexPrintf ("symamd: dense row/col fraction: %f\n",
	    knobs [COLAMD_DENSE_ROW]) ;
    }

    /* === If A is full, convert to a sparse matrix ========================= */

    Ainput = (mxArray *) prhs [0] ;
    if (mxGetNumberOfDimensions (Ainput) != 2)
    {
	mexErrMsgTxt ("symamd: input matrix must be 2-dimensional.") ;
    }
    full = !mxIsSparse (Ainput) ;
    if (full)
    {
	mexCallMATLAB (1, &Ainput, 1, (mxArray **) prhs, "sparse") ;
    }

    /* === Allocate workspace for symamd ==================================== */

    /* get size of matrix */
    n_row = mxGetM (Ainput) ;
    n_col = mxGetN (Ainput) ;
    if (n_col != n_row)
    {
	mexErrMsgTxt ("symamd: matrix must be square.") ;
    }

    A = mxGetIr (Ainput) ;
    p = mxGetJc (Ainput) ;
    perm = (int *) mxCalloc (n_col+1, sizeof (int)) ;

    /* === Order the rows and columns of A (does not destroy A) ============= */

    if (!symamd (n_col, A, p, perm, knobs, stats, &mxCalloc, &mxFree))
    {
	symamd_report (stats) ;
	mexErrMsgTxt ("symamd error!") ;
    }

    if (full)
    {
	mxDestroyArray (Ainput) ;
    }

    /* === Return the permutation vector ==================================== */

    plhs [0] = mxCreateDoubleMatrix (1, n_col, mxREAL) ;
    out_perm = mxGetPr (plhs [0]) ;
    for (i = 0 ; i < n_col ; i++)
    {
	/* symamd is 0-based, but MATLAB expects this to be 1-based */
	out_perm [i] = perm [i] + 1 ;
    }
    mxFree (perm) ;

    /* === Return the stats vector ========================================== */

    /* print stats if spumoni > 0 */
    if (spumoni > 0)
    {
	symamd_report (stats) ;
    }

    if (nlhs == 2)
    {
	plhs [1] = mxCreateDoubleMatrix (1, COLAMD_STATS, mxREAL) ;
	out_stats = mxGetPr (plhs [1]) ;
	for (i = 0 ; i < COLAMD_STATS ; i++)
	{
	    out_stats [i] = stats [i] ;
	}

	/* fix stats (5) and (6), for 1-based information on jumbled matrix. */
	/* note that this correction doesn't occur if symamd returns FALSE */
	out_stats [COLAMD_INFO1] ++ ; 
	out_stats [COLAMD_INFO2] ++ ; 
    }
}
