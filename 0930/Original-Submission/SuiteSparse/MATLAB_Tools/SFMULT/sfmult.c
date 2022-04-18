//==============================================================================
// y = A*x and variants, A is sparse, x is full
//==============================================================================

// y = ytrans (yconj (atrans (aconj (A)) * xtrans (xconj (x))))
//
//	where xtrans(x) is x or x.' and xconj(x) is x or conj(x) and likewise
//	for A and y.  To compute y = x*A, simply flip the use of *trans (see
//	dsmult below).  y = x*A is thus y = (A'*x')'

#include "sfmult.h"

//==============================================================================
//=== sfmult_invalid ===========================================================
//==============================================================================

void sfmult_invalid (void)
{
    mexErrMsgTxt ("Error using ==> sfmult\n"
	"Inner matrix dimensions must agree.") ;
}


//==============================================================================
//=== sfmult_yalloc ============================================================
//==============================================================================

// allocate Y as m-by-n, but do not initialize it

mxArray *sfmult_yalloc	// return Y
(
    Int m,
    Int n,
    int Ycomplex	// true if Y is complex
)
{
    // TODO: guard against integer overflow
    mxArray *Y = mxCreateDoubleMatrix (0, 0, Ycomplex ? mxCOMPLEX : mxREAL) ;
    mxFree (mxGetPr (Y)) ;
    mxSetPr (Y, mxMalloc (m * n * sizeof (double))) ;
    if (Ycomplex)
    {
	mxFree (mxGetPi (Y)) ;
	mxSetPi (Y, mxMalloc (m * n * sizeof (double))) ;
    }
    mxSetM (Y, m) ;
    mxSetN (Y, n) ;
    return (Y) ;
}


//==============================================================================
//=== sfmult_yzero =============================================================
//==============================================================================

// set Y to zero

mxArray *sfmult_yzero (mxArray *Y)
{
    Int n, i ;
    double *Yx, *Yz ;
    n = mxGetNumberOfElements (Y) ;
    Yx = mxGetPr (Y) ;
    for (i = 0 ; i < n ; i++)
    {
	Yx [i] = 0 ;
    }
    if (mxIsComplex (Y))
    {
	Yz = mxGetPi (Y) ;
	for (i = 0 ; i < n ; i++)
	{
	    Yz [i] = 0 ;
	}
    }
    return (Y) ;
}


//==============================================================================
//=== sfmult_walloc ============================================================
//==============================================================================

// Allocate workspace of size k*m

void sfmult_walloc
(
    Int k,
    Int m,
    double **Wx,	// real part (first k*m doubles)
    double **Wz		// imaginary part (next k*m doubles)
)
{
    // TODO Int overflow case
    Int wsize = k*m + 1 ; 
    *Wx = mxMalloc (wsize * sizeof (double)) ;   // TODO more if complex
    *Wz = *Wx + wsize ;
}


//==============================================================================
//=== sfmult ===================================================================
//==============================================================================

mxArray *sfmult		// returns y = A*x or variants
(
    const mxArray *A,
    const mxArray *X,
    int at,		// if true: trans(A)  if false: A
    int ac,		// if true: conj(A)   if false: A. ignored if A real
    int xt,		// if true: trans(x)  if false: x
    int xc,		// if true: conj(x)   if false: x. ignored if x real
    int yt,		// if true: trans(y)  if false: y
    int yc		// if true: conj(y)   if false: y. ignored if y real
)
{
    // TODO error if A not sparse, x sparse
    // TODO error if A not single or double, x not single or double

    if (at)
    {
	if (xt)
	{
	    if (yt)
	    {
		// y = (A'*x')'	    A is m-by-n, x is k-by-m, y is k-by-n
		return (sfmult_AT_XT_YT (A, X, ac, xc, yc)) ;
	    }
	    else
	    {
		// y = A'*x'	    A is m-by-n, x is k-by-m, y is n-by-k
		return (sfmult_AT_XT_YN (A, X, ac, xc, yc)) ;
	    }
	}
	else
	{
	    if (yt)
	    {
		// y = (A'*x)'	    A is m-by-n, x is m-by-k, y is k-by-n
		return (sfmult_AT_XN_YT (A, X, ac, xc, yc)) ;
	    }
	    else
	    {
		// y = A'*x	    A is m-by-n, x is m-by-k, y is n-by-k
		return (sfmult_AT_XN_YN (A, X, ac, xc, yc)) ;
	    }
	}
    }
    else
    {
	if (xt)
	{
	    if (yt)
	    {
		// y = (A*x')'	    A is m-by-n, x is k-by-n, y is k-by-m
		return (sfmult_AN_XT_YT (A, X, ac, xc, yc)) ;
	    }
	    else
	    {
		// y = A*x'	    A is m-by-n, x is k-by-n, y is m-by-k
		return (sfmult_AN_XT_YN (A, X, ac, xc, yc)) ;
	    }
	}
	else
	{
	    if (yt)
	    {
		// y = (A*x)'	    A is m-by-n, x is n-by-k, y is k-by-m
		return (sfmult_AN_XN_YT (A, X, ac, xc, yc)) ;
	    }
	    else
	    {
		// y = A*x	    A is m-by-n, x is n-by-k, y is m-by-k
		return (sfmult_AN_XN_YN (A, X, ac, xc, yc)) ;
	    }
	}
    }
}


//==============================================================================
//=== fsmult ===================================================================
//==============================================================================

mxArray *fsmult		// returns y = x*A or variants
(
    const mxArray *A,
    const mxArray *X,
    int at,		// if true: trans(A)  if false: A
    int ac,		// if true: conj(A)   if false: A. ignored if A real
    int xt,		// if true: trans(x)  if false: x
    int xc,		// if true: conj(x)   if false: x. ignored if x real
    int yt,		// if true: trans(y)  if false: y
    int yc		// if true: conj(y)   if false: y. ignored if y real
)
{
    return (sfmult (A, X, !at, ac, !xt, xc, !yt, yc)) ;
}
