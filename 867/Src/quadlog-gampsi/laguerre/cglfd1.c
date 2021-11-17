/***********************************************************************
This is a hand translation (in order to have clean C code) of Fortran
PROGRAM tglfd1 for testing the C/C++ interface to the log quadrature
package.

************************************************************************
*     (Test Gauss-Laguerre Quadrature w/ Function and Derivative values)
*
*     This is the driver program for subroutine glqfd which
*     returns the nodes and weights necessary for the log-Laguerre
*     quadrature.  The test functions are x^j, j = 0 to n.
*
*     Normal input is on stdin, and normal output on stdout.
*
*     However, there is also an auxiliary file created called simply
*     "output.dat"; it contains data suitable for comparing against
*     similar files from other programs in an integral-independent
*     way.
*
*     Lines in this file are either comment lines (beginning with a
*     sharp sign), or blank or empty lines, or whitespace-separated
*     integral data lines of the form
*
*     np p(1) p(2) ... p(np) neval result opt-relerr opt-abserr
*
*     Here, np is the number of parameters that follow in the p(*)
*     values, neval is the number of function evaluations (0 if
*     unknown or unavailable), result is the floating-point computed
*     value of the integral, opt-relerr is an optional relative error
*     estimate, and opt-abserr is an optional absolute error estimate.
*     If opt-relerr is omitted, then opt-abserr cannot be specified.
*     A zero value for either implies that the value is unknown.
*
*     Normally, comment lines will document what integral is
*     evaluated, and relate any parameters in the integrand to the
*     array elements p(*).  It is acceptable for np to be 0: no p(*)
*     elements are then provided.
*
*     The availability of data files in this standard format makes it
*     relatively easy to compare results of different integrators.  In
*     particular, high-precision values can be computed in symbolic
*     algebra systems, such as Maple, Mathematica, Axiom, Reduce,
*     muPAD, ..., and used to evaluate the accuracy of results from
*     other integrators.
*
************************************************************************

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../gjl.h"

#undef alpha			/* defined by some Compaq/DEC Alpha compilers, sigh... */

#define MAXPTS  1049

static const fortran_double_precision half = 0.5e+00;
static const fortran_double_precision one = 1.0e+00;
static const fortran_double_precision two = 2.0e+00;
static const fortran_double_precision zero = 0.0e+00;
static fortran_integer i0 = 0;
static fortran_integer i1 = 1;
static fortran_integer maxpts = MAXPTS;

#define dfloat Dfloat		/* avoid conflict with ../common/dfloat.f */

static fortran_double_precision
dabs(fortran_double_precision darg)
{
    return ((darg < zero) ? -(darg) : darg);
}

static fortran_double_precision
dlog(fortran_double_precision darg)
{
    return log(darg);
}

static fortran_double_precision
dlog2(fortran_double_precision darg)
{
    return (dlog(darg) / dlog(two));
}

static fortran_double_precision
dfloat(fortran_integer arg)
{
    return ((fortran_double_precision) (arg));
}

static fortran_double_precision
f(fortran_double_precision xarg, fortran_integer narg,
  fortran_double_precision aarg)
{
    aarg = aarg;		/* NO-OP */
    return (pow(xarg, (fortran_double_precision) narg));
}

static fortran_double_precision
fprime(fortran_double_precision xarg, fortran_integer narg,
       fortran_double_precision aarg)
{
    aarg = aarg;		/* NO-OP */
    return (dfloat(narg) * pow(xarg, (fortran_double_precision) (narg - 1)));
}

static fortran_integer
idnint(fortran_double_precision x)
{
    if (x < zero)
	return (fortran_integer)(x - half);
    else if (x >= zero)
	return (fortran_integer)(x + half);
    else				/* must be NaN */
	return (fortran_integer)x;
}

static void
prthdr(FILE * outfile, const char *header)
{
    fortran_double_precision ulp;

    ulp = deps(&one);
    (void)fprintf(outfile, "  Test integral: %s\n", header);
    (void)fprintf(outfile, "\n  1 ulp = %10.2e\n", ulp);
    (void)fprintf(outfile,
		  "\n  Floating-point precision appears to be (1 - log2(ulp) =)%3d bits\n",
		  1 - idnint(dlog2(ulp)));
}

int
main(void)
{
    fortran_double_precision alpha;
    fortran_double_precision b;
    fortran_double_precision c;
    fortran_double_precision deltax[MAXPTS];
    fortran_double_precision dw[MAXPTS];
    fortran_double_precision errbit;
    fortran_double_precision exact;
    fortran_double_precision fn;
    fortran_double_precision relerr;
    fortran_double_precision result;
    fortran_double_precision temp;
    fortran_double_precision ulp;
    fortran_double_precision v[MAXPTS][2];
    fortran_double_precision w[MAXPTS];
    fortran_double_precision x[MAXPTS];
    fortran_double_precision zinit[1];

    FILE *stddat;

    fortran_integer ierr;
    fortran_integer neval;
    fortran_integer np;
    fortran_integer nquad;
    int i;
    int p;
    int pmax;

    /* Initialize all local floating-point variables to NaN, or at
       least an approximation thereto. */
    zinit[0] = getnan();
    dcopy(&maxpts, zinit, &i0, deltax, &i1);
    dcopy(&maxpts, zinit, &i0, dw, &i1);
    dcopy(&maxpts, zinit, &i0, w, &i1);
    dcopy(&maxpts, zinit, &i0, x, &i1);
    alpha = zinit[0];
    b = zinit[0];
    c = zinit[0];
    errbit = zinit[0];
    exact = zinit[0];
    fn = zinit[0];
    relerr = zinit[0];
    result = zinit[0];
    temp = zinit[0];
    ulp = zinit[0];

    stddat = fopen("output.dat", "w");

    if (stddat == (FILE *) NULL)
    {
	(void)fprintf(stderr, "open failure on output.dat\n");
	exit(EXIT_FAILURE);
    }

    (void)fprintf(stddat, "### Numerical integration with glqfd()\n");
    (void)fprintf(stddat, "###\n");

    (void)fprintf(stddat,
		  "### int(x^p * x^alpha * exp(-x) * ln(x), x = 0..infinity)\n");
    (void)fprintf(stddat, "###\n");

    (void)fprintf(stddat,
		  "### Line format: 2 p alpha neval result relerr abserr\n");
    (void)fprintf(stddat, "###\n");

    ulp = deps(&one);
    prthdr(stdout, "int(x^p * x^alpha * exp(-x) * ln(x), x = 0..infinity)");
    while (scanf("%d %le %d\n", &nquad, &alpha, &pmax) == 3)
    {
	if (alpha <= -one)
	{
	    (void)fprintf(stderr, "  Illegal value of alpha: %11.2e <= -1\n",
			  alpha);
	    continue;
	}
	(void)printf("\n  alpha = %10.6f  nquad = %6d\n", alpha, nquad);
	ierr = 0;
	glqfd(x, w, dw, deltax, &alpha, &nquad, &ierr);
	if (ierr != 0)
	{
	    (void)fprintf(stderr, "  ERROR: glqfd() returns ierr = %10d\n", ierr);
	    continue;
	}

	(void)printf("\n     i      x(i)                 w(i)"
		     "                     dw(i)"
		     "                    deltax(i)\n");
	for (i = 0; i < nquad; ++i)

	    (void)printf("  %4d %20.15f%25.15e%25.15e%25.15e\n",
			 i + 1, x[i], w[i], dw[i], deltax[i]);

	/* Evaluate integral as sum over nodes */

	b = dgamma((temp = one + alpha, &temp));
	c = dpsi((temp = one + alpha, &temp));
	(void)printf("\n      p     Quadrature Result"
		     "       Exact Integral      Rel. Error RelE (ULPs)  Err (bits)\n");
	for (p = 0; p <= pmax; ++p)
	{
	    np = p;
	    for (i = 0; i < nquad; ++i)
	    {
		v[i][0] = dw[i] * f(x[i], p, alpha);
		v[i][1] = fprime(x[i], p, alpha) * deltax[i];
	    }

	    neval = 2 * nquad;
	    result = dvsum((fortran_double_precision *)&v[0][0], &neval);

	    /* Analytic calculation of test integral */

	    fn = dfloat(p+1) + alpha;
	    exact = b * c;
	    relerr = (exact - result) / exact;
	    if (relerr == zero)
		errbit = zero;
	    else
		errbit = dlog2(dabs(relerr / ulp));

	    (void)printf("  %5d%24.15e%24.15e%11.2e  %10.2f%10.2f\n",
			 p, result, exact, relerr, relerr / ulp, errbit);
	    (void)fprintf(stddat, "%1d %3d %27.18e %10d %27.18e%10.2e%10.2e\n",
			  2, p, (fn - one), neval, result, relerr,
			  (exact - result));
	    b = b * fn;
	    c = c + one / fn;
	}
    }
    (void)printf("\n  Done\n");
    (void)fclose(stddat);
    return (EXIT_SUCCESS);
}
