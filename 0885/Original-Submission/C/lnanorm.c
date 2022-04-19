/*  This implementation for the log of the standard normal distribution
    is by Jean Marie Linhart
	  StataCorp LP

    It is based on the algorithm for the standard normal distribution
    by W. J. Cody published in
    Algorithm 715, Collected Algorithms from ACM.
       This work published in Transactions on Mathematical Software,
       Vol. 19, No. 1, March, 1993, pp. 22-32.
      
      Includes changes given in remark by Price, TOMS 22 (2)
*/
/*------------------------------------------------------------------

 This function evaluates the log of the normal distribution function:

                                 / x    
                        1       |       -t*t/2
      ln(P(x)) =ln( ----------- |      e       dt )
                    sqrt(2 pi)  |
                                /-oo

   The main computation evaluates near-minimax approximations
   derived from those in "Rational Chebyshev approximations for
   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
   This transportable program uses rational functions that
   theoretically approximate the normal distribution function to
   at least 18 significant decimal digits.  The accuracy achieved
   depends on the arithmetic system, the compiler, the intrinsic
   functions, and proper selection of the machine-dependent
   constants.

 *******************************************************************
 *******************************************************************

 Explanation of machine-dependent constants.  Let

   xmin  = the smallest positive floating-point number.

   Then the following machine-dependent constants must be declared 
   in DATA statements.  IEEE values are provided as a default.

   eps   = argument below which anorm(x) may be represented by
           0.5  and above which  x*x  will not underflow.
           A conservative value is the largest machine number x
           such that   1.0 + x = 1.0   to machine precision.
   Approximate values for some important machines are:

                          xmin        eps 

  CDC 7600      (S.P.)  3.13E-294   7.11E-15
  CRAY-1        (S.P.)  4.58E-246   7.11E-157
  IEEE (IBM/XT,
    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8 
  IEEE (IBM/XT,
    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16
  IBM 195       (D.P.)  5.40D-79    1.39D-17
  VAX D-Format  (D.P.)  2.94D-39    1.39D-17
  VAX G-Format  (D.P.)  5.56D-309   1.11D-16

 *******************************************************************
 *******************************************************************

 Error returns

  The program returns  anorm = 0     for  arg > LNNORM_MAX_X



 Intrinsic functions required are:

     fabs, trunc or floor, exp


  Original normal distribution Algorithm Author: W. J. Cody
          Mathematics and Computer Science Division
          Argonne National Laboratory
          Argonne, IL 60439
  Author for the log of the normal distribution:
	Jean Marie Linhart
	StataCorp LP

  Latest modification: January 3, 2008

 -----------------------------------------------------------------*/
#include "norminc.h"

double lnanorm(double arg)
{
	int i;
	double del, eps=1.11e-16,
		half=0.50e0, one=1.0e0, result, sixten=1.6e1,
		x, xden, xnum, y, xsq,
		zero=0.0e0;
/* ------------------------------------------------------------------
   Mathematical constants

   sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), and
   thrsh is the argument for which anorm = 0.75.
   ------------------------------------------------------------------*/
/* Next two lines represent limits on double precision arithmetic
   for this problem.  One may want to test on one's own machine
   and change -- the chips hold data on a register with more
   precision than double, and this can effect results obtained.

   LNNORM_MAX_X
   LNNORM_MIN_X
	defined in norminc.h
*/
/*------------------------------------------------------------------
   Coefficients for approximation in first interval
  ------------------------------------------------------------------*/
	const static double a[6] = {
		2.2352520354606839287e0,
		1.6102823106855587881e2,
		1.0676894854603709582e3,
		1.8154981253343561249e4,
		6.5682337918207449113e-2
	};

	const static double b[5] = {
		4.7202581904688241870e1,
		9.7609855173777669322e2,
		1.0260932208618978205e4,
		4.5507789335026729956e4 
	};
/*------------------------------------------------------------------
   Coefficients for approximation in second interval
  ------------------------------------------------------------------*/
	const static double c[9] = {
		3.9894151208813466764e-1,
		8.8831497943883759412e0,
		9.3506656132177855979e1,
		5.9727027639480026226e2,
		2.4945375852903726711e3,
		6.8481904505362823326e3,
		1.1602651437647350124e4,
		9.8427148383839780218e3,
		1.0765576773720192317e-8 
	};

	const static double d[8] = {
		2.2266688044328115691e1,
		2.3538790178262499861e2,
		1.5193775994075548050e3,
		6.4855582982667607550e3,
		1.8615571640885098091e4,
		3.4900952721145977266e4,
		3.8912003286093271411e4,
		1.9685429676859990727e4 
	};
/*------------------------------------------------------------------
   Coefficients for approximation in third interval
  ------------------------------------------------------------------*/
	const static double p[6] = {
		2.1589853405795699e-1,
		1.274011611602473639e-1,
		2.2235277870649807e-2,
		1.421619193227893466e-3,
		2.9112874951168792e-5,
		2.307344176494017303e-2
	};

	const static double q[5] = {
		1.28426009614491121e0,
		4.68238212480865118e-1,
		6.59881378689285515e-2,
		3.78239633202758244e-3,
		7.29751555083966205e-5 
	};
/*------------------------------------------------------------------*/
	if (arg > LNNORM_MAX_X) return(0.0e0);
/* Call through to lnnorm to get correct behavior in the lower tail */
	if (arg < LNANORM_MIN_X) return(lnnorm(arg));
/*------------------------------------------------------------------*/
	x = arg;
	y = fabs(x);
	if (y <= thrsh) {
/*------------------------------------------------------------------
   Evaluate  lnanorm  for  |x| <= 0.67448975
  ------------------------------------------------------------------*/
		xsq = zero;
		if (y > eps) xsq = x * x;
		xnum = a[4]*xsq;
		xden = xsq;
		for(i=0; i< 3; i++) {
			xnum = (xnum + a[i]) * xsq;
			xden = (xden + b[i]) * xsq;
		}
		result = x * (xnum + a[3]) / (xden + b[3]);
		result = log(half + result);
	}
/*------------------------------------------------------------------
   Evaluate  lnanorm  for 0.67448975 <= |x| <= sqrt(32)
  ------------------------------------------------------------------*/
	else if (y <= root32) {
			xnum = c[8]*y;
			xden = y;
			for(i=0; i<7; i++) {
				xnum = (xnum + c[i]) * y;
				xden = (xden + d[i]) * y;
			}
			result = (xnum + c[7]) / (xden + d[7]);
			xsq = TRUNC(y*sixten)/sixten;
			del = (y-xsq)*(y+xsq);
			result = exp(-xsq*xsq*half)*exp(-del*half)*result;
			if (x > zero) {
				if (result < .1) { /* close to 0 */
					result = LOG1P(-result);
				}
				else { /* far from 0 */
					result = log(one-result);
				}
			}
			else {
				result = log(result);
			}
	}
/*------------------------------------------------------------------
   Evaluate  lnanorm  for |x| > sqrt(32)
  ------------------------------------------------------------------*/
	else {
		result = zero;
		xsq = one / (x * x);
		xnum = p[5]*xsq;
		xden = xsq;
		for(i=0; i<4; i++) {
			xnum = (xnum + p[i]) * xsq;
			xden = (xden + q[i]) * xsq;
		}
		result = xsq *(xnum + p[4]) / (xden + q[4]);
		result = (sqrpi -  result) / y;
		xsq = TRUNC(x*sixten)/sixten;
		del = (x-xsq)*(x+xsq);
		result = exp(-xsq*xsq*half)*exp(-del*half)*result;
		if (x > zero) {
			result = LOG1P(-result);
		}
		else {
			result = log(result);
		}
	}
	return(result);
}
