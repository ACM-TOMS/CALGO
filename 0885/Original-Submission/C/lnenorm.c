/*  This algorithm for the log of the standard normal distribution is
    written by
	Jean Marie Linhart
	StataCorp LP

    Last modified: October 26, 2007
*/
#include "norminc.h"
#ifdef HAS_ERF

double lnenorm(double x)
{
	double ex, ans=0;
	int pos=1;

/* Next two lines represent limits on double precision arithmetic
   for this problem.  One may want to test on one's own machine
   and change -- the chips hold data on a register with more
   precision than double, and this can effect results obtained.

   LNNORM_MAX_X
   LNNORM_MIN_X
	defined in norminc.h
*/
	if (x > LNNORM_MAX_X) return(0.0e0);
/* Call through to lnnorm to get correct behavior in the lower tail */
	if (x < LNANORM_MIN_X) return(lnnorm(x));

	if (x<0) {
		pos = 0;
		x = -x;
	}

	if (x <= 0.5) {
		ex =0.5*erf(x/sqrt2);
		if (pos) {
			ans = log(0.5 + ex);
		}
		else {
			ans = log(0.5 - ex);
		}
	}
	else {
		if (pos) {
			ex = 0.5*erfc(x/sqrt2);
			if (ex < .1) {
				ans = LOG1P(-ex);
			}
			else {
				ans = log(1.0e0-ex);
			}
		}
		else {
			ex = 0.5*erfc(x/sqrt2);
			ans = log(ex);
		}
	}
	
        return(ans);
}
#endif

