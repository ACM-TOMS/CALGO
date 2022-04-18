/*  This algorithm for the standard normal distribution is written by
	Jean Marie Linhart
	StataCorp LP

    Last modified: October 16, 2007
*/
#include "norminc.h"

#ifdef HAS_ERF

double enorm(double x)
{
	double ex, ans=0;
	int pos=1;

	if (x<0) {
		pos = 0;
		x = -x;
	}

	if (x <= 0.5) {
		ex =0.5*erf(x/sqrt2);
		if (pos) {
			ans = 0.5 + ex;
		}
		else {
			ans = 0.5-ex;
		}
	}
	else {
		if (pos) {
			ex = 0.5*erfc(x/sqrt2);
			ans = 1.0e0-ex;
		}
		else {
			ex = 0.5*erfc(x/sqrt2);
			ans = ex;
		}
	}
	
        return(ans);
}
#endif

