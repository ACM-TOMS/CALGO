/*  This little function for log(1+x) is based on the Taylor series.
    It is meant to be called only with small x.
    Written by:
	Jean Marie Linhart
	StataCorp LP

    Last modified October 26, 2007
*/

#include "math.h"
#include "norminc.h"
/*
 * log1px takes a double and returns a double.
 * It is a Taylor series expansion of log(1+x).
 * x is presumed to be < 1.  As I have called it, x < .1,
 * and so I know the algorithm will terminate quickly.
 * The closer x is to 1, the slower this will be.
 */
double log1px(double x)
{
        int n, sn;
        double xn, ans, oans, term, eps;
        
	if (!(fabs(x) < 1.0e0)) {
		return(0.0e0/0.0e0); /* NaN */
	}
	term = ans= oans = x;
	oans = ans + (double) 1.0e0;
        n = 1;
        sn = 1;
        xn = x;
/* Comparing ans!=oans is done here to insure that this calculation
continues until the accuracy of the machine is reached.  At some point,
the value is not being updated in successive iterations, that is time to
quit. */
        while (ans != oans ) {
		oans = ans;
                sn *= -1;
                xn *= x;
                term= ((double)sn/(double)++n)*xn;
                ans += term;
	}
        return(ans);
}


