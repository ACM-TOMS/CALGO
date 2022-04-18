/* C code for the lnnorm(x) function to return the logarithm of the normal
   distribution.

   It is based on the attached code for the normal distribution, originally
   written at StataCorp LP, and based on Algorithm 304 and the remarks by
   Adams and Holmgren.

   This code is written by Jean Marie Linhart
                           StataCorp LP
                           jlinhart@stata.com

   Last modified January 4, 2008
*/

#include "norminc.h"

/* The following is a modification of Algorithm 304 to give the log of the
 * normal distribution function at x.
 *
 * It takes a double precision argument and returns double precision
 */

double lnnorm(double z)
{
        int upper, lower ;

        double z2, y, s, p1, q1, p2, q2, t, a1, a2 ;
        double n, m ;

	if (z==0.0e0) return(log(0.50e0)) ;

/* Next two lines represent limits on double precision arithmetic
   for this problem.  One may want to test on one's own machine
   and change -- the chips hold data on a register with more
   precision than double, and this can effect results obtained.

   LNNORM_MAX_X
   LNNORM_MIN_X
	defined in norminc.h
*/
	if (z > LNNORM_MAX_X) return(0.0e0); 
        if (z <= LNNORM_MIN_X) return(-0.5e0*z*z);


/*
        In the original algorithm, the logical variable upper served a
        dual purpose.  On input, it indicated whether the user wanted
        F(z) or 1-F(z).  Inside the routine, it combined this with
        information on the positivity of z.
        In this version, the algorithm always returns F(z)
*/

        if (z<0.0e0) {
                z= -z ;
                lower=1 ;
        }
        else    lower=0 ;
        upper = !lower ;

        z2 = z*z ;

/*
        y is the standard normal density function evaluated at z.
*/
        y = exp(-0.5*z2) / sqr2pi ;
        n = y/z ;
/* The original Algorithm 304 had z < 2.32 or z < 3.5 depending on
   whether you were doing upper or lower.  
   This was changed for (slightly) greater accuracy */
        if (!( ( z>2.00e0))) {
/* This is the series representation for F(z) 
   where z lies near center of the distribution */
                z *= y ;
                s=z ;
                t=0.0e0 ;
/* Comparing s!=t is done here to insure that this calculation continues until
the accuracy of the machine is reached.  At some point, the value is not
being updated in successive iterations, that is time to quit. */
                for (n=3.0e0;s!=t;n+=2.0e0) {
                        t=s ;
                        z *= ((z2)/ (n)) ;
                        s+=z ;
                }
                if (lower) return(log(0.50e0-s)) ;
                return(log(0.50e0+s)) ;
        }
/*
	This is the continued fraction representation for Q(z) which
        is used to evaluate F(z) when z lies in one of the tails
        of this distribution.  This section contains the modifications
        suggested by Homgren and Adams.
*/

        a1=2.0e0 ;
        a2=0.0e0 ;
        n=z2+3.0e0 ;
        p1=1.0e0 ;
        q1=z ;
        p2=(n-1.0e0) ;
        q2=n*z ;
        m = p1/q1 ;
        t = p2/q2 ;
        s = m ;                 /* dummy assignment to not stop for     */
/* Comparing  m!=t and s!=t is done here to insure that this calculation
continues until the accuracy of the machine is reached.  At some point,
the value is not being updated in successive iterations,
that is time to quit. */
        for (n+=4.0e0; m!=t && s!=t; n+=4.0e0) {
                a1 -= 8.0 ;
                a2 += a1 ;
                s = a2*p1 + n*p2 ;
                p1=p2 ;
                p2=s ;
                s = a2*q1 + n*q2 ;
                q1=q2 ;
                q2=s ;
                s=m ;
                m=t ;
                if (q2>1.0e30) {
                        p1 /= 1.0e30 ;
                        p2 /= 1.0e30 ;
                        q1 /= 1.0e30 ;
                        q2 /= 1.0e30 ;
                }
 		t = p2/q2;
        }
        t = lower ? log(t) - 0.5*z2 - log(sqr2pi) : LOG1P(-y*t);
        return(t) ;
}

