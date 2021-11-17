/* C code for the norm(x) function to return the norm distribution
   code is based on Algorithm 304 of the CACM

   This code was originally written at StataCorp LP.

*/

#include "norminc.h"

double norm(double z)
{
        int upper, lower ;

        double z2, y, s, p1, q1, p2, q2, t, a1, a2 ;
        double n, m ;

        if (z==0.0e0) return(0.5e0) ;


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
        y = (exp(-0.5*z2) / sqr2pi) ;
        n = (y/z) ;
        if (!(lower || (1.0e0-n)!=1.0e0)) return(1.0e0) ;
        if (!(upper || (n!=0.0e0))) return(0.0e0) ;

/* In the original Algorithm 304 code, this next condition had
   z < 2.32 or z < 3.5 depending on whether you were doing upper or lower.
   This was changed for (slightly) greater accuracy. */
	if (! (z>2.00) ) {
/* This is the series representation for F(z) where z lies near
   the center of the distribution 
*/
                z *= y ;
                s=z ;
                t=0.0e0 ;
/* Comparing s!=t is done here to insure that this calculation continues until
the accuracy of the machine is reached.  At some point, the value is not
being updated in successive iterations, that is time to quit. */
                for (n=3.0e0;s!=t;n+=2.0e0) {
                        t=s ;
                        z *= (z2/n) ;
                        s+=z ;
                }
                if (lower) return(0.5e0-s) ;
                return(0.5e0+s) ;
        }
/*
	This is the continued fraction representation for Q(z) which is
        used to evaluate F(z) when z lies in one of the tails
        of this distribution.  This section contains the modifications
        suggested by Holmgren and Adams.
*/

        a1=2.0e0 ;
        a2=0.0e0 ;
        n=z2+3.0e0 ;
        p1=y ;
        q1=z ;
        p2=(n-1.0e0)*y ;
        q2=n*z ;
        m = (p1/q1) ;
        t = (p2/q2) ;
        if (!lower) {
                m=1.0e0-m ;
                t=1.0e0-t ;
        }
        s = m ;                 /* dummy assignment to not stop for     */
/* Comparing m!=t and s!=t is done here to insure that this calculation
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
                t=lower ? (p2/q2) : 1.0e0-(p2/q2) ;
        }
        return(t) ;
}

