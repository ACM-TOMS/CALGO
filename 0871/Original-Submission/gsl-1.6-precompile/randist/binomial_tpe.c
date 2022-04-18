#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/binomial_tpe.c
 * 
 * Copyright (C) 1996-2003 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

/* The binomial distribution has the form,

   f(x) =  n!/(x!(n-x)!) * p^x (1-p)^(n-x) for integer 0 <= x <= n
        =  0                               otherwise

   This implementation follows the public domain ranlib function
   "ignbin", the bulk of which is the BTPE (Binomial Triangle
   Parallelogram Exponential) algorithm introduced in
   Kachitvichyanukul and Schmeiser[1].  It has been translated to use
   modern C coding standards.

   If n is small and/or p is near 0 or near 1 (specifically, if
   n*min(p,1-p) < SMALL_MEAN), then a different algorithm, called
   BINV, is used which has an average runtime that scales linearly
   with n*min(p,1-p).

   But for larger problems, the BTPE algorithm takes the form of two
   functions b(x) and t(x) -- "bottom" and "top" -- for which b(x) <
   f(x)/f(M) < t(x), with M = floor(n*p+p).  b(x) defines a triangular
   region, and t(x) includes a parallelogram and two tails.  Details
   (including a nice drawing) are in the paper.

   [1] Kachitvichyanukul, V. and Schmeiser, B. W.  Binomial Random
   Variate Generation.  Communications of the ACM, 31, 2 (February,
   1988) 216.

   Note, Bruce Schmeiser (personal communication) points out that if
   you want very fast binomial deviates, and you are happy with
   approximate results, and/or n and n*p are both large, then you can
   just use gaussian estimates: mean=n*p, variance=n*p*(1-p).

   This implementation by James Theiler, April 2003, after obtaining
   permission -- and some good advice -- from Drs. Kachitvichyanukul
   and Schmeiser to use their code as a starting point, and then doing
   a little bit of tweaking.

   Additional polishing for GSL coding standards by Brian Gough.  */

#define SMALL_MEAN 14           /* If n*p < SMALL_MEAN then use BINV
                                   algorithm. The ranlib
                                   implementation used cutoff=30; but
                                   on my computer 14 works better */

#define BINV_CUTOFF 110         /* In BINV, do not permit ix too large */

#define FAR_FROM_MEAN 20        /* If ix-n*p is larger than this, then
                                   use the "squeeze" algorithm.
                                   Ranlib used 20, and this seems to
                                   be the best choice on my machine as
                                   well */

#define LNFACT(x) gsl_sf_lnfact(x)

inline static MpIeee Stirling(MpIeee y1)
{
  MpIeee y2=  y1 * y1;
  MpIeee s= 
    (MpIeee( "13860.0" ) -
     (MpIeee( "462.0" ) - (MpIeee( "132.0" ) - (MpIeee( "99.0" ) - MpIeee( "140.0" ) / y2) / y2) / y2) / y2) / y1 / MpIeee( "166320.0" );
  return s;
}

unsigned int
 gsl_ran_binomial_tpe(const gsl_rng * rng, MpIeee pp, unsigned int  n)
{
  int  ix;                       /* return value */

  const MpIeee p=  (pp > 0.5) ? 1 - pp : pp;    /* choose p=min(pp,1-pp) */
  const MpIeee q=  1 - p;
  const MpIeee s=  p / q;
  const MpIeee np=  n * p;

  if (n == 0)
    return 0;

  /* Inverse cdf logic for small mean (BINV in K+S) */

  if (np < SMALL_MEAN)
    {
      MpIeee f0=  gsl_pow_int (q, n);   /* f(x), starting with x=0 */

      while (1)
        {
          /* This while(1) loop will almost certainly only loop once; but
           * if u=1 to within a few epsilons of machine precision, then it
           * is possible for roundoff to prevent the main loop over ix to
           * achieve its proper value.  following the ranlib implementation,
           * we introduce a check for that situation, and when it occurs,
           * we just try again.
           */

          MpIeee f=  f0;
          MpIeee u=  gsl_rng_uniform (rng);

          for (ix = 0; ix <= BINV_CUTOFF; ++ix)
            {
              if (u < f)
                goto Finish;
              u -= f;
              /* Use recursion f(x+1) = f(x)*[(n-x)/(x+1)]*[p/(1-p)] */
              f *= s * (n - ix) / (ix + MpIeee( "1" ));
            }

          /* It should be the case that the 'goto Finish' was encountered
           * before this point was ever reached.  But if we have reached
           * this point, then roundoff has prevented u from decreasing
           * all the way to zero.  This can happen only if the initial u
           * was very nearly equal to 1, which is a rare situation.  In
           * that rare situation, we just try again.
           *
           * Note, following the ranlib implementation, we loop ix only to
           * a hardcoded value of SMALL_MEAN_LARGE_N=110; we could have
           * looped to n, and 99.99...% of the time it won't matter.  This
           * choice, I think is a little more robust against the rare
           * roundoff error.  If n>LARGE_N, then it is technically
           * possible for ix>LARGE_N, but it is astronomically rare, and
           * if ix is that large, it is more likely due to roundoff than
           * probability, so better to nip it at LARGE_N than to take a
           * chance that roundoff will somehow conspire to produce an even
           * larger (and more improbable) ix.  If n<LARGE_N, then once
           * ix=n, f=0, and the loop will continue until ix=LARGE_N.
           */
        }
    }
  else
    {
      /* For n >= SMALL_MEAN, we invoke the BTPE algorithm */

      int  k;

      MpIeee ffm=  np + p;      /* ffm = n*p+p             */
      int  m=  floor(ffm).toInt();        /* m = int floor[n*p+p]    */
      MpIeee fm=  m;            /* fm = double m;          */
      MpIeee xm=  fm + MpIeee( "0.5" );     /* xm = half integer mean (tip of triangle)  */
      MpIeee npq=  np * q;      /* npq = n*p*q            */

      /* Compute cumulative area of tri, para, exp tails */

      /* p1: radius of triangle region; since height=1, also: area of region */
      /* p2: p1 + area of parallelogram region */
      /* p3: p2 + area of left tail */
      /* p4: p3 + area of right tail */
      /* pi/p4: probability of i'th area (i=1,2,3,4) */

      /* Note: magic numbers 2.195, 4.6, 0.134, 20.5, 15.3 */
      /* These magic numbers are not adjustable...at least not easily! */

      MpIeee p1=  floor (MpIeee( "2.195" ) * sqrt (npq) - MpIeee( "4.6" ) * q) + MpIeee( "0.5" );

      /* xl, xr: left and right edges of triangle */
      MpIeee xl=  xm - p1;
      MpIeee xr=  xm + p1;

      /* Parameter of exponential tails */
      /* Left tail:  t(x) = c*exp(-lambda_l*[xl - (x+0.5)]) */
      /* Right tail: t(x) = c*exp(-lambda_r*[(x+0.5) - xr]) */

      MpIeee c=  MpIeee( "0.134" ) + MpIeee( "20.5" ) / (MpIeee( "15.3" ) + fm);
      MpIeee p2=  p1 * (MpIeee( "1.0" ) + c + c);

      MpIeee al=  (ffm - xl) / (ffm - xl * p);
      MpIeee lambda_l=  al * (MpIeee( "1.0" ) + MpIeee( "0.5" ) * al);
      MpIeee ar=  (xr - ffm) / (xr * q);
      MpIeee lambda_r=  ar * (MpIeee( "1.0" ) + MpIeee( "0.5" ) * ar);
      MpIeee p3=  p2 + c / lambda_l;
      MpIeee p4=  p3 + c / lambda_r;

      MpIeee var;MpIeee  accept;
      MpIeee u;MpIeee  v;              /* random variates */

    TryAgain:

      /* generate random variates, u specifies which region: Tri, Par, Tail */
      u = gsl_rng_uniform (rng) * p4;
      v = gsl_rng_uniform (rng);

      if (u <= p1)
        {
          /* Triangular region */
          ix = (int) MpIeee(xm - p1 * v + u).toInt();
          goto Finish;
        }
      else if (u <= p2)
        {
          /* Parallelogram region */
          MpIeee x=  xl + (u - p1) / c;
          v = v * c + MpIeee( "1.0" ) - fabs (x - xm) / p1;
          if (v > MpIeee( "1.0" ) || v <= MpIeee( "0.0" ))
            goto TryAgain;
          ix = x.toInt();
        }
      else if (u <= p3)
        {
          /* Left tail */
          ix = MpIeee(xl + log (v) / lambda_l).toInt();
          if (ix < 0)
            goto TryAgain;
          v *= ((u - p2) * lambda_l);
        }
      else
        {
          /* Right tail */
          ix = MpIeee(xr - log (v) / lambda_r).toInt();
          if (ix > (double)n)
            goto TryAgain;
          v *= ((u - p3) * lambda_r);
        }

      /* At this point, the goal is to test whether v <= f(x)/f(m) 
       *
       *  v <= f(x)/f(m) = (m!(n-m)! / (x!(n-x)!)) * (p/q)^{x-m}
       *
       */

      /* Here is a direct test using logarithms.  It is a little
       * slower than the various "squeezing" computations below, but
       * if things are working, it should give exactly the same answer
       * (given the same random number seed).  */

#ifdef DIRECT
      var = log (v);

      accept =
        LNFACT (m) + LNFACT (n - m) - LNFACT (ix) - LNFACT (n - ix)
        + (ix - m) * log (p / q);

#else /* SQUEEZE METHOD */

      /* More efficient determination of whether v < f(x)/f(M) */

      k = abs (ix - m);

      if (k <= FAR_FROM_MEAN)
        {
          /* 
           * If ix near m (ie, |ix-m|<FAR_FROM_MEAN), then do
           * explicit evaluation using recursion relation for f(x)
           */
          MpIeee g=  (n + MpIeee( "1" )) * s;
          MpIeee f=  MpIeee( "1.0" );

          var = v;

          if (m < ix)
            {
              int  i;
              for (i = m + 1; i <= ix; i++)
                {
                  f *= (g / i - s);
                }
            }
          else if (m > ix)
            {
              int  i;
              for (i = ix + 1; i <= m; i++)
                {
                  f /= (g / i - s);
                }
            }

          accept = f;
        }
      else
        {
          /* If ix is far from the mean m: k=ABS(ix-m) large */

          var = log (v);

          if (k < npq / 2 - 1)
            {
              /* "Squeeze" using upper and lower bounds on
               * log(f(x)) The squeeze condition was derived
               * under the condition k < npq/2-1 */
              MpIeee amaxp= 
                k / npq * ((k * (k / MpIeee( "3.0" ) + MpIeee( "0.625" )) + (MpIeee( "1.0" ) / MpIeee( "6.0" ))) / npq + MpIeee( "0.5" ));
              MpIeee ynorm=  -(k * k / (MpIeee( "2.0" ) * npq));
              if (var < ynorm - amaxp)
                goto Finish;
              if (var > ynorm + amaxp)
                goto TryAgain;
            }

          /* Now, again: do the test log(v) vs. log f(x)/f(M) */

#if USE_EXACT
          /* This is equivalent to the above, but is a little (~20%) slower */
          /* There are five log's vs three above, maybe that's it? */

          accept = LNFACT (m) + LNFACT (n - m)
            - LNFACT (ix) - LNFACT (n - ix) + (ix - m) * log (p / q);

#else /* USE STIRLING */
          /* The "#define Stirling" above corresponds to the first five
           * terms in asymptoic formula for
           * log Gamma (y) - (y-0.5)log(y) + y - 0.5 log(2*pi);
           * See Abramowitz and Stegun, eq 6.1.40
           */

          /* Note below: two Stirling's are added, and two are
           * subtracted.  In both K+S, and in the ranlib
           * implementation, all four are added.  I (jt) believe that
           * is a mistake -- this has been confirmed by personal
           * correspondence w/ Dr. Kachitvichyanukul.  Note, however,
           * the corrections are so small, that I couldn't find an
           * example where it made a difference that could be
           * observed, let alone tested.  In fact, define'ing Stirling
           * to be zero gave identical results!!  In practice, alv is
           * O(1), ranging 0 to -10 or so, while the Stirling
           * correction is typically O(10^{-5}) ...setting the
           * correction to zero gives about a 2% performance boost;
           * might as well keep it just to be pendantic.  */

          {
            MpIeee x1=  ix + MpIeee( "1.0" );
            MpIeee w1=  n - ix + MpIeee( "1.0" );
            MpIeee f1=  fm + MpIeee( "1.0" );
            MpIeee z1=  n + MpIeee( "1.0" ) - fm;

            accept = xm * log (f1 / x1) + (n - m + MpIeee( "0.5" )) * log (z1 / w1)
              + (ix - m) * log (w1 * p / (x1 * q))
              + Stirling (f1) + Stirling (z1) - Stirling (x1) - Stirling (w1);
          }
#endif
#endif
        }


      if (var <= accept)
        {
          goto Finish;
        }
      else
        {
          goto TryAgain;
        }
    }

Finish:

  return (pp > 0.5) ? n - ix : ix;
}
