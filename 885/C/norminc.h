#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* if your compiler does not have an erf function remove the next
   line
*/
#define HAS_ERF HAS_ERF

/* if your compiler does not have a log1p function
   remove the next line.
*/
#define HAS_LOG1P HAS_LOG1P

/* if your compiler does not have the trunc function, remove
 * the next line.
 */

#define HAS_TRUNC HAS_TRUNC

/* Assuming we are using the GNU glibc timing functions */
/* otherwise remove the next line. */
#define CALC_TIME CALC_TIME

#ifdef CALC_TIME
#include <time.h>
#endif 

#ifndef HAS_LOG1P
#define LOG1P(X) log1px(X)
#else
#define LOG1P(X) log1p(X)
#endif

#ifdef HAS_TRUNC
/* Documentation says that trunc() is defined in math.h (included above),
 * but in reality it is often not there.  So I'm definining it here to
 * make a compiler warning shut up. */
extern double trunc(double x);
#define TRUNC(X) trunc(X)
#else
#define TRUNC(X) ((X) != 0.0e0 ? (((X)/fabs(X))*floor(fabs(X))) : 0.0e0)
#endif
		
#define RELDIF(X,Y) ((X) != 0.0e0 ? fabs((X-Y)/(X)): fabs(X-Y))
/* these are constants relating to double precision arithmetic */
#define LNNORM_MAX_X 38.0
#define LNNORM_MIN_X -1.00e9
#define LNANORM_MIN_X -37.519e0

extern double lnnorm(double x);
extern double norm(double x);
extern double anorm(double x);
extern double lnanorm(double x);

#ifndef HAS_LOG1P
extern double log1px(double x);
#endif

#ifdef HAS_ERF
extern double enorm(double x);
extern double lnenorm(double x);
#endif

/* Some more constants.  If upping precision from double, these next
three will need to be made more precise. */
static double sqrpi = 3.989422804014326779399461e-1;
static double sqr2pi = 2.506628274631000502415765e0;
static double sqrt2 = 1.4142135623730950488e0;
/* these next two are used for comparison, not calculation, they do not 
   need to be as precise as the above */
static double thrsh = .67448975e0; 
static double root32 =   5.65685424949238058e0; 
