#ifndef KERNEL_H

#define KERNEL_H

/* several macro variables */
#define ITMAX 1000000000 /* Maximum allowed number of iterations */
#define DPMIN 1.0e-300 /* Number near the smallest representable double-point number */
#define EPS 2.22e-16 /* Machine epsilon */
#define NITERMAX_ROMBERG 15 /* Maximum allowed number of Romberg iterations */
#define TOL_ROMBERG 0.1 /* Tolerance factor used to stop the Romberg iterations */
#define TOL_DIFF 0.2 /* Tolerance factor used for the approximation of I_{x,y}^{mu,p} using differences */

#endif
