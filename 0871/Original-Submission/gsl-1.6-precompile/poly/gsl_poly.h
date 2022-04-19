#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* poly/gsl_poly.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Brian Gough
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

#ifndef __GSL_POLY_H__
#define __GSL_POLY_H__

#include <stdlib.h>
#include <gsl/gsl_complex.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Evaluate polynomial
 *
 * c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
 *
 * exceptions: none
 */
MpIeee gsl_poly_eval(const MpIeee c[], const int len, const MpIeee x);


#ifdef HAVE_INLINE
extern inline
MpIeee gsl_poly_eval(const MpIeee c[], const int len, const MpIeee x)
{
  int  i;
  MpIeee ans=  c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}
#endif /* HAVE_INLINE */

/* Work with divided-difference polynomials, Abramowitz & Stegun 25.2.26 */

int
 gsl_poly_dd_init(MpIeee dd[], const MpIeee x[], const MpIeee y[],
                  size_t size);

MpIeee gsl_poly_dd_eval(const MpIeee dd[], const MpIeee xa[], const size_t size, const MpIeee x);

#ifdef HAVE_INLINE
extern inline
MpIeee gsl_poly_dd_eval(const MpIeee dd[], const MpIeee xa[], const size_t size, const MpIeee x)
{
  size_t i;
  MpIeee y=  dd[size - 1];
  for (i = size - 1; i--;) y = dd[i] + (x - xa[i]) * y;
  return y;
}
#endif /* HAVE_INLINE */


int
 gsl_poly_dd_taylor(MpIeee c[], MpIeee xp,
                    const MpIeee dd[], const MpIeee x[], size_t size,
                    MpIeee w[]);

/* Solve for real or complex roots of the standard quadratic equation,
 * returning the number of real roots.
 *
 * Roots are returned ordered.
 */
int  gsl_poly_solve_quadratic(MpIeee a, MpIeee b, MpIeee c, 
                              MpIeee * x0, MpIeee * x1);

int 
 gsl_poly_complex_solve_quadratic(MpIeee a, MpIeee b, MpIeee c, 
                                  gsl_complex * z0, gsl_complex * z1);


/* Solve for real roots of the cubic equation
 * x^3 + a x^2 + b x + c = 0, returning the
 * number of real roots.
 *
 * Roots are returned ordered.
 */
int  gsl_poly_solve_cubic(MpIeee a, MpIeee b, MpIeee c, 
                          MpIeee * x0, MpIeee * x1, MpIeee * x2);

int 
 gsl_poly_complex_solve_cubic(MpIeee a, MpIeee b, MpIeee c, 
                              gsl_complex * z0, gsl_complex * z1, 
                              gsl_complex * z2);


/* Solve for the complex roots of a general real polynomial */

typedef struct 
{ 
  size_t nc ;
  MpIeee * matrix; 
} 
gsl_poly_complex_workspace ;

gsl_poly_complex_workspace * gsl_poly_complex_workspace_alloc (size_t n);
void gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w);

int
 gsl_poly_complex_solve(const MpIeee * a, size_t n, 
                        gsl_poly_complex_workspace * w,
                        gsl_complex_packed_ptr z);

__END_DECLS

#endif /* __GSL_POLY_H__ */
