#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_axpy_c.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

{
  INDEX i;
  INDEX ix = OFFSET(N, incX);
  INDEX iy = OFFSET(N, incY);

  const MpIeee alpha_real=  CONST_REAL0(alpha);
  const MpIeee alpha_imag=  CONST_IMAG0(alpha);

  if (fabs(alpha_real) == 0 && fabs(alpha_imag) == 0) {
    return;
  }

  for (i = 0; i < N; i++) {
    const MpIeee x_real=  CONST_REAL(X, ix);
    const MpIeee x_imag=  CONST_IMAG(X, ix);
    REAL(Y, iy) += (alpha_real * x_real - alpha_imag * x_imag);
    IMAG(Y, iy) += (alpha_real * x_imag + alpha_imag * x_real);
    ix += incX;
    iy += incY;
  }
}
