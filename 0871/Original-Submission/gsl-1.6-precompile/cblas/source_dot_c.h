#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_dot_c.h
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
  MpIeee r_real=  0.0;
  MpIeee r_imag=  0.0;
  INDEX i;
  INDEX ix = OFFSET(N, incX);
  INDEX iy = OFFSET(N, incY);
  for (i = 0; i < N; i++) {
    const MpIeee x_real=  CONST_REAL(X, ix);
    const MpIeee x_imag=  CONST_IMAG(X, ix);
    const MpIeee y_real=  CONST_REAL(Y, iy);
    const MpIeee y_imag=  CONST_IMAG(Y, iy);
    r_real += x_real * y_real - CONJ_SIGN * x_imag * y_imag;
    r_imag += x_real * y_imag + CONJ_SIGN * x_imag * y_real;
    ix += incX;
    iy += incY;
  }
  REAL0(result) = r_real;
  IMAG0(result) = r_imag;
}
