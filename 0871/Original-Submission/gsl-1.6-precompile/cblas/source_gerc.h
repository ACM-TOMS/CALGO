#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_gerc.h
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
  INDEX i, j;

  const MpIeee alpha_real=  CONST_REAL0(alpha);
  const MpIeee alpha_imag=  CONST_IMAG0(alpha);

  if (order == CblasRowMajor) {
    INDEX ix = OFFSET(M, incX);
    for (i = 0; i < M; i++) {
      const MpIeee X_real=  CONST_REAL(X, ix);
      const MpIeee X_imag=  CONST_IMAG(X, ix);
      const MpIeee tmp_real=  alpha_real * X_real - alpha_imag * X_imag;
      const MpIeee tmp_imag=  alpha_imag * X_real + alpha_real * X_imag;
      INDEX jy = OFFSET(N, incY);
      for (j = 0; j < N; j++) {
        const MpIeee Y_real=  CONST_REAL(Y, jy);
        const MpIeee Y_imag=  -CONST_IMAG(Y, jy);
        REAL(A, lda * i + j) += Y_real * tmp_real - Y_imag * tmp_imag;
        IMAG(A, lda * i + j) += Y_imag * tmp_real + Y_real * tmp_imag;
        jy += incY;
      }
      ix += incX;
    }
  } else if (order == CblasColMajor) {
    INDEX jy = OFFSET(N, incY);
    for (j = 0; j < N; j++) {
      const MpIeee Y_real=  CONST_REAL(Y, jy);
      const MpIeee Y_imag=  -CONST_IMAG(Y, jy);
      const MpIeee tmp_real=  alpha_real * Y_real - alpha_imag * Y_imag;
      const MpIeee tmp_imag=  alpha_imag * Y_real + alpha_real * Y_imag;
      INDEX ix = OFFSET(M, incX);
      for (i = 0; i < M; i++) {
        const MpIeee X_real=  CONST_REAL(X, ix);
        const MpIeee X_imag=  CONST_IMAG(X, ix);
        REAL(A, i + lda * j) += X_real * tmp_real - X_imag * tmp_imag;
        IMAG(A, i + lda * j) += X_imag * tmp_real + X_real * tmp_imag;
        ix += incX;
      }
      jy += incY;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
