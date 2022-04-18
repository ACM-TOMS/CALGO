#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_gemv_c.h
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
  INDEX lenX, lenY;

  const MpIeee alpha_real=  CONST_REAL0(alpha);
  const MpIeee alpha_imag=  CONST_IMAG0(alpha);

  const MpIeee beta_real=  CONST_REAL0(beta);
  const MpIeee beta_imag=  CONST_IMAG0(beta);

  if (M == 0 || N == 0)
    return;

  if ((alpha_real == 0.0 && alpha_imag == 0.0)
      && (beta_real == 1.0 && beta_imag == 0.0))
    return;

  if (TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }

  /* form  y := beta*y */

  if (beta_real == 0.0 && beta_imag == 0.0) {
    INDEX iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      REAL(Y, iy) = 0.0;
      IMAG(Y, iy) = 0.0;
      iy += incY;
    }
  } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
    INDEX iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      const MpIeee y_real=  REAL(Y, iy);
      const MpIeee y_imag=  IMAG(Y, iy);
      const MpIeee tmpR=  y_real * beta_real - y_imag * beta_imag;
      const MpIeee tmpI=  y_real * beta_imag + y_imag * beta_real;
      REAL(Y, iy) = tmpR;
      IMAG(Y, iy) = tmpI;
      iy += incY;
    }
  }

  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;

  if ((order == CblasRowMajor && TransA == CblasNoTrans)
      || (order == CblasColMajor && TransA == CblasTrans)) {
    /* form  y := alpha*A*x + y */
    INDEX iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      MpIeee dotR=  0.0;
      MpIeee dotI=  0.0;
      INDEX ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        const MpIeee x_real=  CONST_REAL(X, ix);
        const MpIeee x_imag=  CONST_IMAG(X, ix);
        const MpIeee A_real=  CONST_REAL(A, lda * i + j);
        const MpIeee A_imag=  CONST_IMAG(A, lda * i + j);

        dotR += A_real * x_real - A_imag * x_imag;
        dotI += A_real * x_imag + A_imag * x_real;
        ix += incX;
      }

      REAL(Y, iy) += alpha_real * dotR - alpha_imag * dotI;
      IMAG(Y, iy) += alpha_real * dotI + alpha_imag * dotR;
      iy += incY;
    }
  } else if ((order == CblasRowMajor && TransA == CblasTrans)
             || (order == CblasColMajor && TransA == CblasNoTrans)) {
    /* form  y := alpha*A'*x + y */
    INDEX ix = OFFSET(lenX, incX);
    for (j = 0; j < lenX; j++) {
      MpIeee x_real=  CONST_REAL(X, ix);
      MpIeee x_imag=  CONST_IMAG(X, ix);
      MpIeee tmpR=  alpha_real * x_real - alpha_imag * x_imag;
      MpIeee tmpI=  alpha_real * x_imag + alpha_imag * x_real;

      INDEX iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        const MpIeee A_real=  CONST_REAL(A, lda * j + i);
        const MpIeee A_imag=  CONST_IMAG(A, lda * j + i);
        REAL(Y, iy) += A_real * tmpR - A_imag * tmpI;
        IMAG(Y, iy) += A_real * tmpI + A_imag * tmpR;
        iy += incY;
      }
      ix += incX;
    }
  } else if (order == CblasRowMajor && TransA == CblasConjTrans) {
    /* form  y := alpha*A^H*x + y */
    INDEX ix = OFFSET(lenX, incX);
    for (j = 0; j < lenX; j++) {
      MpIeee x_real=  CONST_REAL(X, ix);
      MpIeee x_imag=  CONST_IMAG(X, ix);
      MpIeee tmpR=  alpha_real * x_real - alpha_imag * x_imag;
      MpIeee tmpI=  alpha_real * x_imag + alpha_imag * x_real;

      INDEX iy = OFFSET(lenY, incY);
      for (i = 0; i < lenY; i++) {
        const MpIeee A_real=  CONST_REAL(A, lda * j + i);
        const MpIeee A_imag=  CONST_IMAG(A, lda * j + i);
        REAL(Y, iy) += A_real * tmpR - (-A_imag) * tmpI;
        IMAG(Y, iy) += A_real * tmpI + (-A_imag) * tmpR;
        iy += incY;
      }
      ix += incX;
    }
  } else if (order == CblasColMajor && TransA == CblasConjTrans) {
    /* form  y := alpha*A^H*x + y */
    INDEX iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      MpIeee dotR=  0.0;
      MpIeee dotI=  0.0;
      INDEX ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        const MpIeee x_real=  CONST_REAL(X, ix);
        const MpIeee x_imag=  CONST_IMAG(X, ix);
        const MpIeee A_real=  CONST_REAL(A, lda * i + j);
        const MpIeee A_imag=  CONST_IMAG(A, lda * i + j);

        dotR += A_real * x_real - (-A_imag) * x_imag;
        dotI += A_real * x_imag + (-A_imag) * x_real;
        ix += incX;
      }

      REAL(Y, iy) += alpha_real * dotR - alpha_imag * dotI;
      IMAG(Y, iy) += alpha_real * dotI + alpha_imag * dotR;
      iy += incY;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }

}
