#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_trsv_c.h
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
  const int conj = (TransA == CblasConjTrans) ? -1 : 1;
  const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;
  const int nonunit = (Diag == CblasNonUnit);
  INDEX i, j;
  INDEX ix, jx;

  if (N == 0)
    return;

  /* form  x := inv( A )*x */

  if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasUpper)
      || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasLower)) {

    ix = OFFSET(N, incX) + incX * (N - 1);

    if (nonunit) {
      const MpIeee a_real=  CONST_REAL(A, lda * (N - 1) + (N - 1));
      const MpIeee a_imag=  conj * CONST_IMAG(A, lda * (N - 1) + (N - 1));
      const MpIeee x_real=  REAL(X, ix);
      const MpIeee x_imag=  IMAG(X, ix);
      const MpIeee s=  xhypot(a_real, a_imag);
      const MpIeee b_real=  a_real / s;
      const MpIeee b_imag=  a_imag / s;
      REAL(X, ix) = (x_real * b_real + x_imag * b_imag) / s;
      IMAG(X, ix) = (x_imag * b_real - b_imag * x_real) / s;
    }

    ix -= incX;

    for (i = N - 1; i > 0 && i--;) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * i + j);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * i + j);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }

      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(a_real, a_imag);
        const MpIeee b_real=  a_real / s;
        const MpIeee b_imag=  a_imag / s;
        REAL(X, ix) = (tmp_real * b_real + tmp_imag * b_imag) / s;
        IMAG(X, ix) = (tmp_imag * b_real - tmp_real * b_imag) / s;
      } else {
        REAL(X, ix) = tmp_real;
        IMAG(X, ix) = tmp_imag;
      }
      ix -= incX;
    }

  } else if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasUpper)) {
    /* forward substitution */

    ix = OFFSET(N, incX);

    if (nonunit) {
      const MpIeee a_real=  CONST_REAL(A, lda * 0 + 0);
      const MpIeee a_imag=  conj * CONST_IMAG(A, lda * 0 + 0);
      const MpIeee x_real=  REAL(X, ix);
      const MpIeee x_imag=  IMAG(X, ix);
      const MpIeee s=  xhypot(a_real, a_imag);
      const MpIeee b_real=  a_real / s;
      const MpIeee b_imag=  a_imag / s;
      REAL(X, ix) = (x_real * b_real + x_imag * b_imag) / s;
      IMAG(X, ix) = (x_imag * b_real - b_imag * x_real) / s;
    }

    ix += incX;

    for (i = 1; i < N; i++) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      jx = OFFSET(N, incX);
      for (j = 0; j < i; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * i + j);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * i + j);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }
      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(a_real, a_imag);
        const MpIeee b_real=  a_real / s;
        const MpIeee b_imag=  a_imag / s;
        REAL(X, ix) = (tmp_real * b_real + tmp_imag * b_imag) / s;
        IMAG(X, ix) = (tmp_imag * b_real - tmp_real * b_imag) / s;
      } else {
        REAL(X, ix) = tmp_real;
        IMAG(X, ix) = tmp_imag;
      }
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasUpper)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasLower)) {
    /* form  x := inv( A' )*x */

    /* forward substitution */

    ix = OFFSET(N, incX);

    if (nonunit) {
      const MpIeee a_real=  CONST_REAL(A, lda * 0 + 0);
      const MpIeee a_imag=  conj * CONST_IMAG(A, lda * 0 + 0);
      const MpIeee x_real=  REAL(X, ix);
      const MpIeee x_imag=  IMAG(X, ix);
      const MpIeee s=  xhypot(a_real, a_imag);
      const MpIeee b_real=  a_real / s;
      const MpIeee b_imag=  a_imag / s;
      REAL(X, ix) = (x_real * b_real + x_imag * b_imag) / s;
      IMAG(X, ix) = (x_imag * b_real - b_imag * x_real) / s;
    }

    ix += incX;

    for (i = 1; i < N; i++) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      jx = OFFSET(N, incX);
      for (j = 0; j < i; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * j + i);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * j + i);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }
      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(a_real, a_imag);
        const MpIeee b_real=  a_real / s;
        const MpIeee b_imag=  a_imag / s;
        REAL(X, ix) = (tmp_real * b_real + tmp_imag * b_imag) / s;
        IMAG(X, ix) = (tmp_imag * b_real - tmp_real * b_imag) / s;
      } else {
        REAL(X, ix) = tmp_real;
        IMAG(X, ix) = tmp_imag;
      }
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasUpper)) {

    /* backsubstitution */

    ix = OFFSET(N, incX) + incX * (N - 1);

    if (nonunit) {
      const MpIeee a_real=  CONST_REAL(A, lda * (N - 1) + (N - 1));
      const MpIeee a_imag=  conj * CONST_IMAG(A, lda * (N - 1) + (N - 1));
      const MpIeee x_real=  REAL(X, ix);
      const MpIeee x_imag=  IMAG(X, ix);
      const MpIeee s=  xhypot(a_real, a_imag);
      const MpIeee b_real=  a_real / s;
      const MpIeee b_imag=  a_imag / s;
      REAL(X, ix) = (x_real * b_real + x_imag * b_imag) / s;
      IMAG(X, ix) = (x_imag * b_real - b_imag * x_real) / s;
    }

    ix -= incX;

    for (i = N - 1; i > 0 && i--;) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * j + i);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * j + i);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }

      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(a_real, a_imag);
        const MpIeee b_real=  a_real / s;
        const MpIeee b_imag=  a_imag / s;
        REAL(X, ix) = (tmp_real * b_real + tmp_imag * b_imag) / s;
        IMAG(X, ix) = (tmp_imag * b_real - tmp_real * b_imag) / s;
      } else {
        REAL(X, ix) = tmp_real;
        IMAG(X, ix) = tmp_imag;
      }
      ix -= incX;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
