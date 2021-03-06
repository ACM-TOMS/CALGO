#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_tbsv_c.h
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

  if (N == 0)
    return;

  /* form  x := inv( A )*x */

  if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasUpper)
      || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasLower)) {

    INDEX ix = OFFSET(N, incX) + incX * (N - 1);

    for (i = N; i > 0 && i--;) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      const INDEX j_min = i + 1;
      const INDEX j_max = GSL_MIN(N, i + K + 1);
      INDEX jx = OFFSET(N, incX) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * i + (j - i));
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * i + (j - i));
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }

      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + 0);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + 0);
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

    INDEX ix = OFFSET(N, incX);

    for (i = 0; i < N; i++) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      const INDEX j_min = (K > i ? 0 : i - K);
      const INDEX j_max = i;
      INDEX jx = OFFSET(N, incX) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, lda * i + (K + j - i));
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, lda * i + (K + j - i));
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }
      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, lda * i + K);
        const MpIeee a_imag=  conj * CONST_IMAG(A, lda * i + K);
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

    INDEX ix = OFFSET(N, incX);

    for (i = 0; i < N; i++) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      const INDEX j_min = (K > i ? 0 : i - K);
      const INDEX j_max = i;
      INDEX jx = OFFSET(N, incX) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, (i - j) + lda * j);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, (i - j) + lda * j);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }
      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, 0 + lda * i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, 0 + lda * i);
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

    INDEX ix = OFFSET(N, incX) + incX * (N - 1);

    for (i = N; i > 0 && i--;) {
      MpIeee tmp_real=  REAL(X, ix);
      MpIeee tmp_imag=  IMAG(X, ix);
      const INDEX j_min = i + 1;
      const INDEX j_max = GSL_MIN(N, i + K + 1);
      INDEX jx = OFFSET(N, incX) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const MpIeee Aij_real=  CONST_REAL(A, (K + i - j) + lda * j);
        const MpIeee Aij_imag=  conj * CONST_IMAG(A, (K + i - j) + lda * j);
        const MpIeee x_real=  REAL(X, jx);
        const MpIeee x_imag=  IMAG(X, jx);
        tmp_real -= Aij_real * x_real - Aij_imag * x_imag;
        tmp_imag -= Aij_real * x_imag + Aij_imag * x_real;
        jx += incX;
      }

      if (nonunit) {
        const MpIeee a_real=  CONST_REAL(A, K + lda * i);
        const MpIeee a_imag=  conj * CONST_IMAG(A, K + lda * i);
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
