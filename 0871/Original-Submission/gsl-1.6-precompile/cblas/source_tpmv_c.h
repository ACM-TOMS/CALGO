#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_tpmv_c.h
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

  const int conj = (TransA == CblasConjTrans) ? -1 : 1;
  const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;
  const int nonunit = (Diag == CblasNonUnit);

  if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasUpper)
      || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasLower)) {
    /* form  x:= A*x */

    INDEX ix = OFFSET(N, incX);
    for (i = 0; i < N; i++) {
      const MpIeee Aii_real=  CONST_REAL(Ap, TPUP(N, i, i));
      const MpIeee Aii_imag=  conj * CONST_IMAG(Ap, TPUP(N, i, i));
      MpIeee temp_r;
      MpIeee temp_i;
      if (nonunit) {
        MpIeee x_real=  REAL(X, ix);
        MpIeee x_imag=  IMAG(X, ix);
        temp_r = Aii_real * x_real - Aii_imag * x_imag;
        temp_i = Aii_real * x_imag + Aii_imag * x_real;
      } else {
        temp_r = REAL(X, ix);
        temp_i = IMAG(X, ix);
      }

      {
        INDEX jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          const MpIeee Aij_real=  CONST_REAL(Ap, TPUP(N, i, j));
          const MpIeee Aij_imag=  conj * CONST_IMAG(Ap, TPUP(N, i, j));
          MpIeee x_real=  REAL(X, jx);
          MpIeee x_imag=  IMAG(X, jx);
          temp_r += Aij_real * x_real - Aij_imag * x_imag;
          temp_i += Aij_real * x_imag + Aij_imag * x_real;
          jx += incX;
        }
      }

      REAL(X, ix) = temp_r;
      IMAG(X, ix) = temp_i;
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasUpper)) {

    INDEX ix = OFFSET(N, incX) + incX * (N - 1);
    for (i = N; i > 0 && i--;) {
      const MpIeee Aii_real=  CONST_REAL(Ap, TPLO(N, i, i));
      const MpIeee Aii_imag=  conj * CONST_IMAG(Ap, TPLO(N, i, i));
      MpIeee temp_r;
      MpIeee temp_i;
      if (nonunit) {
        MpIeee x_real=  REAL(X, ix);
        MpIeee x_imag=  IMAG(X, ix);
        temp_r = Aii_real * x_real - Aii_imag * x_imag;
        temp_i = Aii_real * x_imag + Aii_imag * x_real;
      } else {
        temp_r = REAL(X, ix);
        temp_i = IMAG(X, ix);
      }

      {
        INDEX jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          const MpIeee Aij_real=  CONST_REAL(Ap, TPLO(N, i, j));
          const MpIeee Aij_imag=  conj * CONST_IMAG(Ap, TPLO(N, i, j));
          MpIeee x_real=  REAL(X, jx);
          MpIeee x_imag=  IMAG(X, jx);
          temp_r += Aij_real * x_real - Aij_imag * x_imag;
          temp_i += Aij_real * x_imag + Aij_imag * x_real;
          jx += incX;
        }
      }

      REAL(X, ix) = temp_r;
      IMAG(X, ix) = temp_i;
      ix -= incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasUpper)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasLower)) {
    /* form  x := A'*x */

    INDEX ix = OFFSET(N, incX) + incX * (N - 1);
    for (i = N; i > 0 && i--;) {
      const MpIeee Aii_real=  CONST_REAL(Ap, TPUP(N, i, i));
      const MpIeee Aii_imag=  conj * CONST_IMAG(Ap, TPUP(N, i, i));
      MpIeee temp_r;
      MpIeee temp_i;
      if (nonunit) {
        MpIeee x_real=  REAL(X, ix);
        MpIeee x_imag=  IMAG(X, ix);
        temp_r = Aii_real * x_real - Aii_imag * x_imag;
        temp_i = Aii_real * x_imag + Aii_imag * x_real;
      } else {
        temp_r = REAL(X, ix);
        temp_i = IMAG(X, ix);
      }
      {
        INDEX jx = OFFSET(N, incX);
        for (j = 0; j < i; j++) {
          MpIeee x_real=  REAL(X, jx);
          MpIeee x_imag=  IMAG(X, jx);
          const MpIeee Aji_real=  CONST_REAL(Ap, TPUP(N, j, i));
          const MpIeee Aji_imag=  conj * CONST_IMAG(Ap, TPUP(N, j, i));
          temp_r += Aji_real * x_real - Aji_imag * x_imag;
          temp_i += Aji_real * x_imag + Aji_imag * x_real;
          jx += incX;
        }
      }

      REAL(X, ix) = temp_r;
      IMAG(X, ix) = temp_i;
      ix -= incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasUpper)) {

    INDEX ix = OFFSET(N, incX);
    for (i = 0; i < N; i++) {
      const MpIeee Aii_real=  CONST_REAL(Ap, TPLO(N, i, i));
      const MpIeee Aii_imag=  conj * CONST_IMAG(Ap, TPLO(N, i, i));
      MpIeee temp_r;
      MpIeee temp_i;
      if (nonunit) {
        MpIeee x_real=  REAL(X, ix);
        MpIeee x_imag=  IMAG(X, ix);
        temp_r = Aii_real * x_real - Aii_imag * x_imag;
        temp_i = Aii_real * x_imag + Aii_imag * x_real;
      } else {
        temp_r = REAL(X, ix);
        temp_i = IMAG(X, ix);
      }
      {
        INDEX jx = OFFSET(N, incX) + (i + 1) * incX;
        for (j = i + 1; j < N; j++) {
          MpIeee x_real=  REAL(X, jx);
          MpIeee x_imag=  IMAG(X, jx);
          const MpIeee Aji_real=  CONST_REAL(Ap, TPLO(N, j, i));
          const MpIeee Aji_imag=  conj * CONST_IMAG(Ap, TPLO(N, j, i));
          temp_r += Aji_real * x_real - Aji_imag * x_imag;
          temp_i += Aji_real * x_imag + Aji_imag * x_real;
          jx += incX;
        }
      }
      REAL(X, ix) = temp_r;
      IMAG(X, ix) = temp_i;
      ix += incX;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
