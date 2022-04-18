#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_trsm_c.h
 * 
 * Copyright (C) 2001 Brian Gough
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
  INDEX i, j, k;
  INDEX n1, n2;
  const int nonunit = (Diag == CblasNonUnit);
  const int conj = (TransA == CblasConjTrans) ? -1 : 1;
  int side, uplo, trans;

  const MpIeee alpha_real=  CONST_REAL0(alpha);
  const MpIeee alpha_imag=  CONST_IMAG0(alpha);

  if (Order == CblasRowMajor) {
    n1 = M;
    n2 = N;
    side = Side;
    uplo = Uplo;
    trans = TransA;
    trans = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
  } else {
    n1 = N;
    n2 = M;
    side = (Side == CblasLeft) ? CblasRight : CblasLeft;        /* exchanged */
    uplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;      /* exchanged */
    trans = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;       /* same */
  }

  if (side == CblasLeft && uplo == CblasUpper && trans == CblasNoTrans) {

    /* form  B := alpha * inv(TriU(A)) *B */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = n1; i > 0 && i--;) {
      if (nonunit) {
        const MpIeee Aii_real=  CONST_REAL(A, lda * i + i);
        const MpIeee Aii_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(Aii_real, Aii_imag);
        const MpIeee a_real=  Aii_real / s;
        const MpIeee a_imag=  Aii_imag / s;

        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }
      }

      for (k = 0; k < i; k++) {
        const MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        const MpIeee Aki_imag=  conj * CONST_IMAG(A, k * lda + i);
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * k + j) -= Aki_real * Bij_real - Aki_imag * Bij_imag;
          IMAG(B, ldb * k + j) -= Aki_real * Bij_imag + Aki_imag * Bij_real;
        }
      }
    }

  } else if (side == CblasLeft && uplo == CblasUpper && trans == CblasTrans) {

    /* form  B := alpha * inv(TriU(A))' *B */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {

      if (nonunit) {
        const MpIeee Aii_real=  CONST_REAL(A, lda * i + i);
        const MpIeee Aii_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(Aii_real, Aii_imag);
        const MpIeee a_real=  Aii_real / s;
        const MpIeee a_imag=  Aii_imag / s;

        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }
      }

      for (k = i + 1; k < n1; k++) {
        const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
        const MpIeee Aik_imag=  conj * CONST_IMAG(A, i * lda + k);
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * k + j) -= Aik_real * Bij_real - Aik_imag * Bij_imag;
          IMAG(B, ldb * k + j) -= Aik_real * Bij_imag + Aik_imag * Bij_real;
        }
      }
    }

  } else if (side == CblasLeft && uplo == CblasLower && trans == CblasNoTrans) {

    /* form  B := alpha * inv(TriL(A))*B */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {

      if (nonunit) {
        const MpIeee Aii_real=  CONST_REAL(A, lda * i + i);
        const MpIeee Aii_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(Aii_real, Aii_imag);
        const MpIeee a_real=  Aii_real / s;
        const MpIeee a_imag=  Aii_imag / s;

        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }
      }

      for (k = i + 1; k < n1; k++) {
        const MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        const MpIeee Aki_imag=  conj * CONST_IMAG(A, k * lda + i);
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * k + j) -= Aki_real * Bij_real - Aki_imag * Bij_imag;
          IMAG(B, ldb * k + j) -= Aki_real * Bij_imag + Aki_imag * Bij_real;
        }
      }
    }


  } else if (side == CblasLeft && uplo == CblasLower && trans == CblasTrans) {

    /* form  B := alpha * TriL(A)' *B */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = n1; i > 0 && i--;) {
      if (nonunit) {
        const MpIeee Aii_real=  CONST_REAL(A, lda * i + i);
        const MpIeee Aii_imag=  conj * CONST_IMAG(A, lda * i + i);
        const MpIeee s=  xhypot(Aii_real, Aii_imag);
        const MpIeee a_real=  Aii_real / s;
        const MpIeee a_imag=  Aii_imag / s;

        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }
      }

      for (k = 0; k < i; k++) {
        const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
        const MpIeee Aik_imag=  conj * CONST_IMAG(A, i * lda + k);
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * k + j) -= Aik_real * Bij_real - Aik_imag * Bij_imag;
          IMAG(B, ldb * k + j) -= Aik_real * Bij_imag + Aik_imag * Bij_real;
        }
      }
    }

  } else if (side == CblasRight && uplo == CblasUpper && trans == CblasNoTrans) {

    /* form  B := alpha * B * inv(TriU(A)) */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, lda * j + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, lda * j + j);
          const MpIeee s=  xhypot(Ajj_real, Ajj_imag);
          const MpIeee a_real=  Ajj_real / s;
          const MpIeee a_imag=  Ajj_imag / s;
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }

        {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          for (k = j + 1; k < n2; k++) {
            const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
            const MpIeee Ajk_imag=  conj * CONST_IMAG(A, j * lda + k);
            REAL(B, ldb * i + k) -= Ajk_real * Bij_real - Ajk_imag * Bij_imag;
            IMAG(B, ldb * i + k) -= Ajk_real * Bij_imag + Ajk_imag * Bij_real;
          }
        }
      }
    }

  } else if (side == CblasRight && uplo == CblasUpper && trans == CblasTrans) {

    /* form  B := alpha * B * inv(TriU(A))' */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, lda * j + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, lda * j + j);
          const MpIeee s=  xhypot(Ajj_real, Ajj_imag);
          const MpIeee a_real=  Ajj_real / s;
          const MpIeee a_imag=  Ajj_imag / s;
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }

        {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          for (k = 0; k < j; k++) {
            const MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
            const MpIeee Akj_imag=  conj * CONST_IMAG(A, k * lda + j);
            REAL(B, ldb * i + k) -= Akj_real * Bij_real - Akj_imag * Bij_imag;
            IMAG(B, ldb * i + k) -= Akj_real * Bij_imag + Akj_imag * Bij_real;
          }
        }
      }
    }


  } else if (side == CblasRight && uplo == CblasLower && trans == CblasNoTrans) {

    /* form  B := alpha * B * inv(TriL(A)) */

    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, lda * j + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, lda * j + j);
          const MpIeee s=  xhypot(Ajj_real, Ajj_imag);
          const MpIeee a_real=  Ajj_real / s;
          const MpIeee a_imag=  Ajj_imag / s;
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }

        {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          for (k = 0; k < j; k++) {
            const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
            const MpIeee Ajk_imag=  conj * CONST_IMAG(A, j * lda + k);
            REAL(B, ldb * i + k) -= Ajk_real * Bij_real - Ajk_imag * Bij_imag;
            IMAG(B, ldb * i + k) -= Ajk_real * Bij_imag + Ajk_imag * Bij_real;
          }
        }
      }
    }

  } else if (side == CblasRight && uplo == CblasLower && trans == CblasTrans) {

    /* form  B := alpha * B * inv(TriL(A))' */


    if (!(alpha_real == 1.0 && alpha_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = alpha_real * Bij_real - alpha_imag * Bij_imag;
          IMAG(B, ldb * i + j) = alpha_real * Bij_imag + alpha_imag * Bij_real;
        }
      }
    }

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, lda * j + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, lda * j + j);
          const MpIeee s=  xhypot(Ajj_real, Ajj_imag);
          const MpIeee a_real=  Ajj_real / s;
          const MpIeee a_imag=  Ajj_imag / s;
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);
          REAL(B, ldb * i + j) = (Bij_real * a_real + Bij_imag * a_imag) / s;
          IMAG(B, ldb * i + j) = (Bij_imag * a_real - Bij_real * a_imag) / s;
        }

        {
          const MpIeee Bij_real=  REAL(B, ldb * i + j);
          const MpIeee Bij_imag=  IMAG(B, ldb * i + j);

          for (k = j + 1; k < n2; k++) {
            const MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
            const MpIeee Akj_imag=  conj * CONST_IMAG(A, k * lda + j);
            REAL(B, ldb * i + k) -= Akj_real * Bij_real - Akj_imag * Bij_imag;
            IMAG(B, ldb * i + k) -= Akj_real * Bij_imag + Akj_imag * Bij_real;
          }
        }
      }
    }


  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
