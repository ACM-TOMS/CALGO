#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_trmm_c.h
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
    trans = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
  } else {
    n1 = N;
    n2 = M;
    side = (Side == CblasLeft) ? CblasRight : CblasLeft;        /* exchanged */
    uplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;      /* exchanged */
    trans = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;       /* same */
  }

  if (side == CblasLeft && uplo == CblasUpper && trans == CblasNoTrans) {

    /* form  B := alpha * TriU(A)*B */

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        if (nonunit) {
          const MpIeee Aii_real=  CONST_REAL(A, i * lda + i);
          const MpIeee Aii_imag=  conj * CONST_IMAG(A, i * lda + i);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real = Aii_real * Bij_real - Aii_imag * Bij_imag;
          temp_imag = Aii_real * Bij_imag + Aii_imag * Bij_real;
        } else {
          temp_real = REAL(B, i * ldb + j);
          temp_imag = IMAG(B, i * ldb + j);
        }

        for (k = i + 1; k < n1; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  conj * CONST_IMAG(A, i * lda + k);
          const MpIeee Bkj_real=  REAL(B, k * ldb + j);
          const MpIeee Bkj_imag=  IMAG(B, k * ldb + j);
          temp_real += Aik_real * Bkj_real - Aik_imag * Bkj_imag;
          temp_imag += Aik_real * Bkj_imag + Aik_imag * Bkj_real;
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasLeft && uplo == CblasUpper && trans == CblasTrans) {

    /* form  B := alpha * (TriU(A))' *B */

    for (i = n1; i > 0 && i--;) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        for (k = 0; k < i; k++) {
          const MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
          const MpIeee Aki_imag=  conj * CONST_IMAG(A, k * lda + i);
          const MpIeee Bkj_real=  REAL(B, k * ldb + j);
          const MpIeee Bkj_imag=  IMAG(B, k * ldb + j);
          temp_real += Aki_real * Bkj_real - Aki_imag * Bkj_imag;
          temp_imag += Aki_real * Bkj_imag + Aki_imag * Bkj_real;
        }

        if (nonunit) {
          const MpIeee Aii_real=  CONST_REAL(A, i * lda + i);
          const MpIeee Aii_imag=  conj * CONST_IMAG(A, i * lda + i);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real += Aii_real * Bij_real - Aii_imag * Bij_imag;
          temp_imag += Aii_real * Bij_imag + Aii_imag * Bij_real;
        } else {
          temp_real += REAL(B, i * ldb + j);
          temp_imag += IMAG(B, i * ldb + j);
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasLeft && uplo == CblasLower && trans == CblasNoTrans) {

    /* form  B := alpha * TriL(A)*B */


    for (i = n1; i > 0 && i--;) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        for (k = 0; k < i; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  conj * CONST_IMAG(A, i * lda + k);
          const MpIeee Bkj_real=  REAL(B, k * ldb + j);
          const MpIeee Bkj_imag=  IMAG(B, k * ldb + j);
          temp_real += Aik_real * Bkj_real - Aik_imag * Bkj_imag;
          temp_imag += Aik_real * Bkj_imag + Aik_imag * Bkj_real;
        }

        if (nonunit) {
          const MpIeee Aii_real=  CONST_REAL(A, i * lda + i);
          const MpIeee Aii_imag=  conj * CONST_IMAG(A, i * lda + i);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real += Aii_real * Bij_real - Aii_imag * Bij_imag;
          temp_imag += Aii_real * Bij_imag + Aii_imag * Bij_real;
        } else {
          temp_real += REAL(B, i * ldb + j);
          temp_imag += IMAG(B, i * ldb + j);
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }



  } else if (side == CblasLeft && uplo == CblasLower && trans == CblasTrans) {

    /* form  B := alpha * TriL(A)' *B */

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        if (nonunit) {
          const MpIeee Aii_real=  CONST_REAL(A, i * lda + i);
          const MpIeee Aii_imag=  conj * CONST_IMAG(A, i * lda + i);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real = Aii_real * Bij_real - Aii_imag * Bij_imag;
          temp_imag = Aii_real * Bij_imag + Aii_imag * Bij_real;
        } else {
          temp_real = REAL(B, i * ldb + j);
          temp_imag = IMAG(B, i * ldb + j);
        }

        for (k = i + 1; k < n1; k++) {
          const MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
          const MpIeee Aki_imag=  conj * CONST_IMAG(A, k * lda + i);
          const MpIeee Bkj_real=  REAL(B, k * ldb + j);
          const MpIeee Bkj_imag=  IMAG(B, k * ldb + j);
          temp_real += Aki_real * Bkj_real - Aki_imag * Bkj_imag;
          temp_imag += Aki_real * Bkj_imag + Aki_imag * Bkj_real;
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasRight && uplo == CblasUpper && trans == CblasNoTrans) {

    /* form  B := alpha * B * TriU(A) */

    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        for (k = 0; k < j; k++) {
          const MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          const MpIeee Akj_imag=  conj * CONST_IMAG(A, k * lda + j);
          const MpIeee Bik_real=  REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  IMAG(B, i * ldb + k);
          temp_real += Akj_real * Bik_real - Akj_imag * Bik_imag;
          temp_imag += Akj_real * Bik_imag + Akj_imag * Bik_real;
        }

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, j * lda + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, j * lda + j);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real += Ajj_real * Bij_real - Ajj_imag * Bij_imag;
          temp_imag += Ajj_real * Bij_imag + Ajj_imag * Bij_real;
        } else {
          temp_real += REAL(B, i * ldb + j);
          temp_imag += IMAG(B, i * ldb + j);
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasRight && uplo == CblasUpper && trans == CblasTrans) {

    /* form  B := alpha * B * (TriU(A))' */

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, j * lda + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, j * lda + j);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real = Ajj_real * Bij_real - Ajj_imag * Bij_imag;
          temp_imag = Ajj_real * Bij_imag + Ajj_imag * Bij_real;
        } else {
          temp_real = REAL(B, i * ldb + j);
          temp_imag = IMAG(B, i * ldb + j);
        }

        for (k = j + 1; k < n2; k++) {
          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  conj * CONST_IMAG(A, j * lda + k);
          const MpIeee Bik_real=  REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  IMAG(B, i * ldb + k);
          temp_real += Ajk_real * Bik_real - Ajk_imag * Bik_imag;
          temp_imag += Ajk_real * Bik_imag + Ajk_imag * Bik_real;
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasRight && uplo == CblasLower && trans == CblasNoTrans) {

    /* form  B := alpha *B * TriL(A) */

    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, j * lda + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, j * lda + j);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real = Ajj_real * Bij_real - Ajj_imag * Bij_imag;
          temp_imag = Ajj_real * Bij_imag + Ajj_imag * Bij_real;
        } else {
          temp_real = REAL(B, i * ldb + j);
          temp_imag = IMAG(B, i * ldb + j);
        }

        for (k = j + 1; k < n2; k++) {
          const MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          const MpIeee Akj_imag=  conj * CONST_IMAG(A, k * lda + j);
          const MpIeee Bik_real=  REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  IMAG(B, i * ldb + k);
          temp_real += Akj_real * Bik_real - Akj_imag * Bik_imag;
          temp_imag += Akj_real * Bik_imag + Akj_imag * Bik_real;
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (side == CblasRight && uplo == CblasLower && trans == CblasTrans) {

    /* form  B := alpha * B * TriL(A)' */

    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;

        for (k = 0; k < j; k++) {
          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  conj * CONST_IMAG(A, j * lda + k);
          const MpIeee Bik_real=  REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  IMAG(B, i * ldb + k);
          temp_real += Ajk_real * Bik_real - Ajk_imag * Bik_imag;
          temp_imag += Ajk_real * Bik_imag + Ajk_imag * Bik_real;
        }

        if (nonunit) {
          const MpIeee Ajj_real=  CONST_REAL(A, j * lda + j);
          const MpIeee Ajj_imag=  conj * CONST_IMAG(A, j * lda + j);
          const MpIeee Bij_real=  REAL(B, i * ldb + j);
          const MpIeee Bij_imag=  IMAG(B, i * ldb + j);
          temp_real += Ajj_real * Bij_real - Ajj_imag * Bij_imag;
          temp_imag += Ajj_real * Bij_imag + Ajj_imag * Bij_real;
        } else {
          temp_real += REAL(B, i * ldb + j);
          temp_imag += IMAG(B, i * ldb + j);
        }

        REAL(B, ldb * i + j) = alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(B, ldb * i + j) = alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
