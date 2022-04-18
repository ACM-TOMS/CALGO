#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_syr2k_c.h
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
  int uplo, trans;

  const MpIeee alpha_real=  CONST_REAL0(alpha);
  const MpIeee alpha_imag=  CONST_IMAG0(alpha);
  const MpIeee beta_real=  CONST_REAL0(beta);
  const MpIeee beta_imag=  CONST_IMAG0(beta);

  if ((alpha_real == 0.0 && alpha_imag == 0.0)
      && (beta_real == 1.0 && beta_imag == 0.0))
    return;

  if (Order == CblasRowMajor) {
    uplo = Uplo;
    trans = Trans;
  } else {
    uplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;
    trans = (Trans == CblasNoTrans) ? CblasTrans : CblasNoTrans;
  }

  /* form  C := beta*C */

  if (beta_real == 0.0 && beta_imag == 0.0) {
    if (uplo == CblasUpper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          REAL(C, ldc * i + j) = 0.0;
          IMAG(C, ldc * i + j) = 0.0;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          REAL(C, ldc * i + j) = 0.0;
          IMAG(C, ldc * i + j) = 0.0;
        }
      }
    }
  } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
    if (uplo == CblasUpper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          const MpIeee Cij_real=  REAL(C, ldc * i + j);
          const MpIeee Cij_imag=  IMAG(C, ldc * i + j);
          REAL(C, ldc * i + j) = beta_real * Cij_real - beta_imag * Cij_imag;
          IMAG(C, ldc * i + j) = beta_real * Cij_imag + beta_imag * Cij_real;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          const MpIeee Cij_real=  REAL(C, ldc * i + j);
          const MpIeee Cij_imag=  IMAG(C, ldc * i + j);
          REAL(C, ldc * i + j) = beta_real * Cij_real - beta_imag * Cij_imag;
          IMAG(C, ldc * i + j) = beta_real * Cij_imag + beta_imag * Cij_real;
        }
      }
    }
  }


  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;

  if (uplo == CblasUpper && trans == CblasNoTrans) {

    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);
          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  CONST_IMAG(A, j * lda + k);
          const MpIeee Bjk_real=  CONST_REAL(B, j * ldb + k);
          const MpIeee Bjk_imag=  CONST_IMAG(B, j * ldb + k);
          temp_real += ((Aik_real * Bjk_real - Aik_imag * Bjk_imag)
                        + (Bik_real * Ajk_real - Bik_imag * Ajk_imag));
          temp_imag += ((Aik_real * Bjk_imag + Aik_imag * Bjk_real)
                        + (Bik_real * Ajk_imag + Bik_imag * Ajk_real));
        }
        REAL(C, i * ldc + j) += alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(C, i * ldc + j) += alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (uplo == CblasUpper && trans == CblasTrans) {

    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        MpIeee Aki_imag=  CONST_IMAG(A, k * lda + i);
        MpIeee Bki_real=  CONST_REAL(B, k * ldb + i);
        MpIeee Bki_imag=  CONST_IMAG(B, k * ldb + i);
        MpIeee temp1_real=  alpha_real * Aki_real - alpha_imag * Aki_imag;
        MpIeee temp1_imag=  alpha_real * Aki_imag + alpha_imag * Aki_real;
        MpIeee temp2_real=  alpha_real * Bki_real - alpha_imag * Bki_imag;
        MpIeee temp2_imag=  alpha_real * Bki_imag + alpha_imag * Bki_real;
        for (j = i; j < N; j++) {
          MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          MpIeee Akj_imag=  CONST_IMAG(A, k * lda + j);
          MpIeee Bkj_real=  CONST_REAL(B, k * ldb + j);
          MpIeee Bkj_imag=  CONST_IMAG(B, k * ldb + j);
          REAL(C, i * lda + j) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          IMAG(C, i * lda + j) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }
      }
    }

  } else if (uplo == CblasLower && trans == CblasNoTrans) {


    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);
          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  CONST_IMAG(A, j * lda + k);
          const MpIeee Bjk_real=  CONST_REAL(B, j * ldb + k);
          const MpIeee Bjk_imag=  CONST_IMAG(B, j * ldb + k);
          temp_real += ((Aik_real * Bjk_real - Aik_imag * Bjk_imag)
                        + (Bik_real * Ajk_real - Bik_imag * Ajk_imag));
          temp_imag += ((Aik_real * Bjk_imag + Aik_imag * Bjk_real)
                        + (Bik_real * Ajk_imag + Bik_imag * Ajk_real));
        }
        REAL(C, i * ldc + j) += alpha_real * temp_real - alpha_imag * temp_imag;
        IMAG(C, i * ldc + j) += alpha_real * temp_imag + alpha_imag * temp_real;
      }
    }

  } else if (uplo == CblasLower && trans == CblasTrans) {

    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        MpIeee Aki_imag=  CONST_IMAG(A, k * lda + i);
        MpIeee Bki_real=  CONST_REAL(B, k * ldb + i);
        MpIeee Bki_imag=  CONST_IMAG(B, k * ldb + i);
        MpIeee temp1_real=  alpha_real * Aki_real - alpha_imag * Aki_imag;
        MpIeee temp1_imag=  alpha_real * Aki_imag + alpha_imag * Aki_real;
        MpIeee temp2_real=  alpha_real * Bki_real - alpha_imag * Bki_imag;
        MpIeee temp2_imag=  alpha_real * Bki_imag + alpha_imag * Bki_real;
        for (j = 0; j <= i; j++) {
          MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          MpIeee Akj_imag=  CONST_IMAG(A, k * lda + j);
          MpIeee Bkj_real=  CONST_REAL(B, k * ldb + j);
          MpIeee Bkj_imag=  CONST_IMAG(B, k * ldb + j);
          REAL(C, i * lda + j) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          IMAG(C, i * lda + j) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }
      }
    }


  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
