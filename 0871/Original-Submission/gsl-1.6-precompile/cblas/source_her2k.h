#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_her2k_c.h
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
  MpIeee alpha_imag=  CONST_IMAG0(alpha);

  if (beta == 1.0 && ((alpha_real == 0.0 && alpha_imag == 0.0) || K == 0))
    return;

  if (Order == CblasRowMajor) {
    uplo = Uplo;
    trans = Trans;
  } else {
    uplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;
    trans = (Trans == CblasNoTrans) ? CblasConjTrans : CblasNoTrans;
    alpha_imag *= -1;           /* conjugate alpha */
  }

  /* form  C := beta*C */

  if (beta == 0.0) {
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
  } else if (beta != 1.0) {
    if (uplo == CblasUpper) {
      for (i = 0; i < N; i++) {
        REAL(C, ldc * i + i) *= beta;
        IMAG(C, ldc * i + i) = 0.0;
        for (j = i + 1; j < N; j++) {
          REAL(C, ldc * i + j) *= beta;
          IMAG(C, ldc * i + j) *= beta;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
          REAL(C, ldc * i + j) *= beta;
          IMAG(C, ldc * i + j) *= beta;
        }
        REAL(C, ldc * i + i) *= beta;
        IMAG(C, ldc * i + i) = 0.0;
      }
    }
  } else {
    for (i = 0; i < N; i++) {
      IMAG(C, ldc * i + i) = 0.0;
    }
  }

  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;

  if (uplo == CblasUpper && trans == CblasNoTrans) {

    for (i = 0; i < N; i++) {

      /* Cii += alpha Aik conj(Bik) + conj(alpha) Bik conj(Aik) */
      {
        MpIeee temp_real=  0.0;
        /* BASE temp_imag = 0.0; */
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          /* temp1 = alpha * Aik */
          const MpIeee temp1_real=  alpha_real * Aik_real - alpha_imag * Aik_imag;
          const MpIeee temp1_imag=  alpha_real * Aik_imag + alpha_imag * Aik_real;
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);
          temp_real += temp1_real * Bik_real + temp1_imag * Bik_imag;
        }

        REAL(C, i * ldc + i) += 2 * temp_real;
        IMAG(C, i * ldc + i) = 0.0;
      }

      /* Cij += alpha Aik conj(Bjk) + conj(alpha) Bik conj(Ajk) */
      for (j = i + 1; j < N; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          /* temp1 = alpha * Aik */
          const MpIeee temp1_real=  alpha_real * Aik_real - alpha_imag * Aik_imag;
          const MpIeee temp1_imag=  alpha_real * Aik_imag + alpha_imag * Aik_real;
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);

          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  CONST_IMAG(A, j * lda + k);
          /* temp2 = alpha * Ajk */
          const MpIeee temp2_real=  alpha_real * Ajk_real - alpha_imag * Ajk_imag;
          const MpIeee temp2_imag=  alpha_real * Ajk_imag + alpha_imag * Ajk_real;
          const MpIeee Bjk_real=  CONST_REAL(B, j * ldb + k);
          const MpIeee Bjk_imag=  CONST_IMAG(B, j * ldb + k);

          /* Cij += alpha * Aik * conj(Bjk) + conj(alpha) * Bik * conj(Ajk) */
          temp_real += ((temp1_real * Bjk_real + temp1_imag * Bjk_imag)
                        + (Bik_real * temp2_real + Bik_imag * temp2_imag));
          temp_imag += ((temp1_real * (-Bjk_imag) + temp1_imag * Bjk_real)
                        + (Bik_real * (-temp2_imag) + Bik_imag * temp2_real));
        }
        REAL(C, i * ldc + j) += temp_real;
        IMAG(C, i * ldc + j) += temp_imag;
      }
    }

  } else if (uplo == CblasUpper && trans == CblasConjTrans) {

    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        MpIeee Aki_imag=  CONST_IMAG(A, k * lda + i);
        MpIeee Bki_real=  CONST_REAL(B, k * ldb + i);
        MpIeee Bki_imag=  CONST_IMAG(B, k * ldb + i);
        /* temp1 = alpha * conj(Aki) */
        MpIeee temp1_real=  alpha_real * Aki_real - alpha_imag * (-Aki_imag);
        MpIeee temp1_imag=  alpha_real * (-Aki_imag) + alpha_imag * Aki_real;
        /* temp2 = conj(alpha) * conj(Bki) */
        MpIeee temp2_real=  alpha_real * Bki_real - alpha_imag * Bki_imag;
        MpIeee temp2_imag=  -(alpha_real * Bki_imag + alpha_imag * Bki_real);

        /* Cii += alpha * conj(Aki) * Bki + conj(alpha) * conj(Bki) * Aki */
        {
          REAL(C, i * lda + i) += 2 * (temp1_real * Bki_real - temp1_imag * Bki_imag);
          IMAG(C, i * lda + i) = 0.0;
        }

        for (j = i + 1; j < N; j++) {
          MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          MpIeee Akj_imag=  CONST_IMAG(A, k * lda + j);
          MpIeee Bkj_real=  CONST_REAL(B, k * ldb + j);
          MpIeee Bkj_imag=  CONST_IMAG(B, k * ldb + j);
          /* Cij += alpha * conj(Aki) * Bkj + conj(alpha) * conj(Bki) * Akj */
          REAL(C, i * lda + j) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          IMAG(C, i * lda + j) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }
      }
    }

  } else if (uplo == CblasLower && trans == CblasNoTrans) {

    for (i = 0; i < N; i++) {

      /* Cij += alpha Aik conj(Bjk) + conj(alpha) Bik conj(Ajk) */

      for (j = 0; j < i; j++) {
        MpIeee temp_real=  0.0;
        MpIeee temp_imag=  0.0;
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          /* temp1 = alpha * Aik */
          const MpIeee temp1_real=  alpha_real * Aik_real - alpha_imag * Aik_imag;
          const MpIeee temp1_imag=  alpha_real * Aik_imag + alpha_imag * Aik_real;
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);

          const MpIeee Ajk_real=  CONST_REAL(A, j * lda + k);
          const MpIeee Ajk_imag=  CONST_IMAG(A, j * lda + k);
          /* temp2 = alpha * Ajk */
          const MpIeee temp2_real=  alpha_real * Ajk_real - alpha_imag * Ajk_imag;
          const MpIeee temp2_imag=  alpha_real * Ajk_imag + alpha_imag * Ajk_real;
          const MpIeee Bjk_real=  CONST_REAL(B, j * ldb + k);
          const MpIeee Bjk_imag=  CONST_IMAG(B, j * ldb + k);

          /* Cij += alpha * Aik * conj(Bjk) + conj(alpha) * Bik * conj(Ajk) */
          temp_real += ((temp1_real * Bjk_real + temp1_imag * Bjk_imag)
                        + (Bik_real * temp2_real + Bik_imag * temp2_imag));
          temp_imag += ((temp1_real * (-Bjk_imag) + temp1_imag * Bjk_real)
                        + (Bik_real * (-temp2_imag) + Bik_imag * temp2_real));
        }
        REAL(C, i * ldc + j) += temp_real;
        IMAG(C, i * ldc + j) += temp_imag;
      }

      /* Cii += alpha Aik conj(Bik) + conj(alpha) Bik conj(Aik) */
      {
        MpIeee temp_real=  0.0;
        /* BASE temp_imag = 0.0; */
        for (k = 0; k < K; k++) {
          const MpIeee Aik_real=  CONST_REAL(A, i * lda + k);
          const MpIeee Aik_imag=  CONST_IMAG(A, i * lda + k);
          /* temp1 = alpha * Aik */
          const MpIeee temp1_real=  alpha_real * Aik_real - alpha_imag * Aik_imag;
          const MpIeee temp1_imag=  alpha_real * Aik_imag + alpha_imag * Aik_real;
          const MpIeee Bik_real=  CONST_REAL(B, i * ldb + k);
          const MpIeee Bik_imag=  CONST_IMAG(B, i * ldb + k);
          temp_real += temp1_real * Bik_real + temp1_imag * Bik_imag;
        }

        REAL(C, i * ldc + i) += 2 * temp_real;
        IMAG(C, i * ldc + i) = 0.0;
      }
    }

  } else if (uplo == CblasLower && trans == CblasConjTrans) {

    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        MpIeee Aki_real=  CONST_REAL(A, k * lda + i);
        MpIeee Aki_imag=  CONST_IMAG(A, k * lda + i);
        MpIeee Bki_real=  CONST_REAL(B, k * ldb + i);
        MpIeee Bki_imag=  CONST_IMAG(B, k * ldb + i);
        /* temp1 = alpha * conj(Aki) */
        MpIeee temp1_real=  alpha_real * Aki_real - alpha_imag * (-Aki_imag);
        MpIeee temp1_imag=  alpha_real * (-Aki_imag) + alpha_imag * Aki_real;
        /* temp2 = conj(alpha) * conj(Bki) */
        MpIeee temp2_real=  alpha_real * Bki_real - alpha_imag * Bki_imag;
        MpIeee temp2_imag=  -(alpha_real * Bki_imag + alpha_imag * Bki_real);

        for (j = 0; j < i; j++) {
          MpIeee Akj_real=  CONST_REAL(A, k * lda + j);
          MpIeee Akj_imag=  CONST_IMAG(A, k * lda + j);
          MpIeee Bkj_real=  CONST_REAL(B, k * ldb + j);
          MpIeee Bkj_imag=  CONST_IMAG(B, k * ldb + j);
          /* Cij += alpha * conj(Aki) * Bkj + conj(alpha) * conj(Bki) * Akj */
          REAL(C, i * lda + j) += (temp1_real * Bkj_real - temp1_imag * Bkj_imag)
              + (temp2_real * Akj_real - temp2_imag * Akj_imag);
          IMAG(C, i * lda + j) += (temp1_real * Bkj_imag + temp1_imag * Bkj_real)
              + (temp2_real * Akj_imag + temp2_imag * Akj_real);
        }

        /* Cii += alpha * conj(Aki) * Bki + conj(alpha) * conj(Bki) * Aki */
        {
          REAL(C, i * lda + i) += 2 * (temp1_real * Bki_real - temp1_imag * Bki_imag);
          IMAG(C, i * lda + i) = 0.0;
        }
      }
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
