#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_tpsv_r.h
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
  const int nonunit = (Diag == CblasNonUnit);
  const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;

  if (N == 0)
    return;

  /* form  x := inv( A )*x */

  if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasUpper)
      || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasLower)) {
    /* backsubstitution */
    INDEX ix = OFFSET(N, incX) + incX * (N - 1);
    if (nonunit) {
      X[ix] = X[ix] / Ap[TPUP(N, (N - 1), (N - 1))];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      MpIeee tmp=  X[ix];
      INDEX jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const MpIeee Aij=  Ap[TPUP(N, i, j)];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[TPUP(N, i, i)];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasNoTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasTrans && Uplo == CblasUpper)) {

    /* forward substitution */
    INDEX ix = OFFSET(N, incX);
    if (nonunit) {
      X[ix] = X[ix] / Ap[TPLO(N, 0, 0)];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      MpIeee tmp=  X[ix];
      INDEX jx = OFFSET(N, incX);
      for (j = 0; j < i; j++) {
        const MpIeee Aij=  Ap[TPLO(N, i, j)];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[TPLO(N, i, j)];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasUpper)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasLower)) {

    /* form  x := inv( A' )*x */

    /* forward substitution */
    INDEX ix = OFFSET(N, incX);
    if (nonunit) {
      X[ix] = X[ix] / Ap[TPUP(N, 0, 0)];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      MpIeee tmp=  X[ix];
      INDEX jx = OFFSET(N, incX);
      for (j = 0; j < i; j++) {
        const MpIeee Aji=  Ap[TPUP(N, j, i)];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[TPUP(N, i, i)];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans && Uplo == CblasLower)
             || (order == CblasColMajor && Trans == CblasNoTrans && Uplo == CblasUpper)) {

    /* backsubstitution */
    INDEX ix = OFFSET(N, incX) + (N - 1) * incX;
    if (nonunit) {
      X[ix] = X[ix] / Ap[TPLO(N, (N - 1), (N - 1))];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      MpIeee tmp=  X[ix];
      INDEX jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const MpIeee Aji=  Ap[TPLO(N, j, i)];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[TPLO(N, i, i)];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }

}
