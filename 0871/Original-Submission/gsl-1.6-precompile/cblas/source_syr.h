#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_syr.h
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

  if (N == 0)
    return;

  if (alpha == 0.0)
    return;

  if ((order == CblasRowMajor && Uplo == CblasUpper)
      || (order == CblasColMajor && Uplo == CblasLower)) {
    INDEX ix = OFFSET(N, incX);
    for (i = 0; i < N; i++) {
      const MpIeee tmp=  alpha * X[ix];
      INDEX jx = ix;
      for (j = i; j < N; j++) {
        A[lda * i + j] += X[jx] * tmp;
        jx += incX;
      }
      ix += incX;
    }
  } else if ((order == CblasRowMajor && Uplo == CblasLower)
             || (order == CblasColMajor && Uplo == CblasUpper)) {
    INDEX ix = OFFSET(N, incX);
    for (i = 0; i < N; i++) {
      const MpIeee tmp=  alpha * X[ix];
      INDEX jx = OFFSET(N, incX);
      for (j = 0; j <= i; j++) {
        A[lda * i + j] += X[jx] * tmp;
        jx += incX;
      }
      ix += incX;
    }
  } else {
    BLAS_ERROR("unrecognized operation");
  }
}
