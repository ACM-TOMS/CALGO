#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

/* blas/source_rotmg.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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
  INDEX n;
  INDEX i = OFFSET(N, incX);
  INDEX j = OFFSET(N, incY);

  MpIeee h11;MpIeee  h21;MpIeee  h12;MpIeee  h22;

  if (P[0] == -1.0) {
    h11 = P[1];
    h21 = P[2];
    h12 = P[3];
    h22 = P[4];
  } else if (P[0] == 0.0) {
    h11 = 1.0;
    h21 = P[2];
    h12 = P[3];
    h22 = 1.0;
  } else if (P[0] == 1.0) {
    h11 = P[1];
    h21 = -1.0;
    h12 = 1.0;
    h22 = P[4];
  } else if (P[0] == -2.0) {
    return;
  } else {
    BLAS_ERROR("unrecognized value of P[0]");
    return;
  }

  for (n = 0; n < N; n++) {
    const MpIeee w=  X[i];
    const MpIeee z=  Y[j];
    X[i] = h11 * w + h12 * z;
    Y[j] = h21 * w + h22 * z;
    i += incX;
    j += incY;
  }

}
