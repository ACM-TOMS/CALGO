//GENIAL - GENeric Image Array Library
//Copyright (C) 2007  Patrick LAURENT
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#include "signal/dct.h"

#ifdef FFT_PRECOMPILE

DEFINE_DCT(float)
extern "C" void sdct  (      int n, const float  *p, float  *q) { DataVector<float >::self Y(  n,q); real_dct (data_vector(  n,(float  *)p),Y); }
extern "C" void sidct (      int n, const float  *p, float  *q) { DataVector<float >::self Y(  n,q); real_idct(data_vector(  n,(float  *)p),Y); }
extern "C" void sdct2 (int m,int n, const float  *p, float  *q) { DataMatrix<float >::self Y(m,n,q); dct (data_matrix(m,n,(float  *)p),Y); }
extern "C" void sidct2(int m,int n, const float  *p, float  *q) { DataMatrix<float >::self Y(m,n,q); idct(data_matrix(m,n,(float  *)p),Y); }

DEFINE_DCT(double)
extern "C" void ddct  (      int n, const double *p, double *q) { DataVector<double>::self Y(  n,q); real_dct (data_vector(  n,(double *)p),Y); }
extern "C" void didct (      int n, const double *p, double *q) { DataVector<double>::self Y(  n,q); real_idct(data_vector(  n,(double *)p),Y); }
extern "C" void ddct2 (int m,int n, const double *p, double *q) { DataMatrix<double>::self Y(m,n,q); dct      (data_matrix(m,n,(double *)p),Y); }
extern "C" void didct2(int m,int n, const double *p, double *q) { DataMatrix<double>::self Y(m,n,q); idct     (data_matrix(m,n,(double *)p),Y); }

#endif
