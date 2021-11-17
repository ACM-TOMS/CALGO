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


#include "signal/fft.h"

#ifdef FFT_PRECOMPILE

DEFINE_ROW_IFFT(float)
extern "C" void cifft1(int m, int n ,         const float  *p, float  *q) { DataMatrix<complex<float > >::self X(m,n ,(complex<float >*)p); DataMatrix<complex<float > >::self Y(m,n ,(complex<float >*)q); row_ifft(X,Y); }
extern "C" void sifft1(int m, int nx, int ny, const float  *p, float  *q) { DataMatrix<complex<float > >::self X(m,nx,(complex<float >*)p); DataMatrix<        float   >::self Y(m,ny,(        float  *)q); row_ifft(X,Y); }

DEFINE_IFFT2(float)
extern "C" void cifft2(int m, int n ,         const float  *p, float  *q) { DataMatrix<complex<float > >::self X(m,n ,(complex<float >*)p); DataMatrix<complex<float > >::self Y(m,n ,(complex<float >*)q); complex_ifft(X,Y); }
extern "C" void sifft2(int m, int nx, int ny, const float  *p, float  *q) { DataMatrix<complex<float > >::self X(m,nx,(complex<float >*)p); DataMatrix<        float   >::self Y(m,ny,(        float  *)q); real_ifft   (X,Y); }

#endif

