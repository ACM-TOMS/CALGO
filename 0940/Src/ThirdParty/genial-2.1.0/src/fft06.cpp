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

DEFINE_IFFT(float)
extern "C" void cifft(int n, const float  *p, float  *q) { DataVector<complex<float > >::self X(n,(complex<float >*)p); DataVector<complex<float > >::self Y(n,(complex<float >*)q); complex_ifft(X,Y); }
extern "C" void sifft(int n,       float  *p, float  *q) { DataVector<complex<float > >::self X(n,(complex<float >*)p); DataVector<        float   >::self Y(n,(        float  *)q); real_ifft   (X,Y); }

DEFINE_IFFT(double)
extern "C" void zifft(int n, const double *p, double *q) { DataVector<complex<double> >::self X(n,(complex<double>*)p); DataVector<complex<double> >::self Y(n,(complex<double>*)q); complex_ifft(X,Y); }
extern "C" void difft(int n,       double *p, double *q) { DataVector<complex<double> >::self X(n,(complex<double>*)p); DataVector<        double  >::self Y(n,(        double *)q); real_ifft   (X,Y); }

#endif

