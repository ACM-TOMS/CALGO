//GENIAL - GENeric Image Array Library
//Copyright (C) 2006  Patrick LAURENT
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

DEFINE_FFT(double)
extern "C" void zfft  (int n         , const double *p, double *q) { DataVector<complex<double> >::self Y(  n ,(complex<double>*)q); complex_fft  (data_vector(n ,(complex<double>*)p),Y); }
extern "C" void dhfft (int nx, int ny, const double *p, double *q) { DataVector<complex<double> >::self Y(  ny,(complex<double>*)q); half_real_fft(data_vector(nx,(double         *)p),Y); }
extern "C" void dfft  (int n         , const double *p, double *q) { DataVector<complex<double> >::self Y(  n ,(complex<double>*)q); real_fft     (data_vector(n ,(double         *)p),Y); }

#endif
