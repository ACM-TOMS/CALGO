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

DEFINE_ROW_FFT(double)
extern "C" void zfft1 (int m, int n         , const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,n ,(complex<double>*)q); row_fft          (data_matrix(m,n ,(complex<double>*)p),Y); }
extern "C" void dhfft1(int m, int nx, int ny, const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,ny,(complex<double>*)q); half_real_row_fft(data_matrix(m,nx,(double         *)p),Y); }
extern "C" void dfft1 (int m, int n         , const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,n ,(complex<double>*)q); row_fft          (data_matrix(m,n ,(double         *)p),Y); }

DEFINE_FFT2(double)
extern "C" void zfft2 (int m, int n         , const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,n ,(complex<double>*)q); complex_fft  (data_matrix(m,n ,(complex<double>*)p),Y); }
extern "C" void dhfft2(int m, int nx, int ny, const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,ny,(complex<double>*)q); half_real_fft(data_matrix(m,nx,(double         *)p),Y); }
extern "C" void dfft2 (int m, int n         , const double *p, double *q) { DataMatrix<complex<double> >::self Y(m,n ,(complex<double>*)q); real_fft     (data_matrix(m,n ,(double         *)p),Y); }

#endif
