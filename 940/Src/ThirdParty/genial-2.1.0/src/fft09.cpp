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

DEFINE_NORM_FFT(float)
extern "C" void cnfft (int n         , const float  *p, float  *q) { DataVector<float >::self Y(n ,(float  *)q); complex_norm_fft  (data_vector(n ,(complex<float >*)p),Y); }
extern "C" void shnfft(int nx, int ny, const float  *p, float  *q) { DataVector<float >::self Y(ny,(float  *)q); half_real_norm_fft(data_vector(nx,(float          *)p),Y); }
extern "C" void snfft (int n         , const float  *p, float  *q) { DataVector<float >::self Y(n ,(float  *)q); real_norm_fft     (data_vector(n ,(float          *)p),Y); }

#endif
