//GENIAL - GENeric Image & Array Library
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


#include "blas/copy.h"


#ifdef BLAS_PRECOMPILE

DEFINE_BLAS_COPY(float)
DEFINE_BLAS_COPY(double)
DEFINE_BLAS_COPY(complex<float>)
DEFINE_BLAS_COPY(complex<double>)

DEFINE_BLAS_FILL(float)
DEFINE_BLAS_FILL(double)
DEFINE_BLAS_FILL(complex<float>)
DEFINE_BLAS_FILL(complex<double>)

//{unsecret}
extern "C" void BLAS_NAME(sscal)(int n, const float  a, float  *pX, int dx) { return scal(n,a,pX,dx); }
extern "C" void BLAS_NAME(dscal)(int n, const double a, double *pX, int dx) { return scal(n,a,pX,dx); }
extern "C" void BLAS_NAME(cscal)(int n, const void  *a, void   *pX, int dx) { return scal(n,*(complex<float > *)a,(complex<float > *)pX,dx); }
extern "C" void BLAS_NAME(zscal)(int n, const void  *a, void   *pX, int dx) { return scal(n,*(complex<double> *)a,(complex<double> *)pX,dx); }

extern "C" void BLAS_NAME(scopy)(int n, const float  *pX, int dx, float  *pY, int dy) { return copy(n,pX,dx,pY,dy); }
extern "C" void BLAS_NAME(dcopy)(int n, const double *pX, int dx, double *pY, int dy) { return copy(n,pX,dx,pY,dy); }
extern "C" void BLAS_NAME(ccopy)(int n, const void   *pX, int dx, void   *pY, int dy) { return copy(n,(complex<float > *)pX,dx,(complex<float > *)pY,dy); }
extern "C" void BLAS_NAME(zcopy)(int n, const void   *pX, int dx, void   *pY, int dy) { return copy(n,(complex<double> *)pX,dx,(complex<double> *)pY,dy); }

extern "C" void BLAS_NAME(saxpy)(int n, const float  a, const float  *pX, int dx, float  *pY, int dy) { return axpy(n,a,pX,dx,pY,dy); }
extern "C" void BLAS_NAME(daxpy)(int n, const double a, const double *pX, int dx, double *pY, int dy) { return axpy(n,a,pX,dx,pY,dy); }
extern "C" void BLAS_NAME(caxpy)(int n, const void  *a, const void   *pX, int dx, void   *pY, int dy) { return axpy(n,*(complex<float > *)a,(complex<float > *)pX,dx,(complex<float > *)pY,dy); }
extern "C" void BLAS_NAME(zaxpy)(int n, const void  *a, const void   *pX, int dx, void   *pY, int dy) { return axpy(n,*(complex<double> *)a,(complex<double> *)pX,dx,(complex<double> *)pY,dy); }

#endif



