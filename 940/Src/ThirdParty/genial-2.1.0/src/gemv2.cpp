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


#include "blas/gemv.h"


#ifdef BLAS_PRECOMPILE

DEFINE_BLAS_GEMV(double)
DEFINE_BLAS_GEMV(complex<double>)

extern "C" void BLAS_NAME(dgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const double a, const double *pA, int lda, const double *pX, int dx, const double b, double *pY, int dy) { return gemv(Order,TransA,m,n,a,pA,lda,pX,dx,b,pY,dy); }
extern "C" void BLAS_NAME(zgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const void  *a, const void   *pA, int lda, const void   *pX, int dx, const void  *b, void   *pY, int dy) { return gemv(Order,TransA,m,n,*(complex<double> *)a,(complex<double> *)pA,lda,(complex<double> *)pX,dx,*(complex<double> *)b,(complex<double> *)pY,dy); }

#endif



