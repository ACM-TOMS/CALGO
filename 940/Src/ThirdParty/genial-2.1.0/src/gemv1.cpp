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

DEFINE_BLAS_GEMV(float)
DEFINE_BLAS_GEMV(complex<float>)

extern "C" void BLAS_NAME(sgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const float  a, const float  *pA, int lda, const float  *pX, int dx, const float  b, float  *pY, int dy) { return gemv(Order,TransA,m,n,a,pA,lda,pX,dx,b,pY,dy); }
extern "C" void BLAS_NAME(cgemv)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, int m, int n, const void  *a, const void   *pA, int lda, const void   *pX, int dx, const void  *b, void   *pY, int dy) { return gemv(Order,TransA,m,n,*(complex<float > *)a,(complex<float > *)pA,lda,(complex<float > *)pX,dx,*(complex<float > *)b,(complex<float > *)pY,dy); }

#endif



