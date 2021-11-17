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


#include "blas/gemm.h"


#ifdef BLAS_PRECOMPILE

DEFINE_BLAS_GEMM(float)
DEFINE_BLAS_GEMM(complex<float>)

extern "C" void BLAS_NAME(sgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const float  a, const float  *pA, int lda, const float  *pB, int ldb, const float  b, float  *pC, int ldc) { return gemm(Order,TransA,TransB,m,n,k,a,pA,lda,pB,ldb,b,pC,ldc); }
extern "C" void BLAS_NAME(cgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const void  *a, const void   *pA, int lda, const void   *pB, int ldb, const void  *b, void   *pC, int ldc) { return gemm(Order,TransA,TransB,m,n,k,*(complex<float > *)a,(complex<float > *)pA,lda,(complex<float > *)pB,ldb,*(complex<float > *)b,(complex<float > *)pC,ldc); }

#endif



