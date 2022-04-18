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

DEFINE_BLAS_GEMM(double)
DEFINE_BLAS_GEMM(complex<double>)

extern "C" void BLAS_NAME(dgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const double a, const double *pA, int lda, const double *pB, int ldb, const double b, double *pC, int ldc) { return gemm(Order,TransA,TransB,m,n,k,a,pA,lda,pB,ldb,b,pC,ldc); }
extern "C" void BLAS_NAME(zgemm)(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int m, int n, int k, const void  *a, const void   *pA, int lda, const void   *pB, int ldb, const void  *b, void   *pC, int ldc) { return gemm(Order,TransA,TransB,m,n,k,*(complex<double> *)a,(complex<double> *)pA,lda,(complex<double> *)pB,ldb,*(complex<double> *)b,(complex<double> *)pC,ldc); }

#endif



