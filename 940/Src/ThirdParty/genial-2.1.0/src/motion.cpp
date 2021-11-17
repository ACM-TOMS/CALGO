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


#include "image/motion.h"

#ifdef MOTION_ESTIMATION_PRECOMPILE

extern "C" void motion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int m,int n, int *p)
{
  DataMatrix<MatrixIndex<int> >::self M(mx/8,nx/8,(MatrixIndex<int>*)p); 
  M=motion<8,8>(data_matrix(mx,nx,(unsigned char *)px),data_matrix(mx,nx,(unsigned char *)py),m,n);
}

extern "C" void nmotion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int m,int n, int *p)
{
  DataMatrix<MatrixIndex<int> >::self M(mx/8,nx/8,(MatrixIndex<int>*)p); 
  M=neighbourhood_motion<8,8>(data_matrix(mx,nx,(unsigned char *)px),data_matrix(mx,nx,(unsigned char *)py),m,n);
}

extern "C" void qmotion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int l, int m,int n, int m0,int n0, int *p)
{
  DataMatrix<MatrixIndex<int> >::self M(mx/8,nx/8,(MatrixIndex<int>*)p); 
  M=quad_motion<8,8>(data_matrix(mx,nx,(unsigned char *)px),data_matrix(mx,nx,(unsigned char *)py),l,m,n,m0,n0);
}

#endif