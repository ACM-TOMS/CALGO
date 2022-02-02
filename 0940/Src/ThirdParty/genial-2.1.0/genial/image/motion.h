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

#ifndef MOTION_ESTIMATION_H
#define MOTION_ESTIMATION_H

#ifdef __cplusplus

#include "array/matrix.h"
#include "signal/conv.h"

#include "image/draw.h"
//namespace genial
//{

//group=Motion Estimation

#ifndef MOTION_ESTIMATION_LEVEL
#define MOTION_ESTIMATION_LEVEL 32
#endif


////////////////////////////////////////////////////////////////////////////////
template<int M,int N,class G1,class G2>
struct block_dist_function
{
};

#ifdef MMX
template<class G1,class G2>
struct block_dist_function<8,8,G1,G2>
{
  inline PROMOTE2(short,typename Matrix<G1>::value_type) operator()(const Matrix<G1> &X,const Matrix<G2> &Y)
  {
    assert(X.size()==typename Matrix<G1>::size_type(8,8));
    assert(Y.size()==typename Matrix<G2>::size_type(8,8));

    SimdVector<4,short>::self d;
    
    SimdVector<8,unsigned char >::self x,y;
    loadu(x,&X(0,0)); loadu(y,&Y(0,0)); d =sad(x,y);
    loadu(x,&X(1,0)); loadu(y,&Y(1,0)); d+=sad(x,y);
    loadu(x,&X(2,0)); loadu(y,&Y(2,0)); d+=sad(x,y);
    loadu(x,&X(3,0)); loadu(y,&Y(3,0)); d+=sad(x,y);
    loadu(x,&X(4,0)); loadu(y,&Y(4,0)); d+=sad(x,y);
    loadu(x,&X(5,0)); loadu(y,&Y(5,0)); d+=sad(x,y);
    loadu(x,&X(6,0)); loadu(y,&Y(6,0)); d+=sad(x,y);
    loadu(x,&X(7,0)); loadu(y,&Y(7,0)); d+=sad(x,y);
    
    return get_sad(d);
  }
};
#endif

#ifdef SSE2
template<class G1,class G2>
struct block_dist_function<16,16,G1,G2>
{
  inline PROMOTE2(short,typename Matrix<G1>::value_type) operator()(const Matrix<G1> &X,const Matrix<G2> &Y)
  {
    assert(X.size()==typename Matrix<G1>::size_type(16,16));
    assert(Y.size()==typename Matrix<G2>::size_type(16,16));
  
    SimdVector<8,short>::self d;

    SimdVector<16,unsigned char >::self x,y;
    loadu(x,&X( 0,0)); loadu(y,&Y( 0,0)); d =sad(x,y);
    loadu(x,&X( 1,0)); loadu(y,&Y( 1,0)); d+=sad(x,y);
    loadu(x,&X( 2,0)); loadu(y,&Y( 2,0)); d+=sad(x,y);
    loadu(x,&X( 3,0)); loadu(y,&Y( 3,0)); d+=sad(x,y);
    loadu(x,&X( 4,0)); loadu(y,&Y( 4,0)); d+=sad(x,y);
    loadu(x,&X( 5,0)); loadu(y,&Y( 5,0)); d+=sad(x,y);
    loadu(x,&X( 6,0)); loadu(y,&Y( 6,0)); d+=sad(x,y);
    loadu(x,&X( 7,0)); loadu(y,&Y( 7,0)); d+=sad(x,y);
    loadu(x,&X( 8,0)); loadu(y,&Y( 8,0)); d+=sad(x,y);
    loadu(x,&X( 9,0)); loadu(y,&Y( 9,0)); d+=sad(x,y);
    loadu(x,&X(10,0)); loadu(y,&Y(10,0)); d+=sad(x,y);
    loadu(x,&X(11,0)); loadu(y,&Y(11,0)); d+=sad(x,y);
    loadu(x,&X(12,0)); loadu(y,&Y(12,0)); d+=sad(x,y);
    loadu(x,&X(13,0)); loadu(y,&Y(13,0)); d+=sad(x,y);
    loadu(x,&X(14,0)); loadu(y,&Y(14,0)); d+=sad(x,y);
    loadu(x,&X(15,0)); loadu(y,&Y(15,0)); d+=sad(x,y);
    
    return get_sad(d);
  }
};
#endif

template<int M,int N,class G1,class G2> 
inline PROMOTE2(short,typename Matrix<G1>::value_type)
block_dist(const Matrix<G1> &X,const Matrix<G2> &Y)
{
	assert(X.size()==Y.size());
  return block_dist_function<M,N,G1,G2>()(X,Y);
}


#ifdef MMX
template<int H,class G>
inline void block_sad_init(const SimdVector<8,unsigned char>::self &x0,const SimdVector<8,unsigned char>::self &x1, const SimdVector<8,unsigned char>::self &x2,const SimdVector<8,unsigned char>::self &x3, const SimdVector<8,unsigned char>::self &x4,const SimdVector<8,unsigned char>::self &x5, const SimdVector<8,unsigned char>::self &x6,const SimdVector<8,unsigned char>::self &x7,
                                 SimdVector<8,unsigned char>::self &y0,      SimdVector<8,unsigned char>::self &y1,       SimdVector<8,unsigned char>::self &y2,      SimdVector<8,unsigned char>::self &y3,       SimdVector<8,unsigned char>::self &y4,      SimdVector<8,unsigned char>::self &y5,       SimdVector<8,unsigned char>::self &y6,      SimdVector<8,unsigned char>::self &y7,
                           const Matrix<G> &Y, SimdVector<4,short>::self &d, SimdVector<4,short>::self &dm)
{
  typedef SimdVector<4,short>::self array_type;
  array_type s,d2,b;
  loadu(y7,&Y( 7,0)); s=sad(y0,x0);s+=sad(y1,x1);s+=sad(y2,x2);s+=sad(y3,x3);s+=sad(y4,x4);s+=sad(y5,x5);s+=sad(y6,x6);s+=sad(y7,x7); d =s;             
  loadu(y0,&Y( 8,0)); s=sad(y1,x0);s+=sad(y2,x1);s+=sad(y3,x2);s+=sad(y4,x3);s+=sad(y5,x4);s+=sad(y6,x5);s+=sad(y7,x6);s+=sad(y0,x7); d+=shiftr<1>(s); 
  loadu(y1,&Y( 9,0)); s=sad(y2,x0);s+=sad(y3,x1);s+=sad(y4,x2);s+=sad(y5,x3);s+=sad(y6,x4);s+=sad(y7,x5);s+=sad(y0,x6);s+=sad(y1,x7); d+=shiftr<2>(s); 
  loadu(y2,&Y(10,0)); s=sad(y3,x0);s+=sad(y4,x1);s+=sad(y5,x2);s+=sad(y6,x3);s+=sad(y7,x4);s+=sad(y0,x5);s+=sad(y1,x6);s+=sad(y2,x7); d+=shiftr<3>(s); 
  dm=array_type(0,1,2,3);
  
  loadu(y3,&Y(11,0)); s=sad(y4,x0);s+=sad(y5,x1);s+=sad(y6,x2);s+=sad(y7,x3);s+=sad(y0,x4);s+=sad(y1,x5);s+=sad(y2,x6);s+=sad(y3,x7); d2 =s;             
  loadu(y4,&Y(12,0)); s=sad(y5,x0);s+=sad(y6,x1);s+=sad(y7,x2);s+=sad(y0,x3);s+=sad(y1,x4);s+=sad(y2,x5);s+=sad(y3,x6);s+=sad(y4,x7); d2+=shiftr<1>(s); 
  loadu(y5,&Y(13,0)); s=sad(y6,x0);s+=sad(y7,x1);s+=sad(y0,x2);s+=sad(y1,x3);s+=sad(y2,x4);s+=sad(y3,x5);s+=sad(y4,x6);s+=sad(y5,x7); d2+=shiftr<2>(s); 
  loadu(y6,&Y(14,0)); s=sad(y7,x0);s+=sad(y0,x1);s+=sad(y1,x2);s+=sad(y2,x3);s+=sad(y3,x4);s+=sad(y4,x5);s+=sad(y5,x6);s+=sad(y6,x7); d2+=shiftr<3>(s); 
  b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(4,5,6,7));
}

template<int H,int D,class G>
inline void block_sad(const SimdVector<8,unsigned char>::self &x0,const SimdVector<8,unsigned char>::self &x1, const SimdVector<8,unsigned char>::self &x2,const SimdVector<8,unsigned char>::self &x3, const SimdVector<8,unsigned char>::self &x4,const SimdVector<8,unsigned char>::self &x5, const SimdVector<8,unsigned char>::self &x6,const SimdVector<8,unsigned char>::self &x7,
                            SimdVector<8,unsigned char>::self &y0,      SimdVector<8,unsigned char>::self &y1,       SimdVector<8,unsigned char>::self &y2,      SimdVector<8,unsigned char>::self &y3,       SimdVector<8,unsigned char>::self &y4,      SimdVector<8,unsigned char>::self &y5,       SimdVector<8,unsigned char>::self &y6,      SimdVector<8,unsigned char>::self &y7,
                      const Matrix<G> &Y, SimdVector<4,short>::self &d, SimdVector<4,short>::self &dm)
{
  typedef SimdVector<4,short>::self array_type;
  array_type s,d2,b;
  loadu(y7,&Y(D+ 7,0)); s=sad(y0,x0);s+=sad(y1,x1);s+=sad(y2,x2);s+=sad(y3,x3);s+=sad(y4,x4);s+=sad(y5,x5);s+=sad(y6,x6);s+=sad(y7,x7); d2 =s;             
  loadu(y0,&Y(D+ 8,0)); s=sad(y1,x0);s+=sad(y2,x1);s+=sad(y3,x2);s+=sad(y4,x3);s+=sad(y5,x4);s+=sad(y6,x5);s+=sad(y7,x6);s+=sad(y0,x7); d2+=shiftr<1>(s); 
  loadu(y1,&Y(D+ 9,0)); s=sad(y2,x0);s+=sad(y3,x1);s+=sad(y4,x2);s+=sad(y5,x3);s+=sad(y6,x4);s+=sad(y7,x5);s+=sad(y0,x6);s+=sad(y1,x7); d2+=shiftr<2>(s); 
  loadu(y2,&Y(D+10,0)); s=sad(y3,x0);s+=sad(y4,x1);s+=sad(y5,x2);s+=sad(y6,x3);s+=sad(y7,x4);s+=sad(y0,x5);s+=sad(y1,x6);s+=sad(y2,x7); d2+=shiftr<3>(s); 
  if (D  <=H/2) { b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(D  ,D+1,D+2,D+3)); }
           else { b=d>d2; d=min(d,d2); dm=and_not(b,dm)|(b&array_type(D  ,D+1,D+2,D+3)); }

  loadu(y3,&Y(D+11,0)); s=sad(y4,x0);s+=sad(y5,x1);s+=sad(y6,x2);s+=sad(y7,x3);s+=sad(y0,x4);s+=sad(y1,x5);s+=sad(y2,x6);s+=sad(y3,x7); d2 =s;             
  loadu(y4,&Y(D+12,0)); s=sad(y5,x0);s+=sad(y6,x1);s+=sad(y7,x2);s+=sad(y0,x3);s+=sad(y1,x4);s+=sad(y2,x5);s+=sad(y3,x6);s+=sad(y4,x7); d2+=shiftr<1>(s); 
  loadu(y5,&Y(D+13,0)); s=sad(y6,x0);s+=sad(y7,x1);s+=sad(y0,x2);s+=sad(y1,x3);s+=sad(y2,x4);s+=sad(y3,x5);s+=sad(y4,x6);s+=sad(y5,x7); d2+=shiftr<2>(s); 
  loadu(y6,&Y(D+14,0)); s=sad(y7,x0);s+=sad(y0,x1);s+=sad(y1,x2);s+=sad(y2,x3);s+=sad(y3,x4);s+=sad(y4,x5);s+=sad(y5,x6);s+=sad(y6,x7); d2+=shiftr<3>(s); 
  if (D+4<=H/2) { b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(D+4,D+5,D+6,D+7)); }
           else { b=d>d2; d=min(d,d2); dm=and_not(b,dm)|(b&array_type(D+4,D+5,D+6,D+7)); }
}
#endif // MMX

#ifdef SSE2
template<int H,class G>
inline void block_sad_init(const SimdVector<16,unsigned char>::self &x0,const SimdVector<16,unsigned char>::self &x1, const SimdVector<16,unsigned char>::self &x2,const SimdVector<16,unsigned char>::self &x3, const SimdVector<16,unsigned char>::self &x4,const SimdVector<16,unsigned char>::self &x5, const SimdVector<16,unsigned char>::self &x6,const SimdVector<16,unsigned char>::self &x7,
                                 SimdVector<16,unsigned char>::self &y0,      SimdVector<16,unsigned char>::self &y1,       SimdVector<16,unsigned char>::self &y2,      SimdVector<16,unsigned char>::self &y3,       SimdVector<16,unsigned char>::self &y4,      SimdVector<16,unsigned char>::self &y5,       SimdVector<16,unsigned char>::self &y6,      SimdVector<16,unsigned char>::self &y7,
                           const Matrix<G> &Y, SimdVector<8,short>::self &d, SimdVector<8,short>::self &dm)
{
  typedef SimdVector<8,short>::self array_type;
  array_type s,d2,b;
  loadu(y7,&Y( 7,0)); s=sad(y0,x0);s+=sad(y1,x1);s+=sad(y2,x2);s+=sad(y3,x3);s+=sad(y4,x4);s+=sad(y5,x5);s+=sad(y6,x6);s+=sad(y7,x7); d =s;             
  loadu(y0,&Y( 8,0)); s=sad(y1,x0);s+=sad(y2,x1);s+=sad(y3,x2);s+=sad(y4,x3);s+=sad(y5,x4);s+=sad(y6,x5);s+=sad(y7,x6);s+=sad(y0,x7); d+=shiftr<1>(s); 
  loadu(y1,&Y( 9,0)); s=sad(y2,x0);s+=sad(y3,x1);s+=sad(y4,x2);s+=sad(y5,x3);s+=sad(y6,x4);s+=sad(y7,x5);s+=sad(y0,x6);s+=sad(y1,x7); d+=shiftr<2>(s); 
  loadu(y2,&Y(10,0)); s=sad(y3,x0);s+=sad(y4,x1);s+=sad(y5,x2);s+=sad(y6,x3);s+=sad(y7,x4);s+=sad(y0,x5);s+=sad(y1,x6);s+=sad(y2,x7); d+=shiftr<3>(s); 
  dm=array_type(m128s(0,1,2,3,0,1,2,3));
  
  loadu(y3,&Y(11,0)); s=sad(y4,x0);s+=sad(y5,x1);s+=sad(y6,x2);s+=sad(y7,x3);s+=sad(y0,x4);s+=sad(y1,x5);s+=sad(y2,x6);s+=sad(y3,x7); d2 =s;             
  loadu(y4,&Y(12,0)); s=sad(y5,x0);s+=sad(y6,x1);s+=sad(y7,x2);s+=sad(y0,x3);s+=sad(y1,x4);s+=sad(y2,x5);s+=sad(y3,x6);s+=sad(y4,x7); d2+=shiftr<1>(s); 
  loadu(y5,&Y(13,0)); s=sad(y6,x0);s+=sad(y7,x1);s+=sad(y0,x2);s+=sad(y1,x3);s+=sad(y2,x4);s+=sad(y3,x5);s+=sad(y4,x6);s+=sad(y5,x7); d2+=shiftr<2>(s); 
  loadu(y6,&Y(14,0)); s=sad(y7,x0);s+=sad(y0,x1);s+=sad(y1,x2);s+=sad(y2,x3);s+=sad(y3,x4);s+=sad(y4,x5);s+=sad(y5,x6);s+=sad(y6,x7); d2+=shiftr<3>(s); 
  b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(m128s(4,5,6,7,4,5,6,7)));
}

template<int H,int D,class G>
inline void block_sad(const SimdVector<16,unsigned char>::self &x0,const SimdVector<16,unsigned char>::self &x1, const SimdVector<16,unsigned char>::self &x2,const SimdVector<16,unsigned char>::self &x3, const SimdVector<16,unsigned char>::self &x4,const SimdVector<16,unsigned char>::self &x5, const SimdVector<16,unsigned char>::self &x6,const SimdVector<16,unsigned char>::self &x7,
                            SimdVector<16,unsigned char>::self &y0,      SimdVector<16,unsigned char>::self &y1,       SimdVector<16,unsigned char>::self &y2,      SimdVector<16,unsigned char>::self &y3,       SimdVector<16,unsigned char>::self &y4,      SimdVector<16,unsigned char>::self &y5,       SimdVector<16,unsigned char>::self &y6,      SimdVector<16,unsigned char>::self &y7,
                      const Matrix<G> &Y, SimdVector<8,short>::self &d, SimdVector<8,short>::self &dm)
{
  typedef SimdVector<8,short>::self array_type;
  array_type s,d2,b;
  loadu(y7,&Y(D+ 7,0)); s=sad(y0,x0);s+=sad(y1,x1);s+=sad(y2,x2);s+=sad(y3,x3);s+=sad(y4,x4);s+=sad(y5,x5);s+=sad(y6,x6);s+=sad(y7,x7); d2 =s;             
  loadu(y0,&Y(D+ 8,0)); s=sad(y1,x0);s+=sad(y2,x1);s+=sad(y3,x2);s+=sad(y4,x3);s+=sad(y5,x4);s+=sad(y6,x5);s+=sad(y7,x6);s+=sad(y0,x7); d2+=shiftr<1>(s); 
  loadu(y1,&Y(D+ 9,0)); s=sad(y2,x0);s+=sad(y3,x1);s+=sad(y4,x2);s+=sad(y5,x3);s+=sad(y6,x4);s+=sad(y7,x5);s+=sad(y0,x6);s+=sad(y1,x7); d2+=shiftr<2>(s); 
  loadu(y2,&Y(D+10,0)); s=sad(y3,x0);s+=sad(y4,x1);s+=sad(y5,x2);s+=sad(y6,x3);s+=sad(y7,x4);s+=sad(y0,x5);s+=sad(y1,x6);s+=sad(y2,x7); d2+=shiftr<3>(s); 
  if (D  <=H/2) { b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(m128s(D  ,D+1,D+2,D+3,D  ,D+1,D+2,D+3))); }
           else { b=d>d2; d=min(d,d2); dm=and_not(b,dm)|(b&array_type(m128s(D  ,D+1,D+2,D+3,D  ,D+1,D+2,D+3))); }

  loadu(y3,&Y(D+11,0)); s=sad(y4,x0);s+=sad(y5,x1);s+=sad(y6,x2);s+=sad(y7,x3);s+=sad(y0,x4);s+=sad(y1,x5);s+=sad(y2,x6);s+=sad(y3,x7); d2 =s;             
  loadu(y4,&Y(D+12,0)); s=sad(y5,x0);s+=sad(y6,x1);s+=sad(y7,x2);s+=sad(y0,x3);s+=sad(y1,x4);s+=sad(y2,x5);s+=sad(y3,x6);s+=sad(y4,x7); d2+=shiftr<1>(s); 
  loadu(y5,&Y(D+13,0)); s=sad(y6,x0);s+=sad(y7,x1);s+=sad(y0,x2);s+=sad(y1,x3);s+=sad(y2,x4);s+=sad(y3,x5);s+=sad(y4,x6);s+=sad(y5,x7); d2+=shiftr<2>(s); 
  loadu(y6,&Y(D+14,0)); s=sad(y7,x0);s+=sad(y0,x1);s+=sad(y1,x2);s+=sad(y2,x3);s+=sad(y3,x4);s+=sad(y4,x5);s+=sad(y5,x6);s+=sad(y6,x7); d2+=shiftr<3>(s); 
  if (D+4<=H/2) { b=d<d2; d=min(d,d2); dm=(b&dm)|and_not(b,array_type(m128s(D+4,D+5,D+6,D+7,D+4,D+5,D+6,D+7))); }
           else { b=d>d2; d=min(d,d2); dm=and_not(b,dm)|(b&array_type(m128s(D+4,D+5,D+6,D+7,D+4,D+5,D+6,D+7))); }
}
#endif //SSE2

template<int H,int N,class G> 
inline void motion(const Vector<simd_vector_generator<N,unsigned char> > &x0,const Vector<simd_vector_generator<N,unsigned char> > &x1,const Vector<simd_vector_generator<N,unsigned char> > &x2,const Vector<simd_vector_generator<N,unsigned char> > &x3,const Vector<simd_vector_generator<N,unsigned char> > &x4,const Vector<simd_vector_generator<N,unsigned char> > &x5,const Vector<simd_vector_generator<N,unsigned char> > &x6,const Vector<simd_vector_generator<N,unsigned char> > &x7, 
       const Matrix<G> &Y, Vector<simd_vector_generator<N/2,short> > &d, Vector<simd_vector_generator<N/2,short> > &dm)
{
  typename SimdVector<N,unsigned char >::self y0,y1,y2,y3,y4,y5,y6,y7;
  loadu(y0,&Y(0,0)); loadu(y1,&Y(1,0)); loadu(y2,&Y(2,0)); loadu(y3,&Y(3,0)); loadu(y4,&Y(4,0)); loadu(y5,&Y(5,0)); loadu(y6,&Y(6,0));
  if (H>= 8) block_sad_init<H>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=16) block_sad  <H, 8>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=24) block_sad  <H,16>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=32) block_sad  <H,24>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=40) block_sad  <H,32>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=48) block_sad  <H,40>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=56) block_sad  <H,48>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
  if (H>=64) block_sad  <H,56>(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7,Y,d,dm);
}

#ifdef MMX
template<int H,class G1,class G2> 
pair<typename Matrix<G2>::index_type,short> 
mmx_motion_8x8(const Matrix<G1> &X,const Matrix<G2> &Y)
{
  _mm_empty();
  SimdVector<8,unsigned char>::self x0,x1,x2,x3,x4,x5,x6,x7;
  load(x0,&X(0,0)); load(x1,&X(1,0)); load(x2,&X(2,0)); load(x3,&X(3,0)); load(x4,&X(4,0)); load(x5,&X(5,0)); load(x6,&X(6,0)); load(x7,&X(7,0));
  
  int jmax=Y.ncols()-7;
  SimdVector<4,short>::self dn,dm,d2,di,b, d(32767), dj(-jmax/2);
  for(int j=0; j<jmax; ++j, dj+=SimdVector<4,short>::self(1,1,1,1))
  {
    motion<H>(x0,x1,x2,x3,x4,x5,x6,x7,sub<H+7,8>(Y,0,j),d2,di);
    b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  }
  d2=shiftl<2>(d); di=shiftl<2>(dm); dj=shiftl<2>(dn); b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  d2=shiftl<1>(d); di=shiftl<1>(dm); dj=shiftl<1>(dn); b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  _mm_empty();
  return make_pair(typename Matrix<G2>::index_type(dm[0],dn[0]+jmax/2),d[0]);
}
#endif // MMX

#ifdef SSE2
template<int H,class G1,class G2> 
pair<typename Matrix<G2>::index_type,short> 
sse2_motion_8x8(const Matrix<G1> &X,const Matrix<G2> &Y)
{
  SimdVector<16,unsigned char >::self x0,x1,x2,x3,x4,x5,x6,x7;
  loadlh(x0,&X(0,0)); loadlh(x1,&X(1,0)); loadlh(x2,&X(2,0)); loadlh(x3,&X(3,0)); loadlh(x4,&X(4,0)); loadlh(x5,&X(5,0)); loadlh(x6,&X(6,0)); loadlh(x7,&X(7,0));

  int jmax=Y.ncols()-7;
  SimdVector<8,short>::self dn,dm,d2,di,b, d(32767), dj(m128s(-jmax/2)+m128s(0,0,0,0,8,8,8,8));    
  for (int j=0; j<jmax; j+=8,dj+=SimdVector<8,short>::self(m128s(8,8,8,8,8,8,8,8)))
    for (int jmax=j+8; j<jmax; ++j,dj+=SimdVector<8,short>::self(m128s(1,1,1,1,1,1,1,1)))
    {
      motion<H>(x0,x1,x2,x3,x4,x5,x6,x7,sub<H+7,16>(Y,0,j),d2,di);
      b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
    }
  d2=shiftl<4>(d); di=shiftl<4>(dm); dj=shiftl<4>(dn); b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  d2=shiftl<2>(d); di=shiftl<2>(dm); dj=shiftl<2>(dn); b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  d2=shiftl<1>(d); di=shiftl<1>(dm); dj=shiftl<1>(dn); b=(d2<d)|(d2==d&abs(dj)<abs(dn)); d=min(d,d2); dm=(b&di)|and_not(b,dm); dn=(b&dj)|and_not(b,dn); 
  return make_pair(typename Matrix<G2>::index_type(dm[0],dn[0]+jmax/2),d[0]);
}
#endif // SSE2

template<int M,int N,int M2,int N2,class T1, class T2> struct motion_function {};
#if defined(MMX) && !defined(SSE2)
template<int M2> struct motion_function<8,8,M2,0,unsigned char,unsigned char> { template<class G1,class G2> inline pair<typename Matrix<G2>::index_type,short> operator()(const Matrix<G1> &X,const Matrix<G2> &Y) const { return  mmx_motion_8x8<M2>(X,Y); } };
#elif defined(SSE2) && defined(MMX)
template<int M2> struct motion_function<8,8,M2,8,unsigned char,unsigned char> { template<class G1,class G2> inline pair<typename Matrix<G2>::index_type,short> operator()(const Matrix<G1> &X,const Matrix<G2> &Y) const { return  mmx_motion_8x8<M2>(X,Y); } };
template<int M2> struct motion_function<8,8,M2,0,unsigned char,unsigned char> { template<class G1,class G2> inline pair<typename Matrix<G2>::index_type,short> operator()(const Matrix<G1> &X,const Matrix<G2> &Y) const { return sse2_motion_8x8<M2>(X,Y); } };
#elif defined(SSE2) && !defined(MMX)
template<int M2> struct motion_function<8,8,M2,0,unsigned char,unsigned char> { template<class G1,class G2> inline pair<typename Matrix<G2>::index_type,short> operator()(const Matrix<G1> &X,const Matrix<G2> &Y) const { return sse2_motion_8x8<M2>(X,Y); } };
#endif


template<int M,int N,int M2,int N2,class G1,class G2> 
inline pair<typename Matrix<G2>::index_type,short> motion(const Matrix<G1> &X,const Matrix<G2> &Y) 
{
  return motion_function<M,N,M2,N2,typename Matrix<G1>::value_type,typename Matrix<G2>::value_type>()(X,Y); 
}

template<int M,int N,int M2       ,class G1,class G2> 
inline pair<typename Matrix<G2>::index_type,short> 
motion(const Matrix<G1> &X,const Matrix<G2> &Y) 
{
  return motion<M,N,M2,0>(X,Y); 
}


template<int M,int N,class G1,class G2>
inline pair<typename Matrix<G2>::index_type,short> aux_motion(const Matrix<G1> &X,const Matrix<G2> &Y)
{
  assert(X.ncols()<=Y.ncols());
  assert(X.nrows()<=Y.nrows());
  
  switch (Y.nrows()-(M-1))
  {
    case 8:  
    {
#if defined(SSE2) && defined(MMX)
      if (Y.ncols()-(N-1)==8) return mmx_motion_8x8<8>(X,Y);
#endif      
      return motion<M,N, 8>(X,Y);
    }
#if MOTION_ESTIMATION_LEVEL>=16
    case 16: return motion<M,N,16>(X,Y);
#endif
#if MOTION_ESTIMATION_LEVEL>=32
    case 32: return motion<M,N,32>(X,Y);
#endif
#if MOTION_ESTIMATION_LEVEL>=64
    case 64: return motion<M,N,64>(X,Y);
#endif
    default: throw error("Motion estimation error: height of searcharea not allowed!");
  }
}


template<int M,int N,class A1,class G2> pair<typename A1::index_type,short> motion(const Matrix<tiny_sub_matrix_generator      <M,N,const A1> > &X, const Matrix<G2> &Y) { return aux_motion<M,N>(X,Y); }
template<int M,int N,class A1,class G2> pair<typename A1::index_type,short> motion(const Matrix<tiny_sub_dense_matrix_generator<M,N,const A1> > &X, const Matrix<G2> &Y) { return aux_motion<M,N>(X,Y); }


///////////////////////////////////////////////////////////////////////////////////////////////

template<int M,int N,class A1,class A2,class V>
class motion_array_generator : public two_arrays_value_generator<const A1,const A2,V>
{
  public:
    typedef motion_array_generator self;
    typedef two_arrays_value_generator<const A1,const A2,V> base;
    typedef typename base::first_array_type  first_array_type;
    typedef typename base::second_array_type second_array_type;
    typedef typename base::value_type value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::size_type  size_type;
    typedef typename base::index_type index_type;

  private:
    size_type search_size;

  public:
    inline motion_array_generator(first_array_type &x, second_array_type &y,int m,int n       ) : base(x,y), search_size(m,n) { assert(x.size()==y.size()); }
    inline motion_array_generator(first_array_type &x, second_array_type &y,const size_type &s) : base(x,y), search_size(s  ) { assert(x.size()==y.size()); }

    using base::first_array;
    using base::second_array;

    const_reference operator[](const index_type &i) const 
    {
      assert(!( first_array().nrows()%M)); assert(!( first_array().ncols()%N)); assert(!(second_array().nrows()%M)); assert(!(second_array().ncols()%N));
      index_type p1 = i*index_type(M,N)-search_size/2;
      index_type p2= p1 + search_size + (index_type(M,N)-2);  
      check_rect(p1,p2);
      typename TinySubMatrix<M,N,first_array_type>::self X1(sub<M,N>(first_array(),i*index_type(M,N)));
      typename SubArray<second_array_type>::self Y1(sub(second_array(),p1,size_type(p2-p1+1)));
      return value_type(motion(X1,Y1).first+p1-i*index_type(M,N)); 
    }
 
    inline void check_rect(index_type &p1, index_type &p2) const
    {
      int m = p2.i-p1.i; p1.i=__max(p1.i,0); p2.i=p1.i+m; p2.i=__min(p2.i,first_array().row_upper_bound()); p1.i=p2.i-m;
      int n = p2.j-p1.j; p1.j=__max(p1.j,0); p2.j=p1.j+n; p2.j=__min(p2.j,first_array().col_upper_bound()); p1.j=p2.j-n;
    }
    
    size_type size() const { return second_array().size()/size_type(M,N); }
    void resize(const size_type &s) { assert(s==size()); }
};

template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type>
class neighbourhood_motion_array_generator : public two_arrays_value_generator<const A1,const A2,V>
{
  public:
    typedef neighbourhood_motion_array_generator self;
    typedef two_arrays_value_generator<const A1,const A2,V> base;
    typedef typename base::first_array_type  first_array_type;
    typedef typename base::second_array_type second_array_type;    
    typedef typename base::value_type value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::size_type  size_type;
    typedef typename base::index_type index_type;
    
  private:
    size_type search_size;
    mutable typename DenseMatrix<index_type>::self result_matrix;

  public:
    inline neighbourhood_motion_array_generator(first_array_type &x, second_array_type &y,int m,int n       ) : base(x,y), search_size(m,n), result_matrix(x.size()/size_type(M,N),index_type(0,0)) { assert(x.size()==y.size()); }
    inline neighbourhood_motion_array_generator(first_array_type &x, second_array_type &y,const size_type &s) : base(x,y), search_size(s  ), result_matrix(x.size()/size_type(M,N),index_type(0,0)) { assert(x.size()==y.size()); }
     
    using base::first_array;
    using base::second_array;
     
    const_reference operator[](const index_type &i) const 
    {
      assert(!( first_array().nrows()%M)); assert(!( first_array().ncols()%N)); assert(!(second_array().nrows()%M)); assert(!(second_array().ncols()%N));
      index_type p1 = i*index_type(M,N)-search_size/2;
      index_type p2 = p1 + search_size + (index_type(M,N)-2);
      if (!result_matrix.topboundary(i) && !result_matrix.leftboundary(i) && !result_matrix.rightboundary(i))
      {
        index_type prediction = median(result_matrix(i.i,i.j-1),result_matrix(i.i-1,i.j),result_matrix(i.i-1,i.j+1));
        p1 += prediction;
        p2 += prediction;
      }
      check_rect(p1,p2);
      typename SubArray<second_array_type>::self Y1(sub(second_array(),p1,size_type(p2-p1+1)));
      typename TinySubMatrix<M,N,first_array_type>::self X1(sub<M,N>(first_array(),i*index_type(M,N)));
      return result_matrix[i]=value_type(motion(X1,Y1).first+p1-i*index_type(M,N));
    }
    
    inline void check_rect(index_type &p1, index_type &p2) const
    {
      int m = p2.i-p1.i; p1.i=__max(p1.i,0); p2.i=p1.i+m; p2.i=__min(p2.i,first_array().row_upper_bound()); p1.i=p2.i-m;
      int n = p2.j-p1.j; p1.j=__max(p1.j,0); p2.j=p1.j+n; p2.j=__min(p2.j,first_array().col_upper_bound()); p1.j=p2.j-n;
    }

    size_type size() const { return second_array().size()/size_type(M,N); }
    void resize(const size_type &s) { assert(s==size()); }
};

//template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type> struct MotionArray;
//template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type> struct QuadMotionArray;
//template<int M,int N,class G> inline typename     MotionArray<M,N,typename ArrayData<const Matrix<G> >::self>::self motion     (const Matrix<G> &X,const Matrix<G> &Y,      const typename Matrix<G>::size_type &s);
//template<int M,int N,class G> inline typename QuadMotionArray<M,N,typename ArrayData<const Matrix<G> >::self>::self quad_motion(const Matrix<G> &X,const Matrix<G> &Y,int l,const typename Matrix<G>::size_type &s,const typename Matrix<G>::size_type &s0,const typename Matrix<G>::size_type &si);

template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type>
class quad_motion_array_generator : public two_arrays_value_generator<const A1,const A2,V>
{
  public:
    typedef quad_motion_array_generator self;
    typedef two_arrays_value_generator<const A1,const A2,V> base;
    typedef typename base::first_array_type  first_array_type;  
    typedef typename base::second_array_type second_array_type;    
    typedef typename base::value_type value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::size_type  size_type;
    typedef typename base::index_type index_type;
    typedef pair<index_type,index_type> rect_type;   

  private:
    int layer;
    size_type search_size;
    size_type layer_search_size;
    typename DenseMatrix<index_type>::self predictions;    

  public:
    inline quad_motion_array_generator(first_array_type &x, second_array_type &y,int l,int m,int n                                                   ) : base(x,y), layer(l), search_size(m,n), layer_search_size(m ,n ) { init(first_array().size()); }
    inline quad_motion_array_generator(first_array_type &x, second_array_type &y,int l,int m,int n       ,int ml, int nl                             ) : base(x,y), layer(l), search_size(m,n), layer_search_size(ml,nl) { init(first_array().size()); }
    inline quad_motion_array_generator(first_array_type &x, second_array_type &y,int l,const size_type &s                                            ) : base(x,y), layer(l), search_size(s  ), layer_search_size(s    ) { init(first_array().size()); }
    inline quad_motion_array_generator(first_array_type &x, second_array_type &y,int l,const size_type &s,const size_type &sl                        ) : base(x,y), layer(l), search_size(s  ), layer_search_size(sl   ) { init(first_array().size()); }
    inline quad_motion_array_generator(first_array_type &x, second_array_type &y,int l,const size_type &s,const size_type &sl,const size_type &simage) : base(x,y), layer(l), search_size(s  ), layer_search_size(sl   ) { init(simage); }

    using base::first_array;
    using base::second_array;

    void init(const size_type &simage);

    const_reference operator[](const index_type &i) const 
    {
      assert(!( first_array().nrows()%M)); assert(!( first_array().ncols()%N)); assert(!(second_array().nrows()%M)); assert(!(second_array().ncols()%N));
      index_type p1 = i*index_type(M,N)-search_size/2;
      index_type p2 = p1 + search_size + (index_type(M,N)-2);
      
      index_type prediction = (predictions[i/2])*2;
      p1 += prediction;
      p2 += prediction;

      check_rect(p1,p2);
      typename SubArray<second_array_type>::self Y1(sub(second_array(),p1,size_type(p2-p1+1)));
      typename TinySubMatrix<M,N,first_array_type>::self X1(sub<M,N>(first_array(),i*index_type(M,N)));
      return value_type(motion(X1,Y1).first+p1-i*index_type(M,N));
    }
    
    inline void check_rect(index_type &p1, index_type &p2) const
    {
      int m = p2.i-p1.i; p1.i=__max(p1.i,0); p2.i=p1.i+m; p2.i=__min(p2.i,first_array().row_upper_bound()); p1.i=p2.i-m;
      int n = p2.j-p1.j; p1.j=__max(p1.j,0); p2.j=p1.j+n; p2.j=__min(p2.j,first_array().col_upper_bound()); p1.j=p2.j-n;
    }
     
    size_type size() const { return second_array().size()/size_type(M,N); }
    void resize(const size_type &s) { assert(s==size()); }
};


//////////////////////////////////////////////////////////////////////////////////



template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type>
struct MotionArray
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef V value_type;
  typedef motion_array_generator<M,N,first_array_type, second_array_type, value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
};

template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type>
struct NeighbourhoodMotionArray
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef V value_type;
  typedef neighbourhood_motion_array_generator<M,N,first_array_type, second_array_type, value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
};


template<int M,int N,class A1,class A2=A1,class V=typename A1::index_type>
struct QuadMotionArray
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef V value_type;
  typedef quad_motion_array_generator<M,N,first_array_type, second_array_type, value_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};


template<int M,int N,class G> inline typename              MotionArray<M,N,const Matrix<G> >::self motion              (const Matrix<G> &X,const Matrix<G> &Y,      const typename Matrix<G>::size_type &s                                                                                ) { return typename              MotionArray<M,N,const Matrix<G> >::self(X,Y,  s      ); }
template<int M,int N,class G> inline typename NeighbourhoodMotionArray<M,N,const Matrix<G> >::self neighbourhood_motion(const Matrix<G> &X,const Matrix<G> &Y,      const typename Matrix<G>::size_type &s                                                                                ) { return typename NeighbourhoodMotionArray<M,N,const Matrix<G> >::self(X,Y,  s      ); }
template<int M,int N,class G> inline typename          QuadMotionArray<M,N,const Matrix<G> >::self quad_motion         (const Matrix<G> &X,const Matrix<G> &Y,int l,const typename Matrix<G>::size_type &s                                                                                ) { return typename          QuadMotionArray<M,N,const Matrix<G> >::self(X,Y,l,s      ); }
template<int M,int N,class G> inline typename          QuadMotionArray<M,N,const Matrix<G> >::self quad_motion         (const Matrix<G> &X,const Matrix<G> &Y,int l,const typename Matrix<G>::size_type &s,const typename Matrix<G>::size_type &s0                                        ) { return typename          QuadMotionArray<M,N,const Matrix<G> >::self(X,Y,l,s,s0   ); }
template<int M,int N,class G> inline typename          QuadMotionArray<M,N,const Matrix<G> >::self quad_motion         (const Matrix<G> &X,const Matrix<G> &Y,int l,const typename Matrix<G>::size_type &s,const typename Matrix<G>::size_type &s0,const typename Matrix<G>::size_type &si) { return typename          QuadMotionArray<M,N,const Matrix<G> >::self(X,Y,l,s,s0,si); }


template<int M,int N,class A1,class A2,class V>
void quad_motion_array_generator<M,N,A1,A2,V>::init(const typename quad_motion_array_generator<M,N,A1,A2,V>::size_type &simage) 
{
  assert(first_array().size()==second_array().size());
  
  if (layer==0) { predictions.resize((first_array().nrows()+(2*M-1))/(2*M),(first_array().ncols()+(2*N-1))/(2*N)); predictions=index_type(0,0); return; }
  
  TinyVector<7,short>::self filter("-14 0 44 69 44 0 -14");

  if (layer==1) 
  { 
    DenseMatrix<unsigned char>::self Xq((first_array().nrows()+1)/2,16*((first_array().ncols()/2+15)/16)); 
    DenseMatrix<unsigned char>::self Yq(Xq.size());                                
    sample2_conv<7>(first_array() ,filter,Xq);
    sample2_conv<7>(second_array(),filter,Yq);
    predictions=sub(motion<M,N>(Xq,Yq,layer_search_size),index_type(0,0),simage/size_type(2*M,2*N)); 

    // pas véridié, non général à cause de la dimension 8
    //predictions.resize(simage/size_type(2*M,2*N));
    //motion8(Xq.nrows(),Xq.ncols(),&Xq(0,0),&Yq(0,0), layer_search_size.nrows(), layer_search_size.ncols(), (int *)&predictions(0,0));    
  }
  else
  {
    DenseMatrix<unsigned char>::self Xq((first_array().nrows()+1)/2,32*((first_array().ncols()/2+31)/32));
    DenseMatrix<unsigned char>::self Yq(Xq.size());
    sample2_conv<7>(first_array() ,filter,Xq);
    sample2_conv<7>(second_array(),filter,Yq);
    predictions=sub(quad_motion<M,N>(Xq,Yq,layer-1,search_size,layer_search_size,simage/2),index_type(0,0),simage/size_type(2*M,2*N));  
  }   
}
 


//{unsecret}
//Summary: Motion Estimation
//Return: 
//  A matrix of indexes representing the motion estimation vectors of each MxN block of X.
//  Motion vectors are in the range [-m/2 ... m/2-1] for vertical shifting and [-n/2...n/2-1] for horizontal shifting.
//  Positive values mean displacement to the bottom, respectively to the right.
//Arguments: 
//  M - Height of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  N - Width of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  X - The first image - Only images of unsigned char are currently accepted.
//  Y - The second image - Only images of unsigned char are currently accepted.
//  m - Height of search areas, has to be 8, 16, 32 or 64.
//  n - Width of search areas, has to be 8, 16, 32 or 64.
//Remarks:
//  The seach areas are centered at the corresponding searched blocks.
//  Choose the size of search areas carefully in order to have fair computation time so as a fair motion estimation.
//  SSE2 instuctions deal better for 16x16 and bigger search areas than MMX/SSE instructions.
//
//  {image:motion}
//Example:
//  typedef ucharImage::index_type index_type;
//  RGBImage Image0=..., Image1=...;
//  ucharImage X=first(value_cast<YUV>(Image0)); // Isolates the greyscale component (Y)
//  ucharImage Y=first(value_cast<YUV>(Image1));
//  DenseMatrix<index_type>::self M = motion<8,8>(X,Y,16,16); // stores all motion vectors
//  cout << motion<8,8>(X,Y,32,32)(5,5) << endl; // computes only one motion vector, for block at position (5,5)
//See: ^neigbourhood_motion^, ^quad_motion^
template<int M,int N,class G> inline typename MotionArray<M,N,const Matrix<G> >::self 
motion(const Matrix<G> &X,const Matrix<G> &Y,int m,int n) { return motion<M,N>(data(X),data(Y),typename Matrix<G>::size_type(m,n)); }


//{unsecret}
//Summary: Motion Estimation with neigbourhood predictions
//Return: 
//  A matrix of indexes representing the motion estimation vectors of each MxN block of X.
//  Positive values mean displacement to the bottom, respectively to the right.
//Arguments: 
//  M - Height of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  N - Width of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  X - The first image - Only images of unsigned char are currently accepted
//  Y - The second image - Only images of unsigned char are currently accepted
//  m - Height of search areas, has to be 8, 16, 32 or 64
//  n - Width of search areas, has to be 8, 16, 32 or 64
//Remarks:
//  The seach areas are not centered at the corresponding searched blocks, 
//  they are instead translated with previouly computed motion vectors of neighboor blocks (if already computed).
//  Choose the size of search areas carefully in order to have fair computation time so as a fair motion estimation.
//  SSE2 instuctions deal better for 16x16 and bigger search areas than MMX/SSE instructions.
//
//  {image:neighbourhood_motion}
//Example:
//  typedef ucharImage::index_type index_type;
//  RGBImage Image0=..., Image1=...;
//  ucharImage X=first(value_cast<YUV>(Image0)); // Isolates the greyscale component (Y)
//  ucharImage Y=first(value_cast<YUV>(Image1));
//  DenseMatrix<index_type>::self M = neigbourhood_motion<8,8>(X,Y,16,16); // stores all motion vectors
//See: ^motion^, ^quad_motion^
template<int M,int N,class G> inline typename NeighbourhoodMotionArray<M,N,const Matrix<G> >::self  
neighbourhood_motion(const Matrix<G> &X,const Matrix<G> &Y,int m,int n) { return neighbourhood_motion<M,N>(data(X),data(Y),typename Matrix<G>::size_type(m,n)); }


//{unsecret}
//Summary: Motion Estimation with quad tree predictions
//Return: 
//  A matrix of indexes representing the motion estimation vectors of each MxN block of X.
//  Positive values mean displacement to the bottom, respectively to the right.
//Arguments: 
//  M - Height of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  N - Width of searched blocks. Only 8x8 blocks are currently accepted because of some SSE/SSE2 restrictions.
//  X - The first image - Only images of unsigned char are currently accepted
//  Y - The second image - Only images of unsigned char are currently accepted
//  l - Number of layers (image layer excluded)
//      The value 0 makes /quad_motion/ equivalent to the ^motion^ function.
//  m - Height of search areas, has to be 8, 16, 32 or 64
//  n - Width of search areas, has to be 8, 16, 32 or 64
//  m0 - Height of search areas at first layer, has to be 8, 16, 32 or 64
//  n0 - Width of search areas at first layer, has to be 8, 16, 32 or 64
//Remarks:
//  The seach areas are not centered at the corresponding searched blocks,
//  their positions are instead recursively estimated from filtered and halfed images.
//  Choose the size of search areas carefully in order to have fair computation time so as a fair motion estimation.
//  Small areas are best for this recursive algorithm (8x8 or 16x16).
//  SSE2 instuctions deal better for 16x16 search areas than MMX/SSE instructions.
//
//  {image:quad_motion}
//Example:
//  typedef ucharImage::index_type index_type;
//  RGBImage Image0=..., Image1=...;
//  ucharImage X=first(value_cast<YUV>(Image0)); // Isolates the greyscale component (Y)
//  ucharImage Y=first(value_cast<YUV>(Image1));
//  DenseMatrix<index_type>::self M = quad_motion<8,8>(X,Y,2,8,8,16,16); // stores all motion vectors
//See: ^motion^, ^neighbourhood_motion^
template<int M,int N,class G> inline typename QuadMotionArray<M,N,const Matrix<G> >::self 
quad_motion(const Matrix<G> &X,const Matrix<G> &Y,int l,int m,int n               ) { return quad_motion<M,N>(data(X),data(Y),l,typename Matrix<G>::size_type(m,n));                                      }
//{unsecret}
template<int M,int N,class G> inline typename QuadMotionArray<M,N,const Matrix<G> >::self 
quad_motion(const Matrix<G> &X,const Matrix<G> &Y,int l,int m,int n,int m0, int n0) { return quad_motion<M,N>(data(X),data(Y),l,typename Matrix<G>::size_type(m,n),typename Matrix<G>::size_type(m0,n0)); }


//{unsecret}
//Summary: Draws the motion estimation on an image
//Arguments: 
//  X - The image to draw on
//  E - The motion estimation
//  val - The color to draw with
template<class V,class G>
void motion_draw(Image<V> &X, const Matrix<G> &E, const V &val=V())
{
  _mm_empty();
  typedef typename Image<V>::size_type size_type;
  typedef typename Image<V>::index_type index_type;
  size_type d = X.size()/E.size();
  size_type d2 = d/2;
  for (int i=0, imax=E.nrows(); i<imax;++i)
    for (int j=0, jmax=E.ncols(); j<jmax;++j)
    {
      index_type p1 = index_type(i,j)*d+d2;
      index_type p2 = p1+E(i,j);
      draw_line(X,p1,p2,val);
    }
}

//{unsecret}
//Summary: Error (sum of absolute difference) with an estimated image.
//Arguments: 
//  M - Height of the blocks (has currently to be 8)
//  N - Width of the blocks (has currently to be 8)
//  X - The image to compare with.
//  Y - The image on which apply the motion estimation
//  E - The motion estimation
template <int M,int N,class G,class V>
long int motion_dist(const Image<V> &X, const Image<V> &Y, const Matrix<G> &E)
{
  assert(X.size()==Y.size());
  long int dist=0;
  for (int i=0, imax=E.nrows(); i<imax; ++i)
    for (int j=0, jmax=E.ncols(); j<jmax; ++j)
    {
      typename Matrix<G>::index_type p(M*i,N*j);
      dist+=block_dist<M,N>(sub<M,N>(X,p),sub<M,N>(Y,p+E(i,j)));
    }
  _mm_empty();
  return dist;
}

template <class G>
int vector_field_comparison (const Matrix<G> &M1, const Matrix<G> &M2)
{ 
  assert(M1.size()==M2.size());
  int s=0;
  for (int i=0, imax=M1.nrows(); i<imax; ++i)
    for (int j=0, jmax=M1.ncols(); j<jmax; ++j)
      if(M1(i,j)!=M2(i,j)) s+=1;
      
  return s;
}  

#ifdef MOTION_PRECOMPILE

GENIAL_API pair<Matrix<data_matrix_generator<unsigned char> >::index_type,short> motion(const Matrix<tiny_sub_dense_matrix_generator<8,8,Matrix<data_matrix_generator<unsigned char> > > > &X, const Matrix<sub_dense_matrix_generator<Matrix<data_matrix_generator<unsigned char> > > > &Y);

#endif //MOTION_PRECOMPILE

    
//} // namespace genial

#endif //__cplusplus


#ifdef __cplusplus
extern "C" {
#endif

#ifdef MOTION_ESTIMATION_PRECOMPILE

//Group = Motion Estimation C Interface

//{unsecret}
//Summary: Motion Estimation of two greyscale images with 8x8 blocks
//Arguments: 
//  mx - Height of the images
//  nx - Width of the images
//  px - Pointer to the first image
//  py - Pointer to the second image
//  m  - Height of search areas, has to be 8, 16, 32 or 64.
//  n  - Width of search areas, has to be 8, 16, 32 or 64.
//  p  - A matrix of indexes representing the motion estimation vectors of each 8x8 block of X.
//       Motion vectors are in the range [-m/2 ... m/2-1] for vertical shifting and [-n/2...n/2-1] for horizontal shifting.
//       Positive values mean displacement to the bottom, respectively to the right.
//       At least 2*(mx/8)*(nx/8) integers have to be allocated.
//Example:
//  ucharImage Image0, Image1;
//  PGMFile fin0("x.pgm"), fin1("y.pgm");
//  fin0>>Image0; fin1>>Image1;
//  DenseVector<ucharImage::index_type>::self M(X.size()/8);
//  motion8(X.nrows(),X.ncols(), &X(0,0), &Y(0,0), 32,32, (int *)&M(0,0));
//See: ^motion^
GENIAL_API void motion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int m,int n, int *p);

//{unsecret}
//Summary: Motion Estimation of two greyscale images with 8x8 blocks
//Arguments: 
//  mx - Height of the images
//  nx - Width of the images
//  px - Pointer to the first image
//  py - Pointer to the second image
//  m  - Height of search areas, has to be 8, 16, 32 or 64.
//  n  - Width of search areas, has to be 8, 16, 32 or 64.
//  p  - A matrix of indexes representing the motion estimation vectors of each 8x8 block of X.
//       Motion vectors are in the range [-m/2 ... m/2-1] for vertical shifting and [-n/2...n/2-1] for horizontal shifting.
//       Positive values mean displacement to the bottom, respectively to the right.
//       At least 2*(mx/8)*(nx/8) integers have to be allocated.
//Example:
//  ucharImage Image0, Image1;
//  PGMFile fin0("x.pgm"), fin1("y.pgm");
//  fin0>>Image0; fin1>>Image1;
//  DenseVector<ucharImage::index_type>::self M(X.size()/8);
//  nmotion8(X.nrows(),X.ncols(), &X(0,0), &Y(0,0), 16,16, (int *)&M(0,0));
//See: ^neighbourhood_motion^
GENIAL_API void nmotion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int m,int n, int *p);


//{unsecret}
//Summary: Motion Estimation of two greyscale images with 8x8 blocks
//Arguments: 
//  mx - Height of the images
//  nx - Width of the images
//  px - Pointer to the first image
//  py - Pointer to the second image
//  l  - Number of layers (image layer excluded)
//       The value 0 makes /quad_motion/ equivalent to the ^motion^ function.
//  m  - Height of search areas, has to be 8, 16, 32 or 64.
//  n  - Width of search areas, has to be 8, 16, 32 or 64.
//  m0 - Height of search areas at first layer, has to be 8, 16, 32 or 64
//  n0 - Width of search areas at first layer, has to be 8, 16, 32 or 64
//  p  - A matrix of indexes representing the motion estimation vectors of each 8x8 block of X.
//       Motion vectors are in the range [-m/2 ... m/2-1] for vertical shifting and [-n/2...n/2-1] for horizontal shifting.
//       Positive values mean displacement to the bottom, respectively to the right.
//       At least 2*(mx/8)*(nx/8) integers have to be allocated.
//Example:
//  ucharImage Image0, Image1;
//  PGMFile fin0("x.pgm"), fin1("y.pgm");
//  fin0>>Image0; fin1>>Image1;
//  DenseVector<ucharImage::index_type>::self M(X.size()/8);
//  qmotion8(X.nrows(),X.ncols(), &X(0,0), &Y(0,0), 8,8, 16,16, (int *)&M(0,0));
//See: ^quad_motion^
GENIAL_API void qmotion8(int mx,int nx, const unsigned char *px, const unsigned char *py, int l, int m,int n, int m0,int n0, int *p);


#endif //MOTION_ESTIMATION_PRECOMPILE

#ifdef __cplusplus
}
#endif


#endif














