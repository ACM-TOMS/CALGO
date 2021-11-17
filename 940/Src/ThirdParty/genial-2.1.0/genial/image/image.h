//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
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

#ifndef IMAGE_H
#define IMAGE_H

#include "array/matrix.h"
#include "color.h"

//namespace genial
//{

//using namespace genial;


//{unsecret}
//{group:Image Overview}
//Summary: Image container
//Remark: Use it like a dense matrix
template<class V> 
class Image : public Matrix<dense_matrix_generator<V> >
{
  public:
    typedef Image self;
    typedef Matrix<dense_matrix_generator<V> > base;
    
  public:
    Image() : base() {}
    template<class A> Image(const A &a) : base(a) {}
    template<class A,class B> Image(const A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> Image(const A &a, const B &b, const C &c) : base (a,b,c) {}
  
    template<class T> self &operator=(const T &Y) { static_cast<base &>(*this)=Y; return *this; }
};

template<class V1,class V2> struct promotion2_traits<Image<V1>, V2        > { typedef typename Image<PROMOTE2(typename Image<V1>::const_value_type,V2)>::self value_type; };
template<class V1,class V2> struct promotion2_traits<V1,        Image<V2> > { typedef typename Image<PROMOTE2(V1,typename Image<V2>::const_value_type)>::self value_type; };
template<class V1,class V2> struct promotion2_traits<Image<V1>, Image<V2> > { typedef typename Image<PROMOTE2(typename Image<V1>::const_value_type, typename Image<V2>::const_value_type)>::self value_type; };

//Group = Images Type Definitions

//{unsecret}
typedef Image<RGBColor<unsigned char> > RGBImage;
//{unsecret}
typedef Image<RGBA>           RGBAImage;
//{unsecret}
typedef Image<YUV>            YUVImage;
//{unsecret}
typedef Image<YUVA>           YUVAImage;
//{unsecret}
typedef Image<YCbCr>          YCbCrImage;
//{unsecret}
typedef Image<unsigned char>  ucharImage;
//{unsecret}
typedef Image<unsigned short> ushortImage;
//{unsecret}
typedef Image<unsigned char>  charImage;
//{unsecret}
typedef Image<unsigned short> shortImage;
//{unsecret}
typedef Image<float>          floatImage;
//{unsecret}
typedef Image<double>         doubleImage;
//{unsecret}
typedef Image<bool>           boolImage;



//}

#endif

