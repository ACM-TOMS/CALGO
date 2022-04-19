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

#ifndef MATRIXINDEX_H
#define MATRIXINDEX_H

#include <iostream>

//namespace genial
//{

template<class I> class MatrixIndex;
template<class I> class MatrixSize;

template<class V1,class V2> struct promotion2_traits<MatrixIndex<V1>, V2 > { typedef MatrixIndex<PROMOTE2(V1,V2)> value_type; };
template<class V1,class V2> struct promotion2_traits<V1, MatrixIndex<V2> > { typedef MatrixIndex<PROMOTE2(V1,V2)> value_type; };
template<class V1,class V2> struct promotion2_traits<MatrixIndex<V1>, MatrixIndex<V2> > { typedef MatrixIndex<PROMOTE2(V1,V2)> value_type; };

//{unsecret}
//{group:Matrices}
template<class I>
class MatrixIndex
{
  public:
    typedef MatrixIndex self;
    typedef I value_type;
    typedef value_type int_type;
    typedef value_type first_type;
    typedef value_type second_type;

    template<class U> struct rebind { typedef MatrixIndex<U> other; };

  public:
    union { value_type i; value_type m; value_type y; value_type first; };
    union { value_type j; value_type n; value_type x; value_type second; };

  public:
    inline MatrixIndex() : i(), j() { }
    inline explicit MatrixIndex(const value_type &m) : i(m), j(m) {}
    inline MatrixIndex(const value_type &m, const value_type &n) : i(m), j(n) {}
    template<class T> inline explicit MatrixIndex(const complex<T> &z) : y(-z.imag()), x(z.real()) { }

    inline MatrixIndex(const self &r) : i(r.i), j(r.j) { }
    template<class I2> inline MatrixIndex(const MatrixIndex<I2> &x) : i(value_type(x.i)), j(value_type(x.j)) {}

    inline MatrixIndex &operator=(const self &x) { i=x.i; j=x.j; return *this; }
    template<class J> inline MatrixIndex &operator=(const MatrixIndex<J> &x) { i=x.i; j=x.j; return *this; }
    template<class T> inline MatrixIndex &operator=(const complex    <T> &z) { y=-z.imag(); x=z.real(); return *this; }

    //inline operator complex<float>() const { return complex<float>(x,y); } //pas sûr, mais bon pour fourier_function
    //inline operator complex<float>() const { return complex<float>(x,-y); }

    inline void set(value_type m=0, value_type n=0) { i=m; j=n; }

    inline bool operator==(const self &p) const { return i==p.i && j==p.j; }
    inline bool operator!=(const self &p) const { return i!=p.i || j!=p.j; }
};

//Group=Arrays functions

template<class E,class T,class I> inline basic_ostream<E,T>& operator<<(basic_ostream<E,T> &os, const MatrixIndex<I>& p) { return os << "(" << p.i << "," << p.j << ")"; }
template<class E,class T,class I> inline basic_istream<E,T>& operator>>(basic_istream<E,T> &is,       MatrixIndex<I>& p) { E c; is>>c && c=='(' && is>>p.i>>c && c==',' && is>>p.j>>c && c==')'; return is; }

template<class V> inline V       &first (      MatrixIndex<V> &x) { return x.first; }
template<class V> inline const V &first (const MatrixIndex<V> &x) { return x.first; }
template<class V> inline V       &second(      MatrixIndex<V> &x) { return x.second; }
template<class V> inline const V &second(const MatrixIndex<V> &x) { return x.second; }


template<class I> inline MatrixIndex<I> operator-(const MatrixIndex<I> &p) { return MatrixIndex<I>(-p.i,-p.j); }

template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator+(const MatrixIndex<I> &p1, const MatrixIndex<J> &p2) { return MatrixIndex<PROMOTE2(I,J)>(p1.i+p2.i,p1.j+p2.j); }
template<class I,class J> inline MatrixIndex<PROMOTE3(int,I,J)> operator-(const MatrixIndex<I> &p1, const MatrixIndex<J> &p2) { return MatrixIndex<PROMOTE3(int,I,J)>(p1.i-p2.i,p1.j-p2.j); }

template<class I> inline MatrixIndex<I> operator+(const MatrixIndex<I> &p, const I &x) { return MatrixIndex<I>(p.i+x,p.j+x); }
template<class I> inline MatrixIndex<I> operator-(const MatrixIndex<I> &p, const I &x) { return MatrixIndex<I>(p.i-x,p.j-x); }

template<class I> inline MatrixIndex<I> operator+(const I &x, const MatrixIndex<I> &p) { return MatrixIndex<I>(x+p.i,x+p.j); }
template<class I> inline MatrixIndex<I> operator-(const I &x, const MatrixIndex<I> &p) { return MatrixIndex<I>(x-p.i,x-p.j); }

template<class I> inline MatrixIndex<I> operator+(const MatrixIndex<I> &p, const MatrixSize <I> &s) { return MatrixIndex<I>(p.i+s.nrows(),p.j+s.ncols()); }
template<class I> inline MatrixIndex<I> operator-(const MatrixIndex<I> &p, const MatrixSize <I> &s) { return MatrixIndex<I>(p.i-s.nrows(),p.j-s.ncols()); }

template<class I,class J> inline MatrixIndex<I> &operator+=(MatrixIndex<I> &p1, const MatrixIndex<J> &p2) { p1.i+=p2.i; p1.j+=p2.j; return p1; }
template<class I,class J> inline MatrixIndex<I> &operator-=(MatrixIndex<I> &p1, const MatrixIndex<J> &p2) { p1.i-=p2.i; p1.j-=p2.j; return p1; }

template<class I> inline MatrixIndex<I> &operator+=(MatrixIndex<I> &p, const I &i) { p.x-=i; p.y-=i; return *p; }
template<class I> inline MatrixIndex<I> &operator-=(MatrixIndex<I> &p, const I &i) { p.x-=i; p.y-=i; return *p; }

template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator*(const MatrixIndex<I> &x, const MatrixIndex<J> &y) { return MatrixIndex<PROMOTE2(I,J)>(x.i*y.i, x.j*y.j); }
template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator*(const MatrixIndex<I> &p, const MatrixSize <J> &d) { return MatrixIndex<PROMOTE2(I,J)>(p.i*d.nrows(), p.j*d.ncols()); }
template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator*(const MatrixSize <I> &d, const MatrixIndex<J> &p) { return MatrixIndex<PROMOTE2(I,J)>(d.nrows()*p.i, d.ncols()*p.j); }

template<class I> inline MatrixIndex<I> operator*(const MatrixIndex<I> &p, const I &x) { return MatrixIndex<I>(p.i*x,p.j*x); }
template<class I> inline MatrixIndex<I> operator*(const I &x, const MatrixIndex<I> &p) { return MatrixIndex<I>(x*p.i,x*p.j); }

template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator/(const MatrixIndex<I> &p, const J &x) { return MatrixIndex<PROMOTE2(I,J)>(p.i/x,p.j/x); }
template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator/(const MatrixIndex<I> &x, const MatrixIndex<J> &y) { return MatrixIndex<PROMOTE2(I,J)>(x.i/y.i,x.j/y.j); }
template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> operator/(const MatrixIndex<I> &x, const MatrixSize <J> &y) { return MatrixIndex<PROMOTE2(I,J)>(x.i/y.nrows(),x.j/y.ncols()); }
template<class I,class J> inline MatrixIndex<PROMOTE2(I,J)> &operator/=(MatrixIndex<I> &p, const J &x) { p.i/=x; p.j/=x; return p; }

template<class I> inline bool operator>=(const MatrixIndex<I> &x, const I &n) { return x.i>=n && x.j>=n; }
template<class I> inline bool operator< (const MatrixIndex<I> &x, const MatrixIndex<I> &y) { return x.i<y.i       || (x.i==y.i       && x.j<y.j      ); }
template<class I> inline bool operator< (const MatrixIndex<I> &x, const MatrixSize <I> &y) { return x.i<y.nrows() || (x.i==y.nrows() && x.j<y.ncols()); }

template<class I> inline MatrixIndex<I> max(const MatrixIndex<I> &x, const MatrixIndex<I> &y) { return MatrixIndex<I>(max(x.i,y.i),max(x.j,y.j)); }
template<class I> inline MatrixIndex<I> min(const MatrixIndex<I> &x, const MatrixIndex<I> &y) { return MatrixIndex<I>(min(x.i,y.i),min(x.j,y.j)); }

template<class I> inline PROMOTE2(float,I) norm(const MatrixIndex<I> &p) { return sqr(PROMOTE2(float,I)(p.i))+sqr(PROMOTE2(float,I)(p.j)); }
template<class I> inline PROMOTE2(float,I) abs (const MatrixIndex<I> &p) { return sqrt(norm(p)); }
template<class I,class J> inline PROMOTE3(float,I,J) dist(const MatrixIndex<I> &p1,const MatrixIndex<J> &p2) { return abs(p2-p1); }

template<class I> inline MatrixIndex<I> round(const MatrixIndex<I> &p) { return MatrixIndex<I>(round(p.i),round(p.j)); }

template<class I> inline MatrixIndex<I> trn(const MatrixIndex<I> &p) { return MatrixIndex<I>(p.j,p.i); }

template<class I> inline MatrixIndex<I> mod(const MatrixIndex<I> &p, const MatrixSize<I> &s) { return MatrixIndex<I>( mod(p.i,s.nrows()), mod(p.j,s.ncols()) ); }

template<class I> inline MatrixIndex<I> median(MatrixIndex<I> p1, MatrixIndex<I> p2, MatrixIndex<I> p3)
{
  if (p1.i>p2.i) swap(p1.i,p1.i); if (p2.i>p3.i) swap(p2.i,p3.i);
  if (p1.j>p2.j) swap(p1.j,p1.j); if (p2.j>p3.j) swap(p2.j,p3.j);
  return max(p1,p2);
}

//}


#endif

