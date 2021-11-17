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

#ifndef MATRIXSIZE_H
#define MATRIXSIZE_H

//namespace genial
//{

template<class I> class MatrixIndex;

template<class V1,class V2> struct promotion2_traits<MatrixSize<V1>, V2 > { typedef MatrixSize<PROMOTE2(V1,V2)> value_type; };
template<class V1,class V2> struct promotion2_traits<V1, MatrixSize<V2> > { typedef MatrixSize<PROMOTE2(V1,V2)> value_type; };
template<class V1,class V2> struct promotion2_traits<MatrixSize<V1>, MatrixSize<V2> > { typedef MatrixSize<PROMOTE2(V1,V2)> value_type; };

//{unsecret}
//{group:Matrices}
template<class I>
class MatrixSize
{
  public:
    typedef MatrixSize self;
    typedef I value_type;
    typedef value_type int_type;
    typedef MatrixIndex<int_type> index_type;

  public:
    union { value_type i; value_type m; value_type y; value_type first; };
    union { value_type j; value_type n; value_type x; value_type second; };
     
  public:
    inline MatrixSize() : m(0), n(0) {}
    inline MatrixSize(int_type iw, int_type ih) : m(iw), n(ih) {}
    inline MatrixSize(const MatrixIndex<I> &p) : m(p.i), n(p.j) {}

    inline MatrixSize(const self &r) : m(r.m), n(r.n) {}
    template<class S2> inline MatrixSize(const MatrixSize<S2> &r) : m(r.nrows()), n(r.ncols()) {}

    inline operator complex<float>() const { return complex<float>(ncols(),nrows()); }
      
    inline int_type nelms() const { return m*n; }
    inline int_type nrows() const { return m; }
    inline int_type ncols() const { return n; }
    inline int_type height() const { return m; }
    inline int_type width () const { return n; }

    inline void set(const self &r) { m=r.m; n=r.n; }
    inline void set(int_type i, int_type j) { m=i; n=j; }
    inline void setnrows(int_type i) { m=i; }
    inline void setncols(int_type j) { n=j; }

    inline bool operator==(const self &r) const { return (m==r.m) && (n==r.n); }
    inline bool operator!=(const self &r) const { return (m!=r.m) || (n!=r.n); }
    inline bool null() const { return nelms()==0; }

    inline bool withinbounds  (const index_type &p) const { return p.i>=0 && p.j>=0 && p.i<nrows() && p.j<ncols(); }
    inline bool leftboundary  (const index_type &p) const { return p.j==0; }
    inline bool rightboundary (const index_type &p) const { return p.j==ncols()-1; }
    inline bool topboundary   (const index_type &p) const { return p.i==0; }
    inline bool bottomboundary(const index_type &p) const { return p.i==nrows()-1; }
    inline bool boundary      (const index_type &p) const { return leftboundary(p) || rightboundary(p) ||	topboundary(p) || bottomboundary(p); }

    inline index_type upper_bound() const { return index_type(nrows()-1, ncols()-1); }
     
    inline void transpose() { int_type t=m; m=n; n=t; }
};

//Group=Arrays functions

template<class E,class Tr,class I> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const MatrixSize<I>& d) { return os << '(' << d.nrows() << "x" << d.ncols() << ')'; }
template<class E,class Tr,class I> inline basic_istream<E,Tr> &operator>>(basic_istream<E,Tr> &is,       MatrixSize<I>& d) { E c; I m=0,n=0; is>>c && c=='(' && is>>m>>c && c=='x' && is>>n>>c && c==')'; d.set(m,n); return is; }

template<class I> inline MatrixSize<I> operator+(const MatrixSize<I> &p1, const MatrixSize<I> &p2) { return MatrixSize<I>(p1.nrows()+p2.nrows(),p1.ncols()+p2.ncols()); }
template<class I> inline MatrixSize<I> operator+(const MatrixSize<I> &p, const I &x) { return MatrixSize<I>(p.nrows()+x,p.ncols()+x); }
template<class I> inline MatrixSize<I> operator+(const I &x, const MatrixSize<I> &p) { return MatrixSize<I>(x+p.nrows(),x+p.ncols()); }

template<class I> inline MatrixSize<PROMOTE2(I,int)> operator-(const MatrixSize<I> &p1, const MatrixSize<I> &p2) { return MatrixSize<PROMOTE2(I,int)>(p1.nrows()-p2.nrows(),p1.ncols()-p2.ncols()); }
template<class I> inline MatrixSize<PROMOTE2(I,int)> operator-(const MatrixSize<I> &p, const I &x) { return MatrixSize<PROMOTE2(I,int)>(p.nrows()-x,p.ncols()-x); }
template<class I> inline MatrixSize<PROMOTE2(I,int)> operator-(const I &x, const MatrixSize<I> &p) { return MatrixSize<PROMOTE2(I,int)>(x-p.nrows(),x-p.ncols()); }

template<class I> inline MatrixSize<I> operator*(const I &x, const MatrixSize<I> &s) { return MatrixSize<I>(x*s.nrows(),x*s.ncols()); }
template<class I> inline MatrixSize<I> operator*(const MatrixSize<I> &s, const I &x) { return MatrixSize<I>(s.nrows()*x,s.ncols()*x); }

template<class I,class J> inline MatrixSize<PROMOTE2(I,J)> operator*(const MatrixSize<I> &x, const MatrixSize<J> &y) { return MatrixSize<PROMOTE2(I,J)>(x.nrows()*y.nrows(), x.ncols()*y.ncols()); }

template<class I> inline MatrixSize<I> operator/(const MatrixSize<I> &d, const I &x) { return MatrixSize<I>(d.nrows()/x,d.ncols()/x); }
template<class I> inline MatrixSize<I> operator/(const MatrixSize<I> &d1, const MatrixSize<I> &d2) { return MatrixSize<I>(d1.nrows()/d2.nrows(), d1.ncols()/d2.ncols()); }

template<class I> inline double norm(const MatrixSize<I> &s) { return sqr(double(s.nrows()))+sqr(double(s.ncols())); }
template<class I> inline double abs (const MatrixSize<I> &s) { return sqrt(norm(s)); }

template<class I> inline MatrixSize<I> trn(const MatrixSize<I> &s) { return MatrixSize<I>(s.ncols(),s.nrows()); }


template<class I> inline MatrixSize<I> operator-(const MatrixSize<I> &p1, const MatrixIndex<I> &p2) { return MatrixSize<I>(p1.nrows()-p2.i,p1.ncols()-p2.j); }

//}

#endif

