//-*- C++ -*-
//
//  Dense Vector and Matrix
//
//     Copyright (C) 2000 Masakatsu Ito
//     Institute for Molecular Science
//     Myodaiji  Okazaki  444-0867 Japan
//     E-mail mito@ims.ac.jp

//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 

#ifndef VECMAT_H_
#define VECMAT_H_

#include <iostream>
#include <iomanip>
#include <fstream>

#include <assert.h>

#include "mathfunc.h"

#include "xpr1.h"
#include "xprD.h"
#include "xpr2.h"

#ifdef NUMRECIPES
#define FORTARRAY
#endif

#ifdef NUMRECIPES

#ifdef UNDERF77
#define ZTRED2 ztred2_
#define ZTQLI ztqli_
#define DTRED2 dtred2_
#define DTQLI dtqli_
#define TRED2 tred2_
#define TQLI tqli_
#else
#define ZTRED2 ztred2
#define ZTQLI ztqli
#define DTRED2 dtred2
#define DTQLI dtqli
#define TRED2 tred2
#define TQLI tqli
#endif

extern "C" {
  void ZTRED2(Complex*,const int*,const int*,Complex*,Complex*);
  void ZTQLI(Complex*,Complex*,const int*, const int*, Complex*);
  void DTRED2(double*,const int*,const int*,double*,double*);
  void DTQLI(double*,double*,const int*, const int*, double*);
  void TRED2(float*,const int*,const int*,float*,float*);
  void TQLI(float*,float*,const int*, const int*, float*);
}

#endif

static const double EPS = double(1.0e-75);

template <class T> class RangeVec;

template <class T>
class Vec : public Dim1<T,Vec<T> > {
private:
  int sz;
  T* v;

  void error(const char *msg) const {
      std::cerr << "Vec error: " << msg << std::endl;
    exit(1);
  }

public:
  explicit Vec(int sz_) : Dim1<T,Vec<T> >(), sz(sz_), v(new T[sz]) {}
  Vec(const Vec<T>& x) : Dim1<T,Vec<T> >(), sz(x.sz), v(new T[sz]) {
    assignFrom(x);
  }
  Vec() : Dim1<T,Vec<T> >() {}
  void init(int sz_) {
    sz = sz_;
    v = new T[sz];
  }
  ~Vec() { delete [] v; }
  template <class X> Vec<T>& operator=(const Xpr1<T,X>& rhs) {
    return this->assignFrom(rhs);
  }
  template <class V> Vec<T>& operator=(const Dim1<T,V>& rhs) {
    return this->assignFrom(rhs);
  }
  Vec<T>& operator=(const Vec<T>& rhs) { return assignFrom(rhs); }
  Vec<T>& operator=(T rhs) { return assignFrom(rhs); }
  //template <int N,class V> Vec<T>& operator=(const TDim1<N,T,V>& rhs) {
  //  return assignFrom(rhs);
  //}
  template <class Closure> Vec<T>& operator=(const Closure& rhs) {
    rhs.assignTo(v);
    return *this;
  }
//   Vec(const Vec<T>& x) : Dim1<T,Vec<T> >(), sz(x.sz), v(new T[sz]) {
//     (*this) = x;
//   }
  int size() const { return sz; }
  T& operator()(int i) { 
#ifdef RANGECHECK
    if ( i < 0 || sz <= i ) error("out of range");
#endif
    return v[i]; 
  }
  T operator()(int i) const { 
#ifdef RANGECHECK
    if ( i < 0 || sz <= i ) error("out of range");
#endif
    return v[i]; 
  }
  RangeVec<T> operator()(Range sl) const;

  double norm() const { 
    double n = v[0] * v[0];
    for (int i=1; i<sz; i++) n += v[i]*v[i];
    return n;
  }
  double abs() const { return sqrt( norm() ); }

  const T* datPtr() { return v; }
  void write(std::ofstream& str) const { str.write((char*)v, sizeof(T)*sz); }
  void read(std::ifstream& str) { 
    str.read((char*)v, sizeof(T)*sz);
    if (str.bad()) error("Bad in reading.");
    if (str.fail()) error("Failed in reading.");
    if (str.eof()) error("EOF in reading.");
  }
  long seekSize() const { return sizeof(T)*sz; }
};

template <class T>
class RangeVec : public Dim1<T,RangeVec<T> > {
  const Range sl;
  T* const dat;

  void error(const char *msg) const {
    std::cerr << "RangeVec error: " << msg << std::endl;
    exit(1);
  }

public:
  RangeVec(Range sl_,T* dat_) : Dim1<T,RangeVec<T> >(), sl(sl_), dat(dat_) {}
  template <class X> RangeVec<T>& operator=(const Xpr1<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class V> RangeVec<T>& operator=(const Dim1<T,V>& rhs) {
    return assignFrom(rhs);
  }
  //template <int N,class V> RangeVec<T>& operator=(const TDim1<N,T,V>& rhs) {
  //  return assignFrom(rhs);
  //}
  RangeVec<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return sl.size(); }
  T& operator()(int i) { 
#ifdef RANGECHECK
    if ( i < 0 || sl.size() <= i ) error("out of range");
#endif
    return dat[sl(i)]; 
  }
  T operator()(int i) const { 
#ifdef RANGECHECK
    if ( i < 0 || sl.size() <= i ) error("out of range");
#endif
    return dat[sl(i)]; 
  }
};

template <class T>
RangeVec<T> Vec<T>::operator()(Range sl) const {
  return RangeVec<T>(sl,v);
}

template <class T> class RangeRecMat;

template <class T>
class RecMat : public Dim2<T,RecMat<T> > {
private:
  int rsz, csz, tsz;
  T* const dat;
  T** m;

  void error(const char *msg) const {
    std::cerr << "RecMat error: " << msg << std::endl;
    exit(1);
  }

public:
  RecMat(int rsz_, int csz_) : Dim2<T,RecMat<T> >(),
    rsz(rsz_), csz(csz_), tsz(rsz*csz), dat(new T[tsz]), 
#ifdef FORTARRAY
    m(new T*[csz])
#else
    m(new T*[rsz])
#endif
    {
#ifdef FORTARRAY
      for (int i=0; i<csz; i++) m[i] = dat + i*rsz;
#else
      for (int i=0; i<rsz; i++) m[i] = dat + i*csz;
#endif
    }
  RecMat(const RecMat<T>& x) : Dim2<T,RecMat<T> >(),
    rsz(x.rsz), csz(x.csz), tsz(rsz*csz), dat(new T[tsz]), 
#ifdef FORTARRAY
    m(new T*[csz])
#else
    m(new T*[rsz])
#endif
    {
#ifdef FORTARRAY
      for (int i=0; i<csz; i++) m[i] = dat + i*rsz;
#else
      for (int i=0; i<rsz; i++) m[i] = dat + i*csz;
#endif
      assignFrom(x);
    }
  RecMat() : Dim2<T,RecMat<T> >() {}
  void init(int rsz_, int csz_) {
    rsz = rsz_; csz = csz_; tsz = rsz*csz;
    dat = new T[tsz];
#ifdef FORTARRAY
    m = new T*[csz];
    for (int i=0; i<csz; i++) m[i] = dat + i*rsz;
#else
    m = new T*[rsz];
    for (int i=0; i<rsz; i++) m[i] = dat + i*csz;
#endif
  }
  ~RecMat() {
    delete [] m;
    delete [] dat;
  }
  RecMat<T>& operator=(const RecMat<T>& rhs) { return assignFrom(rhs); }
  template <class X> RecMat<T>& operator=(const Xpr2<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> RecMat<T>& operator=(const Dim2<T,M>& rhs) {
    return assignFrom(rhs);
  }
  //template <int Nr,int Nc,class M>
  //  RecMat<T>& operator=(const TDim2<Nr,Nc,T,M>& rhs) { return assignFrom(rhs); }
  RecMat<T>& operator=(T rhs) { return assignFrom(rhs); }
  template <class Closure> RecMat<T>& operator=(const Closure& rhs) {
    rhs.assignTo(m);
    return *this;
  }
  int size() const { return tsz; }
  int rows() const { return rsz; }
  int cols() const { return csz; }
  //T& operator()(int n) { return dat[n]; }
  //T operator()(int n) const { return dat[n]; }
  T& operator()(int i, int j) { 
#ifdef RANGECHECK
    if ( i < 0 || rsz <= i ) error("out of row range");
    if ( j < 0 || csz <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[j][i]; 
#else
    return m[i][j]; 
#endif
  }
  T operator()(int i, int j) const { 
#ifdef RANGECHECK
    if ( i < 0 || rsz <= i ) error("out of row range");
    if ( j < 0 || csz <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[j][i]; 
#else
    return m[i][j]; 
#endif
  }
  //RangeVec<T> operator()(RowRange rsl, int j) const {
  //  return RangeVec<T>(Range(rsl.start()*csz+j,rsl.size()*rsl.step()*csz+j,
		//	     rsl.step()*csz),
		//       v);
  //}
  //RangeVec<T> operator()(int i, ColRange csl) const {
  //  return RangeVec<T>(Range(i*csz+csl.start(),i*csz+csl.size()*csl.step(),
		//	     csl.step()),
		//       v);
  //}
  RangeRecMat<T> operator()(RowRange rsl, ColRange csl) const;
};

template <class T>
class RangeRecMat : public Dim2<T,RangeRecMat<T> > {
private:
  const RowRange rsl;
  const ColRange csl;
  T** const  m;

  void error(const char *msg) const {
    std::cerr << "RangeRecMat error: " << msg << std::endl;
    exit(1);
  }

public:
  RangeRecMat(RowRange rsl_, ColRange csl_, T** m_) :
  Dim2<T,RangeRecMat<T> >(), rsl(rsl_), csl(csl_), m(m_) {}
  template <class X> RangeRecMat<T>& operator=(const Xpr2<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> RangeRecMat<T>& operator=(const Dim2<T,M>& rhs) {
    return assignFrom(rhs);
  }
  //template <int Nr,int Nc,class M>
  //  RangeRecMat<T>& operator=(const TDim2<Nr,Nc,T,M>& rhs) {
  //    return assignFrom(rhs);
  //  }
  RangeRecMat<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return rsl.rows() * csl.cols(); }
  int rows() const { return rsl.rows(); }
  int cols() const { return csl.cols(); }
  T& operator()(int i, int j) { 
#ifdef RANGECHECK
    if ( i < 0 || rsl.rows() <= i ) error("out of row range");
    if ( j < 0 || csl.cols() <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[csl(j)][rsl(i)]; 
#else
    return m[rsl(i)][csl(j)];
#endif
  }
  T operator()(int i, int j) const { 
#ifdef RANGECHECK
    if ( i < 0 || rsl.rows() <= i ) error("out of row range");
    if ( j < 0 || csl.cols() <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[csl(j)][rsl(i)]; 
#else
    return m[rsl(i)][csl(j)];
#endif
  }
};

template <class T>
RangeRecMat<T> RecMat<T>::operator()(RowRange rsl, ColRange csl) const {
  return RangeRecMat<T>(rsl,csl,m);
};


template <class T> class Mat;
template <class T> class DiagRangeMat;

template <class T>
class DiagMat : public DimD<T,DiagMat<T> > {
  friend class Mat<T>;

private:
  int sz;
  T* dat;

  void error(const char *msg) const {
    std::cerr << "DiagMat error: " << msg << std::endl;
    exit(1);
  }

public:
  explicit DiagMat(int sz_) : DimD<T,DiagMat<T> >(), sz(sz_), dat(new T[sz]) {}
  DiagMat(const DiagMat<T>& x) : DimD<T,DiagMat<T> >(), 
    sz(x.sz), dat(new T[sz]) {
      assignFrom(x);
  }
  DiagMat() : DimD<T,DiagMat<T> >() {}
  void init(int sz_) {
    sz = sz_;
    dat = new T[sz];
  }
  ~DiagMat() {
    delete [] dat;
  }
  template <class X> DiagMat<T>& operator=(const XprD<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> DiagMat<T>& operator=(const DimD<T,M>& rhs) {
    return assignFrom(rhs);
  }
  DiagMat<T>& operator=(const DiagMat<T>& rhs) { return assignFrom(rhs); }
  //template <int Nr,int Nc,class M>
  //  DiagMat<T>& operator=(const TDimD<Nr,Nc,T,M>& rhs) { return assignFrom(rhs); }
  DiagMat<T>& operator=(T rhs) { return assignFrom(rhs); }
  template <class Closure> DiagMat<T>& operator=(const Closure& rhs) {
    rhs.assignTo(this);
    return *this;
  }
  int size() const { return sz; }
  int rows() const { return sz; }
  int cols() const { return sz; }
  T& operator()(int n) { 
#ifdef RANGECHECK
    if ( n < 0 || sz <= n ) error("out of range");
#endif
    return dat[n]; 
  }
  T operator()(int n) const { 
#ifdef RANGECHECK
    if ( n < 0 || sz <= n ) error("out of range");
#endif
    return dat[n]; 
  }
  T operator()(int i, int j) const { 
#ifdef RANGECHECK
    if ( i < 0 || sz <= i ) error("out of row range");
    if ( j < 0 || sz <= j ) error("out of column range");
#endif
    return ( i==j ? dat[i] : (T)0 ); 
  }
  DiagRangeMat<T> operator()(DiagRange sl) const;

  const T* datPtr() { return dat; }
  void write(std::ofstream& str) const { str.write((char*)dat, sizeof(T)*sz); }
  void read(std::ifstream& str) { 
    str.read((char*)dat, sizeof(T)*sz);
    if (str.bad()) error("Bad in reading.");
    if (str.fail()) error("Failed in reading.");
    if (str.eof()) error("EOF in reading.");
  }
  long seekSize() const { return sizeof(T)*sz; }
};

template <class T>
class DiagRangeMat : public DimD<T,DiagRangeMat<T> > {
private:
  const DiagRange sl;
  T* const dat;
public:
  DiagRangeMat(DiagRange sl_, T* dat_) : 
    DimD<T,DiagRangeMat<T> >(), sl(sl_), dat(dat_) {}
  template <class X> DiagRangeMat<T>& operator=(const XprD<T,X>& rhs) {
    return assignFrom(rhs);
  }
  template <class M> DiagRangeMat<T>& operator=(const DimD<T,M>& rhs) {
    return assignFrom(rhs);
  }
  //template <int Nr,int Nc,class M>
  //  DiagRangeMat<T>& operator=(const TDimD<Nr,Nc,T,M>& rhs) {
  //    return assignFrom(rhs);
  //  }
  DiagRangeMat<T>& operator=(T rhs) { return assignFrom(rhs); }
  int size() const { return sl.size(); }
  int rows() const { return sl.size(); }
  int cols() const { return sl.size(); }
  T& operator()(int i) { 
#ifdef RANGECHECK
    if ( i < 0 || sl.size() <= i ) error("out of range");
#endif
    return dat[sl(i)]; 
  }
  T operator()(int i) const { 
#ifdef RANGECHECK
    if ( i < 0 || sl.size() <= i ) error("out of range");
#endif
    return dat[sl(i)]; 
  }
  T operator()(int i, int j) const {
#ifdef RANGECHECK
    if ( i < 0 || sl.size() <= i ) error("out of row range");
    if ( j < 0 || sl.size() <= j ) error("out of column range");
#endif
    return ( i==j ? dat[sl(i)] : (T)0 ); 
  }
};

template <class T>
DiagRangeMat<T> DiagMat<T>::operator()(DiagRange sl) const {
  return DiagRangeMat<T>(sl,dat);
};


template <class T>
class Mat : public Dim2<T,Mat<T> > {
private:
  int sz, tsz;
  T* dat;
  T** m;
  
  void error(const char *msg) const {
    std::cerr << "Mat error: " << msg << std::endl;
    exit(1);
  }

public:
  explicit Mat(int sz_) : Dim2<T,Mat<T> >(), sz(sz_), tsz(sz*sz),
    dat(new T[tsz]), m(new T*[sz]) {
      for (int i=0; i<sz; i++) m[i] = dat + i*sz;
    }
  Mat(const Mat<T>& x) : Dim2<T,Mat<T> >(), sz(x.sz), tsz(sz*sz),
    dat(new T[tsz]), m(new T*[sz]) {
      for (int i=0; i<sz; i++) m[i] = dat + i*sz;
      this->assignFrom(x);
    }
  Mat() : Dim2<T,Mat<T> >() {}
  void init(int sz_) {
    sz = sz_; tsz = sz*sz;
    dat = new T[tsz];
    m = new T*[sz];
    for (int i=0; i<sz; i++) m[i] = dat + i*sz;
  }
  ~Mat() {
    delete [] m;
    delete [] dat;
  }
  Mat<T>& operator=(const Mat<T>& rhs) { return assignFrom(rhs); }
  template <class X> Mat<T>& operator=(const Xpr2<T,X>& rhs) {
    return this->assignFrom(rhs);
  }
  template <class M> Mat<T>& operator=(const Dim2<T,M>& rhs) {
    return assignFrom(rhs);
  }
  //template <int Nr,int Nc,class M>
  //  Mat<T>& operator=(const TDim2<Nr,Nc,T,M>& rhs) { return assignFrom(rhs); }
  Mat<T>& operator=(T rhs) { return assignFrom(rhs); }
  template <class Closure> Mat<T>& operator=(const Closure& rhs) {
    rhs.assignTo(m);
    return *this;
  }
  int size() const { return tsz; }
  int rows() const { return sz; }
  int cols() const { return sz; }
  T& operator()(int i, int j) { 
#ifdef RANGECHECK
    if ( i < 0 || sz <= i ) error("out of row range");
    if ( j < 0 || sz <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[j][i]; 
#else
    return m[i][j]; 
#endif
  }
  T operator()(int i, int j) const { 
#ifdef RANGECHECK
    if ( i < 0 || sz <= i ) error("out of row range");
    if ( j < 0 || sz <= j ) error("out of column range");
#endif
#ifdef FORTARRAY
    return m[j][i]; 
#else
    return m[i][j]; 
#endif
  }
  RangeVec<T> operator()(RowRange rsl, int j) const {
    return RangeVec<T>(Range(rsl.start()*sz+j,rsl.size()*rsl.step()*sz+j,
			     rsl.step()*sz),
		       dat);
  }
  RangeVec<T> operator()(int i, ColRange csl) const {
    return RangeVec<T>(Range(i*sz+csl.start(),i*sz+csl.size()*csl.step(),
			     csl.step()),
		       dat);
  }
  RangeRecMat<T> operator()(RowRange rsl, ColRange csl) const;

#ifdef FORTARRAY
  void luDecompose(int* ip) {
    int k;
    for (k=0; k<sz; k++)
      ip[k] = k;

    for (k=0; k<sz; k++) {
      int i,j;
      int mi = k;
      double max = Abs(m[k][ip[mi]]);
      for (i=k+1; i<sz; i++) {
	if (Abs(m[k][ip[i]]) > max) {
	  mi = i;
	  max = Abs(m[k][ip[mi]]);
	}
      }
    
      if (mi != k) {
	int temp = ip[k];
	ip[k] = ip[mi];
	ip[mi] = temp;
      }
    
      if (Abs(m[k][ip[k]]) < EPS) error("singular !");
    
      // Gauss elimination
      m[k][ip[k]] = T(1.0)/m[k][ip[k]];

      for (i=k+1; i<sz; i++) m[k][ip[i]] *= m[k][ip[k]];
      for (j=k+1; j<sz; j++) {
	for (i=k+1; i<sz; i++)
	  m[j][ip[i]] -= m[k][ip[i]]*m[j][ip[k]];
      }
    }
  }
  
  template <class V> void luSubst(const int* ip, 
				  const Dim1<T,V>& b, Vec<T>& ans) const {
    int i,j;
    for (i=0; i<sz; i++) ans(i) = b(ip[i]);
    
    for (j=0; j<sz; j++) {
      for (i=j+1; i<sz; i++) ans(i) -= m[j][ip[i]]*ans(j);
    }
    
    ans(sz-1) *= m[sz-1][ip[sz-1]];
    for (j=sz-1; j>=1; j--) {
      for (i=0; i<=j-1; i++) ans(i) -= m[j][ip[i]]*ans(j);
      ans(j-1) *= m[j-1][ip[j-1]];
    }
  }

#endif

  void house(DiagMat<T>& d, T* e);
  void triQL(DiagMat<T>& d, T* e);

  void diagonalize(Mat<T>& u, DiagMat<T>& d, T* e) const {
    u = (*this);
    u.house(d,e);
    u.triQL(d,e);
  }
};

#if defined(FORTARRAY) && defined(NUMRECIPES)
inline void Mat<double>::house(DiagMat<double>& d, double* e) {
  DTRED2(dat,&sz,&sz,d.dat,e);
}

inline void Mat<double>::triQL(DiagMat<double>& d, double* e) {
  DTQLI(d.dat,e,&sz,&sz,dat);
}
#endif

template <class T>
RangeRecMat<T> Mat<T>::operator()(RowRange rsl, ColRange csl) const {
  return RangeRecMat<T>(rsl,csl,m);
};


template <class T>
class TransMat : public Dim2<T,TransMat<T> > {
private:
  Mat<T>& mat;

  template <class X> Mat<T>& operator=(X& rhs);

public:
  explicit TransMat(Mat<T>& mat_) : mat(mat_) {}
  ~TransMat() {}
  int size() const { return mat.size(); }
  int rows() const { return mat.rows(); }
  int cols() const { return mat.cols(); }
  T& operator()(int i, int j) { return mat(j,i); }
  T operator()(int i, int j) const { return mat(j,i); }
};

#endif // VECMAT_H_
