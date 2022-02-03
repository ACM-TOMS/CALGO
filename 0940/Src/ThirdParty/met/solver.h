//-*- C++ -*-
//
//  Solver Class Declaration
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

#ifndef SOLVER_H_
#define SOLVER_H_

#include "mathfunc.h"
#include "vecmat.h"

template <class T>
class MatLU {
private:
  Mat<T>& mat;
  int* const ip; // for pivot selection
  //T* const y; // work for transSolve()

public:
  // Just allocating memory, DON'T FORGET TO CALL update() before solving.
  MatLU(Mat<T>& mat_) : 
    mat(mat_),
    ip(new int[mat.rows()]) //,y(new T[mat.rows()])
  {}
    //{ reset(); }
  ~MatLU() { 
    //delete [] y;
    delete [] ip; 
  }

  //void reset() { mat.luDecompose(ip); }
  void update() { mat.luDecompose(ip); }

  template <class V>
  void solve(const Dim1<T,V>& b, Vec<T>& x) const { mat.luSubst(ip,b,x); }
  //void transSolve(const Vec<T>& b, Vec<T>& x) const {
  //  m.transLUSubst(ip,y,b,x);
  //}
  //void conjSolve(const Vec<T>& b, Vec<T>& x) const {
  //  m.conjLUSubst(ip,y,b,x);
  //}
};


template <class T>
class MatDiagonalized {
private:
  Mat<T>& mat;
  Mat<T> u;
  TransMat<T> tu;
  DiagMat<T> d, invd;
  Vec<T> tempx;
  T *subdiag;

public:
  // Just allocating memory, DON'T FORGET TO CALL update() before solving.
  MatDiagonalized(Mat<T>& mat_) : 
    mat(mat_), u(mat_.rows()), tu(u), d(mat_.rows()), invd(mat_.rows()),
    tempx(mat_.rows()),
    subdiag(new T[mat_.rows()]) {}

  ~MatDiagonalized() { delete [] subdiag; }
  
  void update() { 
    u = mat;
    u.house(d,subdiag);
    u.triQL(d,subdiag);
  }

  template <class V>
  void solve(const Dim1<T,V>& b, Vec<T>& x) { 
    for (int i=0; i<d.size(); i++) invd(i) = 1.0/d(i);

    // x = u * d^(-1) * transpose(u) * b;
    x = tu * b;
    tempx = invd * x;
    x = u * tempx;
  }

  template <class V>
  void solve(const Dim1<T,V>& b, Vec<T>& x, T smallest_eigenvalue) { 
    for (int i=0; i<d.size(); i++) 
      if ( Abs( d(i) ) < smallest_eigenvalue ) {
	invd(i) = 0.0;
      } else {
	invd(i) = 1.0/d(i);
      }

    // x = u * d^(-1) * transpose(u) * b;
    x = tu * b;
    tempx = invd * x;
    x = u * tempx;
  }
};


#endif // SOLVER_H_


