#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "vecmat.h"
#include "solver.h"

int main(int argc, char * argv[]) {
  Mat<double> a(5), u(5), ud(5), s(5);
  TransMat<double> tu(u);
  DiagMat<double> d(5);
  double subdiag[5];

  Vec<double> b(5), x(5), ans(5), w(5);

  a(0,0) = 1.0; a(0,1) = 2.0; a(0,2) = 3.0; a(0,3) =-2.0; a(0,4) =-1.0;
  a(1,0) = 2.0; a(1,1) = 8.0; a(1,2) =22.0; a(1,3) =10.0; a(1,4) =-0.5;
  a(2,0) = 3.0; a(2,1) =11.0; a(2,2) =35.0; a(2,3) =-8.0; a(2,4) = 7.0;
  a(3,0) =-1.3; a(3,1) = 0.3; a(3,2) =10.0; a(3,3) =15.0; a(3,4) = 0.0;
  a(4,0) =10.0; a(4,1) =-3.0; a(4,2) = 1.0; a(4,3) =25.0; a(4,4) =50.0;  
  
  Range k5(0,5);
  x(k5) = 5.0 - dble(k5)*0.7;

  b = a*x;

  u = a;

  //cout << "a = " << endl;
  //cout << a;

  MatLU<double> asolver(u);
  asolver.update();
  asolver.solve(b,ans);

  //cout << "LU decomposed a = " << endl;
  //cout << u;

  //cout << "x = " << endl;
  //cout << x;
  //cout << "ans = " << endl;
  //cout << ans;

  double err;
  ans -= x;
  err = ans.norm();
  cout << "error of MatLU = " << err << endl;
  assert(err < 1e-10);


  // Symmetrize
  int r,c;
  for (r=0; r<5; r++) {
    for (c=r+1; c<5; c++) a(r,c) = a(c,r) = 0.5*( a(r,c) + a(c,r) );
  }
  b = a*x;

  a.diagonalize(u,d,subdiag);

  ud = u * d;
  s = ud * tu - a;
  err = 0.0;
  for (r=0; r<5; r++) {
    for (c=0; c<5; c++) err += s(r,c)*s(r,c);
  }
  cout << "error of diagonalization = " << err << endl;
  assert(err < 1e-10);

  MatDiagonalized<double> da(a);
  da.update();
  da.solve(b,ans);

  //for (int n=0; n<5; n++) d(n) = 1.0/d(n);
  //ans = tu*b;
  //w = d*ans;
  //ans = u*w;

  ans -= x;
  cout << "error of Similartiy Solving = " << ans.norm() << endl;
  assert(err < 1e-10);

  return 0;
}



