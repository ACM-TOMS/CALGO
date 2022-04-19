#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "vecmat.h"

int main(int argc, char * argv[]) {
  double err, de;
  int r,c,k;

  RecMat<double> a(8,7), b(8,7), ab(8,7), d(7,8);
  Mat<double> ad(8);
  Vec<double> x(7), y(7), w(8);

  RowRange i8(0,8), i7(0,7);
  ColRange j7(0,7), j8(0,8);

  a(i8,j7) = - dble(i8*4) * 2.0 + 3.0 * dble(j7*2);
  b(i8,j7) =  3.0 * dble(i8) + exp( dble(j7*3) * 0.1 );
  d(i7,j8) = 2.0 * dexp(i7) - dble(j8) * 7.0;

  err = 0.0;
  for (r=0; r<8; r++) {
    for (c=0; c<7; c++) {
      de = -double(r*4)*2.0 + 3.0*double(c*2) - a(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (r=0; r<8; r++) {
    for (c=0; c<7; c++) {
      de = 3.0*double(r) + exp( double(c*3) * 0.1) - b(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (r=0; r<7; r++) {
    for (c=0; c<8; c++) {
      de = 2.0 * exp(double(r)) - double(c)*7.0 - d(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);


  ab = a + b;
  ad = a * d;

  err = 0.0;
  for (r=0; r<8; r++) {
    for (c=0; c<7; c++) {
      de = a(r,c) + b(r,c) - ab(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (r=0; r<8; r++) {
    for (c=0; c<8; c++) {
      de = 0.0;
      for (k=0; k<7; k++) de += a(r,k) * d(k,c);
      de -= ad(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);
  
  Range i(0,7);
  x(i) = 7.0 - dble(i);
  y(i) = exp( dble(i) * 0.3 );

  err = 0.0;
  for (r=0; r<7; r++) {
    de = 7.0 - double(r) - x(r);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (r=0; r<7; r++) {
    de = exp( double(r) * 0.3 ) - y(r);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  w = ( a*x + ab*x )*3.0 - b * y;
  
  err = 0.0;
  for (r=0; r<8; r++) {
    de = 0.0;
    for (k=0; k<7; k++) 
      de += ( a(r,k)*x(k) + ab(r,k)*x(k) ) * 3.0 - b(r,k)*y(k);
    de -= w(r);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  return 0;
}
