#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "vecmat.h"

int main(int argc, char * argv[]) {
  double err, de;
  int r,c,k;

  RecMat<double> a(8,7), b(8,7), s(8,7);
  DiagMat<double> d(7), f(7), g(8);
  Vec<double> x(7), y(8), w(8);

  RowRange i8(0,8), i7(0,7);
  ColRange j7(0,7), j8(0,8);

  a(i8,j7) = - dble(i8*4) * 2.0 + 3.0 * dble(j7*2);
  b(i8,j7) =  3.0 * dble(i8) + exp( dble(j7*3) * 0.1 );

  DiagRange k7(0,7), k8(0,8);
  d(k7) = 2.0 * dexp(k7) - dble(k7) * 7.0;
  f(k7) = - dble(k7*3) + 1.0;
  g(k8) = exp( dble(k8*2) * 0.1 ) + 5.0;

  // a, b are the same as those in xpr2Test.C
  // d,f,g are the same as those in xprDTest.C

  //s = - ( 1.5 * a*d - 3.0 * g*b ) * f * 2.0 + a;
  // This cannot be compiled by egcs-1.1.2
  //xpr2.h:127: template instantiation depth exceeds maximum of 17
  //xpr2.h:127:  (use -ftemplate-depth-NN to increase the maximum)

  //s = - ( 1.5 * a*d - 3.0 * g*b ) * f + a; // too long expression to be compiled ??
  s = 1.5 * a*d - 3.0 * g*b;

  err = 0.0;
  for (r=0; r<8; r++) {
    for (c=0; c<7; c++) {
      //de = - ( 1.5 * a(r,c)*d(c) - 3.0 * g(r)*b(r,c) ) * f(c) +
      //	a(r,c) - s(r,c);
      de = 1.5 * a(r,c)*d(c) - 3.0 * g(r)*b(r,c) - s(r,c);
      err += de*de;
    }
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  Range l7(0,7), l8(0,8);
  x(l7) = 7.0 - dble(l7);
  y(l8) = exp( dble(l8) * 0.3 );

  // x,y are the same as those in xpr2Test.C

  //w = g*( 2.0*( a*x ) + g*y ) - 3.0 * (g*y);
  w = 2.0*( a*x ) + g*y;
  
  err = 0.0;
  for (r=0; r<8; r++) {
    de = 0.0;
    //de = g(r) * g(r) * y(r) - 3.0 *g(r)*y(r);
    de = g(r) * y(r);
    for (k=0; k<7; k++) 
      //de += g(r)* ( 2.0*a(r,k)*x(k) );
      de += 2.0*a(r,k)*x(k);
    de -= w(r);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  return 0;
}
