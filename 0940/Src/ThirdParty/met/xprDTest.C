#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "vecmat.h"

int main(int argc, char * argv[]) {
  double err, de;
  int n,k;

  DiagMat<double> a(7), b(7), ab(7), d(7), ad(7);
  Vec<double> x(7), y(7), w(7);

  DiagRange i7(0,7);

  a(i7) = 2.0 * dexp(i7) - dble(i7) * 7.0; 
  b(i7) = - dble(i7*3) + 1.0;
  d(i7) = exp( dble(i7*2) * 0.1 ) + 5.0;

  err = 0.0;
  for (n=0; n<7; n++) {
    de = 2.0*exp( double(n) ) - double(n)*7.0 - a(n);
    err += de*de;
  }
  assert(err < 1e-10);

  err = 0.0;
  for (n=0; n<7; n++) {
    de =  -double(n*3) + 1.0 - b(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (n=0; n<7; n++) {
    de = exp( double(n*2) * 0.1 ) + 5.0 - d(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  ab = a + b;
  ad = a * d;

  err = 0.0;
  for (n=0; n<7; n++) {
    de = a(n) + b(n) - ab(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
    for (n=0; n<7; n++) {
      de = a(n) * d(n) - ad(n);
      err += de*de;
    }
  cout << "error = " << err << endl;
  assert(err < 1e-10);
  
  Range i(0,7);
  x(i) = 7.0 - dble(i);
  y(i) = exp( dble(i) * 0.3 );

  err = 0.0;
  for (n=0; n<7; n++) {
    de = 7.0 - double(n) - x(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  err = 0.0;
  for (n=0; n<7; n++) {
    de = exp( double(n) * 0.3 ) - y(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  w = ( a*x + ab*x )*3.0 - b * y;
  
  err = 0.0;
  for (n=0; n<7; n++) {
    de = ( a(n)*x(n) + ab(n)*x(n) ) * 3.0 - b(n)*y(n) - w(n);
    err += de*de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  return 0;
}
