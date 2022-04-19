#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "vecmat.h"

int main(int argc, char * argv[]) {
  Vec<double> a(7), b(7), c(7), x(7);
  Range i(0,7,1);
  double err, de;
  int j;

  a = 3.0;
  a = ( b + c ) + x;
  a = x + ( b + c);

  a(i) = 1.5*dble(-i + 7);
  cout << "a = " << endl;
  cout << a;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = 1.5*double(-j+7) - a(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  b(i) = dble(i*4)*2.0 + dexp(i);
  //b = -a*3.0;
  cout << "b = " << endl;
  cout << b;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = j*4.0*2.0 + exp(double(j)) - b(j);
    //de = b(j) + a(j)*3.0;
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  c = dble(i*3);
  cout << "c = " << endl;
  cout << c;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = j*3.0 - c(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  x = -a*3.0 + 2.0*b - c; // -2*i
  cout << " x = -a + b - c =" << endl;
  cout << x;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = - a(j)*3.0 + 2.0*b(j) - c(j) - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  x = exp(a);
  cout << "x = exp(a);" << endl;
  cout << x;
  err = 0.0;
  for (j=0; j<7; j++) {
    de = exp(a(j)) - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  Range sl02(0,2);
  x(sl02) = dexp( sl02 );
  cout << "exp( Range(1.0,2.0) ) = " << endl;
  cout << x(sl02);
  err = 0.0;
  for (j=0; j<2; j++) {
    de = exp(double(sl02(j))) - x(j);
    err += de;
  }
  cout << "error = " << err << endl;
  assert(err < 1e-10);

  return 0;
}
