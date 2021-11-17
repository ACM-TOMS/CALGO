#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "random.h"


int main(int argc, char * argv[]) {

  const int seed = 1, n = 1000000, ndiv = 1000;
  int dist[ndiv];
  Random rnd(seed);

  int i;
  for (i=0; i<ndiv; i++) dist[i] = 0;

  double xprev = rnd(1.0)-0.5, s1 = xprev, x0 = xprev;
  double s2 = x0* x0 , r = 0.0, x;
  
  for (i=0; i<n; i++) {
    x = rnd(1.0);
    dist[int(x*ndiv)]++;

    x -= 0.5;
    s1 += x; s2 += x*x; r += xprev * x; xprev = x;
  }
  r += x0*x;

  double dn = double(n),
    corr = ( dn*r - s1*s1 ) / ( dn*s2 - s1*s1 );

  double navr = s1*sqrt(12.0/dn);

  cout << "With 95 % probability, -2 < Normalized average < 2 " << endl;
  cout << "It's " << navr ;
  assert(-2.0 < navr && navr < 2.0);
  cout << " .... OK." << endl;

  double ncorr = ((dn-1.0)*corr+1.0)*sqrt((dn+1.0)/(dn*(dn-3.0)));

  cout << "With 95 % probability, -2 < Normalized Neighbor Corr. < 2 " << endl;
  cout << "It's " << ncorr; 
  assert(-2.0 < ncorr && ncorr < 2.0);
  cout << " ..... OK." << endl;

  int dev = 0;
  for (int i=0; i<ndiv; i++) 
    dev += ( dist[i] - n/ndiv ) * ( dist[i] - n/ndiv );


  cout << "Root mean square deviation from uniform distribution : ";
  cout << sqrt( double(dev)/double(ndiv) ) / (n/ndiv) << endl;

  return 0;
}
