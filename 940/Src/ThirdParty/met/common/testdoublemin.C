#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

#include <math.h>
#include <float.h>

int main(int argc, char * argv[]) {

  for (int i=0; i<1000; i++) {
    double di = double(i),
      e  = exp( -di ),
      l = -log(e*(1.0 - 1.0e-4));
    if ( di > l ) {
      cout << "Limit of Accuracy !" << endl;
      cout << i << " > " << " -log( -exp( -" << i << " )*(1.0 - 1.0e-4) )" << endl;
      break;
    }
  }

  return 0;
}
