#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "xpr1.h"
#include "xpr2.h"
#include "vecmat.h"

#if defined (UNDERF77) || defined(__DECCXX) || defined(__sgi)
#define CALCABC calcabc_
#else
#define CALCABC calcabc
#endif

extern "C" {
  void CALCABC(const int*,const double*, const double*,const double*, const double*, double*);
}

int main(int argc, char * argv[]) {
  if (argc <= 0) {
    cerr << "Usage : " << argv[0] << " vectorSize [vectorSize2 ..]" << endl;
    exit(1);
  }
  
  cout << "#size     c++       c        f77      error      f77/c++   f77/c" << endl;
  
  for (int k=1; k<argc; k++) {
    
    clock_t ctime, cpptime, forttime, alltime = clock();
    double diff, tdiff = 0.0;
    
    int size = atoi(argv[k]);
    
    int k2,j;
    // Memory Allocation 
    Vec<double> a(size), b(size), c(size);
    Range i(0,size,1);
    
    // C++ calculation
    cpptime = clock();
    double sa = 3.0, sb = 2.0;
    for (k2 = 1; k2 < 1000; k2++) {
      a(i) = dble(-i + size);
      b(i) = dble(i*3 + 2) + exp(dble(i)/double(size));
      //c = -a*sa + sb*b;
      c = a + b;
    }
    cpptime = clock() - cpptime;
    
    // Memory Allocation
    double *fa = new double[size], *fb = new double[size], 
    *fc = new double[size];

    // C calc.
    ctime = clock();
    for (k2 =1; k2 < 1000; k2++) {
      for (j=0;j<size; j++) 
        fa[j] = double(-j+size);
      for (j=0;j<size; j++) 
        fb[j] = double(j*3 + 2) + exp(double(j)/double(size));
      //for (j=0;j<size; j++)
      //  fc[j] = - fa[j]*sa + sb*fb[j];
      for (j=0;j<size; j++)
        fc[j] = fa[j] + fb[j];
    }
    ctime = clock() - ctime;

    // Fortran calc.
    forttime = clock();
    for (k2 = 1; k2 < 1000; k2++) {
      CALCABC(&size,fa,fb,&sa,&sb,fc);
    }
    forttime = clock() - forttime;
    
    // Accurracy Check
    for (j=0; j<size; j++) {
      diff = c(j) - fc[j];
      tdiff += diff * diff;
    }
    // Memory Release
    delete [] fc;
    delete [] fb;
    delete [] fa;

    alltime = clock() -alltime;

    cout << setw(5) << atoi(argv[k]);
    cout << setw(9) << setprecision(4) << double(cpptime)/double(CLOCKS_PER_SEC);
    cout << setw(9) << setprecision(4) << double(ctime)/double(CLOCKS_PER_SEC);
    cout << setw(9) << setprecision(4) << double(forttime)/double(CLOCKS_PER_SEC);
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(16) << setprecision(8) << tdiff;
    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(7) << setprecision(3) << double(forttime)/double(cpptime);
    cout << setw(7) << setprecision(3) << double(forttime)/double(ctime) << endl;
    assert(tdiff < 1e-10);

  }

  return 0;
}




