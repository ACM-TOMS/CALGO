#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "vecmat.h"

#if defined (UNDERF77) || defined(__DECCXX) || defined(__sgi)
#define CALCABC calcabc_
#else
#define CALCABC calcabc
#endif

extern "C" {
  void CALCABC(const int*, double*, double*, double*, double*, double*);
}

int main(int argc, char * argv[]) {
  if (argc <= 0) {
    cerr << "Usage : " << argv[0] << " vectorSize [vectorSize2 ..]" << endl;
    exit(1);
  }

  cout << "#size     c++       c        f77      error      f77/c++   f77/c" << endl;

  for (int t=1; t<argc; t++) {
    
    clock_t ctime, cpptime, forttime, alltime = clock();
    double diff, tdiff = 0.0;
    
    int size = atoi(argv[t]);
    
    double err, de;

    // Memory Allocation
    RecMat<double> a(size+1,size), b(size+1,size);
    Vec<double> x(size), y(size), w(size+1);
    
    RowRange r1(0,size+1);
    ColRange c0(0,size);
    Range i0(0,size);

    // C++ calculation
    cpptime = clock();

    a(r1,c0) = - dble(r1*4*7)/double(size) * 2.0 
      + 3.0 * dble(c0*2*7)/double(size);
    b(r1,c0) =  3.0 * dble(r1*7)/double(size) 
      + exp( dble(c0*3*7)/double(size) * 0.1 );

    x(i0) = 7.0 - dble(i0*7)/double(size);
    y(i0) = exp( dble(i0*7)/double(size) * 0.3 );

    w = ( a*x - b*y )*3.0;

    cpptime = clock() - cpptime;

    // Memory Allocation
    double *fa1 = new double[(size+1)*size], *fb1 = new double[(size+1)*size],
      **fa = new double*[size], **fb = new double*[size],
      *fx = new double[size], *fy = new double[size], *fw = new double[size+1];

    int r,c,i,k;
    for (c=0; c<size; c++) {
      fa[c] = fa1 + c * (size+1);
      fb[c] = fb1 + c * (size+1);
    }

    // C calculation
    ctime = clock();

    for (c=0; c<size; c++) {
      for (r=0; r<size+1; r++)
	fa[c][r] = -double(r*4*7)/double(size)*2.0 
	  + 3.0*double(c*2*7)/double(size);
    }
    for (c=0; c<size; c++) {
      for (r=0; r<size+1; r++)
	fb[c][r] = 3.0 * double(r*7)/double(size) 
	  + exp( double(c*3*7)/double(size) * 0.1 );
    }

    for (i=0; i<size; i++) fx[i] = 7.0 - double(i*7)/double(size);
    for (i=0; i<size; i++) fy[i] = exp( double(i*7)/double(size) * 0.3 );

    for (i=0; i<size+1; i++) fw[i] = 0.0;

    for (k=0; k<size; k++) {
      for (i=0; i<size+1; i++)
        fw[i] += ( fa[k][i]*fx[k]  - fb[k][i]*fy[k] )*3.0;
    }
    ctime = clock() - ctime;
  
    // Accuracy Check
    err = 0.0;
//     for (c=0; c<size; c++) {
//       for (r=0; r<size+1; r++) {
// 	de = a(r,c) - fa[c][r];
// 	err += de*de;
//       }
//     }
//     assert(err/size < 1e-10);
//     for (c=0; c<size; c++) {
//       for (r=0; r<size+1; r++) {
// 	de = b(r,c) - fb[c][r];
// 	err += de*de;
//       }
//     }
//     assert(err/size < 1e-10);
//     for (i=0; i<size; i++) {
//       de = x(i) - fx[i];
//       err += de*de;
//     }
//     assert(err/size < 1e-10);
//     for (i=0; i<size; i++) {
//       de = y(i) - fy[i];
//       err += de*de;
//     }
//     assert(err/size < 1e-10);
    for (i=0; i<size+1; i++) {
      de = w(i) - fw[i];
      err += de*de;
    }
    assert(err/size < 1e-7);
    tdiff = err;

    // Fortran calc.
    forttime = clock();
    CALCABC(&size,fa1,fb1,fx,fy,fw);
    forttime = clock() - forttime;

    // Accuracy Check
    err = 0.0;
//     for (c=0; c<size; c++) {
//       for (r=0; r<size+1; r++) {
// 	de = a(r,c) - fa[c][r];
// 	err += de*de;
//       }
//     }
//     assert(err/size < 1e-10);
//     for (c=0; c<size; c++) {
//       for (r=0; r<size+1; r++) {
// 	de = b(r,c) - fb[c][r];
// 	err += de*de;
//       }
//     }
//     assert(err/size < 1e-10);
//     for (i=0; i<size; i++) {
//       de = x(i) - fx[i];
//       err += de*de;
//     }
//     assert(err/size < 1e-10);
//     for (i=0; i<size; i++) {
//       de = y(i) - fy[i];
//       err += de*de;
//     }
//     assert(err/size < 1e-10);
    for (i=0; i<size+1; i++) {
      de = w(i) - fw[i];
      err += de*de;
    }
    assert(err/fw[size] < 1e-10);
    tdiff += err;

    tdiff /= w.abs();

    // Memory Release
    delete [] fw;
    delete [] fy;
    delete [] fx;
    delete [] fb;
    delete [] fa;
    delete [] fb1;
    delete [] fa1;

    alltime = clock() -alltime;

    cout << setw(5) << atoi(argv[t]);
    cout << setw(9) << setprecision(4) << double(cpptime)/double(CLOCKS_PER_SEC);
    cout << setw(9) << setprecision(4) << double(ctime)/double(CLOCKS_PER_SEC);
    cout << setw(9) << setprecision(4) << double(forttime)/double(CLOCKS_PER_SEC);
    cout.setf(ios::scientific,ios::floatfield);
    cout << setw(16) << setprecision(8) << tdiff;
    cout.setf(ios::fixed,ios::floatfield);

    if ( cpptime ) {
      cout << setw(7) << setprecision(3) << double(forttime)/double(cpptime);
    } else {
      cout << setw(7) << "  N/A";
    }
    
    if ( ctime ) {
    cout << setw(7) << setprecision(3) << double(forttime)/double(ctime);
    } else {
      cout << setw(7) << "  N/A";
    }
    cout << endl;
  }


  return 0;
}
