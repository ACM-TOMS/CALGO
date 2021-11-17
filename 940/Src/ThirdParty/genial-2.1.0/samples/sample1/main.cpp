#include "array/vector.h"

int main()
{
  DenseVector<int>::self X(4, 1);

  cout << X << endl;
  cout << 2*X << endl;
}


//#include "array/matrix.h"
//
//int main()
//{
//  DenseMatrix<int>::self X(4,4, 1);
//
//  cout << X << endl;
//  cout << 2*X << endl;
//}


//#define FFT_LEVEL 4
//#include "signal/fft.h"
//
//int main()
//{
//  int n=8;
//  DenseVector<complex<float> >::self X(n,1), Y(n);
//  fft(X,Y);
//  cout << Y << endl;
//}


//#define FFT_LEVEL 4
//#include "signal/fft.h"
//
//int main()
//{
//  int n=8;
//  DenseVector<float>::self X(n,1);
//  DenseVector<complex<float> >::self Y = fft(X);
//  cout << Y << endl;
//}


//#define FFT_LEVEL 4
//#include "signal/dct.h"
//
//int main()
//{
//  int n=8;
//  DenseVector<float>::self X(n,1);
//  cout << idct(dct(X)) << endl;
//}


//#include "array/vector.h"
//#include "signal/conv.h"
//
//int main()
//{
//  DenseVector<float>::self X(4, 1);
//  DenseVector<float>::self Y(3, "0.25 0.5 0.25");
//  DenseVector<float>::self Z;
//  Z=conv(X,Y);
//  cout << Z << endl;
//}


//#include "array/vector.h"
//#include "signal/conv.h"
//
//int main()
//{
//  DenseVector<float>::self X(4, 1);
//  shiftDenseVector<float>::self Y(3, "0.25 0.5 0.25",-1);
//  DenseVector<float>::self Z;
//  Z=sub(conv(X,Y),0,X.size());
//  cout << Z << endl;
//  cout << conv(X,Y) << endl;
//}


//#include "blas/gemm.h"
//
//int main()
//{
//  int m=4,n=4,k=4;
//  DenseMatrix<float>::self A(m,k,1), B(k,n,1), C(m,n,0);
//  //C=mul(A,B);
//  gemm(1,A,B,0,C);
//  cout << C << endl;
//}

//#include "image/image.h"
//#include "image/ppm.h"
//#include "image/pgm.h"
//
//int main()
//{
//  try
//  {
//    PPMFile fin("x.ppm");
//    RGBImage X;
//    fin >> X;
//
//    PGMFile fout("y.pgm");
//    fout << first(value_cast<YUV>(X));
//  } catch (const error &e) { cout << e.what() << endl; return 1; }
//  return 0;
//}


//// Needs the Windows GDI+ Library, to link with gdiplus.lib
//#include "image/image.h"
//#include "image/jpg.h"
//#include "image/ioimage.h"
//
//int main()
//{
//  try
//  {
//    JPGFile fin("x.jpg"), fout("y.jpg",50);
//    RGBImage X;
//    fin >> X;
//    fout << X;
//  } catch (const error &e) { cout << e.what() << endl; return 1; }
//  return 0;
//}
