// Examples for C Interface with pre-compiled library

#define FFT_PRECOMPILE //in case not defined in genial_config.h
#include "signal/fft.h" 
int main()
{
  DenseVector<float>::self X(4096, 1);
  DenseVector<complex<float> >::self Y(X.size());
  
  double t0=get_clock();
  for (int i=0; i<100000; ++i)
    shfft(X.size(),Y.size(),&X[0],(float *)&Y[0]);
    //sfft (X.size(),&X[0],(float *)&Y[0]);
  cout << get_clock()-t0 << "s" << endl;
}

//#define FFT_PRECOMPILE //in case not defined in genial_config.h
//#include "signal/fft.h" 
//int main()
//{
//  DenseVector<complex<float> >::self X(4096, 1);
//  DenseVector<float>::self Y(X.size());
//  
//  double t0=get_clock();
//  for (int i=0; i<100000; ++i)
//    sifft(X.size(),(float *)&X[0],&Y[0]);
//  cout << get_clock()-t0 << "s" << endl;
//}

//#define FFT_PRECOMPILE //in case not defined in genial_config.h
//#include "signal/fft.h" 
//int main()
//{
//  DenseVector<complex<float> >::self X(4096, 1);
//  DenseVector<complex<float> >::self Y(X.size());
//  
//  double t0=get_clock();
//  for (int i=0; i<100000; ++i)
//    cfft(X.size(),(float *)&X[0],(float *)&Y[0]);
//  cout << get_clock()-t0 << "s" << endl;
//}

//#define FFT_PRECOMPILE //in case not defined in genial_config.h
//#include "signal/fft.h" 
//int main()
//{
//  DenseVector<complex<float> >::self X(4096, 1);
//  DenseVector<complex<float> >::self Y(X.size());
//  
//  double t0=get_clock();
//  for (int i=0; i<100000; ++i)
//    cifft(X.size(),(float *)&X[0],(float *)&Y[0]);
//  cout << get_clock()-t0 << "s" << endl;
//}

//#define FFT_PRECOMPILE //in case not defined in genial_config.h
//#include "signal/dct.h" 
//int main()
//{
//  DenseVector<float>::self X(4096, 1);
//  DenseVector<float>::self Y(X.size());
//  
//  double t0=get_clock();
//  for (int i=0; i<100000; ++i)
//    sdct(X.size(),&X[0],(float *)&Y[0]);
//  cout << get_clock()-t0 << "s" << endl;
//}

//#define BLAS_PRECOMPILE //in case not defined in genial_config.h
//#include "blas/gemm.h" 
//int main()
//{
//  int m=200, n=200, k=200;
//  typedef float value_type;
//  
//  DenseMatrix<value_type>::self A(m,k,1), B(k,n,1), C(m,n,0);
//
//  double t0=get_clock();
//  for (int i=0; i<1000; ++i)
//    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1,&A(0,0),k,&B(0,0),n,0,&C(0,0),n);
//  cout << get_clock()-t0 << "s" << endl;
//}


//#define BLAS_PRECOMPILE //in case not defined in genial_config.h
//#include "blas/gemm.h" 
//int main()
//{
//  int m=1000, n=1000, k=1000;
//  typedef complex<float> value_type;
//  
//  DenseMatrix<value_type>::self A(m,k,1), B(k,n,1), C(m,n,0);
//  value_type alpha(1,0), beta(0,0);
//
//  double t0=get_clock();
//  for (int i=0; i<5; ++i)
//    cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,(void*)&alpha,(void*)&A(0,0),k,(void*)&B(0,0),n,(void*)&beta,(void*)&C(0,0),n);
//  cout << get_clock()-t0 << "s" << endl;
//}



//#define MOTION_ESTIMATION_PRECOMPILE //in case not defined in genial_config.h
//#include "image/motion.h" 
//
//int main()
//{
//  int m=720, n=512;
//  
//  ucharImage X(m,n, 0);
//  ucharImage Y(m,n, 0);
//  DenseMatrix<ucharImage::index_type>::self M(X.size()/8);
//
//  double t0=get_clock();
//  for (int i=0; i<100; ++i)
//    motion8(m,n,&X(0,0),&Y(0,0),16,16, (int *)&M(0,0));
//  cout << get_clock()-t0 << "s" << endl;
//}
