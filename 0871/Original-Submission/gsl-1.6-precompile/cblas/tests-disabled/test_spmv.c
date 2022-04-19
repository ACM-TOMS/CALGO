#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>

#include "tests.h"

void
test_spmv (void) {
const MpIeee flteps=  1e-4;const MpIeee  dbleps=  1e-6;
  {
   int order = 101;
   int uplo = 121;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1134)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1135)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1136)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1137)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1138)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1139)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1140)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   MpIeee alpha=  0.1f;
   MpIeee beta=  -0.3f;
   int N = 2;
   MpIeee A[] =  { -0.174f, 0.878f, 0.478f };
   MpIeee X[] =  { 0.503f, 0.313f };
   int incX = -1;
   MpIeee Y[] =  { -0.565f, -0.109f };
   int incY = -1;
   MpIeee y_expected[] =  { 0.221025f, 0.0714172f };
   cblas_sspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "sspmv(case 1141)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1142)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1143)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1144)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1145)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1146)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1147)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1148)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   MpIeee alpha=  -1;
   MpIeee beta=  0.1;
   int N = 2;
   MpIeee A[] =  { -0.181, -0.071, -0.038 };
   MpIeee X[] =  { -0.015, 0.132 };
   int incX = -1;
   MpIeee Y[] =  { -0.449, -0.219 };
   int incY = -1;
   MpIeee y_expected[] =  { -0.036098, 9.27e-04 };
   cblas_dspmv(order, uplo, N, alpha, A, X, incX, beta, Y, incY);
   {
     int i;
     for (i = 0; i < 2; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "dspmv(case 1149)");
     }
   };
  };


}
