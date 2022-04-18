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
test_syr2 (void) {
const MpIeee flteps=  1e-4;const MpIeee  dbleps=  1e-6;
  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0.0f;
   MpIeee A[] =  { 0.862f };
   MpIeee X[] =  { 0.823f };
   int incX = -1;
   MpIeee Y[] =  { 0.699f };
   int incY = -1;
   MpIeee A_expected[] =  { 0.862f };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1434)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0.0f;
   MpIeee A[] =  { 0.862f };
   MpIeee X[] =  { 0.823f };
   int incX = -1;
   MpIeee Y[] =  { 0.699f };
   int incY = -1;
   MpIeee A_expected[] =  { 0.862f };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1435)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0.0f;
   MpIeee A[] =  { 0.862f };
   MpIeee X[] =  { 0.823f };
   int incX = -1;
   MpIeee Y[] =  { 0.699f };
   int incY = -1;
   MpIeee A_expected[] =  { 0.862f };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1436)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0.0f;
   MpIeee A[] =  { 0.862f };
   MpIeee X[] =  { 0.823f };
   int incX = -1;
   MpIeee Y[] =  { 0.699f };
   int incY = -1;
   MpIeee A_expected[] =  { 0.862f };
   cblas_ssyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], flteps, "ssyr2(case 1437)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0;
   MpIeee A[] =  { -0.824 };
   MpIeee X[] =  { 0.684 };
   int incX = -1;
   MpIeee Y[] =  { 0.965 };
   int incY = -1;
   MpIeee A_expected[] =  { -0.824 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1438)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0;
   MpIeee A[] =  { -0.824 };
   MpIeee X[] =  { 0.684 };
   int incX = -1;
   MpIeee Y[] =  { 0.965 };
   int incY = -1;
   MpIeee A_expected[] =  { -0.824 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1439)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0;
   MpIeee A[] =  { -0.824 };
   MpIeee X[] =  { 0.684 };
   int incX = -1;
   MpIeee Y[] =  { 0.965 };
   int incY = -1;
   MpIeee A_expected[] =  { -0.824 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1440)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int N = 1;
   int lda = 1;
   MpIeee alpha=  0;
   MpIeee A[] =  { -0.824 };
   MpIeee X[] =  { 0.684 };
   int incX = -1;
   MpIeee Y[] =  { 0.965 };
   int incY = -1;
   MpIeee A_expected[] =  { -0.824 };
   cblas_dsyr2(order, uplo, N, alpha, X, incX, Y, incY, A, lda);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(A[i], A_expected[i], dbleps, "dsyr2(case 1441)");
     }
   };
  };


}
