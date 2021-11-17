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
test_rot (void) {
const MpIeee flteps=  1e-4;const MpIeee  dbleps=  1e-6;
  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.314f };
   int incX = 1;
   MpIeee Y[] =  { -0.406f };
   int incY = -1;
   MpIeee x_expected[] =  { 0.0f };
   MpIeee y_expected[] =  { 0.0f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 558)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 559)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784f;
   MpIeee s=  0.5f;
   MpIeee X[] =  { -0.314f };
   int incX = 1;
   MpIeee Y[] =  { -0.406f };
   int incY = -1;
   MpIeee x_expected[] =  { -0.474932f };
   MpIeee y_expected[] =  { -0.194606f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 560)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 561)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  -1.0f;
   MpIeee X[] =  { -0.314f };
   int incX = 1;
   MpIeee Y[] =  { -0.406f };
   int incY = -1;
   MpIeee x_expected[] =  { 0.406f };
   MpIeee y_expected[] =  { -0.314f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 562)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 563)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.314f };
   int incX = 1;
   MpIeee Y[] =  { -0.406f };
   int incY = -1;
   MpIeee x_expected[] =  { 0.314f };
   MpIeee y_expected[] =  { 0.406f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 564)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 565)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  0;
   MpIeee X[] =  { -0.493 };
   int incX = 1;
   MpIeee Y[] =  { -0.014 };
   int incY = -1;
   MpIeee x_expected[] =  { 0.0 };
   MpIeee y_expected[] =  { 0.0 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 566)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 567)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784;
   MpIeee s=  0.5;
   MpIeee X[] =  { -0.493 };
   int incX = 1;
   MpIeee Y[] =  { -0.014 };
   int incY = -1;
   MpIeee x_expected[] =  { -0.433950524066 };
   MpIeee y_expected[] =  { 0.234375644347 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 568)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 569)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  -1;
   MpIeee X[] =  { -0.493 };
   int incX = 1;
   MpIeee Y[] =  { -0.014 };
   int incY = -1;
   MpIeee x_expected[] =  { 0.014 };
   MpIeee y_expected[] =  { -0.493 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 570)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 571)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1;
   MpIeee s=  0;
   MpIeee X[] =  { -0.493 };
   int incX = 1;
   MpIeee Y[] =  { -0.014 };
   int incY = -1;
   MpIeee x_expected[] =  { 0.493 };
   MpIeee y_expected[] =  { 0.014 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 572)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 573)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.808f };
   int incX = -1;
   MpIeee Y[] =  { -0.511f };
   int incY = 1;
   MpIeee x_expected[] =  { 0.0f };
   MpIeee y_expected[] =  { 0.0f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 574)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 575)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784f;
   MpIeee s=  0.5f;
   MpIeee X[] =  { -0.808f };
   int incX = -1;
   MpIeee Y[] =  { -0.511f };
   int incY = 1;
   MpIeee x_expected[] =  { -0.955249f };
   MpIeee y_expected[] =  { -0.038539f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 576)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 577)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  -1.0f;
   MpIeee X[] =  { -0.808f };
   int incX = -1;
   MpIeee Y[] =  { -0.511f };
   int incY = 1;
   MpIeee x_expected[] =  { 0.511f };
   MpIeee y_expected[] =  { -0.808f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 578)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 579)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.808f };
   int incX = -1;
   MpIeee Y[] =  { -0.511f };
   int incY = 1;
   MpIeee x_expected[] =  { 0.808f };
   MpIeee y_expected[] =  { 0.511f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 580)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 581)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  0;
   MpIeee X[] =  { -0.176 };
   int incX = -1;
   MpIeee Y[] =  { -0.165 };
   int incY = 1;
   MpIeee x_expected[] =  { 0.0 };
   MpIeee y_expected[] =  { 0.0 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 582)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 583)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784;
   MpIeee s=  0.5;
   MpIeee X[] =  { -0.176 };
   int incX = -1;
   MpIeee Y[] =  { -0.165 };
   int incY = 1;
   MpIeee x_expected[] =  { -0.234920471066 };
   MpIeee y_expected[] =  { -0.0548941916244 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 584)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 585)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  -1;
   MpIeee X[] =  { -0.176 };
   int incX = -1;
   MpIeee Y[] =  { -0.165 };
   int incY = 1;
   MpIeee x_expected[] =  { 0.165 };
   MpIeee y_expected[] =  { -0.176 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 586)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 587)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1;
   MpIeee s=  0;
   MpIeee X[] =  { -0.176 };
   int incX = -1;
   MpIeee Y[] =  { -0.165 };
   int incY = 1;
   MpIeee x_expected[] =  { 0.176 };
   MpIeee y_expected[] =  { 0.165 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 588)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 589)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.201f };
   int incX = -1;
   MpIeee Y[] =  { 0.087f };
   int incY = -1;
   MpIeee x_expected[] =  { 0.0f };
   MpIeee y_expected[] =  { 0.0f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 590)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 591)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784f;
   MpIeee s=  0.5f;
   MpIeee X[] =  { -0.201f };
   int incX = -1;
   MpIeee Y[] =  { 0.087f };
   int incY = -1;
   MpIeee x_expected[] =  { -0.130571f };
   MpIeee y_expected[] =  { 0.175844f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 592)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 593)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.0f;
   MpIeee s=  -1.0f;
   MpIeee X[] =  { -0.201f };
   int incX = -1;
   MpIeee Y[] =  { 0.087f };
   int incY = -1;
   MpIeee x_expected[] =  { -0.087f };
   MpIeee y_expected[] =  { -0.201f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 594)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 595)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1.0f;
   MpIeee s=  0.0f;
   MpIeee X[] =  { -0.201f };
   int incX = -1;
   MpIeee Y[] =  { 0.087f };
   int incY = -1;
   MpIeee x_expected[] =  { 0.201f };
   MpIeee y_expected[] =  { -0.087f };
   cblas_srot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], flteps, "srot(case 596)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], flteps, "srot(case 597)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  0;
   MpIeee X[] =  { -0.464 };
   int incX = -1;
   MpIeee Y[] =  { 0.7 };
   int incY = -1;
   MpIeee x_expected[] =  { 0.0 };
   MpIeee y_expected[] =  { 0.0 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 598)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 599)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0.866025403784;
   MpIeee s=  0.5;
   MpIeee X[] =  { -0.464 };
   int incX = -1;
   MpIeee Y[] =  { 0.7 };
   int incY = -1;
   MpIeee x_expected[] =  { -0.051835787356 };
   MpIeee y_expected[] =  { 0.838217782649 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 600)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 601)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  0;
   MpIeee s=  -1;
   MpIeee X[] =  { -0.464 };
   int incX = -1;
   MpIeee Y[] =  { 0.7 };
   int incY = -1;
   MpIeee x_expected[] =  { -0.7 };
   MpIeee y_expected[] =  { -0.464 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 602)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 603)");
     }
   };
  };


  {
   int N = 1;
   MpIeee c=  -1;
   MpIeee s=  0;
   MpIeee X[] =  { -0.464 };
   int incX = -1;
   MpIeee Y[] =  { 0.7 };
   int incY = -1;
   MpIeee x_expected[] =  { 0.464 };
   MpIeee y_expected[] =  { -0.7 };
   cblas_drot(N, X, incX, Y, incY, c, s);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(X[i], x_expected[i], dbleps, "drot(case 604)");
     }
   };
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(Y[i], y_expected[i], dbleps, "drot(case 605)");
     }
   };
  };


}
