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
test_dot (void) {
const MpIeee flteps=  1e-4;const MpIeee  dbleps=  1e-6;
  {
   int N = 1;
   MpIeee alpha=  0.0f;
   MpIeee X[] =  { 0.733f };
   MpIeee Y[] =  { 0.825f };
   int incX = 1;
   int incY = -1;
   MpIeee expected=  0.604725f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 1)");
  };


  {
   int N = 1;
   MpIeee alpha=  0.1f;
   MpIeee X[] =  { 0.733f };
   MpIeee Y[] =  { 0.825f };
   int incX = 1;
   int incY = -1;
   MpIeee expected=  0.704725f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 2)");
  };


  {
   int N = 1;
   MpIeee alpha=  1.0f;
   MpIeee X[] =  { 0.733f };
   MpIeee Y[] =  { 0.825f };
   int incX = 1;
   int incY = -1;
   MpIeee expected=  1.604725f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 3)");
  };


  {
   int N = 1;
   MpIeee alpha=  0.0f;
   MpIeee X[] =  { -0.812f };
   MpIeee Y[] =  { -0.667f };
   int incX = -1;
   int incY = 1;
   MpIeee expected=  0.541604f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 4)");
  };


  {
   int N = 1;
   MpIeee alpha=  0.1f;
   MpIeee X[] =  { -0.812f };
   MpIeee Y[] =  { -0.667f };
   int incX = -1;
   int incY = 1;
   MpIeee expected=  0.641604f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 5)");
  };


  {
   int N = 1;
   MpIeee alpha=  1.0f;
   MpIeee X[] =  { -0.812f };
   MpIeee Y[] =  { -0.667f };
   int incX = -1;
   int incY = 1;
   MpIeee expected=  1.541604f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 6)");
  };


  {
   int N = 1;
   MpIeee alpha=  0.0f;
   MpIeee X[] =  { 0.481f };
   MpIeee Y[] =  { 0.523f };
   int incX = -1;
   int incY = -1;
   MpIeee expected=  0.251563f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 7)");
  };


  {
   int N = 1;
   MpIeee alpha=  0.1f;
   MpIeee X[] =  { 0.481f };
   MpIeee Y[] =  { 0.523f };
   int incX = -1;
   int incY = -1;
   MpIeee expected=  0.351563f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 8)");
  };


  {
   int N = 1;
   MpIeee alpha=  1.0f;
   MpIeee X[] =  { 0.481f };
   MpIeee Y[] =  { 0.523f };
   int incX = -1;
   int incY = -1;
   MpIeee expected=  1.251563f;
   MpIeee f;
   f = cblas_sdsdot (N, alpha, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdsdot(case 9)");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.785f };
   MpIeee Y[] =  { -0.7f };
   int incX = 1;
   int incY = -1;
   MpIeee expected=  -0.5495f;
   MpIeee f;
   f = cblas_sdot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdot(case 10)");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.79 };
   MpIeee Y[] =  { -0.679 };
   int incX = 1;
   int incY = -1;
   MpIeee expected=  -0.53641;
   MpIeee f;
   f = cblas_ddot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, dbleps, "ddot(case 11)");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.474f, -0.27f };
   MpIeee Y[] =  { -0.144f, -0.392f };
   int incX = 1;
   int incY = -1;
   MpIeee expected[2] =  {-0.174096f, -0.146928f};
   MpIeee f[2];
   cblas_cdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotu(case 12) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotu(case 12) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.474f, -0.27f };
   MpIeee Y[] =  { -0.144f, -0.392f };
   int incX = 1;
   int incY = -1;
   MpIeee expected[2] =  {0.037584f, -0.224688f};
   MpIeee f[2];
   cblas_cdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotc(case 13) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotc(case 13) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.87, -0.631 };
   MpIeee Y[] =  { -0.7, -0.224 };
   int incX = 1;
   int incY = -1;
   MpIeee expected[2] =  {0.467656, 0.63658};
   MpIeee f[2];
   cblas_zdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotu(case 14) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotu(case 14) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.87, -0.631 };
   MpIeee Y[] =  { -0.7, -0.224 };
   int incX = 1;
   int incY = -1;
   MpIeee expected[2] =  {0.750344, -0.24682};
   MpIeee f[2];
   cblas_zdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotc(case 15) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotc(case 15) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.457f };
   MpIeee Y[] =  { 0.839f };
   int incX = -1;
   int incY = 1;
   MpIeee expected=  -0.383423f;
   MpIeee f;
   f = cblas_sdot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdot(case 16)");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.949 };
   MpIeee Y[] =  { -0.873 };
   int incX = -1;
   int incY = 1;
   MpIeee expected=  -0.828477;
   MpIeee f;
   f = cblas_ddot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, dbleps, "ddot(case 17)");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.852f, -0.045f };
   MpIeee Y[] =  { 0.626f, -0.164f };
   int incX = -1;
   int incY = 1;
   MpIeee expected[2] =  {0.525972f, -0.167898f};
   MpIeee f[2];
   cblas_cdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotu(case 18) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotu(case 18) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { 0.852f, -0.045f };
   MpIeee Y[] =  { 0.626f, -0.164f };
   int incX = -1;
   int incY = 1;
   MpIeee expected[2] =  {0.540732f, -0.111558f};
   MpIeee f[2];
   cblas_cdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotc(case 19) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotc(case 19) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.786, -0.341 };
   MpIeee Y[] =  { -0.271, -0.896 };
   int incX = -1;
   int incY = 1;
   MpIeee expected[2] =  {-0.09253, 0.796667};
   MpIeee f[2];
   cblas_zdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotu(case 20) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotu(case 20) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.786, -0.341 };
   MpIeee Y[] =  { -0.271, -0.896 };
   int incX = -1;
   int incY = 1;
   MpIeee expected[2] =  {0.518542, 0.611845};
   MpIeee f[2];
   cblas_zdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotc(case 21) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotc(case 21) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.088f };
   MpIeee Y[] =  { -0.165f };
   int incX = -1;
   int incY = -1;
   MpIeee expected=  0.01452f;
   MpIeee f;
   f = cblas_sdot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, flteps, "sdot(case 22)");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.434 };
   MpIeee Y[] =  { -0.402 };
   int incX = -1;
   int incY = -1;
   MpIeee expected=  0.174468;
   MpIeee f;
   f = cblas_ddot(N, X, incX, Y, incY);
   gsl_test_rel(f, expected, dbleps, "ddot(case 23)");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.347f, 0.899f };
   MpIeee Y[] =  { -0.113f, -0.858f };
   int incX = -1;
   int incY = -1;
   MpIeee expected[2] =  {0.810553f, 0.196139f};
   MpIeee f[2];
   cblas_cdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotu(case 24) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotu(case 24) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.347f, 0.899f };
   MpIeee Y[] =  { -0.113f, -0.858f };
   int incX = -1;
   int incY = -1;
   MpIeee expected[2] =  {-0.732131f, 0.399313f};
   MpIeee f[2];
   cblas_cdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], flteps, "cdotc(case 25) real");
   gsl_test_rel(f[1], expected[1], flteps, "cdotc(case 25) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.897, -0.204 };
   MpIeee Y[] =  { -0.759, 0.557 };
   int incX = -1;
   int incY = -1;
   MpIeee expected[2] =  {0.794451, -0.344793};
   MpIeee f[2];
   cblas_zdotu_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotu(case 26) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotu(case 26) imag");
  };


  {
   int N = 1;
   MpIeee X[] =  { -0.897, -0.204 };
   MpIeee Y[] =  { -0.759, 0.557 };
   int incX = -1;
   int incY = -1;
   MpIeee expected[2] =  {0.567195, -0.654465};
   MpIeee f[2];
   cblas_zdotc_sub(N, X, incX, Y, incY, &f);
   gsl_test_rel(f[0], expected[0], dbleps, "zdotc(case 27) real");
   gsl_test_rel(f[1], expected[1], dbleps, "zdotc(case 27) imag");
  };


}
