// Public-domain work without warranty — see README:LICENSING for details.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int acm409(int m, int n, int k,          // H. Schmitt’s approx()
    double const x[], double const y[],
    double epsh, int * maxit, int ref[],
    double * hmax, double h[], double a[]);

static double const pi=M_PI, eps=1e-12;

typedef double Interval[2];



static void acm409f(int k, int m0, Interval di[], int ndi, double f(double),
    int num, char const * name, int print, FILE * f_h)
// Wrapper around acm409() that allows abscissae x[] and ordinates y[] to be
// expressed instead as a number ndi of domain intervals di[] of ℝ, and a
// continuous function f() on that domain.  m is given as abs(m0); if m0 is
// negative then a supplied linearly-spaced initial reference replaces the
// default (Chebyshev abscissae).  Diagnostics and results are output via num,
// name, print, and f_h.
{
  int i, j, l, m=abs(m0), n=0, N[ndi], ref[1+m+2], maxit=19, ret;
  double d=0, a[1+m], hmax=HUGE_VAL;

  // Determine N[], the number of evaluations of f per interval (less 1):
  for (i=0; i<ndi; ++i) d += di[i][1]-=di[i][0];
  for (i=0; i<ndi; ++i)
    n += 1 + (N[i] = (int)round(di[i][1]/d*(1024<<!k)*(m+1)));

  // Evaluate abscissae x[] and ordinates y[]:
  double x[1+n], y[1+n], h[1+n];
  for (l=i=0; i<ndi; ++i) for (j=0; j<=N[i]; ++j, ++l)
    y[1+l] = f(x[1+l]=di[i][0]+1.*j/N[i]*di[i][1]);

  // Evaluate the initial reference set (used only if m0<0):
  for (i=0; i<=m+1; ++i) ref[1+i]=1+(int)round(1.*i/(m+1)*(n-1));

  ret=acm409(m, n, k, x, y, m0<0? -eps:eps, &maxit, ref, &hmax, h, a);

  // Output summary ...
  static char const * const exits[] = {" - ", "MAX", "SIG", "PAR"};
  printf("%2i  %s %2i  %10.4e  %i  %2i  %-21s",
      num, exits[ret], maxit, hmax, k, m0, name);
  for (i=0; i<ndi; ++i) printf("  %g:%g", di[i][0], di[i][0]+di[i][1]);
  putchar('\n');

  // ... and optionally, other results:
  if (ret<3) {
    double norm = hmax==0? 1 : 1/hmax; // To normalise the error.
    if (f_h)   for (i=1; i<=n; ++i) fprintf(f_h, "%23.17g\n", h[i]*norm);
    if (print) for (i=m; i>=0; --i)  printf(     "%23.17g\n", a[i]);
  }
}



#ifdef __cplusplus                     // C++11 lamda:
  #define F(e) [](double x){return (double)(e);}
#else                                  // gcc lamda:
  #define F(e) ({double f(double x){return e;}f;})
#endif
#define rm_paren(...) __VA_ARGS__
#define len(x) sizeof(x)/sizeof((x)[0])

#define test(k,m,e,di) if ((i+=inc)==arg) {Interval D[]={rm_paren di};\
  acm409f(k,m,D,len(D),F(e),num,#e,inc,f_h);} ++num



int main(int argc, char * argv[])
// Exercises acm409() approximation with a selection of functions.
{
  int    arg = argc<=1? -1 : atoi (argv[1]), i=-1, inc=arg>=0, num=0;
  FILE * f_h = argc<=2?  0 : fopen(argv[2], "w");

  puts(" # exit its    hmax    "
      " k   m  Expression for f(x)    Domain intervals");
  test( 0,  3, (x-1)*(x-2)*(x-3)    ,( {  0, 4} ));
  test( 1,  9, sin(pi/2*x)          ,( {  0, 1} ));
  test( 0, 17, log2(x)              ,( {  1, 2} ));
  test( 0, 31, sin(exp(3*x))        ,( { -1, 1} ));
  test( 0,  6, sqrt(x)              ,( {  1, 2} ));
  test( 0, -6, sqrt(x)              ,( {  1, 2} ));
  test( 0,  8, cos(x)+cos(21*x+1)   ,( { -1, 1} ));
  test( 0, 10, 1-sin(5*fabs(x-.5))  ,( { -1, 1} ));
  test( 0, 15, x>.4                 ,( { -1, .3}, { .5, 1} ));
  test( 0,1+3, asin(x)              ,( { -1,-.9}, { .9, 1} ));
  test( 1,  3, asin(x)              ,( {.9, 1} ));
  test( 1,  7, erf(x)               ,( {eps, 1} ));
  test( 1,  9, sin(pi/2*x)          ,( {eps, 1} ));
  test( 2, 10, fabs(x)              ,( {  0, 1} ));
  test( 2, 14, x<.7                 ,( {  0, .6}, { .8, 1} ));
  test( 2, 12, asin(cos(x*pi*3))    ,( {  0,  1} ));
  return 0;
}
