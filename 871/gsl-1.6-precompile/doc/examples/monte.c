#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

/* Computation of the integral,

      I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */

/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */

MpIeee exact=  MpIeee( "1.3932039296856768591842462603255" );

MpIeee g(MpIeee *k, size_t dim, void *params)
{
  MpIeee A=  MpIeee( "1.0" ) / (M_PI * M_PI * M_PI);
  return A / (MpIeee( "1.0" ) - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

void
display_results (char *title, MpIeee result, MpIeee error)
{
  {cout<<""<< title<<" ==================\n";}
  {cout<<"result = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<< result;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"sigma  = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<< error;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact  = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<< exact;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"error  = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<< result - exact;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" = "<<setiosflags((ios::floatfield))<<setprecision(1)<<
          fabs (result - exact);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" sigma\n";}
}

int
main (void)
{
  MpIeee res;MpIeee  err;

  MpIeee xl[3] =  { MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ) };
  MpIeee xu[3] =  { M_PI, M_PI, M_PI };

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = { &g, 3, 0 };

  size_t calls = 500000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
    gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s, 
                               &res, &err);
    gsl_monte_plain_free (s);

    display_results ("plain", res, err);
  }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    {cout<<"converging...\n";}

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                   &res, &err);
        {cout<<"result = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<<
                "chisq/dof = %.1f\n";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" sigma = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(6)<< res;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" ";}
      }
    while (fabs (s->chisq - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }
  return 0;
}
