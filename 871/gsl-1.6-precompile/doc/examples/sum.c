#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

#define N 20

int
main (void)
{
  MpIeee t[N];
  MpIeee sum_accel;MpIeee  err;
  MpIeee sum=  MpIeee( "0" );
  int  n;
  
  gsl_sum_levin_u_workspace * w 
    = gsl_sum_levin_u_alloc (N);

  const MpIeee zeta_2=  M_PI * M_PI / 6.0;
  
  /* terms for zeta(2) = \sum_{n=1}^{\infty} 1/n^2 */

  for (n = 0; n < N; n++)
    {
      MpIeee np1=  n + MpIeee( "1.0" );
      t[n] = MpIeee( "1.0" ) / (np1 * np1);
      sum += t[n];
    }
  
  gsl_sum_levin_u_accel (t, N, w, &sum_accel, &err);

  {cout<<"term-by-term sum = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< 
          sum;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" using "<< N<<" terms\n";}

  {cout<<"term-by-term sum = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< 
          w->sum_plain;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" using "<< w->terms_used<<" terms\n";}

  {cout<<"exact value      = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< zeta_2;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"accelerated sum  = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< 
          sum_accel;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" using "<< w->terms_used<<" terms\n";}

  {cout<<"estimated error  = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< err;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"actual error     = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<< 
          sum_accel - zeta_2;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  gsl_sum_levin_u_free (w);
  return 0;
}
