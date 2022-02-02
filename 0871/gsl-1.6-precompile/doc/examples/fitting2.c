#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_multifit.h>

int
main (int argc, char **argv)
{
  int  i;int   n;
  MpIeee xi;MpIeee  yi;MpIeee  ei;MpIeee  chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  if (argc != 2)
    {
      fprintf (stderr,"usage: fit n < data\n");
      exit (-1);
    }

  n = atoi (argv[1]);

  X = gsl_matrix_alloc (n, 3);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);

  c = gsl_vector_alloc (3);
  cov = gsl_matrix_alloc (3, 3);

  for (i = 0; i < n; i++)
    {
      int  count=  fscanf (stdin, "%lg %lg %lg",
                          &xi, &yi, &ei);

      if (count != 3)
        {
          fprintf (stderr, "error reading file\n");
          exit (-1);
        }

      {cout<<""<<setiosflags((ios::floatfield))<< xi;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< yi;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" +/- "<<setiosflags((ios::floatfield))<< ei;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
      
      gsl_matrix_set (X, i, 0, 1.0);
      gsl_matrix_set (X, i, 1, xi);
      gsl_matrix_set (X, i, 2, xi*xi);
      
      gsl_vector_set (y, i, yi);
      gsl_vector_set (w, i, 1.0/(ei*ei));
    }

  {
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, 3);
    gsl_multifit_wlinear (X, w, y, c, cov,
                          &chisq, work);
    gsl_multifit_linear_free (work);
  }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    {cout<<"# best fit: Y = "<<setiosflags((ios::floatfield))<< 
            C(0);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" + "<<setiosflags((ios::floatfield))<<, C(1);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" X + "<<setiosflags((ios::floatfield))<<, C(2);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" X^2\n";}

    {cout<<"# covariance matrix:\n";}
    {cout<<"[ "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<
               COV(0,0);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(0,1);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(0,2);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<"  \n";}
    {cout<<"  "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<< 
               COV(1,0);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(1,1);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(1,2);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<"  \n";}
    {cout<<"  "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<< 
               COV(2,0);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(2,1);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<", "<<setiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<setprecision(5)<<, COV(2,2);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::scientific & ios::floatfield))<<" ]\n";}
    {cout<<"# chisq = "<<setiosflags((ios::floatfield))<< chisq;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  }
  return 0;
}
