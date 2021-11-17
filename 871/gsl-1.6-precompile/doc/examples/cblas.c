#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_cblas.h>

int
main (void)
{
  int  lda=  3;

  MpIeee A[] =  { MpIeee( "0.11" ), MpIeee( "0.12" ), MpIeee( "0.13" ),
                MpIeee( "0.21" ), MpIeee( "0.22" ), MpIeee( "0.23" ) };

  int  ldb=  2;
  
  MpIeee B[] =  { MpIeee( "1011" ), MpIeee( "1012" ),
                MpIeee( "1021" ), MpIeee( "1022" ),
                MpIeee( "1031" ), MpIeee( "1032" ) };

  int  ldc=  2;

  MpIeee C[] =  { MpIeee( "0.00" ), MpIeee( "0.00" ),
                MpIeee( "0.00" ), MpIeee( "0.00" ) };

  /* Compute C = A B */

  cblas_sgemm (CblasRowMajor, 
               CblasNoTrans, CblasNoTrans, 2, 2, 3,
               1.0, A, lda, B, ldb, 0.0, C, ldc);

  {cout<<"[ "<<setiosflags((ios::floatfield))<< C[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< C[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"  "<<setiosflags((ios::floatfield))<< C[2];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< C[3];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" ]\n";}

  return 0;  
}
