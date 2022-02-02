/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>
#include <stdio.h>

#ifdef __GNUG__
#include <std.h>
#include <builtin.h>
#else
#include <stdlib.h>
#endif

int const N=3;     // State Space Dimesion of ODE

/**********************************************************
   This orinigal cersion of the right hand side may be used
   for run time comparisons. However, the results are not 
   representative because the function is extremely small
   and is likely to stay in cache if evaluated repeatedly.
*********************************************************/

void tracerhs(double* y, double* yprime)
{
yprime[0] = -sin(y[2]) + 1e8*y[2]*(1-1/y[0]);
yprime[1] = -10*y[0] + 3e7*y[2]*(1-y[1]);
yprime[2] = -yprime[0] - yprime[1];
}

/********************************************************
   This overloaded function generates the tape by evaluating
   the right hand side in terms of adouble variables.
*********************************************************/
void tracerhs(unsigned int tag, double* py, double* pyprime)
{
// py must represent some feasible state space point
adouble y[N],yprime[N];
int i;
trace_on(tag);

for (i=0;i<N;i++)
  y[i] <<= py[i];    //  Initialize and mark independents
yprime[0] = -sin(y[2]) + 1e8*y[2]*(1-1/y[0]);
yprime[1] = -10*y[1] + 3e7*y[2]*(1-y[1]);
yprime[2] = -yprime[0] - yprime[1];
for (i=0;i<N;i++)
   yprime[i] >>= pyprime[i];   // Mark and pass dependents

trace_off();
}

void main()
{
/*******************************************************
   Select problem set up data and generate the tracerhs
*******************************************************/
cout << "highest derivatives =? \n";
int D;
cin >> D;
int i,j,k;
int yes ;
double **z, *zz, *pyp;
z = new double*[N];
zz = new double[N];
pyp = new double[N];
for (i=0;i<N;i++)
  {
  z[i] = new double[D+1];
  if(i==0) {*z[i] = 1;} else *z[i] =0.0;
  zz[i] = *z[i];
  }
for (i=0;i<N;i++)
  zz[i] = *z[i];
double*** B = myalloc(N,N,D);
double*** A = myalloc(N,N,D);
/******************Compute the unit time for comparison**********/
double t0 = myclock();
for(j=0;j<10000;j++)
{ zz[1] = 1.0/j;
tracerhs(zz,pyp);}
double t1 = myclock();
double rtu = 10000.0/(t1-t0);
/******************Do the actual taping ********************/
unsigned short tag=0;
tracerhs(tag,zz,pyp);

/*******************************************************
 Check out forode and jacode for consistency on tape
*******************************************************/

double tau;
cout <<" Enter nonzero scaling paramater \n ";
cin >> tau;
double **w = myalloc(N,D+1);;
short** nonzero;
nonzero = new short*[N];
for (i=0;i<N;i++)
  nonzero[i] = new short[N];
double ten = 200 ;// D < 50 ? 50 : D;
t0 = myclock();
for(j=0;j<ten;j++)
forode(tag,N,tau,D,z);
t1 = myclock();
for(j=0;j<ten;j++)
reverse (tag,N,N,D-1,A,nonzero);
double t2 = myclock();
for(j=0;j<ten;j++)
accode(N,tau,D-1,A,B,nonzero);
double t3 = myclock();
cout <<" Print out the nonzero pattern? \n";
cin >> yes;
if(yes)
{
cout << " 4 = transcend , 3 = rational , 2 = polynomial , 1 = linear , 0 = zero \n";
cout << " negative number k indicate that entries of all B_j with j < -k vanish  \n";
for(i=0;i<N;i++)
   {
   for(j=0;j<N;j++)
      cout << nonzero[i][j] <<"   ";
   cout <<"\n";
   }
}
/****************************************************
  The D+1 columns of z should now be consistent with the
  ODE as represented by the time. Feeding z into forward
  we obtain a coeffient array w, whose columns should
  equal to the shifted and scaled columns of z
*****************************************************/
cout <<"\nCheck that forward reproduces the Taylor series \n";
forward(tag,N,N,D-1,D,z,w);
double err = 0;
for(i=0;i<D;i++)
  {
   for(j=0;j<N;j++)
     {
     double avg = (fabs(w[j][i]*tau)+fabs(z[j][i+1]*(i+1)))/2;
     if (avg < 0.1) avg = 0.1;
     if(avg) err += fabs(w[j][i]*tau-z[j][i+1]*(i+1))/avg;
     }
  }
cout << err << "= total error \n";
t0 /=ten;  t1 /= ten; t2 /= ten; t3 /= ten;
cout << (t1-t0)*rtu<<" units or "<<t1-t0<<" seconds just forode \n";
cout << (t2-t1)*rtu<<" units or "<<t2-t1<<" seconds just reverse \n";
cout << (t3-t2)*rtu<<" units or "<<t3-t2<<" seconds just accode \n";
/*  If desired print out Jacobians of Taylor coeffcients with
   respect to the base point */
cout <<"\nPrint Jacobian of Taylor coefficient vectors? \n";
cin >> yes;
if(yes)
  {
  for(i=0;i<D;i++){
    cout<<"\n\t<-- B("<<i<<")\n";
    for(j=0;j<N;j++){
      for(k=0;k<N;k++)
        cout<<B[j][k][i]<<"  ";
      cout<<"\n";}}
  } 
cout <<"\n Enter increment for differencing, skipped if zero: \n"; 
double h;
cin >> h;
for (i=0;i<N;i++)
   *w[i] = *z[i];
if(h!=0)   
for (i=0;i<D;i++)
  {
  err = 0;
  for(k=0;k<N;k++)
    {
    *w[k] += h;
    forode(tag,N,tau,D,w);
    *w[k] -= h;
    for(j=0;j<N;j++)
      err += B[j][k][i] != 0 ? fabs(1-(w[j][i+1]-z[j][i+1])/h/B[j][k][i])
	    : fabs((w[j][i+1]-z[j][i+1])/h) ;
    };
cout<<"Relative truncation errors in B("<<i<<") ---> "<<err<<"\n"; 
};  
}
