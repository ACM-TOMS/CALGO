/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include <stream.h>
#include <stdio.h>

#ifdef __GNUG__
#include <std.h>
#else
#include <stdlib.h>
#endif

/*  This program can be used to verify the consistency and correctness
of derivatives computed by ADOL-C in its forward and eevere mode.  
The use is required to selct one integer input id. For positive n = id
the monomial x^n is evaluated recursively at x=0.5 and all its nonzero
Taylor coeffcients at this point are evaluated in the forward and
reverse mode. A negative choice of id >= -9 leads to one of nine
identities, whose derivatives should be trivial. These identities
may be used to check the correctness of particular code segments
in the ADOL-C sources forward.c and reverse.c. No timings are
performed in this example program */

#include "adouble.h"    // These includes provide the compiler with
#include "adutils.h"    // definitions and utilities for `adoubles'.

//   The monomial evaluation routine which has been obtained from
//   the original version by retyping all `doubles' as `adoubles'.

adouble power(adouble x, int n)
  {
  adouble z=1;
  if (n>0)
    {
     int nh =n/2;
     z = power(x,nh);
     z *= z;
     if (2*nh != n) z *= x;
     return z;
     }
  else
    {
     if (n==0) 
       return z;
     else return 1.0/power(x,-n);
    }
  }
void main() 
{
  int n,i,id;
  int tag = 0;
  cout << "problem number(-1 .. -10) / degree of monomial =? \n";
  cin >> id;
  n = id >0 ? id : 3;
  double *xp,*yp;
  xp = new double[n+4];
  yp = new double[n+4];
  yp[0]=0;
  xp[0] = 0.5;
  xp[1] = 1.0;
  adouble y,x; 
  int dum=1; 
  trace_on(tag,dum);   // Begin taping all calculations with 'adoubles'
  x <<= xp[0];
  if( id >= 0 )
    { 
    cout << "Evaluate and differentiate recursive power routine \n";
    y = power(x,n);
    }
  else
    {
    cout<<"Check Operations and Functions by Algebraic Identities \n";
    double pi = 2*asin(1.0);
    switch (id)
      {
      case -1 : 
	 cout << "Addition/Subtraction: y = x + x - (2.0/3)*x - x/3 \n";
	 y =  x + x - (2.0/3)*x - x/3 ;
	 break;
      case -2 : 
         cout << "Multiplication/divison:  y = x*x/x \n";
         y = x*x/x;
	 break;
      case -3 :
	 cout << "Square root and power: y = sqrt(pow(x,2)) \n"; 
	 y = sqrt(pow(x,2));
	 break;
      case -4 :
	 cout << "Exponential and log: y = exp(log(log(exp(x)))) \n";
	 y = exp(log(log(exp(x))));
	 break;
      case -5 :
	 cout << "Trig identity: y = x + sin(2*x)-2*cos(x)*sin(x) \n";
	 y =  x + sin(2.0*x)-2.0*cos(x)*sin(x);
	 break;
      case -6 :
         cout << "Check out quadrature macro \n";
	 y = exp(myquad(myquad(exp(x))));
	 break;
      case -7 :
	 cout << "Arcsin: y = sin(asin(acos(cos(x)))) \n";
	 y = sin(asin(acos(cos(x))));
	 break;
      case -8 :
	 cout << "Hyperbolic tangent: y = x + tanh(x)-sinh(x)/cosh(x) \n";
	 y = x + tanh(x)-sinh(x)/cosh(x) ;
	 break;
      case -9 :
	 cout << "Absolute value: y = x + fabs(x) - fabs(-x) \n";
	 y = x + fabs(-x) - fabs(x);
         break;
      case -10 :
	 cout << "atan2: y = atan2(sin(x-0.5+pi),cos(x-0.5+pi)) \n";
	 y = atan2(sin(x),cos(x));
         // y= tan(atan2(x,1.0));
	 break;
      default : cout << " Please select problem number >= -10 \n";
		exit(-1);
      }
    cout << "Round-off error: " << value(y-x)  << " \n";
    }
  y >>= yp[0];
  trace_off();  // The (partial) execution trace is completed.
  int oper,indep,depen,buf_size,maxlive,deaths;
  int tape_stats[11];
  tapestats(tag,tape_stats);
  indep = tape_stats[0];
  depen = tape_stats[1];
  buf_size = tape_stats[4];
  oper = tape_stats[5];
  deaths = tape_stats[3];
  maxlive = tape_stats[2];

  
  cout<<"\n";
  cout<<"independents "<<indep<<"\n";
  cout<<"dependents   "<<depen<<"\n";
  cout<<"operations   "<<oper<<"\n";
  cout<<"buffer size  "<<buf_size<<"\n";
  cout<<"maxlive      "<<maxlive<<"\n";
  cout<<"deaths       "<<deaths<<"\n";
  double *res;
  res = new double[n+2];
  double u[1]; 
  u[0]=1;
 cout << "\nThe two Taylor coefficients in each row should agree\n\n";

 double ***V = (double***)new double**[1];
 V[0] = new double*[1];
 V[0][0] = new double[n+2];
 double **U = new double*[1]; 
 U[0] = new double[1];
 U[0][0] = 1;
 double** xpoint = &xp;
 double** ypoint = &yp;
 double** respoint = &res;
 cout << " \n \t   forward  \t    reverse  \n";
 for( i=0; i < n+2; i++)
   {  
   xp[i+2]=0;    
   forward(tag,depen,indep,i,i+1,xpoint,respoint);
   cout<< i <<"\t"<< res[i] <<"\t \t "<< yp[i] << "\n";
   reverse(tag,depen,indep,i,u,ypoint); // call higher order scalar reverse
   reverse(tag,depen,indep,i,1,U,V);
   yp[i+1] = yp[i]/(i+1);
   if(V[0][0][i] != yp[i])
     cout << i <<"-th component in error "<<V[0][0][i]-yp[i] <<"\n";
   }
 cout << "\nWhen n<0 all rows except the first two should vanish \n";
 }
