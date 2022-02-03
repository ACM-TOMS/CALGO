/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include "adouble.h"
#include "adutils.h"

#ifdef __GNUG__
#include <std.h>
#include <builtin.h>
#else
#include <stdlib.h>
#endif

#include <stdio.h>
#include <stream.h>
#include <time.h>
int n,it;
double** PA;
double pdet(int k, int m) 
  {
    if(m == 0 ) return 1.0 ;
    else
    {
    double* pt = PA[k-1];
    double t=0 ;
    int p =1;
    int s;
    if (k%2) s = 1;
    else s = -1;
    for(int i=0;i<n;i++)
      {
       int p1 = 2*p;
       if ( m%p1 >= p )
	 {
          t += *pt*s*pdet(k-1, m-p);
          s = -s;
	  }
       ++pt;
       p = p1;
     } 
   return t;
    }
  }

adouble** A;
adouble det(int k, int m) 
  {
    if(m == 0 ) return 1.0 ;
    else
    {
    adouble* pt = A[k-1];
    adouble t=0 ;
    int p =1;
    int s;
    if (k%2) s = 1;
    else s = -1;
    for(int i=0;i<n;i++)
      {
       int p1 = 2*p;
       if ( m%p1 >= p )
	 {
          t += *pt*s*det(k-1, m-p);
          s = -s;
	  }
       ++pt;
       p = p1;
     } 
   return t;
    }
  }
void main()
{
   int i;
   int tag = 1;
   printf("order of matrix = ? \n",n);
   scanf("%d",&n);
   A = new adouble*[n];
   PA = new double*[n];
   int n2 =n*n;
   double* a = new double[n2];
   double diag;
   diag = 0;
   int m=1;
   double t00 = myclock();
   trace_on(tag,m);
   int loc =0;
   for (i=0; i<n; i++) 
      {
       m *=2;
       A[i] = new adouble[n];
       PA[i] = new double[n];
//     printf("\n %d  ",i);
       adouble* pt = A[i];
       double* ppt = PA[i];
       for (int j=0;j<n; j++)
	  {
	  *pt++ <<= j/(1.0+i);
	  *ppt++ = j/(1.0+i);
	  a[loc++] = j/(1.0+i);
//	  printf("%f  ",value(A[i][j]));
	  }
       diag += value(A[i][i]);   // val corrected to value 2/23/91
       A[i][i] += 1.0; 
       PA[i][i] += 1.0; 
       }
    diag += 1;
    adouble deter; 
    deter = det(n,m-1);
    double detout=0.0;
    deter >>= detout;
    printf("\n %f =? %f should be the same \n",detout,diag);
    trace_off();
    double t12 = myclock();
    int itu; itu=8-n; itu=itu*itu*itu*itu;
    itu = itu > 0 ? itu : 1;
    double raus;
    for(it = 0; it < itu; it++)
       raus = pdet(n,m-1);
    double t13 = myclock();
    double rtu = itu/(t13-t12);
//    cout << itu <<" reps -> time "<< 1/rtu <<"\n";
//    cout << t0-tm1+t1-t01 <<"= file time \n";
    double* B = new double[n2];
    double* detaut = new double[1];
    double t11 = myclock();
    for(it = 0; it < itu; it++)
    forward(tag,1,n2,0,1,a,detaut);
    double t21 = myclock();
    double u[1];
    u[0] = 1.0;
    for(it = 0; it < itu; it++)
    reverse(tag,1,n2,1,u,B);
    double t31 = myclock();
//  for(i=0;i<n2;i++) cout<<"xbar["<<i<<"][0] --> "<<B[i][0]<<"\n";
    cout << "\n first base? :   \n";
    for (i=0;i<n;i++) 
      {
      adouble sum = 0;
      adouble* pt;
      pt = A[i];
      for (int j=0;j<n;j++)
         sum += (*pt++)*B[j];
      cout << value(sum) <<" ";
      }
      cout << "\n";
      cout << "\n times for "
      <<"\n tracing: \t"<< (t11-t00)*rtu  
      <<" units \t" << (t11-t00) << "    seconds "
      <<"\n forward: \t"<< (t21-t11)*rtu/itu 
      <<" units \t" << (t21-t11)/itu << "    seconds "
      <<"\n reverse: \t"<<(t31-t21)*rtu/itu
      <<" units \t" << (t31-t21)/itu << "    seconds \n";
}
