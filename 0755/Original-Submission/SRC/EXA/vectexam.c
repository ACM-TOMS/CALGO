/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>
#define abs(x) ((x >= 0) ? (x) : -(x))
#define maxabs(x,y) (((x)>abs(y)) ? (x) : abs(y))
#define TAG 1

#ifdef CLOCKS_PER_SEC
#define clockspeed 1.0/CLOCKS_PER_SEC
#else
#define clockspeed 1.e-6    //machine dependent
#endif
#include <time.h>
long t0, t1,t2,t3,t4,t5,t6,t7,t8,t9;
void main() {
  int n,i,oper,indep,depen,buf_size,maxlive,deaths;
  int tape_stats[11];
  cout << "number of independent variables = ?  \n";
  cin >> n;
  double **xp = new double*[n];
  double yp =1;                     // Undifferenciated double code 
  for (i=0;i<n;i++)
    {
    xp[i] = new double[3];
    xp[i][1]=0;
    xp[i][2]=0;
    }
  t0 = clock();
  for(i=0;i<n;i++)
    {
    *xp[i] = (i+1.0)/(2.0+i); 
    yp *= *xp[i]; 
     };                        
  t1 = clock();
  double yout=0;
  int keep=1;
  adouble* x;
  x = new adouble[n];
  trace_on(TAG,keep);
  adouble y;
  y =1;
  for(i=0;i<n;i++) {
    x[i] <<= xp[i][0]; 
    y *= x[i];
  } 
  y >>= yout;
  t2 = clock();
  cout<< yout <<" =? "<<yp<<" function values should be the same \n";
  delete x;
  trace_off();                        // Reading of tape statistics

  tapestats(TAG,tape_stats);
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
  double **r = new double*[1];
  r[0] = new double[1];
  r[0][0] = yp;
  double err;
  int degree;      
  double *z = new double[n];
  double **g = new double*[n];
  double* h = new double[n];
  double *ind = new double[n];
  for (i=0;i<n;i++)
    g[i] = new double[1];
  for(int it=0; it <2; it++) {        // To come back later with perturbed xp
    degree =0;

    t3 = clock();
    forward(TAG,1,n,degree,keep,xp,r);
    t4 = clock();
    reverse(TAG,1,n,degree,1.0,g);        // Reverse sweep to evaluate gradient
    t5 = clock();
    err=0;
    for(i =0;i<n;i++) {               // Compare with deleted product
      err = maxabs(err,xp[i][0]*g[i][0]/r[0][0] - 1.0);
      ind[i] = xp[i][0];
    }
    cout<< err  <<" = maximum relative errors in gradient \n";
    // Combine previous two sweeps in gradient evaluation
    t6 = clock();
    gradient(TAG,n,ind,z);      //last argument lagrange is ommitted
    t7 = clock();
    err = 0;
    // Compare with previous numerical result
    for( i=0; i<n; i++) 
      {
	err =  maxabs(err,g[i][0]/z[i] - 1.0);
//      printf("g[%d][0] = %f  z[%d] = %f \n",i,float(g[i][0]),i,float(z[i]));
      }
    cout << err <<" = gradient error should be exactly zero \n";
    double *tan = new double[n];      // Compute first row of Hessian
    for(i=1;i<n;i++) tan[i] = 0.0 ;
    tan[0]=1.0; 
    t8 = clock();
    hess_vec(TAG,n,ind,tan,h);       // Computes Hessian times direction z.
    t9 = clock();
    err = abs(h[0]);
    //Compare with doubly deleted product
    for(i=1;i<n;i++) err = maxabs(err,xp[0][0]*h[i]/g[i][0]-1.0);
    cout<< err  <<" = maximum relative error in Hessian column \n";
    double h1n = h[n-1];                // Check for symmetry
    xp[0][1]=tan[0]=0;
    xp[n-1][1]=tan[n-1]=1;
    hess_vec(TAG,n,ind,tan,h);       // Computes Hessian times direction z.
    cout<<h1n<<" = "<<h[0]<<" (1,n) and (n,1) entry should be the same\n";
    // Third derivative tensor stuff
    xp[0][1]=1.0;
    degree = 2;
    double **rv;
    rv = new double*[1];
    rv[0] = new double[degree+1];
    double **t;
    t = new double*[n];
    for (i=0;i<n;i++)  t[i] =new double[degree+1];
  }
  if( t1-t0) 
    {double rtu = 1.0/(t1-t0);   
    cout << "\n times for "
	 <<"\n unitime:   " << (t1-t0)*clockspeed << "  seconds "
         <<"\n tracing:   "<< (t2-t1)*rtu << "   units "
         <<"\n forward:   "<< (t4-t3)*rtu  << "   units "
         <<"\n reverse:   "<<(t5-t4)*rtu   << "   units "
         <<"\n gradeval:   "<<(t7-t6)*rtu  << "   units "
         <<"\n hesseval:   "<<(t9-t8)*rtu  << "   units "
         <<"\n";
    }
  else cout << "\n-> zero timing due to small problem  dimension \n";
}
