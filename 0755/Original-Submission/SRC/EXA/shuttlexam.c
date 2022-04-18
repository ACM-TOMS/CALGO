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
#else
#include <stdlib.h>
#endif

#define tag 3
double getunitime(){
double t0 = myclock();
for(int it=0;it<5000;it++)
  {
  double f[7];
  double H,x,l,V,g,A,b,Hp,xp,lp,Vp,gp,Ap,bp,
          r,gr,rho,L,cd,ma,Om,Z;
  // *** Initialization of independent variables
  H = 264039.328;
  x = 177.718047;
  l = 32.0417885;
  V = 24317.0798;
  g = -0.749986488;
  A = 62.7883367;
  b = 41.100771834;
  Hp = -318;
  xp = 0.01;
  lp = 0.1;
  Vp = -3.6;
  gp = 0.001;
  Ap = 0.1;
  bp =0.06*(it+1);
  double ae = 20902900;
  double mu = 0.14E+17;
  bp /= (it+1);
  r = H+ae;
  gr = mu/(r*r);
  rho = .002378*exp(-H/23800.);
  double a = 40;
  double S =2690;
  double crtd = 180/3.14;
  double  cl = .84-.48*(38.-a*crtd)/26.;
  L = .5*rho*cl*S*V*V;
  cd = .78-.58*(38.-a*crtd)/26.;
  ma = 5964.496499824;
  Om = .72921159e-4;
  Z = .5*rho*cd*S*V*V;
  double C0 = 3.974960446019;
  double C1 = -.01448947694635;
  double C2 = -.2156171551995e-4;
  double C3 = -.1089609507291e-7;
  double V0 = 0;
// evaluate the dynamic equations ...
  double sing,cosg,sinA,cosA,sinl,cosl,tanl; 
  sinA =sin(A);
  cosA = cos(A);
  sing =sin(g);
  cosg = cos(g);
  sinl=sin(l);
  cosl=cos(l);
  tanl = sinl/cosl;
  f[0] = V*sing-Hp;
  f[1] = V*cosg*sinA/(r*cosl)-xp;
  f[2] = V*cosg*cosA/r-lp;
  f[3] = -Z/ma-gr*sing-Om*Om*r*cosl
         *(sinl*cosA*cosg-cosl*sing)-Vp;
  f[4] = L*cos(b)/(ma*V)+cosl/V*(V*V/r-gr)
         +2*Om*cosl*sinA
	 +Om*Om*r*cosl/V*(sinl*cosA*sing+cosl*cosg)
	 -gp;
  f[5] = L*sin(b)/(ma*V*cosg)+V/r*cosg*sinA*tanl
         -2*Om*(cosl*cosA*sing/cosg-sinl)
	 +Om*Om*r*cosl*sinl*sinA/(V*cosg)-Ap;
  f[6] = Z/ma - (C0+(V-V0)*(C1+(V-V0)*(C2+(V-V0)*C3)));
  }
  double t1 = myclock();
  double ti = (t1-t0)/5000.0;
  return ti ;
}

void main() {
  double rtu = 1.0/getunitime();
  int i,j,k,deg;
  adouble f[7];
  adouble H,x,l,V,g,A,b,Hp,xp,lp,Vp,gp,Ap,bp,
          r,gr,rho,L,cd,ma,Om,Z;
  cout<<"\nenter the degree: "; cin>>deg; cout<<"\n";
  trace_on(tag);
  // *** Initialization of independent variables
  H <<= 264039.328;
  x <<= 177.718047;
  l <<= 32.0417885;
  V <<= 24317.0798;
  g <<= -0.749986488;
  A <<= 62.7883367;
  b <<= 41.100771834;
  Hp <<= -318;
  xp <<= 0.01;
  lp <<= 0.1;
  Vp <<= -3.6;
  gp <<= 0.001;
  Ap <<= 0.1;
  bp <<=0.06;
  double ae = 20902900;
  double mu = 0.14E+17;
  r = H+ae;
  gr = mu/(r*r);
  rho = .002378*exp(-H/23800.);
  double a = 40;
  double S =2690;
  double crtd = 180/3.14;
  double  cl = .84-.48*(38.-a*crtd)/26.;
  L = .5*rho*cl*S*V*V;
  cd = .78-.58*(38.-a*crtd)/26.;
  ma = 5964.496499824;
  Om = .72921159e-4;
  Z = .5*rho*cd*S*V*V;
  double C0 = 3.974960446019;
  double C1 = -.01448947694635;
  double C2 = -.2156171551995e-4;
  double C3 = -.1089609507291e-7;
  double V0 = 0;
// evaluate the dynamic equations ...
  adouble sing,cosg,sinA,cosA,sinl,cosl,tanl; 
  sinA =sin(A);
  cosA = cos(A);
  sing =sin(g);
  cosg = cos(g);
  sinl=sin(l);
  cosl=cos(l);
  tanl = sinl/cosl;
  f[0] = V*sing-Hp;
  f[1] = V*cosg*sinA/(r*cosl)-xp;
  f[2] = V*cosg*cosA/r-lp;
  f[3] = -Z/ma-gr*sing-Om*Om*r*cosl
         *(sinl*cosA*cosg-cosl*sing)-Vp;
  f[4] = L*cos(b)/(ma*V)+cosl/V*(V*V/r-gr)
         +2*Om*cosl*sinA
	 +Om*Om*r*cosl/V*(sinl*cosA*sing+cosl*cosg)
	 -gp;
  f[5] = L*sin(b)/(ma*V*cosg)+V/r*cosg*sinA*tanl
         -2*Om*(cosl*cosA*sing/cosg-sinl)
	 +Om*Om*r*cosl*sinl*sinA/(V*cosg)-Ap;
  f[6] = Z/ma - (C0+(V-V0)*(C1+(V-V0)*(C2+(V-V0)*C3)));
  // *** pass final values of active variables to passive variables
  double res[7];
  for(i=0;i<7;i++)
    f[i] >>= res[i];
  trace_off ();                  
  // Set up correspondence between user's Variables and program Variables
  const int m = 7;
  const int n = 14;
  double **X = myalloc(n,deg+1);     // independent variable values
  double **Y = myalloc(m,deg+1);    // output Taylor coefficients
  double **Y0 = myalloc(m,deg+1);   // output Taylor coefficients
  X[0][0] = value(H);
  X[1][0] = value(x);
  X[2][0] = value(l);
  X[3][0] = value(V);
  X[4][0] = value(g);
  X[5][0] = value(A);
  X[6][0] = value(b);
  X[7][0] = value(Hp);
  X[8][0] = value(xp);
  X[9][0] = value(lp);
  X[10][0] = value(Vp);
  X[11][0] = value(gp);
  X[12][0] = value(Ap);
  X[13][0] = value(bp);
  srand(53);
  for(i=1;i<=deg;i++)
    for(j=0;j<n;j++)
      X[j][i] = sin(rand());
  double *u = new double[m];        // weighting vector for scalar reverse 
  for (i=0;i<m;i++)
    u[i] = 0.0;
  double **S_rev = myalloc(n,deg+1); // result adjoints from scalar reverse
  double ***V_rev = myalloc(m,n,deg+1);
  double ***W_rev = myalloc(m,n,deg+1);
  int rep, repu;
  repu = (int)(100.0/(deg+1));
  double epselon = 1.0e-3;
 // time divided differencing using forward
   double t0 = myclock();
   forward (tag,m,n,deg,0,X,Y);
   for (i=0;i<n;i++) {             
      X[i][0] *= 1+epselon;
      for (rep =0;rep < repu ;rep ++)
          forward (tag,m,n,deg,0,X,Y0);  // may wish to store columns of matrix :
      double recx = 1.0/(epselon*X[i][0]);
      for(k=0;k<=deg;k++)
	  for(j=0;j<m;j++)
              W_rev[j][i][k] = (Y0[j][k]-Y[j][k])*recx ;
      X[i][0] /= 1+epselon; }         // (Y[1..m][deg] - yofx[1..m]) / epselon
    double t1 = myclock();
    forward(tag,m,n,deg,deg+1,X,Y);       // in preparation for reverse sweeps
// time single vector reverse call
    double t2 = myclock();
    short** nonzero = new short*[7];
    for (int ii=0;ii<7;ii++)
       nonzero[ii] = new short[14];
    for (rep =0;rep < repu ;rep ++)
    reverse (tag,m,n,deg,V_rev,nonzero);// matrix is stored as :
    double t3  = myclock();
    cout <<" Print out the nonzero pattern? \n";
    int yes; cin >> yes;
    if(yes)
    {
  cout<<" 4 = transcend , 4 = rational , 2 = polynomial , 1 = linear , 0 = zero \n";
    for(i=0;i<7;i++)
      {
      for(j=0;j<14;j++)
         cout << nonzero[i][j] <<"   ";
       cout <<"\n";
       }
     }
// time multiple calls to scalar reverse
    float err =0;
    cout << "Print comparison between reverse and differences? \n"; 
    cin >> yes;
    double t4 = myclock();
    for (i=0;i<m;i++) {
      u[i] = 1;                     
    for (rep =0;rep < repu ;rep ++)
      reverse (tag,m,n,deg,u,S_rev);// rows of matrix are stored as :
      u[i] = 0;                    // S_rev[1..n][deg]
      for ( j=0;j<n;j++) 
	 for (int k=0;k<=deg;k++)
	   {
	   double lerr= fabs(S_rev[j][k]-V_rev[i][j][k]);
	   if (lerr != 0.0)
	   err += lerr/ (fabs(S_rev[j][k])+fabs(V_rev[i][j][k]));
	   if(yes && (S_rev[j][k] != 0 || W_rev[i][j][k] != 0) )
 	   cout << S_rev[j][k] << " =?= " << W_rev[i][j][k] << "\n";
	   }
    }
    double t5 = myclock();
    cout<<"Jacobian discrepancy  : "<<err<<"  between reverse modes\n";
    cout << "forward-eval:  " << rtu*(t1-t0)/(repu*n) << 
    " units or  "<<(t1-t0)/(repu*n) <<"  seconds \n";
    cout << "forward-diff:  " << rtu*(t1-t0)/repu <<
    " units or  "<< (t1-t0)/repu <<"  seconds \n";
    cout << "vect-reverse:  " << rtu*(t3-t2)/repu <<
    " units or  "<< (t3-t2)/repu <<"  seconds \n";
    cout << "scal-reverse:  " << rtu*(t5-t4)/repu <<
    " units or  " << (t5-t4)/repu <<"  seconds \n";
}
