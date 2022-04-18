/* CACM Alg 280 */
#include <stdio.h>
#include <math.h>
#define nmax 100
double pown(a,n)
  double a; int n;
{
  double x; int i;
  if(n == 0) return(1);
  if(n>0)
     {
       x=a;
       for(i=1;i<n;i++) x=x*a;
     }
  return(x);
}
void gregory(n,r,t,w)
  int n,r;
  double w[],t[];
/* 

  Computes the abscissas and weights of the Gregory quadrature rule with r
  differences:
  
    \int_{t_0}^{t_n} f(t) dt  \approx
  	  h \left( \frac{1}{2} f_0 + f_1 + \cdots +
            f_{n-1} + \frac{1}{2} f_n \right ) -
  	  \frac{h}{12}( \nabla f_n - \delta f_0) -
  	  \frac{h}{24}( \nabla^2 f_n + \delta^2 f_0) - \cdots -
  	  h c_{r+1}^{*} (\nabla^r f_n + \delta^r f_0)
  = \sum_{j=0}^{n} w_j f(t_j),
  
  where h = (t_n - t_0)/n, and the c_j^* are given in Henrici (1964). The
  number r must be an integer from 0 to n, the number of subdivisions. The
  left and right endpoints must be in t(0) and t(n) respectively. The
  abscissas are returned in t(0) to t(n) and the corresponding weights in
  w(0) to w(n).
  
  If r=0 the Gregory rule is the same as the repeated trapezoid rule, and if
  r=n the same as the Newton-Cotes rule (closed type). The order p of the
  quadrature rule is p = r+1 for r odd and p = r+2 for r even. For n >= 9
  and large r some of the weights can be negative.
  
  For n<= 32 and r<= 24, the numerical integration of powers (less than r)
  of x on the interval [0,1] gave 9 significant digits correct in an 11
  digit mantissa.
  
  Refs:
    Hildebrand, F. B. Introduction to Numerical Analysis. McGraw-Hill, New
    York, 1956, p. 155.
  
    Henrici, Peter. Elements of Numerical Analysis. Wiley, New York, 1964,
    p. 252.

*/
{
  int i,j; double h,cj,c[nmax+1],b[nmax]; 
  b[0]=1; c[0]=1; c[1]=-0.5; b[n]=0;
  h=(t[n]-t[0])/n;
  w[0]=w[n]=0.5;
  for(i=1;i<n;i++)
    {
      w[i]=1; t[i]=i*h+t[0]; b[i]=0;
    }
  if(r>n) r=n;
  for(j=1;j<=r;j++)
  {
    cj=0.5*c[j];
    for(i=j;i>=1;i--) b[i]=b[i]-b[i-1];
    for(i=3;i<=j+2;i++) cj=cj+c[j+2-i]/i;
    c[j+1]=-cj;
    for(i=0;i<=n;i++) w[i]=w[i]-cj*(b[n-i]+b[i]);
  }
  for(i=0;i<=n;i++) w[i]=w[i]*h;
}
main()
{
  double t[nmax],w[nmax],a=-1.0,b=1.0,I,In;
  int i,n,r,p;
  n=7; r=n;
  t[0]=a; t[n]=b;    /* limits of integration */
/* Generate the weights */
  gregory(n,r,t,w);
  for(i=0;i<=n;i++) printf("%f %f\n",w[i],t[i]);
/* Check - integrate x^p */
  for(p=0;p<=r+4;p++)
  {
    I=(pown(b,p+1)-pown(a,p+1))/(p+1);
    In=0;
    for(i=0;i<=n;i++) In=In+w[i]*pown(t[i],p);
    printf("%d %d %d %f %f %e\n",n,r,p,I,In,I-In);
  }
}
