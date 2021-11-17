#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<signal.h>
#include <complex.h>
#include <unistd.h>
#include <float.h>
static int imaxarg1,imaxarg2;
#define Sqr(x) ((x)*(x))
#ifndef CMPLXL
#define CMPLXL(x,y) ((x)+I*(y))
#endif
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define SIGNL(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
void balancel(long double a[4][4], int n)
{
  const long double RADIX=FLT_RADIX;// numeric_limits<Doub>::radix;
  int i, j;
  long double scale[4]={1.0,1.0,1.0,1.0};
  int done=0;
  long double r, c, g, f, s, sqrdx=RADIX*RADIX;
  while (!done) 
    {
      done=1;
      for (i=0;i<n;i++) 
	{
	  //Calculate row and column norms.
	  //If both are nonzero,
	  //find the integer power of the machine radix that comes closest to balancing the matrix.
	  r=0.0;
	  c=0.0;
	  for (j=0;j<n;j++)
	    if (j != i) 
	      {
		c += fabsl(a[j][i]);
		r += fabsl(a[i][j]);
	      }
	  if (c != 0.0 && r != 0.0) 
	    {
	      g=r/RADIX;
	      f=1.0;
	      s=c+r;
	      while (c<g) {
		f *= RADIX;
		c *= sqrdx;
	      }
	      g=r*RADIX;
	      while (c>g) 
		{
		  f /= RADIX;
		  c /= sqrdx; 
		}
	      if ((c+r)/f < 0.95*s) 
		{
		  done=0;
		  g=1.0/f;
		  scale[i] *= f;
		  for (j=0;j<n;j++) a[i][j] *= g; //Apply similarity transformation
		  for (j=0;j<n;j++) a[j][i] *= f;
		}
	    }
	}
    }
}
void hqrl(long double a[4][4], complex long double wri[4], int *ok, int n)
{
  int nn,m,l,k,j,its,i,mmin;
  long double z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0, anorm=0.0;
  const int MAXITS = 480;
  const long double EPS= 1.92592994438723585305597794258492732E-0034;
  for (i=0;i<n;i++)
    //Compute matrix norm for possible use in lo- cating single small sub diagonal element.
    for (j=IMAX(i-1,0);j<n;j++)
      anorm += fabsl(a[i][j]);
  nn=n-1;
  t=0.0;
  *ok = 1;
  //Gets changed only by an exceptional shift.
  while (nn >= 0) 
    {
      //Begin search for next eigenvalue.
      its=0;
      do 
	{
	  for (l=nn;l>0;l--)
	    {
	      //Begin iteration: look for single small sub di- agonal element.
	      s=fabsl(a[l-1][l-1])+fabsl(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	
	      if (fabsl(a[l][l-1]) <= EPS*s)
		{
		  a[l][l-1] = 0.0;
		  break;
		}
	    }
	  x=a[nn][nn];
	  if (l == nn)
	    {
	      //One root found.  
	      wri[nn--]=x+t;
	    } 
	  else
	    {
	      y=a[nn-1][nn-1];
	      w=a[nn][nn-1]*a[nn-1][nn];
	      if (l == nn-1)
		{
		  //Two roots found...
		  p=0.5*(y-x);
		  q=p*p+w;
		  z=sqrtl(fabsl(q));
		  x += t;
		  if (q >= 0.0)
		    {
		      //...a real pair.
		      z=p+SIGNL(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=CMPLXL(x+p,-z);
		      wri[nn-1]=conjl(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == MAXITS)// il valore era 30 ma l'ho aumentato a 200 altrimenti non ce la fa...boh
		    {
		      printf("LONG Too many iterations in hqr");
		      *ok = 0;
		      return;
		    }
		  /* 
		   * condition from GSL routine
		   */
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=fabsl(a[nn][nn-1])+fabsl(a[nn-1][nn-2]);
		      y=x=0.75*s;
		      w = -0.4375*s*s;
		    }
		  ++its;
		  for (m=nn-2;m>=l;m--)
		    {
		      //Form shift and then look for 2 consecutive small sub- diagonal elements.
		      z=a[m][m];
		      r=x-z;
		      s=y-z;
		      p=(r*s-w)/a[m+1][m]+a[m][m+1];
		      //Equation (W ebnote 16.21).
		      q=a[m+1][m+1]-z-r-s;
		      r=a[m+2][m+1];
		      s=fabsl(p)+fabsl(q)+fabsl(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=fabsl(a[m][m-1])*(fabsl(q)+fabsl(r));
		      v=fabsl(p)*(fabsl(a[m-1][m-1])+fabsl(z)+fabsl(a[m+1][m+1]));
		      /*  if (a1 * (fabs (q) + fabs (r)) <= GSL_DBL_EPSILON * fabs (p) * (a2 + a3))
			  break; */
 
		      if (u <= EPS*v)
			break;
		      //Equation (W ebnote 16.24).
		    }
		  for (i=m;i<nn-1;i++)
		    {
		      a[i+2][i]=0.0;
		      if (i != m) a[i+2][i-1]=0.0;
		    }
		  for (k=m;k<nn;k++)
		    {
		      //Double QR step on rows l to nn and columns m to nn .
		      if (k != m) 
			{
			  p=a[k][k-1];
			  //Begin setup of Householder vector.
			  q=a[k+1][k-1];
			  r=0.0;
			  if (k+1 != nn) 
			    r=a[k+2][k-1];
			  if ((x=fabsl(p)+fabsl(q)+fabsl(r)) != 0.0)
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGNL(sqrtl(p*p+q*q+r*r),p)) != 0.0)
			{
			  if (k == m) 
			    {
			      if (l != m)
				a[k][k-1] = -a[k][k-1];
			    } 
			  else
			    a[k][k-1] = -s*x;
			  p += s;
			  //Equations (Webnote 16.22).
			  x=p/s;
			  y=q/s;
			  z=r/s;
			  q /= p;
			  r /= p;
			  for (j=k;j<nn+1;j++)
			    {
			      //Row mo di cation.
			      p=a[k][j]+q*a[k+1][j];
			      if (k+1 != nn)
				{
				  p += r*a[k+2][j];
				  a[k+2][j] -= p*z;
				}
			      a[k+1][j] -= p*y;
			      a[k][j] -= p*x;
			    }
			  mmin = nn < k+3 ? nn : k+3;
			  for (i=l;i<mmin+1;i++)
			    {
			      //Column modification.
			      p=x*a[i][k]+y*a[i][k+1];
			      if (k+1 != nn) {
				p += z*a[i][k+2];
				a[i][k+2]
				  -= p*r;
			      }
			      a[i][k+1] -= p*q;
			      a[i][k] -= p;
			    }
			}
		    }
		}
	    }
	} 
      while (l+1 < nn);
    }
}
void QRfactorizationl(long double hess[4][4], complex long double sol[4], int *ok, int n)
{
  /* pagina 615 Num. Rec. */  
  balancel(hess, n);
  hqrl(hess, sol, ok, n);
}
void solve_numrecl(long double *coeff, complex long double *csol, int m, int *ok)
{
  /* Find all the roots of a polynomial with real coefficients, a4*x^4+a3*x^3+a2*x^2+a1*x+a0, 
   * given the coefficients coeff[0..m]. The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector csol[0..m-1], sorted in descending order by their real parts.*/
  /* pagina 497 Num. Rec. */
  long double hess[4][4]; 
  int j, k;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorizationl(hess, csol, ok, m);
}
