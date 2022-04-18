#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include <unistd.h>
#include <float.h>
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void balance(double a[4][4])
{
  const double RADIX=FLT_RADIX;
  int i, j;
  double scale[4]={1.0,1.0,1.0,1.0};
  int done=0;
  double r, c, g, f, s, sqrdx=RADIX*RADIX;
  const int n=4;
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
		c += fabs(a[j][i]);
		r += fabs(a[i][j]);
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

void hqr(double a[4][4], complex double wri[4], int *ok)
{
  int nn,m,l,k,j,its,i,mmin;
  double z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0, anorm=0.0;
  const double EPS=2.2204460492503131E-16;
  const int n=4;
  for (i=0;i<n;i++)
    //Compute matrix no rm for possible use in lo- cating single small sub diagonal element.
    for (j=IMAX(i-1,0);j<n;j++)
      anorm += fabs(a[i][j]);
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
	      s=fabs(a[l-1][l-1])+fabs(a[l][l]);
	      if (s == 0.0)
		s=anorm;
	      if (fabs(a[l][l-1]) + s == s)
	       	{
	  	  a[l][l-1]=0.0;
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
		  z=sqrt(fabs(q));
		  x += t;
		  if (q >= 0.0)
		    {
		      //...a real pair.
		      z=p+SIGN(z,p);
		      wri[nn-1]=wri[nn]=x+z;
		      if (z != 0.0)
			wri[nn]=x-w/z;
		    } 
		  else
		    {
		      //...a complex pair.
		      wri[nn]=CMPLX(x+p,-z);
		      wri[nn-1]=conj(wri[nn]);
		    }
		  nn -= 2;
		} 
	      else
		{
		  //No roots found.  Continue iteration.
		  if (its == 480)
		    {
		      printf("Too many iterations in hqr");
		      *ok=0;
		      return;
		      //exit(-1);
		    }
		  if (its % 10 == 0 && its > 0)
		    {
		      //Form exceptional shift.
		      t += x;
		      for (i=0;i<nn+1;i++)
			a[i][i] -= x;
		      s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
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
		      s=fabs(p)+fabs(q)+fabs(r);
		      //Scale to prevent over flow or under flow.
		      p /= s;
		      q /= s;
		      r /= s;
		      if (m == l) 
			break;
		      u=fabs(a[m][m-1])*(fabs(q)+ fabs(r));
		      v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
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
			  if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0)
			    {
			      p /= x;
			      //Scale to prevent over flow or under flow.
			      q /= x;
			      r /= x;
			    }
			}
		      if ((s=SIGN(sqrt(p*p+q*q+ r*r),p)) != 0.0)
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
			      //Row modification.
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
			      p=x*a[i][k]+y*a[i][k+1 ];
			      if (k+1 != nn) {
				p += z*a[i][k+2];
				a[i][k+2] -= p*r;
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
void QRfactorization( double hess[4][4], complex double sol[4], int *ok)
{
  /* pag. 615 Num. Rec. */  
  balance(hess);
  hqr(hess, sol, ok);
}
void solve_numrec (double coeff[5], complex double csol[4], int *ok)
{
  /* Find all the roots of a polynomial with real coefficients, 
   * coeff[4]*x^4+coeff[3]*x^3+coeff[2]*x^2+coeff[1]*x+coeff[0], 
   * The method is to construct an upper Hessenberg matrix whose 
   * eigenvalues are the desired roots and then use the routine Unsymmeig. The roots are returned 
   * in the complex vector rt[0..m-1], sorted in descending order by their real parts.*/
  /* pag. 497 Num. Rec. */
  const int m=4;
  double hess[4][4];
  int j, k;
  for (k=0;k<m;k++) { //Construct the matrix.
    hess[0][k] = -coeff[m-k-1]/coeff[m];
    for (j=1;j<m;j++) hess[j][k]=0.0;
    if (k != m-1) hess[k+1][k]=1.0;
  }
  QRfactorization(hess, csol, ok);
}
