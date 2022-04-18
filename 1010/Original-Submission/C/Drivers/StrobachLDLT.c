#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include <float.h>
#define Sqr(x) ((x)*(x))
extern void solve_cubic_analytic(double *coeff, complex double sol[3]);
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
void solve_quadratic_cmplx(double coeff[3], complex double *sol)
{
  /* numeric error safe version of solve_quadratic from Numerical Recipe */
  double delta, a, b, c, q;
  complex double cq;
  a = coeff[2];
  b = coeff[1];
  c = coeff[0];
  delta = Sqr(b) - 4.0*a*c;
  if (delta > 0.0)
    {
      q = -0.5*(b+copysign(1.0,b)*sqrt(delta));
      sol[0] = q/a;
      sol[1] = c/q;
    } 
  else if (delta == 0)
    {
      sol[0] = sol[1] = -b/(2.0*a);
    }
  else
    {
      cq = -0.5*(b+copysign(1.0,b)*csqrt(delta));
      sol[0] = cq/a;
      sol[1] = c/cq; 
    }
}

void two_lin_eqs(double fmat[2][2],double evec[2], double dalf[2])
{ 
  double r,x1,x2,c,s;

  r=sqrt(Sqr(fmat[0][0])+Sqr(fmat[1][0]));
  if(r>0.0)
    {
      c=fmat[0][0]/r;
      s=fmat[1][0]/r;
    }
  else
    {
      c=1.0;
      s=0.0;
    }

  x1=c*fmat[0][0]+s*fmat[1][0];
  x2=-s*fmat[0][0]+c*fmat[1][0];
  fmat[0][0]=x1;
  fmat[1][0]=x2;

  x1=c*fmat[0][1]+s*fmat[1][1];
  x2=-s*fmat[0][1]+c*fmat[1][1];
  fmat[0][1]=x1;
  fmat[1][1]=x2;      

  x1=c*evec[0]+s*evec[1];
  x2=-s*evec[0]+c*evec[1];
  evec[0]=x1;
  evec[1]=x2; 

  if(fmat[1][1]==0.0)
    dalf[1]=0.0;
  else
    dalf[1]=evec[1]/fmat[1][1];

  if(fmat[0][0]==0.0)
    dalf[0]=0.0;
  else
    dalf[0]=(evec[0]-fmat[0][1]*dalf[1])/fmat[0][0];
}
void quadratic(double aa,double bb,double cc, double dd, double a, double b, complex double roots[2])
{ 
  double diskr,div,zmax,zmin,evec[2], fmat[2][2],dpar[2],at,bt,err,errt;
  int iter; 
  //-------------------------------- parameter backward correction step:       
  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)
  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23)

  if(err!=0.0)    
    {      
      for (iter=0; iter < 8; iter++) 
	{ 
	  fmat[0][0]=-b*aa+2*a*b;                         // equation (5.19)
	  fmat[0][1]=bb-2*b-a*aa+a*a;                    // equation (5.19)
	  fmat[1][0]=b*b-dd;                             // equation (5.19)
	  fmat[1][1]=cc-2*b*aa+2*b*a;                     // equation (5.19)

	  evec[0]=-evec[0];
	  evec[1]=-evec[1];

	  two_lin_eqs(fmat,evec,dpar);              // equation (5.20)

	  at=a;
	  bt=b;

	  a=a+dpar[0];                                   // equation (5.21)
	  b=b+dpar[1];                                   // equation (5.22)

	  evec[0]=bb*b-b*b-a*b*aa+a*a*b-dd;            // equation (5.18)
	  evec[1]=cc*b-b*b*aa+b*b*a-a*dd;              // equation (5.18)

	  errt=err;
	  err=fabs(evec[0])+fabs(evec[1]);                 // equation (5.23) 

	  if(err==0.0)
	    break;	// terminate

	  if(err>=errt)   // terminate without parameter update
	    {
	      a=at;
	      b=bt;  
	      break;
	    }
	}
    } 
  //---------------------------------------- solve a quadratic equation:     

  diskr=a*a-4*b;   

  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      roots[0]=CMPLX(-a/2,sqrt(-diskr)/2);
      roots[1]=CMPLX(-a/2,-sqrt(-diskr)/2);      
    }   
} 

void  ac_fit(double a, double b, double c, double *AQ, double *BQ, double *CQ, double *DQ)
{
  double cmat[3][2],x1,x2,r,cr,sr; 
  double aq, bq, cq, dq;
  int j;

  aq = *AQ;
  bq = *BQ;
  cq = *CQ;
  dq = *DQ;
  if (fabs(aq) > fabs(cq))
    {                    // equation (5.12)
      cmat[0][0]=aq;
      cmat[0][1]=b-bq-dq;      
      cmat[1][0]=bq;
      cmat[1][1]=c-aq*dq;      
      cmat[2][0]=1.0;
      cmat[2][1]=a-aq;
   
      x1=cmat[0][0];
      x2=cmat[2][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2)>fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}     
     
      for (j=0; j < 2; j++)
	{
	  x1=cr*cmat[0][j]+sr*cmat[2][j];
	  x2=-sr*cmat[0][j]+cr*cmat[2][j];
	  cmat[0][j]=x1;
	  cmat[2][j]=x2;
	}
    
      x1=cmat[0][0];
      x2=cmat[1][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}
      for (j=0; j < 2; j++) 
	{
	  x1=cr*cmat[0][j]+sr*cmat[1][j];
	  x2=-sr*cmat[0][j]+cr*cmat[1][j];
	  cmat[0][j]=x1;
	  cmat[1][j]=x2;
	}
      cq=cmat[0][1]/cmat[0][0];
    } 
  else                                          // equation (5.13)
    {
      cmat[0][0]=cq;
      cmat[0][1]=b-bq-dq;      
      cmat[1][0]=dq;
      cmat[1][1]=c-bq*cq;      
      cmat[2][0]=1.0;
      cmat[2][1]=a-cq;

      x1=cmat[0][0];
      x2=cmat[2][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1;
	    } 
	}     

      for (j=0; j < 2; j++)
	{
	  x1=cr*cmat[0][j]+sr*cmat[2][j];
	  x2=-sr*cmat[0][j]+cr*cmat[2][j];
	  cmat[0][j]=x1;
	  cmat[2][j]=x2;
	}

      x1=cmat[0][0];
      x2=cmat[1][0];
      r=sqrt(x1*x1+x2*x2);
      if(r==0.0)
	{
	  cr=1.0;
	  sr=0.0;
	}
      else
	{
	  if(fabs(x2) > fabs(x1))
	    {
	      sr=x2/r;
	      cr=sr*x1/x2;
	    }
	  else
	    {
	      cr=x1/r;
	      sr=cr*x2/x1; 
	    }
	}   

      for (j=0; j < 2; j++) 
	{
	  x1=cr*cmat[0][j]+sr*cmat[1][j];
	  x2=-sr*cmat[0][j]+cr*cmat[1][j];
	  cmat[0][j]=x1;
	  cmat[1][j]=x2;
	}
      aq=cmat[0][1]/cmat[0][0];
    } 
  *AQ=aq;
  *BQ=bq;
  *CQ=cq;
  *DQ=dq;
  //endif   

  //      return
  //      end
}
void r_quadratic(double a,double b,double c, double *xmax,double *xmin,int *iflag) // real solution: iflag=1
                                                    // no real sol:   iflag=0
{
  double diskr,nenn;

  diskr=b*b-4*a*c;

  if(diskr < 0.0)
    {
      *iflag=0;
      *xmax=-b/(2*a);
      *xmin=-b/(2*a);
    }
  else
    {
      *iflag=1;

      if(b > 0.0)
	nenn=-b-sqrt(diskr);
      else
	nenn=-b+sqrt(diskr);

      if(nenn==0.0 || a==0.0)
	{	
	  *xmax=0.0;
	  *xmin=0.0;
	}
      else
	{
	  *xmax=nenn/(2*a);
	  *xmin=2*c/nenn;
	}
    }
}
//====================================================================

void depressed_cubic_root(double g, double h, double *root)
{
  double x0=0.0,xt0,xr,xmax,xmin,a,b,c,xabs,oldabs;
  int iter, iflag;    
  xabs=0.0;

  if(h==0.0)
    {
      x0=0.0;
    }
  else
    {
      //-------------------------------------------------------------------- 
      if(h<0.0)
	{

	  xr=sqrt(-h);                                  // equation (3.20)

	  if(g > xr)
	    x0=-h/g;                                       // equation (3.26)
	  else if(g < -xr*xr)
	    x0=sqrt(-g);                                  // equation (3.22)
	  else
	    x0=xr;                                         // equation (3.24)

	}
      else if(h > 0.0)
	{
	  xr=sqrt(h);                                   // equation (3.21)

	  if(g > xr)
	    x0=-h/g;                                       // equation (3.26)
	  else if(g < -xr*xr)
	    x0=-sqrt(-g);                                 // equation (3.23)
	  else
	    x0=-xr;                                        // equation (3.25)

	}

      //-------------------------------      
      for (iter=0; iter < 8; iter++)	
	{
	  a=x0*x0;
	  b=-h;
	  c=g*x0*x0+2*h*x0;

	  r_quadratic(a,b,c,&xmax,&xmin,&iflag);      // equation (3.10)

	  x0=xmin;                                       // equation (3.11)
	  xt0=x0;

	  a=2*x0;
	  b=g-x0*x0;
	  c=h;

	  r_quadratic(a,b,c,&xmax,&xmin,&iflag);       // equation (3.15)

	  if(h < 0.0)
	    {     // take the positive solution
	      if(xmax>=0.0)
		x0=xmax;                                       // equation (3.16)
	      else
		x0=xmin;                                       // equation (3.16)

	    }

	  if(h > 0.0)
	    {     
	      // take the negative solution 
	      if(xmax <= 0.0)
		x0=xmax;                                       // equation (3.17)
	      else
		x0=xmin;                                       // equation (3.17)

	    }       

	  oldabs=xabs;
	  xabs=fabs(xt0-x0);
	  if(iter > 0 && xabs == oldabs)   // terminate
	    { 
	      break; //goto 100
	    }
	  if(iter > 0  && xabs == 0.0)      // terminate
	    { 
	      //printf("qui2\n");
	      //fine = 1;
	      break;
	      //goto 100
	    }     

	}//enddo

      *root=x0;

      //--------------------------------------------------------------------      
    }
  //--------------------------------------------------------------------      
} 

void  cubic_B_shift(double a, double b, double c, double d, double *phi0)
{
  double g,h,gg,hh,aq,bq,cq,dq,s,diskr;

  // !------------------------------------------------------- the B-shift:

  diskr=9*a*a-24*b;                    //         ! equation (3.58)
      
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2*b/(3*a+diskr);                            // equation (3.58)
      else
	s=-2*b/(3*a-diskr);                            // equation (3.58)
	  
    }
  else
    {      
      s=-a/4;                                    // equations (3.59-60)
    }
      
  // !--------------------------- the shift transformation (Horner forms):

  aq=a+4*s;                                      // equation (3.45)
  bq=b+3*s*(a+2*s);                              // equation (3.46)
  cq=c+s*(2*b+s*(3*a+4*s));                      // equation (3.47)
  dq=d+s*(c+s*(b+s*(a+s)));                      // equation (3.48)     

  gg=bq*bq/9;
  hh=aq*cq;      
  g=hh-4*dq-3*gg;                                // equation (3.49)
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq;         // equation (3.50)

//------------------ call the parabola/reciprocal intersection method:      
      
  depressed_cubic_root(g,h,phi0);
 
// -------------------------------------------------------------------- 
      
} 

/*
   Q(z) = z^4 + a*z^3 + b*z^2 + c*z + d 
   = (z-z_1)*(z-z_2)*(z-z_3)*(z-z_4)

  coefficient: a  b  c  d
                value: 0  0  0  0  (all coefficients zero)
                 value: 0  0  0  x  (only d=nonzero)
                 value: 0  0  x  0
                 value: 0  x  0  0
                 value: x  0  0  0
*/
void CLDLT_quartic(double coeff[5], complex double roots[4])      
{
  double a,b,c,d,phi0,aq,bq,cq,dq,d2,l1,l2,l3;
  double ssd, gamma,del2, cbq[3], cubc[4];
  complex double sbq[2], acx,bcx,cdiskr,zx1,zx2,zxmax,zxmin, qroots[2];
  complex long double rri, rmri;
  //----------------------------- calculate the antidiagonal shift phi0:

  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  /* special cases to handle */
  if (a==0 && b==0 && c==0 && d==0)
    {
      roots[0]=roots[1]=roots[2]=roots[3]=0;
      return;
    }
  else if (a==0 && b==0 && c==0)
    {
      if (d < 0.0)
	{
	  ssd = pow(-d,0.25);
	  roots[0]=ssd+I*0.0;
	  roots[1]=-ssd+I*0.0;
	  roots[2]=0.0+I*ssd;
	  roots[3]=0.0-I*ssd;
	}
      else
	{
	  ssd = pow(d,0.25);
	  rri = cosl(M_PI/4.0)+sinl(M_PI/4.0)*I;
	  rmri=I*(cosl(M_PI/4.0)+sinl(M_PI/4.0)*I);
	  roots[0]=rri*ssd;
	  roots[1]=-rri*ssd;
	  roots[2]=rmri*ssd;
	  roots[3]=-rmri*ssd;
	}
      return;
    }
  else if (a==0 && b==0 && d==0) 
    {
      cubc[3]=1.0;
      cubc[2]=0.0;
      cubc[1]=0.0;
      cubc[0]=c;
      solve_cubic_analytic(cubc, roots);
      roots[3]=0.0+I*0.0;
      return;
    }
  else if (b==0 && c==0 && d==0)
    {
      roots[0]=roots[1]=roots[2]=0.0+I*0.0;
      roots[3]=-a+I*0.0;
      return;
    }
  cubic_B_shift(a,b,c,d,&phi0);     

  //------------------------------------------ compute LDLT parameters:      

  l1=a/2;                                         // equation (4.2)
  l3=b/6+phi0/2;                                  // equation (4.3)
  del2=c-a*l3;                                    // equation (4.10) 

  if(d<=0.0)
    {
      if(del2==0.0)
	{
	  l2=0.0;
	}
      else
	{
	  l2=2*(d-l3*l3)/del2;                           // equation (4.12)
	}
    }  
  else
    {
      if(b-l1*l1-2*l3==0.0)
	{
	  l2=0.0;
	}
      else
	{
	  l2=del2/(2*(b-l1*l1-2*l3));              // equation (4.11)
	}
    }
  if(l2==0.0)
    d2=2*b/3-phi0-l1*l1;
  else
    d2=del2/(2*l2);                                // equation (4.15) 
  if(a==0.0 && c==0.0)  // handle a bi-quadratic equation
    {
      /* solve directly */
      cbq[2]=1.0;
      cbq[1]=b;	
      cbq[0]=d;
      solve_quadratic_cmplx(cbq,sbq);
      roots[0] = csqrt(sbq[0]);
      roots[1] = -csqrt(sbq[0]);
      roots[2] = csqrt(sbq[1]);
      roots[3] = -csqrt(sbq[1]);
      return;
    }

  //-------- decide whether a real domain decomposition (equation (1.3))
  //- or a complex domain decomposition (equation (1.4)) should be used:
  //--------------------------------------------------------------------      

  if(d2<=0.0) 
    {
      // assume a real quadratic decomposition (1.3)
      gamma=sqrt(-d2);                               // equation (2.8) 
      //printf("gamma=%.16G\n", gamma);

      aq=l1+gamma;                                    // equation (5.3)
      bq=l3+gamma*l2;                                 //equation (5.5)

      cq=l1-gamma;                                    // equation (5.4)
      dq=l3-gamma*l2;                                 // equation (5.6)

      if(fabs(dq) < fabs(bq))
	dq=d/bq;                                        // equation (5.9)
      else if(fabs(dq) > fabs(bq))
	bq=d/dq;                                       // equation (5.10)

      //------------------------ aq and cq coefficients backward correction:

      ac_fit(a,b,c,&aq,&bq,&cq,&dq);
      //------------------------------------ roots of the quadratic factors:   

      quadratic(a,b,c,d,aq,bq,qroots);
      roots[0]=qroots[0];
      roots[1]=qroots[1];        
      quadratic(a,b,c,d,cq,dq,qroots);
      roots[2]=qroots[0];
      roots[3]=qroots[1];               

      //--------------------------------------------------------------------      
    }
  else    // (d2.gt.0.0) assume a complex quadratic decomposition: 
    {
      // --------------------------------------------------------------------       

      gamma=sqrt(d2);                                // equation (2.8)

      acx=CMPLX(l1,gamma);                          // equation (2.28)
      bcx=CMPLX(l3,gamma*l2);                       // equation (2.29)

      cdiskr=acx*acx/4-bcx;               

      zx1=-acx/2+csqrt(cdiskr);                       // equation (2.31)
      zx2=-acx/2-csqrt(cdiskr);                       // equation (2.32)

      if(cabs(zx1) > cabs(zx2))                  // equation (2.33)
	zxmax=zx1;
      else
	zxmax=zx2;
	     

      zxmin=bcx/zxmax;                               // equation (2.34)

      roots[0]=zxmin;
      roots[1]=conj(zxmin);
      roots[2]=zxmax;
      roots[3]=conj(zxmax);

      //--------------------------------------------------------------------

    }       
  //--------------------------------------------------------------------      
}
