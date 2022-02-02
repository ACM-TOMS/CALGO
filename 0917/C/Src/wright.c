#include "wright.h"

/**********************************************************************/
/* wrightomega is the simple routine for evaluating the wright omega  */
/* function.                                                          */
/*                                                                    */
/* Calling:                                                           */
/*    w = wrightomega(z)                                              */
/*                                                                    */
/* Input:                                                             */
/*   z  --  double complex                                            */
/*          Value to evaluate Wrightomega(z) at.                      */
/*                                                                    */
/* Output:                                                            */
/*   w  --  double complex                                            */
/*          Value of Wrightomega(z)                                   */
/*                                                                    */
/**********************************************************************/


/**********************************************************************/
/* wrightomega_ext is the extended routine for evaluating the wright  */
/* omega function with the option of extracting the last update step, */
/* the penultimate residual and the condition number estimate.        */
/*                                                                    */
/* Calling:                                                           */
/*   success = wrightomega_ext(z,w,e,r,cond);                         */
/*                                                                    */
/* Input:                                                             */
/*   z  --  double complex                                            */
/*          Value to evaluate Wrightomega(z) at.                      */
/*                                                                    */
/*   w  --  double complex*                                           */
/*          Pointer to return value of Wrightomega(z).                */
/*                                                                    */
/*   e  --  double complex*                                           */
/*          Pointer to the last update step in the iterative scheme.  */
/*                                                                    */
/*   r  --  double complex*                                           */
/*          Pointer to penultimate residual r_k = z - w_k - log(w_k)  */
/*                                                                    */
/*   cond  --  double complex*                                        */
/*         Pointer to the condition number estimate. If NULL the      */
/*         condition number is not calculated.                        */
/*                                                                    */
/* Output: returns 0 on sucessful exit.                               */
/**********************************************************************/

int
wrightomega_ext(double complex z,double complex *w,\
                double complex *e,double complex *r,    \
                double complex *cond)
{
  double pi=M_PI,s=1.0;
  double x,y,ympi,yppi,near;
  double complex pz,wp1,t;


  /* extract real and imaginary parts of z */
  x=creal(z);
  y=cimag(z);
  
  /* compute if we are near the branch cuts */
  ympi=y-pi;
  yppi=y+pi;
  near=0.1e-1;

  /* Test for floating point exceptions */
  /*****************************/
  /* NaN output for NaN input  */
  /*****************************/
  if(isnan(x) || isnan(y))
    {
      *w = 0.0/0.0 + 0.0/0.0*I;
      *e = 0.0 + 0.0*I;
      *r = 0.0 + 0.0*I;
      return 0;
    }
  /*********************************/
  /* Signed zeros between branches */
  /*********************************/
  else if(isinf(x) && (x < 0.0) && (-pi < y) && (y<= pi))
    {
      if( fabs(y) <= pi/2.0)
        {
          *w = +0.0;
        }
      else
        {
          *w = -0.0;
        }
      
      if(y>=0)
        {
          *w += 0.0*I;
        }
      else
        {
          *w += -1.0*0.0*I;
        }

      *e = 0.0 + 0.0*I;
      *r = 0.0 + 0.0*I;
      return 0;
    }
  /**************************/
  /* Asymptotic for large z */
  /**************************/
  else if(isinf(x) || isinf(y))
    {
      *w = x + I*y;
      *e = 0.0 + 0.0*I;
      *r = 0.0 + 0.0*I;
      return 0;
    }

  /******************************************/
  /* Test If exactly on the singular points */
  /******************************************/
  if((x==-1.0) && (fabs(y)==pi))
    {
      *w = -1.0 + 0.0*I;
      *e = 0.0 + 0.0*I;
      *r = 0.0 + 0.0*I;
      return 0;
    }


  /* Choose approximation based on region */
  /**********************************/
  /* Region 1: upper branch point   */
  /* Series about z=-1+Pi*I         */
  /**********************************/
  if ((-2.0<x && x<=1.0 && 1.0<y && y< 2.0*pi))
    {
      pz=conj(csqrt(conj(2*(z+1-I*pi))));
      *w=-1.0+(I+(1.0/3.0+(-1.0/36.0*I+(1.0/270.0+1.0/4320.0*I*pz)*pz)*pz)*pz)*pz;
    }
  /**********************************/
  /* Region 2: lower branch point   */
  /* Series about z=-1-Pi*I         */
  /**********************************/
  else if ((-2.0<x && x<=1.0 && -2.0*pi<y && y<-1.0))
    {
      pz=conj(csqrt(conj(2*(z+1.0+I*pi))));
      *w=-1.0+(-I+(1.0/3.0+(1.0/36.0*I+(1.0/270.0-1.0/4320.0*I*pz)*pz)*pz)*pz)*pz;
    }
  /*********************************/
  /* Region 3: between branch cuts */
  /* Series: About -infinity       */
  /*********************************/
  else if (x <= -2.0 && -pi < y && y <= pi)
    {
      pz=cexp(z);
      *w=(1.0+(-1.0+(3.0/2.0+(-8.0/3.0+125.0/24.0*pz)*pz)*pz)*pz)*pz;
    }
  /**********************/
  /* Region 4: Mushroom */
  /* Series about z=1   */
  /**********************/
  else if (((-2.0 < x) && (x<=1.0) && (-1.0 <= y) && (y <= 1.0))
           || ((-2.0 < x) && (x-0.1e1)*(x-0.1e1)+y*y<=pi*pi))
    {
      pz=z-1.0;
      *w=1.0/2.0+1.0/2.0*z+(1.0/16.0+(-1.0/192.0+(-1.0/3072.0+13.0/61440.0*pz)*pz)*pz)*pz*pz;
    }
  /*************************/
  /* Region 5: Top wing    */
  /* Negative log series   */
  /*************************/
  else if (x<=-0.105e1 && pi<y && y-pi<=-0.75e0*(x+0.1e1))
    {
      t=z-I*pi;
      pz=clog(-t);
      *w=((1.0+(-3.0/2.0+1.0/3.0*pz)*pz)*pz+((-1.0+1.0/2.0*pz)*pz+(pz+(-pz+t)*t)*t)*t)/(t*t*t);
    }
  /***************************/
  /* Region 6: Bottom wing   */
  /* Negative log esries     */
  /***************************/
  else if (x<=-0.105e1 && 0.75e0*(x+0.1e1)< y+pi && y+pi<=0.0e0)
    {
      t=z+I*pi;
      pz=clog(-t);
      *w=((1.0+(-3.0/2.0+1.0/3.0*pz)*pz)*pz+((-1.0+1.0/2.0*pz)*pz+(pz+(-pz+t)*t)*t)*t)/(t*t*t);
    }
  /************************************/
  /* Region 7: Everywhere else        */
  /* Series solution about infinity   */
  /************************************/
  else
    {
      pz=clog(z);
      *w=((1.0+(-3.0/2.0+1.0/3.0*pz)*pz)*pz+((-1.0+1.0/2.0*pz)*pz+(pz+(-pz+z)*z)*z)*z)/(z*z*z);
    }

  /**********************************/
  /* Regularize if near branch cuts */
  /**********************************/
  if (x <= -0.1e1 + near && (fabs(ympi) <= near || fabs(yppi) <= near)) 
    { 
      s = -1.0;
      if (fabs(ympi) <= near)
        {
          /* Recompute ympi with directed rounding */
          fesetround(FE_UPWARD);
          ympi = y-pi;
          
          if( ympi <= 0.0)
            {
              fesetround(FE_DOWNWARD);
              ympi = y-pi;
            }
          
          z = x + I*ympi;

          /* Return rounding to default */
          fesetround(FE_TONEAREST);
        }
      else
        {
          /* Recompute yppi with directed rounding */
          fesetround(FE_UPWARD);
          yppi = y + pi;
          
          if( yppi <= 0.0)
            {
              fesetround(FE_DOWNWARD);
              yppi = y + pi;
            }
          
          z = x + I*yppi;
          /* Return rounding to default */
          fesetround(FE_TONEAREST);
        }
    }
  
  /*****************/
  /* Iteration one */
  /*****************/
  *w=s**w;
  *r=z-s**w-clog(*w);
  wp1=s**w+1.0;
  *e=*r/wp1*(2.0*wp1*(wp1+2.0/3.0**r)-*r)/(2.0*wp1*(wp1+2.0/3.0**r)-2.0**r);
  *w=*w*(1.0+*e);
  
  /*****************/
  /* Iteration two */
  /*****************/
  if(cabs((2.0**w**w-8.0**w-1.0)*pow(cabs(*r),4.0)) >= TWOITERTOL*72.0*pow(cabs(wp1),6.0) )
    {
      *r=z-s**w-clog(*w);
      wp1=s**w+1.0;
      *e=*r/wp1*(2.0*wp1*(wp1+2.0/3.0**r)-*r)/(2.0*wp1*(wp1+2.0/3.0**r)-2.0**r);
      *w=*w*(1.0+*e);
    }

  /***********************/
  /* Undo regularization */
  /***********************/
  *w=s**w;

  /***************************************************/
  /*  Provide condition number estimate if requested */
  /***************************************************/
  if(cond != NULL)
    {
      *cond = z/(1.0+*w);
    }
      
  return 0;
}

double complex
wrightomega(double complex z)
{
  double complex w,e,r;
  wrightomega_ext(z,&w,&e,&r,NULL);
  return w;
}
