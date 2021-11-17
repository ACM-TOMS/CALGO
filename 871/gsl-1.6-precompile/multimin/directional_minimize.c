#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multimin/directional_minimize.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static void
take_step (const gsl_vector * x, const gsl_vector * p,
           MpIeee step, MpIeee lambda, gsl_vector * x1, gsl_vector * dx)
{
  gsl_vector_set_zero (dx);
  gsl_blas_daxpy (-step * lambda, p, dx);

  gsl_vector_memcpy (x1, x);
  gsl_blas_daxpy (1.0, dx, x1);
}

static void 
intermediate_point (gsl_multimin_function_fdf * fdf,
                    const gsl_vector * x, const gsl_vector * p,
                    MpIeee lambda, 
                    MpIeee pg,
                    MpIeee stepa, MpIeee stepc,
                    MpIeee fa, MpIeee fc,
                    gsl_vector * x1, gsl_vector * dx, gsl_vector * gradient,
                    MpIeee * step, MpIeee * f)
{
  MpIeee stepb;MpIeee  fb;

trial:
  {
    MpIeee u=  fabs (pg * lambda * stepc);
    stepb = MpIeee( "0.5" ) * stepc * u / ((fc - fa) + u);
  }

  take_step (x, p, stepb, lambda, x1, dx);

  fb = GSL_MULTIMIN_FN_EVAL_F (fdf, x1);

#ifdef DEBUG
  {cout<<"trying stepb = "<<setiosflags((ios::floatfield))<< stepb;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"  fb = "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(18)<< fb;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
#endif

  if (fb >= fa  && stepb > MpIeee( "0.0" ))
    {
      /* downhill step failed, reduce step-size and try again */
      fc = fb;
      stepc = stepb;
      goto trial;
    }
#ifdef DEBUG
  {cout<<"ok!\n";}
#endif

  *step = stepb;
  *f = fb;
  GSL_MULTIMIN_FN_EVAL_DF(fdf, x1, gradient);
}

static void
minimize (gsl_multimin_function_fdf * fdf,
          const gsl_vector * x, const gsl_vector * p,
          MpIeee lambda,
          MpIeee stepa, MpIeee stepb, MpIeee stepc,
          MpIeee fa, MpIeee fb, MpIeee fc, MpIeee tol,
          gsl_vector * x1, gsl_vector * dx1, 
          gsl_vector * x2, gsl_vector * dx2, gsl_vector * gradient,          
          MpIeee * step, MpIeee * f, MpIeee * gnorm)
{
  /* Starting at (x0, f0) move along the direction p to find a minimum
     f(x0 - lambda * p), returning the new point x1 = x0-lambda*p,
     f1=f(x1) and g1 = grad(f) at x1.  */

  MpIeee u=  stepb;
  MpIeee v=  stepa;
  MpIeee w=  stepc;
  MpIeee fu=  fb;
  MpIeee fv=  fa;
  MpIeee fw=  fc;

  MpIeee old2=  fabs(w - v);
  MpIeee old1=  fabs(v - u);

  MpIeee stepm;MpIeee  fm;MpIeee  pg;MpIeee  gnorm1;

  MpIeee iter=  MpIeee( "0" );

  gsl_vector_memcpy (x2, x1);
  gsl_vector_memcpy (dx2, dx1);

  *f = fb;
  *step = stepb;
  *gnorm = gsl_blas_dnrm2 (gradient);

mid_trial:

  iter++;

  if (iter > MpIeee( "10" ))
    {
      return;  /* MAX ITERATIONS */
    }

  {
    MpIeee dw=  w - u;
    MpIeee dv=  v - u;
    MpIeee du=  MpIeee( "0.0" );

    MpIeee e1=  ((fv - fu) * dw * dw + (fu - fw) * dv * dv);
    MpIeee e2=  MpIeee( "2.0" ) * ((fv - fu) * dw + (fu - fw) * dv);

    if (e2 != MpIeee( "0.0" ))
      {
        du = e1 / e2;
      }

    if (du > MpIeee( "0" ) && du < (stepc - stepb) && fabs(du) < MpIeee( "0.5" ) * old2)
      {
        stepm = u + du;
      }
    else if (du < MpIeee( "0" ) && du > (stepa - stepb) && fabs(du) < MpIeee( "0.5" ) * old2)
      {
        stepm = u + du;
      }
    else if ((stepc - stepb) > (stepb - stepa))
      {
        stepm = MpIeee( "0.38" ) * (stepc - stepb) + stepb;
      }
    else
      {
        stepm = stepb - MpIeee( "0.38" ) * (stepb - stepa);
      }
  }

  take_step (x, p, stepm, lambda, x1, dx1);

  fm = GSL_MULTIMIN_FN_EVAL_F (fdf, x1);

#ifdef DEBUG
  {cout<<"trying stepm = "<<setiosflags((ios::floatfield))<< stepm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"  fm = "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(18)<< fm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
#endif

  if (fm > fb)
    {
      if (fm < fv)
        {
          w = v;
          v = stepm;
          fw = fv;
          fv = fm;
        }
      else if (fm < fw)
        {
          w = stepm;
          fw = fm;
        }

      if (stepm < stepb)
        {
          stepa = stepm;
          fa = fm;
        }
      else
        {
          stepc = stepm;
          fc = fm;
        }
      goto mid_trial;
    }
  else if (fm <= fb)
    {
      old2 = old1;
      old1 = fabs(u - stepm);
      w = v;
      v = u;
      u = stepm;
      fw = fv;
      fv = fu;
      fu = fm;

      gsl_vector_memcpy (x2, x1);
      gsl_vector_memcpy (dx2, dx1);

      GSL_MULTIMIN_FN_EVAL_DF (fdf, x1, gradient);
      gsl_blas_ddot (p, gradient, &pg);
      gnorm1 = gsl_blas_dnrm2 (gradient);

#ifdef DEBUG
      {cout<<"p: ";} gsl_vector_fprintf(stdout, p, "%g");
      {cout<<"g: ";} gsl_vector_fprintf(stdout, gradient, "%g");
      {cout<<"gnorm: "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(18)<< gnorm1;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
      {cout<<"pg: "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(18)<< pg;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
      {cout<<"orth: "<<setiosflags((ios::floatfield))<< fabs (pg * lambda/ gnorm1);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif
      *f = fm;
      *step = stepm;
      *gnorm = gnorm1;

      if (fabs (pg * lambda / gnorm1) < tol)
        {
#ifdef DEBUG
          {cout<<"ok!\n";}
#endif
          return; /* SUCCESS */
        }

      if (stepm < stepb)
        {
          stepc = stepb;
          fc = fb;
          stepb = stepm;
          fb = fm;
        }
      else
        {
          stepa = stepb;
          fa = fb;
          stepb = stepm;
          fb = fm;
        }
      goto mid_trial;
    }
}
