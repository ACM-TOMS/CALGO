#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multiroots/dogleg.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include "enorm.c"

static void compute_diag (const gsl_matrix * J, gsl_vector * diag);
static void update_diag (const gsl_matrix * J, gsl_vector * diag);
static MpIeee compute_delta(gsl_vector * diag, gsl_vector * x);
static void compute_df (const gsl_vector * f_trial, const gsl_vector * f, gsl_vector * df);
static void compute_wv (const gsl_vector * qtdf, const gsl_vector *rdx, const gsl_vector *dx, const gsl_vector *diag, MpIeee pnorm, gsl_vector * w, gsl_vector * v);

static MpIeee scaled_enorm(const gsl_vector * d, const gsl_vector * f);

static MpIeee scaled_enorm(const gsl_vector * d, const gsl_vector * f) {
  MpIeee e2=  MpIeee( "0" ) ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    MpIeee fi=  gsl_vector_get(f, i);
    MpIeee di=  gsl_vector_get(d, i);
    MpIeee u=  di * fi;
    e2 += u * u ;
  }
  return sqrt(e2);
}

static MpIeee enorm_sum(const gsl_vector * a, const gsl_vector * b);

static MpIeee enorm_sum(const gsl_vector * a, const gsl_vector * b) {
  MpIeee e2=  MpIeee( "0" ) ;
  size_t i, n = a->size ;
  for (i = 0; i < n ; i++) {
    MpIeee ai=  gsl_vector_get(a, i);
    MpIeee bi=  gsl_vector_get(b, i);
    MpIeee u=  ai + bi;
    e2 += u * u ;
  }
  return sqrt(e2);
}

static void
compute_wv (const gsl_vector * qtdf, const gsl_vector *rdx, const gsl_vector *dx, const gsl_vector *diag, MpIeee pnorm, gsl_vector * w, gsl_vector * v)
{
  size_t i, n = qtdf->size;

  for (i = 0; i < n; i++)
    {
      MpIeee qtdfi=  gsl_vector_get (qtdf, i);
      MpIeee rdxi=  gsl_vector_get (rdx, i);
      MpIeee dxi=  gsl_vector_get (dx, i);
      MpIeee diagi=  gsl_vector_get (diag, i);

      gsl_vector_set (w, i, (qtdfi - rdxi) / pnorm);
      gsl_vector_set (v, i, diagi * diagi * dxi / pnorm);
    }
}


static void
compute_df (const gsl_vector * f_trial, const gsl_vector * f, gsl_vector * df)
{
  size_t i, n = f->size;

  for (i = 0; i < n; i++)
    {
      MpIeee dfi=  gsl_vector_get (f_trial, i) - gsl_vector_get (f, i);
      gsl_vector_set (df, i, dfi);
    }
}

static void
compute_diag (const gsl_matrix * J, gsl_vector * diag)
{
  size_t i, j, n = diag->size;

  for (j = 0; j < n; j++)
    {
      MpIeee sum=  MpIeee( "0" );
      for (i = 0; i < n; i++)
        {
          MpIeee Jij=  gsl_matrix_get (J, i, j);
          sum += Jij * Jij;
        }
      if (sum == MpIeee( "0" ))
        sum = MpIeee( "1.0" );

      gsl_vector_set (diag, j, sqrt (sum));
    }
}

static void
update_diag (const gsl_matrix * J, gsl_vector * diag)
{
  size_t i, j, n = diag->size;

  for (j = 0; j < n; j++)
    {
      MpIeee cnorm;MpIeee  diagj;MpIeee  sum=  MpIeee( "0" );
      for (i = 0; i < n; i++)
        {
          MpIeee Jij=  gsl_matrix_get (J, i, j);
          sum += Jij * Jij;
        }
      if (sum == MpIeee( "0" ))
        sum = MpIeee( "1.0" );

      cnorm = sqrt (sum);
      diagj = gsl_vector_get (diag, j);

      if (cnorm > diagj)
        gsl_vector_set (diag, j, cnorm);
    }
}

static MpIeee compute_delta(gsl_vector * diag, gsl_vector * x)
{
  MpIeee Dx=  scaled_enorm (diag, x);
  MpIeee factor=  MpIeee( "100" );

  return (Dx > MpIeee( "0" )) ? factor * Dx : factor;
}

static MpIeee compute_actual_reduction(MpIeee fnorm, MpIeee fnorm1)
{
  MpIeee actred;

  if (fnorm1 < fnorm)
    {
      MpIeee u=  fnorm1 / fnorm;
      actred = MpIeee( "1" ) - u * u;
    }
  else
    {
      actred = -MpIeee( "1" );
    }

  return actred;
}

static MpIeee compute_predicted_reduction(MpIeee fnorm, MpIeee fnorm1)
{
  MpIeee prered;

  if (fnorm1 < fnorm)
    {
      MpIeee u=  fnorm1 / fnorm;
      prered = MpIeee( "1" ) - u * u;
    }
  else
    {
      prered = MpIeee( "0" );
    }

  return prered;
}

static void 
compute_qtf (const gsl_matrix * q, const gsl_vector * f, gsl_vector * qtf)
{
  size_t i, j, N = f->size ;

  for (j = 0; j < N; j++)
    {
      MpIeee sum=  MpIeee( "0" );
      for (i = 0; i < N; i++)
        sum += gsl_matrix_get (q, i, j) * gsl_vector_get (f, i);

      gsl_vector_set (qtf, j, sum);
    }
}

static void 
compute_rdx (const gsl_matrix * r, const gsl_vector * dx, gsl_vector * rdx)
{
  size_t i, j, N = dx->size ;

  for (i = 0; i < N; i++)
    {
      MpIeee sum=  MpIeee( "0" );

      for (j = i; j < N; j++)
        {
          sum += gsl_matrix_get (r, i, j) * gsl_vector_get (dx, j);
        }

      gsl_vector_set (rdx, i, sum);
    }
}


static void
compute_trial_step (gsl_vector *x, gsl_vector * dx, gsl_vector * x_trial)
{
  size_t i, N = x->size;

  for (i = 0; i < N; i++)
    {
      MpIeee pi=  gsl_vector_get (dx, i);
      MpIeee xi=  gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + pi);
    }
}

static int
 newton_direction(const gsl_matrix * r, const gsl_vector * qtf, gsl_vector * p)
{
  const size_t N = r->size2;
  size_t i;
  int  status;

  status = gsl_linalg_R_solve (r, qtf, p);

#ifdef DEBUG
  {cout<<"rsolve status = "<< status<<"\n";}
#endif

  for (i = 0; i < N; i++)
    {
      MpIeee pi=  gsl_vector_get (p, i);
      gsl_vector_set (p, i, -pi);
    }

  return status;
}

static void
gradient_direction (const gsl_matrix * r, const gsl_vector * qtf,
                    const gsl_vector * diag, gsl_vector * g)
{
  const size_t M = r->size1;
  const size_t N = r->size2;

  size_t i, j;

  for (j = 0; j < M; j++)
    {
      MpIeee sum=  MpIeee( "0" );
      MpIeee dj;

      for (i = 0; i < N; i++)
        {
          sum += gsl_matrix_get (r, i, j) * gsl_vector_get (qtf, i);
        }

      dj = gsl_vector_get (diag, j);
      gsl_vector_set (g, j, -sum / dj);
    }
}

static void
minimum_step (MpIeee gnorm, const gsl_vector * diag, gsl_vector * g)
{
  const size_t N = g->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee gi=  gsl_vector_get (g, i);
      MpIeee di=  gsl_vector_get (diag, i);
      gsl_vector_set (g, i, (gi / gnorm) / di);
    }
}

static void
compute_Rg (const gsl_matrix * r, const gsl_vector * gradient, gsl_vector * Rg)
{
  const size_t N = r->size2;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      MpIeee sum=  MpIeee( "0" );

      for (j = i; j < N; j++)
        {
          MpIeee gj=  gsl_vector_get (gradient, j);
          MpIeee rij=  gsl_matrix_get (r, i, j);
          sum += rij * gj;
        }

      gsl_vector_set (Rg, i, sum);
    }
}

static void
scaled_addition (MpIeee alpha, gsl_vector * newton, MpIeee beta, gsl_vector * gradient, gsl_vector * p)
{
  const size_t N = p->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      MpIeee ni=  gsl_vector_get (newton, i);
      MpIeee gi=  gsl_vector_get (gradient, i);
      gsl_vector_set (p, i, alpha * ni + beta * gi);
    }
}

static int
 dogleg(const gsl_matrix * r, const gsl_vector * qtf,
        const gsl_vector * diag, MpIeee delta,
        gsl_vector * newton, gsl_vector * gradient, gsl_vector * p)
{
  MpIeee qnorm;MpIeee  gnorm;MpIeee  sgnorm;MpIeee  bnorm;MpIeee  temp;

  newton_direction (r, qtf, newton);

#ifdef DEBUG
  {cout<<"newton = ";} gsl_vector_fprintf(stdout, newton, "%g"); {cout<<"\n";}
#endif

  qnorm = scaled_enorm (diag, newton);

  if (qnorm <= delta)
    {
      gsl_vector_memcpy (p, newton);
#ifdef DEBUG
      {cout<<"took newton (qnorm = "<<setiosflags((ios::floatfield))<< qnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"  <=   delta = "<<setiosflags((ios::floatfield))<< delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<")\n";}
#endif
      return GSL_SUCCESS;
    }

  gradient_direction (r, qtf, diag, gradient);

#ifdef DEBUG
  {cout<<"grad = ";} gsl_vector_fprintf(stdout, gradient, "%g"); {cout<<"\n";}
#endif

  gnorm = enorm (gradient);

  if (gnorm == MpIeee( "0" ))
    {
      MpIeee alpha=  delta / qnorm;
      MpIeee beta=  MpIeee( "0" );
      scaled_addition (alpha, newton, beta, gradient, p);
#ifdef DEBUG
      {cout<<"took scaled newton because gnorm = 0\n";}
#endif
      return GSL_SUCCESS;
    }

  minimum_step (gnorm, diag, gradient);

  compute_Rg (r, gradient, p);  /* Use p as temporary space to compute Rg */

#ifdef DEBUG
  {cout<<"mingrad = ";} gsl_vector_fprintf(stdout, gradient, "%g"); {cout<<"\n";}
  {cout<<"Rg = ";} gsl_vector_fprintf(stdout, p, "%g"); {cout<<"\n";}
#endif

  temp = enorm (p);
  sgnorm = (gnorm / temp) / temp;

  if (sgnorm > delta)
    {
      MpIeee alpha=  MpIeee( "0" );
      MpIeee beta=  delta;
      scaled_addition (alpha, newton, beta, gradient, p);
#ifdef DEBUG
      {cout<<"took gradient\n";}
#endif
      return GSL_SUCCESS;
    }

  bnorm = enorm (qtf);

  {
    MpIeee bg=  bnorm / gnorm;
    MpIeee bq=  bnorm / qnorm;
    MpIeee dq=  delta / qnorm;
    MpIeee dq2=  dq * dq;
    MpIeee sd=  sgnorm / delta;
    MpIeee sd2=  sd * sd;

    MpIeee t1=  bg * bq * sd;
    MpIeee u=  t1 - dq;
    MpIeee t2=  t1 - dq * sd2 + sqrt (u * u + (MpIeee( "1" )-dq2) * (MpIeee( "1" ) - sd2));

    MpIeee alpha=  dq * (MpIeee( "1" ) - sd2) / t2;
    MpIeee beta=  (MpIeee( "1" ) - alpha) * sgnorm;

#ifdef DEBUG
    {cout<<"bnorm = "<<setiosflags((ios::floatfield))<< bnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    {cout<<"gnorm = "<<setiosflags((ios::floatfield))<< gnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    {cout<<"qnorm = "<<setiosflags((ios::floatfield))<< qnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    {cout<<"delta = "<<setiosflags((ios::floatfield))<< delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    {cout<<"alpha = "<<setiosflags((ios::floatfield))<< alpha;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"   beta = "<<setiosflags((ios::floatfield))<< beta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    {cout<<"took scaled combination of newton and gradient\n";}
#endif

    scaled_addition (alpha, newton, beta, gradient, p);
  }

  return GSL_SUCCESS;
}
