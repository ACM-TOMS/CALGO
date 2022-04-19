#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multifit/lmpar.c
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

#include <gsl/gsl_permute_vector_double.h>

#include "qrsolv.c"

static size_t
count_nsing (const gsl_matrix * r)
{
  /* Count the number of nonsingular entries. Returns the index of the
     first entry which is singular. */

  size_t n = r->size2;
  size_t i;

  for (i = 0; i < n; i++)
    {
      MpIeee rii=  gsl_matrix_get (r, i, i);

      if (rii == MpIeee( "0" ))
        {
          break;
        }
    }

  return i;
}


static void
compute_newton_direction (const gsl_matrix * r, const gsl_permutation * perm,
                          const gsl_vector * qtf, gsl_vector * x)
{

  /* Compute and store in x the Gauss-Newton direction. If the
     Jacobian is rank-deficient then obtain a least squares
     solution. */

  const size_t n = r->size2;
  size_t i, j, nsing;

  for (i = 0 ; i < n ; i++)
    {
      MpIeee qtfi=  gsl_vector_get (qtf, i);
      gsl_vector_set (x, i,  qtfi);
    }

  nsing = count_nsing (r);

#ifdef DEBUG
  {cout<<"nsing = "<< nsing<<"\n";}
  {cout<<"r = ";} gsl_matrix_fprintf(stdout, r, "%g"); {cout<<"\n";}
  {cout<<"qtf = ";} gsl_vector_fprintf(stdout, x, "%g"); {cout<<"\n";}
#endif

  for (i = nsing; i < n; i++)
    {
      gsl_vector_set (x, i, 0.0);
    }

  if (nsing > 0)
    {
      for (j = nsing; j > 0 && j--;)
        {
          MpIeee rjj=  gsl_matrix_get (r, j, j);
          MpIeee temp=  gsl_vector_get (x, j) / rjj;
          
          gsl_vector_set (x, j, temp);
          
          for (i = 0; i < j; i++)
            {
              MpIeee rij=  gsl_matrix_get (r, i, j);
              MpIeee xi=  gsl_vector_get (x, i);
              gsl_vector_set (x, i, xi - rij * temp);
            }
        }
    }

  gsl_permute_vector_inverse (perm, x);
}

static void
compute_newton_correction (const gsl_matrix * r, const gsl_vector * sdiag,
                           const gsl_permutation * p, gsl_vector * x,
                           MpIeee dxnorm,
                           const gsl_vector * diag, gsl_vector * w)
{
  size_t n = r->size2;
  size_t i, j;

  for (i = 0; i < n; i++)
    {
      size_t pi = gsl_permutation_get (p, i);

      MpIeee dpi=  gsl_vector_get (diag, pi);
      MpIeee xpi=  gsl_vector_get (x, pi);

      gsl_vector_set (w, i, dpi * (dpi * xpi) / dxnorm);
    }

  for (j = 0; j < n; j++)
    {
      MpIeee sj=  gsl_vector_get (sdiag, j);
      MpIeee wj=  gsl_vector_get (w, j);

      MpIeee tj=  wj / sj;

      gsl_vector_set (w, j, tj);

      for (i = j + 1; i < n; i++)
        {
          MpIeee rij=  gsl_matrix_get (r, i, j);
          MpIeee wi=  gsl_vector_get (w, i);

          gsl_vector_set (w, i, wi - rij * tj);
        }
    }
}

static void
compute_newton_bound (const gsl_matrix * r, const gsl_vector * x, 
                      MpIeee dxnorm,  const gsl_permutation * perm, 
                      const gsl_vector * diag, gsl_vector * w)
{
  /* If the jacobian is not rank-deficient then the Newton step
     provides a lower bound for the zero of the function. Otherwise
     set this bound to zero. */

  size_t n = r->size2;

  size_t i, j;

  size_t nsing = count_nsing (r);

  if (nsing < n)
    {
      gsl_vector_set_zero (w);
      return;
    }

  for (i = 0; i < n; i++)
    {
      size_t pi = gsl_permutation_get (perm, i);

      MpIeee dpi=  gsl_vector_get (diag, pi);
      MpIeee xpi=  gsl_vector_get (x, pi);

      gsl_vector_set (w, i, dpi * (dpi * xpi / dxnorm));
    }

  for (j = 0; j < n; j++)
    {
      MpIeee sum=  MpIeee( "0" );

      for (i = 0; i < j; i++)
        {
          sum += gsl_matrix_get (r, i, j) * gsl_vector_get (w, i);
        }

      {
        MpIeee rjj=  gsl_matrix_get (r, j, j);
        MpIeee wj=  gsl_vector_get (w, j);

        gsl_vector_set (w, j, (wj - sum) / rjj);
      }
    }
}

static void
compute_gradient_direction (const gsl_matrix * r, const gsl_permutation * p,
                            const gsl_vector * qtf, const gsl_vector * diag,
                            gsl_vector * g)
{
  const size_t n = r->size2;

  size_t i, j;

  for (j = 0; j < n; j++)
    {
      MpIeee sum=  MpIeee( "0" );

      for (i = 0; i <= j; i++)
        {
          sum += gsl_matrix_get (r, i, j) * gsl_vector_get (qtf, i);
        }

      {
        size_t pj = gsl_permutation_get (p, j);
        MpIeee dpj=  gsl_vector_get (diag, pj);

        gsl_vector_set (g, j, sum / dpj);
      }
    }
}

static int
 lmpar(gsl_matrix * r, const gsl_permutation * perm, const gsl_vector * qtf,
       const gsl_vector * diag, MpIeee delta, MpIeee * par_inout,
       gsl_vector * newton, gsl_vector * gradient, gsl_vector * sdiag, 
       gsl_vector * x, gsl_vector * w)
{
  MpIeee dxnorm;MpIeee  gnorm;MpIeee  fp;MpIeee  fp_old;MpIeee  par_lower;MpIeee  par_upper;MpIeee  par_c;

  MpIeee par=  *par_inout;

  size_t iter = 0;

#ifdef DEBUG
  {cout<<"ENTERING lmpar\n";}
#endif


  compute_newton_direction (r, perm, qtf, newton);

#ifdef DEBUG
  {cout<<"newton = ";}
  gsl_vector_fprintf (stdout, newton, "%g");
  {cout<<"\n";}

  {cout<<"diag = ";}
  gsl_vector_fprintf (stdout, diag, "%g");
  {cout<<"\n";}
#endif

  /* Evaluate the function at the origin and test for acceptance of
     the Gauss-Newton direction. */

  dxnorm = scaled_enorm (diag, newton);

  fp = dxnorm - delta;

#ifdef DEBUG
  {cout<<"dxnorm = "<<setiosflags((ios::floatfield))<< dxnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", delta = "<<setiosflags((ios::floatfield))<< delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", fp = "<<setiosflags((ios::floatfield))<< fp;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  if (fp <= MpIeee( "0.1" ) * delta)
    {
      gsl_vector_memcpy (x, newton);
#ifdef DEBUG
      {cout<<"took newton (fp = "<<setiosflags((ios::floatfield))<< fp;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", delta = "<<setiosflags((ios::floatfield))<< delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<")\n";}
#endif

      *par_inout = MpIeee( "0" );

      return GSL_SUCCESS;
    }

  compute_newton_bound (r, newton, dxnorm, perm, diag, w);

  {
    MpIeee wnorm=  enorm (w);
    MpIeee phider=  wnorm * wnorm;

    /* w == zero if r rank-deficient, 
       then set lower bound to zero form MINPACK, lmder.f 
       Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no) */
    if ( wnorm > MpIeee( "0" ) )
      par_lower = fp / (delta * phider);
    else
      par_lower = MpIeee( "0.0" );
  }

#ifdef DEBUG
  {cout<<"par       = "<<setiosflags((ios::floatfield))<< par      ;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"par_lower = "<<setiosflags((ios::floatfield))<< par_lower;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  compute_gradient_direction (r, perm, qtf, diag, gradient);

  gnorm = enorm (gradient);

#ifdef DEBUG
  {cout<<"gradient = ";} gsl_vector_fprintf(stdout, gradient, "%g"); {cout<<"\n";}
  {cout<<"gnorm = "<<setiosflags((ios::floatfield))<< gnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  par_upper =  gnorm / delta;

  if (par_upper == MpIeee( "0" ))
    {
      par_upper = GSL_DBL_MIN / GSL_MIN_DBL(delta, MpIeee( "0.1" ));
    }

#ifdef DEBUG
  {cout<<"par_upper = "<<setiosflags((ios::floatfield))<< par_upper;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  if (par > par_upper)
    {
#ifdef DEBUG
  {cout<<"set par to par_upper\n";}
#endif

      par = par_upper;
    }
  else if (par < par_lower)
    {
#ifdef DEBUG
  {cout<<"set par to par_lower\n";}
#endif

      par = par_lower;
    }

  if (par == MpIeee( "0" ))
    {
      par = gnorm / dxnorm;
#ifdef DEBUG
      {cout<<"set par to gnorm/dxnorm = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

    }

  /* Beginning of iteration */

iteration:

  iter++;

#ifdef DEBUG
  {cout<<"lmpar iteration = "<< iter<<"\n";}
#endif

#ifdef BRIANSFIX
  /* Seems like this is described in the paper but not in the MINPACK code */

  if (par < par_lower || par > par_upper) 
    {
      par = GSL_MAX_DBL (MpIeee( "0.001" ) * par_upper, sqrt(par_lower * par_upper));
    }
#endif

  /* Evaluate the function at the current value of par */

  if (par == MpIeee( "0" ))
    {
      par = GSL_MAX_DBL (MpIeee( "0.001" ) * par_upper, GSL_DBL_MIN);
#ifdef DEBUG
      {cout<<"par = 0, set par to  = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

    }

  /* Compute the least squares solution of [ R P x - Q^T f, sqrt(par) D x]
     for A = Q R P^T */

#ifdef DEBUG
  {cout<<"calling qrsolv with par = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  {
    MpIeee sqrt_par=  sqrt(par);

    qrsolv (r, perm, sqrt_par, diag, qtf, x, sdiag, w);
  }

  dxnorm = scaled_enorm (diag, x);

  fp_old = fp;

  fp = dxnorm - delta;

#ifdef DEBUG
  {cout<<"After qrsolv dxnorm = "<<setiosflags((ios::floatfield))<< dxnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", delta = "<<setiosflags((ios::floatfield))<< delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", fp = "<<setiosflags((ios::floatfield))<< fp;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"sdiag = ";} gsl_vector_fprintf(stdout, sdiag, "%g"); {cout<<"\n";}
  {cout<<"x = ";} gsl_vector_fprintf(stdout, x, "%g"); {cout<<"\n";}
  {cout<<"r = ";} gsl_matrix_fprintf(stdout, r, "%g"); {cout<<"\nXXX\n";}
#endif

  /* If the function is small enough, accept the current value of par */

  if (fabs (fp) <= 0.1 * delta)
    goto line220;

  if (par_lower == MpIeee( "0" ) && fp <= fp_old && fp_old < MpIeee( "0" ))
    goto line220;

  /* Check for maximum number of iterations */

  if (iter == 10)
    goto line220;

  /* Compute the Newton correction */

  compute_newton_correction (r, sdiag, perm, x, dxnorm, diag, w);

#ifdef DEBUG
  {cout<<"newton_correction = ";}
  gsl_vector_fprintf(stdout, w, "%g"); {cout<<"\n";}
#endif

  {
    MpIeee wnorm=  enorm (w);
    par_c = fp / (delta * wnorm * wnorm);
  }

#ifdef DEBUG
  {cout<<"fp = "<<setiosflags((ios::floatfield))<< fp;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"par_lower = "<<setiosflags((ios::floatfield))<< par_lower;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"par_upper = "<<setiosflags((ios::floatfield))<< par_upper;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"par_c = "<<setiosflags((ios::floatfield))<< par_c;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif


  /* Depending on the sign of the function, update par_lower or par_upper */

  if (fp > MpIeee( "0" ))
    {
      if (par > par_lower)
        {
          par_lower = par;
#ifdef DEBUG
      {cout<<"fp > 0: set par_lower = par = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

        }
    }
  else if (fp < MpIeee( "0" ))
    {
      if (par < par_upper)
        {
#ifdef DEBUG
      {cout<<"fp < 0: set par_upper = par = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif
          par_upper = par;
        }
    }

  /* Compute an improved estimate for par */

#ifdef DEBUG
      {cout<<"improved estimate par = MAX("<<setiosflags((ios::floatfield))<< par_lower;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< par+par_c;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<") \n";}
#endif

  par = GSL_MAX_DBL (par_lower, par + par_c);

#ifdef DEBUG
      {cout<<"improved estimate par = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" \n";}
#endif


  goto iteration;

line220:

#ifdef DEBUG
  {cout<<"LEAVING lmpar, par = "<<setiosflags((ios::floatfield))<< par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  *par_inout = par;

  return GSL_SUCCESS;
}



