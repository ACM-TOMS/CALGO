#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static void
chop_small_elements (gsl_vector * d, gsl_vector * f)
{
  const size_t N = d->size;
  MpIeee d_i=  gsl_vector_get (d, 0);

  size_t i;

  for (i = 0; i < N - 1; i++)
    {
      MpIeee f_i=  gsl_vector_get (f, i);
      MpIeee d_ip1=  gsl_vector_get (d, i + 1);

      if (fabs (f_i) < GSL_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1)))
        {
          gsl_vector_set (f, i, 0.0);
        }
      d_i = d_ip1;
    }

}

static MpIeee trailing_eigenvalue(const gsl_vector * d, const gsl_vector * f)
{
  const size_t n = d->size;

  MpIeee da=  gsl_vector_get (d, n - 2);
  MpIeee db=  gsl_vector_get (d, n - 1);
  MpIeee fa=  (n > MpIeee( "2" )) ? gsl_vector_get (f, n - 3) : MpIeee( "0.0" );
  MpIeee fb=  gsl_vector_get (f, n - 2);

  MpIeee ta=  da * da + fa * fa;
  MpIeee tb=  db * db + fb * fb;
  MpIeee tab=  da * fb;

  MpIeee dt=  (ta - tb) / MpIeee( "2.0" );

  MpIeee mu;

  if (dt >= MpIeee( "0" ))
    {
      mu = tb - (tab * tab) / (dt + hypot (dt, tab));
    }
  else 
    {
      mu = tb + (tab * tab) / ((-dt) + hypot (dt, tab));
    }

  return mu;
}

static void
create_schur (MpIeee d0, MpIeee f0, MpIeee d1, MpIeee * c, MpIeee * s)
{
  MpIeee apq=  MpIeee( "2.0" ) * d0 * f0;
  
  if (apq != MpIeee( "0.0" ))
    {
      MpIeee t;
      MpIeee tau=  (f0*f0 + (d1 + d0)*(d1 - d0)) / apq;
      
      if (tau >= MpIeee( "0.0" ))
        {
          t = MpIeee( "1.0" )/(tau + hypot(MpIeee( "1.0" ), tau));
        }
      else
        {
          t = -MpIeee( "1.0" )/(-tau + hypot(MpIeee( "1.0" ), tau));
        }

      *c = MpIeee( "1.0" ) / hypot(MpIeee( "1.0" ), t);
      *s = t * (*c);
    }
  else
    {
      *c = MpIeee( "1.0" );
      *s = MpIeee( "0.0" );
    }
}

static void
svd2 (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V)
{
  size_t i;
  MpIeee c;MpIeee  s;MpIeee  a11;MpIeee  a12;MpIeee  a21;MpIeee  a22;

  const size_t M = U->size1;
  const size_t N = V->size1;

  MpIeee d0=  gsl_vector_get (d, 0);
  MpIeee f0=  gsl_vector_get (f, 0);
  
  MpIeee d1=  gsl_vector_get (d, 1);

  if (d0 == MpIeee( "0.0" ))
    {
      /* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

      create_givens (f0, d1, &c, &s);

      /* compute B <= G^T B X,  where X = [0,1;1,0] */

      gsl_vector_set (d, 0, c * f0 - s * d1);
      gsl_vector_set (f, 0, s * f0 + c * d1);
      gsl_vector_set (d, 1, 0.0);
      
      /* Compute U <= U G */

      for (i = 0; i < M; i++)
        {
          MpIeee Uip=  gsl_matrix_get (U, i, 0);
          MpIeee Uiq=  gsl_matrix_get (U, i, 1);
          gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
        }

      /* Compute V <= V X */

      gsl_matrix_swap_columns (V, 0, 1);

      return;
    }
  else if (d1 == MpIeee( "0.0" ))
    {
      /* Eliminate off-diagonal element in [d0,f0;0,0] */

      create_givens (d0, f0, &c, &s);

      /* compute B <= B G */

      gsl_vector_set (d, 0, d0 * c - f0 * s);
      gsl_vector_set (f, 0, 0.0);

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
        {
          MpIeee Vip=  gsl_matrix_get (V, i, 0);
          MpIeee Viq=  gsl_matrix_get (V, i, 1);
          gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
          gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
        }

      return;
    }
  else
    {
      /* Make columns orthogonal, A = [d0, f0; 0, d1] * G */
      
      create_schur (d0, f0, d1, &c, &s);
      
      /* compute B <= B G */
      
      a11 = c * d0 - s * f0;
      a21 = - s * d1;
      
      a12 = s * d0 + c * f0;
      a22 = c * d1;
      
      /* Compute V <= V G */
      
      for (i = 0; i < N; i++)
        {
          MpIeee Vip=  gsl_matrix_get (V, i, 0);
          MpIeee Viq=  gsl_matrix_get (V, i, 1);
          gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
          gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
        }
      
      /* Eliminate off-diagonal elements, bring column with largest
         norm to first column */
      
      if (hypot(a11, a21) < hypot(a12,a22))
        {
          MpIeee t1;MpIeee  t2;

          /* B <= B X */

          t1 = a11; a11 = a12; a12 = t1;
          t2 = a21; a21 = a22; a22 = t2;

          /* V <= V X */

          gsl_matrix_swap_columns(V, 0, 1);
        } 

      create_givens (a11, a21, &c, &s);
      
      /* compute B <= G^T B */
      
      gsl_vector_set (d, 0, c * a11 - s * a21);
      gsl_vector_set (f, 0, c * a12 - s * a22);
      gsl_vector_set (d, 1, s * a12 + c * a22);
      
      /* Compute U <= U G */
      
      for (i = 0; i < M; i++)
        {
          MpIeee Uip=  gsl_matrix_get (U, i, 0);
          MpIeee Uiq=  gsl_matrix_get (U, i, 1);
          gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
        }

      return;
    }
}


static void
chase_out_intermediate_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * U, size_t k0)
{
  const size_t M = U->size1;
  const size_t n = d->size;
  MpIeee c;MpIeee  s;
  MpIeee x;MpIeee  y;
  size_t i, k;

  x = gsl_vector_get (f, k0);
  y = gsl_vector_get (d, k0+1);

  for (k = k0; k < n - 1; k++)
    {
      create_givens (y, -x, &c, &s);
      
      /* Compute U <= U G */
      
      for (i = 0; i < M; i++)
        {
          MpIeee Uip=  gsl_matrix_get (U, i, k0);
          MpIeee Uiq=  gsl_matrix_get (U, i, k + 1);
          gsl_matrix_set (U, i, k0, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
        }
      
      /* compute B <= G^T B */
      
      gsl_vector_set (d, k + 1, s * x + c * y);

      if (k == k0)
        gsl_vector_set (f, k, c * x - s * y );

      if (k < n - 2) 
        {
          MpIeee z=  gsl_vector_get (f, k + 1);
          gsl_vector_set (f, k + 1, c * z); 

          x = -s * z ;
          y = gsl_vector_get (d, k + 2); 
        }
    }
}

static void
chase_out_trailing_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * V)
{
  const size_t N = V->size1;
  const size_t n = d->size;
  MpIeee c;MpIeee  s;
  MpIeee x;MpIeee  y;
  size_t i, k;

  x = gsl_vector_get (d, n - 2);
  y = gsl_vector_get (f, n - 2);

  for (k = n - 1; k > 0 && k--;)
    {
      create_givens (x, y, &c, &s);

      /* Compute V <= V G where G = [c, s ; -s, c] */
      
      for (i = 0; i < N; i++)
        {
          MpIeee Vip=  gsl_matrix_get (V, i, k);
          MpIeee Viq=  gsl_matrix_get (V, i, n - 1);
          gsl_matrix_set (V, i, k, c * Vip - s * Viq);
          gsl_matrix_set (V, i, n - 1, s * Vip + c * Viq);
        }

      /* compute B <= B G */
      
      gsl_vector_set (d, k, c * x - s * y);

      if (k == n - 2)
        gsl_vector_set (f, k, s * x + c * y );

      if (k > 0) 
        {
          MpIeee z=  gsl_vector_get (f, k - 1);
          gsl_vector_set (f, k - 1, c * z); 

          x = gsl_vector_get (d, k - 1); 
          y = s * z ;
        }
    }
}

static void
qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V)
{
  const size_t M = U->size1;
  const size_t N = V->size1;
  const size_t n = d->size;
  MpIeee y;MpIeee  z;
  MpIeee ak;MpIeee  bk;MpIeee  zk;MpIeee  ap;MpIeee  bp;MpIeee  aq;MpIeee  bq;
  size_t i, k;

  if (n == 1)
    return;  /* shouldn't happen */

  /* Compute 2x2 svd directly */

  if (n == 2)
    {
      svd2 (d, f, U, V);
      return;
    }

  /* Chase out any zeroes on the diagonal */

  for (i = 0; i < n - 1; i++)
    {
      MpIeee d_i=  gsl_vector_get (d, i);
      
      if (d_i == MpIeee( "0.0" ))
        {
          chase_out_intermediate_zero (d, f, U, i);
          return;
        }
    }

  /* Chase out any zero at the end of the diagonal */

  {
    MpIeee d_nm1=  gsl_vector_get (d, n - 1);

    if (d_nm1 == MpIeee( "0.0" )) 
      {
        chase_out_trailing_zero (d, f, V);
        return;
      }
  }


  /* Apply QR reduction steps to the diagonal and offdiagonal */

  {
    MpIeee d0=  gsl_vector_get (d, 0);
    MpIeee f0=  gsl_vector_get (f, 0);
    
    MpIeee d1=  gsl_vector_get (d, 1);
    MpIeee f1=  gsl_vector_get (f, 1);
    
    {
      MpIeee mu=  trailing_eigenvalue (d, f);
    
      y = d0 * d0 - mu;
      z = d0 * f0;
    }
    
    /* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    
    ak = MpIeee( "0" );
    bk = MpIeee( "0" );
    
    ap = d0;
    bp = f0;
    
    aq = d1;
    bq = f1;
  }

  for (k = 0; k < n - 1; k++)
    {
      MpIeee c;MpIeee  s;
      create_givens (y, z, &c, &s);

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
        {
          MpIeee Vip=  gsl_matrix_get (V, i, k);
          MpIeee Viq=  gsl_matrix_get (V, i, k + 1);
          gsl_matrix_set (V, i, k, c * Vip - s * Viq);
          gsl_matrix_set (V, i, k + 1, s * Vip + c * Viq);
        }

      /* compute B <= B G */

      {
        MpIeee bk1=  c * bk - s * z;

        MpIeee ap1=  c * ap - s * bp;
        MpIeee bp1=  s * ap + c * bp;
        MpIeee zp1=  -s * aq;

        MpIeee aq1=  c * aq;

        if (k > 0)
          {
            gsl_vector_set (f, k - 1, bk1);
          }

        ak = ap1;
        bk = bp1;
        zk = zp1;

        ap = aq1;

        if (k < n - 2)
          {
            bp = gsl_vector_get (f, k + 1);
          }
        else
          {
            bp = MpIeee( "0.0" );
          }

        y = ak;
        z = zk;
      }

      create_givens (y, z, &c, &s);

      /* Compute U <= U G */

      for (i = 0; i < M; i++)
        {
          MpIeee Uip=  gsl_matrix_get (U, i, k);
          MpIeee Uiq=  gsl_matrix_get (U, i, k + 1);
          gsl_matrix_set (U, i, k, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
        }

      /* compute B <= G^T B */

      {
        MpIeee ak1=  c * ak - s * zk;
        MpIeee bk1=  c * bk - s * ap;
        MpIeee zk1=  -s * bp;

        MpIeee ap1=  s * bk + c * ap;
        MpIeee bp1=  c * bp;

        gsl_vector_set (d, k, ak1);

        ak = ak1;
        bk = bk1;
        zk = zk1;

        ap = ap1;
        bp = bp1;

        if (k < n - 2)
          {
            aq = gsl_vector_get (d, k + 2);
          }
        else
          {
            aq = MpIeee( "0.0" );
          }

        y = bk;
        z = zk;
      }
    }

  gsl_vector_set (f, n - 2, bk);
  gsl_vector_set (d, n - 1, ap);
}


