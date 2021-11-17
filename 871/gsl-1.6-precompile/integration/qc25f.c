#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/qc25f.c
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

struct fn_fourier_params
{
  gsl_function *function;
  MpIeee omega;
};

static MpIeee fn_sin(MpIeee t, void *params);
static MpIeee fn_cos(MpIeee t, void *params);

static void
qc25f (gsl_function * f, MpIeee a, MpIeee b, 
       gsl_integration_qawo_table * wf, size_t level,
       MpIeee *result, MpIeee *abserr, MpIeee *resabs, MpIeee *resasc);

static void
qc25f (gsl_function * f, MpIeee a, MpIeee b, 
       gsl_integration_qawo_table * wf, size_t level,
       MpIeee *result, MpIeee *abserr, MpIeee *resabs, MpIeee *resasc)
{
  const MpIeee center=  0.5 * (a + b);
  const MpIeee half_length=  0.5 * (b - a);
  const MpIeee omega=  wf->omega ;
  
  const MpIeee par=  omega * half_length;

  if (fabs (par) < 2)
    {
      gsl_function weighted_function;
      struct fn_fourier_params fn_params;

      fn_params.function = f;
      fn_params.omega = omega;

      if (wf->sine == GSL_INTEG_SINE) 
        {
          weighted_function.function = &fn_sin;
        }
      else
        {
          weighted_function.function = &fn_cos;
        }

      weighted_function.params = &fn_params;

      gsl_integration_qk15 (&weighted_function, a, b, result, abserr,
                            resabs, resasc);
      
      return;
    }
  else
    {
      MpIeee *moment;
      MpIeee cheb12[13];MpIeee  cheb24[25];
      MpIeee result_abs;MpIeee  res12_cos;MpIeee  res12_sin;MpIeee  res24_cos;MpIeee  res24_sin;
      MpIeee est_cos;MpIeee  est_sin;
      MpIeee c;MpIeee  s;
      size_t i;

      gsl_integration_qcheb (f, a, b, cheb12, cheb24);

      if (level >= wf->n)
        {
          /* table overflow should not happen, check before calling */
          GSL_ERROR_VOID("table overflow in internal function", GSL_ESANITY);
        }

      /* obtain moments from the table */

      moment = &(wf->chebmo[  25 * level ]); // was eerst wf->chebmo + 25* level ???

      res12_cos = cheb12[12] * moment[12];
      res12_sin = MpIeee( "0" ) ;

      for (i = 0; i < 6; i++)
        {
          size_t k = 10 - 2 * i;
          res12_cos += cheb12[k] * moment[k];
          res12_sin += cheb12[k + 1] * moment[k + 1];
        }

      res24_cos = cheb24[24] * moment[24];
      res24_sin = MpIeee( "0" ) ;

      result_abs = fabs(cheb24[24]) ;

      for (i = 0; i < 12; i++)
        {
          size_t k = 22 - 2 * i;
          res24_cos += cheb24[k] * moment[k];
          res24_sin += cheb24[k + 1] * moment[k + 1];
          result_abs += fabs(cheb24[k]) + fabs(cheb24[k+1]);
        }

      est_cos = fabs(res24_cos - res12_cos);
      est_sin = fabs(res24_sin - res12_sin);

      c = half_length * cos(center * omega);
      s = half_length * sin(center * omega);

      if (wf->sine == GSL_INTEG_SINE)
        {
          *result = c * res24_sin + s * res24_cos;
          *abserr = fabs(c * est_sin) + fabs(s * est_cos);
        }
      else
        {
          *result = c * res24_cos - s * res24_sin;
          *abserr = fabs(c * est_cos) + fabs(s * est_sin);
        }
      
      *resabs = result_abs * half_length;
      *resasc = GSL_DBL_MAX;

      return;
    }
}

static MpIeee fn_sin(MpIeee x, void *params)
{
  struct fn_fourier_params *p = (struct fn_fourier_params *) params;
  gsl_function *f = p->function;
  MpIeee w=  p->omega;
  MpIeee wx=  w * x;
  MpIeee sinwx=  sin(wx) ;
  return GSL_FN_EVAL (f, x) * sinwx;
}

static MpIeee fn_cos(MpIeee x, void *params)
{
  struct fn_fourier_params *p = (struct fn_fourier_params *) params;
  gsl_function *f = p->function;
  MpIeee w=  p->omega;
  MpIeee wx=  w * x;
  MpIeee coswx=  cos(wx) ;
  return GSL_FN_EVAL (f, x) * coswx ;
}

