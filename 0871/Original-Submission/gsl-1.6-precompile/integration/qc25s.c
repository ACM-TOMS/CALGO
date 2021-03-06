#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/qc25s.c
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

struct fn_qaws_params
{
  gsl_function *function;
  MpIeee a;
  MpIeee b;
  gsl_integration_qaws_table *table;
};

static MpIeee fn_qaws(MpIeee t, void *params);
static MpIeee fn_qaws_L(MpIeee x, void *params);
static MpIeee fn_qaws_R(MpIeee x, void *params);

static void
compute_result (const MpIeee * r, const MpIeee * cheb12, const MpIeee * cheb24,
                MpIeee * result12, MpIeee * result24);


static void
qc25s (gsl_function * f, MpIeee a, MpIeee b, MpIeee a1, MpIeee b1,
       gsl_integration_qaws_table * t,
       MpIeee *result, MpIeee *abserr, int  *err_reliable);

static void
qc25s (gsl_function * f, MpIeee a, MpIeee b, MpIeee a1, MpIeee b1,
       gsl_integration_qaws_table * t,
       MpIeee *result, MpIeee *abserr, int  *err_reliable)
{
  gsl_function weighted_function;
  struct fn_qaws_params fn_params;
  
  fn_params.function = f;
  fn_params.a = a;
  fn_params.b = b;
  fn_params.table = t;

  weighted_function.params = &fn_params;
    
  if (a1 == a && (t->alpha != MpIeee( "0.0" ) || t->mu != MpIeee( "0" )))
    {
      MpIeee cheb12[13];MpIeee  cheb24[25];

      MpIeee factor=  pow(MpIeee( "0.5" ) * (b1 - a1), t->alpha + MpIeee( "1.0" ));

      weighted_function.function = &fn_qaws_R;

      gsl_integration_qcheb (&weighted_function, a1, b1, cheb12, cheb24);

      if (t->mu == 0)
        {
          MpIeee res12=  MpIeee( "0" );MpIeee  res24=  MpIeee( "0" );
          MpIeee u=  factor;

          compute_result (t->ri, cheb12, cheb24, &res12, &res24);

          *result = u * res24;
          *abserr = fabs(u * (res24 - res12));
        }
      else 
        {
          MpIeee res12a=  MpIeee( "0" );MpIeee  res24a=  MpIeee( "0" );
          MpIeee res12b=  MpIeee( "0" );MpIeee  res24b=  MpIeee( "0" );

          MpIeee u=  factor * log(b1 - a1);
          MpIeee v=  factor;

          compute_result (t->ri, cheb12, cheb24, &res12a, &res24a);
          compute_result (t->rg, cheb12, cheb24, &res12b, &res24b);

          *result = u * res24a + v * res24b;
          *abserr = fabs(u * (res24a - res12a)) + fabs(v * (res24b - res12b));
        }

      *err_reliable = 0;

      return;
    }
  else if (b1 == b && (t->beta != MpIeee( "0.0" ) || t->nu != MpIeee( "0" )))
    {
      MpIeee cheb12[13];MpIeee  cheb24[25];
      MpIeee factor=  pow(MpIeee( "0.5" ) * (b1 - a1), t->beta + MpIeee( "1.0" ));

      weighted_function.function = &fn_qaws_L;

      gsl_integration_qcheb (&weighted_function, a1, b1, cheb12, cheb24);

      if (t->nu == 0)
        {
          MpIeee res12=  MpIeee( "0" );MpIeee  res24=  MpIeee( "0" );
          MpIeee u=  factor;

          compute_result (t->rj, cheb12, cheb24, &res12, &res24);

          *result = u * res24;
          *abserr = fabs(u * (res24 - res12));
        }
      else 
        {
          MpIeee res12a=  MpIeee( "0" );MpIeee  res24a=  MpIeee( "0" );
          MpIeee res12b=  MpIeee( "0" );MpIeee  res24b=  MpIeee( "0" );

          MpIeee u=  factor * log(b1 - a1);
          MpIeee v=  factor;

          compute_result (t->rj, cheb12, cheb24, &res12a, &res24a);
          compute_result (t->rh, cheb12, cheb24, &res12b, &res24b);

          *result = u * res24a + v * res24b;
          *abserr = fabs(u * (res24a - res12a)) + fabs(v * (res24b - res12b));
        }

      *err_reliable = 0;

      return;
    }
  else
    {
      MpIeee resabs;MpIeee  resasc;

      weighted_function.function = &fn_qaws;
  
      gsl_integration_qk15 (&weighted_function, a1, b1, result, abserr,
                            &resabs, &resasc);

      if (*abserr == resasc)
        {
          *err_reliable = 0;
        }
      else 
        {
          *err_reliable = 1;
        }

      return;
    }

}

static MpIeee fn_qaws(MpIeee x, void *params)
{
  struct fn_qaws_params *p = (struct fn_qaws_params *) params;
  gsl_function *f = p->function;
  gsl_integration_qaws_table *t = p->table;

  MpIeee factor=  MpIeee( "1.0" );
  
  if (t->alpha != 0.0)
    factor *= pow(x - p->a, t->alpha);

  if (t->beta != 0.0)
    factor *= pow(p->b - x, t->beta);

  if (t->mu == 1)
    factor *= log(x - p->a);

  if (t->nu == 1)
    factor *= log(p->b - x);

  return factor * GSL_FN_EVAL (f, x);
}

static MpIeee fn_qaws_L(MpIeee x, void *params)
{
  struct fn_qaws_params *p = (struct fn_qaws_params *) params;
  gsl_function *f = p->function;
  gsl_integration_qaws_table *t = p->table;

  MpIeee factor=  MpIeee( "1.0" );
  
  if (t->alpha != 0.0)
    factor *= pow(x - p->a, t->alpha);

  if (t->mu == 1)
    factor *= log(x - p->a);

  return factor * GSL_FN_EVAL (f, x);
}

static MpIeee fn_qaws_R(MpIeee x, void *params)
{
  struct fn_qaws_params *p = (struct fn_qaws_params *) params;
  gsl_function *f = p->function;
  gsl_integration_qaws_table *t = p->table;

  MpIeee factor=  MpIeee( "1.0" );
  
  if (t->beta != 0.0)
    factor *= pow(p->b - x, t->beta);

  if (t->nu == 1)
    factor *= log(p->b - x);

  return factor * GSL_FN_EVAL (f, x);
}


static void
compute_result (const MpIeee * r, const MpIeee * cheb12, const MpIeee * cheb24,
                MpIeee * result12, MpIeee * result24)
{
  size_t i;
  MpIeee res12=  MpIeee( "0" );
  MpIeee res24=  MpIeee( "0" );
  
  for (i = 0; i < 13; i++)
    {
      res12 += r[i] * cheb12[i];
    }
  
  for (i = 0; i < 25; i++)
    {
      res24 += r[i] * cheb24[i];
    }
  
  *result12 = res12;
  *result24 = res24;
}
