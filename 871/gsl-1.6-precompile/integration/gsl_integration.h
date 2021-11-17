#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/gsl_integration.h
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

#ifndef __GSL_INTEGRATION_H__
#define __GSL_INTEGRATION_H__
#include <stdlib.h>
#include <gsl/gsl_math.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


//added for MpIeee
MpIeee myEpsilon();
MpIeee myMax();




/* Workspace for adaptive integrators */

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    MpIeee *alist;
    MpIeee *blist;
    MpIeee *rlist;
    MpIeee *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (const size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);


/* Workspace for QAWS integrator */

typedef struct
{
  MpIeee alpha;
  MpIeee beta;
  int  mu;
  int  nu;
  MpIeee ri[25];
  MpIeee rj[25];
  MpIeee rg[25];
  MpIeee rh[25];
}
gsl_integration_qaws_table;

gsl_integration_qaws_table * 
gsl_integration_qaws_table_alloc (MpIeee alpha, MpIeee beta, int  mu, int  nu);

int
 gsl_integration_qaws_table_set(gsl_integration_qaws_table * t,
                                MpIeee alpha, MpIeee beta, int  mu, int  nu);

void
gsl_integration_qaws_table_free (gsl_integration_qaws_table * t);

/* Workspace for QAWO integrator */

enum gsl_integration_qawo_enum { GSL_INTEG_COSINE, GSL_INTEG_SINE };

typedef struct
{
  size_t n;
  MpIeee omega;
  MpIeee L;
  MpIeee par;
  enum gsl_integration_qawo_enum sine;
  MpIeee *chebmo;
}
gsl_integration_qawo_table;

gsl_integration_qawo_table * 
gsl_integration_qawo_table_alloc (MpIeee omega, MpIeee L, 
                                  enum gsl_integration_qawo_enum sine,
                                  size_t n);

int
 gsl_integration_qawo_table_set(gsl_integration_qawo_table * t,
                                MpIeee omega, MpIeee L,
                                enum gsl_integration_qawo_enum sine);

int
 gsl_integration_qawo_table_set_length(gsl_integration_qawo_table * t,
                                       MpIeee L);

void
gsl_integration_qawo_table_free (gsl_integration_qawo_table * t);


/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
                                   MpIeee a, MpIeee b,
                                   MpIeee *result, MpIeee *abserr,
                                   MpIeee *defabs, MpIeee *resabs);

void gsl_integration_qk15 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qk21 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qk31 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qk41 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qk51 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qk61 (const gsl_function * f, MpIeee a, MpIeee b,
                           MpIeee *result, MpIeee *abserr,
                           MpIeee *resabs, MpIeee *resasc);

void gsl_integration_qcheb (gsl_function * f, MpIeee a, MpIeee b, 
                            MpIeee *cheb12, MpIeee *cheb24);

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
  {
    GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,      /* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
  };

void 
gsl_integration_qk (const int n, const MpIeee xgk[], 
                    const MpIeee wg[], const MpIeee wgk[],
                    MpIeee fv1[], MpIeee fv2[],
                    const gsl_function *f, MpIeee a, MpIeee b,
                    MpIeee * result, MpIeee * abserr, 
                    MpIeee * resabs, MpIeee * resasc);


int  gsl_integration_qng(const gsl_function * f,
                         MpIeee a, MpIeee b,
                         MpIeee epsabs, MpIeee epsrel,
                         MpIeee *result, MpIeee *abserr,
                         size_t * neval);

int  gsl_integration_qag(const gsl_function * f,
                         MpIeee a, MpIeee b,
                         MpIeee epsabs, MpIeee epsrel, size_t limit,
                         int  key,
                         gsl_integration_workspace * workspace,
                         MpIeee *result, MpIeee *abserr);

int  gsl_integration_qagi(gsl_function * f,
                          MpIeee epsabs, MpIeee epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          MpIeee *result, MpIeee *abserr);

int  gsl_integration_qagiu(gsl_function * f,
                           MpIeee a,
                           MpIeee epsabs, MpIeee epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           MpIeee *result, MpIeee *abserr);

int  gsl_integration_qagil(gsl_function * f,
                           MpIeee b,
                           MpIeee epsabs, MpIeee epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           MpIeee *result, MpIeee *abserr);


int  gsl_integration_qags(const gsl_function * f,
                          MpIeee a, MpIeee b,
                          MpIeee epsabs, MpIeee epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          MpIeee *result, MpIeee *abserr);

int  gsl_integration_qagp(const gsl_function * f,
                          MpIeee *pts, size_t npts,
                          MpIeee epsabs, MpIeee epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          MpIeee *result, MpIeee *abserr);

int  gsl_integration_qawc(gsl_function *f,
                          const MpIeee a, const MpIeee b, const MpIeee c,
                          const MpIeee epsabs, const MpIeee epsrel, const size_t limit,
                          gsl_integration_workspace * workspace,
                          MpIeee * result, MpIeee * abserr);

int  gsl_integration_qaws(gsl_function * f,
                          const MpIeee a, const MpIeee b,
                          gsl_integration_qaws_table * t,
                          const MpIeee epsabs, const MpIeee epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          MpIeee *result, MpIeee *abserr);

int  gsl_integration_qawo(gsl_function * f,
                          const MpIeee a,
                          const MpIeee epsabs, const MpIeee epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_qawo_table * wf,
                          MpIeee *result, MpIeee *abserr);

int  gsl_integration_qawf(gsl_function * f,
                          const MpIeee a,
                          const MpIeee epsabs,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_workspace * cycle_workspace,
                          gsl_integration_qawo_table * wf,
                          MpIeee *result, MpIeee *abserr);

__END_DECLS

#endif /* __GSL_INTEGRATION_H__ */
