#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multiroots/test_funcs.c
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#include "test_funcs.h"

/* For information on testing see the following paper,

   J.J More, B.S. Garbow, K.E. Hillstrom, "Testing Unconstrained
   Optimization Software", ACM Transactions on Mathematical Software,
   Vol 7, No 1, (1981) p 17-41

   */

/* Rosenbrock Function */

gsl_multiroot_function_fdf rosenbrock =
{&rosenbrock_f,
 &rosenbrock_df,
 &rosenbrock_fdf,
 2, 0};

void
rosenbrock_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, -1.2);
  gsl_vector_set (x, 1, 1.0);
}

int
 rosenbrock_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee y0=  MpIeee( "1" ) - x0;
  MpIeee y1=  MpIeee( "10" ) * (x1 - x0 * x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 rosenbrock_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));

  MpIeee df00=  -MpIeee( "1" );
  MpIeee df01=  MpIeee( "0" );
  MpIeee df10=  -MpIeee( "20" ) * x0;
  MpIeee df11=  MpIeee( "10" );

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 rosenbrock_fdf(const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * df)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, df);

  return GSL_SUCCESS;
}


/* Freudenstein and Roth function */

gsl_multiroot_function_fdf roth =
{&roth_f,
 &roth_df,
 &roth_fdf,
 2, 0};

void
roth_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 4.5);  /* changed from the value in the paper */
  gsl_vector_set (x, 1, 3.5);  /* otherwise the problem is too hard */
}

int
 roth_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee y0=  -MpIeee( "13.0" ) + x0 + ((MpIeee( "5.0" ) - x1)*x1 - MpIeee( "2.0" ))*x1;
  MpIeee y1=  -MpIeee( "29.0" ) + x0 + ((x1 + MpIeee( "1.0" ))*x1 - MpIeee( "14.0" ))*x1;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 roth_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee df00=  MpIeee( "1" );
  MpIeee df01=  -MpIeee( "3" ) * x1 * x1 + MpIeee( "10" ) * x1 - MpIeee( "2" );
  MpIeee df10=  MpIeee( "1" );
  MpIeee df11=  MpIeee( "3" ) * x1 * x1 + MpIeee( "2" ) * x1 - MpIeee( "14" );

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 roth_fdf(const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * df)
{
  roth_f (x, params, f);
  roth_df (x, params, df);

  return GSL_SUCCESS;
}



/* Powell badly scaled function */

gsl_multiroot_function_fdf powellscal =
{&powellscal_f,
 &powellscal_df,
 &powellscal_fdf,
 2, 0};

void
powellscal_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 0.0);
  gsl_vector_set (x, 1, 1.0);
}

int
 powellscal_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee y0=  MpIeee( "10000.0" ) * x0 * x1 - MpIeee( "1.0" );
  MpIeee y1=  exp (-x0) + exp (-x1) - MpIeee( "1.0001" );

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 powellscal_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee df00=  MpIeee( "10000.0" ) * x1;MpIeee  df01=  MpIeee( "10000.0" ) * x0;
  MpIeee df10=  -exp (-x0);MpIeee  df11=  -exp (-x1);

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 powellscal_fdf(const gsl_vector * x, void *params,
                  gsl_vector * f, gsl_matrix * df)
{
  powellscal_f (x, params, f);
  powellscal_df (x, params, df);

  return GSL_SUCCESS;
}


/* Brown badly scaled function */

gsl_multiroot_function_fdf brownscal =
{&brownscal_f,
 &brownscal_df,
 &brownscal_fdf,
 2, 0};

void
brownscal_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 1.0);
  gsl_vector_set (x, 1, 1.0);
}

int
 brownscal_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee y0=  x0 - MpIeee( "1" )e6;
  MpIeee y1=  x0 * x1 - MpIeee( "2" );

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 brownscal_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee df00=  MpIeee( "1.0" );MpIeee  df01=  MpIeee( "0.0" );
  MpIeee df10=  x1;MpIeee  df11=  x0;

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 brownscal_fdf(const gsl_vector * x, void *params,
                  gsl_vector * f, gsl_matrix * df)
{
  brownscal_f (x, params, f);
  brownscal_df (x, params, df);

  return GSL_SUCCESS;
}


/* Powell Singular Function */

gsl_multiroot_function_fdf powellsing =
{&powellsing_f,
 &powellsing_df,
 &powellsing_fdf,
 4, 0};

void
powellsing_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 3.0);
  gsl_vector_set (x, 1, -1.0);
  gsl_vector_set (x, 2, 0.0);
  gsl_vector_set (x, 3, 1.0);
}

int
 powellsing_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));

  MpIeee y0=  x0 + MpIeee( "10" ) * x1;
  MpIeee y1=  sqrt (MpIeee( "5.0" )) * (x2 - x3);
  MpIeee y2=  pow (x1 - MpIeee( "2" ) * x2, MpIeee( "2.0" ));
  MpIeee y3=  sqrt (MpIeee( "10.0" )) * pow (x0 - x3, MpIeee( "2.0" ));

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 powellsing_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));

  MpIeee df00=  MpIeee( "1" );MpIeee  df01=  MpIeee( "10" );MpIeee  df02=  MpIeee( "0" );MpIeee  df03=  MpIeee( "0" );
  MpIeee df10=  MpIeee( "0" );MpIeee  df11=  MpIeee( "0" );MpIeee  df12=  sqrt (MpIeee( "5.0" ));MpIeee  df13=  -df12;
  MpIeee df20=  MpIeee( "0" );MpIeee  df21=  MpIeee( "2" ) * (x1 - MpIeee( "2" ) * x2);MpIeee  df22=  -MpIeee( "2" ) * df21;MpIeee  df23=  MpIeee( "0" );
  MpIeee df30=  MpIeee( "2" ) * sqrt (MpIeee( "10.0" )) * (x0 - x3);MpIeee  df31=  MpIeee( "0" );MpIeee  df32=  MpIeee( "0" );MpIeee  df33=  -df30;

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 0, 2, df02);
  gsl_matrix_set (df, 0, 3, df03);

  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);
  gsl_matrix_set (df, 1, 2, df12);
  gsl_matrix_set (df, 1, 3, df13);

  gsl_matrix_set (df, 2, 0, df20);
  gsl_matrix_set (df, 2, 1, df21);
  gsl_matrix_set (df, 2, 2, df22);
  gsl_matrix_set (df, 2, 3, df23);

  gsl_matrix_set (df, 3, 0, df30);
  gsl_matrix_set (df, 3, 1, df31);
  gsl_matrix_set (df, 3, 2, df32);
  gsl_matrix_set (df, 3, 3, df33);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 powellsing_fdf(const gsl_vector * x, void *params,
                    gsl_vector * f, gsl_matrix * df)
{
  powellsing_f (x, params, f);
  powellsing_df (x, params, df);

  return GSL_SUCCESS;
}


/* Wood function */

gsl_multiroot_function_fdf wood =
{&wood_f,
 &wood_df,
 &wood_fdf,
 4, 0};

void
wood_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, -3.0);
  gsl_vector_set (x, 1, -1.0);
  gsl_vector_set (x, 2, -3.0);
  gsl_vector_set (x, 3, -1.0);
}

int
 wood_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));

  MpIeee t1=  x1 - x0 * x0;
  MpIeee t2=  x3 - x2 * x2;

  MpIeee y0=  -MpIeee( "200.0" ) * x0 * t1 - (MpIeee( "1" ) - x0);
  MpIeee y1=  MpIeee( "200.0" ) * t1 + MpIeee( "20.2" ) * (x1 - MpIeee( "1" )) + MpIeee( "19.8" ) * (x3 - MpIeee( "1" ));
  MpIeee y2=  -MpIeee( "180.0" ) * x2 * t2 - (MpIeee( "1" ) - x2);
  MpIeee y3=  MpIeee( "180.0" ) * t2 + MpIeee( "20.2" ) * (x3 - MpIeee( "1" )) + MpIeee( "19.8" ) * (x1 - MpIeee( "1" ));

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 wood_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));

  MpIeee t1=  x1 - MpIeee( "3" ) * x0 * x0;
  MpIeee t2=  x3 - MpIeee( "3" ) * x2 * x2;

  MpIeee df00=  -MpIeee( "200.0" ) * t1 + MpIeee( "1" );MpIeee  df01=  -MpIeee( "200.0" ) * x0;MpIeee  df02=  MpIeee( "0" );MpIeee  df03=  MpIeee( "0" );
  MpIeee df10=  -MpIeee( "400.0" )*x0;MpIeee  df11=  MpIeee( "200.0" ) + MpIeee( "20.2" );MpIeee  df12=  MpIeee( "0" );MpIeee  df13=  MpIeee( "19.8" );
  MpIeee df20=  MpIeee( "0" );MpIeee  df21=  MpIeee( "0" );MpIeee  df22=  -MpIeee( "180.0" ) * t2 + MpIeee( "1" );MpIeee  df23=  -MpIeee( "180.0" ) * x2;
  MpIeee df30=  MpIeee( "0" );MpIeee  df31=  MpIeee( "19.8" );MpIeee  df32=  -MpIeee( "2" ) * MpIeee( "180" ) * x2;MpIeee  df33=  MpIeee( "180.0" ) + MpIeee( "20.2" );

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 0, 2, df02);
  gsl_matrix_set (df, 0, 3, df03);

  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);
  gsl_matrix_set (df, 1, 2, df12);
  gsl_matrix_set (df, 1, 3, df13);

  gsl_matrix_set (df, 2, 0, df20);
  gsl_matrix_set (df, 2, 1, df21);
  gsl_matrix_set (df, 2, 2, df22);
  gsl_matrix_set (df, 2, 3, df23);

  gsl_matrix_set (df, 3, 0, df30);
  gsl_matrix_set (df, 3, 1, df31);
  gsl_matrix_set (df, 3, 2, df32);
  gsl_matrix_set (df, 3, 3, df33);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 wood_fdf(const gsl_vector * x, void *params,
                    gsl_vector * f, gsl_matrix * df)
{
  wood_f (x, params, f);
  wood_df (x, params, df);

  return GSL_SUCCESS;
}


/* Helical Valley Function */

gsl_multiroot_function_fdf helical =
{&helical_f,
 &helical_df,
 &helical_fdf,
 3, 0};

void
helical_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, -1.0);
  gsl_vector_set (x, 1, 0.0);
  gsl_vector_set (x, 2, 0.0);
}

int
 helical_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));

  MpIeee t1;MpIeee  t2;
  MpIeee y0;MpIeee  y1;MpIeee  y2;

  if (x0 > MpIeee( "0" )) 
    {
      t1 = atan(x1/x0) / (MpIeee( "2.0" ) * M_PI);
    }
  else if (x0 < MpIeee( "0" ))
    {
      t1 = MpIeee( "0.5" ) + atan(x1/x0) / (MpIeee( "2.0" ) * M_PI);
    }
  else
    {
      t1 = MpIeee( "0.25" ) * (x1 > MpIeee( "0" ) ? +MpIeee( "1" ) : -MpIeee( "1" ));
    }

  t2 = sqrt(x0*x0 + x1*x1) ;
  
  y0 = MpIeee( "10" ) * (x2 - MpIeee( "10" ) * t1);
  y1 = MpIeee( "10" ) * (t2 - MpIeee( "1" ));
  y2 = x2 ;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 helical_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));

  MpIeee t=  x0 * x0 + x1 * x1 ;
  MpIeee t1=  MpIeee( "2" ) * M_PI * t ;
  MpIeee t2=  sqrt(t) ;

  MpIeee df00=  MpIeee( "100" )*x1/t1;MpIeee  df01=  -MpIeee( "100.0" ) * x0/t1;MpIeee  df02=  MpIeee( "10.0" );
  MpIeee df10=  MpIeee( "10" )*x0/t2;MpIeee  df11=  MpIeee( "10" )*x1/t2;MpIeee  df12=  MpIeee( "0" );
  MpIeee df20=  MpIeee( "0" );MpIeee  df21=  MpIeee( "0" );MpIeee  df22=  MpIeee( "1.0" );

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 0, 2, df02);

  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);
  gsl_matrix_set (df, 1, 2, df12);

  gsl_matrix_set (df, 2, 0, df20);
  gsl_matrix_set (df, 2, 1, df21);
  gsl_matrix_set (df, 2, 2, df22);

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 helical_fdf(const gsl_vector * x, void *params,
                    gsl_vector * f, gsl_matrix * df)
{
  helical_f (x, params, f);
  helical_df (x, params, df);

  return GSL_SUCCESS;
}


/* Discrete Boundary Value Function */

#define N 10

gsl_multiroot_function_fdf dbv =
{&dbv_f,
 &dbv_df,
 &dbv_fdf,
 N, 0};

void
dbv_initpt (gsl_vector * x)
{
  size_t i;
  MpIeee h=  MpIeee( "1.0" ) / (N + MpIeee( "1.0" ));

  for (i = 0; i < N; i++)
    {
      MpIeee t=  (i + MpIeee( "1" )) * h;
      MpIeee z=  t * (t - MpIeee( "1" ));
      gsl_vector_set (x, i, z);
    }
}

int
 dbv_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i;

  MpIeee h=  MpIeee( "1.0" ) / (N + MpIeee( "1.0" ));

  for (i = 0; i < N; i++)
    {
      MpIeee z;MpIeee  ti=  (i + MpIeee( "1" )) * h;
      MpIeee xi=  MpIeee( "0" );MpIeee  xim1=  MpIeee( "0" );MpIeee  xip1=  MpIeee( "0" );

      xi = gsl_vector_get (x, i);
      
      if (i > 0)
        xim1 = gsl_vector_get (x, i - MpIeee( "1" ));

      if (i < N - 1)
        xip1 = gsl_vector_get (x, i + MpIeee( "1" ));

      z = MpIeee( "2" ) * xi - xim1 - xip1 + h * h * pow(xi + ti + MpIeee( "1" ), MpIeee( "3.0" )) / MpIeee( "2.0" );

      gsl_vector_set (f, i, z);

    }

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 dbv_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  size_t i, j;

  MpIeee h=  MpIeee( "1.0" ) / (N + MpIeee( "1.0" ));

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      gsl_matrix_set (df, i, j, 0.0);

  for (i = 0; i < N; i++)
    {
      MpIeee dz_dxi;MpIeee  ti=  (i + MpIeee( "1" )) * h;

      MpIeee xi=  gsl_vector_get (x, i);
      
      dz_dxi = MpIeee( "2.0" ) + (MpIeee( "3.0" ) / MpIeee( "2.0" )) * h * h * pow(xi + ti + MpIeee( "1" ), MpIeee( "2.0" )) ;
      
      gsl_matrix_set (df, i, i, dz_dxi);

      if (i > 0)
        gsl_matrix_set (df, i, i-1, -1.0);

      if (i < N - 1)
        gsl_matrix_set (df, i, i+1, -1.0);

    }

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 dbv_fdf(const gsl_vector * x, void *params,
                    gsl_vector * f, gsl_matrix * df)
{
  dbv_f (x, params, f);
  dbv_df (x, params, df);

  return GSL_SUCCESS;
}

/* Trigonometric Function */

gsl_multiroot_function_fdf trig =
{&trig_f,
 &trig_df,
 &trig_fdf,
 N, 0};

void
trig_initpt (gsl_vector * x)
{
  size_t i;

  for (i = 0; i < N; i++)       /* choose an initial point which converges */
    {
      gsl_vector_set (x, i, 0.05);   
    }
}

int
 trig_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i;
  MpIeee sum=  MpIeee( "0" );

  for (i = 0; i < N; i++)
    {
      sum += cos(gsl_vector_get(x,i));
    }

  for (i = 0; i < N; i++)
    {
      MpIeee xi=  gsl_vector_get (x,i);
      MpIeee z=  N - sum + (i + MpIeee( "1" )) * (MpIeee( "1" ) - cos(xi)) - sin(xi);

      gsl_vector_set (f, i, z);
    }

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 trig_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          MpIeee dz;
          MpIeee xi=  gsl_vector_get(x, i);
          MpIeee xj=  gsl_vector_get(x, j);

          if (j == i)
            dz = sin(xi) + (i + MpIeee( "1" )) * sin(xi) - cos(xi);
          else
            dz = sin(xj);
          
          gsl_matrix_set(df, i, j, dz);
        }
    }

  params = 0;                   /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int
 trig_fdf(const gsl_vector * x, void *params,
                    gsl_vector * f, gsl_matrix * df)
{
  trig_f (x, params, f);
  trig_df (x, params, df);

  return GSL_SUCCESS;
}
