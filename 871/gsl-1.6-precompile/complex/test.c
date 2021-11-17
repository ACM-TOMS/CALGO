#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* complex/test.c
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
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

struct f
{
  char *name;
  MpIeee(*f) (gsl_complex z);
  MpIeee x;
  MpIeee y;
  MpIeee fx;
  MpIeee fy;
};

struct fz
{
  char *name;
  gsl_complex (*f) (gsl_complex z);
  MpIeee x;
  MpIeee y;
  MpIeee fx;
  MpIeee fy;
};

struct freal
{
  char *name;
  gsl_complex (*f) (MpIeee x);
  MpIeee x;
  MpIeee fx;
  MpIeee fy;
};


#define FN(x) "gsl_complex_" #x, gsl_complex_ ## x
#define ARG(x,y) x, y
#define RES(x,y) x, y

struct f list[] =
{
#include "results1.h"
  {"", 0, 0, 0, 0, 0}
};


struct fz listz[] =
{
#include "results.h"
  {"", 0, 0, 0, 0, 0}
};

struct freal listreal[] =
{
#include "results_real.h"
  {"", 0, 0, 0, 0}
};


int
main (void)
{
  size_t i = 0;

  gsl_ieee_env_setup();


  for (i = 0 ; i < 10; i++) 
    {
      MpIeee r=  (i - MpIeee( "5.0" )) * MpIeee( "0.3" ) ;
      MpIeee t=  MpIeee( "2.0" ) * M_PI * i / MpIeee( "5" ) ;
      MpIeee x=  r * cos(t);MpIeee  y=  r * sin(t) ;
      gsl_complex z = gsl_complex_polar (r, t) ;
      gsl_test_rel (GSL_REAL(z), x, 10 * GSL_DBL_EPSILON, "gsl_complex_polar real part at (r=%g,t=%g)", r, t);
      
      gsl_test_rel (GSL_IMAG(z), y, 10 * GSL_DBL_EPSILON, "gsl_complex_polar imag part at (r=%g,t=%g)", r, t);
    }
    
    i = 0;

  while (list[i].f)
    {
      struct f t = list[i];
      gsl_complex z = gsl_complex_rect (t.x, t.y);
      MpIeee f=  (t.f) (z);
      gsl_test_rel (f, t.fx, 10 * GSL_DBL_EPSILON, "%s at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  i = 0;

  while (listz[i].f)
    {
      struct fz t = listz[i];
      gsl_complex z = gsl_complex_rect (t.x, t.y);
      gsl_complex fz = (t.f) (z);
      MpIeee fx=  GSL_REAL (fz);MpIeee  fy=  GSL_IMAG (fz);

#ifdef DEBUG
      {cout<<"x = ";} gsl_ieee_fprintf_double (stdout, &t.x); {cout<<"\n";}
      {cout<<"y = ";} gsl_ieee_fprintf_double (stdout, &t.y); {cout<<"\n";}
      {cout<<"fx = ";} gsl_ieee_fprintf_double (stdout, &fx); {cout<<"\n";}
      {cout<<"ex = ";} gsl_ieee_fprintf_double (stdout, &t.fx); {cout<<"\n";}
      {cout<<"fy = ";} gsl_ieee_fprintf_double (stdout, &fy); {cout<<"\n";}
      {cout<<"ey = ";} gsl_ieee_fprintf_double (stdout, &t.fy); {cout<<"\n";}
#endif

      gsl_test_rel (fx, t.fx, 10 * GSL_DBL_EPSILON, "%s real part at (%g,%g)", t.name, t.x, t.y);
      gsl_test_rel (fy, t.fy, 10 * GSL_DBL_EPSILON, "%s imag part at (%g,%g)", t.name, t.x, t.y);
      i++;
    }


  i = 0;

  while (listreal[i].f)
    {
      struct freal t = listreal[i];
      gsl_complex fz = (t.f) (t.x);
      MpIeee fx=  GSL_REAL (fz);MpIeee  fy=  GSL_IMAG (fz);

#ifdef DEBUG
      {cout<<"x = ";} gsl_ieee_fprintf_double (stdout, &t.x); {cout<<"\n";}
      {cout<<"fx = ";} gsl_ieee_fprintf_double (stdout, &fx); {cout<<"\n";}
      {cout<<"ex = ";} gsl_ieee_fprintf_double (stdout, &t.fx); {cout<<"\n";}
      {cout<<"fy = ";} gsl_ieee_fprintf_double (stdout, &fy); {cout<<"\n";}
      {cout<<"ey = ";} gsl_ieee_fprintf_double (stdout, &t.fy); {cout<<"\n";}
#endif

      gsl_test_rel (fx, t.fx, 10 * GSL_DBL_EPSILON, "%s real part at (%g,0)", t.name, t.x);
      gsl_test_rel (fy, t.fy, 10 * GSL_DBL_EPSILON, "%s imag part at (%g,0)", t.name, t.x);
      i++;
    }

  exit (gsl_test_summary ());
}
