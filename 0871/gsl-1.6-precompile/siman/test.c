#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* siman/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
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
#include <string.h>
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>
#include <stdio.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 1000      /* how many iterations for each T? */
#define STEP_SIZE 1.0           /* max step size in random walk */
#define K 1.0                   /* Boltzmann constant */
#define T_INITIAL 0.008         /* initial temperature */
#define MU_T 1.003              /* damping factor for temperature */
#define T_MIN 2.0e-6

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                             K, T_INITIAL, MU_T, T_MIN};

inline MpIeee square(MpIeee x) ;
inline MpIeee square(MpIeee x) { return x * x ; } 

MpIeee E1(void *xp);
MpIeee M1(void *xp, void *yp);
void S1(const gsl_rng * r, void *xp, MpIeee step_size);
void P1(void *xp);

/* now some functions to test in one dimension */
MpIeee E1(void *xp)
{
  MpIeee x=  * ((MpIeee *) xp);

  return exp(-square(x-MpIeee( "1" )))*sin(MpIeee( "8" )*x) - exp(-square(x-MpIeee( "1000" )))*MpIeee( "0.89" );
}

MpIeee M1(void *xp, void *yp)
{
  MpIeee x=  *((MpIeee *) xp);
  MpIeee y=  *((MpIeee *) yp);

  return fabs(x - y);
}

void S1(const gsl_rng * r, void *xp, MpIeee step_size)
{
  MpIeee old_x=  *((MpIeee *) xp);
  MpIeee new_x;

  new_x = gsl_rng_uniform(r)*MpIeee( "2" )*step_size - step_size + old_x;

  memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp)
{
  {cout<<" "<<setiosflags((ios::floatfield))<<setw(12)<< *((double *);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" ";}
}

int main(void)
{
  MpIeee x_min=  MpIeee( "1.36312999455315182" ) ;
  MpIeee x;

  gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;

  gsl_ieee_env_setup ();

  /* The function tested here has multiple mimima. 
     The global minimum is at    x = 1.36312999, (f = -0.87287)
     There is a local minimum at x = 0.60146196, (f = -0.84893) */

  x = -MpIeee( "10.0" ) ;
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
                  sizeof(MpIeee), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=-10") ;

  x = +MpIeee( "10.0" ) ;
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
                  sizeof(MpIeee), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=10") ;

  /* Start at the false minimum */

  x = +MpIeee( "0.6" ) ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
                  sizeof(MpIeee), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.6") ;

  x = +MpIeee( "0.5" ) ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
                  sizeof(MpIeee), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.5") ;

  x = +MpIeee( "0.4" ) ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
                  sizeof(MpIeee), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.4") ;

  gsl_rng_free(r);
  exit (gsl_test_summary ());

#ifdef JUNK 
  x0.D1 = 12.0;
  {cout<<"#one dimensional problem, x0 = "<<setiosflags((ios::fixed & ios::floatfield))<< x0.D1;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  gsl_siman_Usolve(r, &x0, test_E_1D, test_step_1D, distance_1D,
                   print_pos_1D, params);


  x0.D2[0] = 12.0;
  x0.D2[1] = 5.5;
  {cout<<"#two dimensional problem, (x0,y0) = ("<<setiosflags((ios::fixed & ios::floatfield))<<
         x0.D2[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<","<<setiosflags((ios::fixed & ios::floatfield))<< x0.D2[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<")\n";}
  gsl_siman_Usolve(r, &x0, test_E_2D, test_step_2D, distance_2D,
                   print_pos_2D, params); 

  x0.D3[0] = 12.2;
  x0.D3[1] = 5.5;
  x0.D3[2] = -15.5;
  {cout<<"#three dimensional problem, (x0,y0,z0) = ("<<setiosflags((ios::fixed & ios::floatfield))<<
         x0.D3[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<","<<setiosflags((ios::fixed & ios::floatfield))<< x0.D3[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<","<<setiosflags((ios::fixed & ios::floatfield))<< x0.D3[2];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<")\n";}
  gsl_siman_Usolve(r, &x0, test_E_3D, test_step_3D, distance_3D, 
                   print_pos_3D, params); 

  x0.D2[0] = 12.2;
  x0.D2[1] = 5.5;

  gsl_siman_solve(r, &x0, test_E_2D, test_step_2D, distance_2D, print_pos_2D, params);
  
  x0.D3[0] = 12.2;
  x0.D3[1] = 5.5;
  x0.D3[2] = -15.5;

  gsl_siman_solve(r, &x0, test_E_3D, test_step_3D, distance_3D, print_pos_3D, params);

  return 0;
#endif
}



