#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* siman/siman.c
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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>

/* implementation of a basic simulated annealing algorithm */

void 
gsl_siman_solve (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                 gsl_siman_step_t take_step,
                 gsl_siman_metric_t distance,
                 gsl_siman_print_t print_position,
                 gsl_siman_copy_t copyfunc,
                 gsl_siman_copy_construct_t copy_constructor,
                 gsl_siman_destroy_t destructor,
                 size_t element_size,
                 gsl_siman_params_t params)
{
  void *x, *new_x, *best_x;
  MpIeee E;MpIeee  new_E;MpIeee  best_E;
  int  i;int   done;
  MpIeee T;
  int  n_evals=  1;int   n_iter=  0;int   n_accepts;int   n_rejects;int   n_eless;

  /* this function requires that either the dynamic functions (copy,
     copy_constructor and destrcutor) are passed, or that an element
     size is given */
  assert((copyfunc != NULL && copy_constructor != NULL && destructor != NULL)
         || (element_size != 0));

  distance = 0 ; /* This parameter is not currently used */
  E = Ef(x0_p);

  if (copyfunc) {
    x = copy_constructor(x0_p);
    new_x = copy_constructor(x0_p);
    best_x = copy_constructor(x0_p);
  } else {
    x = (void *) malloc (element_size);
    memcpy (x, x0_p, element_size);
    new_x = (void *) malloc (element_size);
    best_x =  (void *) malloc (element_size);
    memcpy (best_x, x0_p, element_size);
  }

  best_E = E;

  T = params.t_initial;
  done = 0;

  if (print_position) {
    {cout<<"#-iter  #-evals   temperature     position   energy\n";}
  }

  while (!done) {

    n_accepts = 0;
    n_rejects = 0;
    n_eless = 0;
    for (i = 0; i < params.iters_fixed_T; ++i) {
      if (copyfunc) {
        copyfunc(x, new_x);
      } else {
        memcpy (new_x, x, element_size);
      }

      take_step (r, new_x, params.step_size);
      new_E = Ef (new_x);

      if(new_E <= best_E){
        if (copyfunc) {
          copyfunc(new_x,best_x);
        } else {
          memcpy (best_x, new_x, element_size);
        }
        best_E=new_E;
      }

      ++n_evals;                /* keep track of Ef() evaluations */
      /* now take the crucial step: see if the new point is accepted
         or not, as determined by the boltzman probability */
      if (new_E < E) {
        /* yay! take a step */
        if (copyfunc) {
          copyfunc(new_x, x);
        } else {
          memcpy (x, new_x, element_size);
        }
        E = new_E;
        ++n_eless;
      } else if (gsl_rng_uniform(r) < exp (-(new_E - E)/(params.k * T)) ) {
        /* yay! take a step */
        if (copyfunc) {
          copyfunc(new_x, x);
        } else {
          memcpy(x, new_x, element_size);
        }
        E = new_E;
        ++n_accepts;
      } else {
        ++n_rejects;
      }
    }

    if (print_position) {
      /* see if we need to print stuff as we go */
      /*       printf("%5d %12g %5d %3d %3d %3d", n_iter, T, n_evals, */
      /*           100*n_eless/n_steps, 100*n_accepts/n_steps, */
      /*           100*n_rejects/n_steps); */
      {cout<<""<<setw(5)<< n_iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"   "<<setw(7)<< n_evals;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"  "<<setiosflags((ios::floatfield))<<setw(12)<< T;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"";}
      print_position (x);
      {cout<<"  "<<setiosflags((ios::floatfield))<<setw(12)<< E;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }

    /* apply the cooling schedule to the temperature */
    /* FIXME: I should also introduce a cooling schedule for the iters */
    T /= params.mu_t;
    ++n_iter;
    if (T < params.t_min) {
      done = 1;
    }
  }

  /* at the end, copy the result onto the initial point, so we pass it
     back to the caller */
  if (copyfunc) {
    copyfunc(best_x, x0_p);
  } else {
    memcpy (x0_p, best_x, element_size);
  }

  if (copyfunc) {
    destructor(x);
    destructor(new_x);
    destructor(best_x);
  } else {
    free (x);
    free (new_x);
    free (best_x);
  }
}

/* implementation of a simulated annealing algorithm with many tries */

void 
gsl_siman_solve_many (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                      gsl_siman_step_t take_step,
                      gsl_siman_metric_t distance,
                      gsl_siman_print_t print_position,
                      size_t element_size,
                      gsl_siman_params_t params)
{
  /* the new set of trial points, and their energies and probabilities */
  void *x, *new_x;
  MpIeee *energies;MpIeee  *probs;MpIeee  *sum_probs;
  MpIeee Ex;                    /* energy of the chosen point */
  MpIeee T;                     /* the temperature */
  int  i;int   done;
  MpIeee u;                     /* throw the die to choose a new "x" */
  int  n_iter;

  if (print_position) {
    {cout<<"#-iter    temperature       position";}
    {cout<<"         delta_pos        energy\n";}
  }

  x = (void *) malloc (params.n_tries * element_size);
  new_x = (void *) malloc (params.n_tries * element_size);
  energies = (MpIeee *) malloc (params.n_tries * sizeof (MpIeee));
  probs = (MpIeee *) malloc (params.n_tries * sizeof (MpIeee));
  sum_probs = (MpIeee *) malloc (params.n_tries * sizeof (MpIeee));

  T = params.t_initial;
/*    memcpy (x, x0_p, element_size); */
  memcpy (x, x0_p, element_size);
  done = 0;

  n_iter = 0;
  while (!done)
    {
      Ex = Ef (x);
      for (i = 0; i < params.n_tries - 1; ++i)
        {                       /* only go to N_TRIES-2 */
          /* center the new_x[] around x, then pass it to take_step() */
          sum_probs[i] = MpIeee( "0" );
          memcpy ((char *)new_x + i * element_size, x, element_size);
          take_step (r, (char *)new_x + i * element_size, params.step_size);
          energies[i] = Ef ((char *)new_x + i * element_size);
          probs[i] = exp (-(energies[i] - Ex) / (params.k * T));
        }
      /* now add in the old value of "x", so it is a contendor */
      memcpy ((char *)new_x + (params.n_tries - 1) * element_size, x, element_size);
      energies[params.n_tries - 1] = Ex;
      probs[params.n_tries - 1] = exp (-(energies[i] - Ex) / (params.k * T));

      /* now throw biased die to see which new_x[i] we choose */
      sum_probs[0] = probs[0];
      for (i = 1; i < params.n_tries; ++i)
        {
          sum_probs[i] = sum_probs[i - 1] + probs[i];
        }
      u = gsl_rng_uniform (r) * sum_probs[params.n_tries - 1];
      for (i = 0; i < params.n_tries; ++i)
        {
          if (u < sum_probs[i])
            {
              memcpy (x, (char *)new_x + i * element_size, element_size);
              break;
            }
        }
      if (print_position)
        {
          {cout<<""<<setw(5)<< n_iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"\t"<<setiosflags((ios::floatfield))<<setw(12)<< T;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\t";}
          print_position (x);
          {cout<<"\t"<<setiosflags((ios::floatfield))<<setw(12)<< distance (x, x0_p);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\t"<<setiosflags((ios::floatfield))<<setw(12)<< Ex;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
        }
      T /= params.mu_t;
      ++n_iter;
      if (T < params.t_min)
        {
          done = 1;
        }
    }

  /* now return the value via x0_p */
  memcpy (x0_p, x, element_size);

  /*  printf("the result is: %g (E=%g)\n", x, Ex); */

  free (x);
  free (new_x);
  free (energies);
  free (probs);
  free (sum_probs);
}
