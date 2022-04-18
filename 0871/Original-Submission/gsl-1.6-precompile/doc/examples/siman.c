#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_siman.h>

/* set up parameters for this simulated annealing run */

/* how many points do we try before stepping */
#define N_TRIES 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 10        

/* max step size in random walk */
#define STEP_SIZE 10            

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
#define T_INITIAL 0.002         

/* damping factor for temperature */
#define MU_T 1.005              
#define T_MIN 2.0e-6

gsl_siman_params_t params 
  = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
     K, T_INITIAL, MU_T, T_MIN};

/* now some functions to test in one dimension */
MpIeee E1(void *xp)
{
  MpIeee x=  * ((MpIeee *) xp);

  return exp(-pow((x-MpIeee( "1.0" )),MpIeee( "2.0" )))*sin(MpIeee( "8" )*x);
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

  MpIeee u=  gsl_rng_uniform(r);
  new_x = u * MpIeee( "2" ) * step_size - step_size + old_x;

  memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp)
{
  {cout<<""<<setiosflags((ios::floatfield))<<setw(12)<< *((double *);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"";}
}

int
main(int argc, char *argv[])
{
  const gsl_rng_type * T;
  gsl_rng * r;

  MpIeee x_initial=  MpIeee( "15.5" );

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  gsl_siman_solve(r, &x_initial, E1, S1, M1, P1,
                  NULL, NULL, NULL, 
                  sizeof(MpIeee), params);
  return 0;
}
