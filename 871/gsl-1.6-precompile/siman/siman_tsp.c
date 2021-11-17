#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* siman/siman_tsp.c
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
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

/* set up parameters for this simulated annealing run */

#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 2000      /* how many iterations for each T? */
#define STEP_SIZE 1.0           /* max step size in random walk */
#define K 1.0                   /* Boltzmann constant */
#define T_INITIAL 5000.0        /* initial temperature */
#define MU_T 1.002              /* damping factor for temperature */
#define T_MIN 5.0e-1

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                             K, T_INITIAL, MU_T, T_MIN};

struct s_tsp_city {
  const char * name;
  MpIeee lat;MpIeee  longitude;        /* coordinates */
};
typedef struct s_tsp_city Stsp_city;

void prepare_distance_matrix(void);
void exhaustive_search(void);
void print_distance_matrix(void);
MpIeee city_distance(Stsp_city c1, Stsp_city c2);
MpIeee Etsp(void *xp);
MpIeee Mtsp(void *xp, void *yp);
void Stsp(const gsl_rng * r, void *xp, MpIeee step_size);
void Ptsp(void *xp);

/* in this table, latitude and longitude are obtained from the US
   Census Bureau, at http://www.census.gov/cgi-bin/gazetteer */

Stsp_city cities[] = {{"Santa Fe",    35.68,   105.95},
                      {"Phoenix",     33.54,   112.07},
                      {"Albuquerque", 35.12,   106.62},
                      {"Clovis",      34.41,   103.20},
                      {"Durango",     37.29,   107.87},
                      {"Dallas",      32.79,    96.77},
                      {"Tesuque",     35.77,   105.92},
                      {"Grants",      35.15,   107.84},
                      {"Los Alamos",  35.89,   106.28},
                      {"Las Cruces",  32.34,   106.76},
                      {"Cortez",      37.35,   108.58},
                      {"Gallup",      35.52,   108.74}};

#define N_CITIES (sizeof(cities)/sizeof(Stsp_city))

MpIeee distance_matrix[N_CITIES][N_CITIES];

/* distance between two cities */
MpIeee city_distance(Stsp_city c1, Stsp_city c2)
{
  const MpIeee earth_radius=  6375.000; /* 6000KM approximately */
  /* sin and cos of lat and long; must convert to radians */
  MpIeee sla1=  sin(c1.lat*M_PI/MpIeee( "180" ));MpIeee  cla1=  cos(c1.lat*M_PI/MpIeee( "180" ));MpIeee 
    slo1=  sin(c1.longitude*M_PI/MpIeee( "180" ));MpIeee  clo1=  cos(c1.longitude*M_PI/MpIeee( "180" ));
  MpIeee sla2=  sin(c2.lat*M_PI/MpIeee( "180" ));MpIeee  cla2=  cos(c2.lat*M_PI/MpIeee( "180" ));MpIeee 
    slo2=  sin(c2.longitude*M_PI/MpIeee( "180" ));MpIeee  clo2=  cos(c2.longitude*M_PI/MpIeee( "180" ));

  MpIeee x1=  cla1*clo1;
  MpIeee x2=  cla2*clo2;

  MpIeee y1=  cla1*slo1;
  MpIeee y2=  cla2*slo2;

  MpIeee z1=  sla1;
  MpIeee z2=  sla2;

  MpIeee dot_product=  x1*x2 + y1*y2 + z1*z2;

  MpIeee angle=  acos(dot_product);

  /* distance is the angle (in radians) times the earth radius */
  return angle*earth_radius;
}

/* energy for the travelling salesman problem */
MpIeee Etsp(void *xp)
{
  /* an array of N_CITIES integers describing the order */
  int  *route=  (int *) xp;
  MpIeee E=  MpIeee( "0" );
  unsigned int  i;

  for (i = 0; i < N_CITIES; ++i) {
    /* use the distance_matrix to optimize this calculation; it had
       better be allocated!! */
    E += distance_matrix[route[i]][route[(i + 1) % N_CITIES]];
  }

  return E;
}

MpIeee Mtsp(void *xp, void *yp)
{
  int  *route1=  (int *) xp;int   *route2=  (int *) yp;
  MpIeee distance=  MpIeee( "0" );
  unsigned int  i;

  for (i = 0; i < N_CITIES; ++i) {
    distance += ((route1[i] == route2[i]) ? MpIeee( "0" ) : MpIeee( "1" ));
  }

  return distance;
}

/* take a step through the TSP space */
void Stsp(const gsl_rng * r, void *xp, MpIeee step_size)
{
  int  x1;int   x2;int   dummy;
  int  *route=  (int *) xp;

  step_size = MpIeee( "0" ) ; /* prevent warnings about unused parameter */

  /* pick the two cities to swap in the matrix; we leave the first
     city fixed */
  x1 = (gsl_rng_get (r) % (N_CITIES-1)) + 1;
  do {
    x2 = (gsl_rng_get (r) % (N_CITIES-1)) + 1;
  } while (x2 == x1);

  dummy = route[x1];
  route[x1] = route[x2];
  route[x2] = dummy;
}

void Ptsp(void *xp)
{
  unsigned int  i;
  int  *route=  (int *) xp;
  {cout<<"  [";}
  for (i = 0; i < N_CITIES; ++i) {
    {cout<<" "<< route[i]<<" ";}
  }
  {cout<<"]  ";}
}

int main(void)
{
  int  x_initial[N_CITIES];
  unsigned int  i;

  const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;

  gsl_ieee_env_setup ();

  prepare_distance_matrix();

  /* set up a trivial initial route */
  {cout<<"# initial order of cities:\n";}
  for (i = 0; i < N_CITIES; ++i) {
    {cout<<"# \""<< cities[i].name<<"\"\n";}
    x_initial[i] = i;
  }

  {cout<<"# distance matrix is:\n";}
  print_distance_matrix();

  {cout<<"# initial coordinates of cities (longitude and latitude)\n";}
  /* this can be plotted with */
  /* ./siman_tsp > hhh ; grep city_coord hhh | awk '{print $2 "   " $3}' | xyplot -ps -d "xy" > c.eps */
  for (i = 0; i < N_CITIES+1; ++i) {
    {cout<<"###initial_city_coord: "<<setiosflags((ios::floatfield))<<
           -cities[x_initial[i % N_CITIES]].longitude;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<<
           cities[x_initial[i % N_CITIES]].lat;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" \""<<
           cities[x_initial[i % N_CITIES]].name<<"\"\n";}
  }

/*   exhaustive_search(); */

  gsl_siman_solve(r, x_initial, Etsp, Stsp, Mtsp, Ptsp, NULL, NULL, NULL,
                  N_CITIES*sizeof(int), params);

  {cout<<"# final order of cities:\n";}
  for (i = 0; i < N_CITIES; ++i) {
    {cout<<"# \""<< cities[x_initial[i]].name<<"\"\n";}
  }

  {cout<<"# final coordinates of cities (longitude and latitude)\n";}
  /* this can be plotted with */
  /* ./siman_tsp > hhh ; grep city_coord hhh | awk '{print $2 "   " $3}' | xyplot -ps -d "xy" > c.eps */
  for (i = 0; i < N_CITIES+1; ++i) {
    {cout<<"###final_city_coord: "<<setiosflags((ios::floatfield))<<
           -cities[x_initial[i % N_CITIES]].longitude;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<<
           cities[x_initial[i % N_CITIES]].lat;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<
           cities[x_initial[i % N_CITIES]].name<<"\n";}
  }

  {cout<<"# ";}
  fflush(stdout);
#if 0
  system("date");
#endif /* 0 */
  fflush(stdout);

  return 0;
}

void prepare_distance_matrix()
{
  unsigned int  i;int   j;
  MpIeee dist;

  for (i = 0; i < N_CITIES; ++i) {
    for (j = 0; j < N_CITIES; ++j) {
      if (i == j) {
        dist = MpIeee( "0" );
      } else {
        dist = city_distance(cities[i], cities[j]);
      }
      distance_matrix[i][j] = dist;
    }
  }
}

void print_distance_matrix()
{
  unsigned int  i;int   j;

  for (i = 0; i < N_CITIES; ++i) {
    {cout<<"# ";}
    for (j = 0; j < N_CITIES; ++j) {
      {cout<<""<<setiosflags((ios::fixed & ios::floatfield))<<setw(15)<<setprecision(8)<< distance_matrix[i][j];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"   ";}
    }
    {cout<<"\n";}
  }
}

/* [only works for 12] search the entire space for solutions */
static MpIeee best_E=  MpIeee( "1.0e100" );MpIeee  second_E=  MpIeee( "1.0e100" );MpIeee  third_E=  MpIeee( "1.0e100" );
static int  best_route[N_CITIES];
static int  second_route[N_CITIES];
static int  third_route[N_CITIES];
static void do_all_perms(int  *route, int  n);

void exhaustive_search()
{
  static int  initial_route[N_CITIES] =  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  {cout<<"\n# ";}
  fflush(stdout);
#if 0
  system("date");
#endif
  fflush(stdout);
  do_all_perms(initial_route, 1);
  {cout<<"\n# ";}
  fflush(stdout);
#if 0
  system("date");
#endif /* 0 */
  fflush(stdout);

  {cout<<"# exhaustive best route: ";}
  Ptsp(best_route);
  {cout<<"\n# its energy is: "<<setiosflags((ios::floatfield))<< best_E;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  {cout<<"# exhaustive second_best route: ";}
  Ptsp(second_route);
  {cout<<"\n# its energy is: "<<setiosflags((ios::floatfield))<< second_E;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  {cout<<"# exhaustive third_best route: ";}
  Ptsp(third_route);
  {cout<<"\n# its energy is: "<<setiosflags((ios::floatfield))<< third_E;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
}

/* James Theiler's recursive algorithm for generating all routes */
static void do_all_perms(int  *route, int  n)
{
  if (n == (N_CITIES-1)) {
    /* do it! calculate the energy/cost for that route */
    MpIeee E;
    E = Etsp(route);            /* TSP energy function */
    /* now save the best 3 energies and routes */
    if (E < best_E) {
      third_E = second_E;
      memcpy(third_route, second_route, N_CITIES*sizeof(*route));
      second_E = best_E;
      memcpy(second_route, best_route, N_CITIES*sizeof(*route));
      best_E = E;
      memcpy(best_route, route, N_CITIES*sizeof(*route));
    } else if (E < second_E) {
      third_E = second_E;
      memcpy(third_route, second_route, N_CITIES*sizeof(*route));
      second_E = E;
      memcpy(second_route, route, N_CITIES*sizeof(*route));
    } else if (E < third_E) {
      third_E = E;
      memcpy(route, third_route, N_CITIES*sizeof(*route));
    }
  } else {
    int  new_route[N_CITIES];
    unsigned int  j;
    int  swap_tmp;
    memcpy(new_route, route, N_CITIES*sizeof(*route));
    for (j = n; j < N_CITIES; ++j) {
      swap_tmp = new_route[j];
      new_route[j] = new_route[n];
      new_route[n] = swap_tmp;
      do_all_perms(new_route, n+1);
    }
  }
}
