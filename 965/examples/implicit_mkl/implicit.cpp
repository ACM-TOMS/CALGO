#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "ridc.h"
#include "implicit.h"

int main(int argc, char *argv[]) {
  int order, nt;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <order> <nt>  >  output_file\n");
    fflush(stdout);
    exit(1);
  } else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
  }

  int neq = 2;
  int ti = 0;
  int tf = 1;
  double dt = (double)(tf - ti)/nt; // compute dt
  
  // initialize ODE variable
  ImplicitMKL *ode = new ImplicitMKL(neq,nt,ti,tf,dt);

  sol = new double[neq];
  // specify initial condition
  for (int i =0; i<neq; i++) {
    sol[i]=1.0;
  }


  // call ridc 
  ridc_be(ode, order, sol);
  
  // output solution to screen
  for (int i = 0; i < neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;

  delete ode;
}
