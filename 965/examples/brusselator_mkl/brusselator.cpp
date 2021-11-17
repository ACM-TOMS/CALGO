#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK
#include <cmath>

#include "ridc.h"
#include "brusselator.h"

int main(int argc, char *argv[]) {
  int order, nt,neq;
  double *sol;

  if (argc != 4) {
    printf("usage: <executable> <order> <nt> <neq> >  output_file\n");
    fflush(stdout);
    exit(1);
  } else {
    order = atoi(argv[1]); // order of method
    nt = atoi(argv[2]); // number of time steps
    neq = atoi(argv[3]); // number of equations
  }

  int ti = 0;
  int tf = 10.0;
  double dt = (double)(tf - ti)/nt; // compute dt
  
  // initialize ODE variable
  Brusselator_MKL *ode = new Brusselator_MKL(neq,nt,ti,tf,dt);

  sol = new double[neq];
  // specify initial condition

  double xi;
  int Nx=neq/2;
  double dx = 1.0/(Nx+1);
  
  for (int i =0; i<Nx; i++) {
    xi = (i+1)*dx;
    sol[i]=1.0 + sin(2*3.14159265359*xi);
    sol[Nx+i] = 3.0;
  }


  // call ridc 
  ridc_be(ode,order, sol);

  // output solution to screen
  for (int i = 0; i < neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;

  delete ode;
}
