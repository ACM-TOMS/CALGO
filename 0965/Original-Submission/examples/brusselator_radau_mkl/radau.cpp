#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK
#include <cmath>
#include "ode.h"


/** @brief This is the the main function for the brusselator_radau_mkl example
 *
 * This will pass user given options along with some standard options 
 * for this type of problem in to the PARAMETER struct and start 
 * the solving the problem
 */

int main(int argc, char *argv[]) {
  int nt,neq;
  double *sol;

  if (argc != 3) {
    printf("usage: <executable> <nt> <Nx> >  output_file\n");
    fflush(stdout);
    exit(1);
  }
  else {
    nt = atoi(argv[1]); // number of time steps
    neq = atoi(argv[2]); // number of spatial intervals
  }


  
  // initialize parameters for problem
  PARAMETER param;
  param.ti = 0; // initial time
  param.tf = 10.0; // final time
  param.dt = (double)(param.tf - param.ti)/nt; // compute dt
  param.neq = neq; // number of equations
  param.nt = nt; // store number of time steps


  sol = new double[param.neq];
  // specify initial condition

  
  double xi;
  int Nx=param.neq/2;
  double dx = 1.0/(Nx+1);
 
  for (int i =0; i<Nx; i++) {
    xi = (i+1)*dx;
    sol[i]=1.0 + sin(2*3.14159265359*xi);
    sol[Nx+i] = 3.0;
  }

  // Specify Butcher Tableau for RADAU IIA method, 5th order
  BUTCHER rk;
  
  rk.S = 3; // number of stages
  rk.c  = new double[rk.S];
  rk.b = new double[rk.S];
  rk.A = new double*[rk.S];
  for (int s = 0; s<rk.S; s++) {
    rk.A[s] = new double[rk.S];
  }

  rk.c[0] = 0.155051025721682;
  rk.c[1] = 0.644948974278318;
  rk.c[2] = 1.000000000000000;
  rk.A[0][0] = 0.196815477223660;
  rk.A[1][0] = 0.394424314739087;
  rk.A[2][0] = 0.376403062700467;
  rk.A[0][1] = -0.065535425850198;
  rk.A[1][1] = 0.292073411665228;
  rk.A[2][1] = 0.512485826188422;
  rk.A[0][2] = 0.023770974348220;
  rk.A[1][2] = -0.041548752125998;
  rk.A[2][2] = 0.111111111111111;
  rk.b[0] = rk.A[2][0];
  rk.b[1] = rk.A[2][1];
  rk.b[2] = rk.A[2][2];

  double * unew = new double[param.neq];
  
  for (int ii = 0; ii < param.nt; ii++) {
    // advance using RADAU IIA 5th order method
    double t = ii*param.dt;
    step(t,sol,param,unew,rk);

    for (int i =0; i<param.neq; i++) {
      sol[i] = unew[i];
    }
  }


  // output solution to screen
  for (int i = 0; i < param.neq; i++)
    printf("%14.12f\n", sol[i]);

  delete [] sol;

  // delete memory for radau method
  delete [] rk.c;
  delete [] rk.b;
  for (int s = 0; s < rk.S; s++) {
    delete [] rk.A[s];
  }
  delete [] rk.A;
  
}
