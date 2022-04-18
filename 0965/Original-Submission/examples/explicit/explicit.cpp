#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

#include "ridc.h"

using namespace std;

class ExplicitOde : public ODE {
public:
  ExplicitOde(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt) {
    neq = my_neq;
    nt = my_nt;
    ti = my_ti;
    tf = my_tf;
    dt = my_dt;
  }
  
  void rhs(double t, double *u, double *f) {
    for (int i =0; i<neq; i++) {
      f[i]=-(i+1)*t*u[i];
    }    
  }

  void step(double t, double * u, double * unew) {
    double* fold = new double[neq];
    rhs(t,u,fold);
  
    for (int i = 0; i < neq; i++)  {
      unew[i] = u[i] + dt*(fold[i]);
    }
    delete [] fold;
  }
};


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
  ExplicitOde *ode = new ExplicitOde(neq,nt,ti,tf,dt);

  sol = new double[neq];
  // specify initial condition
  for (int i =0; i<neq; i++) {
    sol[i]=1.0;
  }

  // call ridc 
  ridc_fe(ode,order, sol);

  // output solution to screen
  for (int i = 0; i < neq; i++)
    printf("%14.12f\n", sol[i]);
  delete [] sol;
  delete ode;
}
