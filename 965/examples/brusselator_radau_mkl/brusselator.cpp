#include <stdlib.h>
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK
#include <cmath>
#include "ode.h"


void rhs(double t, double *u, PARAMETER param, double *f) {
  double A = 1.0;
  double B = 3.0;
  double alpha = 0.02;
  double dx = 1.0/(param.neq+1);
  int Nx = param.neq/2;

  f[0] = A+u[0]*u[0]*u[Nx] -(B+1.0)*u[0] +alpha/dx/dx*(u[1]-2*u[0]+1.0);
  f[Nx-1] = A+u[Nx-1]*u[Nx-1]*u[2*Nx-1] -(B+1.0)*u[Nx-1]
    + alpha/dx/dx*(1.0-2*u[Nx-1]+u[Nx-2]);

  f[Nx] = B*u[0]- u[0]*u[0]*u[Nx] + alpha/dx/dx*(u[Nx+1]-2*u[Nx]+3.0);
  f[2*Nx-1] = B*u[Nx-1]- u[Nx-1]*u[Nx-1]*u[2*Nx-1]
    + alpha/dx/dx*(3.0-2*u[2*Nx-1]+u[2*Nx-2]);

  for (int i=1; i<Nx-1; i++) {
    f[i] = A + u[i]*u[i]*u[Nx+i] -(B+1.0)*u[i] +alpha/dx/dx*(u[i+1]-2*u[i]+u[i-1]);
    f[Nx+i] = B*u[i]- u[i]*u[i]*u[Nx+i] + alpha/dx/dx*(u[Nx+i+1]-2*u[Nx+i]+u[Nx+i-1]);
  }
}



void newt(double t, double *uprev, double * Kguess,
	  double *g, PARAMETER param, BUTCHER rk) {

  // note: creation and subsequent deletion of temp variables may lead
  // to an unfair benchmark comparison.
  double * temp = new double[param.neq];
  double * ftemp = new double[param.neq];

  int start_ind;
  for (int ss =0; ss<rk.S; ss++) {

    // initialize
    for (int i =0; i<param.neq; i++) {
      temp[i] = uprev[i];
    }

    // compute temp = y + a_ij*K_j
    for (int s = 0; s<rk.S; s++) {
      double fact = rk.A[ss][s]*param.dt;
      start_ind = s*param.neq;
      for (int i = 0; i<param.neq; i++) {
	temp[i] = temp[i] + fact*Kguess[start_ind+i];
      }
    }

    rhs(t+rk.c[ss]*param.dt,temp,param,ftemp);
    start_ind = ss*param.neq;
    for (int i =0; i<param.neq; i++) {
      g[start_ind+i] = Kguess[start_ind+i]-ftemp[i];
    }
  }

  delete [] temp;
  delete [] ftemp;

}


void jac(double t, double *uprev, double *Kguess,
	 double *J, PARAMETER param, BUTCHER rk) {


  double d = 1e-5; // numerical jacobian approximation

  double *g1;
  double *g2;

  // need to allocate memory
  g1 = new double[rk.S*param.neq];
  g2 = new double[rk.S*param.neq];

  // note, instead of using newt, use rhs for efficiency
  newt(t,uprev,Kguess,g1,param,rk);

  for (int i =0; i<param.neq*rk.S; i++) {
    Kguess[i] = Kguess[i] + d;
    newt(t,uprev,Kguess,g2,param,rk);
    Kguess[i] = Kguess[i] - d;

    for (int j = 0; j < param.neq*rk.S; j++) {
      J[j*param.neq*rk.S + i] = (g2[j]-g1[j])/d;
    }
  }
  

  
  // need to delete memory
  delete [] g1;
  delete [] g2;
}



void step(double t, double* uold, PARAMETER param, double* unew, BUTCHER rk) {

  double NEWTON_TOL = 1.0e-12;
  int NEWTON_MAXSTEP = 10;

  typedef double *pdoub;

  pdoub J;
  pdoub stepsize;

  int neq = param.neq*rk.S;
  stepsize = new double[neq];

  J = new double[neq*neq];

  // initialize memory for stage weights
  double * K = new double [neq];
  double * Ktemp = new double[param.neq];
  
  // try and find a good initial guess 
  double * ftemp = new double[param.neq];
  rhs(t,uold,param,ftemp);
  
  double * temp = new double[param.neq];

  for (int s = 0; s<rk.S; s++){
    
    for (int i =0; i< param.neq; i++) {
      temp[i] = uold[i] + rk.c[s]*param.dt*ftemp[i];
    }

    
    rhs(t+rk.c[s]*param.dt,temp,param,Ktemp);

    for (int i = 0; i < param.neq; i++) {
      K[s*param.neq+i] = Ktemp[i];
    }
  }
    
  delete [] temp;
  delete [] ftemp;
  delete [] Ktemp;

  
  double maxstepsize;

  int counter=0;
  int * pivot = new int[neq];

  //iterate until stages K_i, i=0..s-1 converge
  while (1) {

    // finding zero of nonlinear system.  residual stored in stepsized
    newt(t,uold,K,stepsize,param,rk);


    
    // compute numerical Jacobian
    jac(t,uold,K,J,param,rk);

    int nrhs = 1;
    
    // blas linear solve
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, neq,
		  nrhs, J, neq,pivot,stepsize,1);

    // check for convergence
    maxstepsize = 0.0;
    for (int i = 0; i<neq; i++) {
      if (std::abs(stepsize[i])>maxstepsize) {
	maxstepsize = std::abs(stepsize[i]);
      }
    }

    
    // updated Kguess using stepsize
    for (int i =0; i<neq; i++) {
      K[i] = K[i] - stepsize[i];
    }   

    // if update sufficiently small enough
    if ( maxstepsize < NEWTON_TOL) {
      break;
    }

    counter++;
    //error if too many steps
    if (counter > NEWTON_MAXSTEP) {
      fprintf(stderr,"max newton iterations reached\n");
      exit(42);
    }
  } // end newton iteration


  
  // form solution from stage weights
  for (int i =0; i<param.neq; i++) {
    unew[i] = uold[i];
  }
  
  for (int s = 0; s < rk.S; s++) {
    for (int i=0; i<param.neq; i++) {
      unew[i] = unew[i] + param.dt*rk.b[s]*K[s*param.neq+i];
    }
  }

  
  delete [] pivot;
  delete [] stepsize;
  delete [] J;
  delete [] K;

}
