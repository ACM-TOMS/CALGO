#include "ridc.h"
#include <stdio.h>
#include "mkl.h" // MKL Standard Lib
#include "mkl_lapacke.h" // MKL LAPACK

#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_

using namespace std;

class ImplicitMKL : public ODE {
public:
  ImplicitMKL(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt) {
    neq = my_neq; 
    nt = my_nt; 
    ti = my_ti;
    tf = my_tf;
    dt = my_dt;
  }
  
  void rhs(double t, double *u, double *f) {
    /** user implemented rhs function, u'=rhs(t,u) 
       @return (by reference) f: rhs(t,u)
       @param t current time step
       @param u solution u at time t
       @param f rhs(t,u)
    */
    for (int i =0; i<neq; i++) {
      f[i]=-(i+1)*t*u[i];
    }    
  }

  void step(double t, double * u, double * unew) {
  /** user implemented step function, for advancing the solution from t to t+dt 
     @return (by reference) unew: solution at time t+dt
     @param t current time step
     @param u solution u at time t
     @param unew solution at time t+dt
  */
    
    double tnew = t + dt;

    double NEWTON_TOL = 1.0e-14;
    int NEWTON_MAXSTEP = 1000;

    double * J;
    double * stepsize;
    double * uguess;
    
    stepsize = new double[neq];
    uguess = new double[neq];

    J = new double[neq*neq];

    // set initial condition
    for (int i=0; i<neq; i++) {
      uguess[i] = u[i];
    }

    double maxstepsize;
    
    int counter=0;
    int * pivot = new int[neq];

    while (1) {

      // create rhs
      newt(tnew,u,uguess,stepsize);

      // compute numerical Jacobian
      jac(tnew,uguess,J);

      // blas linear solve
      LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) neq,
		    (lapack_int) 1, J, neq,pivot,stepsize,1);
      
      // check for convergence
      maxstepsize = 0.0;
      for (int i = 0; i<neq; i++) {
	uguess[i] = uguess[i] - stepsize[i];
	if (std::abs(stepsize[i])>maxstepsize) {
	  maxstepsize = std::abs(stepsize[i]);
	}
      }

      // if update sufficiently small enough
      if ( maxstepsize < NEWTON_TOL) {
	break;
      }

      counter++;
      //error if too many steps
      if (counter > NEWTON_MAXSTEP) {
	fprintf(stderr,"max Newton iterations reached\n");
	exit(42);
      }
    } // end Newton iteration
    
    for (int i=0; i<neq; i++) {
      unew[i] = uguess[i];
    }

    delete [] pivot;
    delete [] uguess;
    delete [] stepsize;
    delete [] J;
    
  }

 private:
  void newt(double t, double *uprev, double *uguess,
	    double *g){
    /**  Helper function for computing the next Newton step
       for solving a system of equations
       
       @return (by reference) g, how far from zero we are
       @param t current time step
       @param uguess current solution guess
       @param uprev solution at previous time step
       @param g how far from zero we are, returned by reference
    */
    rhs(t,uguess,g);
    for (int i =0; i<neq; i++) {
      g[i] = uguess[i]-uprev[i]-dt*g[i];
    }
  }
    
  void jac(double t, double *u, double *J){
    /** Helper function for computing the jacobian matrix (using
       finite differences) for advancing the solution from time t(n)
       to t(n+1) using an implicit Euler step on a system of equations
       
       @return (by reference) J the Jacobian for the Newton step
       @param t current time step
       @param u function value at the current time step
       @param J Jacobian, returned by reference
    */
    double d = 1e-5; // numerical Jacobian approximation
    double *u1;
    double *f1;
    double *f;
    
    
    // need to allocate memory
    u1 = new double[neq];
    f1 = new double[neq];
    f = new double[neq];
    
    // note, instead of using newt, use rhs for efficiency
    rhs(t,u,f);
    
    for (int i = 0; i<neq; i++) {
      for (int j = 0; j<neq; j++) {
	u1[j] = u[j];
      }
      u1[i] = u1[i] + d;
      
      rhs(t,u1,f1);
      
      for (int j = 0; j<neq; j++) {
	J[j*neq+i] = -dt*(f1[j]-f[j])/d;
      }
      J[i*neq+i] = 1.0+J[i*neq+i];
    }
    
    // need to delete memory
    delete [] u1;
    delete [] f1;
    delete [] f;
  }
  
};






#endif // _IMPLICIT_H_
