/**
 * 2014--6--18
 * ode.h
 */
#ifndef _ODE_H_
#define _ODE_H_

/**
 * Requires at least
 * dt - delta t
 * neq - number of equations
 */
struct PARAMETER {
    int neq;
    int nt;
    double ti;
    double tf;
    double dt;
};

struct BUTCHER {
  int S;
  double * b;
  double * c;
  double ** A;
  
};

void rhs(double t, double *u, PARAMETER param, double* f);
/**< rhs function, u'=rhs(t,u) 
    @return (by reference) f rhs(t,u)
    @param t current time step
    @param u solution u at time t
    @param f rhs(t,u)
    @param param structure containing number of equations, number of time steps, initial and final time, time step    
*/

void newt(double t, double *uprev, double * Kguess,
	  double *g, PARAMETER param, BUTCHER rk);
/**< Helper function for advancing the solution from time t(n) to t(n+1) using an implicit RK step on a non linear system using a Newton step.
   @return (by reference) g distance from the root, returned by reference
   @param param structure containing number of equations, number of time steps, initial and final time, time step
   @param t current time step
   @param uprev function value at the previous (known) time step
   @param Kguess iterated guess for the stage values
   @param rk Butcher Tableau coefficients
   @param g distance from the root, returned by reference
*/

void jac(double t, double *uprev, double * Kguess,
	 double *J, PARAMETER param, BUTCHER rk);
/**< Helper function for computing the Jacobian matrix (using finite differences) for advancing the solution from time t(n) to t(n+1) using an implicit RK step on a system of equations
    @return (by reference) J the Jacobian for the Newton step
    @param param structure containing number of equations, number of time steps, initial and final time, time step
    @param t current time step
    @param uprev function value at the previous (known) time step
    @param Kguess iterated guess for the stage values
    @param J Jacobian, returned by reference
    @param rk Butcher Tableau coefficients
*/

void step(double t, double* u, PARAMETER param, double* unew, BUTCHER rk);
/** rhs function, u'=rhs(t,u) 
    @return (by reference) unew solution at time t+dt
    @param t current time step
    @param u solution u at time t
    @param unew solution at time t+dt, returned by reference
    @param param structure containing number of equations, number of time steps, initial and final time, time step    
    @param rk Butcher Tableau coefficients
*/

#endif // _ODE_H_
