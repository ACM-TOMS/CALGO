/** 
    @file ridc.h
    @author Ong:Benjamin
    @version Revision 0.2
    @brief header file containing explanation of functions for the RIDC integrator
    @date 2015-09-04
*/


#ifndef _RIDC_H_
#define _RIDC_H_

#include <omp.h>
#include <cmath>
#include <algorithm>

class ODE
{
 public:
  /** number of equations */
  int neq;

  /** number of time steps */
  int nt;

  /** initial time */
  double ti;

  /** final time */
  double tf;

  /** time step */
  double dt;  
  
  virtual void rhs(double t, double *u, double *f) = 0;
  /**< user implemented rhs function, u'=rhs(t,u) 
     @return (by reference) f: rhs(t,u)
     @param t current time step
     @param u solution u at time t
     @param f rhs(t,u)
  */
  
  virtual void step(double t, double *u, double *unew) = 0;
  /**< user implemented step function, for advancing the solution from t to t+dt 
     @return (by reference) unew: solution at time t+dt
     @param t current time step
     @param u solution u at time t
     @param unew solution at time t+dt
  */
  
};


void ridc_fe(ODE *ode, int order, double *sol);
/**< Main explicit ridc loop that initializes variables, integrates
   solution from ti to tf by bootstrapping the step function.
   @return (by reference) sol, the solution at the final time, param.tf
   @param ode abstract class containing parameters and step/rhs functions
   @param order order of the RIDC method (predictor + number of correctors)
   @param sol initial condition of the IVP
*/

void ridc_be(ODE *ode, int order,  double *sol);
/**< Main implicit ridc loop that initializes variables, integrates
   solution from ti to tf by bootstrapping the step function.

   @return (by reference) sol, the solution at the final time, param.tf
   @param ode abstract class containing parameters and step/rhs functions
   @param order order of the RIDC method (predictor + number of correctors)
   @param sol initial condition of the IVP
*/


void lagrange_coeff(double *x, int Nx, int i, double *L);
/**< RIDC helper function -- generates the coefficients for the
   lagrange interpolatory polynomials.

   @return (by reference) L: coefficients for the Lagrange
   intepolatory polynomial.  L is a vector of elements such that p(x)
   = L(0) + L(1)*x + L(2)*x^2 + ...
   @param x quadrature nodes
   @param i the i^{th} Lagrange polynomial
   @param Nx number of quadrature nodes
   @param L coefficients, returned by reference
*/


double get_quad_weight(double *L, int Nx, double a, double b);
/**< RIDC helper function -- generates quadrature weight,
   int(L_{n,i}(x),x=a..b)

   @return quadrature weights

   @param a range of integration
   @param b range of integration
   @param Nx number of quadrature nodes
   @param L coefficients for Lagrange poly, L[0] + L[1]x + L[2]x^2 + ...
*/


void integration_matrices(int Nx, double **S);
/**< RIDC helper function -- constructions the integration matrix
   using get_quad_weight

   @return (by reference) the integration matrix S

   @param Nx number of quadrature nodes
   @param S integration matrix (by reference)

*/

void init_unif_nodes(double *x, int Nx, double a, double b);
/**< RIDC helper function -- initializes uniformly spaced quadrature nodes

   @return (by reference) x: uniformly spaced quadrature nodes

   @param Nx number of quadrature nodes
   @param a range of integration
   @param b range of integration
   @param x quadrature node location (returned by reference)

*/

void corr_fe(ODE * ode,
	     double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             double * unew);
/**< RIDC helper function - solves error equation, updating the
   solution from time t to time t+param.dt.

   @return (by reference) unew: solution at time level t + param.dt

   @param ode abstract class containing parameters and step/rhs functions
   @param uold solution at time level t
   @param fprev matrix containing derivative information from previous steps, previous level
   @param S integration matrix (quadrature weights)
   @param index decides which quadrature weights to use
   @param level determines size of quadrature stencil
   @param t current time iterate
   @param unew solution at the new time level, passed by reference

*/

void corr_be(ODE * ode,
	     double * uold,
             double ** fprev,
             double ** S,
             int index, int level,
             double t,
             double * unew);
/**< RIDC helper function - solves error equation, updating the
   solution from time t to time t+param.dt.

   @return (by reference) unew: solution at time level t + param.dt

   @param ode abstract class containing parameters and step/rhs functions
   @param uold solution at time level t
   @param fprev matrix containing derivative information from previous steps, previous level
   @param S integration matrix (quadrature weights)
   @param index decides which quadrature weights to use
   @param level determines size of quadrature stencil
   @param t current time iterate
   @param unew solution at the new time level, passed by reference

*/

#endif // _RIDC_H_
