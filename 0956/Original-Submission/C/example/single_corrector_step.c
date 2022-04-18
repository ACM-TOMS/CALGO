#include "KS_example.h"
int
single_corrector_step (int N_real, double *Z, double *T)
/**********************************************************************/
/* This function computes a single Newton corrector step for the BVP  */
/*                                                                    */
/* -c \psi' + \psi \psi' + \psi'' + \nu\psi^{(iv)} + A\sin(\psi) = 0, */
/* \psi(0) = \psi(2\pi).                                              */
/*                                                                    */
/* The BVP is discretized with a trigonometric polynomial basis, i.e.,*/
/* \psi(\xi) = \sum_{k=0}^{N_{grid}-1} X_{k} \exp(i k \xi)            */
/* and using collocation on a uniform grid to determine the           */
/* coefficients X_k. That is, the collocation points are taken to be  */
/* \xi_\ell=2\pi\ell/N_{grid}. The resulting nonlinear vector field   */
/* can be evaluated efficiently using discrete Fourier transforms     */
/* (hence the use of complex FFTs from the GSL library).              */
/*                                                                    */
/* This BVP arises in seeking travelling wave solutions of the form   */
/* \psi(x-ct)for the time-dependent PDE                               */
/*                                                                    */
/* u_t + u u_x + u_{xx} + \nu u_{xxxx} + A \sin(u) = 0, 0<x<1,        */
/* u(0,t) = u(2*\pi,t) for all t                                      */
/*                                                                    */
/* (a modified Kuramoto-Shivashinsky equation with periodic BCs).     */
/*                                                                    */
/* The input vector Z consists of N_{real}=2 N_{grid} + 4 components: */
/* in order, the components of Z are the real and imaginary parts of  */
/* X_{k}, the Fourier coefficients in alternating order, the wave     */
/* speed c (real and imaginary parts) and the viscosity \nu (real and */
/* imaginary parts).                                                  */
/*                                                                    */
/* To clarify, the components of Z are as follows:                    */
/*                                                                    */
/* Z[0] = Re(X[0]) <- 0th Fourier coeff.                              */
/* Z[1] = Im(X[0]) <- 0th Fourier coeff.                              */
/* Z[2] = Re(X[1]) <- 1st Fourier coeff.                              */
/* Z[3] = Im(X[1]) <- 1st Fourier coeff.                              */
/*   ...                                                              */
/* Z[2*N_grid-2] = Re(X[N_grid-1]) <- (N_grid-1)st Fourier coeff.     */
/* Z[2*N_grid-1] = Im(X[N_grid-1]) <- (N_grid-1)st Fourier coeff.     */
/* Z[2*N_grid]   = Re(c)   <- wavespeed                               */
/* Z[2*N_grid+1] = Im(c)   <- wavespeed                               */
/* Z[2*N_grid+2] = Re(nu)  <- viscosity                               */
/* Z[2*N_grid+3] = Im(nu)  <- viscosity                               */
/*                                                                    */
/* The input array T represents the tangent vector to the curve on    */
/* which Z lies. The array of doubles T is also ordered using         */
/* compressed complex storage. The result of the computation modifies */
/* the array Z in place, i.e., Z is overwritten upon return from      */
/* invoking this function.                                            */
/**********************************************************************/
{
  bool is_valid;
  int N_grid, k, status;
  complex c, nu;
  /* N_real = number of REAL variables input    */
  /* N_grid = number of collocation points used */
  N_grid = (N_real / 2) - 2;

  for (k = 0; k < N_grid; k++)
    X_cplx[k] = Z[2*k] + 1.0I * Z[2*k + 1];
  c = Z[N_real - 4] + 1.0I*Z[N_real-3];
  nu = Z[N_real - 2] + 1.0I*Z[N_real-1];

  is_valid = compute_Jacobian (N_real, c, nu, T);
  if (!is_valid) {
	fprintf (stderr, "%s: Error\n", __func__);
	return -1;
  }
  /* Compute the right-hand side and pack into a complex vector */
  status = compute_residual ( N_real, Z, real_workspace );
  if (status!=0) {
	fprintf (stderr, "%s: Error in computing residual\n", __func__);
	fprintf (stderr, "%s: error status = %i\n", __func__, status);
	return status;
  }

  for (k=0; k<N_grid; k++)
    RHS_cplx[k] = real_workspace[2*k] + 1.0I * real_workspace[2*k+1];
  RHS_cplx[N_grid] = 0.0;
  RHS_cplx[N_grid+1] = 0.0;

  /* Solve linear system Jac_cplx*X=RHS_cplx; solution overwrites RHS_cplx */
  status = clapack_zgesv( CblasColMajor, N_grid+2, 1, Jac_cplx, N_grid+2,
                        ipiv, RHS_cplx, N_grid+2);
  if (status!= 0) {
	fprintf (stderr, "%s: Error in clapack_zgesv\n", __func__);
	fprintf (stderr, "%s: error status = %i\n", __func__, status);
	return -1;
  }
  /* Filter out the part of the Newton update that
     gives complex-valued solutions or parameters
     due to accumulated numerical noise */
  RHS_cplx[0] = creal(RHS_cplx[0]);
  RHS_cplx[N_grid/2] = creal(RHS_cplx[N_grid/2]);
  RHS_cplx[N_grid] = creal(RHS_cplx[N_grid]);
  RHS_cplx[N_grid+1] = creal(RHS_cplx[N_grid+1]);
  for (k=1; k<N_grid/2; k++) {
    RHS_cplx[k] = 0.5*(RHS_cplx[k]+conj(RHS_cplx[N_grid-k]));
    RHS_cplx[N_grid-k] = conj(RHS_cplx[k]);
  };
  /* Increment array Z using RHS computed; unwrap real & imaginary parts
     (notice minus sign) */
  for (k=0; k<N_grid+2; k++) {
    Z[2*k] -= creal(RHS_cplx[k]);
    Z[2*k+1] -= cimag(RHS_cplx[k]);
  }
  return 0;
}
