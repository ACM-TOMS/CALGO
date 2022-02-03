#include "KS_example.h"
/**********************************************************************/
int
compute_residual (int N_real, const double *Z, double *Resdl)
/**********************************************************************/
/* This function computes a particular nonlinear function of Z and    */
/* returns the result in the array Resdl. The function arises from    */
/* a pseudo-spectral discretisation of the one-dimensional BVP        */
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
/* The computed vector field is returned in the array of doubles      */
/* Resdl that has 2*N_{grid}+3 components. The array Resdl also uses  */
/* a compressed complex storage scheme (i.e., alternating real and    */
/* imaginary parts) to represent a complex vector using an array of   */
/* doubles. The last three components of are set to zero.             */
/**********************************************************************/
{
  int N_grid, k;
  complex c, nu;
  bool is_valid;

  /* N_real = number of REAL variables input    */
  /* N_grid = number of collocation points used */
  N_grid = (N_real / 2) - 2;
  c  = Z[N_real - 4] + 1.0I*Z[N_real-3];
  nu = Z[N_real - 2] + 1.0I*Z[N_real-1];
  /* Accumulate linear terms first */
  for (k = 0; k < N_grid; k++) {
    X_cplx[k] = Z[2*k] + 1.0I * Z[2*k + 1];
    RHS_cplx[k] = Dmatrix2[k] + nu * Dmatrix4[k];
    RHS_cplx[k] *= X_cplx[k];
    DX_cplx[k] = Dmatrix[k] * X_cplx[k];
    RHS_cplx[k] -= c * DX_cplx[k];
  }
  /* Now use DFTs to compute nonlinear terms of function. */
  /* Equivalent to X_cplx = ifft(X_cplx) in Matlab... */
  is_valid = fft_wrapper (false, N_grid, X_cplx);
  if (!is_valid) {
    fprintf (stderr, "%s: Error.\n", __func__);
	return (-1);
  }
  /* Equivalent to DX_cplx = ifft(DX_cplx) in Matlab... */
  is_valid = fft_wrapper (false, N_grid, DX_cplx);
  if (!is_valid) {
    fprintf (stderr, "%s: Error.\n", __func__);
	return (-1);
  }

  for (k = 0; k < N_grid; k++) {
    /* Equivalent to DX_cplx = ifft(DX_cplx) .* ifft(X_cplx) */
    DX_cplx[k] *= X_cplx[k];
    /* Equivalent to X_cplx = sin( ifft(X_cplx) ) */
    X_cplx[k] = csin (X_cplx[k]);
  }

  /* Equivalent to X_cplx = fft(X_cplx) in Matlab */
  is_valid = fft_wrapper (true, N_grid, X_cplx);
  if (!is_valid) {
	  fprintf (stderr, "%s: Error.\n", __func__);
	  return (-1);
  }
  /* Equivalent to DX_cplx = fft(DX_cplx) in Matlab */
  is_valid = fft_wrapper (true, N_grid, DX_cplx);
  if (!is_valid) {
	  fprintf (stderr, "%s: Error.\n", __func__);
	  return (-1);
  }
  /* Finally, accumulate nonlinear terms into residual vector & unpack
     real/imaginary components into vector Resdl of *doubles*. */
  for (k = 0; k < N_grid; k++) {
    RHS_cplx[k] += DX_cplx[k];
    RHS_cplx[k] += Aval * X_cplx[k];
    Resdl[2*k] = creal (RHS_cplx[k]);
    Resdl[2*k + 1] = cimag (RHS_cplx[k]);
  }
  Resdl[N_real-4] = 0.0;
  Resdl[N_real-3] = 0.0;
  Resdl[N_real-2] = 0.0;
  return 0;
}
/******************************************************************************/
