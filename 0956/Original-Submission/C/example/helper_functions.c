#include "KS_example.h"
/**********************************************************************/
bool setup_globals(int N_dim) {
  /* N_dim = number of REAL variables input */
  /* N_grid = number of collocation points  */
  int N_grid = (N_dim / 2) - 2;
  if (N_grid == 2048)
    Aval = -8.089150779042718; /* for N=2048 test case */
  else if (N_grid == 1024)
    Aval = -9.135914993072223; /* for N=1024 test case */
  else {
    fprintf (stderr, "%s: N_grid = %i\n", __func__, N_grid);
    fprintf (stderr, "%s: must have N_grid = 1024 or 2048\n", __func__);
    return false;
  }
  /* Setting up static arrays and data structures */
  Dmatrix = malloc (N_grid * sizeof (*Dmatrix));
  Dmatrix2 = malloc (N_grid * sizeof (*Dmatrix2));
  Dmatrix4 = malloc (N_grid * sizeof (*Dmatrix4));
  set_differentiation_matrices (N_grid);
  
  real_workspace = malloc (N_dim * sizeof (*real_workspace));
  X_cplx = malloc (N_grid * sizeof (*X_cplx));
  DX_cplx = malloc (N_grid * sizeof (*DX_cplx));
  Jac_cplx = malloc (pow(N_grid+2,2) * sizeof (*Jac_cplx));
  RHS_cplx = malloc ((N_grid+2) * sizeof (*RHS_cplx));
  ipiv = malloc ((N_grid+2) * sizeof(*ipiv));

  /* Structures here required for complex FFT in GSL */
  wavetable = gsl_fft_complex_wavetable_alloc (N_grid);
  workspace = gsl_fft_complex_workspace_alloc (N_grid);
  bool success = (Dmatrix!=NULL)   && (Dmatrix2!=NULL) &&
                 (Dmatrix4!=NULL)  && (DX_cplx!=NULL)  &&
                 (Jac_cplx!=NULL)  && (RHS_cplx!=NULL) &&
                 (RHS_cplx!=NULL)  && (ipiv!=NULL) &&
                 (wavetable!=NULL) && (workspace!=NULL);  
  if (!success) {
     printf ("%s: Memory allocation failure.\n", __func__);
     teardown_globals (); /* Release memory allocated successfully */
     return false;
  }
  return true;
}
/**********************************************************************/
void teardown_globals () {
  /* Cleaning up static arrays and data structures */
  if (Dmatrix != NULL)
    free (Dmatrix);
  if (Dmatrix2 != NULL)
    free (Dmatrix2);
  if (Dmatrix4 != NULL)
    free (Dmatrix4);
  if (real_workspace!=NULL)
    free (real_workspace);
  if (X_cplx != NULL)
    free (X_cplx);
  if (DX_cplx != NULL)
    free (DX_cplx);
  if (Jac_cplx != NULL)
    free (Jac_cplx);
  if (RHS_cplx != NULL)
    free (RHS_cplx);
  if (ipiv != NULL)
    free (ipiv);
  /* Data structures required for complex FFT in GSL */
  if (wavetable!=NULL)
    gsl_fft_complex_wavetable_free (wavetable);
  if (workspace!=NULL)
    gsl_fft_complex_workspace_free (workspace);
  return;
}
/**********************************************************************/
void
set_differentiation_matrices (int N_grid)
/**********************************************************************/
/* Helper routine for compute_residual and for single_corrector_step. */
/* Observe that differentiation matrices are diagonal and hence are   */
/* stored as 1D arrays. Also notice that the array Dmatrix is purely  */
/* imaginary, so Dmatrix2==(Dmatrix)^2 and Dmatrix4==(Dmatrix)^4 are  */
/* both real matrices (and hence are declared as doubles).            */
/* Note: it is assumed that memory allocation occurs elsewhere.       */
/**********************************************************************/
{
  int k;
  Dmatrix[0] = 0.0;
  Dmatrix2[0] = 0.0;
  Dmatrix4[0] = 0.0;
  for (k = 1; k < N_grid / 2; k++) {
    Dmatrix[k] = 1.0I*k;
    Dmatrix[N_grid - k] = -Dmatrix[k];
    Dmatrix2[k] = -pow (k, 2);
    Dmatrix2[N_grid - k] = -pow (k, 2);
    Dmatrix4[k] = pow (k, 4);
    Dmatrix4[N_grid - k] = pow (k, 4);
  }

  Dmatrix[N_grid / 2] = 0.0;
  Dmatrix2[N_grid / 2] = -pow (N_grid, 2) / 4.0;
  Dmatrix4[N_grid / 2] = pow (N_grid, 4) / 16.0;
  return;
}
/**********************************************************************/
bool compute_Jacobian (int N_dim, complex c, complex nu, double * T) {
  int N_grid, k, ell, index;
  bool success;
  /* N_dim = number of REAL variables input */
  /* N_grid = number of collocation points  */
  N_grid = (N_dim / 2) - 2;
  /* Resetting Jacobian entries to 0.0 */
  for (k=0; k<N_grid+2; k++)
    for (ell=0; ell<N_grid+2; ell++)
      Jac_cplx[k+(N_grid+2)*ell] = 0.0;

  for (k=0; k<N_grid; k++) {
    /* Setting diagonal elements */
    index = k + (N_grid+2) * k;
    Jac_cplx[index] = (-c) * Dmatrix[k];
    Jac_cplx[index] += Dmatrix2[k];
    Jac_cplx[index] += nu*Dmatrix4[k];
    /* Derivative wrt wavespeed c */
    index = k + (N_grid+2) * N_grid;
    Jac_cplx[index] = -Dmatrix[k] * X_cplx[k];
    /* Derivative wrt viscosity nu */
    index = k + (N_grid+2) * (N_grid+1);
    Jac_cplx[index] = Dmatrix4[k] * X_cplx[k];
    /* 2nd to last row = (D*X)' */
    index = N_grid + (N_grid+2) * k;
    Jac_cplx[index] = conj( Dmatrix[k] * X_cplx[k] );
    /* last row = T' */
    index = (N_grid+1) + (N_grid+2) * k;
    Jac_cplx[index] = T[2*k] - 1.I*T[2*k+1];
  }

  /* Complete last two entries of last row */
  index = (N_grid+1) + (N_grid+2) * N_grid;
  Jac_cplx[index] = T[N_dim-4]-1.0I*T[N_dim-3];
  index = (N_grid+1) + (N_grid+2) * (N_grid+1);
  Jac_cplx[index] = T[N_dim-2]-1.0I*T[N_dim-1];

  /* Accumulate advective terms in Jacobian matrix */
  for (k=0; k<=N_grid-1; k++)
    for (ell=0; ell<=k; ell++) {
      index = k + (N_grid+2) * ell;
      Jac_cplx[index] += Dmatrix[k-ell] * X_cplx[k-ell] / N_grid;
      Jac_cplx[index] += Dmatrix[ell] * X_cplx[k-ell] / N_grid;
    }

  for (k=0; k<=N_grid-2; k++)
    for (ell=k+1; ell<=N_grid-1; ell++) {
      index = k + (N_grid+2) * ell;
      Jac_cplx[index] += Dmatrix[k-ell+N_grid] *
                         X_cplx[k-ell+N_grid] / N_grid;
      Jac_cplx[index] += Dmatrix[ell] * X_cplx[k-ell+N_grid] / N_grid;
    }

  /* Use DFTs to compute nonlinear terms of function. */
  /* Equivalent to "X_cplx = Aval*cos(ifft(X_cplx))" in Matlab... */
  success = fft_wrapper (false, N_grid, X_cplx);
  if (!success)
    return false;
  for (k=0; k<N_grid; k++)
    X_cplx[k] = Aval*ccos( X_cplx[k] );

  /* Equivalent to X_cplx = fft(X_cplx) in Matlab... */
  success = fft_wrapper (true, N_grid, X_cplx);
  if (!success)
    return false;

  /* Accumulate nonlinear trigonometric derivative terms in Jacobian matrix */
  for (k=0; k<=N_grid-1; k++)
    for (ell=0; ell<=k; ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell] / N_grid;

  for (k=0; k<=N_grid-2; k++)
    for (ell=k+1; ell<=N_grid-1; ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell+N_grid] / N_grid;
  return true;
}
/**********************************************************************/
bool fft_wrapper (bool forward, int N, complex * Y) {
  int k, status;
  /* Unwrap array of complex values into array of doubles */
  for (k=0; k<N; k++) {
    real_workspace[2*k]   = creal (Y[k]);
    real_workspace[2*k+1] = cimag (Y[k]);
  }
  if (forward)
    status = gsl_fft_complex_forward (real_workspace, 1, N, 
                                      wavetable, workspace);
  else
    status = gsl_fft_complex_inverse (real_workspace, 1, N,
                                      wavetable, workspace);
  /* Check error status from call to GSL FFT routines */
  if (status) {
    if (status == GSL_EINVAL) {
       fprintf (stderr, "%s: invalid match, N=%d\n", 
                __func__, N);
    } else if (status == GSL_EDOM) {
       fprintf (stderr, "%s: invalid argument, N=%d\n", 
                __func__, N);
    } else {
       fprintf (stderr, "%s: Failed, gsl_errno=%d\n", 
                __func__, status);
    }
    return false;
  }
  /* Repack computed double values into complex array */
  for (k=0; k<N; k++)
    Y[k] = real_workspace[2*k] + 1.0I*real_workspace[2*k+1];
  return true;
}
