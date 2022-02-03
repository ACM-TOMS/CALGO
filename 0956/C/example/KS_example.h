#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <clapack.h>
#include "pampac.h"

extern bool setup_globals (int);
extern void teardown_globals ();
extern void set_differentiation_matrices (int);
extern bool compute_Jacobian (int, complex, complex, double*);
extern bool fft_wrapper (bool, int, complex* );

/* Note: These variables are all *globally* defined. Avoid name collision. */
double Aval, *Dmatrix2, *Dmatrix4, *real_workspace;
complex *Dmatrix, *X_cplx, *DX_cplx, *RHS_cplx, *Jac_cplx;
int *ipiv;
gsl_fft_complex_wavetable * wavetable;
gsl_fft_complex_workspace * workspace;

