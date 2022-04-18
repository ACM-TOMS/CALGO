#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static int
 iterate(void *vstate, gsl_multifit_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx, int  scale)
{
  lmder_state_t *state = (lmder_state_t *) vstate;

  gsl_matrix *r = state->r;
  gsl_vector *tau = state->tau;
  gsl_vector *diag = state->diag;
  gsl_vector *qtf = state->qtf;
  gsl_vector *x_trial = state->x_trial;
  gsl_vector *f_trial = state->f_trial;
  gsl_vector *rptdx = state->rptdx;
  gsl_vector *newton = state->newton;
  gsl_vector *gradient = state->gradient;
  gsl_vector *sdiag = state->sdiag;
  gsl_vector *w = state->w;
  gsl_vector *work1 = state->work1;
  gsl_permutation *perm = state->perm;

  MpIeee prered;MpIeee  actred;
  MpIeee pnorm;MpIeee  fnorm1;MpIeee  fnorm1p;MpIeee  gnorm;
  MpIeee ratio;
  MpIeee dirder;

  int  iter=  0;

  MpIeee p1=  MpIeee( "0.1" );MpIeee  p25=  MpIeee( "0.25" );MpIeee  p5=  MpIeee( "0.5" );MpIeee  p75=  MpIeee( "0.75" );MpIeee  p0001=  MpIeee( "0.0001" );

  if (state->fnorm == 0.0) 
    {
      return GSL_SUCCESS;
    }

  /* Compute qtf = Q^T f */

  gsl_vector_memcpy (qtf, f);
  gsl_linalg_QR_QTvec (r, tau, qtf);

  /* Compute norm of scaled gradient */

  compute_gradient_direction (r, perm, qtf, diag, gradient);

  { 
    size_t iamax = gsl_blas_idamax (gradient);

    gnorm = fabs(gsl_vector_get (gradient, iamax) / state->fnorm);
  }

  /* Determine the Levenberg-Marquardt parameter */

lm_iteration:
  
  iter++ ;

  {
    int  status=  lmpar (r, perm, qtf, diag, state->delta, &(state->par), newton, gradient, sdiag, dx, w);
    if (status)
      return status;
  }

  /* Take a trial step */

  gsl_vector_scale (dx, -1.0); /* reverse the step to go downhill */

  compute_trial_step (x, dx, state->x_trial);

  pnorm = scaled_enorm (diag, dx);

  if (state->iter == 1)
    {
      if (pnorm < state->delta)
        {
#ifdef DEBUG
          {cout<<"set delta = pnorm = "<<setiosflags((ios::floatfield))<<  pnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif
          state->delta = pnorm;
        }
    }

  /* Evaluate function at x + p */
  /* return immediately if evaluation raised error */
  {
    int  status=  GSL_MULTIFIT_FN_EVAL_F (fdf, x_trial, f_trial);
    if (status)
      return status;
  }

  fnorm1 = enorm (f_trial);

  /* Compute the scaled actual reduction */

  actred = compute_actual_reduction (state->fnorm, fnorm1);

#ifdef DEBUG
  {cout<<"lmiterate: fnorm = "<<setiosflags((ios::floatfield))<< state->fnorm;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" fnorm1 = "<<setiosflags((ios::floatfield))<< fnorm1;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"  actred = "<<setiosflags((ios::floatfield))<< actred;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif

  /* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */

  compute_rptdx (r, perm, dx, rptdx);

  fnorm1p = enorm (rptdx);

  /* Compute the scaled predicted reduction = |J dx|^2 + 2 par |D dx|^2 */

  { 
    MpIeee t1=  fnorm1p / state->fnorm;
    MpIeee t2=  (sqrt(state->par) * pnorm) / state->fnorm;
    
    prered = t1 * t1 + t2 * t2 / p5;
    dirder = -(t1 * t1 + t2 * t2);
  }

  /* compute the ratio of the actual to predicted reduction */

  if (prered > MpIeee( "0" ))
    {
      ratio = actred / prered;
    }
  else
    {
      ratio = MpIeee( "0" );
    }

#ifdef DEBUG
  {cout<<"lmiterate: prered = "<<setiosflags((ios::floatfield))<< prered;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" dirder = "<<setiosflags((ios::floatfield))<< dirder;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" ratio = "<<setiosflags((ios::floatfield))<<ratio;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif


  /* update the step bound */

  if (ratio > p25)
    {
#ifdef DEBUG
      {cout<<"ratio > p25\n";}
#endif
      if (state->par == 0 || ratio >= p75)
        {
          state->delta = pnorm / p5;
          state->par *= p5;
#ifdef DEBUG
          {cout<<"updated step bounds: delta = "<<setiosflags((ios::floatfield))<< state->delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", par = "<<setiosflags((ios::floatfield))<< state->par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif
        }
    }
  else
    {
      MpIeee temp=  (actred >= MpIeee( "0" )) ? p5 : p5*dirder / (dirder + p5 * actred);

#ifdef DEBUG
      {cout<<"ratio < p25\n";}
#endif

      if (p1 * fnorm1 >= state->fnorm || temp < p1 ) 
        {
          temp = p1;
        }

      state->delta = temp * GSL_MIN_DBL (state->delta, pnorm/p1);

      state->par /= temp;
#ifdef DEBUG
      {cout<<"updated step bounds: delta = "<<setiosflags((ios::floatfield))<< state->delta;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", par = "<<setiosflags((ios::floatfield))<< state->par;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
#endif
    }


  /* test for successful iteration, termination and stringent tolerances */

  if (ratio >= p0001)
    {
      gsl_vector_memcpy (x, x_trial);
      gsl_vector_memcpy (f, f_trial);

      /* return immediately if evaluation raised error */
      {
        int  status=  GSL_MULTIFIT_FN_EVAL_DF (fdf, x_trial, J);
        if (status)
          return status;
      }

      /* wa2_j  = diag_j * x_j */
      state->xnorm = scaled_enorm(diag, x);
      state->fnorm = fnorm1;
      state->iter++;

      /* Rescale if necessary */

      if (scale)
        {
          update_diag (J, diag);
        }

      {
        int  signum;
        gsl_matrix_memcpy (r, J);
        gsl_linalg_QRPT_decomp (r, tau, perm, &signum, work1);
      }
      
      return GSL_SUCCESS;
    }
  else if (fabs(actred) <= GSL_DBL_EPSILON  && prered <= GSL_DBL_EPSILON 
           && p5 * ratio <= MpIeee( "1.0" ))
    {
      return GSL_ETOLF ;
    }
  else if (state->delta <= GSL_DBL_EPSILON * state->xnorm)
    {
      return GSL_ETOLX;
    }
  else if (gnorm <= GSL_DBL_EPSILON)
    {
      return GSL_ETOLG;
    }
  else if (iter < 10)
    {
      /* Repeat inner loop if unsuccessful */
      goto lm_iteration;
    }

  return GSL_CONTINUE;
}
