// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

// Track using predictor methods with error estimates
int newton_error_d(double errorEstimate, vec_d newParam, vec_d newPoint, tracker_config_t *T, FILE *OUT, comp_d t1, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int newton_error_mp(mpf_t errorEstimate, vec_mp newParam, vec_mp newPoint, tracker_config_t *T, FILE *OUT, comp_mp t1, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int step_error_d(double *predictorError, int consecSuccess, point_data_d *P, comp_d newTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Take a step using error estimating predictor method    *
*  run newton if we are in the endgame or if we are about to    *
*  increase the step size                                       *
\***************************************************************/
{ // assume odePredictor == 3, 4, 5, 6, 7, or 8
  int    retVal = 0;
  vec_d  newParam, newPoint;
  comp_d dT;

  // initialize
  init_vec_d(newParam, 0);
  init_vec_d(newPoint, 0);

  // find dT = newTime - P->time
  sub_d(dT, newTime, P->time);

  // predictor step
  if (T->odePredictor == 3)
  { // use Euler-Heun method
    retVal = heun_error_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 4)
  { // use Runge-Kutta-Norsett (RKN34)
    retVal = rkn34_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 5)
  { // use Runge-Kutta-Fehlberg (RKF45)
    retVal = rkf45_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 6)
  { // use Runge-Kutta-Cash-Karp (RKCK45)
    retVal = rkck45_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 7)
  { // use Runge-Kutta-Dormand-Prince (RKDP56)
    retVal = rkdp56_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 8)
  { // use Runge-Kutta-Verner (RKV67)
    retVal = rkv67_d(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }

  // check for prediction success
  if (retVal)
  { // failure - clear memory and return
    clear_vec_d(newParam);
    clear_vec_d(newPoint);

    return retVal;
  }

  // corrector step - use if endgame or if predictor error is larger than tolerance or if we are about to increase the step size
  if (T->endgameSwitch || *predictorError >= T->currentNewtonTol || consecSuccess >= T->cSecInc - 1)
  { // run newton iterations
    retVal = newton_error_d(*predictorError, newParam, newPoint, T, OUT, newTime, e, ED, eval_func);

    // check for newton success
    if (retVal)
    { // if the newPoint is at infinity, we want to copy it so that we can see the point so that we could adjust MAXNORM appropriately, if possible
      if (retVal == retVal_going_to_infinity) 
      {
        point_cp_d(P->point, newPoint);
        set_d(P->time, newTime);
      }

      // clear
      clear_vec_d(newParam);
      clear_vec_d(newPoint);

      return retVal;
    }
  }
  else
  { // setup residual
    T->latest_newton_residual_d = *predictorError;
  }

  // successful step - copy it over
  point_cp_d(P->point, newPoint);
  set_d(P->time, newTime);

  // clear
  clear_vec_d(newParam);
  clear_vec_d(newPoint);

  return 0;
}

int newton_error_d(double errorEstimate, vec_d newParam, vec_d newPoint, tracker_config_t *T, FILE *OUT, comp_d t1, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int its = 0, retVal = 0, num_digits = 16;
  double size_of_newPoint, norm_J, norm_J_inv;

  while (its < T->maxNewtonIts)
  { // do a newton iteration
    retVal = newton_iteration_d(&T->latest_newton_residual_d, T->MPType == 2, &norm_J, &norm_J_inv, newPoint, t1, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_error_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_error_d - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_error_d - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_error_d - the step size will be decreased.\n");
        return retVal;
      }
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "residual = %e\n", T->latest_newton_residual_d);
      if (T->screenOut)
        printf("residual = %e\n", T->latest_newton_residual_d);
    }

    // check to see if the residual is small enough to call newton a success
    if (T->latest_newton_residual_d < T->currentNewtonTol)
      return 0;

    // find the size of newPoint
    size_of_newPoint = infNormVec_d(newPoint);

    // do AMP criterion checks if needed
    if (T->MPType == 2)
    { // check criterion B, if possible
      if ((T->maxNewtonIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, T->maxNewtonIts, its))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_error_d sees that AMP Criterion B has been violated!\n");
        fprintf(OUT, "NOTE: Newton_error_d sees that AMP Criterion B has been violated - the precision will be increased.\n");

        return retVal_higher_prec_needed;
      }

      // check criterion C
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_error_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: Newton_error_d sees that AMP Criterion C has been violated - the precision will be increased.\n");

        return retVal_higher_prec_needed;
      }
    }

    // check to see if the newPoint is too large
    if ((T->endgameSwitch) && (size_of_newPoint > T->goingToInfinity))
      return retVal_going_to_infinity;

    its++;
  }

  // after doing the max number of newton iterations, check to see if we are near a solution and return the appropriate value
  if (T->latest_newton_residual_d > T->currentNewtonTol)
    return retVal_Failed_to_converge;
  else
    return 0;
}

// MULTI PRECISION

int step_error_mp(mpf_t predictorError, int consecSuccess, point_data_mp *P, comp_mp newTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Take a step using error estimating predictor method    *
*  run newton if we are in the endgame or if we are about to    *
*  increase the step size                                       *
\***************************************************************/
{ // assume odePredictor == 3, 4, 5, 6, 7, or 8
  int     retVal = 0;
  comp_mp dT;
  vec_mp  newParam, newPoint;

  // initialize
  init_mp(dT);
  init_vec_mp(newParam, 0);
  init_vec_mp(newPoint, 0);

  // find dT = newTime - P->time
  sub_mp(dT, newTime, P->time);

  // predictor step
  if (T->odePredictor == 3)
  { // use Euler-Heun method
    retVal = heun_error_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 4)
  { // use Runge-Kutta-Norsett (RKN34)
    retVal = rkn34_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 5)
  { // use Runge-Kutta-Fehlberg (RKF45)
    retVal = rkf45_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 6)
  { // use Runge-Kutta-Cash-Karp (RKCK45)
    retVal = rkck45_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 7)
  { // use Runge-Kutta-Dormand-Prince (RKDP56)
    retVal = rkdp56_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 8)
  { // use Runge-Kutta-Verner (RKV67)
    retVal = rkv67_mp(predictorError, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }

  // check for prediction success
  if (retVal)
  { // failure - clear memory and return
    clear_mp(dT);
    clear_vec_mp(newParam);
    clear_vec_mp(newPoint);

    return retVal;
  }

  // corrector step - use if endgame or if predictor error is larger than tolerance or if we are about to increase the step size
  if (T->endgameSwitch || mpf_get_d(predictorError) >= T->currentNewtonTol || consecSuccess >= T->cSecInc - 1)
  { // run newton iterations
    retVal = newton_error_mp(predictorError, newParam, newPoint, T, OUT, newTime, e, ED, eval_func);

    // check for newton success
    if (retVal)
    { // if the newPoint is at infinity, we want to copy it so that we can see the point so that we could adjust MAXNORM appropriately, if possible
      if (retVal == retVal_going_to_infinity)
      {
        point_cp_mp(P->point, newPoint);
        set_mp(P->time, newTime);
      }

      // clear
      clear_vec_mp(newParam);
      clear_vec_mp(newPoint);

      return retVal;
    }
  }
  else
  { // setup residual
    mpf_set(T->latest_newton_residual_mp, predictorError);
  }

  // successful step - copy it over
  point_cp_mp(P->point, newPoint);
  set_mp(P->time, newTime);

  // clear
  clear_mp(dT);
  clear_vec_mp(newParam);
  clear_vec_mp(newPoint);

  return 0;
}

int newton_error_mp(mpf_t errorEstimate, vec_mp newParam, vec_mp newPoint, tracker_config_t *T, FILE *OUT, comp_mp t1, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int its = 0, retVal = 0, num_digits = prec_to_digits(T->Precision);
  double size_of_newPoint, norm_J, norm_J_inv;

  while (its < T->maxNewtonIts)
  { // do a newton iteration
    retVal = newton_iteration_mp(T->latest_newton_residual_mp, T->MPType == 2, &norm_J, &norm_J_inv, newPoint, t1, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_error_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_error_mp - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_error_mp - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_error_mp - the step size will be decreased.\n");
        return retVal;
      }
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "residual = %e\n", mpf_get_d(T->latest_newton_residual_mp));
      if (T->screenOut)
        printf("residual = %e\n", mpf_get_d(T->latest_newton_residual_mp));
    }

    // find double approximation
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);

    // check to see if the residual is small enough to call newton a success
    if (T->latest_newton_residual_d < T->currentNewtonTol)
      return 0;

    // find the size of newPoint
    size_of_newPoint = infNormVec_mp(newPoint);

    // do AMP criterion checks if needed
    if (T->MPType == 2)
    { // check criterion B, if possible
      if ((T->maxNewtonIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, T->maxNewtonIts, its))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_error_mp sees that AMP Criterion B has been violated!\n");
        fprintf(OUT, "NOTE: Newton_error_mp sees that AMP Criterion B has been violated - the precision will be increased.\n");

        return retVal_higher_prec_needed;
      }

      // check criterion C
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_error_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: Newton_error_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");

        return retVal_higher_prec_needed;
      }
    }

    // check to see if the newPoint is too large
    if ((T->endgameSwitch) && (size_of_newPoint > T->goingToInfinity))
      return retVal_going_to_infinity;

    its++;
  }

  // after doing the max number of newton iterations, check to see if we are near a solution and return the appropriate value
  if (T->latest_newton_residual_d > T->currentNewtonTol)
    return retVal_Failed_to_converge;
  else
    return 0;
}


