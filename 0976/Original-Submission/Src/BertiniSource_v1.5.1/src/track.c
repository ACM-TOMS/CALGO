// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

/************************************************************/
int track_d(point_data_d *Final, point_data_d *Start, comp_d FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/************************************************************/
/* Track the point on the given path, Start, to FinalT. The */
/* result is returned in Final.                             */
/* Return value:                                            */
/*    0 on success, nonzero otherwise.                      */
/* Assuming Start was a valid point on the path, then Final */
/* will also be a valid point on the path regardless of     */
/* the return value. It just may not have made it all the   */
/* way to the desired stop time.                            */
/************************************************************/
{
  double absValOfDistanceLeft, tempD, predictorError = 0;
  int cont, retVal = 0, consecSuccess, numSteps;
  comp_d nextTime, distanceLeft;
  point_data_d Current;
  eval_struct_d e;

  // initialize
  init_point_data_d(&Current, Start->point->size);
  init_eval_struct_d(e, 0, 0, 0);

  // copy Start to Current
  point_cp_d(Current.point, Start->point);
  set_d(Current.time, Start->time);

  numSteps = consecSuccess = 0;

  // if this is the first step of the path, initialize the current step size as the maximum step size - otherwise, it is already set!
  if (T->first_step_of_path)
  {
    T->currentStepSize = T->maxStepSize;
    T->first_step_of_path = 0;
  }

  // main tracking loop
  cont = 1;
  while (cont)
  { // print out the current time
    if (T->outputLevel > 0)
    {
      fprintf(OUT, "Time = %.15e\n", Current.time->r);
      if (T->screenOut)
        printf("Time = %.15e\n", Current.time->r);
    }

    // find the distance to the end = end time - current time
    sub_d(distanceLeft, FinalT, Current.time);
    absValOfDistanceLeft = d_abs_d(distanceLeft);

    // check to see how close we are to the target
    if (FinalT->r == Current.time->r && FinalT->i == Current.time->i)
    { // we are exactly at the end!
      retVal = cont = 0;
      break;
    }
    else if (absValOfDistanceLeft < T->currentStepSize)
    { // we can take a smaller step than our current step size to hit the target time
      set_d(nextTime, FinalT);
    }
    else
    { // the distance left is larger than the current step size, so we head in the correct direction will a step size of the current step size
      tempD = T->currentStepSize / absValOfDistanceLeft;
      mul_rdouble_d(nextTime, distanceLeft, tempD); // this will put deltaT = T->currentStepSize * sign(distanceleft)
      add_d(nextTime, Current.time, nextTime); // nextTime = current_time + deltaT
    }

    // take a step
    if (T->odePredictor == 3 || T->odePredictor == 4 || T->odePredictor == 5 || T->odePredictor == 6 || T->odePredictor == 7 || T->odePredictor == 8)
    { // use a ode predictor that has an error estimate with it
      retVal = step_error_d(&predictorError, consecSuccess, &Current, nextTime, T, OUT, &e, ED, eval_func);
    }
    else 
    { // use a standard ode predictor
      retVal = step_d(&Current, nextTime, T, OUT, &e, ED, eval_func);
    }
    numSteps++;

    // complete the output
    if (T->outputLevel > 1)
    {
      fprintf(OUT, "\n");
      if (T->screenOut)
        printf("\n");
    }

    // check to see if the step was successful
    if (retVal == 0)
    { // success!!

      // increment the number of consecutive successful steps
      consecSuccess++;
      // check to see if we should change the step size
      step_size_increase_d(predictorError, &consecSuccess, T);
    }
    else
    { // bad step - reset the number of consecutive successes
      consecSuccess = 0;
      // shrink the step size
      T->currentStepSize *= T->step_fail_factor;

      // do some checks
      if (T->currentStepSize < T->minStepSize)
      { // step size is too small, so we give up
        retVal = retVal_step_size_too_small;
        cont = 0;
        break;
      }
      else if (retVal == retVal_going_to_infinity)
      { // path is going to infinity
        retVal = retVal_going_to_infinity;
        cont = 0;
        break;
      }
    }

    // check to see if we have taken too many steps
    if (numSteps > T->maxNumSteps)
    {
      retVal = retVal_too_many_steps;
      cont = 0;
      break;
    }
  }

  // copy over the current information
  point_cp_d(Final->point, Current.point);
  set_d(Final->time, Current.time);

  // display how many steps it took
  if (T->outputLevel >= 0)
  {
    fprintf(OUT, "numSteps = %d\n", numSteps);
    if (T->screenOut)
      printf("numSteps = %d\n", numSteps);
  }

  // clear
  clear_point_data_d(&Current);
  clear_eval_struct_d(e);

  return retVal;
}
 
void step_size_increase_d(double predictorError, int *consecSuccess, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: adjust current step size                               *
\***************************************************************/
{
  // see if we can adjust the step size
  if (*consecSuccess >= T->cSecInc)
  { // determine which ode predictor we are using
    if (T->odePredictor >= 3)
    { // Runge-Kutta with error control
      double tempD = T->step_success_factor * T->currentStepSize;
      // the new step size is the minimum of the max step size and tempD
      T->currentStepSize = MIN(T->maxStepSize, tempD);
    }
    else // T->odePredictor < 3
    { // find the new step size if we increase it
      double tempD = T->step_success_factor * T->currentStepSize;
      // the new step size is the minimum of the max step size and tempD
      T->currentStepSize = MIN(T->maxStepSize, tempD);
    }

    // reset the number of consective successful steps
    *consecSuccess = 0;
  }

  return;
}

/************************************************************/
int step_d(point_data_d *P, comp_d newTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/

/************************************************************/
/* Use Euler and Newton to track the given point from       */
/* (P->time)  to  newTime                                   */
/************************************************************/
{ 
  int    retVal;
  vec_d  newParam, newPoint;
  comp_d dT;

  // initialize
  init_vec_d(newParam, 0);
  init_vec_d(newPoint, 0);

  // find dT = newTime - P->time
  sub_d(dT, newTime, P->time);

  if (T->odePredictor == -1)
  { // use constant method
    retVal = constant_d(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 1)
  { // use heun method
    retVal = heun_d(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 2)
  { // use standard Runge-Kutta predictor (RK4)
    retVal = runge_kutta_d(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else
  { // use euler predictor
    retVal = euler_d(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }

  if (retVal)
  { // clear
    clear_vec_d(newParam);
    clear_vec_d(newPoint);

    return retVal;
  }

  retVal = newton_d(newParam, newPoint, T, OUT, newTime, e, ED, eval_func);
  if (retVal)
  {
    if (retVal == retVal_going_to_infinity) // if the newPoint is at infinity, we want to copy it so that we can see the point so that we could adjust MAXNORM appropriately, if possible
    {
      point_cp_d(P->point, newPoint);
      set_d(P->time, newTime);
    }

    // clear
    clear_vec_d(newParam);
    clear_vec_d(newPoint);

    return retVal;
  }
  
  point_cp_d(P->point, newPoint);
//  point_cp_d(P->param, newParam);
  set_d(P->time, newTime);

  // clear
  clear_vec_d(newParam);
  clear_vec_d(newPoint);

  return 0;
}
    
int euler_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  // Get an estimate via Euler's method.
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
    }
  }

  // find Y = -dH/dt = -dH/dp * dp/dt
  rows = e->funcVals->size; // == e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_d(&e->funcVals->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_d(&e->funcVals->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_d(&e->funcVals->coord[i], &e->funcVals->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_d(e->parVals, e->Jv, e->funcVals, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(e->parVals, e->Jv, e->funcVals);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_d(e->parVals, e->Jv, e->funcVals);
  }  
  // matrixSolve calculated dX = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal) 
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in euler_d - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in euler_d - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in euler_d - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in euler_d - the step size will be decreased.\n");
    }
    return retVal;
  }

  // check to make sure that rules A & C are satisfied for AMP
  if (T->MPType == 2)
  {
    retVal = 0;
    // check criterion A from the AMP paper.
    if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
    {
      if (T->screenOut)
        fprintf(stderr, "Euler_d sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: Euler_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "Euler_d sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: Euler_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    if (retVal)
      return retVal;
  }

  // find newPoint = P->point + dX * dT
  rows = e->parVals->size;
  change_size_vec_d(newPoint, rows);
  newPoint->size = rows;
  for (i = 0; i < rows; i++)
  {
    mul_d(&newPoint->coord[i], dT, &e->parVals->coord[i]);
    add_d(&newPoint->coord[i], &newPoint->coord[i], &P->point->coord[i]);
  }

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("dX = ");  printVec_d(stdout, 10, e->parVals);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  return retVal;
}

int constant_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, retVal = 0;

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // evaluate for output
    eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
    }
  }

  // newPoint = P->point 
  vec_cp_d(newPoint,P->point);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // dX = 0
    change_size_vec_d(e->parVals, P->point->size);
    for (i = 0; i < P->point->size; i++)
      set_zero_d(&e->parVals->coord[i]);

    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("dX = ");  printVec_d(stdout, 10, e->parVals);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  return retVal;
}

int newton_d(vec_d newParam, vec_d newPoint, tracker_config_t *T, FILE *OUT, comp_d t1, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
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
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_d - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_d - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_d - the step size will be decreased.\n");
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
          fprintf(stderr, "Newton_d sees that AMP Criterion B has been violated!\n");
        fprintf(OUT, "NOTE: Newton_d sees that AMP Criterion B has been violated - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }

      // check criterion C
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: Newton_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
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

int newton_iteration_d(double *residual, int return_norm, double *norm_J, double *norm_J_inv, vec_d P, comp_d time, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find norm_j & norm_J_inv *
* NOTES: performs one newton iteration                          *
\***************************************************************/
{
  int retVal;
  double cond_num;

  // evaluate the function

  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED, eval_func);

  if (return_norm) 
  { // find the norm & norm_J_inv by doing matrixSolve & cond_num together
    retVal = matrixSolve_cond_num_norms_d(e->parVals, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_d(e->parVals, e->Jv, e->funcVals);
  }
  // matrixSolve calculated dX = (dH/dX)^(-1) * H(x,t).

  // check for matrixSolve failure
  if (!retVal)
  { // update P = P - dX
    vec_sub_d(P, P, e->parVals);
    // find the size of dX
    *residual = infNormVec_d(e->parVals);
  }

  return retVal;
}

int newton_residual_d(vec_d dX, double *residual, int return_norm, double *norm_J, double *norm_J_inv, vec_d P, comp_d time, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find norm_j & norm_J_inv *
* NOTES: finds the newton residual dX but does not update P     *
\***************************************************************/
{
  int retVal;
  double cond_num;

  // evaluate the function
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED, eval_func);

  if (return_norm)
  { // find the norm & norm_J_inv by doing matrixSolve & cond_num together
    retVal = matrixSolve_cond_num_norms_d(dX, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_d(dX, e->Jv, e->funcVals);
  }
  // matrixSolve calculated dX = (dH/dX)^(-1) * H(x,t).

  // check for matrixSolve failure
  *residual = 0;
  if (!retVal)
  { // find the size of dX
    *residual = infNormVec_d(dX);
  }

  return retVal;
}

int refine_d(point_data_d *P, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: newton_d makes sure that we have convergence or else P *
* is not updated by step_d. Refine_d updates P only if it       *
* converges.                                                    *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, maxIts = 20, num_digits = 16;
  double residual_old = 1e300, correction_tolerance = 1e-14;
  double size_of_newPoint, norm_J, norm_J_inv;

  point_data_d P_temp;
  eval_struct_d e;

  // initialize
  init_point_data_d(&P_temp, 0);
  init_eval_struct_d(e, 0, 0, 0);

  // copy P to P_temp
  point_data_cp_d(&P_temp, P);

  while (cont)  // We want to keep refining until we hit the desired residual or the residual levels off.
  { // do a newton iteration
    retVal = newton_iteration_d(&T->latest_newton_residual_d, T->MPType == 2, &norm_J, &norm_J_inv, P_temp.point, P_temp.time, &e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in refine_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in refine_d - the precision will be increased.\n");
 
        // clear
        clear_point_data_d(&P_temp);
        clear_eval_struct_d(e);

       return retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve failed in refine_d.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve failed in refine_d.");

        // clear
        clear_point_data_d(&P_temp);
        clear_eval_struct_d(e);

        return retVal;
      }
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "refining residual = %e\n", T->latest_newton_residual_d);
      if (T->screenOut)
        printf("refining residual = %e\n", T->latest_newton_residual_d);
    }

    // check to see if we need to do more iterations
    if (T->latest_newton_residual_d < correction_tolerance)
    { // we have converged
      cont = 0;
    }
    else if (0.95 * T->latest_newton_residual_d > residual_old)
    { // quit since residual is worse than previous one
      cont = 0;
    }
    else
    { // update residual_old
      residual_old = T->latest_newton_residual_d;

      // do AMP criterion checks if needed
      if (T->MPType == 2)
      {
        // check criterion B, if possible - notice that we are using maxIts instead of maxNewtonIts
        if ((maxIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, maxIts, its))
        {
          if (T->screenOut)
            fprintf(stderr, "NOTE: Refine_d sees that AMP Criterion B has been violated!\n");
          fprintf(OUT, "NOTE: Refine_d sees that AMP Criterion B has been violated.\n");

          // clear
          clear_point_data_d(&P_temp);
          clear_eval_struct_d(e);

          return retVal_higher_prec_needed;
        }

        // find the size of newPoint
        size_of_newPoint = infNormVec_d(P_temp.point);

        // check criterion C
        if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
        {
          if (T->screenOut)
            fprintf(stderr, "NOTE: Refine_d sees that AMP Criterion C has been violated!\n");
          fprintf(OUT, "NOTE: Refine_d sees that AMP Criterion C has been violated.\n");

          // clear
          clear_point_data_d(&P_temp);
          clear_eval_struct_d(e);

          return retVal_higher_prec_needed;
        }
      }

      its++; // increment the number of iterations completed
      // check to see if we have done too many iterations
      if (its > maxIts)
      {
        cont = 0;
      }
    }
  }

  // check for convergence - as long as it is better than either the correction tolerance (based on precision) or the final tolerance, we can call refine a success
  if (T->latest_newton_residual_d > correction_tolerance && T->latest_newton_residual_d > T->final_tolerance)
  {
    if (T->MPType == 2)
    { // AMP needs to use higher precision
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      retVal = retVal_refining_failed;
    }
  }
  else
  { // update P since converged
    point_data_cp_d(P, &P_temp);
    retVal = 0;
  }

  // clear
  clear_point_data_d(&P_temp);
  clear_eval_struct_d(e);

  return retVal;
}

/************************************************************/
int track_mp(point_data_mp *Final, point_data_mp *Start, comp_mp FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/************************************************************/
/* Track the point on the given path, Start, to FinalT. The */
/* result is returned in Final.                             */
/* Return value:                                            */
/*    0 on success, nonzero otherwise.                      */
/* Assuming Start was a valid point on the path, then Final */
/* will also be a valid point on the path regardless of     */
/* the return value. It just may not have made it all the   */
/* way to the desired stop time.                            */
/************************************************************/
{
  int cont, retVal = 0, consecSuccess, numSteps;

  mpf_t predictorError, currentStepSize, absValOfDistanceLeft, tempMPF;
  comp_mp nextTime, distanceLeft;
  point_data_mp Current;
  eval_struct_mp e;

  // initialize the MP data types
  mpf_init(predictorError); mpf_init(currentStepSize); mpf_init(absValOfDistanceLeft); mpf_init(tempMPF);
  init_mp(nextTime); init_mp(distanceLeft);
  init_point_data_mp(&Current, Start->point->size);
  init_eval_struct_mp(e, 0, 0, 0);

  // copy Start to Current
  point_cp_mp(Current.point, Start->point);
  set_mp(Current.time, Start->time);

  numSteps = consecSuccess = 0;

  // if this is the first step of the path, initialize the current step size as the maximum step size - otherwise, it is already set!
  if (T->first_step_of_path)
  {
    T->currentStepSize = T->maxStepSize;
    T->first_step_of_path = 0;
  }

  // start with the current step size
  mpf_set_d(currentStepSize, T->currentStepSize);

  // main tracking loop
  cont = 1;
  while (cont)
  { // print out the current time
    if (T->outputLevel > 0)
    {
      fprintf(OUT, "Time = "); print_mp(OUT, 0, Current.time); fprintf(OUT, "\n");
      if (T->screenOut)
      {
        fprintf(stdout, "Time = "); print_mp(stdout, 0, Current.time); fprintf(stdout, "\n");
      }
    }

    // find the distance to the end = end time - current time
    sub_mp(distanceLeft, FinalT, Current.time);
    mpf_mul(tempMPF, distanceLeft->r, distanceLeft->r);
    mpf_mul(absValOfDistanceLeft, distanceLeft->i, distanceLeft->i);
    mpf_add(absValOfDistanceLeft, absValOfDistanceLeft, tempMPF);
    mpf_sqrt(absValOfDistanceLeft, absValOfDistanceLeft);

    // check to see how close we are to the target
    if (mpfr_equal_p(FinalT->r, Current.time->r) && mpfr_equal_p(FinalT->i, Current.time->i))
    { // we are exactly at the end!
      retVal = cont = 0;
      break;
    }
    else if (mpfr_less_p(absValOfDistanceLeft, currentStepSize))
    { // we can take a smaller step than our current step size to hit the target time
      set_mp(nextTime, FinalT);
    }
    else
    { // the distance left is larger than the current step size, so we head in the correct direction will a step size of the current step size
      mpf_div(tempMPF, currentStepSize, absValOfDistanceLeft);
      // this will put deltaT = currentStepSize * sign(distanceLeft)
      mpf_mul(nextTime->r, distanceLeft->r, tempMPF);
      mpf_mul(nextTime->i, distanceLeft->i, tempMPF);
      // nextTime = current_time + deltaT
      add_mp(nextTime, nextTime, Current.time);
    }

    // take a step
    if (T->odePredictor == 3 || T->odePredictor == 4 || T->odePredictor == 5 || T->odePredictor == 6 || T->odePredictor == 7 || T->odePredictor == 8)
    { // use a ode predictor that has an error estimate with it
      retVal = step_error_mp(predictorError, consecSuccess, &Current, nextTime, T, OUT, &e, ED, eval_func);
    }
    else
    { // use a standard ode predictor
      retVal = step_mp(&Current, nextTime, T, OUT, &e, ED, eval_func);
    }
    numSteps++;

    // complete the output
    if (T->outputLevel > 1)
    {
      fprintf(OUT, "\n");
      if (T->screenOut)
        printf("\n");
    }

    // check to see if the step was successful
    if (retVal == 0)
    { // success!!
      consecSuccess++;

      // check to see if we should change the step size
      step_size_increase_mp(predictorError, currentStepSize, &consecSuccess, T);
    }
    else
    { // bad step - reset the number of consecutive successes
      consecSuccess = 0;
      // shrink the step size
      mpf_set_d(tempMPF, T->step_fail_factor);
      mpf_mul(currentStepSize, currentStepSize, tempMPF);

      // do some checks
      if (mpf_get_d(currentStepSize) < T->minStepSize)
      { // step size is too small, so we give up
        retVal = retVal_step_size_too_small;
        cont = 0;
        break;
      }
      else if (retVal == retVal_going_to_infinity)
      { // path is going to infinity
        retVal = retVal_going_to_infinity;
        cont = 0;
        break;
      }
    }

    // check to see if we have taken too many steps
    if (numSteps > T->maxNumSteps)
    {
      retVal = retVal_too_many_steps;
      cont = 0;
      break;
    }
  }

  // copy over the current information
  point_cp_mp(Final->point, Current.point);
  set_mp(Final->time, Current.time);
  T->currentStepSize = mpf_get_d(currentStepSize);

  // display how many steps it took
  if (T->outputLevel >= 0)
  {
    fprintf(OUT, "numSteps = %d\n", numSteps);
    if (T->screenOut)
      printf("numSteps = %d\n", numSteps);
  }

  // clear the MP data types
  mpf_clear(predictorError); mpf_clear(currentStepSize); mpf_clear(absValOfDistanceLeft); mpf_clear(tempMPF);
  clear_mp(nextTime); clear_mp(distanceLeft);
  clear_point_data_mp(&Current);
  clear_eval_struct_mp(e);

  return retVal;
}

void step_size_increase_mp(mpf_t predictorError, mpf_t currentStepSize, int *consecSuccess, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: adjust current step size                               *
\***************************************************************/
{
  mpf_t tempMPF;
  mpf_init(tempMPF);

  // see if we can adjust the step size
  if (*consecSuccess >= T->cSecInc)
  { // determine which ode predictor we are using
    if (T->odePredictor >= 3)
    { // Runge-Kutta with error control
      mpf_set_d(tempMPF, T->step_success_factor);
      mpf_mul(tempMPF, tempMPF, currentStepSize);

      // the new step size is the minimum of the max step size and tempD
      if (T->maxStepSize < mpf_get_d(tempMPF))
        mpf_set_d(currentStepSize, T->maxStepSize);
      else
        mpf_set(currentStepSize, tempMPF);
    }
    else // T->odePredictor < 3
    { // find the new step size if we increase it
      mpf_set_d(tempMPF, T->step_success_factor);
      mpf_mul(tempMPF, tempMPF, currentStepSize);

      // the new step size is the minimum of the max step size and tempD
      if (T->maxStepSize < mpf_get_d(tempMPF))
        mpf_set_d(currentStepSize, T->maxStepSize);
      else
        mpf_set(currentStepSize, tempMPF);
    }

    // reset the number of consective successful steps
    *consecSuccess = 0;
  }

  mpf_clear(tempMPF);

  return;
}
 
/************************************************************/
int step_mp(point_data_mp *P, comp_mp newTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/************************************************************/
/* Use Euler and Newton to track the given point from       */
/* (P->time)  to  newTime.                           */
/************************************************************/
{
  int retVal;

  vec_mp  newParam, newPoint;
  comp_mp dT;

  // initialize MP
  init_mp(dT);
  init_vec_mp(newParam, 0);
  init_vec_mp(newPoint, 0);

  // find dT = newTime - P->time
  sub_mp(dT, newTime, P->time);

  if (T->odePredictor == -1)
  { // use constant method
    retVal = constant_mp(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  } 
  else if (T->odePredictor == 1)
  { // use heun method
    retVal = heun_mp(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else if (T->odePredictor == 2)
  { // use standard Runge-Kutta predictor (RK4)
    retVal = runge_kutta_mp(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }
  else
  { // use euler predictor
    retVal = euler_mp(newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);
  }

  if (retVal)
  { // clear MP and return
    clear_mp(dT);
    clear_vec_mp(newParam);
    clear_vec_mp(newPoint);
    return retVal;
  }

  retVal = newton_mp(newParam, newPoint, P, T, OUT, newTime, e, ED, eval_func);
  if (retVal)
  {
    if (retVal == retVal_going_to_infinity) // if the newPoint is at infinity, we want to copy it so that we can see the point so that we could adjust MAXNORM appropriately, if possible
    {
      point_cp_mp(P->point, newPoint);
      set_mp(P->time, newTime);
    }
    // clear MP
    clear_mp(dT);
    clear_vec_mp(newParam);
    clear_vec_mp(newPoint);
    return retVal;
  }
  
  point_cp_mp(P->point, newPoint);
  //point_cp_mp(P->param, newParam);
  set_mp(P->time, newTime);

  // clear MP
  clear_mp(dT);
  clear_vec_mp(newParam);
  clear_vec_mp(newPoint);

  return 0;
}

int euler_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0;
  int num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;

  // Get an estimate using Euler's method. 
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
    }
  }

  // find Y = -dH/dt = -dH/dp * dp/dt
  rows = e->funcVals->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_mp(&e->funcVals->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_mp(&e->funcVals->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_mp(&e->funcVals->coord[i], &e->funcVals->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_mp(e->parVals, e->Jv, e->funcVals, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_mp(e->parVals, e->Jv, e->funcVals);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_mp(e->parVals, e->Jv, e->funcVals);
  }
  // matrixSolve calculated dX = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work     
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in euler_mp - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in euler_mp - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in euler_mp - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in euler_mp - the step size will be decreased.\n");
    }

    return retVal;
  }

  // check to make sure that rules A & C are satisfied for AMP
  if (T->MPType == 2)
  {
    retVal = 0;
    // check criterion A from the AMP paper.
    if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
    {
      if (T->screenOut)
        fprintf(stderr, "Euler_mp sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: Euler_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "Euler_mp sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: Euler_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    if (retVal)
      return retVal;
  }

  // find newPoint = P->point + dX * dT
  rows = e->parVals->size;
  change_size_point_mp(newPoint, rows);
  newPoint->size = rows;
  for (i = 0; i < rows; i++)
  {
    mul_mp(&newPoint->coord[i], dT, &e->parVals->coord[i]);
    add_mp(&newPoint->coord[i], &newPoint->coord[i], &P->point->coord[i]);
  }

  // Spit out the rest of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("dX = ");  printVec_mp(stdout, 0, e->parVals);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }
 
  return retVal; 
}  
  
int constant_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, retVal = 0;

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // evaluate for output
    eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
    }
  }

  // newPoint = P->point
  vec_cp_mp(newPoint, P->point);

  // Spit out the rest of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // dX == 0
    change_size_vec_mp(e->parVals, P->point->size);
    for (i = 0; i < P->point->size; i++)
      set_zero_mp(&e->parVals->coord[i]);

    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("dX = ");  printVec_mp(stdout, 0, e->parVals);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  return retVal;
}
 
int newton_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp t1, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int its = 0, retVal = 0, num_digits = (int) floor(T->Precision * log10(2.0));
  double size_of_newPoint, norm_J, norm_J_inv;

  while (its < T->maxNewtonIts)
  { // do a newton iteration
    retVal = newton_iteration_mp(T->latest_newton_residual_mp, T->MPType == 2, &norm_J, &norm_J_inv, newPoint, t1, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_mp - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in newton_mp - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in newton_mp - the step size will be decreased.\n");
        return retVal;
      }
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "residual = "); mpf_out_str(OUT, 10, 6, T->latest_newton_residual_mp); fprintf(OUT, "\n");
      if (T->screenOut)
      {
        printf("residual = "); mpf_out_str(stdout, 10, 6, T->latest_newton_residual_mp); printf("\n");
      }
    }

    // find double approximation
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);

    // check to see if the residual is small enough to call newton a success
    if (T->latest_newton_residual_d < T->currentNewtonTol)
    { // we have converged
      return 0;
    }

    // find the size of newPoint
    size_of_newPoint = infNormVec_mp(newPoint);

    // do AMP criterion checks if needed
    if (T->MPType == 2) 
    { // check criterion B, if possible
      if ((T->maxNewtonIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, T->maxNewtonIts, its))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_mp sees that AMP Criterion B has been violated!\n");
        fprintf(OUT, "NOTE: Newton_mp sees that AMP Criterion B has been violated - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }

      // check criterion C
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
      {
        if (T->screenOut)
          fprintf(stderr, "Newton_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: Newton_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        return retVal_higher_prec_needed;
      }
    }

    // check to see if the newPoint is too large
    if (T->endgameSwitch && size_of_newPoint > T->goingToInfinity)
      return retVal_going_to_infinity;

    its++;
  }

  // after doing the max number of newton iterations, check to see if we are near a solution and return the appropriate value
  if (T->latest_newton_residual_d > T->currentNewtonTol)
    return retVal_Failed_to_converge;
  else
    return 0;
}

int newton_iteration_mp(mpf_t residual, int return_norm, double *norm_J, double *norm_J_inv, vec_mp P, comp_mp time, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find norm_j & norm_J_inv *
* NOTES: performs one newton iteration                          *
\***************************************************************/
{
  int retVal;
  double cond_num;

  // evaluate the function
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED, eval_func);

  if (return_norm)
  { // find the norm & norm_J_inv by doing matrixSolve & cond_num together
    retVal = matrixSolve_cond_num_norms_mp(e->parVals, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_mp(e->parVals, e->Jv, e->funcVals);
  }
  // matrixSolve calculated dX = (dH/dX)^(-1) * H(x,t).

  // check for matrixSolve failure
  if (!retVal)
  { // matrixSolve success
    // update P = P - dX
    vec_sub_mp(P, P, e->parVals);
    // find the size of dX
    infNormVec_mp2(residual, e->parVals);
  }

  return retVal;
}

int newton_residual_mp(vec_mp dX, mpf_t residual, int return_norm, double *norm_J, double *norm_J_inv, vec_mp P, comp_mp time, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find norm_j & norm_J_inv *
* NOTES: finds the newton residual dX but does not update P     *
\***************************************************************/
{
  int retVal;
  double cond_num;

  // evaluate the function
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED, eval_func);

  if (return_norm)
  { // find the norm & norm_J_inv by doing matrixSolve & cond_num together
    retVal = matrixSolve_cond_num_norms_mp(dX, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_mp(dX, e->Jv, e->funcVals);
  }
  // matrixSolve calculated dX = (dH/dX)^(-1) * H(x,t).

  // check for matrixSolve failure
  if (!retVal)
  { // find the size of dX
    infNormVec_mp2(residual, dX);
  }

  return retVal;
}

int refine_mp(point_data_mp *P, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: newton_mp makes sure that we have convergence or else P*
* is not updated by step_mp. Refine_mp updates P only if it     *
* converges.                                                    *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, maxIts = 20;
  int num_digits = (int) floor(T->Precision * log10(2.0));
  int num_digits_for_correction = - (int) floor(T->Precision * log10(2.0) - 3.5);
  double residual_old = 1e308, correction_tolerance = num_digits_for_correction < -308 ? 1e-308 : pow(10, num_digits_for_correction);
  double size_of_newPoint, norm_J, norm_J_inv;

  point_data_mp P_temp;
  eval_struct_mp e;

  // initialize MP
  init_point_data_mp(&P_temp, 0);
  init_eval_struct_mp(e, 0, 0, 0);

  // copy P to P_temp
  point_data_cp_mp(&P_temp, P);

  while (cont)  //We want to keep refining until we hit the desired residual or the residual levels off.
  { // do a newton iteration
    retVal = newton_iteration_mp(T->latest_newton_residual_mp, T->MPType == 2, &norm_J, &norm_J_inv, P_temp.point, P_temp.time, &e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failured)
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in refine_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in refine_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in refine_mp.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in refine_mp.\n");
      }
      // clear MP
      clear_point_data_mp(&P_temp);
      clear_eval_struct_mp(e);

      // return error
      return retVal;
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "refining residual = "); mpf_out_str(OUT, 10, 7, T->latest_newton_residual_mp); fprintf(OUT, "\n");
      if (T->screenOut)
      {
        printf("refining residual = "); mpf_out_str(stdout, 10, 7, T->latest_newton_residual_mp); printf("\n");
      }
    }

    // find double approximation
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);

    // check to see if we need to do more iterations
    if (T->latest_newton_residual_d < correction_tolerance)
    { // we have converged
      cont = 0;
    }
    else if (0.95 * T->latest_newton_residual_d > residual_old)
    { // quit since residual is worse than previous one
      cont = 0;
    }
    else
    { // update residual_old
      residual_old = T->latest_newton_residual_d;

      // do AMP criterion checks if needed
      if (T->MPType == 2)
      {
        retVal = 0;
        // check criterion B, if possible - notice that we are using maxIts instead of maxNewtonIts
        if ((maxIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, residual_old, maxIts, its))
        {
          if (T->screenOut)
            fprintf(stderr, "NOTE: Refine_mp sees that AMP Criterion B has been violated!\n");
          fprintf(OUT, "NOTE: Refine_mp sees that AMP Criterion B has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }

        // find size_of_newPoint
        size_of_newPoint = infNormVec_mp(P_temp.point);

        // check criterion C
        if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, size_of_newPoint))
        {
          if (T->screenOut)
            fprintf(stderr, "NOTE: Refine_mp sees that AMP Criterion C has been violated!\n");
          fprintf(OUT, "NOTE: Refine_mp sees that AMP Criterion C has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }
        if (retVal)
        { // clear MP
          clear_point_data_mp(&P_temp);
          clear_eval_struct_mp(e);
          // return error
          return retVal;
        }
      }

      its++; // increment the number of iterations completed
      // check to see if we have done too many iterations
      if (its > maxIts)
      {
        cont = 0;
      }
    }
  }

  // check for convergence - as long as it is better than either the correction tolerance (based on precision) or the final tolerance, we can call refine a success
  if (T->latest_newton_residual_d > correction_tolerance && T->latest_newton_residual_d > T->final_tolerance)
  {
    if (T->MPType == 2)
    { // AMP needs to use higher precision
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      retVal = retVal_refining_failed;
    }
  }
  else
  { // update P since converged
    point_data_cp_mp(P, &P_temp);
    retVal = 0;
  }

  // clear MP
  clear_point_data_mp(&P_temp);
  clear_eval_struct_mp(e);

  return retVal;
}

