// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// this file contains different ODE methods

int heun_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Heun's method - Kincaid/Cheney pg 540-541         *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d  newTime;
  point_d Y, dX1, dX2;

  init_point_d(Y, 0); init_point_d(dX1, 0); init_point_d(dX2, 0);

  // find f(x,t) = -(Jv^-1)*(Jp*dp/dt) = dX

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_d(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_d(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_d(&Y->coord[i], &Y->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_d(dX1, e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX1, e->Jv, Y);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_d(dX1, e->Jv, Y);
  }
  // matrixSolve calculated dX1 = f(x,t) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_d - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_d - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

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
        fprintf(stderr, "heun_d sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: heun_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "heun_d sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: heun_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    if (retVal)
    {
      clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

      return retVal;
    }
  }

  // find Y = x + dT * dX1
  increase_size_vec_d(Y, dX1->size);
  rows = Y->size = dX1->size;
  for (i = 0; i < rows; i++)
  {
    mul_d(&Y->coord[i], dT, &dX1->coord[i]);
    add_d(&Y->coord[i], &Y->coord[i], &P->point->coord[i]);
  }

  // find newTime = t + dT
  add_d(newTime, dT, P->time);

  // do the second function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_d(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_d(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_d(&Y->coord[i], &Y->coord[i]);
  }

  // just do matrixSolve
  retVal = matrixSolve_d(dX2, e->Jv, Y);

  // matrixSolve calculated dX2 = f(x + dT*dX1,t + dT) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_d - the precision will be increased.\n");
      if (T->screenOut)         
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_d - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

    return retVal;
  }

  // find newPoint = P->point + 1/2*dT*(dX1 + dX2)
  increase_size_vec_d(newPoint, dX2->size);
  mul_rdouble_d(newTime, dT, 0.5);
  rows = newPoint->size = dX2->size;
  for (i = 0; i < rows; i++)
  {
    add_d(&dX1->coord[i], &dX1->coord[i], &dX2->coord[i]);
    mul_d(&dX1->coord[i], &dX1->coord[i], newTime);
    add_d(&newPoint->coord[i], &dX1->coord[i], &P->point->coord[i]);
  }

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX1);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX1);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

  return retVal;
}

int heun_error_d(double *error, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Heun's method with error approximation            *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d newTime;
  point_d Y, dX1, dX2;

  init_point_d(Y, 0); init_point_d(dX1, 0); init_point_d(dX2, 0);

  // find f(x,t) = -(Jv^-1)*(Jp*dp/dt) = dX

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_d(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_d(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_d(&Y->coord[i], &Y->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_d(dX1, e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX1, e->Jv, Y);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_d(dX1, e->Jv, Y);
  }
  // matrixSolve calculated dX1 = f(x,t) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_error_d - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_error_d - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

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
        fprintf(stderr, "heun_error_d sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: heun_error_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "heun_error_d sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: heun_error_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    if (retVal)
    {
      clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

      return retVal;
    }
  }

  // setup dX1 = dT * dX1
  // find Y = x + dX1
  increase_size_vec_d(Y, dX1->size);
  rows = Y->size = dX1->size;
  for (i = 0; i < rows; i++)
  {
    mul_d(&dX1->coord[i], dT, &dX1->coord[i]);
    add_d(&Y->coord[i], &P->point->coord[i], &dX1->coord[i]);
  }

  // find newTime = t + dT
  add_d(newTime, dT, P->time);

  // do the second function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_d(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_d(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_d(&Y->coord[i], &Y->coord[i]);
  }

  // just do matrixSolve
  retVal = matrixSolve_d(dX2, e->Jv, Y);

  // matrixSolve calculated dX2 = f(x + dT*dX1,t + dT) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_error_d - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_error_d - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

    return retVal;
  }

  // setup dX2 = dT*dX2
  // compute error: 1/2*(dX1 - dX2)
  // find newPoint = P->point + 1/2*(dX1 + dX2)
  increase_size_vec_d(Y, dX2->size);
  increase_size_vec_d(newPoint, dX2->size);
  rows = newPoint->size = Y->size = dX2->size;
  for (i = 0; i < rows; i++)
  {
    mul_d(&dX2->coord[i], &dX2->coord[i], dT);

    sub_d(&Y->coord[i], &dX1->coord[i], &dX2->coord[i]);
    add_d(&dX1->coord[i], &dX1->coord[i], &dX2->coord[i]);

    mul_rdouble_d(&Y->coord[i], &Y->coord[i], 0.5);
    mul_rdouble_d(&dX1->coord[i], &dX1->coord[i], 0.5);

    add_d(&newPoint->coord[i], &dX1->coord[i], &P->point->coord[i]);
  }
  *error = infNormVec_d(Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX1);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX1);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y); clear_point_d(dX1); clear_point_d(dX2);

  return retVal;
}

int runge_kutta_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the Runge Kutta method - Kincaid/Cheney pg 541    *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;
 
  comp_d  newTime;
  point_d Y, dX[4];

  init_point_d(Y, 0);
  for (i = 0; i < 4; i++)
    init_point_d(dX[i], 0);
 
  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 4; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }
 
    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_d(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y
 
    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in runge_kutta_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in runge_kutta_d - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      clear_point_d(Y);
      for (i = 0; i < 4; i++)
        clear_point_d(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "runge_kutta_d sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: runge_kutta_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "runge_kutta_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: runge_kutta_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      {
        clear_point_d(Y);
        for (i = 0; i < 4; i++)
          clear_point_d(dX[i]);

        return retVal;
      }
    }
 
    // calculate dX[k] = dT * dX[k]
    vec_mulcomp_d(dX[k], dX[k], dT);

    if (k < 3)
    { // another function evaluation is needed

      // find Y = x + a * dX[k], where a = 1/2 when k = 0 || k == 1 and a = 1 when k == 2
      increase_size_point_d(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      if (k == 0 || k == 1)
      { // find Y
        for (i = 0; i < rows; i++)
        {
          mul_rdouble_d(&Y->coord[i], &dX[k]->coord[i], 1.0 / 2);
          add_d(&Y->coord[i], &Y->coord[i], &P->point->coord[i]);
        }

        // find newTime = t + a * dT
        mul_rdouble_d(newTime, dT, 1.0 / 2);
        add_d(newTime, newTime, P->time);
      }
      else // k == 2
      { // find Y
        for (i = 0; i < rows; i++)
        {
          add_d(&Y->coord[i], &dX[k]->coord[i], &P->point->coord[i]);
        }
      
        // find newTime = t + dT
        add_d(newTime, dT, P->time);
      }

      // do the next function evaluation
      eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
  }

  // find newPoint = P->point + 1/6*(dX[0] + dX[3] + 2*(dX[1] + dX[2]))
  change_size_vec_d(newPoint, dX[0]->size);
  rows = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    add_d(&dX[1]->coord[i], &dX[1]->coord[i], &dX[2]->coord[i]);
    mul_rdouble_d(&dX[1]->coord[i], &dX[1]->coord[i], 2.0);
    add_d(&dX[0]->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
    add_d(&dX[0]->coord[i], &dX[0]->coord[i], &dX[3]->coord[i]);
    mul_rdouble_d(&dX[0]->coord[i], &dX[0]->coord[i], 1.0 / 6);
    add_d(&newPoint->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
  }

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX[0]);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y);
  for (i = 0; i < 4; i++)
    clear_point_d(dX[i]);
 
  return retVal;
}

int rkn34_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKN34 method - Norsett method                 *
\***************************************************************/
{
  int rV;
  double a[5] = {25.0/162, 32.0/135, 256.0/567, 0, 11.0/70};
  double a_minus_b[5] = {-41.0/4050, -244.0/1755, 256.0/567, -448.0/975, 11.0/70};
  double c[5] = {0, 3.0/8, 9.0/16, 25.0/32, 1};
  double d[5][4] = {{0, 0, 0, 0},
                    {3.0/8, 0, 0, 0},
                    {0, 9.0/16, 0, 0},
                    {-125.0/672, 325.0/336, 0, 0},
                    {371.0/891, -200.0/297, 1120.0/891, 0}};

  // run the Runge-Kutta34 method
  rV = rk34_methods_d(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rk34_methods_d(double *err, double a[5], double a_minus_b[5], double c[5], double d[5][4], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK34 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d  newTime;
  point_d Y, dX[5];

  init_point_d(Y, 0);
  for (i = 0; i < 5; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 5; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_d(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk34_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk34_d - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      clear_point_d(Y);
      for (i = 0; i < 5; i++)
        clear_point_d(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk34_d sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk34_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk34_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk34_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      {
        clear_point_d(Y);
        for (i = 0; i < 5; i++)
          clear_point_d(dX[i]);

        return retVal;
      }
    }

    if (k < 4)
    { // another function evaluation is needed
      // compute dX[k] = dX[k] * dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_d(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_d(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_d(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rdouble_d(newTime, dT, c[k+1]);
      add_d(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_d(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_d(Y, dX[0]->size);
  increase_size_vec_d(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rdouble_d(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 5; k++)
    {
      mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_d(&Y->coord[i], &Y->coord[i], newTime);

      mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
      add_d(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_d(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  *err = infNormVec_d(Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX[0]);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y);
  for (i = 0; i < 5; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rkf45_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKF45 method - Kincaid/Cheney pg 544-545      *
\***************************************************************/
{
  int rV;
  double a[6] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
  double a_minus_b[6] = {1.0/360, 0, -128.0/4275, -2197.0/75240, 1.0/50, 2.0/55};
  double c[6] = {0, 1.0/4, 3.0/8, 12.0/13, 1, 1.0/2};
  double d[6][5] = {{0, 0, 0, 0, 0}, 
                    {1.0/4, 0, 0, 0, 0}, 
                    {3.0/32, 9.0/32, 0, 0, 0}, 
                    {1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0}, 
                    {439.0/216, -8, 3680.0/513, -845.0/4104, 0}, 
                    {-8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40}};

  // run the Runge-Kutta45 method
  rV = rk45_methods_d(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkck45_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKCK45 method - Cash/Karp method              *
\***************************************************************/
{
  int rV;
  double a[6] = {37.0/378, 0, 250.0/621, 125.0/594, 0, 512.0/1771};
  double a_minus_b[6] = {-277.0/64512, 0, 6925.0/370944, -6925.0/202752, -277.0/14336, 277.0/7084};
  double c[6] = {0, 1.0/5, 3.0/10, 3.0/5, 1, 7.0/8};
  double d[6][5] = {{0, 0, 0, 0, 0},
                    {1.0/5, 0, 0, 0, 0},
                    {3.0/40, 9.0/40, 0, 0, 0},
                    {3.0/10, -9.0/10, 6.0/5, 0, 0},
                    {-11.0/54, 5.0/2, -70.0/27, 35.0/27, 0},
                    {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096}};

  // run the Runge-Kutta45 method
  rV = rk45_methods_d(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rk45_methods_d(double *err, double a[6], double a_minus_b[6], double c[6], double d[6][5], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d  newTime;
  point_d Y, dX[6];

  init_point_d(Y, 0);
  for (i = 0; i < 6; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 6; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_d(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk45_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk45_d - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      clear_point_d(Y);
      for (i = 0; i < 6; i++)
        clear_point_d(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk45_d sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk45_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
  
      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk45_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk45_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      {
        clear_point_d(Y);
        for (i = 0; i < 6; i++)
          clear_point_d(dX[i]);

        return retVal;
      }
    }

    if (k < 5)
    { // another function evaluation is needed
      // compute dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_d(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_d(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_d(&Y->coord[i], &Y->coord[i], newTime);
        }
      }
 
      // find newTime = t + c[k+1]*dT
      mul_rdouble_d(newTime, dT, c[k+1]);
      add_d(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_d(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_d(Y, dX[0]->size);
  increase_size_vec_d(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rdouble_d(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 6; k++)
    {
      mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_d(&Y->coord[i], &Y->coord[i], newTime);

      mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
      add_d(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_d(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  *err = infNormVec_d(Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX[0]);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y);
  for (i = 0; i < 6; i++)
    clear_point_d(dX[i]);

  return retVal;
} 

int rkdp56_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKDP56 method - Dormand/Prince method         *
\***************************************************************/
{
  int rV;
  double a[8] = {821.0/10800, 0, 19683.0/71825, 175273.0/912600, 395.0/3672, 785.0/2704, 3.0/50, 0};
  double a_minus_b[8] = {13.0/2400, 0, -19683.0/618800, 2401.0/31200, -65.0/816, 15.0/416, 521.0/5600, -1.0/10};
  double c[8] = {0, 1.0/10, 2.0/9, 3.0/7, 3.0/5, 4.0/5, 1, 1};
  double d[8][7] = {{0, 0, 0, 0, 0, 0, 0},
                    {1.0/10, 0, 0, 0, 0, 0, 0},
                    {-2.0/81, 20.0/81, 0, 0, 0, 0, 0},
                    {615.0/1372, -270.0/343, 1053.0/1372, 0, 0, 0, 0},
                    {3243.0/5500, -54.0/55, 50949.0/71500, 4998.0/17875, 0, 0, 0},
                    {-26492.0/37125, 72.0/55, 2808.0/23375, -24206.0/37125, 338.0/459, 0, 0},
                    {5561.0/2376, -35.0/11, -24117.0/31603, 899983.0/200772, -5225.0/1836, 3925.0/4056, 0},
                    {465467.0/266112, -2945.0/1232, -5610201.0/14158144, 10513573.0/3212352, -424325.0/205632, 376225.0/454272, 0}};

  // run the Runge-Kutta56 method
  rV = rk56_methods_d(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rk56_methods_d(double *err, double a[8], double a_minus_b[8], double c[8], double d[8][7], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK56 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d  newTime;
  point_d Y, dX[8];

  init_point_d(Y, 0);
  for (i = 0; i < 8; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 8; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_d(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk56_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk56_d - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      clear_point_d(Y);
      for (i = 0; i < 8; i++)
        clear_point_d(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk56_d sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk56_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk56_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk56_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      {
        clear_point_d(Y);
        for (i = 0; i < 8; i++)
          clear_point_d(dX[i]);

        return retVal;
      }
    }

    if (k < 7)
    { // another function evaluation is needed
      // calculate dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_d(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_d(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_d(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rdouble_d(newTime, dT, c[k+1]);
      add_d(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_d(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..7, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..7, a_k * dX[k])
  increase_size_vec_d(Y, dX[0]->size);
  increase_size_vec_d(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rdouble_d(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 8; k++)
    {
      mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_d(&Y->coord[i], &Y->coord[i], newTime);

      mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
      add_d(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_d(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  *err = infNormVec_d(Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX[0]);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y);
  for (i = 0; i < 8; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rkv67_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKV67 method - Verner method                  *
\***************************************************************/
{
  int rV;
  double a[10], a_minus_b[10], c[10], d[10][9];

a[0] = 0.0471556184862722217043176510883817567956932673866498263717863225066525727846555841973346401995827631062036755590558062534375047415234260388676624358999573840784918991178510655224244139232889913788333823544199823935451154623648793450503723002665841976598046959664007442181911565685230727013689593392982;
a[1] = 0;
a[2] = 0;
a[3] = 0.257505642984341518959643610103768758098639571258883586561584933570028952897166612637862059798925621609396907594245725006327476849765387793842037461290921370430264152893947313918692476051185173305706391851963836018457313459645894963869257643822167873501457534217624142946250030154297352910582419606143;
a[4] = 0.262166539774126204771386309576452771112914046736971319203584006643876812614319887681344771902494347386463073904093142446268038453678784950161124772440805335094258253338958758076230419920570973073193926469534220996564683552221244076754949769609944998538014278529773548771150093620971914747932175414550;
a[5] = 0.152160926567385574032313319916511753552259518225131674367096139045695165543077407261967511436494354991858803879708884258567824761923036313221569840599729116798542290424111444095404797846598555249340492159421091497761496624438637702252254645961603915777481188523759942071760307205888241753387037725744;
a[6] = 0.493996917003248424690717589322787684429585308571357948791678035642674335879997785662081156612226972186240516630330279108098440721497194764675468796945045160954090374742547119955506272386697362858820688622832619000271491070870242751247344166824766097061388163627211199461983342148211784138952128223367;
a[7] = -0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817;
a[8] = 0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153;
a[9] = 0;

a_minus_b[0] = 0.00254701187993104541699947511358977898137393503651067684217948288343659637043927033143148723131125158014121649938900878979858451840352977282793140381518240008963362033067932018366864275478963166454694523116781829305324081676250402401051739684741551799350892643994927331531661450056821725901262420206562;
a_minus_b[1] = 0;
a_minus_b[2] = 0;
a_minus_b[3] = -0.00965839487279574909126661599061503187520108884519528045821452710507663251123547574163069610539919389963989845238921904972854937606654792208299092082065749981679263365462373361333668366236519725184736326271819113885999268066494138424342590942206813226759750600367884910177105022385487684166900829606156;
a_minus_b[4] = 0.0420647097563969027734147319113774614805948824559621257010925807491537576600651305254295976225991591885211868879798925804616780192287010286942291836640369052613351689178787207331605464157207419903756767065167872280544436640338531219590767268952538968560279938316779313158803412904221702386669877100151;
a_minus_b[5] = -0.0666822437469301090659987634347776289054953577328478620986708910886522673765780883582819579617699092650805636589852395875931279767916652676612036293647748766874420794927278867725689176383726462247671413895793673537650626399046027049371138790153665954418043148360367277640961073763921508851946773219422;
a_minus_b[6] = 0.265009746462128136352900200346432447893361385631363781535494557795019991145447029415600259733188324877703429783214090406452549289091431065680284500120422909432522977485881928123107977267794033274265297233044517672944951473291445755463803660526233224001811988404937384560301642803203641418032914674657;
a_minus_b[7] = -0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817;
a_minus_b[8] = 0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153;
a_minus_b[9] = -0.0202951846633562822276705479381043035855420443667990862261517658249536055689205887319585505702055732014823934917746960666918489454776188166903872302377514709236100830696726470857731850092255075866785330602598147948274804639773599730786794693464008286038012269720794348562965112960546349366261206596328;

c[0] = 0;
c[1] = 0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[2] = 0.108888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889;
c[3] = 0.163333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
c[4] = 0.455500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[5] = 0.609509448997838131708700442148602494963796759286605673642939302386787496187674982310265299374033791406090583319455261759148494730544940521238957202640385592141084318199590688285015525265373754263537303527537491030484005373807797476632985866792541212386248529601147679387596870467556824874147722141014;
c[6] = 0.884000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[7] = 0.925000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[8] = 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[9] = 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;

d[0][0] = 0;
d[0][1] = 0;
d[0][2] = 0;
d[0][3] = 0;
d[0][4] = 0;
d[0][5] = 0;
d[0][6] = 0;
d[0][7] = 0;
d[0][8] = 0;
d[1][0] = 0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
d[1][1] = 0;
d[1][2] = 0;
d[1][3] = 0;
d[1][4] = 0;
d[1][5] = 0;
d[1][6] = 0;
d[1][7] = 0;
d[1][8] = 0;
d[2][0] = -1.07679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012346;
d[2][1] = 1.18567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901235;
d[2][2] = 0;
d[2][3] = 0;
d[2][4] = 0;
d[2][5] = 0;
d[2][6] = 0;
d[2][7] = 0;
d[2][8] = 0;
d[3][0] = 0.0408333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
d[3][1] = 0;
d[3][2] = 0.122500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
d[3][3] = 0;
d[3][4] = 0;
d[3][5] = 0;
d[3][6] = 0;
d[3][7] = 0;
d[3][8] = 0;
d[4][0] = 0.638913923625572678050812161599333610995418575593502707205331112036651395251978342357351103706788837984173261141191170345689296126613910870470637234485630987088713036234902124114952103290295710120783007080383173677634319033735943356934610578925447730112453144523115368596418159100374843815077051228655;
d[4][1] = 0;
d[4][2] = -2.45567263822365680966264056643065389421074552269887546855476884631403581840899625156184922948771345272802998750520616409829237817576009995835068721366097459391920033319450229071220324864639733444398167430237401082882132444814660558100791336942940441482715535193669304456476468138275718450645564348188;
d[4][3] = 2.27225871459808413161182840483132028321532694710537276134943773427738442315701790920449812578092461474385672636401499375260308204914618908788004997917534360683048729695960016659725114535610162432319866722199083715118700541441066222407330279050395668471470220741357767596834652228238234069137859225323;
d[4][4] = 0;
d[4][5] = 0;
d[4][6] = 0;
d[4][7] = 0;
d[4][8] = 0;
d[5][0] = -2.66157737501875713111925929786181811927886110621279491101469713930509612941385509594430505557023368554958603063493231787317799674230845336779694159618434634321299954680555583166102818777787674050555772208977921672782109407296127769462859055244305774353176625390759049478562797393162386458986585502767;
d[5][1] = 0;
d[5][2] = 10.8045138864561376956539665536553283848216695313861604686310006136290059392424850728726033575720412843056269529183992807963144686325598710039491665420757292749882852097214396218692807538243757762105076420989844983852423950323217516089733038440975804125387258140136771071745292560309970797058687391147;
d[5][3] = -8.35391465739619941196804854781929169154088257610338993266161781790756725509495082412344107605516870213559043057090330898930263408054993414287063406173197791833013415175411049756858203990308618552762172172286334472369559218815026430538538071155831209364354686322679692202799997090849267285704422656095;
d[5][4] = 0.820487594956656979142041734174383920961870910216630048688253645970444941453995829505408073427394894785640091606891607825314656920843457027957366318480980578695932807037817395645344999121960904086209105241195554096758296602597587867673653286696330637022835832721857989026695559276676282615189064614936;
d[5][5] = 0;
d[5][6] = 0;
d[5][7] = 0;
d[5][8] = 0;
d[6][0] = 6.06774143469677099271836018387727671467926640176887306359909941550772632909985821015131404165728812111281168328257663180016064215895470690198142759429883346899999316325970007167545285470047913584083460611017718063130166260596491115835032527423915877478243890066807867827599446375139672648735819080763;
d[6][1] = 0;
d[6][2] = -24.7112736359110857973420348529074600180341771629677540813665554596587326431239932816061210269509909363811041380357177346279533645394890272496798875375316795980503341386187230715331278943723231323700151408836595698414095201199019330839199397683748689498911238776565489210858226525290942123615012204069;
d[6][3] = 20.4275179307888939404577311174834661269660931649405558961215868278847940059501212187832475382874361032894160252654902848069119955759637823202167209261918250178828181639508851014288273511443537985849961840706061832183847204574232612353992610081640468028660624990850592258752345152578286375382830911574;
d[6][4] = -1.90615797881664715062409678435275701087897468130154449732332455905089888161653727987321175840553252820469326122182450914811549396345902609080269405892288783291419739936108774599029695502963878984287699984458255527075364785177970366374833241148008090451090875831710847373226588365269422264785656055595;
d[6][5] = 1.00617224924206801479004033589947418726779227755986961896919377531711118969055113254477120541179924018356969070947532716899622076802956411828443307596390894408172021076922564441914464355712898778706135054745876126247678490829346435391868589745174427675353123622051949066685955717256307098371649899780;
d[6][6] = 0;
d[6][7] = 0;
d[6][8] = 0;
d[7][0] = 12.0546700762532029950910945289277831164804252410675022456504679687895689139960146039443288283750336379411096035825954447760936334776440213949443362851020478648773614432833479739455614899658192055871179070263606189309615200393669693029537872090663785823585220176083264357375138384872002392464453156991;
d[7][1] = 0;
d[7][2] = -49.7547849504689893280725761533144475832181722583426930533724074758711662694085370746679713708359769175525414749792979166191877040702950632035050470020509416758662857711589354027529455771937047867243267859099031357708288506308006555337876940221477891866864320545092351815961637286145785810654620897866;
d[7][3] = 41.1428886386046766325969841671015735420861993168074196817767197622564349824524439840170095382691179158858647537995469065281803314037669775100891174117560375843405111217539168488777686205987779081576532479718384543237590330774192826000651008343724280693109569745879303050600351583118738080778029866863;
d[7][4] = -4.46176014997400418564191160348481537505077261504527943671205073698503173782183498817306004035577444077774338030344144153116606679374731727848251310165306184178863494208558304215340239111889056504039756992455997594732226828711730693128980332843828825010470766013799305694315017127922117864625331835752;
d[7][5] = 2.04233482223917495982171707770860854373769387362950704675226264453594041136972385119392178374708596422006064736887487517564055580246077638154888337104289952170304921640066775478538257502794816329229535477961389665728962587873947154164966547885933962957589901056567981878726157221164947050255341552053;
d[7][6] = -0.0983484366540610737953080169387022440353735581164564840949921627257463005878103763142287391994861597167501494682778683295607498198293948045947769641969814532660010681934141327023647172799499252723421539433498581938590600776077609795910561717120688444542382881147083210454966691169237581150863097617698;
d[7][7] = 0;
d[7][8] = 0;
d[8][0] = 10.1381465228818078764184514198168903076855566197446625878810302234504456896166594014779323382898853196615353532636712517123358396606673215965993494780266380948135104983713885563789845413826484979295440534003324926587827813571056844823251010555486301677111922636549074237942689192388927633686568909689;
d[8][1] = 0;
d[8][2] = -42.6411360317175021462284600673663573062458353006064366534390081892169164380334190644420973788613121461309320350305871333287662229028588588586066541960004850314875558898834391865978268773660940651739460069750441404794804377117425966233263555407232236441737361866796074494069051528927970522898826975093;
d[8][3] = 35.7638400399225700713502117802316005403429826863865820113565153536300404598612532898965354818454694648874725072682749267373157568233054310410187975820486383990095600494571284321105493325839468611507648538420067322279253828587432646374273647616387422018633170284509612151810792245462886190112774102403;
d[8][4] = -4.34802284039290765334037029690824594371028341639116518865331287370628888237874447914620105044959090476557742632948537833967528093187874711244268408364492971056602670317493656831158047324970779777483136577206179452089410068354309842483825320643538174054369045758176371157191719967939556428314730251243;
d[8][5] = 2.00986226837703589544194359301182755477072972921068058748766475060301022678114409555285903080395717143894554579171389259003194150520251106072972691640555512448806839231389179903663905816067082512966845731905978567106788706342683552320444934837670908237905142357004848267401763746008483191555662392919;
d[8][6] = 0.348749046033827240595382285305314587914041325248182841863580609691330339196441065124674708479092059422639244419745826821459351853291206642381263602952792439394331494334732664490055557979195367934536312102247333268705820911625949222101271883714392689492479451747299183798095407000894534464887545496646;
d[8][7] = -0.271439005104831284237158714091029740757191643592506186496469874451621395043334308463703130107500964514083189383333386192701386007728864369679799299788209315651887841418765697106821139490659689195736303916540408826107333795616038816893578302119868756728613523161845144468638835673968132187348470613286;
d[8][8] = 0;
d[9][0] = -45.0300720342986771243532240507376963515145444719937644824079824393198401727221975824532835519581554198813324085760134180570787821182548227893314627083718036426491592921518503051297688605921310892888780618959209224202919510149533695706542210627124350751677569019920750022536127060269294301531394376504;
d[9][1] = 0;
d[9][2] = 187.327243765458884075241820615420199738351561133306193730472042811504816158488076206432002582649614412920275002256672174701510739653751114685455028887061290838474955260652565481309225324569424103179332013504226474976360325989130909739828895637852868788385760452234355087997545517424616009154634705758;
d[9][3] = -154.028823693501869059672862103451040258175735549044835668465248022257814538167848348770086411169018627838612812103724152543810567982425588659759737860803806260371145723240966294468826535262918904290878915967941233375788211216428688383225716305408979344377644058147926817639142528931780694303713181806;
d[9][4] = 18.5646530634753623385949233295843913676513806454050558851549018133209261097165524415954345328693490705442847387683039834385800075030766667530814210089033384064285859281150173647602629309861028182523307187208812695652827444150045993819313815469829521981879407525662683342644692190497470756749789063511;
d[9][5] = -7.14180967929507885492542049682355119282094179375667055538129282740352869535208365846699842107061792731728890413284789902812520734947617347293091211988712457899696810183968746594422973631393050228738091580973160939356720771138507764901376187221457776450991563372884707638237400874442565117025303276579;
d[9][6] = 1.30880857816137862511476270600769669650828003608402109062757866415544113803750094166293126867882849157267438378760931148892381029332880348348566279309810523711373192846492121947333687661345357443547516144848602064800429953863162648113342205550017119748161538906822547401311450722877269079749204011337;
d[9][7] = 0;
d[9][8] = 0;

  // run the Runge-Kutta56 method
  rV = rk67_methods_d(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rk67_methods_d(double *err, double a[10], double a_minus_b[10], double c[10], double d[10][9], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK67 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, num_digits = 16;
  double cond_num, norm_J, norm_J_inv;

  comp_d  newTime;
  point_d Y, dX[10];

  init_point_d(Y, 0);
  for (i = 0; i < 10; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 10; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_d(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk67_d - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk67_d - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      clear_point_d(Y);
      for (i = 0; i < 10; i++)
        clear_point_d(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk67_d sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk67_d sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_d(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk67_d sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk67_d sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      {
        clear_point_d(Y);
        for (i = 0; i < 10; i++)
          clear_point_d(dX[i]);

        return retVal;
      }
    }

    if (k < 9)
    { // another function evaluation is needed
      // calculate dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_d(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_d(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_d(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rdouble_d(newTime, dT, c[k+1]);
      add_d(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_d(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..7, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..7, a_k * dX[k])
  increase_size_vec_d(Y, dX[0]->size);
  increase_size_vec_d(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rdouble_d(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 10; k++)
    {
      mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_d(&Y->coord[i], &Y->coord[i], newTime);

      mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
      add_d(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_d(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  *err = infNormVec_d(Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_d(OUT, 10, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
      printf("dX = ");  printVec_d(stdout, 10, dX[0]);
      printf("newPoint = ");  printVec_d(stdout, 10, newPoint);
    }
  }

  clear_point_d(Y);
  for (i = 0; i < 10; i++)
    clear_point_d(dX[i]);

  return retVal;
}

///////////////// MP VERSIONS /////////////////////////////////

int heun_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Heun's method - Kincaid/Cheney pg 540-541         *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0;
  int num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;

  comp_mp newTime;
  vec_mp Y, dX1, dX2;

  // initialize MP
  init_mp2(newTime, T->Precision);
  init_vec_mp2(Y, 0, T->Precision); init_vec_mp2(dX1, 0, T->Precision); init_vec_mp2(dX2, 0, T->Precision);

  // find dX1 = f(x,t) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_mp(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_mp(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_mp(&Y->coord[i], &Y->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_mp(dX1, e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX1, e->Jv, Y);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_mp(dX1, e->Jv, Y);
  }
  // matrixSolve calculated dX1 = f(x, t) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_mp - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_mp - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    // clear MP
    clear_mp(newTime);
    clear_vec_mp(Y); clear_vec_mp(dX1); clear_vec_mp(dX2); 

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
        fprintf(stderr, "heun_mp sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: heun_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "heun_mp sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: heun_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    if (retVal)
    { // clear MP
      clear_mp(newTime);
      clear_vec_mp(Y); clear_vec_mp(dX1); clear_vec_mp(dX2);

      return retVal;
    }
  }

  // find Y = x + dT * dX
  increase_size_vec_mp(Y, dX1->size);
  rows = Y->size = dX1->size;
  for (i = 0; i < rows; i++)
  {
    mul_mp(&Y->coord[i], dT, &dX1->coord[i]);
    add_mp(&Y->coord[i], &Y->coord[i], &P->point->coord[i]);
  }

  // find newTime = t + dT
  add_mp(newTime, dT, P->time);

  // do the second function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_mp(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_mp(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_mp(&Y->coord[i], &Y->coord[i]);
  }

  // just do matrixSolve
  retVal = matrixSolve_mp(dX2, e->Jv, Y);

  // matrixSolve calculated dX2 = f(x + dT*dX1,t + dT) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_mp - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_mp - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    // clear MP
    clear_mp(newTime);
    clear_vec_mp(Y); clear_vec_mp(dX1); clear_vec_mp(dX2);

    return retVal;
  }

  // find newPoint = P->point + 1/2*dT*(dX1 + dX2)
  increase_size_vec_mp(newPoint, dX1->size);
  div_ui_mp(newTime, dT, 2);
  rows = newPoint->size = dX1->size;
  for (i = 0; i < rows; i++)
  {
    add_mp(&dX1->coord[i], &dX1->coord[i], &dX2->coord[i]);
    mul_mp(&dX1->coord[i], &dX1->coord[i], newTime);
    add_mp(&newPoint->coord[i], &dX1->coord[i], &P->point->coord[i]);
  }

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX1);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX1);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  // clear MP
  clear_mp(newTime);
  clear_vec_mp(Y); clear_vec_mp(dX1); clear_vec_mp(dX2);

  return retVal;
}

int runge_kutta_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the Runge Kutta method - Kincaid/Cheney pg 541    *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;
  int num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;
 
  comp_mp newTime;
  vec_mp  Y, dX[4];
 
  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 4; i++)
    init_point_mp(dX[i], 0);

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 4; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_mp(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y
 
    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in runge_kutta_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in runge_kutta_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }
   
      // clear MP
      clear_mp(newTime);
      clear_point_mp(Y);
      for (i = 0; i < 4; i++)
        clear_point_mp(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "runge_kutta_mp sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: runge_kutta_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
 
      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "runge_kutta_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: runge_kutta_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
   
      if (retVal)
      { // clear MP
        clear_mp(newTime);
        clear_point_mp(Y);
        for (i = 0; i < 4; i++)
          clear_point_mp(dX[i]);

        return retVal;
      }
    }

    // calculate dX[k] = dT * dX[k]
    vec_mulcomp_mp(dX[k], dX[k], dT);

    if (k < 3)
    { // another function evaluation is needed

      // find Y = x + a * dX[k], where a = 1/2 when k = 0 || k == 1 and a = 1 when k == 2
      increase_size_vec_mp(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      if (k == 0 || k == 1)
      { // find Y
        for (i = 0; i < rows; i++)
        {
          mpf_div_ui(Y->coord[i].r, dX[k]->coord[i].r, 2);
          mpf_div_ui(Y->coord[i].i, dX[k]->coord[i].i, 2);
          add_mp(&Y->coord[i], &Y->coord[i], &P->point->coord[i]);
        }
 
        // find newTime = t + dT / 2
        mpf_div_ui(newTime->r, dT->r, 2);
        mpf_div_ui(newTime->i, dT->i, 2);
        add_mp(newTime, newTime, P->time);
      }
      else // k == 2
      { // find Y
        for (i = 0; i < rows; i++)
        {
          add_mp(&Y->coord[i], &dX[k]->coord[i], &P->point->coord[i]);
        }

        // find newTime = t + dT
        add_mp(newTime, dT, P->time);
      }
      
      // do the second function evaluation
      eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
  }

  // find newPoint = P->point + 1/6*(dX[0] + dX[3] + 2*(dX[1] + dX[2]))
  change_size_vec_mp(newPoint, dX[0]->size);
  rows = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    add_mp(&dX[1]->coord[i], &dX[1]->coord[i], &dX[2]->coord[i]);
    mul_ui_mp(&dX[1]->coord[i], &dX[1]->coord[i], 2);
    add_mp(&dX[0]->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
    add_mp(&dX[0]->coord[i], &dX[0]->coord[i], &dX[3]->coord[i]);
    div_ui_mp(&dX[0]->coord[i], &dX[0]->coord[i], 6);
    add_mp(&newPoint->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
  }

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX[0]);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }
 
  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 4; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int heun_error_mp(mpf_t error, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Heun's method with error approximation            *
\***************************************************************/
{
  int i, j, rows, cols, retVal = 0, num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;

  comp_mp newTime;
  point_mp Y, dX1, dX2;

  init_mp(newTime);
  init_point_mp(Y, 0); init_point_mp(dX1, 0); init_point_mp(dX2, 0);

  // find f(x,t) = -(Jv^-1)*(Jp*dp/dt) = dX

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_mp(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_mp(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_mp(&Y->coord[i], &Y->coord[i]);
  }

  if (T->MPType == 2) // using AMP
  {
    if (T->steps_since_last_CN > -1)
    { // do matrixSolve & cond_num together
      retVal = matrixSolve_cond_num_norms_mp(dX1, e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX1, e->Jv, Y);
      T->steps_since_last_CN++;
    }
  }
  else
  { // just do matrixSolve
    retVal = matrixSolve_mp(dX1, e->Jv, Y);
  }
  // matrixSolve calculated dX1 = f(x,t) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_error_mp - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_error_mp - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_mp(newTime);
    clear_point_mp(Y); clear_point_mp(dX1); clear_point_mp(dX2);

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
        fprintf(stderr, "heun_error_mp sees that AMP Criterion A has been violated!\n");
      fprintf(OUT, "NOTE: heun_error_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }

    // check criterion C from the AMP paper.
    if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
    {
      if (T->screenOut)
        fprintf(stderr, "heun_error_mp sees that AMP Criterion C has been violated!\n");
      fprintf(OUT, "NOTE: heun_error_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    if (retVal)
    {
      clear_mp(newTime);
      clear_point_mp(Y); clear_point_mp(dX1); clear_point_mp(dX2);

      return retVal;
    }
  }

  // setup dX1 = dT * dX1
  // find Y = x + dX1
  increase_size_vec_mp(Y, dX1->size);
  rows = Y->size = dX1->size;
  for (i = 0; i < rows; i++)
  {
    mul_mp(&dX1->coord[i], dT, &dX1->coord[i]);
    add_mp(&Y->coord[i], &P->point->coord[i], &dX1->coord[i]);
  }

  // find newTime = t + dT
  add_mp(newTime, dT, P->time);

  // do the second function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);

  // find Y = -dH/dt = -dH/dp * dp/dt
  increase_size_vec_mp(Y, e->Jp->rows);
  rows = Y->size = e->Jp->rows;
  cols = e->Jp->cols;
  for (i = 0; i < rows; i++)
  {
    set_zero_mp(&Y->coord[i]);
    for (j = 0; j < cols; j++)
    {
      sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
    }
    neg_mp(&Y->coord[i], &Y->coord[i]);
  }

  // just do matrixSolve
  retVal = matrixSolve_mp(dX2, e->Jv, Y);

  // matrixSolve calculated dX2 = f(x + dT*dX1,t + dT) = (dH/dX)^(-1) * Y

  // check for matrixSolve failures
  if (retVal)
  {
    if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
    {
      fprintf(OUT, "NOTE: matrixSolve has failed in heun_error_mp - the precision will be increased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed in heun_error_mp - the precision will be increased.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    {
      fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
    }

    clear_mp(newTime);
    clear_point_mp(Y); clear_point_mp(dX1); clear_point_mp(dX2);

    return retVal;
  }

  // setup dX2 = dT*dX2
  // compute error: 1/2*(dX1 - dX2)
  // find newPoint = P->point + 1/2*(dX1 + dX2)
  increase_size_vec_mp(Y, dX2->size);
  increase_size_vec_mp(newPoint, dX2->size);
  rows = newPoint->size = Y->size = dX2->size;
  for (i = 0; i < rows; i++)
  {
    mul_mp(&dX2->coord[i], &dX2->coord[i], dT);

    sub_mp(&Y->coord[i], &dX1->coord[i], &dX2->coord[i]);
    add_mp(&dX1->coord[i], &dX1->coord[i], &dX2->coord[i]);

    div_ui_mp(&Y->coord[i], &Y->coord[i], 2);
    div_ui_mp(&dX1->coord[i], &dX1->coord[i], 2);

    add_mp(&newPoint->coord[i], &dX1->coord[i], &P->point->coord[i]);
  }
  infNormVec_mp2(error, Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX1);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX1);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  clear_mp(newTime);
  clear_point_mp(Y); clear_point_mp(dX1); clear_point_mp(dX2);

  return retVal;
}

int rkn34_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKN34 method - Norsett method                 *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[5] = {25, 32, 256, 0, 11};
  int a_denom[5] = {162, 135, 567, 1, 70};
  int a_minus_b_numer[5] = {-41, -244, 256, -448, 11};
  int a_minus_b_denom[5] = {4050, 1755, 567, 975, 70};
  int c_numer[5] = {0, 3, 9, 25, 1};
  int c_denom[5] = {1, 8, 16, 32, 1};
  int d_numer[5][4] = {{0, 0, 0, 0},
                       {3, 0, 0, 0},
                       {0, 9, 0, 0},
                       {-125, 325, 0, 0},
                       {371, -200, 1120, 0}};
  int d_denom[5][4] = {{1, 1, 1, 1},
                       {8, 1, 1, 1},
                       {1, 16, 1, 1},
                       {672, 336, 1, 1},
                       {891, 297, 891, 1}};

  mpf_t a[5], a_minus_b[5], c[5], d[5][4];
  for (i = 0; i < 5; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 4; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta34 method
  rV = rk34_methods_mp(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 5; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 4; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk34_methods_mp(mpf_t err, mpf_t a[5], mpf_t a_minus_b[5], mpf_t c[5], mpf_t d[5][4], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK34 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;
  int num_digits = (int) floor(T->Precision * log10(2.0));
  double cond_num, norm_J, norm_J_inv;

  comp_mp newTime;
  point_mp Y, dX[5];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 5; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 5; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_mp(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk34_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk34_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      // clear MP
      clear_mp(newTime);
      clear_point_mp(Y);
      for (i = 0; i < 5; i++)
        clear_point_mp(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk34_mp sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk34_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk34_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk34_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      { // clear MP
        clear_mp(newTime);
        clear_point_mp(Y);
        for (i = 0; i < 5; i++)
          clear_point_mp(dX[i]);

        return retVal;
      }
    }

    if (k < 4)
    { // another function evaluation is needed
 
      // find dX[k] = dX[k] * dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_mp(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_mp(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_mp(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rmpf_mp(newTime, dT, c[k+1]);
      add_mp(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_mp(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_mp(Y, dX[0]->size);
  increase_size_vec_mp(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rmpf_mp(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 5; k++)
    {
      mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_mp(&Y->coord[i], &Y->coord[i], newTime);

      mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
      add_mp(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_mp(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  infNormVec_mp2(err, Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX[0]);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 5; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkf45_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKF45 method - Kincaid/Cheney pg 544-545      *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[6] = {16, 0, 6656, 28561, -9, 2};
  int a_denom[6] = {135, 1, 12825, 56430, 50, 55};
  int a_minus_b_numer[6] = {1, 0, -128, -2197, 1, 2};
  int a_minus_b_denom[6] = {360, 1, 4275, 75240, 50, 55};
  int c_numer[6] = {0, 1, 3, 12, 1, 1};
  int c_denom[6] = {1, 4, 8, 13, 1, 2};
  int d_numer[6][5] = {{0, 0, 0, 0, 0},
                       {1, 0, 0, 0, 0},
                       {3, 9, 0, 0, 0}, 
                       {1932, -7200, 7296, 0, 0},
                       {439, -8, 3680, -845, 0},
                       {-8, 2, -3544, 1859, -11}};
  int d_denom[6][5] = {{1, 1, 1, 1, 1},
                       {4, 1, 1, 1, 1},
                       {32, 32, 1, 1, 1},
                       {2197, 2197, 2197, 1, 1},
                       {216, 1, 513, 4104, 1},
                       {27, 1, 2565, 4104, 40}};

  mpf_t a[6], a_minus_b[6], c[6], d[6][5];
  for (i = 0; i < 6; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 5; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta45 method
  rV = rk45_methods_mp(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 6; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]); 
    mpf_clear(c[i]); 
    for (j = 0; j < 5; j++)
      mpf_clear(d[i][j]);
  } 

  return rV;
}

int rkck45_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKCK45 method - Cash/Karp method              *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[6] = {37, 0, 250, 125, 0, 512};
  int a_denom[6] = {378, 1, 621, 594, 1, 1771};
  int a_minus_b_numer[6] = {-277, 0, 6925, -6925, -277, 277};
  int a_minus_b_denom[6] = {64512, 1, 370944, 202752, 14336, 7084};
  int c_numer[6] = {0, 1, 3, 3, 1, 7};
  int c_denom[6] = {1, 5, 10, 5, 1, 8};
  int d_numer[6][5] = {{0, 0, 0, 0, 0},
                       {1, 0, 0, 0, 0},
                       {3, 9, 0, 0, 0},
                       {3, -9, 6, 0, 0},
                       {-11, 5, -70, 35, 0},
                       {1631, 175, 575, 44275, 253}};
  int d_denom[6][5] = {{1, 1, 1, 1, 1},
                       {5, 1, 1, 1, 1},
                       {40, 40, 1, 1, 1},
                       {10, 10, 5, 1, 1},
                       {54, 2, 27, 27, 1},
                       {55296, 512, 13824, 110592, 4096}};

  mpf_t a[6], a_minus_b[6], c[6], d[6][5];
  for (i = 0; i < 6; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 5; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta45 method
  rV = rk45_methods_mp(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 6; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 5; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk45_methods_mp(mpf_t err, mpf_t a[6], mpf_t a_minus_b[6], mpf_t c[6], mpf_t d[6][5], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;
  int num_digits = (int) floor(T->Precision * log10(2.0));
  double cond_num, norm_J, norm_J_inv;

  comp_mp newTime;
  point_mp Y, dX[6];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 6; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 6; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_mp(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk45_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk45_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      // clear MP
      clear_mp(newTime);
      clear_point_mp(Y);
      for (i = 0; i < 6; i++)
        clear_point_mp(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk45_mp sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk45_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk45_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk45_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      { // clear MP
        clear_mp(newTime);
        clear_point_mp(Y);
        for (i = 0; i < 6; i++)
          clear_point_mp(dX[i]);

        return retVal;
      }
    }

    if (k < 5)
    { // another function evaluation is needed
      // calculate dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_mp(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_mp(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_mp(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rmpf_mp(newTime, dT, c[k+1]);
      add_mp(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_mp(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_mp(Y, dX[0]->size);
  increase_size_vec_mp(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rmpf_mp(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 6; k++)
    {
      mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_mp(&Y->coord[i], &Y->coord[i], newTime);

      mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
      add_mp(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_mp(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  infNormVec_mp2(err, Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {       
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX[0]);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 6; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkdp56_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKDP56 method - Dormand/Prince method         *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[8] = {821, 0, 19683, 175273, 395, 785, 3, 0};
  int a_denom[8] = {10800, 1, 71825, 912600, 3672, 2704, 50, 1};
  int a_minus_b_numer[8] = {13, 0, -19683, 2401, -65, 15, 521, -1};
  int a_minus_b_denom[8] = {2400, 1, 618800, 31200, 816, 416, 5600, 10};
  int c_numer[8] = {0, 1, 2, 3, 3, 4, 1, 1};
  int c_denom[8] = {1, 10, 9, 7, 5, 5, 1, 1};
  int d_numer[8][7] = {{0, 0, 0, 0, 0, 0, 0},
                          {1, 0, 0, 0, 0, 0, 0},
                          {-2, 20, 0, 0, 0, 0, 0},
                          {615, -270, 1053, 0, 0, 0, 0},
                          {3243, -54, 50949, 4998, 0, 0, 0},
                          {-26492, 72, 2808, -24206, 338, 0, 0},
                          {5561, -35, -24117, 899983, -5225, 3925, 0},
                          {465467, -2945, -5610201, 10513573, -424325, 376225, 0}};
  int d_denom[8][7] = {{1, 1, 1, 1, 1, 1, 1},
                       {10, 1, 1, 1, 1, 1, 1},
                       {81, 81, 1, 1, 1, 1, 1},
                       {1372, 343, 1372, 1, 1, 1, 1},
                       {5500, 55, 71500, 17875, 1, 1, 1},
                       {37125, 55, 23375, 37125, 459, 1, 1},
                       {2376, 11, 31603, 200772, 1836, 4056, 1},
                       {266112, 1232, 14158144, 3212352, 205632, 454272, 1}};

  mpf_t a[8], a_minus_b[8], c[8], d[8][7];
  for (i = 0; i < 8; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 7; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta56 method
  rV = rk56_methods_mp(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 8; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 7; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk56_methods_mp(mpf_t err, mpf_t a[8], mpf_t a_minus_b[8], mpf_t c[8], mpf_t d[8][7], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;
  int num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;

  comp_mp newTime;
  point_mp Y, dX[8];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 8; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 8; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_mp(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk56_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk56_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      // clear MP
      clear_mp(newTime);
      clear_point_mp(Y);
      for (i = 0; i < 8; i++)
        clear_point_mp(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk56_mp sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk56_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk56_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk56_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      { // clear MP
        clear_mp(newTime);
        clear_point_mp(Y);
        for (i = 0; i < 8; i++)
          clear_point_mp(dX[i]);

        return retVal;
      }
    }

    if (k < 7)
    { // another function evaluation is needed
      // calculate dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_mp(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_mp(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_mp(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rmpf_mp(newTime, dT, c[k+1]);
      add_mp(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_mp(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_mp(Y, dX[0]->size);
  increase_size_vec_mp(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rmpf_mp(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 8; k++)
    {
      mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_mp(&Y->coord[i], &Y->coord[i], newTime);

      mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
      add_mp(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_mp(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  infNormVec_mp2(err, Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX[0]);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 8; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkv67_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKV67 method - Verner method                  *
\***************************************************************/
{
  int i, j, rV;

  mpf_t a[10], a_minus_b[10], c[10], d[10][9];
  for (i = 0; i < 10; i++)
  {
    mpf_init(a[i]); mpf_init(a_minus_b[i]); mpf_init(c[i]); 
    for (j = 0; j < 9; j++)
      mpf_init(d[i][j]);
  }

mpf_set_str(a[0],"0.0471556184862722217043176510883817567956932673866498263717863225066525727846555841973346401995827631062036755590558062534375047415234260388676624358999573840784918991178510655224244139232889913788333823544199823935451154623648793450503723002665841976598046959664007442181911565685230727013689593392982",10);
mpf_set_str(a[1],"0",10);
mpf_set_str(a[2],"0",10);
mpf_set_str(a[3],"0.257505642984341518959643610103768758098639571258883586561584933570028952897166612637862059798925621609396907594245725006327476849765387793842037461290921370430264152893947313918692476051185173305706391851963836018457313459645894963869257643822167873501457534217624142946250030154297352910582419606143",10);
mpf_set_str(a[4],"0.262166539774126204771386309576452771112914046736971319203584006643876812614319887681344771902494347386463073904093142446268038453678784950161124772440805335094258253338958758076230419920570973073193926469534220996564683552221244076754949769609944998538014278529773548771150093620971914747932175414550",10);
mpf_set_str(a[5],"0.152160926567385574032313319916511753552259518225131674367096139045695165543077407261967511436494354991858803879708884258567824761923036313221569840599729116798542290424111444095404797846598555249340492159421091497761496624438637702252254645961603915777481188523759942071760307205888241753387037725744",10);
mpf_set_str(a[6],"0.493996917003248424690717589322787684429585308571357948791678035642674335879997785662081156612226972186240516630330279108098440721497194764675468796945045160954090374742547119955506272386697362858820688622832619000271491070870242751247344166824766097061388163627211199461983342148211784138952128223367",10);
mpf_set_str(a[7],"-0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817",10);
mpf_set_str(a[8],"0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153",10);
mpf_set_str(a[9],"0",10);

mpf_set_str(a_minus_b[0],"0.00254701187993104541699947511358977898137393503651067684217948288343659637043927033143148723131125158014121649938900878979858451840352977282793140381518240008963362033067932018366864275478963166454694523116781829305324081676250402401051739684741551799350892643994927331531661450056821725901262420206562",10);
mpf_set_str(a_minus_b[1],"0",10);
mpf_set_str(a_minus_b[2],"0",10);
mpf_set_str(a_minus_b[3],"-0.00965839487279574909126661599061503187520108884519528045821452710507663251123547574163069610539919389963989845238921904972854937606654792208299092082065749981679263365462373361333668366236519725184736326271819113885999268066494138424342590942206813226759750600367884910177105022385487684166900829606156",10);
mpf_set_str(a_minus_b[4],"0.0420647097563969027734147319113774614805948824559621257010925807491537576600651305254295976225991591885211868879798925804616780192287010286942291836640369052613351689178787207331605464157207419903756767065167872280544436640338531219590767268952538968560279938316779313158803412904221702386669877100151",10);
mpf_set_str(a_minus_b[5],"-0.0666822437469301090659987634347776289054953577328478620986708910886522673765780883582819579617699092650805636589852395875931279767916652676612036293647748766874420794927278867725689176383726462247671413895793673537650626399046027049371138790153665954418043148360367277640961073763921508851946773219422",10);
mpf_set_str(a_minus_b[6],"0.265009746462128136352900200346432447893361385631363781535494557795019991145447029415600259733188324877703429783214090406452549289091431065680284500120422909432522977485881928123107977267794033274265297233044517672944951473291445755463803660526233224001811988404937384560301642803203641418032914674657",10);
mpf_set_str(a_minus_b[7],"-0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817",10);
mpf_set_str(a_minus_b[8],"0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153",10);
mpf_set_str(a_minus_b[9],"-0.0202951846633562822276705479381043035855420443667990862261517658249536055689205887319585505702055732014823934917746960666918489454776188166903872302377514709236100830696726470857731850092255075866785330602598147948274804639773599730786794693464008286038012269720794348562965112960546349366261206596328",10);

mpf_set_str(c[0],"0",10);
mpf_set_str(c[1],"0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[2],"0.108888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889",10);
mpf_set_str(c[3],"0.163333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333",10);
mpf_set_str(c[4],"0.455500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[5],"0.609509448997838131708700442148602494963796759286605673642939302386787496187674982310265299374033791406090583319455261759148494730544940521238957202640385592141084318199590688285015525265373754263537303527537491030484005373807797476632985866792541212386248529601147679387596870467556824874147722141014",10);
mpf_set_str(c[6],"0.884000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[7],"0.925000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[8],"1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[9],"1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);

mpf_set_str(d[0][0],"0",10);
mpf_set_str(d[0][1],"0",10);
mpf_set_str(d[0][2],"0",10);
mpf_set_str(d[0][3],"0",10);
mpf_set_str(d[0][4],"0",10);
mpf_set_str(d[0][5],"0",10);
mpf_set_str(d[0][6],"0",10);
mpf_set_str(d[0][7],"0",10);
mpf_set_str(d[0][8],"0",10);
mpf_set_str(d[1][0],"0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(d[1][1],"0",10);
mpf_set_str(d[1][2],"0",10);
mpf_set_str(d[1][3],"0",10);
mpf_set_str(d[1][4],"0",10);
mpf_set_str(d[1][5],"0",10);
mpf_set_str(d[1][6],"0",10);
mpf_set_str(d[1][7],"0",10);
mpf_set_str(d[1][8],"0",10);
mpf_set_str(d[2][0],"-1.07679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012346",10);
mpf_set_str(d[2][1],"1.18567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901235",10);
mpf_set_str(d[2][2],"0",10);
mpf_set_str(d[2][3],"0",10);
mpf_set_str(d[2][4],"0",10);
mpf_set_str(d[2][5],"0",10);
mpf_set_str(d[2][6],"0",10);
mpf_set_str(d[2][7],"0",10);
mpf_set_str(d[2][8],"0",10);
mpf_set_str(d[3][0],"0.0408333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333",10);
mpf_set_str(d[3][1],"0",10);
mpf_set_str(d[3][2],"0.122500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(d[3][3],"0",10);
mpf_set_str(d[3][4],"0",10);
mpf_set_str(d[3][5],"0",10);
mpf_set_str(d[3][6],"0",10);
mpf_set_str(d[3][7],"0",10);
mpf_set_str(d[3][8],"0",10);
mpf_set_str(d[4][0],"0.638913923625572678050812161599333610995418575593502707205331112036651395251978342357351103706788837984173261141191170345689296126613910870470637234485630987088713036234902124114952103290295710120783007080383173677634319033735943356934610578925447730112453144523115368596418159100374843815077051228655",10);
mpf_set_str(d[4][1],"0",10);
mpf_set_str(d[4][2],"-2.45567263822365680966264056643065389421074552269887546855476884631403581840899625156184922948771345272802998750520616409829237817576009995835068721366097459391920033319450229071220324864639733444398167430237401082882132444814660558100791336942940441482715535193669304456476468138275718450645564348188",10);
mpf_set_str(d[4][3],"2.27225871459808413161182840483132028321532694710537276134943773427738442315701790920449812578092461474385672636401499375260308204914618908788004997917534360683048729695960016659725114535610162432319866722199083715118700541441066222407330279050395668471470220741357767596834652228238234069137859225323",10);
mpf_set_str(d[4][4],"0",10);
mpf_set_str(d[4][5],"0",10);
mpf_set_str(d[4][6],"0",10);
mpf_set_str(d[4][7],"0",10);
mpf_set_str(d[4][8],"0",10);
mpf_set_str(d[5][0],"-2.66157737501875713111925929786181811927886110621279491101469713930509612941385509594430505557023368554958603063493231787317799674230845336779694159618434634321299954680555583166102818777787674050555772208977921672782109407296127769462859055244305774353176625390759049478562797393162386458986585502767",10);
mpf_set_str(d[5][1],"0",10);
mpf_set_str(d[5][2],"10.8045138864561376956539665536553283848216695313861604686310006136290059392424850728726033575720412843056269529183992807963144686325598710039491665420757292749882852097214396218692807538243757762105076420989844983852423950323217516089733038440975804125387258140136771071745292560309970797058687391147",10);
mpf_set_str(d[5][3],"-8.35391465739619941196804854781929169154088257610338993266161781790756725509495082412344107605516870213559043057090330898930263408054993414287063406173197791833013415175411049756858203990308618552762172172286334472369559218815026430538538071155831209364354686322679692202799997090849267285704422656095",10);
mpf_set_str(d[5][4],"0.820487594956656979142041734174383920961870910216630048688253645970444941453995829505408073427394894785640091606891607825314656920843457027957366318480980578695932807037817395645344999121960904086209105241195554096758296602597587867673653286696330637022835832721857989026695559276676282615189064614936",10);
mpf_set_str(d[5][5],"0",10);
mpf_set_str(d[5][6],"0",10);
mpf_set_str(d[5][7],"0",10);
mpf_set_str(d[5][8],"0",10);
mpf_set_str(d[6][0],"6.06774143469677099271836018387727671467926640176887306359909941550772632909985821015131404165728812111281168328257663180016064215895470690198142759429883346899999316325970007167545285470047913584083460611017718063130166260596491115835032527423915877478243890066807867827599446375139672648735819080763",10);
mpf_set_str(d[6][1],"0",10);
mpf_set_str(d[6][2],"-24.7112736359110857973420348529074600180341771629677540813665554596587326431239932816061210269509909363811041380357177346279533645394890272496798875375316795980503341386187230715331278943723231323700151408836595698414095201199019330839199397683748689498911238776565489210858226525290942123615012204069",10);
mpf_set_str(d[6][3],"20.4275179307888939404577311174834661269660931649405558961215868278847940059501212187832475382874361032894160252654902848069119955759637823202167209261918250178828181639508851014288273511443537985849961840706061832183847204574232612353992610081640468028660624990850592258752345152578286375382830911574",10);
mpf_set_str(d[6][4],"-1.90615797881664715062409678435275701087897468130154449732332455905089888161653727987321175840553252820469326122182450914811549396345902609080269405892288783291419739936108774599029695502963878984287699984458255527075364785177970366374833241148008090451090875831710847373226588365269422264785656055595",10);
mpf_set_str(d[6][5],"1.00617224924206801479004033589947418726779227755986961896919377531711118969055113254477120541179924018356969070947532716899622076802956411828443307596390894408172021076922564441914464355712898778706135054745876126247678490829346435391868589745174427675353123622051949066685955717256307098371649899780",10);
mpf_set_str(d[6][6],"0",10);
mpf_set_str(d[6][7],"0",10);
mpf_set_str(d[6][8],"0",10);
mpf_set_str(d[7][0],"12.0546700762532029950910945289277831164804252410675022456504679687895689139960146039443288283750336379411096035825954447760936334776440213949443362851020478648773614432833479739455614899658192055871179070263606189309615200393669693029537872090663785823585220176083264357375138384872002392464453156991",10);
mpf_set_str(d[7][1],"0",10);
mpf_set_str(d[7][2],"-49.7547849504689893280725761533144475832181722583426930533724074758711662694085370746679713708359769175525414749792979166191877040702950632035050470020509416758662857711589354027529455771937047867243267859099031357708288506308006555337876940221477891866864320545092351815961637286145785810654620897866",10);
mpf_set_str(d[7][3],"41.1428886386046766325969841671015735420861993168074196817767197622564349824524439840170095382691179158858647537995469065281803314037669775100891174117560375843405111217539168488777686205987779081576532479718384543237590330774192826000651008343724280693109569745879303050600351583118738080778029866863",10);
mpf_set_str(d[7][4],"-4.46176014997400418564191160348481537505077261504527943671205073698503173782183498817306004035577444077774338030344144153116606679374731727848251310165306184178863494208558304215340239111889056504039756992455997594732226828711730693128980332843828825010470766013799305694315017127922117864625331835752",10);
mpf_set_str(d[7][5],"2.04233482223917495982171707770860854373769387362950704675226264453594041136972385119392178374708596422006064736887487517564055580246077638154888337104289952170304921640066775478538257502794816329229535477961389665728962587873947154164966547885933962957589901056567981878726157221164947050255341552053",10);
mpf_set_str(d[7][6],"-0.0983484366540610737953080169387022440353735581164564840949921627257463005878103763142287391994861597167501494682778683295607498198293948045947769641969814532660010681934141327023647172799499252723421539433498581938590600776077609795910561717120688444542382881147083210454966691169237581150863097617698",10);
mpf_set_str(d[7][7],"0",10);
mpf_set_str(d[7][8],"0",10);
mpf_set_str(d[8][0],"10.1381465228818078764184514198168903076855566197446625878810302234504456896166594014779323382898853196615353532636712517123358396606673215965993494780266380948135104983713885563789845413826484979295440534003324926587827813571056844823251010555486301677111922636549074237942689192388927633686568909689",10);
mpf_set_str(d[8][1],"0",10);
mpf_set_str(d[8][2],"-42.6411360317175021462284600673663573062458353006064366534390081892169164380334190644420973788613121461309320350305871333287662229028588588586066541960004850314875558898834391865978268773660940651739460069750441404794804377117425966233263555407232236441737361866796074494069051528927970522898826975093",10);
mpf_set_str(d[8][3],"35.7638400399225700713502117802316005403429826863865820113565153536300404598612532898965354818454694648874725072682749267373157568233054310410187975820486383990095600494571284321105493325839468611507648538420067322279253828587432646374273647616387422018633170284509612151810792245462886190112774102403",10);
mpf_set_str(d[8][4],"-4.34802284039290765334037029690824594371028341639116518865331287370628888237874447914620105044959090476557742632948537833967528093187874711244268408364492971056602670317493656831158047324970779777483136577206179452089410068354309842483825320643538174054369045758176371157191719967939556428314730251243",10);
mpf_set_str(d[8][5],"2.00986226837703589544194359301182755477072972921068058748766475060301022678114409555285903080395717143894554579171389259003194150520251106072972691640555512448806839231389179903663905816067082512966845731905978567106788706342683552320444934837670908237905142357004848267401763746008483191555662392919",10);
mpf_set_str(d[8][6],"0.348749046033827240595382285305314587914041325248182841863580609691330339196441065124674708479092059422639244419745826821459351853291206642381263602952792439394331494334732664490055557979195367934536312102247333268705820911625949222101271883714392689492479451747299183798095407000894534464887545496646",10);
mpf_set_str(d[8][7],"-0.271439005104831284237158714091029740757191643592506186496469874451621395043334308463703130107500964514083189383333386192701386007728864369679799299788209315651887841418765697106821139490659689195736303916540408826107333795616038816893578302119868756728613523161845144468638835673968132187348470613286",10);
mpf_set_str(d[8][8],"0",10);
mpf_set_str(d[9][0],"-45.0300720342986771243532240507376963515145444719937644824079824393198401727221975824532835519581554198813324085760134180570787821182548227893314627083718036426491592921518503051297688605921310892888780618959209224202919510149533695706542210627124350751677569019920750022536127060269294301531394376504",10);
mpf_set_str(d[9][1],"0",10);
mpf_set_str(d[9][2],"187.327243765458884075241820615420199738351561133306193730472042811504816158488076206432002582649614412920275002256672174701510739653751114685455028887061290838474955260652565481309225324569424103179332013504226474976360325989130909739828895637852868788385760452234355087997545517424616009154634705758",10);
mpf_set_str(d[9][3],"-154.028823693501869059672862103451040258175735549044835668465248022257814538167848348770086411169018627838612812103724152543810567982425588659759737860803806260371145723240966294468826535262918904290878915967941233375788211216428688383225716305408979344377644058147926817639142528931780694303713181806",10);
mpf_set_str(d[9][4],"18.5646530634753623385949233295843913676513806454050558851549018133209261097165524415954345328693490705442847387683039834385800075030766667530814210089033384064285859281150173647602629309861028182523307187208812695652827444150045993819313815469829521981879407525662683342644692190497470756749789063511",10);
mpf_set_str(d[9][5],"-7.14180967929507885492542049682355119282094179375667055538129282740352869535208365846699842107061792731728890413284789902812520734947617347293091211988712457899696810183968746594422973631393050228738091580973160939356720771138507764901376187221457776450991563372884707638237400874442565117025303276579",10);
mpf_set_str(d[9][6],"1.30880857816137862511476270600769669650828003608402109062757866415544113803750094166293126867882849157267438378760931148892381029332880348348566279309810523711373192846492121947333687661345357443547516144848602064800429953863162648113342205550017119748161538906822547401311450722877269079749204011337",10);
mpf_set_str(d[9][7],"0",10);
mpf_set_str(d[9][8],"0",10);

  // run the Runge-Kutta67 method
  rV = rk67_methods_mp(err, a, a_minus_b, c, d, newParam, newPoint, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 10; i++)
  {
    mpf_clear(a[i]); mpf_clear(a_minus_b[i]); mpf_clear(c[i]);
    for (j = 0; j < 9; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk67_methods_mp(mpf_t err, mpf_t a[10], mpf_t a_minus_b[10], mpf_t c[10], mpf_t d[10][9], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;
  int num_digits = prec_to_digits(T->Precision);
  double cond_num, norm_J, norm_J_inv;
  
  comp_mp newTime;
  point_mp Y, dX[10];
  
  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 10; i++)
    init_point_mp(dX[i], 0);
    
  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)
  
  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);
  
  for (k = 0; k < 10; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }
    
    if (k == 0 && T->MPType == 2) // using AMP
    {
      if (T->steps_since_last_CN > -1)
      { // do matrixSolve & cond_num together
        retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, &norm_J, &norm_J_inv);
        T->steps_since_last_CN = 0;
        T->latest_cond_num_exp = ceil(log10(cond_num));
      }
      else
      { // just do matrixSolve
        retVal = matrixSolve_mp(dX[k], e->Jv, Y);
        T->steps_since_last_CN++;
      }
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      if (T->MPType == 2)  // adaptive precision - should trigger an increase in precision to make matrixSolve work
      {
        fprintf(OUT, "NOTE: matrixSolve has failed in rk67_mp - the precision will be increased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed in rk67_mp - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      else
      {
        fprintf(OUT, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
        if (T->screenOut)
          fprintf(stderr, "NOTE: matrixSolve has failed - the step size will be decreased.\n");
      }

      // clear MP
      clear_mp(newTime);
      clear_point_mp(Y);
      for (i = 0; i < 10; i++)
        clear_point_mp(dX[i]);

      return retVal;
    }

    // check to make sure that rules A & C are satisfied for AMP on the first function evaluation
    if (k == 0 && T->MPType == 2)
    {
      retVal = 0;
      // check criterion A from the AMP paper.
      if (num_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
      {
        if (T->screenOut)
          fprintf(stderr, "rk67_mp sees that AMP Criterion A has been violated!\n");
        fprintf(OUT, "NOTE: rk67_mp sees that AMP Criterion A has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }

      // check criterion C from the AMP paper.
      if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, infNormVec_mp(P->point)))
      {
        if (T->screenOut)
          fprintf(stderr, "rk67_mp sees that AMP Criterion C has been violated!\n");
        fprintf(OUT, "NOTE: rk67_mp sees that AMP Criterion C has been violated - the precision will be increased.\n");
        retVal = retVal_higher_prec_needed;
      }
      if (retVal)
      { // clear MP
        clear_mp(newTime);
        clear_point_mp(Y);
        for (i = 0; i < 10; i++)
          clear_point_mp(dX[i]);

        return retVal;
      }
    }

    if (k < 9)
    { // another function evaluation is needed
      // calculate dX[k] *= dT
      // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
      increase_size_vec_mp(Y, dX[k]->size);
      rows = Y->size = dX[k]->size;
      for (i = 0; i < rows; i++)
      {
        mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
        set_mp(&Y->coord[i], &P->point->coord[i]);
        for (j = 0; j <= k; j++)
        {
          mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
          add_mp(&Y->coord[i], &Y->coord[i], newTime);
        }
      }

      // find newTime = t + c[k+1]*dT
      mul_rmpf_mp(newTime, dT, c[k+1]);
      add_mp(newTime, newTime, P->time);

      // do the next - for k+1 - function evaluation
      eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
    }
    else
    { // calculate dX[k] = dT * dX[k]
      vec_mulcomp_mp(dX[k], dX[k], dT);
    }
  }

  // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k])
  // find newPoint = P->point + sum(k = 0..5, a_k * dX[k])
  increase_size_vec_mp(Y, dX[0]->size);
  increase_size_vec_mp(newPoint, dX[0]->size);
  rows = Y->size = newPoint->size = dX[0]->size;
  for (i = 0; i < rows; i++)
  {
    mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
    mul_rmpf_mp(&dX[0]->coord[i], &dX[0]->coord[i], a[0]);
    for (k = 1; k < 10; k++)
    {
      mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
      add_mp(&Y->coord[i], &Y->coord[i], newTime);

      mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
      add_mp(&dX[0]->coord[i], &dX[0]->coord[i], newTime);
    }
    add_mp(&newPoint->coord[i], &P->point->coord[i], &dX[0]->coord[i]);
  }
  infNormVec_mp2(err, Y);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 0, P->point);
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 0, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 0, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 0, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 0, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 0, e->Jp);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 0, dX[0]);
    fprintf(OUT, "newPoint = ");  printVec_mp(OUT, 0, newPoint);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 0, P->point);
      printf("H_of_X = ");  printPoint_mp(stdout, 0, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 0, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 0, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 0, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 0, e->Jp);
      printf("dX = ");  printVec_mp(stdout, 0, dX[0]);
      printf("newPoint = ");  printVec_mp(stdout, 0, newPoint);
    }
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 10; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

