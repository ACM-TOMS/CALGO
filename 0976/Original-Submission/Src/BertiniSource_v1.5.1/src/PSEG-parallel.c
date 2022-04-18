// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

// This file contains the functions needed to run the power series endgame
void change_size_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int num_new_points);

void refine_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int find_cycle_num_loop_d(PSEG_samples_struct_d *S, int start, int finish, int max_cycle_num, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int setup_PSEG_prediction_d(int setup_dX_Z, int cycle_num, PSEG_samples_struct_d *S, _comp_d *s, _point_d *dZ_ds, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int prediction_PSEG_samples_struct_d(vec_d approx, comp_d finalT, int cycle_num, int setup_dX_Z, PSEG_samples_struct_d *S, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int PSEG_samples_prediction_d(point_d approx, comp_d finalTime, int *cycle_num, PSEG_samples_struct_d *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

void change_size_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int num_new_points);

void refine_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int find_cycle_num_loop_mp(PSEG_samples_struct_mp *S, int start, int finish, int max_cycle_num, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int setup_PSEG_prediction_mp(int setup_dX_Z, int cycle_num, PSEG_samples_struct_mp *S, _comp_mp *s, _point_mp *dZ_ds, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int prediction_PSEG_samples_struct_mp(vec_mp approx, comp_mp finalT, int cycle_num, int setup_dX_Z, PSEG_samples_struct_mp *S, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int PSEG_samples_prediction_mp(point_mp approx, comp_mp finalTime_mp, int *cycle_num, PSEG_samples_struct_mp *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));


int PSEG_d(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does power series endgame in fixed double precision    *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_d lastSample;

  init_point_data_d(&lastSample, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  retVal = PSEG_d2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  clear_point_data_d(&lastSample);

  return retVal;
}

int PSEG_rank_d(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Power series endgame using double precision and   *
* uses Cauchy endgame for rank determination                    *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_d finalTime;
  point_data_d lastSample;

  init_point_data_d(&lastSample, 0);

  // setup finalTime
  set_double_d(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // perform the Power Series endgame - returns the last sample point so that we can use it for rank determination
  retVal = PSEG_d2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_d(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &lastSample, cycle_num, samples_per_loop, T, OUT, ED, eval_func);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx
    point_cp_d(last_approx, Final->point);
    for (i = 0; i < last_approx->size; i++)
    {
      get_comp_rand_d(finalTime);
      mul_rdouble_d(finalTime, finalTime, T->final_tolerance);
      add_d(&last_approx->coord[i], &last_approx->coord[i], finalTime);
    }

    // assume failure since we could not converge
    if (rankType == 0)
      *rankDef = -1;
    else
      *corank = 1;
  }

  clear_point_data_d(&lastSample);

  return retVal;
}

int PSEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *lastSample, int *cycle_num, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does power series endgame in fixed double precision &  *
* returns the last sample point                                 *
\***************************************************************/
{
  int i, retVal, numSamples = T->num_PSEG_sample_points;
  comp_d endTime;
  PSEG_samples_struct_d PSEG_samples;

  // initialize PSEG_samples
  init_PSEG_samples_struct_d(&PSEG_samples, numSamples);
  // copy Start to samples[0]
  point_data_cp_d(&PSEG_samples.samples[0], Start);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT 
  set_double_d(endTime, T->endgameBoundary, 0.0);
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // refine the start point and then track
  refine_PSEG_samples_struct_d(&PSEG_samples, 0, 1, T, OUT, ED, eval_func);
  retVal = track_d(&PSEG_samples.samples[0], &PSEG_samples.samples[0], endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < PSEG_samples.samples[0].point->size; i++)
    fprintf(midOUT, "%.15e %.15e\n", PSEG_samples.samples[0].point->coord[i].r, PSEG_samples.samples[0].point->coord[i].i);
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Power series endgame never started!\n");

    // copy the last sample to Final
    point_data_cp_d(Final, &PSEG_samples.samples[0]);
    Final->cycle_num = *cycle_num = 0;
    T->t_val_at_latest_sample_point = PSEG_samples.samples[0].time->r;
  }
  else
  { // success - so run the actual endgame now that we are at the endgame boundary

    // we have 1 sample point
    PSEG_samples.num_samples = 1;

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // run the actual PSEG endgame
    retVal = PSEG_struct_d(Final, last_approx, &PSEG_samples, T, OUT, ED, eval_func, find_dehom);

    // setup lastSample & cycle_num
    *cycle_num = Final->cycle_num;
    point_data_cp_d(lastSample, &PSEG_samples.samples[PSEG_samples.num_samples - 1]);
  }

  // free the memory
  clear_PSEG_samples_struct_d(&PSEG_samples);

  return retVal;
}

void cubicInterp_d(point_d Res, _comp_d *T, _point_d *Y, _point_d *dHdT, comp_d T0)
/***************************************************************\
* USAGE: does a cubic interpolation - i.e. uses Y[0] & Y[1] with*
* T[0] & T[1] and dHdT[0] & dHdT[1] to do a cubic approximation *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = Y[0].size;
  comp_d a, b, c, d, h, recip_h_sqr, diff_0_0, diff_0_1, sqr0;

  // set the size
  change_size_point_d(Res, size);
  Res->size = size;

  // find h = T[1] - T[0]
  sub_d(h, &T[1], &T[0]);

  // make sure |h| > 0
  if (d_abs_d(h) == 0)
  {
    printf("ERROR: The times in cubicInterp_d are equal!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find 1 / h^2
  sqr_d(recip_h_sqr, h);
  recip_d(recip_h_sqr, recip_h_sqr);

  // find diff_0_0
  sub_d(diff_0_0, T0, &T[0]);
  // find diff_0_1
  sub_d(diff_0_1, T0, &T[1]);
  // find sqr0
  sqr_d(sqr0, diff_0_0);

  // find each coordinate
  for (i = 0; i < size; i++)
  { // find a = Y[0]
    set_d(a, &Y[0].coord[i]);
    // find b = dH/dt[0]
    set_d(b, &dHdT[0].coord[i]);

    // find c = (Y[1] - a - b*h) / h^2
    mul_d(c, b, h);
    sub_d(c, &Y[1].coord[i], c);
    sub_d(c, c, a);
    mul_d(c, c, recip_h_sqr);

    // find d = (dHdT[1] - b - 2*c*h) / h^2
    mul_d(d, c, h);
    mul_rdouble_d(d, d, 2.0);
    sub_d(d, &dHdT[1].coord[i], d);
    sub_d(d, d, b);
    mul_d(d, d, recip_h_sqr);

    // find Res = a + b(T0 - T[0]) + c(T0 - T[0])^2 + d(T0 - T[0])^2(T0 - T[1])
    mul_d(&Res->coord[i], sqr0, diff_0_1);
    mul_d(&Res->coord[i], &Res->coord[i], d);
    sum_mul_d(&Res->coord[i], c, sqr0);
    sum_mul_d(&Res->coord[i], b, diff_0_0);
    add_d(&Res->coord[i], &Res->coord[i], a);
  }

  return;
}

void cubicInterp_mp(point_mp Res, _comp_mp *T, _point_mp *Y, _point_mp *dHdT, comp_mp T0)
/***************************************************************\
* USAGE: does a cubic interpolation - i.e. uses Y[0] & Y[1] with*
* T[0] & T[1] and dHdT[0] & dHdT[1] to do a cubic approximation *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = Y[0].size;
  comp_mp a, b, c, d, h, recip_h_sqr, diff_0_0, diff_0_1, sqr0;

  // initialize the MP
  init_mp(a); init_mp(b); init_mp(c); init_mp(d); init_mp(h);
  init_mp(recip_h_sqr); init_mp(diff_0_0); init_mp(diff_0_1); init_mp(sqr0);

  // set the size
  change_size_point_mp(Res, size);
  Res->size = size;

  // find h = T[1] - T[0]
  sub_mp(h, &T[1], &T[0]);

  // make sure |h| > 0
  if (d_abs_mp(h) == 0)
  {
    printf("ERROR: The times in cubicInterp_mp are equal!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find 1 / h^2
  sqr_mp(recip_h_sqr, h);
  recip_mp(recip_h_sqr, recip_h_sqr);

  // find diff_0_0
  sub_mp(diff_0_0, T0, &T[0]);
  // find diff_0_1
  sub_mp(diff_0_1, T0, &T[1]);
  // find sqr0
  sqr_mp(sqr0, diff_0_0);

  // find each coordinate
  for (i = 0; i < size; i++)
  {
    // find a = Y[0]
    set_mp(a, &Y[0].coord[i]);
    // find b = dH/dt[0]
    set_mp(b, &dHdT[0].coord[i]);

    // find c = (Y[1] - a - b*h) / h^2
    mul_mp(c, b, h);
    sub_mp(c, &Y[1].coord[i], c);
    sub_mp(c, c, a);
    mul_mp(c, c, recip_h_sqr);

    // find d = (dHdT[1] - b - 2*c*h) / h^2
    mul_mp(d, c, h);
    mpf_mul_ui(d->r, d->r, 2);
    mpf_mul_ui(d->i, d->i, 2);
    sub_mp(d, &dHdT[1].coord[i], d);
    sub_mp(d, d, b);
    mul_mp(d, d, recip_h_sqr);

    // find Res = a + b(T0 - T[0]) + c(T0 - T[0])^2 + d(T0 - T[0])^2(T0 - T[1])
    mul_mp(&Res->coord[i], sqr0, diff_0_1);
    mul_mp(&Res->coord[i], &Res->coord[i], d);
    sum_mul_mp(&Res->coord[i], c, sqr0);
    sum_mul_mp(&Res->coord[i], b, diff_0_0);
    add_mp(&Res->coord[i], &Res->coord[i], a);
  }

  // release the MP
  clear_mp(a); clear_mp(b); clear_mp(c); clear_mp(d); clear_mp(h);
  clear_mp(recip_h_sqr); clear_mp(diff_0_0); clear_mp(diff_0_1); clear_mp(sqr0);

  return;
}

void init_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int num_points)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initiailizes S                                         *
\***************************************************************/
{
  int i;

  S->samples = (point_data_d *)bmalloc(num_points * sizeof(point_data_d));
  S->dX = (vec_d *)bmalloc(num_points * sizeof(vec_d));
  S->Z_rev = (_point_d *)bmalloc(num_points * sizeof(_point_d));
  for (i = 0; i < num_points; i++)
  {
    init_point_data_d(&S->samples[i], 0);
    init_vec_d(S->dX[i], 0);
    init_point_d(&S->Z_rev[i], 0);
    S->samples[i].point->size = S->dX[i]->size = S->Z_rev[i].size = 0;
    S->samples[i].cycle_num = 0;
  }
  S->num_samples = 0;
  S->mem_count = num_points;
 
  return;
}

void change_size_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int num_new_points)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: changes the size of S                                  *
\***************************************************************/
{
  int i;

  if (num_new_points != S->mem_count)
  {
    S->samples = (point_data_d *)brealloc(S->samples, num_new_points * sizeof(point_data_d));
    S->dX = (vec_d *)brealloc(S->dX, num_new_points * sizeof(vec_d));
    S->Z_rev = (_point_d *)brealloc(S->Z_rev, num_new_points * sizeof(_point_d));
    for (i = S->num_samples; i < num_new_points; i++)
    {
      init_point_data_d(&S->samples[i], 0);
      init_vec_d(S->dX[i], 0);
      init_point_d(&S->Z_rev[i], 0);
      S->samples[i].point->size = S->dX[i]->size = S->Z_rev[i].size = 0;
      S->samples[i].cycle_num = 0;
    }
    S->mem_count = num_new_points;
  }

  return;
}

void clear_PSEG_samples_struct_d(PSEG_samples_struct_d *S)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears S                                               *
\***************************************************************/
{
  int i;

  for (i = S->mem_count - 1; i >= 0; i--)
  {
    clear_point_data_d(&S->samples[i]);
    clear_vec_d(S->dX[i]);
    clear_point_d(&S->Z_rev[i]);
  }
  free(S->samples);
  free(S->dX);
  free(S->Z_rev);
  S->num_samples = 0;
  S->mem_count = 0;

  return;
}

int setup_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int start, int finish, int MPType, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup dX                                               *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;

  eval_struct_d e;
  init_eval_struct_d(e, 0, 0, 0);

  for (i = start; i < finish; i++)
  { // evaluate the function
    eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples[i].point, S->samples[i].time, ED, eval_func);

    // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(e.funcVals, e.Jp->rows);
    rows = e.funcVals->size = e.Jp->rows;
    cols = e.Jp->cols;
    for (k = 0; k < rows; k++)
    {
      set_zero_d(&e.funcVals->coord[k]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&e.funcVals->coord[k], &e.Jp->entry[k][j], &e.parDer->coord[j]);
      }
      neg_d(&e.funcVals->coord[k], &e.funcVals->coord[k]);
    }

    // find dX = (dH/dX)^(-1) * (-dH/dt)
    retVal = matrixSolve_d(S->dX[i], e.Jv, e.funcVals);
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve has failed.\n");

      clear_eval_struct_d(e);
      // return error
      if (MPType == 2)
        return retVal_higher_prec_needed;
      else
        return retVal;
    }
  }

  clear_eval_struct_d(e);

  return retVal;
}

void init_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int num_points, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initiailizes S                                         *
\***************************************************************/
{
  int i;

  S->samples = (point_data_mp *)bmalloc(num_points * sizeof(point_data_mp));
  S->dX = (vec_mp *)bmalloc(num_points * sizeof(vec_mp));
  S->Z_rev = (_point_mp *)bmalloc(num_points * sizeof(_point_mp));
  for (i = 0; i < num_points; i++)
  {
    init_point_data_mp2(&S->samples[i], 0, prec);
    init_vec_mp2(S->dX[i], 0, prec);
    init_point_mp2(&S->Z_rev[i], 0, prec);
    S->samples[i].point->size = S->dX[i]->size = S->Z_rev[i].size = 0;
    S->samples[i].cycle_num = 0; 
  }
  S->curr_prec = prec;
  S->num_samples = 0;
  S->mem_count = num_points;

  return;
}

void change_size_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int num_new_points)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initiailizes S                                         *
\***************************************************************/
{
  int i, old_size = S->num_samples;

  if (S->mem_count != num_new_points)
  {
    point_data_mp *old_samples = S->samples; // store the old samples
    vec_mp *old_dX = S->dX;
    _point_mp *old_Z = S->Z_rev;

    // allocate new memory
    S->samples = NULL;
    S->dX = NULL;
    S->Z_rev = NULL;
    S->samples = (point_data_mp *)bmalloc(num_new_points * sizeof(point_data_mp));
    S->dX = (vec_mp *)bmalloc(num_new_points * sizeof(vec_mp));
    S->Z_rev = (_point_mp *)bmalloc(num_new_points * sizeof(_point_mp));
    // initialize new memory and copy over old memory
    for (i = num_new_points - 1; i >= 0; i--)
    {
      init_point_data_mp2(&S->samples[i], 0, S->curr_prec);
      init_vec_mp2(S->dX[i], 0, S->curr_prec);
      init_point_mp2(&S->Z_rev[i], 0, S->curr_prec);

      if (i < old_size)
      { // copy over the old data - Z does not need copied over
        point_data_cp_mp(&S->samples[i], &old_samples[i]);
        point_cp_mp(S->dX[i], old_dX[i]);
      } 
    }
    S->mem_count = num_new_points;

    // clear old memory - this is outside the other loop since num_new_points could be < old_size
    for (i = old_size - 1; i >= 0; i--)
    {
      clear_point_data_mp(&old_samples[i]);
      clear_vec_mp(old_dX[i]);
      clear_point_mp(&old_Z[i]);
    }
    free(old_samples);
    free(old_dX);
    free(old_Z);
  }

  return;
}

void clear_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears S                                               *
\***************************************************************/
{
  int i;

  for (i = S->mem_count - 1; i >= 0; i--)
  {
    clear_point_data_mp(&S->samples[i]);
    clear_vec_mp(S->dX[i]);
    clear_point_mp(&S->Z_rev[i]);
  } 
  free(S->samples);
  free(S->dX);
  free(S->Z_rev);
  S->curr_prec = 0;
  S->num_samples = 0;
  S->mem_count = 0;           

  return;
}

int setup_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int start, int finish, int MPType, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup dX                                               *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0;

  if (start < finish)
  { // only initialize if we actually use it
    eval_struct_mp e;
    init_eval_struct_mp(e, 0, 0, 0);

    for (i = start; i < finish; i++)
    { // evaluate the function
      eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples[i].point, S->samples[i].time, ED, eval_func);

      // find Y = -dH/dt = -dH/dp * dp/dt
      increase_size_vec_mp(e.funcVals, e.Jp->rows);
      rows = e.funcVals->size = e.Jp->rows;
      cols = e.Jp->cols;
      for (k = 0; k < rows; k++)
      {
        set_zero_mp(&e.funcVals->coord[k]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_mp(&e.funcVals->coord[k], &e.Jp->entry[k][j], &e.parDer->coord[j]);
        }
        neg_mp(&e.funcVals->coord[k], &e.funcVals->coord[k]);
      }

      // find dX = (dH/dX)^(-1) * (-dH/dt)
      retVal = matrixSolve_mp(S->dX[i], e.Jv, e.funcVals);
      if (retVal)
      {
        fprintf(OUT, "NOTE: matrixSolve has failed.\n");
        // clear e
        clear_eval_struct_mp(e);

        // return error
        if (MPType == 2)
          return retVal_higher_prec_needed;
        else
          return retVal;
      }
    }

    // clear if it was initialized
    clear_eval_struct_mp(e);
  }

  return retVal;
}

double compare_approximations_d(point_d approx0, point_d approx1)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compares the approximations and returns the error      *
\***************************************************************/
{
  int i, size;
  double error, tempD;
  comp_d tempComp;

  // find the error
  size = approx1->size;
  error = 0;
  for (i = 0; i < size; i++)
  {
    sub_d(tempComp, &approx0->coord[i], &approx1->coord[i]);
    tempD = tempComp->r * tempComp->r + tempComp->i * tempComp->i;
    if (tempD > error)
      error = tempD;
  }
  error = sqrt(error);

  return error;
}

int PSEG_samples_prediction_d(point_d approx, comp_d finalTime, int *cycle_num, PSEG_samples_struct_d *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sets up and finds the approximation based on start &   *
* finish                                                        *
\***************************************************************/
{
  int retVal = 0;

  // error checking
  if (T->num_PSEG_sample_points < 2)
  {
    printf("ERROR: You must use atleast 2 samples points for a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (finish - start != T->num_PSEG_sample_points + 1)
  {
    printf("ERROR: The number of samples points passed to the predictor does not match the number expected!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (start < 0 || finish > PSEG_samples->num_samples)
  {
    printf("ERROR: Incorrect settings for power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // refine the sample points
  refine_PSEG_samples_struct_d(PSEG_samples, start, finish, T, OUT, ED, eval_func);

  // find an approximation based on the last 3 samples
  *cycle_num = find_cycle_num_d(PSEG_samples->samples, finish - 3, T->power_series_sample_factor);
  // use this approximation as a way to determine where to stop the exhaustive search for finding the best cycle number
  *cycle_num = find_cycle_num_loop_d(PSEG_samples, start, finish, MAX(*cycle_num * 5, T->cycle_num_max), OUT, ED, eval_func);

  if (*cycle_num > 0)
  { // the structures are setup properly - find the prediction

    // store the cycle number to [start]
    PSEG_samples->samples[start].cycle_num = *cycle_num;

    // find the approximation
    retVal = prediction_PSEG_samples_struct_d(approx, finalTime, *cycle_num, 0, PSEG_samples, start, finish, OUT, ED, eval_func);
  }
  else
  { // there was an error
    retVal = *cycle_num;
  }

  return retVal;
}

void refine_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: refines the samples points from start to < finish      *
\***************************************************************/
{
  int i;
  eval_struct_d e;
  init_eval_struct_d(e, 0, 0, 0);

  for (i = start; i < finish; i++)
  { // refine
    refine_d_basic(&S->samples[i], T, OUT, &e, ED, eval_func);
  }

  clear_eval_struct_d(e);

  return;
}

int find_cycle_num_d(point_data_d *Pts, int start, double sample_factor)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximates the cycle number using the sample factor  *
* and the 3 sample points {start, start+1, start_2}             *
\***************************************************************/
{
  int i, size, retVal;
  double c;
  comp_d temp1, temp2, sum1, sum2, rand;

  // find the size
  size = Pts[start].point->size;

static vec_d randVec;
static int init = 0;
if (!init) {init_vec_d(randVec,size); make_vec_random_d(randVec,size); }

  // find a random sum of the subtractions
  set_zero_d(sum1);
  set_zero_d(sum2);
  for (i = 0; i < size; i++)
  { // temp1 = [start+1] - [start]
    sub_d(temp1, &Pts[start+1].point->coord[i], &Pts[start].point->coord[i]);
    // temp2 = [start+2] - [start+1]
    sub_d(temp2, &Pts[start+2].point->coord[i], &Pts[start+1].point->coord[i]);

    // find random number
//    get_comp_rand_d(rand);
set_d(rand,&randVec->coord[i]);

    // add to sum1 & sum2
    sum_mul_d(sum1, rand, temp1);
    sum_mul_d(sum2, rand, temp2);
  }

  // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1
  if (d_abs_d(sum1) > 0 && d_abs_d(sum2) > 0)
  { // sum2 / sum1
    div_d(temp1, sum2, sum1);

    // find the norm of sum2 / sum1
    c = d_abs_d(temp1);

    // find the approximate cycle number
    c = log(sample_factor) / log(c);

    // find the cycle number - a positive integer
    if (c < 1)
      retVal = 1;
    else
    { // round c
      retVal = (int) floor(c + 0.5);
    }
  }
  else
  { // cycle number is 1
    retVal = 1;
  }

  // double check retVal >= 1
  if (retVal < 1)
    retVal = 1;

  return retVal;
}

int find_cycle_num_loop_d(PSEG_samples_struct_d *S, int start, int finish, int max_cycle_num, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - error, otherwise best_cycle_num           *
* NOTES: finds the number between 1 and max_cycle_num that best *
* uses {start,..,finish-1} to approximate finish                *
\***************************************************************/
{
  int c, j, size, retVal, best_c = 1, n = finish - start - 1;
  double min_error = 1e300, error = 0, tempD;
  point_d approx;

  init_point_d(approx, 0);

  // error checking
  if (start < 0 || finish > S->num_samples || max_cycle_num < 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup s & dZ_ds
  _comp_d *s = (_comp_d *)bmalloc((n + 1) * sizeof(_comp_d));
  _point_d *dZ_ds = (_point_d *)bmalloc((n + 1) * sizeof(_point_d));

  for (j = 0; j <= n; j++)
    init_point_d(&dZ_ds[j], 0);

  for (c = 1; c <= max_cycle_num; c++)
  { // setup the structures to do an interpolation
    retVal = setup_PSEG_prediction_d((c == 1), c, S, s, dZ_ds, start, finish, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // clear structurs
      clear_point_d(approx);
      for (j = n; j >= 0; j--)
        clear_point_d(&dZ_ds[j]);
      free(dZ_ds);
      free(s);
      // return error
      return -1;
    }

    // interpolate from 1, .. n
    if (n == 2)
      cubicInterp_d(approx, &s[1], &S->Z_rev[1], &dZ_ds[1], &s[0]);
    else
      hermiteInterpCW_d(approx, &s[1], &S->Z_rev[1], &dZ_ds[1], &s[0], n);

    // find the error
    error = 0;
    size = approx->size;
    for (j = 0; j < size; j++)
    {
      sub_d(&approx->coord[j], &approx->coord[j], &S->samples[finish-1].point->coord[j]);
      tempD = approx->coord[j].r * approx->coord[j].r + approx->coord[j].i * approx->coord[j].i;
      if (tempD > error)
        error = tempD;
    }

    // check to see if this is the minimum error
    if (min_error > error)
    {
      min_error = error;
      best_c = c;
    }
  }

  // clear structures
  for (j = n; j >= 0; j--)
    clear_point_d(&dZ_ds[j]);
  free(dZ_ds);
  free(s);
  clear_point_d(approx);

  return best_c;
}

int setup_PSEG_prediction_d(int setup_dX_Z, int cycle_num, PSEG_samples_struct_d *S, _comp_d *s, _point_d *dZ_ds, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE: setup the structures needed to do a prediction         *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - everything is setup, otherwise, some error *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0, i, j, k, rows, cols, size, n = finish - start - 1;
  double tempD;
  comp_d factor, tempComp, c_times_t;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (setup_dX_Z)
  { // setup dX & Z_rev
    eval_struct_d e;
    init_eval_struct_d(e, 0, 0, 0);

    for (i = start; i < finish; i++)
    { // evaluate the function
      eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples[i].point, S->samples[i].time, ED, eval_func);

      // find Y = -dH/dt = -dH/dp * dp/dt
      increase_size_vec_d(e.funcVals, e.Jp->rows);
      rows = e.funcVals->size = e.Jp->rows;
      cols = e.Jp->cols;
      for (k = 0; k < rows; k++)
      {
        set_zero_d(&e.funcVals->coord[k]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_d(&e.funcVals->coord[k], &e.Jp->entry[k][j], &e.parDer->coord[j]);
        }
        neg_d(&e.funcVals->coord[k], &e.funcVals->coord[k]);
      }

      // find dX = (dH/dX)^(-1) * (-dH/dt)
      retVal = matrixSolve_d(S->dX[i], e.Jv, e.funcVals);
      if (retVal)
      {
        fprintf(OUT, "NOTE: matrixSolve has failed.\n");
        clear_eval_struct_d(e);
        // return error
        return retVal;
      }

      // setup Z_rev
      k = finish - i - 1;
      increase_size_vec_d(&S->Z_rev[k], S->samples[i].point->size);
      size = S->Z_rev[k].size = S->samples[i].point->size;
      for (j = 0; j < size; j++)
      {
        set_d(&S->Z_rev[k].coord[j], &S->samples[i].point->coord[j]);
      }
    }

    clear_eval_struct_d(e);
  }

  // setup dZ_ds & s
  tempD = 1 / (double) cycle_num;
  // since normalize time, set s[n] == 1
  pow_rdouble_d(factor, S->samples[start].time, tempD);
  recip_d(factor, factor);
  set_double_d(&s[n], 1, 0);
  // find dZ_ds[n]
  mul_rdouble_d(c_times_t, S->samples[start].time, cycle_num);
  vec_mulcomp_d(&dZ_ds[n], S->dX[start], c_times_t);

  // loop over i = n - 1, .., 0 to find the dZ_ds[i] & s[i]
  for (i = n - 1; i >= 0; i--)
  { // find j
    j = finish - i - 1;
    // s_d[i] = t^(1/c) * factor
    pow_rdouble_d(&s[i], S->samples[j].time, tempD);
    mul_d(&s[i], &s[i], factor);
    // dZ_ds[i]
    k = cycle_num - 1;
    pow_rdouble_d(tempComp, &s[i], k);
    mul_d(tempComp, tempComp, c_times_t);
    vec_mulcomp_d(&dZ_ds[i], S->dX[j], tempComp);
  }

  return retVal;
}

int prediction_PSEG_samples_struct_d(vec_d approx, comp_d finalT, int cycle_num, int setup_dX_Z, PSEG_samples_struct_d *S, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE: does a power series prediction based on samples to     *
* predict the value at finalT                                   *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - have an approx, otherwise, some error      *
* NOTES:                                                        *
\***************************************************************/
{
  int j, retVal = 0, n = finish - start - 1;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup s & dZ_ds
  _comp_d *s = (_comp_d *)bmalloc((n + 1) * sizeof(_comp_d));
  _point_d *dZ_ds = (_point_d *)bmalloc((n + 1) * sizeof(_point_d));

  for (j = 0; j <= n; j++)
    init_point_d(&dZ_ds[j], 0);

  // setup the structures
  retVal = setup_PSEG_prediction_d(setup_dX_Z, cycle_num, S, s, dZ_ds, start, finish, OUT, ED, eval_func);

  // check for error
  if (!retVal)
  { // interpolate
    if (n + 1 == 2)
      cubicInterp_d(approx, &s[0], &S->Z_rev[0], &dZ_ds[0], finalT);
    else
      hermiteInterpCW_d(approx, &s[0], &S->Z_rev[0], &dZ_ds[0], finalT, n + 1);
  }

  // free the memory
  for (j = n; j >= 0; j--)
    clear_point_d(&dZ_ds[j]);
  free(dZ_ds);
  free(s);

  return retVal;
}

int PSEG_struct_d(point_data_d *Final, point_d last_approx, PSEG_samples_struct_d *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual PSEG endgame in double precision       *
* storing the samples along the way to samples                  *
\***************************************************************/
{
  int i, dehom_prec, cycle_num = 0, retVal = 0, start, finish, numSamples = T->num_PSEG_sample_points;
  double approx_err, norm_last = 0, norm_new = 0;
  comp_d endTime, finalTime;
  point_d new_approx, dehom;

  // initialize new_approx & dehom
  init_point_d(new_approx, 0);
  init_point_d(dehom, 0);

  // initialize cycle number
  Final->cycle_num = 0;

  // make sure that there is atleast 1 good sample point
  if (PSEG_samples->num_samples <= 0)
  {
    printf("WARNING: The number of good sample points must be positive in PSEG_struct_d!!\n");
    // clear new_approx & dehom
    clear_point_d(new_approx);
    clear_point_d(dehom);

    return retVal_PSEG_failed;
  }

  // setup start and finish
  finish = PSEG_samples->num_samples - 1;
  start = finish - numSamples;
  if (start < 0)
    start = 0;

  // setup finalTime
  set_double_d(finalTime, T->targetT, 0.0);

  // make sure there are enough locations to find the initial set of sample points
  if (PSEG_samples->mem_count <= numSamples)
  { // increase the memory
    change_size_PSEG_samples_struct_d(PSEG_samples, numSamples + 1);
  }

  // find the initial set of samples points to create the first approximation
  for (i = finish + 1; i <= numSamples; i++)
  { // find the end time for the next samples
    mul_rdouble_d(endTime, PSEG_samples->samples[i-1].time, T->power_series_sample_factor);

    // refine the sample point
    refine_PSEG_samples_struct_d(PSEG_samples, i-1, i, T, OUT, ED, eval_func);
    // track to the next sample point
    retVal = track_d(&PSEG_samples->samples[i], &PSEG_samples->samples[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      point_data_cp_d(Final, &PSEG_samples->samples[i]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples->samples[i].time->r;

      // clear new_approx & dehom
      clear_point_d(new_approx);
      clear_point_d(dehom);

      // return error
      return retVal;
    }
    else
    { // update the number of good samples
      PSEG_samples->num_samples++;
    }
  }
  // so, we have enough sample points for the initial approximation

  // find the start & end
  finish = PSEG_samples->num_samples;
  start = finish - numSamples - 1;

  // find the first approximation
  retVal = PSEG_samples_prediction_d(last_approx, finalTime, &cycle_num, PSEG_samples, start, finish, T, OUT, ED, eval_func);
  // find the dehom norm of last_approx, if needed

  // check for error
  if (retVal)
  { // total failure since we cannot find the first approximation - setup Final and return error
    point_data_cp_d(Final, &PSEG_samples->samples[finish-1]);
    Final->cycle_num = 0;
    T->t_val_at_latest_sample_point = PSEG_samples->samples[finish-1].time->r;

    // clear new_approx & dehom
    clear_point_d(new_approx);
    clear_point_d(dehom);

    // return error
    return retVal;
  }
  else if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom, NULL, &dehom_prec, last_approx, NULL, 52, ED, NULL);
    norm_last = infNormVec_d(dehom);
  }

  // now that we have the first approximation, we can loop until we have convergence
  while (1)
  { // find the index of the next sample and the current precision
    i = PSEG_samples->num_samples;
    // setup endTime
    mul_rdouble_d(endTime, PSEG_samples->samples[i-1].time, T->power_series_sample_factor);

    // make sure there is enough memory
    if (PSEG_samples->mem_count <= PSEG_samples->num_samples)
    { // increase the memory
      change_size_PSEG_samples_struct_d(PSEG_samples, 2 * PSEG_samples->mem_count);
    }

    // refine the sample point
    refine_PSEG_samples_struct_d(PSEG_samples, i-1, i, T, OUT, ED, eval_func);
    // track to the next sample point
    retVal = track_d(&PSEG_samples->samples[i], &PSEG_samples->samples[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      point_data_cp_d(Final, &PSEG_samples->samples[i]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples->samples[i].time->r;

      // break out of the loop to return error
      break;
    }
 
    // update the number of good samples
    PSEG_samples->num_samples++;

    // find the start & end
    finish = PSEG_samples->num_samples;
    start = finish - numSamples - 1;

    // find the next approximation
    retVal = PSEG_samples_prediction_d(new_approx, finalTime, &cycle_num, PSEG_samples, start, finish, T, OUT, ED, eval_func);

    // check for error
    if (retVal)
    { // total failure since we cannot find the first approximation - setup Final and return error
      point_data_cp_d(Final, &PSEG_samples->samples[finish-1]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples->samples[finish-1].time->r;

      // break out of the loop to return error
      break;
    }
    else if (T->securityLevel <= 0)
    { // find dehom & its norm
      find_dehom(dehom, NULL, &dehom_prec, new_approx, NULL, 52, ED, NULL);
      norm_new = infNormVec_d(dehom);
    }
    
    // compare the approximations
    T->error_at_latest_sample_point = approx_err = compare_approximations_d(new_approx, last_approx);

    // check for convergence
    if (approx_err < T->final_tolerance)
    { // accept the approximation

      // setup Final
      point_cp_d(Final->point, new_approx);
      set_d(Final->time, finalTime);
      Final->cycle_num = cycle_num;

      T->t_val_at_latest_sample_point = PSEG_samples->samples[finish-1].time->r;

      // break out of loop
      retVal = 0;
      break;
    }
    else if (PSEG_samples->samples[finish-1].time->r < T->minTrackT)
    { // too close to t = 0 but we do not have the correct tolerance - so we exit

      // setup Final
      point_cp_d(Final->point, new_approx);
      set_d(Final->time, finalTime);
      Final->cycle_num = cycle_num;

      T->t_val_at_latest_sample_point = PSEG_samples->samples[finish-1].time->r;

      // break out of loop
      retVal = retVal_EG_failed_to_converge;
      break;
    }
    else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
    { // we are too large 
      point_data_cp_d(Final, &PSEG_samples->samples[finish-1]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples->samples[finish-1].time->r;

      // break out of the loop to return error
      retVal = retVal_security_max;
      break;
    }
    else
    { // simply copy new_approx to last_approx and try again
      point_cp_d(last_approx, new_approx);
      norm_last = norm_new;
    }
  }

  // clear new_approx & dehom
  clear_point_d(new_approx);
  clear_point_d(dehom);

  return retVal;
}

///////////// MAIN MP VERSIONS //////////////////

int PSEG_mp(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does power series endgame in fixed multi precision     *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_mp lastSample;

  init_point_data_mp(&lastSample, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  retVal = PSEG_mp2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  // clear
  clear_point_data_mp(&lastSample);

  return retVal;
}

int PSEG_rank_mp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Power series endgame using fixed multi precision  *
* and uses Cauchy endgame for rank determination                *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_mp finalTime;
  point_data_mp lastSample;

  init_mp(finalTime);
  init_point_data_mp(&lastSample, 0);

  // setup finalTime
  set_double_mp(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // perform the Power Series endgame - returns that last sample so that it can be used for rank determination
  retVal = PSEG_mp2(pathNum, Final, last_approx, &lastSample, &cycle_num, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_mp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &lastSample, cycle_num, samples_per_loop, T, OUT, ED, eval_func);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx 
    point_cp_mp(last_approx, Final->point);
    for (i = 0; i < last_approx->size; i++)
    {
      get_comp_rand_mp(finalTime);
      mul_rdouble_mp(finalTime, finalTime, T->final_tolerance);
      add_mp(&last_approx->coord[i], &last_approx->coord[i], finalTime);
    }
  }
 
  // clear
  clear_mp(finalTime);
  clear_point_data_mp(&lastSample);

  return retVal;
} 

int PSEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *lastSample, int *cycle_num, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does power series endgame in fixed multi precision &   *
* returns the last sample point                                 *
\***************************************************************/
{
  int i, retVal, numSamples = T->num_PSEG_sample_points;
  comp_mp endTime;
  PSEG_samples_struct_mp PSEG_samples;

  // initialize endTime
  init_mp(endTime);

  // initialize PSEG_samples
  init_PSEG_samples_struct_mp(&PSEG_samples, numSamples, T->Precision);

  // copy Start to samples[0]
  point_data_cp_mp(&PSEG_samples.samples[0], Start);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_mp(endTime, T->endgameBoundary, 0.0);
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // refine the start point and then track
  refine_PSEG_samples_struct_mp(&PSEG_samples, 0, 1, T, OUT, ED, eval_func);
  retVal = track_mp(&PSEG_samples.samples[0], &PSEG_samples.samples[0], endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < PSEG_samples.samples[0].point->size; i++)
  {
    print_mp(midOUT, 0, &PSEG_samples.samples[0].point->coord[i]);
    fprintf(midOUT, "\n");
  }
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Power series endgame never started!\n");

    // copy the last sample to Final
    point_data_cp_mp(Final, &PSEG_samples.samples[0]);
    Final->cycle_num = *cycle_num = 0;
    T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples.samples[0].time->r);
  }
  else
  { // success - so run the actual endgame now that we are at the endgame boundary

    // we have 1 sample point
    PSEG_samples.num_samples = 1;

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // run the actual PSEG endgame
    retVal = PSEG_struct_mp(Final, last_approx, &PSEG_samples, T, OUT, ED, eval_func, find_dehom);

    // setup lastSample & cycle_num
    *cycle_num = Final->cycle_num;
    point_data_cp_mp(lastSample, &PSEG_samples.samples[PSEG_samples.num_samples - 1]);
  }

  // free the memory
  clear_mp(endTime);
  clear_PSEG_samples_struct_mp(&PSEG_samples);

  return retVal;
}

double compare_approximations_mp(point_mp approx0, point_mp approx1)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compares the approximations and returns the error      *
\***************************************************************/
{
  int i, size;
  double error_d;
  mpf_t error_mp;
  comp_mp tempComp;

  mpf_init(error_mp);
  init_mp(tempComp);

  mpf_set_ui(error_mp, 0);

  // find the error
  size = approx1->size;
  for (i = 0; i < size; i++)
  {
    sub_mp(tempComp, &approx0->coord[i], &approx1->coord[i]);
    mpf_mul(tempComp->r, tempComp->r, tempComp->r);
    mpf_mul(tempComp->i, tempComp->i, tempComp->i);
    mpf_add(tempComp->r, tempComp->r, tempComp->i);
    if (mpf_cmp(tempComp->r, error_mp) > 0)
      mpf_set(error_mp, tempComp->r);
  }
  mpf_sqrt(error_mp, error_mp);
  // set error
  error_d = mpf_get_d(error_mp);

  // clear MP
  mpf_clear(error_mp);
  clear_mp(tempComp);

  return error_d;
}

int PSEG_samples_prediction_mp(point_mp approx, comp_mp finalTime, int *cycle_num, PSEG_samples_struct_mp *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sets up and finds the approximation based on start &   *
* finish                                                        *
\***************************************************************/
{
  int retVal = 0;

  // error checking
  if (T->num_PSEG_sample_points < 2)
  {
    printf("ERROR: You must use atleast 2 samples points for a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (finish - start != T->num_PSEG_sample_points + 1)
  {
    printf("ERROR: The number of samples points passed to the predictor does not match the number expected!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (start < 0 || finish > PSEG_samples->num_samples)
  {
    printf("ERROR: Incorrect settings for power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // refine the sample points
  refine_PSEG_samples_struct_mp(PSEG_samples, start, finish, T, OUT, ED, eval_func);

  // find an approximation based on the last 3 samples
  *cycle_num = find_cycle_num_mp(PSEG_samples->samples, finish - 3, T->power_series_sample_factor);
  // use this approximation as a way to determine where to stop the exhaustive search for finding the best cycle number
  *cycle_num = find_cycle_num_loop_mp(PSEG_samples, start, finish, MAX(*cycle_num * 5, T->cycle_num_max), OUT, ED, eval_func);

  if (*cycle_num > 0)
  { // the structures are setup properly - find the prediction

    // store the cycle number to [start]
    PSEG_samples->samples[start].cycle_num = *cycle_num;

    // find the approximation
    retVal = prediction_PSEG_samples_struct_mp(approx, finalTime, *cycle_num, 0, PSEG_samples, start, finish, OUT, ED, eval_func);
  }
  else
  { // there was an error
    retVal = *cycle_num;
  }

  return retVal;
}

void refine_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: refines the samples points from start to < finish      *
\***************************************************************/
{
  int i;
  eval_struct_mp e;
  init_eval_struct_mp(e, 0, 0, 0);

  for (i = start; i < finish; i++)
  { // refine
    refine_mp_basic(&S->samples[i], T, OUT, &e, ED, eval_func);
  }

  clear_eval_struct_mp(e);

  return;
}

int find_cycle_num_mp(point_data_mp *Pts, int start, double sample_factor)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximates the cycle number using the sample factor  *
* and the 3 sample points {start, start+1, start_2}             *
\***************************************************************/
{
  int i, size, retVal;
  double c;
  comp_mp temp1, temp2, sum1, sum2, rand;

  // initialize
  init_mp(temp1); init_mp(temp2); init_mp(sum1); init_mp(sum2); init_mp(rand);

  // find the size
  size = Pts[start].point->size;

  // find a random sum of the subtractions
  set_zero_mp(sum1);
  set_zero_mp(sum2);
  for (i = 0; i < size; i++)
  { // temp1 = [start+1] - [start]
    sub_mp(temp1, &Pts[start+1].point->coord[i], &Pts[start].point->coord[i]);
    // temp2 = [start+2] - [start+1]
    sub_mp(temp2, &Pts[start+2].point->coord[i], &Pts[start+1].point->coord[i]);

    // find random number
    get_comp_rand_mp(rand);

    // add to sum1 & sum2
    sum_mul_mp(sum1, rand, temp1);
    sum_mul_mp(sum2, rand, temp2);
  }

  // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1
  mpf_abs_mp(temp1->r, sum1);
  mpf_abs_mp(temp1->i, sum2);
  if (mpf_cmp_ui(temp1->r, 0) > 0 && mpf_cmp_ui(temp1->i, 0) > 0)
  { // sum2 / sum1
    div_mp(temp1, sum2, sum1);

    // find the norm of sum2 / sum1
    c = d_abs_mp(temp1);

    // find the approximate cycle number
    c = log(sample_factor) / log(c);

    // find the cycle number - a positive integer
    if (c < 1)
      retVal = 1;
    else
    { // round c
      retVal = (int) floor(c + 0.5);
    }
  }
  else
  { // cycle number is 1
    retVal = 1;
  }

  // double check retVal >= 1
  if (retVal < 1)
    retVal = 1;

  clear_mp(temp1); clear_mp(temp2); clear_mp(sum1); clear_mp(sum2); clear_mp(rand);

  return retVal;
}

int find_cycle_num_loop_mp(PSEG_samples_struct_mp *S, int start, int finish, int max_cycle_num, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - error, otherwise best_cycle_num           *
* NOTES: finds the number between 1 and max_cycle_num that best *
* uses {start,..,finish-1} to approximate finish                *
\***************************************************************/
{
  int c, j, size, retVal, best_c = 1, n = finish - start - 1;
  mpf_t min_error, error;
  comp_mp tempComp;
  point_mp approx;

  // error checking
  if (start < 0 || finish > S->num_samples || max_cycle_num < 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize
  mpf_init(min_error); mpf_init(error);
  init_mp(tempComp);
  init_point_mp(approx, 0);

  mpf_set_d(min_error, 1e300);

  // setup s & dZ_ds
  _comp_mp *s = (_comp_mp *)bmalloc((n + 1) * sizeof(_comp_mp));
  _point_mp *dZ_ds = (_point_mp *)bmalloc((n + 1) * sizeof(_point_mp));
  for (j = 0; j <= n; j++)
  {
    init_mp(&s[j]);
    init_point_mp(&dZ_ds[j], 0);
  }

  for (c = 1; c <= max_cycle_num; c++)
  { // setup the structures to do an interpolation
    retVal = setup_PSEG_prediction_mp((c == 1), c, S, s, dZ_ds, start, finish, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // clear structures
      for (j = n; j >= 0; j--)
      {
        clear_mp(&s[j]);
        clear_point_mp(&dZ_ds[j]);
      }
      mpf_clear(min_error); mpf_clear(error);
      clear_mp(tempComp);
      clear_point_mp(approx);
      free(s);
      free(dZ_ds);
      // return error
      return -1;
    }

    // interpolate from 1, .. n
    if (n == 2)
      cubicInterp_mp(approx, &s[1], &S->Z_rev[1], &dZ_ds[1], &s[0]);
    else
      hermiteInterpCW_mp(approx, &s[1], &S->Z_rev[1], &dZ_ds[1], &s[0], n);

    // find the error
    mpf_set_ui(error, 0);
    size = approx->size;
    for (j = 0; j < size; j++)
    {
      sub_mp(&approx->coord[j], &approx->coord[j], &S->samples[finish-1].point->coord[j]);
      mpf_mul(tempComp->r, approx->coord[j].r, approx->coord[j].r);
      mpf_mul(tempComp->i, approx->coord[j].i, approx->coord[j].i);
      mpf_add(tempComp->r, tempComp->r, tempComp->i);
      if (mpf_cmp(tempComp->r, error) > 0)
        mpf_set(error, tempComp->r);
    }
    // check to see if this is the minimum error
    if (mpf_cmp(min_error, error) > 0)
    {
      mpf_set(min_error, error);
      best_c = c;
    }
  }

  // clear structures
  for (j = n; j >= 0; j--)
  {
    clear_mp(&s[j]);
    clear_point_mp(&dZ_ds[j]);
  }
  mpf_clear(min_error); mpf_clear(error);
  clear_mp(tempComp);
  clear_point_mp(approx);
  free(s);
  free(dZ_ds);

  return best_c;
}

int setup_PSEG_prediction_mp(int setup_dX_Z, int cycle_num, PSEG_samples_struct_mp *S, _comp_mp *s, _point_mp *dZ_ds, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE: setup the structures needed to do a prediction         *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - everything is setup, otherwise, some error *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0, i, j, k, rows, cols, size, n = finish - start - 1;
  mpf_t tempMPF;
  comp_mp factor, tempComp, c_times_t;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize
  mpf_init(tempMPF);
  init_mp(factor); init_mp(tempComp); init_mp(c_times_t);
  
  if (setup_dX_Z)
  { // setup dX & Z_rev
    eval_struct_mp e;
    init_eval_struct_mp(e, 0, 0, 0);

    for (i = start; i < finish; i++)
    { // evaluate the function
      eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples[i].point, S->samples[i].time, ED, eval_func);

      // find Y = -dH/dt = -dH/dp * dp/dt
      increase_size_vec_mp(e.funcVals, e.Jp->rows);
      rows = e.funcVals->size = e.Jp->rows;
      cols = e.Jp->cols;
      for (k = 0; k < rows; k++)
      {
        set_zero_mp(&e.funcVals->coord[k]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_mp(&e.funcVals->coord[k], &e.Jp->entry[k][j], &e.parDer->coord[j]);
        }
        neg_mp(&e.funcVals->coord[k], &e.funcVals->coord[k]);
      }

      // find dX = (dH/dX)^(-1) * (-dH/dt)
      retVal = matrixSolve_mp(S->dX[i], e.Jv, e.funcVals);
      if (retVal)
      {
        fprintf(OUT, "NOTE: matrixSolve has failed.\n");

        // clear memory
        mpf_clear(tempMPF);
        clear_mp(factor); clear_mp(tempComp); clear_mp(c_times_t);
        clear_eval_struct_mp(e);

        // return error
        return retVal;
      }

      // setup Z_rev
      k = finish - i - 1;
      increase_size_vec_mp(&S->Z_rev[k], S->samples[i].point->size);
      size = S->Z_rev[k].size = S->samples[i].point->size;
      for (j = 0; j < size; j++)
      {
        set_mp(&S->Z_rev[k].coord[j], &S->samples[i].point->coord[j]);
      }
    }
  }

  // setup dZ_ds & s
  mpf_set_ui(tempMPF, 1);
  mpf_div_ui(tempMPF, tempMPF, cycle_num);
  // since normalize time, set s[n] == 1
  pow_rmpf_mp(factor, S->samples[start].time, tempMPF);
  recip_mp(factor, factor);
  mpf_set_ui(s[n].r, 1); mpf_set_ui(s[n].i, 0);
  // find dZ_ds[n]
  mpf_mul_ui(c_times_t->r, S->samples[start].time->r, cycle_num);
  mpf_mul_ui(c_times_t->i, S->samples[start].time->i, cycle_num);
  vec_mulcomp_mp(&dZ_ds[n], S->dX[start], c_times_t);

  // loop over i = 0, .., n - 1 to find the dZ_ds[i] & s[i]
  for (i = 0; i < n; i++)
  { // find j
    j = finish - i - 1;
    // s_mp[i] = t^(1/c) * factor
    pow_rmpf_mp(&s[i], S->samples[j].time, tempMPF);
    mul_mp(&s[i], &s[i], factor);

    // dZ_ds_mp[i]
    k = cycle_num - 1;
    pow_rdouble_mp(tempComp, &s[i], k);
    mul_mp(tempComp, tempComp, c_times_t);
    vec_mulcomp_mp(&dZ_ds[i], S->dX[j], tempComp);
  }

  mpf_clear(tempMPF);
  clear_mp(factor); clear_mp(tempComp); clear_mp(c_times_t);

  return retVal;
}

int prediction_PSEG_samples_struct_mp(vec_mp approx, comp_mp finalT, int cycle_num, int setup_dX_Z, PSEG_samples_struct_mp *S, int start, int finish, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE: does a power series prediction based on samples to     *
* predict the value at finalT                                   *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - have an approx, otherwise, some error      *
* NOTES:                                                        *
\***************************************************************/
{
  int i, retVal = 0, n = finish - start - 1;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup s & dZ_ds
  _comp_mp *s = (_comp_mp *)bmalloc((n + 1) * sizeof(_comp_mp));
  _point_mp *dZ_ds = (_point_mp *)bmalloc((n + 1) * sizeof(_point_mp));

  // initialize
  for (i = 0; i <= n; i++)
  {
    init_mp(&s[i]);
    init_point_mp(&dZ_ds[i], 0);
  }

  // setup the structures
  retVal = setup_PSEG_prediction_mp(setup_dX_Z, cycle_num, S, s, dZ_ds, start, finish, OUT, ED, eval_func);

  // check for error
  if (!retVal)
  { // interpolate
    if (n + 1 == 2)
      cubicInterp_mp(approx, &s[0], &S->Z_rev[0], &dZ_ds[0], finalT);
    else
      hermiteInterpCW_mp(approx, &s[0], &S->Z_rev[0], &dZ_ds[0], finalT, n + 1);
  }

  // free the memory
  for (i = n; i >= 0; i--)
  {
    clear_mp(&s[i]);
    clear_point_mp(&dZ_ds[i]);
  }
  free(s);
  free(dZ_ds);

  return retVal;
}

int PSEG_struct_mp(point_data_mp *Final, point_mp last_approx, PSEG_samples_struct_mp *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual PSEG endgame in double precision       *
* storing the samples along the way to samples                  *
\***************************************************************/
{
  int i, dehom_prec, cycle_num = 0, retVal = 0, start, finish, numSamples = T->num_PSEG_sample_points;
  double approx_err, norm_last = 0, norm_new = 0;
  comp_mp endTime, finalTime;   
  point_mp new_approx, dehom;

  // initialize MP
  init_mp(endTime); init_mp(finalTime);
  init_point_mp(new_approx, 0); init_point_mp(dehom, 0);

  // initialize cycle number
  Final->cycle_num = 0;

  // make sure that there is atleast 1 good sample point
  if (PSEG_samples->num_samples <= 0)
  {
    printf("WARNING: The number of good sample points must be positive in PSEG_struct_d!!\n");

    // clear MP
    clear_mp(endTime); clear_mp(finalTime);
    clear_point_mp(new_approx); clear_point_mp(dehom);

    return retVal_PSEG_failed;
  }

  // setup start and finish
  finish = PSEG_samples->num_samples - 1;
  start = finish - numSamples;
  if (start < 0)
    start = 0;

  // setup finalTime
  set_double_mp(finalTime, T->targetT, 0);

  // make sure there are enough locations to find the initial set of sample points
  if (PSEG_samples->mem_count <= numSamples)
  { // increase the memory
    change_size_PSEG_samples_struct_mp(PSEG_samples, numSamples + 1);
  }

  // find the initial set of samples points to create the first approximation
  for (i = finish + 1; i <= numSamples; i++)
  { // find the end time for the next samples
    mul_rdouble_mp(endTime, PSEG_samples->samples[i-1].time, T->power_series_sample_factor);

    // refine the sample point
    refine_PSEG_samples_struct_mp(PSEG_samples, i-1, i, T, OUT, ED, eval_func);
    // track to the next sample point
    retVal = track_mp(&PSEG_samples->samples[i], &PSEG_samples->samples[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      point_data_cp_mp(Final, &PSEG_samples->samples[i]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[i].time->r);

      // clear MP
      clear_mp(endTime); clear_mp(finalTime);
      clear_point_mp(new_approx); clear_point_mp(dehom);

      // return error
      return retVal;
    }
    else
    { // update the number of good samples
      PSEG_samples->num_samples++;
    }
  }
  // so, we have enough sample points for the initial approximation

  // find the start & end
  finish = PSEG_samples->num_samples;
  start = finish - numSamples - 1;

  // find the first approximation
  retVal = PSEG_samples_prediction_mp(last_approx, finalTime, &cycle_num, PSEG_samples, start, finish, T, OUT, ED, eval_func);

  // check for error
  if (retVal)
  { // total failure since we cannot find the first approximation - setup Final and return error
    point_data_cp_mp(Final, &PSEG_samples->samples[finish-1]);
    Final->cycle_num = 0;
    T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[finish-1].time->r);

    // clear MP
    clear_mp(endTime); clear_mp(finalTime);
    clear_point_mp(new_approx); clear_point_mp(dehom);

    // return error
    return retVal;
  }
  else if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(NULL, dehom, &dehom_prec, NULL, last_approx, T->Precision, NULL, ED);
    norm_last = infNormVec_mp(dehom);
  }

  // now that we have the first approximation, we can loop until we have convergence
  while (1)
  { // find the index of the next sample and the current precision
    i = PSEG_samples->num_samples;
    // setup endTime
    mul_rdouble_mp(endTime, PSEG_samples->samples[i-1].time, T->power_series_sample_factor);

    // make sure there is enough memory
    if (PSEG_samples->mem_count <= PSEG_samples->num_samples)
    { // increase the memory
      change_size_PSEG_samples_struct_mp(PSEG_samples, 2 * PSEG_samples->mem_count);
    }

    // refine the sample point
    refine_PSEG_samples_struct_mp(PSEG_samples, i-1, i, T, OUT, ED, eval_func);
    // track to the next sample point
    retVal = track_mp(&PSEG_samples->samples[i], &PSEG_samples->samples[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      point_data_cp_mp(Final, &PSEG_samples->samples[i]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[i].time->r);

      // break out of the loop to return error
      break;
    }

    // update the number of good samples
    PSEG_samples->num_samples++;

    // find the start & end
    finish = PSEG_samples->num_samples;
    start = finish - numSamples - 1;

    // find the next approximation
    retVal = PSEG_samples_prediction_mp(new_approx, finalTime, &cycle_num, PSEG_samples, start, finish, T, OUT, ED, eval_func);

    // check for error
    if (retVal)
    { // total failure since we cannot find the first approximation - setup Final and return error
      point_data_cp_mp(Final, &PSEG_samples->samples[finish-1]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[finish-1].time->r);

      // break out of the loop to return error
      break;
    }
    else if (T->securityLevel <= 0)
    { // find dehom & its norm
      find_dehom(NULL, dehom, &dehom_prec, NULL, new_approx, T->Precision, NULL, ED);
      norm_new = infNormVec_mp(dehom);
    }

    // compare the approximations
    T->error_at_latest_sample_point = approx_err = compare_approximations_mp(new_approx, last_approx);

    // check for convergence
    if (approx_err < T->final_tolerance)
    { // accept the approximation

      // setup Final
      point_cp_mp(Final->point, new_approx);
      set_mp(Final->time, finalTime);
      Final->cycle_num = cycle_num;

      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[finish-1].time->r);

      // break out of loop
      retVal = 0;
      break;
    }
    else if (mpf_get_d(PSEG_samples->samples[finish-1].time->r) < T->minTrackT)
    { // too close to t = 0 but we do not have the correct tolerance - so we exit

      // setup Final
      point_cp_mp(Final->point, new_approx);
      set_mp(Final->time, finalTime);
      Final->cycle_num = cycle_num;

      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[finish-1].time->r);

      // break out of loop
      retVal = retVal_EG_failed_to_converge;
      break;
    }
    else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
    { // we are too large
      point_data_cp_mp(Final, &PSEG_samples->samples[finish-1]);
      Final->cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples[finish-1].time->r);

      // break out of the loop to return error
      retVal = retVal_security_max;
      break;
    }
    else
    { // simply copy new_approx to last_approx and try again
      point_cp_mp(last_approx, new_approx);
      norm_last = norm_new;
    }
  }

  // clear MP
  clear_mp(endTime); clear_mp(finalTime);
  clear_point_mp(new_approx); clear_point_mp(dehom);

  return retVal;
}

