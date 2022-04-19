// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

// This file contains the power series endgame using adaptive multiprecision tracking

int PSEG_struct_amp(int *final_prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, PSEG_samples_struct_amp *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int prediction_PSEG_samples_struct_amp(vec_d approx_d, vec_mp approx_mp, int *approx_prec, comp_d finalT_d, comp_mp finalT_mp, int currPrec, int cycle_num, int setup_dX_Z, PSEG_samples_struct_amp *S, int start, int finish, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

double compare_approximations_amp(point_d approx0_d, point_mp approx0_mp, int approx0_prec, point_d approx1_d, point_mp approx1_mp, int approx1_prec);
int find_cycle_num_amp(PSEG_samples_struct_amp *S, int start, int currPrec, double sample_factor);
int find_cycle_num_loop_amp(PSEG_samples_struct_amp *S, int start, int finish, int currPrec, int max_cycle_num, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int setup_PSEG_prediction_amp(int setup_dX_Z, int currPrec, int cycle_num, PSEG_samples_struct_amp *S, _comp_d *s_d, _point_d *dZ_ds_d, _comp_mp *s_mp, _point_mp *dZ_ds_mp, int start, int finish, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int PSEG_samples_prediction_amp(point_d approx_d, point_mp approx_mp, int *approx_prec, comp_d finalTime_d, comp_mp finalTime_mp, int currPrec, int *cycle_num, PSEG_samples_struct_amp *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void refine_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S, int start, int finish, int prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void init_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S, int num_points, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initiailizes S                                         *
\***************************************************************/
{
  int i;

  // verify prec >= 64 since this is the precision that MP values will be initialized to
  prec = MAX(prec, 64);

  S->max_prec = (int *)bmalloc(num_points * sizeof(int));
  S->samples_d  = (point_data_d *)bmalloc(num_points * sizeof(point_data_d));
  S->samples_mp = (point_data_mp *)bmalloc(num_points * sizeof(point_data_mp));
  S->dX_d  = (vec_d *)bmalloc(num_points * sizeof(vec_d));
  S->dX_mp = (vec_mp *)bmalloc(num_points * sizeof(vec_mp));
  S->Z_rev_d  = (_point_d *)bmalloc(num_points * sizeof(_point_d));
  S->Z_rev_mp = (_point_mp *)bmalloc(num_points * sizeof(_point_mp));
  for (i = 0; i < num_points; i++)
  { // initialize
    init_point_data_d(&S->samples_d[i], 0);
    init_point_data_mp2(&S->samples_mp[i], 0, prec);
    init_vec_d(S->dX_d[i], 0);
    init_vec_mp2(S->dX_mp[i], 0, prec);
    init_point_d(&S->Z_rev_d[i], 0);
    init_point_mp2(&S->Z_rev_mp[i], 0, prec);

    // initialize max_prec
    S->max_prec[i] = prec;
  }
  S->num_samples = 0;
  S->mem_count = num_points;

  return;
}

void change_size_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S, int num_new_points, int currPrec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: changes the size of S                                  *
\***************************************************************/
{
  int i, old_sample_size = S->num_samples, newPrec = MAX(currPrec, 64);

  if (S->mem_count != num_new_points)
  { // store the old MP samples
    point_data_mp *old_samples_mp = S->samples_mp;
    vec_mp *old_dX_mp = S->dX_mp;
    _point_mp *old_Z_mp = S->Z_rev_mp;

    // setup max_prec
    S->max_prec = (int *)brealloc(S->max_prec, num_new_points * sizeof(int));

    // setup the double precision memory
    S->samples_d = (point_data_d *)brealloc(S->samples_d, num_new_points * sizeof(point_data_d));
    S->dX_d = (vec_d *)brealloc(S->dX_d, num_new_points * sizeof(vec_d));
    S->Z_rev_d = (_point_d *)brealloc(S->Z_rev_d, num_new_points * sizeof(_point_d));
    for (i = S->num_samples; i < num_new_points; i++)
    { // initialize the new points
      init_point_data_d(&S->samples_d[i], 0);
      init_vec_d(S->dX_d[i], 0);
      init_point_d(&S->Z_rev_d[i], 0);
      S->samples_d[i].point->size = S->dX_d[i]->size = S->Z_rev_d[i].size = 0;
      S->samples_d[i].cycle_num = 0;
      // set the precision for the new points
      S->max_prec[i] = newPrec;
    }

    // setup the multi precision memory
    S->samples_mp = NULL;
    S->dX_mp = NULL;
    S->Z_rev_mp = NULL;
    S->samples_mp = (point_data_mp *)bmalloc(num_new_points * sizeof(point_data_mp));
    S->dX_mp = (vec_mp *)bmalloc(num_new_points * sizeof(vec_mp));
    S->Z_rev_mp = (_point_mp *)bmalloc(num_new_points * sizeof(_point_mp));
    // initialize new memory and copy over old memory
    for (i = num_new_points - 1; i >= 0; i--)
    {
      init_point_data_mp2(&S->samples_mp[i], 0, S->max_prec[i]);
      init_vec_mp2(S->dX_mp[i], 0, S->max_prec[i]);
      init_point_mp2(&S->Z_rev_mp[i], 0, S->max_prec[i]);

      if (i < old_sample_size)
      { // copy over the old data - Z does not need copied over
        point_data_cp_mp(&S->samples_mp[i], &old_samples_mp[i]);
        point_cp_mp(S->dX_mp[i], old_dX_mp[i]);
      }
    }
    S->mem_count = num_new_points;

    // clear all of the old MP memory
    for (i = old_sample_size - 1; i >= 0; i--)
    {
      clear_point_data_mp(&old_samples_mp[i]);
      clear_vec_mp(old_dX_mp[i]);
      clear_point_mp(&old_Z_mp[i]);
    }
    free(old_samples_mp);
    free(old_dX_mp);
    free(old_Z_mp);
  }

  return;
}

void clear_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears S                                               *
\***************************************************************/
{
  int i;

  for (i = S->mem_count - 1; i >= 0; i--)
  { // clear 
    clear_point_data_d(&S->samples_d[i]);
    clear_point_data_mp(&S->samples_mp[i]);
    clear_vec_d(S->dX_d[i]);
    clear_vec_mp(S->dX_mp[i]);
    clear_point_d(&S->Z_rev_d[i]);
    clear_point_mp(&S->Z_rev_mp[i]);
  }
  // free structures
  free(S->samples_d);
  free(S->samples_mp);
  free(S->dX_d);
  free(S->dX_mp);
  free(S->Z_rev_d);
  free(S->Z_rev_mp);
  free(S->max_prec);
  S->num_samples = 0;
  S->mem_count = 0;

  return;
}

int PSEG_amp(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: prec - the ending precision - < 64 - Final_d   *
*                                             >= 64 - Final_mp  *
* time_first_increase - the time of the last increase to MP     *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, cycle_num;
  point_data_d lastSample_d;
  point_data_mp lastSample_mp;

  init_point_data_d(&lastSample_d, 0);
  init_point_data_mp(&lastSample_mp, 0);

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  retVal = PSEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &lastSample_d, &lastSample_mp, &cycle_num, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  // clear 
  clear_point_data_d(&lastSample_d);
  clear_point_data_mp(&lastSample_mp);

  return retVal;
}

int PSEG_rank_amp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Power series endgame using adaptive precision and *
* uses Cauchy endgame for rank determination                    *
\***************************************************************/
{
  int i, retVal, cycle_num, samples_per_loop = T->num_PSEG_sample_points + 1;
  comp_d finalTime_d;
  comp_mp finalTime_mp;
  point_data_d lastSample_d;
  point_data_mp lastSample_mp;

  init_mp(finalTime_mp);
  init_point_data_d(&lastSample_d, 0);
  init_point_data_mp(&lastSample_mp, 0);

  // setup finalTime
  set_double_d(finalTime_d, T->targetT, 0);
  d_to_mp(finalTime_mp, finalTime_d);

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  // perform the Power Series endgame - return that last sample so that it can be used for rank determination
  retVal = PSEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &lastSample_d, &lastSample_mp, &cycle_num, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_amp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, *last_approx_prec, finalTime_d, finalTime_mp, &lastSample_d, &lastSample_mp, *prec, cycle_num, samples_per_loop, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;

    // setup last_approx 
    *last_approx_prec = *prec;
    if (*prec < 64)
    { // setup _d
      point_cp_d(last_approx_d, Final_d->point);
      for (i = 0; i < last_approx_d->size; i++)
      {
        get_comp_rand_d(finalTime_d);
        mul_rdouble_d(finalTime_d, finalTime_d, T->final_tolerance);
        add_d(&last_approx_d->coord[i], &last_approx_d->coord[i], finalTime_d);
      }
    }
    else
    { // setup _mp
      setprec_point_mp(last_approx_mp, *last_approx_prec);
      point_cp_mp(last_approx_mp, Final_mp->point);
      for (i = 0; i < last_approx_mp->size; i++)
      {
        get_comp_rand_mp(finalTime_mp);
        mul_rdouble_mp(finalTime_mp, finalTime_mp, T->final_tolerance);
        add_mp(&last_approx_mp->coord[i], &last_approx_mp->coord[i], finalTime_mp);
      }
    }
  }

  // clear
  clear_mp(finalTime_mp);
  clear_point_data_d(&lastSample_d);
  clear_point_data_mp(&lastSample_mp);

  return retVal;
}

int PSEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *lastSample_d, point_data_mp *lastSample_mp, int *cycle_num, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: prec - the ending precision - < 64 - Final_d   *
*                                             >= 64 - Final_mp  *
* time_first_increase - the time of the last increase to MP     *
* NOTES: returns the last sample point                          *
\***************************************************************/
{
  int i, retVal, numSamples = T->num_PSEG_sample_points;
  comp_d endTime_d;
  comp_mp endTime_mp;

  PSEG_samples_struct_amp PSEG_samples;

  // start off with atleast 64 bit precision
  T->Precision = MAX(64, prec_in);
  initMP(T->Precision);
  change_prec(ED_mp, T->Precision);

  // initialize PSEG_samples
  init_PSEG_samples_struct_amp(&PSEG_samples, numSamples + 1, T->Precision);

  // initialize MP data types
  init_mp(endTime_mp);

  // initialze prec & time_first_increase
  *prec = prec_in; // start off at the intial precision
  *time_first_increase = 0; // no increase so far!
  Final_d->cycle_num = Final_mp->cycle_num = 0;

  // setup to track to the endgame boundary
  set_double_d(endTime_d, T->endgameBoundary, 0.0);
  set_double_mp(endTime_mp, T->endgameBoundary, 0.0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0.0;

  // setup last_approx in case of failure
  *last_approx_prec = prec_in;
  if (prec_in < 64)
  {
    point_cp_d(last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(last_approx_mp, prec_in);
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  // try to track to the endgame boundary - this will handle everything
  retVal = AMPtrack(&PSEG_samples.samples_d[0], &PSEG_samples.samples_mp[0], prec, time_first_increase, Start_d, Start_mp, prec_in, endTime_d, endTime_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // store the precision of the first sample point
  PSEG_samples.max_prec[0] = *prec;

  // so, either made it to the endgame boundary or path failure - first print to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  if (*prec == 52)
  { // print to midOUT using samples_d
    for (i = 0; i < PSEG_samples.samples_d[0].point->size; i++)
      fprintf(midOUT, "%.15e %.15e\n", PSEG_samples.samples_d[0].point->coord[i].r, PSEG_samples.samples_d[0].point->coord[i].i);
    fprintf(midOUT, "\n");
  }
  else
  { // print to midOUT using samples_mp
    for (i = 0; i < PSEG_samples.samples_mp[0].point->size; i++)
    {
      print_mp(midOUT, 0, &PSEG_samples.samples_mp[0].point->coord[i]);
      fprintf(midOUT, "\n");
    }
    fprintf(midOUT, "\n");
  }
  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check to see if we made it to the endgame boundary successfully
  if (retVal)
  { // had path failure
    fprintf(OUT, "NOTE: Power series endgame never started!\n");

    // copy the last sample to Final
    if (*prec == 52)
    {
      point_data_cp_d(Final_d, &PSEG_samples.samples_d[0]);
      Final_d->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples.samples_d[0].time->r;
    }
    else
    {
      setprec_point_mp(Final_mp->point, *prec);
      setprec_mp(Final_mp->time, *prec);
      point_data_cp_mp(Final_mp, &PSEG_samples.samples_mp[0]);
      Final_mp->cycle_num = *cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples.samples_mp[0].time->r);
    }
  }
  else
  { // success - setup to run the acutal endgame!

    // we have 1 sample point
    PSEG_samples.num_samples = 1;

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;    

    // run the actual PSEG
    retVal = PSEG_struct_amp(prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &PSEG_samples, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    // setup lastSample & cycle_num
    if (*prec < 64)
    { // use _d
      *cycle_num = Final_d->cycle_num;

      // setup lastSample_d
      if (PSEG_samples.max_prec[PSEG_samples.num_samples - 1] < 64)
      { // use samples_d
        point_data_cp_d(lastSample_d, &PSEG_samples.samples_d[PSEG_samples.num_samples - 1]);
      }
      else
      { // use samples_mp
        convert_point_data_mp_to_d(lastSample_d, &PSEG_samples.samples_mp[PSEG_samples.num_samples - 1]); 
      }
    }
    else
    { // use _mp
      *cycle_num = Final_mp->cycle_num;

      // setup lastSample_mp
      setprec_point_mp(lastSample_mp->point, *prec);
      setprec_mp(lastSample_mp->time, *prec);

      if (PSEG_samples.max_prec[PSEG_samples.num_samples - 1] < 64)
      { // use samples_d
        convert_point_data_d_to_mp(lastSample_mp, &PSEG_samples.samples_d[PSEG_samples.num_samples - 1]);
      }
      else
      { // use samples_mp
        point_data_cp_mp(lastSample_mp, &PSEG_samples.samples_mp[PSEG_samples.num_samples - 1]);
      }
    }
  }

  // clear the MP data types
  clear_mp(endTime_mp);
  clear_PSEG_samples_struct_amp(&PSEG_samples);

  return retVal;
}

int PSEG_struct_amp(int *final_prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, PSEG_samples_struct_amp *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: prec - the ending precision - < 64 - Final_d   *
*                                             >= 64 - Final_mp  *
* time_first_increase - the time of the last increase to MP     *
* NOTES: does the actual AMP PSEG                               *
\***************************************************************/
{
  int i, dehom_prec = 52, currPrec, new_approx_prec, cycle_num = 0, retVal = 0, start, finish, numSamples = T->num_PSEG_sample_points;
  double approx_err, norm_last = 0, norm_new = 0;
  comp_d endTime_d, finalTime_d;
  comp_mp endTime_mp, finalTime_mp;
  point_d new_approx_d, dehom_d;
  point_mp new_approx_mp, dehom_mp;

  // initialize cycle number
  Final_d->cycle_num = Final_mp->cycle_num = 0;

  // make sure there is atleast 1 sample point
  if (PSEG_samples->num_samples <= 0)
  {
    printf("WARNING: The number of good sample points must be positive in PSEG_struct_amp!!\n");
    return retVal_PSEG_failed;
  }

  // initialize MP
  init_mp(endTime_mp); init_mp(finalTime_mp);
  init_point_d(new_approx_d, 0); init_point_d(dehom_d, 0);
  init_point_mp(new_approx_mp, 0); init_point_mp(dehom_mp, 0);

  // setup start and finish
  finish = PSEG_samples->num_samples - 1;
  start = finish - numSamples;
  if (start < 0)
    start = 0;

  // find the current precision of the last sample point
  currPrec = PSEG_samples->max_prec[finish];

  // setup finalTime
  set_double_d(finalTime_d, T->targetT, 0);
  d_to_mp(finalTime_mp, finalTime_d);

  // make sure there are enough locations to find the initial set of sample points
  if (PSEG_samples->mem_count <= numSamples)
  { // increase the memory
    change_size_PSEG_samples_struct_amp(PSEG_samples, numSamples + 1, currPrec);
  }

  // find the initial set of sample points to create the first approximation
  for (i = finish + 1; i <= numSamples; i++)
  { // find the end time for the next sample
    currPrec = PSEG_samples->max_prec[i-1];
    if (currPrec < 64)
    { // setup endTime_d
      mul_rdouble_d(endTime_d, PSEG_samples->samples_d[i-1].time, T->power_series_sample_factor);
    }
    else
    { // setup endTime_mp
      setprec_mp(endTime_mp, currPrec);
      mul_rdouble_mp(endTime_mp, PSEG_samples->samples_mp[i-1].time, T->power_series_sample_factor);
    }

    // track to the next sample point
    retVal = AMPtrack(&PSEG_samples->samples_d[i], &PSEG_samples->samples_mp[i], &PSEG_samples->max_prec[i], time_first_increase, &PSEG_samples->samples_d[i-1], &PSEG_samples->samples_mp[i-1], currPrec, endTime_d, endTime_mp, currPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      *final_prec = PSEG_samples->max_prec[i];
      if (*final_prec < 64)
      { // setup Final_d
        point_data_cp_d(Final_d, &PSEG_samples->samples_d[i]);
        Final_d->cycle_num = 0;
        T->t_val_at_latest_sample_point = PSEG_samples->samples_d[i].time->r;
      }
      else
      { // setup Final_mp
        point_data_cp_mp(Final_mp, &PSEG_samples->samples_mp[i]);
        Final_mp->cycle_num = 0;
        T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[i].time->r);
      }

      // clear MP
      clear_mp(endTime_mp); clear_mp(finalTime_mp);
      clear_point_d(new_approx_d); clear_point_d(dehom_d);
      clear_point_mp(new_approx_mp); clear_point_mp(dehom_mp);

      // return error
      return retVal;
    }
    else
    { // update the number of good samples
      PSEG_samples->num_samples++;
    }
  }
  // so, we have enough sample points for the initial approximation

  // find the start & end and the current precision
  finish = PSEG_samples->num_samples;
  start = finish - numSamples - 1;
  currPrec = PSEG_samples->max_prec[finish - 1];

  // find the first approximation
  retVal = PSEG_samples_prediction_amp(last_approx_d, last_approx_mp, last_approx_prec, finalTime_d, finalTime_mp, currPrec, &cycle_num, PSEG_samples, start, finish, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // check for errors
  if (retVal)
  { // total failure since we cannot find the first approximation - setup Final, clear MP and return error
    *final_prec = PSEG_samples->max_prec[finish - 1];
    if (*final_prec < 64)
    { // setup Final_d
      point_data_cp_d(Final_d, &PSEG_samples->samples_d[finish - 1]);
      Final_d->cycle_num = 0;
      T->t_val_at_latest_sample_point = PSEG_samples->samples_d[finish - 1].time->r;
    }
    else
    { // setup Final_mp
      point_data_cp_mp(Final_mp, &PSEG_samples->samples_mp[finish - 1]);
      Final_mp->cycle_num = 0;
      T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[finish - 1].time->r);
    }

    // clear MP
    clear_mp(endTime_mp); clear_mp(finalTime_mp);
    clear_point_d(new_approx_d); clear_point_d(dehom_d);
    clear_point_mp(new_approx_mp); clear_point_mp(dehom_mp);

    // return error
    return retVal;
  }
  else if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom_d, dehom_mp, &dehom_prec, last_approx_d, last_approx_mp, *last_approx_prec, ED_d, ED_mp);
    if (dehom_prec < 64)
      norm_last = infNormVec_d(dehom_d);
    else
      norm_last = infNormVec_mp(dehom_mp);
  }

  // now that we have the first approximation, we can loop until we have convergence
  while (1)
  { // find the index of the next sample and the current precision
    i = PSEG_samples->num_samples;
    currPrec = PSEG_samples->max_prec[i-1];
    // setup endTime
    if (currPrec < 64)
    { // setup endTime_d
      mul_rdouble_d(endTime_d, PSEG_samples->samples_d[i-1].time, T->power_series_sample_factor);
    }
    else
    { // setup endTime_mp
      setprec_mp(endTime_mp, currPrec);
      mul_rdouble_mp(endTime_mp, PSEG_samples->samples_mp[i-1].time, T->power_series_sample_factor);
    }

    // make sure there is enough memory
    if (PSEG_samples->mem_count <= PSEG_samples->num_samples)
    { // increase the memory
      change_size_PSEG_samples_struct_amp(PSEG_samples, 2 * PSEG_samples->mem_count, currPrec);
    }

    // track to the next sample point
    retVal = AMPtrack(&PSEG_samples->samples_d[i], &PSEG_samples->samples_mp[i], &PSEG_samples->max_prec[i], time_first_increase, &PSEG_samples->samples_d[i-1], &PSEG_samples->samples_mp[i-1], currPrec, endTime_d, endTime_mp, currPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for errors
    if (retVal)
    { // tracking did not make it to the next sample point - copy last sample to Final and return path failure
      *final_prec = PSEG_samples->max_prec[i];
      if (*final_prec < 64)
      { // setup Final_d
        point_data_cp_d(Final_d, &PSEG_samples->samples_d[i]);
        Final_d->cycle_num = 0;
        T->t_val_at_latest_sample_point = PSEG_samples->samples_d[i].time->r;
      }
      else
      { // setup Final_mp
        point_data_cp_mp(Final_mp, &PSEG_samples->samples_mp[i]);
        Final_mp->cycle_num = 0;
        T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[i].time->r);
      }

      // break out of the loop to return error
      break;
    }

    // update the number of good samples
    PSEG_samples->num_samples++;

    // find the start & end and the current precision
    finish = PSEG_samples->num_samples;
    start = finish - numSamples - 1;
    currPrec = PSEG_samples->max_prec[finish - 1];

    // find the next approximation
    retVal = PSEG_samples_prediction_amp(new_approx_d, new_approx_mp, &new_approx_prec, finalTime_d, finalTime_mp, currPrec, &cycle_num, PSEG_samples, start, finish, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for error
    if (retVal)
    { // total failure since we cannot find the next approximation - setup Final and break out of the loop to return error
      *final_prec = PSEG_samples->max_prec[finish - 1];
      if (*final_prec < 64)
      { // setup Final_d
        point_data_cp_d(Final_d, &PSEG_samples->samples_d[finish - 1]);
        Final_d->cycle_num = 0;
        T->t_val_at_latest_sample_point = PSEG_samples->samples_d[finish - 1].time->r;
      }
      else
      { // setup Final_mp
        point_data_cp_mp(Final_mp, &PSEG_samples->samples_mp[finish - 1]);
        Final_mp->cycle_num = 0;
        T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[finish - 1].time->r);
      }

      break;
    }
    else if (T->securityLevel <= 0)
    { // find dehom & its norm
      find_dehom(dehom_d, dehom_mp, &dehom_prec, new_approx_d, new_approx_mp, new_approx_prec, ED_d, ED_mp);
      if (dehom_prec < 64)
        norm_new = infNormVec_d(dehom_d);
      else
        norm_new = infNormVec_mp(dehom_mp);
    }

    // compare the approximations
    T->error_at_latest_sample_point = approx_err = compare_approximations_amp(new_approx_d, new_approx_mp, new_approx_prec, last_approx_d, last_approx_mp, *last_approx_prec);

    // check for convergence
    if (approx_err < T->final_tolerance)
    { // accept the approximation

      // setup Final
      *final_prec = new_approx_prec;
      if (new_approx_prec < 64)
      { // setup Final_d
        point_cp_d(Final_d->point, new_approx_d);
        set_d(Final_d->time, finalTime_d);
        Final_d->cycle_num = cycle_num;
      }
      else
      { // setup Final_mp
        setprec_point_mp(Final_mp->point, *final_prec);
        setprec_mp(Final_mp->time, *final_prec);

        point_cp_mp(Final_mp->point, new_approx_mp);
        set_mp(Final_mp->time, finalTime_mp);
        Final_mp->cycle_num = cycle_num;
      }
      if (PSEG_samples->max_prec[finish-1] < 64)
        T->t_val_at_latest_sample_point = PSEG_samples->samples_d[finish-1].time->r;
      else
        T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[finish-1].time->r);
      // break out of loop
      retVal = 0;
      break;
    }
    else
    { // check to see if we are too close to t = 0
      retVal = 0;
      if (PSEG_samples->max_prec[finish-1] < 64)
      { // check double precision
        if (PSEG_samples->samples_d[finish-1].time->r < T->minTrackT)
          retVal = 1;
      }
      else
      { // check multi precision
        if (mpf_get_d(PSEG_samples->samples_mp[finish-1].time->r) < T->minTrackT)
          retVal = 1;
      }

      if (retVal)
      { // too close to t = 0
        // setup Final
        *final_prec = new_approx_prec;
        if (new_approx_prec < 64)
        { // setup Final_d
          point_cp_d(Final_d->point, new_approx_d);
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = cycle_num;
        }
        else
        { // setup Final_mp
          setprec_point_mp(Final_mp->point, *final_prec);
          setprec_mp(Final_mp->time, *final_prec);

          point_cp_mp(Final_mp->point, new_approx_mp);
          set_mp(Final_mp->time, finalTime_mp);
          Final_mp->cycle_num = cycle_num;
        }
        if (PSEG_samples->max_prec[finish-1] < 64)
          T->t_val_at_latest_sample_point = PSEG_samples->samples_d[finish-1].time->r;
        else
          T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[finish-1].time->r);
        // break out of loop
        retVal = retVal_EG_failed_to_converge;
        break;
      }
      else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
      { // we are too large
        *final_prec = PSEG_samples->max_prec[i];
        if (*final_prec < 64)
        { // setup Final_d
          point_data_cp_d(Final_d, &PSEG_samples->samples_d[i]);
          Final_d->cycle_num = 0;
          T->t_val_at_latest_sample_point = PSEG_samples->samples_d[i].time->r;
        }
        else
        { // setup Final_mp
          point_data_cp_mp(Final_mp, &PSEG_samples->samples_mp[i]);
          Final_mp->cycle_num = 0;
          T->t_val_at_latest_sample_point = mpf_get_d(PSEG_samples->samples_mp[i].time->r);
        }

        // break out of the loop to return error
        retVal = retVal_security_max;
        break;
      }
      else
      { // simply copy new_approx to last_approx and try again
        *last_approx_prec = new_approx_prec;
        if (new_approx_prec < 64)
        { // copy to last_approx_d
          point_cp_d(last_approx_d, new_approx_d);
        }
        else
        { // copy to last_approx_mp
          setprec_point_mp(last_approx_mp, *last_approx_prec);
          point_cp_mp(last_approx_mp, new_approx_mp);
        }
        norm_last = norm_new;
      }
    }
  }

  // clear MP
  clear_mp(endTime_mp); clear_mp(finalTime_mp);
  clear_point_d(new_approx_d); clear_point_d(dehom_d);
  clear_point_mp(new_approx_mp); clear_point_mp(dehom_mp);

  return retVal;
}

double compare_approximations_amp(point_d approx0_d, point_mp approx0_mp, int approx0_prec, point_d approx1_d, point_mp approx1_mp, int approx1_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compares the approximations and returns the error      *
\***************************************************************/
{
  int i, size;
  double error;

  if (approx0_prec <= approx1_prec)
  { // compare using the precision of approx1
    if (approx1_prec < 64)
    { // both in double precision
      double tempD;
      comp_d tempComp;

      // find the error
      size = approx1_d->size;
      error = 0;
      for (i = 0; i < size; i++)
      {
        sub_d(tempComp, &approx0_d->coord[i], &approx1_d->coord[i]);
        tempD = tempComp->r * tempComp->r + tempComp->i * tempComp->i;
        if (tempD > error)
          error = tempD;
      }
      error = sqrt(error);
    }
    else
    { // since approx1 is using MP, move approx0 to the same precision and then find error in MP
      mpf_t error_mp;
      comp_mp tempCoord;

      mpf_init2(error_mp, approx1_prec);
      init_mp2(tempCoord, approx1_prec);

      // find the error
      size = approx1_mp->size;
      mpf_set_ui(error_mp, 0);
      for (i = 0; i < size; i++)
      { // setup tempCoord
        if (approx0_prec < 64)
        { // move _d to _mp
          d_to_mp(tempCoord, &approx0_d->coord[i]);
        }
        else
        { // just use set
          set_mp(tempCoord, &approx0_mp->coord[i]);
        }

        // find the error in this coordinate
        sub_mp(tempCoord, tempCoord, &approx1_mp->coord[i]);
        mpf_mul(tempCoord->r, tempCoord->r, tempCoord->r);
        mpf_mul(tempCoord->i, tempCoord->i, tempCoord->i);  
        mpf_add(tempCoord->r, tempCoord->r, tempCoord->i);
        if (mpf_cmp(tempCoord->r, error_mp) > 0)
          mpf_set(error_mp, tempCoord->r);
      }
      mpf_sqrt(error_mp, error_mp);
      // set error
      error = mpf_get_d(error_mp);

      mpf_clear(error_mp);
      clear_mp(tempCoord);
    }
  }
  else // approx0_prec > approx1_prec
  { // so we know that approx0 is in MP
    mpf_t error_mp;
    comp_mp tempCoord;

    mpf_init2(error_mp, approx0_prec);
    init_mp2(tempCoord, approx0_prec);

    // find the error
    size = approx0_mp->size;
    mpf_set_ui(error_mp, 0);
    for (i = 0; i < size; i++)
    { // setup tempCoord
      if (approx1_prec < 64)
      { // move _d to _mp
        d_to_mp(tempCoord, &approx1_d->coord[i]);
      }
      else
      { // just use set
        set_mp(tempCoord, &approx1_mp->coord[i]);
      }

      // find the error in this coordinate
      sub_mp(tempCoord, tempCoord, &approx0_mp->coord[i]);
      mpf_mul(tempCoord->r, tempCoord->r, tempCoord->r);
      mpf_mul(tempCoord->i, tempCoord->i, tempCoord->i);
      mpf_add(tempCoord->r, tempCoord->r, tempCoord->i);
      if (mpf_cmp(tempCoord->r, error_mp) > 0)
        mpf_set(error_mp, tempCoord->r);
    }
    mpf_sqrt(error_mp, error_mp);
    // set error
    error = mpf_get_d(error_mp);

    mpf_clear(error_mp);
    clear_mp(tempCoord);
  }
  
  return error;
}

int PSEG_samples_prediction_amp(point_d approx_d, point_mp approx_mp, int *approx_prec, comp_d finalTime_d, comp_mp finalTime_mp, int currPrec, int *cycle_num, PSEG_samples_struct_amp *PSEG_samples, int start, int finish, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
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

  // refine the sample points to the current precision
  refine_PSEG_samples_struct_amp(PSEG_samples, start, finish, currPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

  // find an approximation based on the last 3 samples
  *cycle_num = find_cycle_num_amp(PSEG_samples, finish - 3, currPrec, T->power_series_sample_factor);

  // use this approximation as a way to determine where to stop the exhaustive search for finding the best cycle number
  *cycle_num = find_cycle_num_loop_amp(PSEG_samples, start, finish, currPrec, MAX(*cycle_num * 5, T->cycle_num_max), OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

  if (*cycle_num > 0)
  { // the structures are setup properly - find the prediction

    // store the cycle number to [start]
    if (currPrec < 64)
      PSEG_samples->samples_d[start].cycle_num = *cycle_num;
    else
      PSEG_samples->samples_mp[start].cycle_num = *cycle_num;

    // find the approximation
    retVal = prediction_PSEG_samples_struct_amp(approx_d, approx_mp, approx_prec, finalTime_d, finalTime_mp, currPrec, *cycle_num, 0, PSEG_samples, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);
  }
  else
  { // there was an error
    retVal = *cycle_num;
  }

  // loop to make sure we find an approximation or we run out of precision
  while (retVal)
  { // increase the precision
    if (currPrec < 64)
    { // increase to 64-bit precision
      currPrec = T->Precision = 64;
    }
    else
    { // increase by 32 bits
      currPrec = T->Precision = currPrec + 32;
    }
    // see if the precision is too high
    if (currPrec > T->AMP_max_prec)
    { // exit loop immediately!
      retVal = retVal_max_prec_reached;
      break;
    }

    // set the new precision
    initMP(T->Precision);
    change_prec(ED_mp, T->Precision);

    // refine the sample points to the current precision
    refine_PSEG_samples_struct_amp(PSEG_samples, start, finish, currPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

    // find an approximation based on the last 3 samples
    *cycle_num = find_cycle_num_amp(PSEG_samples, finish - 3, currPrec, T->power_series_sample_factor);
    // use this approximation as a way to determine where to stop the exhaustive search for finding the best cycle number
    *cycle_num = find_cycle_num_loop_amp(PSEG_samples, start, finish, currPrec, MAX(*cycle_num * 5, T->cycle_num_max), OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

    if (*cycle_num > 0)
    { // the structures are setup properly - find the prediction

      // store the cycle number to [start] - must be MP!
      PSEG_samples->samples_mp[start].cycle_num = *cycle_num;

      // find the approximation
      retVal = prediction_PSEG_samples_struct_amp(approx_d, approx_mp, approx_prec, finalTime_d, finalTime_mp, currPrec, *cycle_num, 0, PSEG_samples, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);
    }
    else
    { // there was an error
      retVal = *cycle_num;
    }
  }

  return retVal;
}

int prediction_PSEG_samples_struct_amp(vec_d approx_d, vec_mp approx_mp, int *approx_prec, comp_d finalT_d, comp_mp finalT_mp, int currPrec, int cycle_num, int setup_dX_Z, PSEG_samples_struct_amp *S, int start, int finish, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE: does a power series prediction based on samples to     *
* predict the value at finalT                                   *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - have an approx, otherwise, some error      *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0, i, n = finish - start - 1;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // determine whether to do prediction in double or multi precision
  if (currPrec < 64)
  { // do prediction in double precision
    _comp_d *s = (_comp_d *)bmalloc((n + 1) * sizeof(_comp_d));
    _point_d *dZ_ds = (_point_d *)bmalloc((n + 1) * sizeof(_point_d));

    for (i = 0; i <= n; i++)
      init_point_d(&dZ_ds[i], 0);

    // setup the structures
    retVal = setup_PSEG_prediction_amp(setup_dX_Z, currPrec, cycle_num, S, s, dZ_ds, NULL, NULL, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

    // check for error
    if (!retVal)
    { // interpolate
      *approx_prec = currPrec;
      if (n + 1 == 2)
        cubicInterp_d(approx_d, &s[0], &S->Z_rev_d[0], &dZ_ds[0], finalT_d);
      else
        hermiteInterpCW_d(approx_d, &s[0], &S->Z_rev_d[0], &dZ_ds[0], finalT_d, n + 1);
    }

    // free the memory
    for (i = 0; i <= n; i++)
      clear_point_d(&dZ_ds[i]);
    free(dZ_ds);
    free(s);
  }
  else
  { // do the prediction in multi precision
    _comp_mp *s = (_comp_mp *)bmalloc((n + 1) * sizeof(_comp_mp));
    _point_mp *dZ_ds = (_point_mp *)bmalloc((n + 1) * sizeof(_point_mp));

    // initialize
    for (i = 0; i <= n; i++)
    {
      init_mp2(&s[i], currPrec);
      init_point_mp2(&dZ_ds[i], 0, currPrec);
    }

    // setup the precision of approx_mp
    change_prec_vec_mp(approx_mp, currPrec);

    // setup the structures
    retVal = setup_PSEG_prediction_amp(setup_dX_Z, currPrec, cycle_num, S, NULL, NULL, s, dZ_ds, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

    // check for error
    if (!retVal)
    { // interpolate
      *approx_prec = currPrec;
      if (n + 1 == 2)
        cubicInterp_mp(approx_mp, &s[0], &S->Z_rev_mp[0], &dZ_ds[0], finalT_mp);
      else
        hermiteInterpCW_mp(approx_mp, &s[0], &S->Z_rev_mp[0], &dZ_ds[0], finalT_mp, n + 1);
    }

    // free the memory
    for (i = n; i >= 0; i--)
    {
      clear_mp(&s[i]);
      clear_point_mp(&dZ_ds[i]);
    }
    free(s);
    free(dZ_ds);
  }

  return retVal;
}

int find_cycle_num_amp(PSEG_samples_struct_amp *S, int start, int currPrec, double sample_factor)
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

  // error checking
  if (start + 2 > S->num_samples)
  {
    printf("ERROR: There needs to be atleast 3 good samples to provide a cycle number approximation!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // determine whether to compute in double or multi precision
  if (currPrec < 64)
  { // use double precision
    comp_d temp1_d, temp2_d, sum1_d, sum2_d, rand_d;

    // find the size
    size = S->samples_d[start].point->size;

    // find a random sum of the subtractions
    set_zero_d(sum1_d);
    set_zero_d(sum2_d);
    for (i = 0; i < size; i++)
    { // temp1 = [start+1] - [start]
      sub_d(temp1_d, &S->samples_d[start+1].point->coord[i], &S->samples_d[start].point->coord[i]);
      // temp2 = [start+2] - [start+1]
      sub_d(temp2_d, &S->samples_d[start+2].point->coord[i], &S->samples_d[start+1].point->coord[i]);

      // find random number
      get_comp_rand_d(rand_d);

      // add to sum1 & sum2
      sum_mul_d(sum1_d, rand_d, temp1_d);
      sum_mul_d(sum2_d, rand_d, temp2_d);
    }

    // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1 
    if (d_abs_d(sum1_d) > 0 && d_abs_d(sum2_d) > 0)
    { // sum2 / sum1
      div_d(temp1_d, sum2_d, sum1_d);

      // find the norm of sum2 / sum1
      c = d_abs_d(temp1_d);

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
  } 
  else
  { // use multi precision
    comp_mp temp1_mp, temp2_mp, sum1_mp, sum2_mp, rand_mp;
    // initialize
    init_mp(temp1_mp); init_mp(temp2_mp); init_mp(sum1_mp); init_mp(sum2_mp); init_mp(rand_mp);

    // find the size
    size = S->samples_mp[start].point->size;

    // find a random sum of the subtractions
    set_zero_mp(sum1_mp);
    set_zero_mp(sum2_mp);
    for (i = 0; i < size; i++)
    { // temp1 = [start+1] - [start]
      sub_mp(temp1_mp, &S->samples_mp[start+1].point->coord[i], &S->samples_mp[start].point->coord[i]);
      // temp2 = [start+2] - [start+1]
      sub_mp(temp2_mp, &S->samples_mp[start+2].point->coord[i], &S->samples_mp[start+1].point->coord[i]);

      // find a random number
      get_comp_rand_mp(rand_mp);

      // add to sum1 & sum2
      sum_mul_mp(sum1_mp, rand_mp, temp1_mp);
      sum_mul_mp(sum2_mp, rand_mp, temp2_mp);
    }

    // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1
    mpf_abs_mp(temp1_mp->r, sum1_mp);
    mpf_abs_mp(temp1_mp->i, sum2_mp);
    if (mpf_cmp_ui(temp1_mp->r, 0) > 0 && mpf_cmp_ui(temp1_mp->i, 0) > 0)
    { // sum2 / sum1
      div_mp(temp1_mp, sum2_mp, sum1_mp);

      // find the norm of sum2 / sum1
      c = d_abs_mp(temp1_mp);

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

    // clear
    clear_mp(temp1_mp); clear_mp(temp2_mp); clear_mp(sum1_mp); clear_mp(sum2_mp); clear_mp(rand_mp);
  }
  // double check retVal >= 1
  if (retVal < 1)
    retVal = 1;

  return retVal;
}

void refine_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S, int start, int finish, int prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: refines the samples points from start to < finish      *
* to prec - assume precision is already set appropriately in ED *
\***************************************************************/
{
  int i;

  // we do not want AMP checks done in refine - just refine the points the best that they can in the precision given!

  // refine the sample points from start to finish to the current precision
  if (prec < 64)
  { // make sure all of the sample points that we care about have a double precision refined value
    eval_struct_d e;
    init_eval_struct_d(e, 0, 0, 0);
    for (i = start; i < finish; i++)
      if (S->max_prec[i] < 64)
      { // sample point is already in double precision - just refine
        refine_d_basic(&S->samples_d[i], T, OUT, &e, ED_d, eval_func_d);
      }
      else
      { // copy MP to D and refine it
        convert_point_data_mp_to_d(&S->samples_d[i], &S->samples_mp[i]);
        refine_d_basic(&S->samples_d[i], T, OUT, &e, ED_d, eval_func_d);
      }
    clear_eval_struct_d(e);
  }
  else
  { // make sure all of the sample points that we care about have a 'prec' multi precision refined value
    eval_struct_mp e;
    init_eval_struct_mp2(e, 0, 0, 0, prec);

    for (i = start; i < finish; i++)
      if (S->max_prec[i] < 64)
      { // copy D to MP and refine it
        setprec_point_data_mp(&S->samples_mp[i], prec);
        convert_point_data_d_to_mp(&S->samples_mp[i], &S->samples_d[i]);
        refine_mp_basic(&S->samples_mp[i], T, OUT, &e, ED_mp, eval_func_mp);
        // set new maximum precision
        S->max_prec[i] = prec;
      }
      else if (S->max_prec[i] < prec)
      { // increase the precision and refine it
        change_prec_point_data_mp2(&S->samples_mp[i], prec);
        refine_mp_basic(&S->samples_mp[i], T, OUT, &e, ED_mp, eval_func_mp);
        // set new maximum precision
        S->max_prec[i] = prec;
      }
      else
      { // sample point has atleast 'prec' precision - just refine
        refine_mp_basic(&S->samples_mp[i], T, OUT, &e, ED_mp, eval_func_mp);
      }
    clear_eval_struct_mp(e);
  }

  return;
}

int setup_PSEG_prediction_amp(int setup_dX_Z, int currPrec, int cycle_num, PSEG_samples_struct_amp *S, _comp_d *s_d, _point_d *dZ_ds_d, _comp_mp *s_mp, _point_mp *dZ_ds_mp, int start, int finish, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE: setup the structures needed to do a prediction         *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - everything is setup, otherwise, some error *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0, i, j, k, rows, cols, size, n = finish - start - 1;

  // error checking
  if (start < 0 || finish > S->num_samples || n < 2)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // determine whether to setup in double or mutli precision
  if (currPrec < 64)
  { // setup in double precision
    double tempD;
    comp_d factor, tempComp, c_times_t;

    if (setup_dX_Z)
    { // setup dX_d & Z_rev_d
      eval_struct_d e;
      init_eval_struct_d(e, 0, 0, 0);

      for (i = start; i < finish; i++)
      { // evaluate the function
        eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples_d[i].point, S->samples_d[i].time, ED_d, eval_func_d);

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
        retVal = matrixSolve_d(S->dX_d[i], e.Jv, e.funcVals);
        if (retVal)
        {
          fprintf(OUT, "NOTE: matrixSolve has failed.\n");

          // return error
          return retVal;
        }

        // setup Z_rev_d
        k = finish - i - 1;
        increase_size_vec_d(&S->Z_rev_d[k], S->samples_d[i].point->size);
        size = S->Z_rev_d[k].size = S->samples_d[i].point->size;
        for (j = 0; j < size; j++)
        {
          set_d(&S->Z_rev_d[k].coord[j], &S->samples_d[i].point->coord[j]);
        }
      }
      clear_eval_struct_d(e);
    }

    // setup dZ_ds_d & s_d
    tempD = 1 / (double) cycle_num;
    // since normalize time, set s[n] == 1
    pow_rdouble_d(factor, S->samples_d[start].time, tempD);
    recip_d(factor, factor);
    set_double_d(&s_d[n], 1, 0);
    // find dZ_ds[n]
    mul_rdouble_d(c_times_t, S->samples_d[start].time, cycle_num);
    vec_mulcomp_d(&dZ_ds_d[n], S->dX_d[start], c_times_t);

    // loop over i = n - 1, .., 0 to find the dZ_ds[i] & s[i]
    for (i = n - 1; i >= 0; i--)
    { // find j
      j = finish - i - 1;
      // s_d[i] = t^(1/c) * factor
      pow_rdouble_d(&s_d[i], S->samples_d[j].time, tempD);
      mul_d(&s_d[i], &s_d[i], factor);
      // dZ_ds_d[i]
      k = cycle_num - 1;
      pow_rdouble_d(tempComp, &s_d[i], k);
      mul_d(tempComp, tempComp, c_times_t);
      vec_mulcomp_d(&dZ_ds_d[i], S->dX_d[j], tempComp);
    }
  }
  else
  { // setup in multi precision
    mpf_t tempMPF;
    comp_mp factor, tempComp, c_times_t;

    // initialize
    mpf_init(tempMPF);
    init_mp(factor); init_mp(tempComp); init_mp(c_times_t);

    if (setup_dX_Z)
    { // setup dX_mp & Z_rev_mp
      eval_struct_mp e;
      init_eval_struct_mp(e, 0, 0, 0);

      for (i = start; i < finish; i++)
      { // evaluate the function
        eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, S->samples_mp[i].point, S->samples_mp[i].time, ED_mp, eval_func_mp);

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
        retVal = matrixSolve_mp(S->dX_mp[i], e.Jv, e.funcVals);
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

        // setup Z_rev_mp
        k = finish - i - 1;
        increase_size_vec_mp(&S->Z_rev_mp[k], S->samples_mp[i].point->size);
        size = S->Z_rev_mp[k].size = S->samples_mp[i].point->size;
        for (j = 0; j < size; j++)
        {
          set_mp(&S->Z_rev_mp[k].coord[j], &S->samples_mp[i].point->coord[j]);
        }
      }
      clear_eval_struct_mp(e);
    }

    // setup dZ_ds_mp & s_mp
    mpf_set_ui(tempMPF, 1);
    mpf_div_ui(tempMPF, tempMPF, cycle_num);
    // since normalize time, set s[n] == 1
    pow_rmpf_mp(factor, S->samples_mp[start].time, tempMPF);
    recip_mp(factor, factor);
    mpf_set_ui(s_mp[n].r, 1); mpf_set_ui(s_mp[n].i, 0);
    // find dZ_ds[n]
    mpf_mul_ui(c_times_t->r, S->samples_mp[start].time->r, cycle_num);
    mpf_mul_ui(c_times_t->i, S->samples_mp[start].time->i, cycle_num);
    vec_mulcomp_mp(&dZ_ds_mp[n], S->dX_mp[start], c_times_t);

    // loop over i = 0, .., n - 1 to find the dZ_ds[i] & s[i]
    for (i = 0; i < n; i++)
    { // find j
      j = finish - i - 1;
      // s_mp[i] = t^(1/c) * factor
      pow_rmpf_mp(&s_mp[i], S->samples_mp[j].time, tempMPF);
      mul_mp(&s_mp[i], &s_mp[i], factor);

      // dZ_ds_mp[i]
      k = cycle_num - 1;
      pow_rdouble_mp(tempComp, &s_mp[i], k);
      mul_mp(tempComp, tempComp, c_times_t);
      vec_mulcomp_mp(&dZ_ds_mp[i], S->dX_mp[j], tempComp);
    }

    // free the memory
    mpf_clear(tempMPF);
    clear_mp(factor); clear_mp(tempComp); clear_mp(c_times_t);
  }
 
  return retVal;
}

int find_cycle_num_loop_amp(PSEG_samples_struct_amp *S, int start, int finish, int currPrec, int max_cycle_num, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - error, otherwise best_cycle_num           *
* NOTES: finds the number between 1 and max_cycle_num that best *
* uses {start,..,finish-1} to approximate finish                *
\***************************************************************/
{
  int c, j, size, retVal, best_c = 1, n = finish - start - 1;

  // error checking
  if (start < 0 || finish > S->num_samples || max_cycle_num < 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series prediction!\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (currPrec < 64)
  { // do in double precision
    double min_error = 1e300, error = 0, tempD;
    point_d approx;
    _comp_d *s = (_comp_d *)bmalloc((n + 1) * sizeof(_comp_d));
    _point_d *dZ_ds = (_point_d *)bmalloc((n + 1) * sizeof(_point_d));

    init_point_d(approx, 0);
    for (j = 0; j <= n; j++)
      init_point_d(&dZ_ds[j], 0);

    for (c = 1; c <= max_cycle_num; c++)
    { // setup the structures to do an interpolation
      retVal = setup_PSEG_prediction_amp((c == 1), currPrec, c, S, s, dZ_ds, NULL, NULL, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

      // check for errors
      if (retVal)
      { // clear structurs
        for (j = n; j >= 0; j--)
          clear_point_d(&dZ_ds[j]);
        free(dZ_ds);
        free(s);
        // return error
        return -1;
      }
        
      // interpolate from 1, .. n
      if (n == 2)
        cubicInterp_d(approx, &s[1], &S->Z_rev_d[1], &dZ_ds[1], &s[0]);
      else
        hermiteInterpCW_d(approx, &s[1], &S->Z_rev_d[1], &dZ_ds[1], &s[0], n);

      // find the error
      error = 0;
      size = approx->size;
      for (j = 0; j < size; j++)
      {
        sub_d(&approx->coord[j], &approx->coord[j], &S->samples_d[finish-1].point->coord[j]);
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
  }
  else
  { // do in multi precision
    mpf_t min_error, error;
    comp_mp tempComp;
    point_mp approx;
    _comp_mp *s = (_comp_mp *)bmalloc((n + 1) * sizeof(_comp_mp));
    _point_mp *dZ_ds = (_point_mp *)bmalloc((n + 1) * sizeof(_point_mp));

    mpf_init(min_error); mpf_init(error); 
    init_mp(tempComp);
    init_point_mp(approx, 0);
    for (j = 0; j <= n; j++)
    {
      init_mp(&s[j]);
      init_point_mp(&dZ_ds[j], 0);
    }

    // initialize min_error
    mpf_set_d(min_error, 1e300);

    for (c = 1; c <= max_cycle_num; c++)
    { // setup the structures to do an interpolation
      retVal = setup_PSEG_prediction_amp((c == 1), currPrec, c, S, NULL, NULL, s, dZ_ds, start, finish, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp);

      // check for errors
      if (retVal)
      { // clear MP
        mpf_clear(min_error); mpf_clear(error);
        clear_mp(tempComp);
        clear_point_mp(approx);
        for (j = 0; j <= n; j++)
        {
          clear_mp(&s[j]);
          clear_point_mp(&dZ_ds[j]);
        }
        free(s);
        free(dZ_ds);
        // return error
        return -1;
      }

      // interpolate from 1, .. n
      if (n == 2)
        cubicInterp_mp(approx, &s[1], &S->Z_rev_mp[1], &dZ_ds[1], &s[0]);
      else
        hermiteInterpCW_mp(approx, &s[1], &S->Z_rev_mp[1], &dZ_ds[1], &s[0], n);

      // find the error
      mpf_set_ui(error, 0);
      size = approx->size;
      for (j = 0; j < size; j++)
      {
        sub_mp(&approx->coord[j], &approx->coord[j], &S->samples_mp[finish-1].point->coord[j]);
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

    // clear MP
    mpf_clear(min_error); mpf_clear(error);
    clear_mp(tempComp);
    clear_point_mp(approx);
    for (j = 0; j <= n; j++)
    {
      clear_mp(&s[j]);
      clear_point_mp(&dZ_ds[j]);
    }
    free(s);
    free(dZ_ds);
  }

  return best_c;
}




