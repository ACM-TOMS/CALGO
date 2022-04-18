// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

#define NUM_AGREEMENT 3 // number of cycle agreement for 'pre-Cauchy endgame'
#define TIME_CUTOFF 1e-10 // time cutoff for cycle agreement for 'pre-Cauchy endgame'
#define FAIL_SAFE_MAX 250 // maximum number of loops before giving up

/////////////////////////////// CAUCHY ENDGAME //////////////////////////

// _d functions
int CauchyEG_preEG_d(point_d preEG_approx, int *cycle_num, point_data_d *Final, comp_d finalTime, point_data_d *Start, int num_agreement, double timeCutoff, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int CauchyEG_actual_d(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, point_d preEG_approx, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

// _mp functions
int CauchyEG_preEG_mp(point_mp preEG_approx, int *cycle_num, point_data_mp *Final, comp_mp finalTime, point_data_mp *Start, int num_agreement, double timeCutoff, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int CauchyEG_actual_mp(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, point_mp preEG_approx, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// _amp functions
int CauchyEG_preEG_amp(point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, int *cycle_num, point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, int num_agreement, double timeCutoff, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int CauchyEG_actual_amp(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, point_d preEG_approx_d, point_mp preEG_approx_mp, int preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

/////// Double Precision ////////

int CauchyEG_d(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in double precision       *
\***************************************************************/
{
  int i, retVal, cycle, samples_per_loop;
  comp_d finalTime;
  point_data_d *endSamples = NULL;

  // setup finalTime
  set_double_d(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  retVal = CauchyEG_d2(pathNum, Final, last_approx, &endSamples, &cycle, &samples_per_loop, finalTime, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  // clear endSamples since is it not being used
  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // clear
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_d(&endSamples[i]);
    free(endSamples);
  }

  return retVal;
}

int CauchyEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d **endSamples, int *cycle, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endSamples, cycle, samples_per_loop            *
* NOTES: does Cauchy integral endgame in double precision       *
\***************************************************************/
{
  int i, retVal = 0;
  comp_d endTime;

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_d(endTime, T->endgameBoundary, 0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0;

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // track to the endgame boundary
  retVal = track_d(Final, Start, endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < Final->point->size; i++)
    fprintf(midOUT, "%.15e %.15e\n", Final->point->coord[i].r, Final->point->coord[i].i);
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check for success
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Cauchy Integral Endgame never started!\n");

    T->t_val_at_latest_sample_point = Final->time->r;
  }
  else
  { // success - so run the acutal endgame

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // run the endgame
    retVal = CauchyEG_main_d2(Final, last_approx, endSamples, cycle, samples_per_loop, finalTime, Final, T, OUT, ED, eval_func, find_dehom);
  }

  return retVal;
}

int CauchyEG_main_d(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in double precision            *
\***************************************************************/
{
  int retVal = 0, num_agreement = NUM_AGREEMENT;
  double timeCutoff = TIME_CUTOFF;
  eval_struct_d e;

  init_eval_struct_d(e, 0, 0, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // perform the 'pre-Cauchy EG' - use cycle number approximations to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_d(last_approx, cycle_num, Final, finalTime, Start, num_agreement, timeCutoff, T, OUT, &e, ED, eval_func);

  if (retVal == 0) 
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_d(Final, last_approx, endLoopSamples, cycle_num, samples_per_loop, last_approx, finalTime, Final, T, OUT, &e, ED, eval_func);
  }
 
  clear_eval_struct_d(e);
 
  return retVal;
}

int CauchyEG_preEG_d(point_d preEG_approx, int *cycle_num, point_data_d *Final, comp_d finalTime, point_data_d *Start, int num_agreement, double timeCutoff, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  *
* can efficiently start using the actual Cauchy EG.             *
* Uses timeCutoff as a fail-safe stopping criterion.            *
\***************************************************************/
{
  int i, last_sample = 0, retVal = 0, agreement = 0, *cycle = (int *)bmalloc(num_agreement * sizeof(int));
  point_data_d *samples_d = (point_data_d *)bmalloc(3 * sizeof(point_data_d)); // only 3 sample points are needed for a cycle number approximation
  comp_d endTime;

  // initialize samples_d
  for (i = 0; i < 3; i++)
    init_point_data_d(&samples_d[i], 0);

  // initialize cycle
  *cycle_num = Final->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    cycle[i] = -1;

  // setup the first sample
  point_data_cp_d(&samples_d[0], Start);
  // refine the sample point
  refine_d_basic(&samples_d[0], T, OUT, e_d, ED, eval_func);

  // find the first 3 samples
  for (i = 1; i < 3 && retVal == 0; i++)
  { // find the end time for the next sample
    mul_rdouble_d(endTime, samples_d[i-1].time, T->power_series_sample_factor);

    // track to the next sample point
    retVal = track_d(&samples_d[i], &samples_d[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (!retVal)
    { // refine the sample point
      refine_d_basic(&samples_d[i], T, OUT, e_d, ED, eval_func);
    }
    else
    { // store where this error occured
      last_sample = i;
    }
  }

  // see if we should continue
  if (!retVal)
  { // store where the last sample will be
    last_sample = 2;

    // find the first approximation for the cycle number
    cycle[0] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);

    // loop to find the other cycle numbers
    i = 1;
    while (retVal == 0 && i < num_agreement && d_abs_d(samples_d[2].time) > timeCutoff)
    { // update samples_d
      point_data_cp_d(&samples_d[0], &samples_d[1]);
      point_data_cp_d(&samples_d[1], &samples_d[2]);

      // find the end time for the next sample
      mul_rdouble_d(endTime, samples_d[1].time, T->power_series_sample_factor);

      // track to the next sample point
      retVal = track_d(&samples_d[2], &samples_d[1], endTime, T, OUT, ED, eval_func);

      // check for errors
      if (!retVal)
      { // refine the sample point
        refine_d_basic(&samples_d[2], T, OUT, e_d, ED, eval_func);

        // find the next approximation for the cycle number
        cycle[i] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);
      }

      // increment i;
      i++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = 1;
      for (i = 1; i < num_agreement && agreement; i++)
        if (cycle[0] != cycle[i])
          agreement = 0;

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && d_abs_d(samples_d[2].time) > timeCutoff)
      { // update cycle
        for (i = 1; i < num_agreement; i++)
          cycle[i-1] = cycle[i];

        // update samples_d
        for (i = 1; i < 3; i++)
          point_data_cp_d(&samples_d[i-1], &samples_d[i]);

        // find the end time for the next sample
        mul_rdouble_d(endTime, samples_d[1].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_d(&samples_d[2], &samples_d[1], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_d_basic(&samples_d[2], T, OUT, e_d, ED, eval_func);

          // find the next approximation for the cycle number
          cycle[num_agreement - 1] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);

          // check for agreement
          agreement = 1;
          for (i = 1; i < num_agreement && agreement; i++)
            if (cycle[0] != cycle[i])
              agreement = 0;
        }
      }
    }
  }
 
  // check for errors
  if (retVal == 0)
  { // find the power series approximation described by the 3 sample points
    find_Cauchy_preEG_approx_d(preEG_approx, finalTime, samples_d, 3, cycle[num_agreement - 1], OUT, e_d, ED, eval_func);
  }

  // update Final
  point_data_cp_d(Final, &samples_d[last_sample]);
  *cycle_num = Final->cycle_num = cycle[num_agreement - 1];
  T->t_val_at_latest_sample_point = samples_d[last_sample].time->r;

  // clear memory
  for (i = 2; i >= 0; i--)
    clear_point_data_d(&samples_d[i]);
  free(samples_d);
  free(cycle);

  return retVal;
}

void find_Cauchy_preEG_approx_d(point_d preEG_approx, comp_d finalTime, point_data_d *samples_d, int num_samples, int cycle_num, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the approximation for Cauchy preEG               *
\***************************************************************/
{ // error checking
  if (cycle_num <= 0 || num_samples <= 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series predicition!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // cycle_num >= 1 && num_samples >= 2
  int i, j, k, rows, cols, retVal = 0;
  double tempD = 1 / (double) cycle_num;
  comp_d factor, tempComp, c_times_t;

  // allocate s, Z & dZ_ds
  _comp_d *s = (_comp_d *)bmalloc(num_samples * sizeof(_comp_d));
  _point_d *Z = (_point_d *)bmalloc(num_samples * sizeof(_point_d)), *dZ_ds = (_point_d *)bmalloc(num_samples * sizeof(_point_d));

  // initialize
  for (i = 0; i < num_samples; i++)
  {
    init_point_d(&Z[i], 0);
    init_point_d(&dZ_ds[i], 0);
  }

  // normalize time (s[n-1] = 1)
  pow_rdouble_d(factor, samples_d[0].time, tempD);
  recip_d(factor, factor);
  set_one_d(&s[num_samples - 1]);
 
  // setup c_times_t
  mul_rdouble_d(c_times_t, samples_d[0].time, cycle_num);

  // setup the structures
  for (i = 0; i < num_samples && retVal == 0; i++)
  { // evaluate the function
    eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, samples_d[i].point, samples_d[i].time, ED, eval_func);

    // find Y = -dH/dt = -dH/dp * dp/dt
    rows = e_d->funcVals->size; // == Jp->rows
    cols = e_d->Jp->cols;
    for (k = 0; k < rows; k++)
    {
      set_zero_d(&e_d->funcVals->coord[k]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&e_d->funcVals->coord[k], &e_d->Jp->entry[k][j], &e_d->parDer->coord[j]);
      }
      neg_d(&e_d->funcVals->coord[k], &e_d->funcVals->coord[k]);
    }

    // find dZ = (dH/dZ)^-1 * (-dH/dt)
    retVal = matrixSolve_d(e_d->parVals, e_d->Jv, e_d->funcVals);
   
    // check for errors
    if (!retVal)
    { // find the index
      k = num_samples - i - 1;

      // setup the structures - note that s[n-1] == 1 is already setup
      if (k != num_samples - 1)
      { // s[k] = t^(1/c) * factor
        pow_rdouble_d(&s[k], samples_d[i].time, tempD);
        mul_d(&s[k], &s[k], factor);
      
        // setup tempComp = s[k]^(c-1) * c * t
        j = cycle_num - 1;
        pow_rdouble_d(tempComp, &s[k], j);
        mul_d(tempComp, tempComp, c_times_t);
      }
      else
      { // setup tempComp
        set_d(tempComp, c_times_t);
      }

      // setup Z[k] & dZ_ds[k]
      rows = samples_d[i].point->size;
      change_size_point_d(&Z[k], rows);
      change_size_point_d(&dZ_ds[k], rows);
      Z[k].size = dZ_ds[k].size = rows;
      for (j = 0; j < rows; j++)
      {
        set_d(&Z[k].coord[j], &samples_d[i].point->coord[j]);
        mul_d(&dZ_ds[k].coord[j], &e_d->parVals->coord[j], tempComp);
      }
    }
  }

  // check for errors
  if (!retVal)   
  { // find the approximation
    if (num_samples == 2)
      cubicInterp_d(preEG_approx, s, Z, dZ_ds, finalTime);
    else
      hermiteInterpCW_d(preEG_approx, s, Z, dZ_ds, finalTime, num_samples);
  }
  else
  { // just copy the next to last sample point and call it good enough - Cauchy will take over from the last sample point and do its thing!
    point_cp_d(preEG_approx, samples_d[num_samples - 2].point);
  }

  // clear
  for (i = num_samples - 1; i >= 0; i--)
  {
    clear_point_d(&Z[i]);
    clear_point_d(&dZ_ds[i]);
  }
  free(s);
  free(Z);
  free(dZ_ds);

  return;
}

int CauchyEG_actual_d(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, point_d preEG_approx, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual Cauchy integral endgame in double prec *
\***************************************************************/
{
  int i, retVal = 0, curr_cycle_size = 1;
  double closed_loop_tol, min_closed_loop_tol = 1e-10, max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult), approx_err;
  comp_d endTime;
  point_d new_approx;
  point_data_d *Samples = NULL;

  init_point_d(new_approx, 0);

  // setup the number of samples per loop
  *samples_per_loop = T->num_PSEG_sample_points + 1;

  // initially allocate for Samples
  Samples = (point_data_d *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));
  for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    init_point_data_d(&Samples[i], 0);

  // setup last_approx to be preEG_approx
  point_cp_d(last_approx, preEG_approx);

  // find the tolerance to close the first loop
  closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, Start, NULL, 52, T, e_d, NULL, ED, NULL, eval_func, NULL, NULL);

  // find the first set of samples
  retVal = find_Cauchy_samples_d(closed_loop_tol, &Samples, cycle_num, *samples_per_loop, &curr_cycle_size, Start, T, OUT, e_d, ED, eval_func);

  // see if we can continue
  if (!retVal)
  { // find new_approx
    find_Cauchy_approx_d(new_approx, Samples, *cycle_num, *samples_per_loop);

    // main loop - keep tracking until we either have failure or we converge
    do 
    { // compare the approximations
      T->error_at_latest_sample_point = approx_err = compare_approximations_d(new_approx, last_approx);
fprintf(OUT, "cycle: %d error: %e\n", *cycle_num, approx_err);

      // check for convergence
      if (approx_err < T->final_tolerance)
      { // we accept our approximation

        // setup Final
        point_cp_d(Final->point, new_approx);
        set_d(Final->time, finalTime);
        Final->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = Samples[0].time->r;

        retVal = 0;
        break;
      }
      else if (d_abs_d(Samples[0].time) < T->minTrackT)
      { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

        // setup Final
        point_cp_d(Final->point, new_approx);
        set_d(Final->time, Samples[0].time);
        Final->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = Samples[0].time->r;

        retVal = retVal_EG_failed_to_converge;
        break;
      }
      else
      { // copy new_approx to last_approx and find the next approximation
        point_cp_d(last_approx, new_approx);

        // find end time for the next sample point
        mul_rdouble_d(endTime, Samples[0].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_d(&Samples[0], &Samples[0], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_d_basic(&Samples[0], T, OUT, e_d, ED, eval_func);

          // find the tolerance to close the next loop
          closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &Samples[0], NULL, 52, T, e_d, NULL, ED, NULL, eval_func, NULL, NULL);

          // find the next set of samples
          retVal = find_Cauchy_samples_d(closed_loop_tol, &Samples, cycle_num, *samples_per_loop, &curr_cycle_size, &Samples[0], T, OUT, e_d, ED, eval_func);

          // check for errors
          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_d(new_approx, Samples, *cycle_num, *samples_per_loop);
          }
        }
      }
    } while (!retVal);
  }

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    point_cp_d(Final->point, Samples[0].point);
    set_d(Final->time, Samples[0].time);
    Final->cycle_num = *cycle_num;

    T->t_val_at_latest_sample_point = Samples[0].time->r;
  }
  else
  { // allocate endLoopSamples
    *endLoopSamples = (point_data_d *)bmalloc(*samples_per_loop * (*cycle_num) * sizeof(point_data_d));

    // setup endLoopSamples
    for (i = *samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
    {
      init_point_data_d(&(*endLoopSamples)[i], Samples[i].point->size);
      point_data_cp_d(&(*endLoopSamples)[i], &Samples[i]);
    }
  }

  // release memory
  for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_d(&Samples[i]);
  free(Samples);
  clear_point_d(new_approx);

  return retVal;
}

int circle_track_d(point_data_d Final[], int Final_digits[], point_data_d *Start, int M, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks around a circle to find the M+1 sample points   *
\***************************************************************/
{
  int i, cont, retVal = 0, consecSuccess, numSteps = 0;
  double absValOfDistanceLeft, tempD, radius, next_step_angle = 0, next_stop_angle = 0, curr_angle = -2 * M_PI;
  comp_d nextTime;

  // error checking - want M >= 3 so that we do not have to predict directly across the origin
  if (M < 3)
  {
    printf("ERROR: The number of sample points (%d) for circle tracking must be >= 3!\n", M);
    bexit(ERROR_CONFIGURATION);
  }

  // find the radius of the circle that we are going to track around
  radius = d_abs_d(Start->time);

  // error checking
  if (radius <= 0)
  {
    printf("ERROR: The radius of the circle (%e) needs to be positive!\n", radius);
    bexit(ERROR_CONFIGURATION);
  }

  // if this is the first step of the path, initialize the current step size as the maximum step size - otherwise, it is already set!
  if (T->first_step_of_path)
  {
    T->currentStepSize = T->maxStepSize;
    T->first_step_of_path = 0;
  }

  // sample point [0] is Start
  point_data_cp_d(&Final[0], Start);
  refine_d_basic(&Final[0], T, OUT, e, ED, eval_func);
  Final_digits[0] = residual_to_digits_d(T->latest_newton_residual_d);

  // loop over each sample point
  consecSuccess = 0;
  for (i = 1; i <= M && retVal == 0; i++)
  { // setup the starting point
    point_cp_d(Final[i].point, Final[i-1].point);
    set_d(Final[i].time, Final[i-1].time);

    // reset the number of steps
    numSteps = 0;

    // calculate the next angle to stop at = -2 PI + 2 PI i / M
    next_stop_angle = 2 * M_PI * (i * 1.0 / M - 1);

    // track to the next sample point
    cont = 1;
    while (cont)
    { // print out the current time
      if (T->outputLevel > 0)
      {
        fprintf(OUT, "Time = %.15e %.15e\n", Final[i].time->r, Final[i].time->i);
        if (T->screenOut)
          printf("Time = %.15e %.15e\n", Final[i].time->r, Final[i].time->i);
      }

      // find the distance to the end = radius * (next_stop_angle - curr_angle)
      absValOfDistanceLeft = (next_stop_angle - curr_angle) * radius;

      // check to see how close we are to the target
      if (next_stop_angle == curr_angle)
      { // we are exactly at the end!
        cont = 0;
      }
      else if (absValOfDistanceLeft < T->currentStepSize)
      { // we can take a smaller step than our current step size to hit the target time
        next_step_angle = next_stop_angle;
      }
      else
      { // the distance left is larger than the current step size, so we head in the correct direction will a step size of the current step size
        tempD = T->currentStepSize / radius; // calculate the change in angle that we can do
        next_step_angle = curr_angle + tempD; // calculate the next step angle
      }

      // check to see if we should continue
      if (cont)
      { // setup nextTime
        set_double_d(nextTime, radius * cos(next_step_angle), radius * sin(next_step_angle));

        // take a step
        retVal = step_d(&Final[i], nextTime, T, OUT, e, ED, eval_func);
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
          // update curr_angle
          curr_angle = next_step_angle;

          // check to see if we increase the step size
          if (consecSuccess >= T->cSecInc)
          { // find the new step size if we increase it
            tempD = T->step_success_factor * T->currentStepSize;
            // the new step size is the minimum of the max step size and tempD
            T->currentStepSize = MIN(T->maxStepSize, tempD);
            consecSuccess = 0;
          }
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
          }
          else if (retVal == retVal_going_to_infinity)
          { // path is going to infinity
            retVal = retVal_going_to_infinity;
            cont = 0;
          }
          else if (retVal == retVal_higher_prec_needed)
          { // higher precision is needed
            retVal = retVal_higher_prec_needed;
            cont = 0;
          }
        }

        // check to see if we have taken too many steps
        if (numSteps > T->maxNumSteps)
        {
          retVal = retVal_too_many_steps;
          cont = 0;
        }
      }
    }

    // display how many steps it took
    if (T->outputLevel >= 0)
    {
      fprintf(OUT, "numSteps = %d\n", numSteps);
      if (T->screenOut)
        printf("numSteps = %d\n", numSteps);
    }

    // refine the sample, if it was good
    if (retVal == 0)
    {
      refine_d_basic(&Final[i], T, OUT, e, ED, eval_func);
      Final_digits[i] = residual_to_digits_d(T->latest_newton_residual_d);
    }
  }

  return retVal;
}

int find_Cauchy_samples_d(double closed_loop_max_error, point_data_d **Samples, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_d *Start, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Populates Samples and finds the cycle number to close  *
* the loop                                                      *
\***************************************************************/
{
  int i, tempInt, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  comp_d time;

  // save the correct time
  set_d(time, Start->time);

  // setup Samples[0] - circle_track_d will refine this sample
  point_data_cp_d(&(*Samples)[0], Start);

  // loop until get back to the original point
  *cycle_num = 0;
  while (1)
  { // go around the origin once and refine the samples as they are found
    retVal = circle_track_d(&(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], &(*Samples)[*cycle_num * num_samples], num_samples, T, OUT, e_d, ED, eval_func);

    // increment the cycle_num
    (*cycle_num)++;
    
    if (retVal) 
    { // return error
      if (retVal != retVal_higher_prec_needed)
        fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
      break;
    }
    else
    { // check to see if we have closed the loop
      refine_d_basic(&(*Samples)[0], T, OUT, e_d, ED, eval_func);
      refine_d_basic(&(*Samples)[*cycle_num * num_samples], T, OUT, e_d, ED, eval_func);
      if (check_closed_loop_d(closed_loop_max_error, &(*Samples)[0], &s_digits[0], &(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time, e_d, T, OUT, ED, eval_func))
      { // error is small enough - exit loop with success
        retVal = 0;
        break;
      }
      else if (*cycle_num > fail_safe_cycle_max)
      { // too many iterations 
        fprintf(OUT, "ERROR: Cycle number too high to detect!\n");
        retVal = retVal_cycle_num_too_high;
        break;
      }
      else if (*cycle_num >= *curr_cycle_size)
      { // need to increase the size of memory
        s_digits = (int *)brealloc(s_digits, (2 * (*curr_cycle_size) * num_samples + 1) * sizeof(int));

        point_data_d *temp_ptr = *Samples;
        *Samples = NULL;
        tempInt = *curr_cycle_size * num_samples;
        *Samples = (point_data_d *)bmalloc((2 * (*curr_cycle_size) * num_samples + 1) * sizeof(point_data_d));
        for (i = 2 * (*curr_cycle_size) * num_samples; i >= 0; i--)
        {
          init_point_data_d(&(*Samples)[i], 0);

          if (i <= tempInt)
          { // copy temp_ptr to Samples
            point_data_cp_d(&(*Samples)[i], &temp_ptr[i]);
            clear_point_data_d(&temp_ptr[i]);
          }
        }
        free(temp_ptr);
        temp_ptr = NULL;
        // update curr_cycle_size
        *curr_cycle_size = 2 * (*curr_cycle_size);
      }
    }
  }

  // release memory
  free(s_digits);

  return retVal;
}

void find_Cauchy_approx_d(point_d approx, point_data_d Samples[], int cycle_num, int num_samples)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Finds the Cauchy integral approximation using Trapezoid*
* rule based on cycle_num & num_samples                         *
\***************************************************************/
{
  int i, j, size, sample_num = cycle_num * num_samples; // the first and last sample are equal so we skip the last one

  // set size
  size = Samples[0].point->size;
  change_size_point_d(approx, size);
  approx->size = size;

  for (i = 0; i < size; i++)
  { // initialize to 0
    set_zero_d(&approx->coord[i]);

    // sum over all of the samples
    for (j = 0; j < sample_num; j++)
    {
      add_d(&approx->coord[i], &approx->coord[i], &Samples[j].point->coord[i]);
    }

    // divide by number of samples added together
    approx->coord[i].r /= sample_num;
    approx->coord[i].i /= sample_num;
  }

  return;
}

int check_closed_loop_d(double closed_loop_max_error, point_data_d *Pt0, int *digitsCorrect0, point_data_d *Pt1, int *digitsCorrect1, comp_d time, eval_struct_d *e_d, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - Pt0 == Pt1, 0 - otherwise                  *
* NOTES: check to see if we have closed the loop for Cauchy EG  *
\***************************************************************/
{
  int retVal = 0, min_digits = 0;

  // setup min_digits
  min_digits = ceil(-log10(closed_loop_max_error) + 0.5); // minimum number of digits that we need to have correct for each point

  // check to see if we have the correct number of digits
  if (*digitsCorrect0 < min_digits)
  { // need to refine Pt0
    refine_digits_d(T->outputLevel, min_digits, &T->latest_newton_residual_d, *digitsCorrect0, Pt0, Pt0, time, OUT, e_d, ED, eval_func);

    // update digitsCorrect0
    *digitsCorrect0 = min_digits;
  }

  if (*digitsCorrect1 < min_digits)
  { // need to refine Pt1
    refine_digits_d(T->outputLevel, min_digits, &T->latest_newton_residual_d, *digitsCorrect0, Pt0, Pt0, time, OUT, e_d, ED, eval_func);

    // update digitsCorrect1
    *digitsCorrect1 = min_digits;
  }

  // check to see if we have closed the loop
  retVal = isSamePoint(Pt0->point, NULL, 52, Pt1->point, NULL, 52, closed_loop_max_error);

  return retVal;
}

////// Multi Precision ///////

int CauchyEG_mp(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in fixed multi precision  *
\***************************************************************/
{
  int i, retVal, cycle, samples_per_loop;
  comp_mp finalTime;
  point_data_mp *endSamples = NULL;

  init_mp(finalTime);

  // setup finalTime
  set_double_mp(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  retVal = CauchyEG_mp2(pathNum, Final, last_approx, &endSamples, &cycle, &samples_per_loop, finalTime, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  // clear endSamples since is it not being used
  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // clear
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_mp(&endSamples[i]);
    free(endSamples);
  }

  clear_mp(finalTime);

  return retVal;
}

int CauchyEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp **endSamples, int *cycle, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endSamples, cycle, samples_per_loop            *
* NOTES: does Cauchy integral endgame in fixed multi precision  *
\***************************************************************/
{ 
  int i, retVal = 0;
  comp_mp endTime;
  
  // set everything to the current precision
  initMP(T->Precision);
  
  // initialize MP
  init_mp(endTime);
  
  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_mp(endTime, T->endgameBoundary, 0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0;
  
  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // track to the endgame boundary
  retVal = track_mp(Final, Start, endTime, T, OUT, ED, eval_func);
  
  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < Final->point->size; i++)
  {
    print_mp(midOUT, 0, &Final->point->coord[i]);
    fprintf(midOUT, "\n");
  } 
  fprintf(midOUT, "\n");
    
  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check for success
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Cauchy Integral Endgame never started!\n");

    T->t_val_at_latest_sample_point = mpf_get_d(Final->time->r);
  }
  else
  { // success - so run the acutal endgame

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // run the endgame
    retVal = CauchyEG_main_mp2(Final, last_approx, endSamples, cycle, samples_per_loop, finalTime, Final, T, OUT, ED, eval_func, find_dehom);
  }

  // clear MP
  clear_mp(endTime);

  return retVal;
}

int CauchyEG_main_mp(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in multi precision             *
\***************************************************************/
{
  int retVal = 0, num_agreement = NUM_AGREEMENT;
  double timeCutoff = TIME_CUTOFF;
  eval_struct_mp e;

  // initialize
  init_eval_struct_mp(e, 0, 0, 0);

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // perform the 'pre-Cauchy EG' - use cycle number approximations to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_mp(last_approx, cycle_num, Final, finalTime, Start, num_agreement, timeCutoff, T, OUT, &e, ED, eval_func);

  if (retVal == 0)
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_mp(Final, last_approx, endLoopSamples, cycle_num, samples_per_loop, last_approx, finalTime, Final, T, OUT, &e, ED, eval_func);
  }

  // clear
  clear_eval_struct_mp(e);

  return retVal;
}

int CauchyEG_preEG_mp(point_mp preEG_approx, int *cycle_num, point_data_mp *Final, comp_mp finalTime, point_data_mp *Start, int num_agreement, double timeCutoff, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  * 
* can efficiently start using the actual Cauchy EG.             * 
* Uses timeCutoff as a fail-safe stopping criterion.            * 
\***************************************************************/
{
  int i, last_sample = 0, retVal = 0, agreement = 0, *cycle = (int *)bmalloc(num_agreement * sizeof(int));
  point_data_mp *samples = (point_data_mp *)bmalloc(3 * sizeof(point_data_mp)); // only 3 sample points are needed for a cycle number approximation
  comp_mp endTime;

  init_mp(endTime);
  for (i = 0; i < 3; i++)
    init_point_data_mp(&samples[i], 0);

  // initialize cycle
  *cycle_num = Final->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    cycle[i] = -1;

  // setup the first sample
  point_data_cp_mp(&samples[0], Start);
  // refine the sample point
  refine_mp_basic(&samples[0], T, OUT, e_mp, ED, eval_func);

  // find the first 3 samples
  for (i = 1; i < 3 && retVal == 0; i++)
  { // find the end time for the next sample
    mul_rdouble_mp(endTime, samples[i-1].time, T->power_series_sample_factor);

    // track to the next sample point
    retVal = track_mp(&samples[i], &samples[i-1], endTime, T, OUT, ED, eval_func);

    // check for errors
    if (!retVal)
    { // refine the sample point
      refine_mp_basic(&samples[i], T, OUT, e_mp, ED, eval_func);
    }
    else
    { // store where this error occured
      last_sample = i;
    }
  }

  // see if we should continue
  if (!retVal)
  { // store where the last sample will be
    last_sample = 2;

    // find the first approximation for the cycle number
    cycle[0] = find_cycle_num_mp(samples, 0, T->power_series_sample_factor);

    // loop to find the other cycle numbers
    i = 1;
    while (retVal == 0 && i < num_agreement && d_abs_mp(samples[2].time) > timeCutoff)
    { // update samples
      point_data_cp_mp(&samples[0], &samples[1]);
      point_data_cp_mp(&samples[1], &samples[2]);

      // find the end time for the next sample
      mul_rdouble_mp(endTime, samples[1].time, T->power_series_sample_factor);

      // track to the next sample point
      retVal = track_mp(&samples[2], &samples[1], endTime, T, OUT, ED, eval_func);

      // check for errors
      if (!retVal)
      { // refine the sample point
        refine_mp_basic(&samples[2], T, OUT, e_mp, ED, eval_func);

        // find the next approximation for the cycle number
        cycle[i] = find_cycle_num_mp(samples, 0, T->power_series_sample_factor);
      }

      // increment i;
      i++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = 1;
      for (i = 1; i < num_agreement && agreement; i++)
        if (cycle[0] != cycle[i])
          agreement = 0;

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && d_abs_mp(samples[2].time) > timeCutoff)
      { // update cycle
        for (i = 1; i < num_agreement; i++)
          cycle[i-1] = cycle[i];

        // update samples
        for (i = 1; i < 3; i++)
          point_data_cp_mp(&samples[i-1], &samples[i]);

        // find the end time for the next sample
        mul_rdouble_mp(endTime, samples[1].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_mp(&samples[2], &samples[1], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_mp_basic(&samples[2], T, OUT, e_mp, ED, eval_func);

          // find the next approximation for the cycle number
          cycle[num_agreement - 1] = find_cycle_num_mp(samples, 0, T->power_series_sample_factor);

          // check for agreement
          agreement = 1;
          for (i = 1; i < num_agreement && agreement; i++)
            if (cycle[0] != cycle[i])
              agreement = 0;
        }
      }
    }
  }

  // check for errors
  if (retVal == 0)
  { // find the power series approximation described by the 3 sample points
    find_Cauchy_preEG_approx_mp(preEG_approx, finalTime, samples, 3, cycle[num_agreement - 1], OUT, e_mp, ED, eval_func);
  }

  // update Final
  point_data_cp_mp(Final, &samples[last_sample]);
  *cycle_num = Final->cycle_num = cycle[num_agreement - 1];
  T->t_val_at_latest_sample_point = mpf_get_d(samples[last_sample].time->r);

  // clear memory
  clear_mp(endTime);
  for (i = 2; i >= 0; i--)
    clear_point_data_mp(&samples[i]);
  free(samples);
  free(cycle);

  return retVal;
}

void find_Cauchy_preEG_approx_mp(point_mp preEG_approx, comp_mp finalTime, point_data_mp *samples_mp, int num_samples, int cycle_num, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the approximation for Cauchy preEG               *
\***************************************************************/
{ // error checking
  if (cycle_num <= 0 || num_samples <= 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series predicition!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // cycle_num >= 1 && num_samples >= 2
  int i, j, k, rows, cols, retVal = 0;
  mpf_t tempMPF;
  comp_mp factor, tempComp, c_times_t;

  // allocate s, Z & dZ_ds
  _comp_mp *s = (_comp_mp *)bmalloc(num_samples * sizeof(_comp_mp));
  _point_mp *Z = (_point_mp *)bmalloc(num_samples * sizeof(_point_mp)), *dZ_ds = (_point_mp *)bmalloc(num_samples * sizeof(_point_mp));

  // initialize
  mpf_init(tempMPF);
  init_mp(factor); init_mp(tempComp); init_mp(c_times_t);
  for (i = 0; i < num_samples; i++)
  {
    init_mp(&s[i]);
    init_point_mp(&Z[i], 0); 
    init_point_mp(&dZ_ds[i], 0);
  }

  // normalize time (s[n-1] = 1)
  mpf_set_ui(tempMPF, cycle_num);
  mpf_ui_div(tempMPF, 1, tempMPF);
  pow_rmpf_mp(factor, samples_mp[0].time, tempMPF);
  recip_mp(factor, factor);
  set_one_mp(&s[num_samples - 1]);

  // setup c_times_t
  mpf_mul_ui(c_times_t->r, samples_mp[0].time->r, cycle_num);
  mpf_mul_ui(c_times_t->i, samples_mp[0].time->i, cycle_num);

  // setup the structures
  for (i = 0; i < num_samples && retVal == 0; i++)
  { // evaluate the function
    eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, samples_mp[i].point, samples_mp[i].time, ED, eval_func);

    // find Y = -dH/dt = -dH/dp * dp/dt
    rows = e_mp->funcVals->size; // == Jp->rows;
    cols = e_mp->Jp->cols;
    for (k = 0; k < rows; k++)
    {
      set_zero_mp(&e_mp->funcVals->coord[k]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&e_mp->funcVals->coord[k], &e_mp->Jp->entry[k][j], &e_mp->parDer->coord[j]);
      }
      neg_mp(&e_mp->funcVals->coord[k], &e_mp->funcVals->coord[k]);
    }

    // find dZ = (dH/dZ)^-1 * (-dH/dt)
    retVal = matrixSolve_mp(e_mp->parVals, e_mp->Jv, e_mp->funcVals);

    // check for errors
    if (!retVal)
    { // find the index
      k = num_samples - i - 1;

      // setup the structures - note that s[n-1] == 1 is already setup
      if (k != num_samples - 1)
      { // s[k] = t^(1/c) * factor
        pow_rmpf_mp(&s[k], samples_mp[i].time, tempMPF);
        mul_mp(&s[k], &s[k], factor);

        // setup tempComp = s[k]^(c-1) * c * t
        j = cycle_num - 1;
        pow_rdouble_mp(tempComp, &s[k], j);
        mul_mp(tempComp, tempComp, c_times_t);
      }
      else
      { // setup tempComp
        set_mp(tempComp, c_times_t);
      }

      // setup Z[k] & dZ_ds[k]
      rows = samples_mp[i].point->size;
      change_size_point_mp(&Z[k], rows);
      change_size_point_mp(&dZ_ds[k], rows);
      Z[k].size = dZ_ds[k].size = rows;
      for (j = 0; j < rows; j++)
      {
        set_mp(&Z[k].coord[j], &samples_mp[i].point->coord[j]);
        mul_mp(&dZ_ds[k].coord[j], &e_mp->parVals->coord[j], tempComp);
      }
    }
  }

  // check for errors
  if (!retVal)
  { // find the approximation
    if (num_samples == 2)
      cubicInterp_mp(preEG_approx, s, Z, dZ_ds, finalTime);
    else
      hermiteInterpCW_mp(preEG_approx, s, Z, dZ_ds, finalTime, num_samples);
  }
  else
  { // just copy the next to last sample point and call it good enough - Cauchy will take over from the last sample point and do its thing!
    point_cp_mp(preEG_approx, samples_mp[num_samples - 2].point);
  }

  // clear
  mpf_clear(tempMPF);
  clear_mp(factor);
  clear_mp(tempComp);
  clear_mp(c_times_t);
  for (i = num_samples - 1; i >= 0; i--)
  {
    clear_mp(&s[i]);
    clear_point_mp(&Z[i]);
    clear_point_mp(&dZ_ds[i]);
  }
  free(s);
  free(Z);
  free(dZ_ds);

  return;
}

int CauchyEG_actual_mp(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, point_mp preEG_approx, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual Cauchy integral endgame in multi prec  *
\***************************************************************/
{
  int i, retVal = 0, curr_cycle_size = 1;
  double closed_loop_tol, min_closed_loop_tol = MAX(1e-150, pow(10, -(prec_to_digits(T->Precision) - 5))), max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  mpf_t approx_err;
  comp_mp endTime;
  point_mp new_approx;
  point_data_mp *Samples = NULL;
  
  // initialize MP
  mpf_init(approx_err);
  init_mp(endTime); 
  init_point_mp(new_approx, 0);

  // setup the number of samples per loop
  *samples_per_loop = T->num_PSEG_sample_points + 1;
  
  // initially allocate for Samples
  Samples = (point_data_mp *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    init_point_data_mp(&Samples[i], 0);

  // setup last_approx to be preEG_approx
  point_cp_mp(last_approx, preEG_approx);

  // find the tolerance to close the first loop
  closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, Start, T->Precision, T, NULL, e_mp, NULL, ED, NULL, eval_func, NULL);

  // find the first set of samples
  retVal = find_Cauchy_samples_mp(closed_loop_tol, &Samples, cycle_num, *samples_per_loop, &curr_cycle_size, Start, T, OUT, e_mp, ED, eval_func);

  // see if we can continue
  if (!retVal)
  { // find new_approx
    find_Cauchy_approx_mp(new_approx, Samples, *cycle_num, *samples_per_loop);

    // main loop - keep tracking until we either have failure or we converge
    do
    { // compare the approximations
      findDiff_point(approx_err, NULL, new_approx, T->Precision, NULL, last_approx, T->Precision);
      T->error_at_latest_sample_point = mpf_get_d(approx_err);
fprintf(OUT, "cycle: %d error: %e\n", *cycle_num, mpf_get_d(approx_err));

      // check for convergence
      if (T->error_at_latest_sample_point < T->final_tolerance)
      { // we accept our approximation

        // setup Final
        point_cp_mp(Final->point, new_approx);
        set_mp(Final->time, finalTime);
        Final->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = mpf_get_d(Samples[0].time->r);

        retVal = 0;
        break;
      }
      else if (d_abs_mp(Samples[0].time) < T->minTrackT)
      { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

        // setup Final
        point_cp_mp(Final->point, new_approx);
        set_mp(Final->time, Samples[0].time);
        Final->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = mpf_get_d(Samples[0].time->r);

        retVal = retVal_EG_failed_to_converge;
        break;
      }
      else
      { // copy new_approx to last_approx
        point_cp_mp(last_approx, new_approx);

        // find end time for the next sample point
        mul_rdouble_mp(endTime, Samples[0].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_mp(&Samples[0], &Samples[0], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_mp_basic(&Samples[0], T, OUT, e_mp, ED, eval_func);

          // find the tolerance to close the next loop
          closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, &Samples[0], T->Precision, T, NULL, e_mp, NULL, ED, NULL, eval_func, NULL);

          // find the next set of samples
          retVal = find_Cauchy_samples_mp(closed_loop_tol, &Samples, cycle_num, *samples_per_loop, &curr_cycle_size, &Samples[0], T, OUT, e_mp, ED, eval_func);

          // check for errors
          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_mp(new_approx, Samples, *cycle_num, *samples_per_loop);
          }
        }
      }
    } while (!retVal);
  }

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    point_cp_mp(Final->point, Samples[0].point);
    set_mp(Final->time, Samples[0].time);
    Final->cycle_num = *cycle_num;

    T->t_val_at_latest_sample_point = mpf_get_d(Samples[0].time->r);
  }
  else
  { // allocate endLoopSamples
    *endLoopSamples = (point_data_mp *)bmalloc(*samples_per_loop * (*cycle_num) * sizeof(point_data_mp));

    // setup endLoopSamples
    for (i = *samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
    {
      init_point_data_mp(&(*endLoopSamples)[i], Samples[i].point->size);
      point_data_cp_mp(&(*endLoopSamples)[i], &Samples[i]);
    }
  }

  // clear memory
  mpf_clear(approx_err);
  clear_mp(endTime); 
  clear_point_mp(new_approx); 
  for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_mp(&Samples[i]);
  free(Samples);

  return retVal;
}

int circle_track_mp(point_data_mp Final[], int Final_digits[], point_data_mp *Start, int M, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks around a circle to find the M+1 sample points   *
\***************************************************************/
{
  int i, cont, retVal = 0, consecSuccess, numSteps = 0;
  mpf_t absValOfDistanceLeft, currentStepSize, tempMPF, radius, next_step_angle, next_stop_angle, curr_angle, two_pi;
  comp_mp nextTime;

  // error checking - want M >= 3 so that we do not have to predict directly across the origin
  if (M < 3)
  {
    printf("ERROR: The number of sample points (%d) for circle tracking must be >= 3!\n", M);
    bexit(ERROR_CONFIGURATION);
  }

  // initialize radius
  mpf_init(radius);

  // find the radius of the circle that we are going to track around
  mpf_abs_mp(radius, Start->time);

  // error checking
  if (mpf_cmp_ui(radius, 0) <= 0)
  {
    printf("ERROR: The radius of the circle (%e) needs to be positive!\n", mpf_get_d(radius));
    mpf_clear(radius);
    bexit(ERROR_CONFIGURATION);
  }

  // initialize other MP datatypes
  mpf_init(absValOfDistanceLeft); mpf_init(currentStepSize); mpf_init(tempMPF);
  mpf_init(next_step_angle); mpf_init(next_stop_angle); mpf_init(curr_angle); mpf_init(two_pi);
  init_mp(nextTime);

  // setup two_pi
  mpfr_const_pi(two_pi, __gmp_default_rounding_mode);
  mpf_mul_ui(two_pi, two_pi, 2);

  // setup curr_angle = -2 * PI
  mpf_neg(curr_angle, two_pi);

  // if this is the first step of the path, initialize the current step size as the maximum step size - otherwise, it is already set!
  if (T->first_step_of_path)
  {
    T->currentStepSize = T->maxStepSize;
    T->first_step_of_path = 0;
  }

  // sample point [0] is Start
  point_data_cp_mp(&Final[0], Start);
  refine_mp_basic(&Final[0], T, OUT, e, ED, eval_func);
  Final_digits[0] = residual_to_digits_mp(T->latest_newton_residual_mp, T->Precision);

  // start with the current step size
  mpf_set_d(currentStepSize, T->currentStepSize);

  // loop over each sample point
  consecSuccess = 0;
  for (i = 1; i <= M && retVal == 0; i++)
  { // setup the starting point
    point_cp_mp(Final[i].point, Final[i-1].point);
    set_mp(Final[i].time, Final[i-1].time);

    // reset the number of steps
    numSteps = 0;

    // calculate the next angle to stop at = -2 PI + 2 PI i / M = 2 PI (i / M - 1)
    mpf_set_ui(next_stop_angle, i);
    mpf_div_ui(next_stop_angle, next_stop_angle, M);
    mpf_sub_ui(next_stop_angle, next_stop_angle, 1);
    mpf_mul(next_stop_angle, next_stop_angle, two_pi);

    // track to the next sample point
    cont = 1;
    while (cont)
    { // print out the current time
      if (T->outputLevel > 0)
      {
        fprintf(OUT, "Time = "); print_mp(OUT, 0, Final[i].time); fprintf(OUT, "\n");
        if (T->screenOut)
        {
          fprintf(stdout, "Time = "); print_mp(stdout, 0, Final[i].time); fprintf(stdout, "\n");
        }
      }

      // find the distance to the end = radius * (next_stop_angle - curr_angle)
      mpf_sub(absValOfDistanceLeft, next_stop_angle, curr_angle);
      mpf_mul(absValOfDistanceLeft, absValOfDistanceLeft, radius);

      // check to see how close we are to the target
      if (mpfr_equal_p(next_stop_angle, curr_angle))
      { // we are exactly at the end!
        cont = 0;
      }
      else if (mpfr_less_p(absValOfDistanceLeft, currentStepSize))
      { // we can take a smaller step than our current step size to hit the target time
        mpf_set(next_step_angle, next_stop_angle);
      }
      else
      { // the distance left is larger than the current step size, so we head in the correct direction will a step size of the current step size
        mpf_div(tempMPF, currentStepSize, radius); // calculate the change in angle that we can do
        mpf_add(next_step_angle, curr_angle, tempMPF); // calculate the next step angle
      }

      // check to see if we should continue
      if (cont)
      { // setup nextTime
        mpfr_cos(nextTime->r, next_step_angle, __gmp_default_rounding_mode);
        mpfr_sin(nextTime->i, next_step_angle, __gmp_default_rounding_mode);
        mul_rmpf_mp(nextTime, nextTime, radius);

        // take a step
        retVal = step_mp(&Final[i], nextTime, T, OUT, e, ED, eval_func);
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
          // update curr_angle
          mpf_set(curr_angle, next_step_angle);

          // check to see if we increase the step size
          if (consecSuccess >= T->cSecInc)
          { // find the new step size if we increase it
            mpf_set_d(tempMPF, T->step_success_factor);
            mpf_mul(tempMPF, currentStepSize, tempMPF);
            // the new step size is the minimum of the max step size and tempD
            if (T->maxStepSize < mpf_get_d(tempMPF))
              mpf_set_d(currentStepSize, T->maxStepSize);
            else
              mpf_set(currentStepSize, tempMPF);

            consecSuccess = 0;
          }
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
          }
          else if (retVal == retVal_going_to_infinity)
          { // path is going to infinity
            retVal = retVal_going_to_infinity;
            cont = 0;
          }
          else if (retVal == retVal_higher_prec_needed)
          { // higher precision is needed
            retVal = retVal_higher_prec_needed;
            cont = 0;
          }
        }

        // check to see if we have taken too many steps
        if (numSteps > T->maxNumSteps)
        {
          retVal = retVal_too_many_steps;
          cont = 0;
        }
      }
    }

    // display how many steps it took
    if (T->outputLevel >= 0)
    {
      fprintf(OUT, "numSteps = %d\n", numSteps);
      if (T->screenOut)
        printf("numSteps = %d\n", numSteps);
    }

    // refine the sample, if it was good
    if (retVal == 0)
    { 
      refine_mp_basic(&Final[i], T, OUT, e, ED, eval_func);
      Final_digits[i] = residual_to_digits_mp(T->latest_newton_residual_mp, T->Precision);
    }
  }

  // clear the MP datatypes
  mpf_clear(absValOfDistanceLeft); mpf_clear(currentStepSize); mpf_clear(tempMPF); mpf_clear(radius);
  mpf_clear(next_step_angle); mpf_clear(next_stop_angle); mpf_clear(curr_angle); mpf_clear(two_pi);
  clear_mp(nextTime);

  // free the cache associated with calculating PI
  mpfr_free_cache();

  return retVal;
}

int find_Cauchy_samples_mp(double closed_loop_max_error, point_data_mp **Samples, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_mp *Start, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Populates Samples and finds the cycle number to close  *
* the loop                                                      *
\***************************************************************/
{
  int j, tempInt, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  comp_mp time;

  // save the correct time
  init_mp(time);
  set_mp(time, Start->time);
    
  // setup Samples[0] - circle_track_mp will refine this sample
  point_data_cp_mp(&(*Samples)[0], Start);

  // loop until get back to the original point
  *cycle_num = 0; 
  while (1)
  { // go around the origin once and refine the samples as they are found
    retVal = circle_track_mp(&(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], &(*Samples)[*cycle_num * num_samples], num_samples, T, OUT, e_mp, ED, eval_func);

    // increment the cycle_num
    (*cycle_num)++;

    if (retVal)
    { // return error
      fprintf(OUT, "Note: The tracking failed while going around the origin!\n");
      break;
    }
    else
    { // check to see if we have closed the loop
      if (check_closed_loop_mp(closed_loop_max_error, &(*Samples)[0], &s_digits[0], &(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time, e_mp, T, OUT, ED, eval_func))
      { // error is small enough - exit loop with success
        retVal = 0;
        break;
      }
      else if (*cycle_num > fail_safe_cycle_max)
      { // too many iterations
        fprintf(OUT, "ERROR: Cycle number too high to detect!\n");
        retVal = retVal_cycle_num_too_high;
        break;
      }
      else if (*cycle_num >= *curr_cycle_size)
      { // need to increase the size of memory
        s_digits = (int *)brealloc(s_digits, (2 * (*curr_cycle_size) * num_samples + 1) * sizeof(int));

        point_data_mp *temp_ptr = *Samples; 
        *Samples = NULL;
        tempInt = (*curr_cycle_size) * num_samples;
        *Samples = (point_data_mp *)bmalloc((2 * (*curr_cycle_size) * num_samples + 1) * sizeof(point_data_mp));
        for (j = 2 * (*curr_cycle_size) * num_samples; j >= 0; j--)
        {
          init_point_data_mp(&(*Samples)[j], 0);

          if (j <= tempInt)
          {
            point_data_cp_mp(&(*Samples)[j], &temp_ptr[j]);
            clear_point_data_mp(&temp_ptr[j]);
          }
        }
        free(temp_ptr);
        temp_ptr = NULL;

        // update curr_cycle_size
        *curr_cycle_size = 2 * (*curr_cycle_size);
      }
    }
  }

  // release memory
  free(s_digits);
  clear_mp(time);

  return retVal;
}

void find_Cauchy_approx_mp(point_mp approx, point_data_mp Samples[], int cycle_num, int num_samples)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Finds the Cauchy integral approximation using Trapezoid*
* rule based on cycle_num & num_samples                         *
\***************************************************************/
{
  int i, j, size, sample_num = cycle_num * num_samples; // the first and last sample are equal so we skip the last one

  // set size
  size = Samples[0].point->size;
  change_size_point_mp(approx, size);
  approx->size = size;

  for (i = 0; i < size; i++)
  { // initialize to 0
    set_zero_mp(&approx->coord[i]);

    // sum over all of the samples
    for (j = 0; j < sample_num; j++)
    {
      add_mp(&approx->coord[i], &approx->coord[i], &Samples[j].point->coord[i]);
    }

    // divide by number of samples added together
    mpf_div_ui(approx->coord[i].r, approx->coord[i].r, sample_num);
    mpf_div_ui(approx->coord[i].i, approx->coord[i].i, sample_num);
  }

  return;
}

int check_closed_loop_mp(double closed_loop_max_error, point_data_mp *Pt0, int *digitsCorrect0, point_data_mp *Pt1, int *digitsCorrect1, comp_mp time, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - Pt0 == Pt1, 0 - otherwise                  *
* NOTES: check to see if we have closed the loop for Cauchy EG  *
\***************************************************************/
{
  int retVal = 0, min_digits = 0;

  // setup min_digits
  min_digits = ceil(-log10(closed_loop_max_error) + 0.5); // minimum number of digits that we need to have correct for each point

  // check to see if we have the correct number of digits
  if (*digitsCorrect0 < min_digits)
  { // need to refine Pt0
    refine_digits_mp(T->outputLevel, min_digits, T->latest_newton_residual_mp, *digitsCorrect0, Pt0, Pt0, T->Precision, time, OUT, e_mp, ED, eval_func);

    // update digitsCorrect0
    *digitsCorrect0 = min_digits;
  }

  if (*digitsCorrect1 < min_digits)
  { // need to refine Pt1
    refine_digits_mp(T->outputLevel, min_digits, T->latest_newton_residual_mp, *digitsCorrect1, Pt1, Pt1, T->Precision, time, OUT, e_mp, ED, eval_func);

    // update digitsCorrect1
    *digitsCorrect1 = min_digits;
  }

  // see if we have closed the loop
  retVal = isSamePoint(NULL, Pt0->point, T->Precision, NULL, Pt1->point, T->Precision, closed_loop_max_error);

  return retVal;
}

////// Adaptive Multi Precision ///////

int CauchyEG_amp(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in adaptive precision     *
\***************************************************************/
{
  int i, retVal, cycle, samples_per_loop;
  comp_d finalTime_d;
  comp_mp finalTime_mp;
  point_data_d *endSamples_d = NULL;
  point_data_mp *endSamples_mp = NULL;

  init_mp(finalTime_mp);

  // setup finalTime
  set_double_d(finalTime_d, T->targetT, 0);
  d_to_mp(finalTime_mp, finalTime_d);

  retVal = CauchyEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &endSamples_d, &endSamples_mp, &cycle, &samples_per_loop, finalTime_d, finalTime_mp, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  // clear endSamples since is it not being used
  if ((retVal == 0 || retVal == retVal_EG_failed_to_converge) && *prec < 64)
  { // clear _d
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_d(&endSamples_d[i]);
    free(endSamples_d);
  }
  else if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // clear _mp
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_mp(&endSamples_mp[i]);
    free(endSamples_mp);
  }

  clear_mp(finalTime_mp);

  return retVal;
}

int CauchyEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endSamples_d, point_data_mp **endSamples_mp, int *cycle, int *samples_per_loop, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *)) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endSamples, cycle, samples_per_loop            *
* NOTES: does Cauchy integral endgame in adaptive precision     *
\***************************************************************/
{
  int i, retVal;
  comp_d endTime_d;
  comp_mp endTime_mp;

  // start off with atleast 64 bit precision
  T->Precision = MAX(64, prec_in);
  initMP(T->Precision);
  change_prec(ED_mp, T->Precision);

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

  // try to track to the endgame boundary - this will handle everything
  retVal = AMPtrack(Final_d, Final_mp, prec, time_first_increase, Start_d, Start_mp, prec_in, endTime_d, endTime_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // so, either made it to the endgame boundary or path failure - first print to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  if (*prec == 52)
  { // print to midOUT using Final_d
    for (i = 0; i < Final_d->point->size; i++)
      fprintf(midOUT, "%.15e %.15e\n", Final_d->point->coord[i].r, Final_d->point->coord[i].i);
    fprintf(midOUT, "\n");
  }
  else
  { // print to midOUT using Final_mp
    for (i = 0; i < Final_mp->point->size; i++)
    {
      print_mp(midOUT, 0, &Final_mp->point->coord[i]);
      fprintf(midOUT, "\n");
    }
    fprintf(midOUT, "\n");
  }
  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check for success
  if (retVal)
  { // failure
    fprintf(OUT, "NOTE: Cauchy Integral Endgame never started!\n");

    if (*prec < 64)
      T->t_val_at_latest_sample_point = Final_d->time->r;
    else
      T->t_val_at_latest_sample_point = mpf_get_d(Final_mp->time->r);
  }
  else
  { // success - so run the acutal endgame

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // run the endgame
    retVal = CauchyEG_main_amp2(prec, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, endSamples_d, endSamples_mp, cycle, samples_per_loop, time_first_increase, finalTime_d, finalTime_mp, Final_d, Final_mp, *prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
  }

  clear_mp(endTime_mp);

  return retVal;
}

int CauchyEG_main_amp(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, double *time_first_increase, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in double precision            *
\***************************************************************/
{
  int retVal = 0, num_agreement = NUM_AGREEMENT;
  double timeCutoff = TIME_CUTOFF;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  // initialize
  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // perform the 'pre-Cauchy EG' - use cycle number approximations to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_amp(last_approx_d, last_approx_mp, last_approx_prec, cycle_num, Final_d, Final_mp, prec_out, finalTime_d, finalTime_mp, Start_d, Start_mp, prec_in, num_agreement, timeCutoff, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  if (retVal == 0)
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_amp(prec_out, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, endLoopSamples_d, endLoopSamples_mp, cycle_num, samples_per_loop, last_approx_d, last_approx_mp, *last_approx_prec, finalTime_d, finalTime_mp, Final_d, Final_mp, *prec_out, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }

  // clear
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  return retVal;
}

int CauchyEG_preEG_amp(point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, int *cycle_num, point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, int num_agreement, double timeCutoff, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  *
* can efficiently starting using the actual Cauchy EG.  Used    *
* timeCutoff as a fail-safe stopping criterion.                 *
\***************************************************************/
{ 
  int i, curr_sample, last_sample = 0, retVal = 0, agreement = 0, *cycle = (int *)bmalloc(num_agreement * sizeof(int)), *samples_prec = (int *)bmalloc(3 * sizeof(int));
  point_data_d *samples_d = (point_data_d *)bmalloc(3 * sizeof(point_data_d)); // only 3 sample points are needed for a cycle number approximation
  point_data_mp *samples_mp = (point_data_mp *)bmalloc(3 * sizeof(point_data_mp));
  comp_d endTime_d;
  comp_mp endTime_mp;

  init_mp(endTime_mp);
  for (i = 0; i < 3; i++)
  {
    init_point_data_d(&samples_d[i], 0);
    init_point_data_mp(&samples_mp[i], 0);
  }

  // initialize cycle
  *cycle_num = Final_d->cycle_num = Final_mp->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    cycle[i] = -1;

  // setup the first sample and refine it
  samples_prec[0] = prec_in;
  if (prec_in < 64)
  { // use _d
    point_data_cp_d(&samples_d[0], Start_d);
    refine_d_basic(&samples_d[0], T, OUT, e_d, ED_d, eval_func_d);
  }
  else
  { // use _mp
    initMP(prec_in);
    change_prec(ED_mp, prec_in);
    setprec_eval_struct_mp(*e_mp, prec_in);
    setprec_point_mp(samples_mp[0].point, samples_prec[0]);
    setprec_mp(samples_mp[0].time, samples_prec[0]);

    point_data_cp_mp(&samples_mp[0], Start_mp);
    refine_mp_basic(&samples_mp[0], T, OUT, e_mp, ED_mp, eval_func_mp);
  }

  // find the first 3 samples
  for (i = 1; i < 3 && retVal == 0; i++)
  { // find the end time for the next sample
    if (samples_prec[i-1] < 64)
    { // setup endTime_d
      mul_rdouble_d(endTime_d, samples_d[i-1].time, T->power_series_sample_factor);
    }
    else
    { // setup endTime_mp
      setprec_mp(endTime_mp, samples_prec[i-1]);
      mul_rdouble_mp(endTime_mp, samples_mp[i-1].time, T->power_series_sample_factor);
    }

    // track to the next sample point
    retVal = AMPtrack(&samples_d[i], &samples_mp[i], &samples_prec[i], time_first_increase, &samples_d[i-1], &samples_mp[i-1], samples_prec[i-1], endTime_d, endTime_mp, samples_prec[i-1], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for errors
    if (retVal)
    { // store where this error occured
      last_sample = i;
    }
  }

  // see if we should continue
  if (!retVal)
  { // store where the last sample will be
    last_sample = 2;

    // set all sample points to the current precision and find the first approximation for the cycle number
    if (samples_prec[2] < 64)
    { // convert all to double and refine
      for (i = 0; i < 2; i++)
      {
        if (samples_prec[i] != samples_prec[2])
        { // create a double precision copy
          convert_point_data_mp_to_d(&samples_d[i], &samples_mp[i]);
        }
        // refine this sample
        refine_d_basic(&samples_d[i], T, OUT, e_d, ED_d, eval_func_d);
      }

      // find approximation for the cycle number using double precision
      cycle[0] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);
    }
    else
    { // convert all to multi precision and refine
      initMP(samples_prec[2]);
      change_prec(ED_mp, samples_prec[2]);
      setprec_eval_struct_mp(*e_mp, samples_prec[2]);
      for (i = 0; i < 2; i++)
      {
        if (samples_prec[i] < 64)
        { // create a multi precison copy
          setprec_point_mp(samples_mp[i].point, samples_prec[2]);
          setprec_mp(samples_mp[i].time, samples_prec[2]);

          convert_point_data_d_to_mp(&samples_mp[i], &samples_d[i]);
        }
        else if (samples_prec[i] < samples_prec[2])
        { // increase the precision
          change_prec_point_mp(samples_mp[i].point, samples_prec[2]);
          change_prec_mp(samples_mp[i].time, samples_prec[2]);
        }
        // refine this sample
        refine_mp_basic(&samples_mp[i], T, OUT, e_mp, ED_mp, eval_func_mp);
      }

      // find approximation for the cycle number using mulit precision
      cycle[0] = find_cycle_num_mp(samples_mp, 0, T->power_series_sample_factor);
    }

    // loop to find the other cycle numbers
    curr_sample = 1;
    while (retVal == 0 && curr_sample < num_agreement && (samples_prec[2] < 64 ? d_abs_d(samples_d[2].time) > timeCutoff : d_abs_mp(samples_mp[2].time) > timeCutoff))
    { // update samples
      for (i = 1; i < 3; i++)
      { // update precision
        samples_prec[i-1] = samples_prec[i];
        // copy samples
        if (samples_prec[i] < 64)
        { // copy _d
          point_data_cp_d(&samples_d[i-1], &samples_d[i]);
        }
        else
        { // copy _mp
          setprec_point_mp(samples_mp[i-1].point, samples_prec[i]);
          setprec_mp(samples_mp[i-1].time, samples_prec[i]);

          point_data_cp_mp(&samples_mp[i-1], &samples_mp[i]);
        }
      }

      // find the end time for the next sample
      if (samples_prec[1] < 64)
      { // setup endTime_d
        mul_rdouble_d(endTime_d, samples_d[1].time, T->power_series_sample_factor);
      }
      else
      { // setup endTime_mp
        setprec_mp(endTime_mp, samples_prec[1]);
        mul_rdouble_mp(endTime_mp, samples_mp[1].time, T->power_series_sample_factor);
      }

      // track to the next sample point
      retVal = AMPtrack(&samples_d[2], &samples_mp[2], &samples_prec[2], time_first_increase, &samples_d[1], &samples_mp[1], samples_prec[1], endTime_d, endTime_mp, samples_prec[1], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

      // check for errors
      if (!retVal)
      { // set all sample points to the current precision and find the next approximation for the cycle number
        if (samples_prec[2] < 64)
        { // convert all to double and refine
          for (i = 0; i < 2; i++)
          {
            if (samples_prec[i] != samples_prec[2])
            { // create a double precision copy
              convert_point_data_mp_to_d(&samples_d[i], &samples_mp[i]);
            }
            // refine this sample
            refine_d_basic(&samples_d[i], T, OUT, e_d, ED_d, eval_func_d);
          }

          // find approximation for the cycle number using double precision
          cycle[curr_sample] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);
        }
        else
        { // convert all to multi precision and refine
          initMP(samples_prec[2]);
          change_prec(ED_mp, samples_prec[2]);
          setprec_eval_struct_mp(*e_mp, samples_prec[2]);
          for (i = 0; i < 2; i++)
          {
            if (samples_prec[i] < 64)
            { // create a multi precison copy
              setprec_point_mp(samples_mp[i].point, samples_prec[2]);
              setprec_mp(samples_mp[i].time, samples_prec[2]);

              convert_point_data_d_to_mp(&samples_mp[i], &samples_d[i]);
            }
            else if (samples_prec[i] < samples_prec[2])
            { // increase the precision
              change_prec_point_mp(samples_mp[i].point, samples_prec[2]);
              change_prec_mp(samples_mp[i].time, samples_prec[2]);
            }
            // refine this sample
            refine_mp_basic(&samples_mp[i], T, OUT, e_mp, ED_mp, eval_func_mp);
          }

          // find approximation for the cycle number using mulit precision
          cycle[curr_sample] = find_cycle_num_mp(samples_mp, 0, T->power_series_sample_factor);
        }
      }

      // increment curr_sample
      curr_sample++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = 1;
      for (i = 1; i < num_agreement && agreement; i++)
        if (cycle[0] != cycle[i])
          agreement = 0;

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && (samples_prec[2] < 64 ? d_abs_d(samples_d[2].time) > timeCutoff : d_abs_mp(samples_mp[2].time) > timeCutoff))
      { // update cycle
        for (i = 1; i < num_agreement; i++)
          cycle[i-1] = cycle[i];

        // update samples
        for (i = 1; i < 3; i++)
        { // update precision
          samples_prec[i-1] = samples_prec[i];
          // copy samples
          if (samples_prec[i] < 64)
          { // copy _d
            point_data_cp_d(&samples_d[i-1], &samples_d[i]);
          }
          else
          { // copy _mp
            setprec_point_mp(samples_mp[i-1].point, samples_prec[i]);
            setprec_mp(samples_mp[i-1].time, samples_prec[i]);

            point_data_cp_mp(&samples_mp[i-1], &samples_mp[i]);
          }
        }

        // find the end time for the next sample
        if (samples_prec[1] < 64)
        { // setup endTime_d
          mul_rdouble_d(endTime_d, samples_d[1].time, T->power_series_sample_factor);
        }
        else
        { // setup endTime_mp
          setprec_mp(endTime_mp, samples_prec[1]);
          mul_rdouble_mp(endTime_mp, samples_mp[1].time, T->power_series_sample_factor);
        }

        // track to the next sample point
        retVal = AMPtrack(&samples_d[2], &samples_mp[2], &samples_prec[2], time_first_increase, &samples_d[1], &samples_mp[1], samples_prec[1], endTime_d, endTime_mp, samples_prec[1], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // check for errors
        if (!retVal)
        { // set all sample points to the current precision and find the next approximation for the cycle number
          if (samples_prec[2] < 64)
          { // convert all to double and refine
            for (i = 0; i < 2; i++)
            {
              if (samples_prec[i] != samples_prec[2])
              { // create a double precision copy
                convert_point_data_mp_to_d(&samples_d[i], &samples_mp[i]);
              }
              // refine this sample
              refine_d_basic(&samples_d[i], T, OUT, e_d, ED_d, eval_func_d);
            }

            // find approximation for the cycle number using double precision
            cycle[num_agreement - 1] = find_cycle_num_d(samples_d, 0, T->power_series_sample_factor);
          }
          else
          { // convert all to multi precision and refine
            initMP(samples_prec[2]);
            change_prec(ED_mp, samples_prec[2]);
            setprec_eval_struct_mp(*e_mp, samples_prec[2]);
            for (i = 0; i < 2; i++)
            {
              if (samples_prec[i] < 64)
              { // create a multi precison copy
                setprec_point_mp(samples_mp[i].point, samples_prec[2]);
                setprec_mp(samples_mp[i].time, samples_prec[2]);

                convert_point_data_d_to_mp(&samples_mp[i], &samples_d[i]);
              }
              else if (samples_prec[i] < samples_prec[2])
              { // increase the precision
                change_prec_point_mp(samples_mp[i].point, samples_prec[2]);
                change_prec_mp(samples_mp[i].time, samples_prec[2]);
              }
              // refine this sample
              refine_mp_basic(&samples_mp[i], T, OUT, e_mp, ED_mp, eval_func_mp);
            }

            // find approximation for the cycle number using mulit precision
            cycle[num_agreement - 1] = find_cycle_num_mp(samples_mp, 0, T->power_series_sample_factor);
          }

          // check for agreement
          agreement = 1;
          for (i = 1; i < num_agreement && agreement; i++)
            if (cycle[0] != cycle[i])
              agreement = 0;
        }
      }
    }
  }


  // check for errors
  if (retVal == 0)
  { // find the power series approximation described by the 3 sample points
    find_Cauchy_preEG_approx_amp(preEG_approx_d, preEG_approx_mp, preEG_prec, finalTime_d, finalTime_mp, samples_d, samples_mp, samples_prec, 3, cycle[num_agreement - 1], OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }

  // update Final
  *prec_out = samples_prec[last_sample];
  if (*prec_out < 64)
  { // copy _d
    point_data_cp_d(Final_d, &samples_d[last_sample]);
    *cycle_num = Final_d->cycle_num = cycle[num_agreement - 1];
    T->t_val_at_latest_sample_point = samples_d[last_sample].time->r;
  }
  else
  { // copy _mp
    setprec_point_mp(Final_mp->point, *prec_out);
    setprec_mp(Final_mp->time, *prec_out);

    point_data_cp_mp(Final_mp, &samples_mp[last_sample]);
    *cycle_num = Final_mp->cycle_num = cycle[num_agreement - 1];
    T->t_val_at_latest_sample_point = mpf_get_d(samples_mp[last_sample].time->r);
  }

  // clear memory
  clear_mp(endTime_mp);
  for (i = 2; i >= 0; i--)
  {
    clear_point_data_d(&samples_d[i]);
    clear_point_data_mp(&samples_mp[i]);
  }
  free(samples_mp);
  free(samples_d);
  free(samples_prec);
  free(cycle);

  return retVal;
}

void find_Cauchy_preEG_approx_amp(point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *samples_d, point_data_mp *samples_mp, int *samples_prec, int num_samples, int cycle_num, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the approximation for Cauchy preEG               *
\***************************************************************/
{ // error checking
  if (cycle_num <= 0 || num_samples <= 1)
  {
    printf("ERROR: Incorrect settings when trying to do a power series predicition!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // cycle_num >= 1 && num_samples >= 2
  int i;

  // find the maximum precision
  *preEG_prec = samples_prec[0];
  for (i = 1; i < num_samples; i++)
    if (*preEG_prec < samples_prec[i])
      *preEG_prec = samples_prec[i];

  if (*preEG_prec < 64)
  { // try to find the approximation using double precision
    find_Cauchy_preEG_approx_d(preEG_approx_d, finalTime_d, samples_d, num_samples, cycle_num, OUT, e_d, ED_d, eval_func_d);
  }
  else
  { // set the precision
    initMP(*preEG_prec);
    change_prec(ED_mp, *preEG_prec);
    setprec_eval_struct_mp(*e_mp, *preEG_prec);
    setprec_point_mp(preEG_approx_mp, *preEG_prec);

    setprec_mp(finalTime_mp, *preEG_prec);
    d_to_mp(finalTime_mp, finalTime_d);

    // convert all samples to the precision 
    for (i = 0; i < num_samples; i++)
      if (samples_prec[i] < 64)
      { // convert to MP
        setprec_point_data_mp(&samples_mp[i], *preEG_prec);
        convert_point_data_d_to_mp(&samples_mp[i], &samples_d[i]);
      }
      else if (samples_prec[i] < *preEG_prec)
      { // increase precision
        change_prec_point_mp(samples_mp[i].point, *preEG_prec);
        change_prec_mp2(samples_mp[i].time, *preEG_prec);
      }

    // try to find the approximation using the precision
    find_Cauchy_preEG_approx_mp(preEG_approx_mp, finalTime_mp, samples_mp, num_samples, cycle_num, OUT, e_mp, ED_mp, eval_func_mp);
  }

  return;
}

void newError(mpf_t error, point_d new_approx_d, point_mp new_approx_mp, int *new_approx_prec, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, double tol)
{
int i, prec = MAX(*new_approx_prec, *last_approx_prec);
if (prec < 64)
{ // _d
int size = new_approx_d->size;
double mult = 1.0 / 999.0, err[2] = {0,0};
point_d diff, diff2, test_approx;
init_point_d(diff, size); init_point_d(diff2,size);init_point_d(test_approx, size);
diff->size = diff2->size = test_approx->size = size;

for (i = 0; i < size; i++)
{
 set_zero_d(&diff->coord[i]);
 sub_d(&diff->coord[i], &new_approx_d->coord[i], &last_approx_d->coord[i]);
 mul_rdouble_d(&test_approx->coord[i], &diff->coord[i], mult);
 sub_d(&test_approx->coord[i], &new_approx_d->coord[i], &test_approx->coord[i]);
 sub_d(&diff2->coord[i], &test_approx->coord[i], &new_approx_d->coord[i]);
}
twoNormVec_d(diff, &err[0]); twoNormVec_d(diff2, &err[1]);
printf("reg error: %e new error: %e\n", err[0], err[1]);

mpf_set_prec(error, 64); mpf_set_d(error, err[1]);

//if (err[1] < tol) point_cp_d(new_approx_d, test_approx);
}
else
{
point_mp approx1, approx2;
init_point_mp2(approx1,0,prec); init_point_mp2(approx2,0,prec);

if (*new_approx_prec < 64) {point_d_to_mp(approx1,new_approx_d);}
else {point_cp_mp(approx1,new_approx_mp);}

if (*last_approx_prec < 64) {point_d_to_mp(approx2,last_approx_d);}
else {point_cp_mp(approx2,last_approx_mp);}

int size = approx1->size;
mpf_t mult;
mpf_init2(mult, prec); mpf_set_ui(mult, 999);
mpf_ui_div(mult, 1, mult);
comp_mp err[2];
point_mp diff, diff2, test_approx;
init_point_mp2(diff, size, prec); init_point_mp2(diff2,size, prec);init_point_mp2(test_approx, size, prec);
diff->size = diff2->size = test_approx->size = size;
init_mp2(err[0],prec); init_mp2(err[1],prec);

for (i = 0; i < size; i++)
{
 set_zero_mp(&diff->coord[i]);
 sub_mp(&diff->coord[i], &approx1->coord[i], &approx2->coord[i]);
 mul_rmpf_mp(&test_approx->coord[i], &diff->coord[i], mult);
 sub_mp(&test_approx->coord[i], &approx1->coord[i], &test_approx->coord[i]);
 sub_mp(&diff2->coord[i], &test_approx->coord[i], &approx1->coord[i]);
}
twoNormVec_mp(diff, err[0]); twoNormVec_mp(diff2, err[1]);
printf("reg error: %e new error: %e\n", mpf_get_d(err[0]->r), mpf_get_d(err[1]->r));

mpf_set_prec(error, prec); mpf_set(error, err[1]->r);

//if (mpf_get_d(err[1]->r) < tol) point_cp_mp(new_approx_mp, test_approx);
}
return;
}

void testCauchy(int new_point, point_d average_d, point_mp average_mp, point_data_d *Samples_d, point_data_mp *Samples_mp, int prec, int cycle_num, int samples_per_loop, double sample_factor)
{ // find min & max
int i, j, size;
//int k, l;
double time, norm, max, min, max_diff = 0, total;

static int count = 0;
static point_data_mp *s_mp = NULL;

if (new_point)
{ for (j = 0; j < count; j++) clear_point_data_mp(&s_mp[j]); 
  free(s_mp); s_mp = NULL;
  count = 0;
}

if (prec < 64)
{ comp_d diff;
  size = Samples_d[0].point->size; 

//printf("z1="); printVec_Matlab_d(stdout,0,Samples_d[0].point);
//printf("z2="); printVec_Matlab_d(stdout,0,Samples_d[2].point);
//printf("z3="); if (cycle_num > 1) printVec_Matlab_d(stdout,0,Samples_d[4].point); else printVec_Matlab_d(stdout,0,Samples_d[0].point);

point_d testPt, test_d[5]; 
//double err_d[4];
init_point_d(testPt, size);
for (i = 0; i < 5; i++) {init_point_d(test_d[i], size); test_d[i]->size = size; for (j = 0; j < size; j++) set_zero_d(&test_d[i]->coord[j]); }
/*
for (i = 0; i < 5; i++)
{
  if (i == 0) l = 1; else if (i == 1) l = 2; else if (i == 2) l = 4; else if (i == 3) l = 8; else l = 16;

  for (j = 0; j < size; j++)
  {
   for (k = 0; k < cycle_num * samples_per_loop; k++)
    if (!(k % l))
    {
      add_d(&test_d[i]->coord[j], &test_d[i]->coord[j], &Samples_d[k].point->coord[j]);
    }
    mul_rdouble_d(&test_d[i]->coord[j], &test_d[i]->coord[j], 1.0 / (samples_per_loop * cycle_num / l));
  }
}
vec_sub_d(testPt, test_d[0], test_d[1]); twoNormVec_d(testPt, &err_d[0]);
vec_sub_d(testPt, test_d[0], test_d[2]); twoNormVec_d(testPt, &err_d[1]);
vec_sub_d(testPt, test_d[0], test_d[3]); twoNormVec_d(testPt, &err_d[2]);
vec_sub_d(testPt, test_d[0], test_d[4]); twoNormVec_d(testPt, &err_d[3]);
printf("err: %e %e %e %e div: %e %e %e\n", err_d[0], err_d[1], err_d[2], err_d[3], err_d[0]/err_d[1], err_d[1]/err_d[2], err_d[2]/err_d[3]);
*/
clear_point_d(testPt);
for (i = 0; i < 5; i++) clear_point_d(test_d[i]);


s_mp = (point_data_mp *)brealloc(s_mp, (count + 1) * sizeof(point_data_mp));
init_point_data_mp(&s_mp[count], 0);
convert_point_data_d_to_mp(&s_mp[count], &Samples_d[0]);
count++;

  time = Samples_d[0].time->r;
  for (j = 0; j < size; j++)
  { max = 0; min = 10000; total = 0;
    for (i = 0; i < cycle_num * samples_per_loop; i++)
    {
      sub_d(diff, &Samples_d[i].point->coord[j], &average_d->coord[j]);
      total += norm = d_abs_d(diff);
      if (norm > max) max = norm;
      if (norm < min) min = norm;
      norm = max - min; if (norm > max_diff) max_diff = norm;
    }
    printf("coord: %d max: %e min: %e ratio: %e average: %e\n", j, max, min, min/max, total / (cycle_num * samples_per_loop));
  }
  printf("time: %e cycle: %d max_diff: %e\n", time, cycle_num, max_diff);
}
else
{ comp_mp diff; init_mp2(diff, prec);
  size = Samples_mp[0].point->size;

//printf("z1="); printVec_Matlab_mp(stdout,0,Samples_mp[0].point);
//printf("z2="); printVec_Matlab_mp(stdout,0,Samples_mp[2].point);
//printf("z3="); if (cycle_num > 1) printVec_Matlab_mp(stdout,0,Samples_mp[4].point); else printVec_Matlab_mp(stdout,0,Samples_mp[0].point);

point_mp testPt, *test_mp = (point_mp *)bmalloc(5 * sizeof(point_mp)); comp_mp err_mp[4];
init_point_mp2(testPt, size, prec);
for (i = 0; i < 5; i++) {init_mp2(err_mp[i], prec); init_point_mp2(test_mp[i], size, prec); test_mp[i]->size = size; for (j = 0; j < size; j++) {set_zero_mp(&test_mp[i]->coord[j]); }}
/*
for (i = 0; i < 5; i++)
{
  if (i == 0) l = 1; else if (i == 1) l = 2; else if (i == 2) l = 4; else if (i == 3) l = 8; else l = 16;

  for (j = 0; j < size; j++)
  {
   for (k = 0; k < cycle_num * samples_per_loop; k++)
    if (!(k % l))
    {
      add_mp(&test_mp[i]->coord[j], &test_mp[i]->coord[j], &Samples_mp[k].point->coord[j]);
    }
    mul_rdouble_mp(&test_mp[i]->coord[j], &test_mp[i]->coord[j], 1.0 / (samples_per_loop * cycle_num / l));
  }
}
vec_sub_mp(testPt, test_mp[0], test_mp[1]); twoNormVec_mp(testPt, err_mp[0]);
vec_sub_mp(testPt, test_mp[0], test_mp[2]); twoNormVec_mp(testPt, err_mp[1]);
vec_sub_mp(testPt, test_mp[0], test_mp[3]); twoNormVec_mp(testPt, err_mp[2]);
vec_sub_mp(testPt, test_mp[0], test_mp[4]); twoNormVec_mp(testPt, err_mp[3]);
printf("err: %e %e %e %e", mpf_get_d(err_mp[0]->r), mpf_get_d(err_mp[1]->r), mpf_get_d(err_mp[2]->r), mpf_get_d(err_mp[3]->r));
div_mp(err_mp[0], err_mp[0], err_mp[1]);
div_mp(err_mp[1], err_mp[1], err_mp[2]);
div_mp(err_mp[2], err_mp[2], err_mp[3]);
printf(" div: %e %e %e\n", mpf_get_d(err_mp[0]->r), mpf_get_d(err_mp[1]->r), mpf_get_d(err_mp[2]->r));
*/
clear_point_mp(testPt);
for (i = 0; i < 5; i++) { clear_mp(err_mp[i]); clear_point_mp(test_mp[i]); }
free(test_mp);


s_mp = (point_data_mp *)brealloc(s_mp, (count + 1) * sizeof(point_data_mp));
init_point_data_mp2(&s_mp[count], 0, prec);
point_data_cp_mp(&s_mp[count], &Samples_mp[0]);
count++;

  time = mpf_get_d(Samples_mp[0].time->r);
  for (j = 0; j < size; j++)
  { max = 0; min = 10000; total = 0;
    for (i = 0; i < cycle_num * samples_per_loop; i++)
    {
      sub_mp(diff, &Samples_mp[i].point->coord[j], &average_mp->coord[j]);
      total += norm = d_abs_mp(diff);
      if (norm > max) max = norm;
      if (norm < min) min = norm;
      norm = max - min; if (norm > max_diff) max_diff = norm;
    }
    printf("coord: %d max: %e min: %e ratio: %e average: %e\n", j, max, min, min/max, total / (cycle_num * samples_per_loop));
  }
  printf("time: %e cycle: %d max_diff: %e\n", time, cycle_num, max_diff);
  clear_mp(diff);
}

if (count >= 3)
{ 
  comp_mp diff1, diff2;
  init_mp(diff1); init_mp(diff2);
//printf("z1="); printVec_Matlab_mp(stdout,0,s_mp[count-3].point);
//printf("z2="); printVec_Matlab_mp(stdout,0,s_mp[count-2].point);
//printf("z3="); printVec_Matlab_mp(stdout,0,s_mp[count-1].point);
  for (j = 0; j < size; j++)
  { // compute k for each coord
    sub_mp(diff2, &s_mp[count - 2].point->coord[j], &s_mp[count - 3].point->coord[j]); 
    sub_mp(diff1, &s_mp[count - 1].point->coord[j], &s_mp[count - 2].point->coord[j]); 
    div_mp(diff2, diff1, diff2); 
    norm = d_abs_mp(diff2);
    norm = log10(norm) / log10(sample_factor);
printf("coord: %d k/c: %e k: %e round(k): %d\n", j, norm, norm * cycle_num, (int) (norm * cycle_num + 0.5));
  }
  clear_mp(diff1); clear_mp(diff2);
}

  return;
}

int CauchyEG_actual_amp(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, point_d preEG_approx_d, point_mp preEG_approx_mp, int preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                * 
* NOTES: does the actual Cauchy integral endgame using AMP      *
\***************************************************************/
{
  int i, Samples_prec = 52, retVal = 0, curr_cycle_size = 1;
  int new_approx_prec = 52;
  double closed_loop_tol, min_closed_loop_tol = 1e-150, max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  mpf_t approx_err;
  comp_d endTime_d;
  comp_mp endTime_mp;
  point_d new_approx_d;
  point_mp new_approx_mp;
  point_data_d *Samples_d = NULL;
  point_data_mp *Samples_mp = NULL;

  // initialize MP
  mpf_init(approx_err);
  init_mp(endTime_mp); 
  init_point_d(new_approx_d, 0);
  init_point_mp(new_approx_mp, 0);

  // setup the number of samples per loop
  *samples_per_loop = T->num_PSEG_sample_points + 1;

  // initially allocate for Samples
  Samples_d = (point_data_d *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));
  Samples_mp = (point_data_mp *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = 0; i <= *samples_per_loop * curr_cycle_size; i++)
  {
    init_point_data_d(&Samples_d[i], 0);
    init_point_data_mp(&Samples_mp[i], 0);
  }

  // setup last_approx to be preEG_approx
  *last_approx_prec = preEG_prec;
  if (preEG_prec < 64)
  { // setup last_approx_d
    point_cp_d(last_approx_d, preEG_approx_d);
  }
  else
  { // setup last_approx_mp
    change_prec_point_mp(last_approx_mp, preEG_prec);
    point_cp_mp(last_approx_mp, preEG_approx_mp);
  }

  // find the tolerance to close the first loop
  closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, Start_d, Start_mp, prec_in, T, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // find the first set of samples
  retVal = find_Cauchy_samples_amp(closed_loop_tol, &Samples_d, &Samples_mp, &Samples_prec, cycle_num, *samples_per_loop, &curr_cycle_size, Start_d, Start_mp, prec_in, time_first_increase, T, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // see if we can continue
  if (!retVal)
  { // find new_approx
    find_Cauchy_approx_amp(new_approx_d, new_approx_mp, &new_approx_prec, Samples_d, Samples_mp, Samples_prec, *cycle_num, *samples_per_loop);

    // main loop - keep tracking until we either have failure or we converge
    do
    { // compare the approximations
      findDiff_point(approx_err, new_approx_d, new_approx_mp, new_approx_prec, last_approx_d, last_approx_mp, *last_approx_prec);
      T->error_at_latest_sample_point = mpf_get_d(approx_err);
fprintf(OUT, "cycle: %d error: %e\n", *cycle_num, mpf_get_d(approx_err));

      // check for convergence
      if (T->error_at_latest_sample_point < T->final_tolerance)
      { // we accept our approximation

        // setup Final
        *prec_out = new_approx_prec;
        if (new_approx_prec < 64)
        { // setup Final_d
          point_cp_d(Final_d->point, new_approx_d);
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = *cycle_num;

          T->t_val_at_latest_sample_point = Samples_d[0].time->r;
        }
        else
        { // setup Final_mp
          setprec_point_mp(Final_mp->point, *prec_out);
          setprec_mp(Final_mp->time, *prec_out);

          point_cp_mp(Final_mp->point, new_approx_mp);
          d_to_mp(Final_mp->time, finalTime_d);
          Final_mp->cycle_num = *cycle_num;

          T->t_val_at_latest_sample_point = mpf_get_d(Samples_mp[0].time->r);
        }
        retVal = 0;
        break;
      }
      else if ((Samples_prec < 64 ? d_abs_d(Samples_d[0].time) : d_abs_mp(Samples_mp[0].time)) < T->minTrackT)
      { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

        // setup Final
        *prec_out = new_approx_prec;
        if (new_approx_prec < 64)
        { // setup Final_d
          point_cp_d(Final_d->point, new_approx_d);
          set_d(Final_d->time, Samples_d[0].time);
          Final_d->cycle_num = *cycle_num;

          T->t_val_at_latest_sample_point = Samples_d[0].time->r;
        }
        else
        { // setup Final_mp
          setprec_point_mp(Final_mp->point, *prec_out);
          setprec_mp(Final_mp->time, *prec_out);

          point_cp_mp(Final_mp->point, new_approx_mp);
          set_mp(Final_mp->time, Samples_mp[0].time);
          Final_mp->cycle_num = *cycle_num;

          T->t_val_at_latest_sample_point = mpf_get_d(Samples_mp[0].time->r);
        }
        retVal = retVal_EG_failed_to_converge;
        break;
      }
      else
      { // copy new_approx to last_approx
        *last_approx_prec = new_approx_prec;
        if (new_approx_prec < 64)
        { // copy _d
          point_cp_d(last_approx_d, new_approx_d);
        }
        else
        { // copy _mp
          setprec_point_mp(last_approx_mp, new_approx_prec);
          point_cp_mp(last_approx_mp, new_approx_mp);
        }

        // find the end time for the next sample
        if (Samples_prec < 64)
        { // setup endTime_d
          mul_rdouble_d(endTime_d, Samples_d[0].time, T->power_series_sample_factor);
        }
        else
        { // setup endTime_mp
          setprec_mp(endTime_mp, Samples_prec);
          mul_rdouble_mp(endTime_mp, Samples_mp[0].time, T->power_series_sample_factor);
        }

        // track to the next sample point
        retVal = AMPtrack(&Samples_d[0], &Samples_mp[0], &Samples_prec, time_first_increase, &Samples_d[0], &Samples_mp[0], Samples_prec, endTime_d, endTime_mp, Samples_prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // check for errors
        if (!retVal)
        { // refine the sample point
          if (Samples_prec < 64)
          { // refine using _d
            refine_d_basic(&Samples_d[0], T, OUT, e_d, ED_d, eval_func_d);
          }
          else
          { // refine using _mp
            initMP(Samples_prec);
            change_prec(ED_mp, Samples_prec);
            setprec_eval_struct_mp(*e_mp, Samples_prec);

            refine_mp_basic(&Samples_mp[0], T, OUT, e_mp, ED_mp, eval_func_mp);
          }

          // find the tolerance to close the next loop
          closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &Samples_d[0], &Samples_mp[0], Samples_prec, T, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

          // find the next set of samples
          retVal = find_Cauchy_samples_amp(closed_loop_tol, &Samples_d, &Samples_mp, &Samples_prec, cycle_num, *samples_per_loop, &curr_cycle_size, &Samples_d[0], &Samples_mp[0], Samples_prec, time_first_increase, T, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

          // check for errors
          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_amp(new_approx_d, new_approx_mp, &new_approx_prec, Samples_d, Samples_mp, Samples_prec, *cycle_num, *samples_per_loop);
          }
        }
      }
    } while (!retVal);
  }

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    *prec_out = Samples_prec;
    if (Samples_prec < 64)
    { // setup Final_d
      point_cp_d(Final_d->point, Samples_d[0].point);
      set_d(Final_d->time, Samples_d[0].time);
      Final_d->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = Samples_d[0].time->r;
    }
    else
    { // setup Final_mp
      setprec_point_mp(Final_mp->point, *prec_out);
      setprec_mp(Final_mp->time, *prec_out);

      point_cp_mp(Final_mp->point, Samples_mp[0].point);
      set_mp(Final_mp->time, Samples_mp[0].time);
      Final_mp->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = mpf_get_d(Samples_mp[0].time->r);
    }
  }
  else
  { // setup endLoopSamples
    if (*prec_out < 64)
    { // allocate endLoopSamples_d
      *endLoopSamples_d = (point_data_d *)bmalloc((*samples_per_loop * (*cycle_num) + 1) * sizeof(point_data_d));

      // setup endLoopSamples
      for (i = *samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
      {
        init_point_data_d(&(*endLoopSamples_d)[i], Samples_d[i].point->size);
        point_data_cp_d(&(*endLoopSamples_d)[i], &Samples_d[i]);
      }

      // NULL out endLoopSamples_mp
      endLoopSamples_mp = NULL;
    }
    else
    { // allocate endLoopSamples_mp
      *endLoopSamples_mp = (point_data_mp *)bmalloc((*samples_per_loop * (*cycle_num) + 1) * sizeof(point_data_mp));

      // setup endLoopSamples_mp
      for (i = *samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
      {
        init_point_data_mp(&(*endLoopSamples_mp)[i], Samples_mp[i].point->size);
        point_data_cp_mp(&(*endLoopSamples_mp)[i], &Samples_mp[i]);
      }

      // NULL out endLoopSamples_d
      endLoopSamples_d = NULL;
    }
  }

  // clear memory
  mpf_clear(approx_err);
  clear_mp(endTime_mp);
  clear_point_d(new_approx_d);
  clear_point_mp(new_approx_mp);
  for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
  {
    clear_point_data_d(&Samples_d[i]);
    clear_point_data_mp(&Samples_mp[i]);
  }
  free(Samples_d);
  free(Samples_mp);

  return retVal;
}

int circle_track_amp(int M, point_data_d Final_d[], point_data_mp Final_mp[], int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks around a circle to find M+1 sample points       *
* This will take care of all of the changes in precision and    *
* will always either complete tracking or the path tracking will*
* have completely failed                                        *
\***************************************************************/
{
  int i, retVal = 0;
  double radius, next_step_angle;
  comp_d endTime;

  // error checking - want M >= 3 so that we do not have to predict directly across the origin
  if (M < 3)
  {
    printf("ERROR: The number of sample points (%d) for circle tracking must be >= 3!\n", M);
    bexit(ERROR_CONFIGURATION);
  }

  // setup radius
  if (prec_in < 64)
  { // time in _d
    radius = d_abs_d(Start_d->time);
  }
  else
  { // time in _mp
    radius = d_abs_mp(Start_mp->time);
  }

  // setup local copies to store the sample points
  int *samples_prec = (int *)bmalloc((M + 1) * sizeof(int));
  point_data_d *samples_d = (point_data_d *)bmalloc((M + 1) * sizeof(point_data_d));
  point_data_mp *samples_mp = (point_data_mp *)bmalloc((M + 1) * sizeof(point_data_mp));

  // initialize the first sample
  init_point_data_d(&samples_d[0], 0);
  init_point_data_mp2(&samples_mp[0], 0, MAX(64, prec_in));
  // initialize other sample points
  for (i = 1; i <= M; i++)
  {
    init_point_data_d(&samples_d[i], 0);
    init_point_data_mp(&samples_mp[i], 0);
  }

  // setup samples point [0] - which is Start
  samples_prec[0] = prec_in;
  if (prec_in < 64)
  { // copy _d
    point_data_cp_d(&samples_d[0], Start_d);
  }
  else
  { // copy _mp
    point_data_cp_mp(&samples_mp[0], Start_mp);
  }

  // loop over each sample point
  for (i = 1; i <= M && retVal == 0; i++)
  { // setup endTime
    next_step_angle = 2 * M_PI * (i * 1.0 / M - 1); // next stop at -2 PI + 2 PI i / M
    set_double_d(endTime, radius * cos(next_step_angle), radius * sin(next_step_angle));

    // track from i-1 to i
    retVal = AMPtrack(&samples_d[i], &samples_mp[i], &samples_prec[i], time_first_increase, &samples_d[i-1], &samples_mp[i-1], samples_prec[i-1], endTime, NULL, 52, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }

  if (!retVal)
  { // all were found successfully so we find the maximum precision and set all sample points to that precision
    *prec_out = 0;
    for (i = 0; i <= M; i++)
      if (samples_prec[i] > *prec_out)
        *prec_out = samples_prec[i];

    if (*prec_out < 64)
    { // all samples in double precision - setup Final_d
      for (i = 0; i <= M; i++)
      {
        point_data_cp_d(&Final_d[i], &samples_d[i]);
      }
    }
    else
    { // set all samples to 'prec_out' precision - setup Final_mp
      for (i = 0; i <= M; i++)
      {
        setprec_point_data_mp(&Final_mp[i], *prec_out);

        if (samples_prec[i] < 64)
        {
          convert_point_data_d_to_mp(&Final_mp[i], &samples_d[i]);
        }
        else
        {
          point_data_cp_mp(&Final_mp[i], &samples_mp[i]);
        }
      }
    }
  }

  // clear memory
  for (i = M; i >= 0; i--)
  {
    clear_point_data_d(&samples_d[i]);
    clear_point_data_mp(&samples_mp[i]);
  }
  free(samples_d);
  free(samples_mp);
  free(samples_prec);

  return retVal;
}

void coord_check_amp(double closed_loop_max_error, point_data_d *Pt0_d, point_data_mp *Pt0_mp, int prec0, point_data_d *Pt1_d, point_data_mp *Pt1_mp, int prec1)
{ // check coordinate by coordinate
  int i, size, pr = 0, max_prec = MAX(prec0, prec1);
  int *closed = NULL;

  if (max_prec < 64)
  { // both _d
    comp_d diff;
    size = Pt0_d->point->size;
closed = (int *)bmalloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    { // check to see if loop closed
      sub_d(diff, &Pt0_d->point->coord[i], &Pt1_d->point->coord[i]);
      if (d_abs_d(diff) < closed_loop_max_error)
        pr = closed[i] = 1;
      else
        closed[i] = 0;
    }
  }
  else
  {
    comp_mp diff;
    point_mp pt0, pt1;
    init_mp2(diff, max_prec); init_point_mp2(pt0, 0, max_prec); init_point_mp2(pt1, 0, max_prec);

    if (prec0 < 64) { point_d_to_mp(pt0, Pt0_d->point); } else { point_cp_mp(pt0, Pt0_mp->point); }
    if (prec1 < 64) { point_d_to_mp(pt1, Pt1_d->point); } else { point_cp_mp(pt1, Pt1_mp->point); }

    size = pt0->size;
closed = (int *)bmalloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    { // check to see if loop closed
      sub_mp(diff, &pt0->coord[i], &pt1->coord[i]);
      if (d_abs_mp(diff) < closed_loop_max_error)
        pr = closed[i] = 1;
      else
        closed[i] = 0;
    }
    clear_mp(diff); clear_point_mp(pt0); clear_point_mp(pt1);
  }

if (pr)
{
  printf("closed coords: ");
  for (i = 0; i < size; i++)
    printf("%d ", closed[i]);
  printf("\n");
}

  free(closed);

  return;
}

int find_Cauchy_samples_amp(double closed_loop_max_error, point_data_d **Samples_d, point_data_mp **Samples_mp, int *Samples_prec, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Populates Samples and finds the cycle number to close  *
* the loop                                                      *
\***************************************************************/
{
  int j, tempInt, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_prec = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int)), *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  comp_d time;

  // compute the correct time using double precision
  if (prec_in < 64)
  { // use _d
    set_d(time, Start_d->time);
  }
  else
  { // use _mp
    mp_to_d(time, Start_mp->time);
  }

  // setup Samples[0] by refining Start - make sure that refining works!
  s_digits[0] = ceil(-log10(closed_loop_max_error));
  retVal = refine_digits_amp(T->outputLevel, s_digits[0], &T->latest_newton_residual_d, T->latest_newton_residual_mp, floor(-log10(T->currentNewtonTol)), &(*Samples_d)[0], &(*Samples_mp)[0], &s_prec[0], Start_d, Start_mp, prec_in, time, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // make sure that the refine was successful and then loop until back to the original point
  *cycle_num = 0;
  while (retVal == 0)
  { // go around the origin once
    retVal = circle_track_amp(num_samples, &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], &tempInt, time_first_increase, &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], s_prec[*cycle_num * num_samples], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // save the precision & number of digits correct for these samples
    for (j = 0; j <= num_samples; j++)
    {
      s_prec[*cycle_num * num_samples + j] = tempInt;
      s_digits[*cycle_num * num_samples + j] = floor(-log10(T->currentNewtonTol));
    }

    // increment the cycle_num
    (*cycle_num)++;

    if (retVal)
    { // return error
      fprintf(OUT, "Note: The tracking failed while going around the origin!\n");
      break;
    }
    else
    { // check to see if we have closed the loop
      retVal = check_closed_loop_amp(closed_loop_max_error, &(*Samples_d)[0], &(*Samples_mp)[0], &s_prec[0], &s_digits[0], &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], &s_prec[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time, e_d, e_mp, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

      if (retVal)
      { // error is small enough - exit loop with success
        retVal = 0;
        break;
      }
      else if (*cycle_num > fail_safe_cycle_max)
      { // too many iterations
        fprintf(OUT, "ERROR: Cycle number too high to detect!\n");
        retVal = retVal_cycle_num_too_high;
        break;
      }
      else if (*cycle_num >= *curr_cycle_size)
      { // need to increase the size of memory

        // _prec is easy
        s_prec = (int *)brealloc(s_prec, (2 * (*curr_cycle_size) * num_samples + 1) * sizeof(int));

        // _digits is easy
        s_digits = (int *)brealloc(s_digits, (2 * (*curr_cycle_size) * num_samples + 1) * sizeof(int));

        // _d
        point_data_d *temp_ptr_d = *Samples_d;
        *Samples_d = NULL;
        tempInt = *curr_cycle_size * num_samples;
        *Samples_d = (point_data_d *)bmalloc((2 * (*curr_cycle_size) * num_samples + 1) * sizeof(point_data_d));
        for (j = 2 * (*curr_cycle_size) * num_samples; j >= 0; j--)
        {
          init_point_data_d(&(*Samples_d)[j], 0);

          if (j <= tempInt)
          { // copy temp_ptr_d to Samples
            point_data_cp_d(&(*Samples_d)[j], &temp_ptr_d[j]);
            clear_point_data_d(&temp_ptr_d[j]);
          }
        }
        free(temp_ptr_d);
        temp_ptr_d = NULL;

        // _mp
        point_data_mp *temp_ptr_mp = *Samples_mp;
        *Samples_mp = NULL;
        tempInt = (*curr_cycle_size) * num_samples;
        *Samples_mp = (point_data_mp *)bmalloc((2 * (*curr_cycle_size) * num_samples + 1) * sizeof(point_data_mp));
        for (j = 2 * (*curr_cycle_size) * num_samples; j >= 0; j--)
          if (j <= tempInt)
          { // copy over the orig value to the correct precision
            init_point_data_mp2(&(*Samples_mp)[j], temp_ptr_mp[j].point->size, s_prec[j]);
            point_data_cp_mp(&(*Samples_mp)[j], &temp_ptr_mp[j]);

            // clear temp_ptr_mp
            clear_point_data_mp(&temp_ptr_mp[j]);
          }
          else
          { // just initialize
            init_point_data_mp(&(*Samples_mp)[j], 0);
          }
        free(temp_ptr_mp);
        temp_ptr_mp = NULL;

        // update curr_cycle_size
        *curr_cycle_size = 2 * (*curr_cycle_size);
      }
    }
  }

  if (retVal == 0)
  { // convert all to the maximum precision used and then refine them
    *Samples_prec = 0;
    for (j = 0; j <= *cycle_num * num_samples; j++)
      if (s_prec[j] > *Samples_prec)
        *Samples_prec = s_prec[j];

    if (*Samples_prec < 64)
    { // all are in double precision
      for (j = 0; j <= *cycle_num * num_samples; j++)
      { // refine _d
        refine_d_basic(&(*Samples_d)[j], T, OUT, e_d, ED_d, eval_func_d);
      }
    }
    else
    { // convert all to the correct precision and refine
      initMP(*Samples_prec);
      change_prec(ED_mp, *Samples_prec);
      setprec_eval_struct_mp(*e_mp, *Samples_prec);

      for (j = 0; j <= *cycle_num * num_samples; j++)
      { // convert it to the correct precision
        if (s_prec[j] < 64)
        { // convert _d to _mp
          setprec_point_mp((*Samples_mp)[j].point, *Samples_prec);
          setprec_mp((*Samples_mp)[j].time, *Samples_prec);

          convert_point_data_d_to_mp(&(*Samples_mp)[j], &(*Samples_d)[j]);
        }
        else if (s_prec[j] < *Samples_prec)
        { // increase precision
          change_prec_point_mp((*Samples_mp)[j].point, *Samples_prec);
          change_prec_mp((*Samples_mp)[j].time, *Samples_prec);
        }

        // refine _mp
        refine_mp_basic(&(*Samples_mp)[j], T, OUT, e_mp, ED_mp, eval_func_mp);
      }
    }
  }
  else
  { // setup Samples_prec to the precision of the first one
    *Samples_prec = prec_in;
  }

  free(s_prec);
  free(s_digits);

  return retVal;
}

void find_Cauchy_approx_amp(point_d approx_d, point_mp approx_mp, int *approx_prec, point_data_d Samples_d[], point_data_mp Samples_mp[], int prec_in, int cycle_num, int num_samples)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Finds the Cauchy integral approximation using Trapezoid*
* rule based on cycle_num & num_samples                         *
\***************************************************************/
{
  int i, j, size, sample_num = cycle_num * num_samples; // the first and last sample are equal so we skip the last one

  // set the precision
  *approx_prec = prec_in;

  if (prec_in < 64)
  { // caluclate using _d

    // set size
    size = Samples_d[0].point->size;
    change_size_point_d(approx_d, size);
    approx_d->size = size;

    for (i = 0; i < size; i++)
    { // initialize to 0
      set_zero_d(&approx_d->coord[i]);

      // sum over all of the samples
      for (j = 0; j < sample_num; j++)
      {
        add_d(&approx_d->coord[i], &approx_d->coord[i], &Samples_d[j].point->coord[i]);
      }

      // divide by number of samples added together
      approx_d->coord[i].r /= sample_num;
      approx_d->coord[i].i /= sample_num;
    }
  }
  else
  { // calculate using _mp

    // set size
    size = Samples_mp[0].point->size;
    change_size_point_mp(approx_mp, size);
    approx_mp->size = size;

    // set the precision on approx_mp
    setprec_point_mp(approx_mp, prec_in);

    for (i = 0; i < size; i++)
    { // initialize to 0
      set_zero_mp(&approx_mp->coord[i]);

      // sum over all of the samples
      for (j = 0; j < sample_num; j++)
      {
        add_mp(&approx_mp->coord[i], &approx_mp->coord[i], &Samples_mp[j].point->coord[i]);
      }

      // divide by number of samples added together
      mpf_div_ui(approx_mp->coord[i].r, approx_mp->coord[i].r, sample_num);
      mpf_div_ui(approx_mp->coord[i].i, approx_mp->coord[i].i, sample_num);
    }
  }

  return;
}

int check_closed_loop_amp(double closed_loop_max_error, point_data_d *Pt0_d, point_data_mp *Pt0_mp, int *prec0, int *digitsCorrect0, point_data_d *Pt1_d, point_data_mp *Pt1_mp, int *prec1, int *digitsCorrect1, comp_d time, eval_struct_d *e_d, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - Pt0 == Pt1, 0 - otherwise                  *
* NOTES: check to see if we have closed the loop for Cauchy EG  *
\***************************************************************/
{
  int retVal = 0, min_digits = 0;

  // setup min_digits
  min_digits = ceil(-log10(closed_loop_max_error) + 2.5); // minimum number of digits that we need to have correct for each point

  // check to see if we have the correct number of digits
  if (*digitsCorrect0 < min_digits)
  { // need to refine Pt0
    refine_digits_amp(T->outputLevel, min_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, *digitsCorrect0, Pt0_d, Pt0_mp, prec0, Pt0_d, Pt0_mp, *prec0, time, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // update digitsCorrect0
    *digitsCorrect0 = min_digits;
  }

  if (*digitsCorrect1 < min_digits)
  { // need to refine Pt0
    refine_digits_amp(T->outputLevel, min_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, *digitsCorrect1, Pt1_d, Pt1_mp, prec1, Pt1_d, Pt1_mp, *prec1, time, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // update digitsCorrect1
    *digitsCorrect1 = min_digits;
  }

  // see if we have closed the loop
  retVal = isSamePoint(Pt0_d->point, Pt0_mp->point, *prec0, Pt1_d->point, Pt1_mp->point, *prec1, closed_loop_max_error);
  
  return retVal;
}

double find_closed_loop_tol(double min_tol, double max_tol, point_data_d *x_d, point_data_mp *x_mp, int x_prec, tracker_config_t *T, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: tol                                            *
* NOTES: tol - a tolerance that if two end loop points are      *
* within, then we know that we have closed the loop             *
\***************************************************************/
{
  int n = 0, d_max = MAX(T->AMP_bound_on_degree, 2); // d_max needs to be atleast 2
  double tol = 0, norm_x = 0, L = 0, N = 0, M = 0, K = T->AMP_bound_on_abs_vals_of_coeffs;

  // error checking on min_tol & max_tol
  if (max_tol < min_tol)
    max_tol = min_tol;

  // calculate tol
  if (x_prec < 64)
  { // use double
    double min_sing_d = 0, error_tol_d = 1e-13;

    // setup n
    n = x_d->point->size;

    // setup N
    if (n <= 1)
      N = d_max;
    else
      N = combination(d_max + n - 1, n - 1);

    // setup M
    M = d_max * (d_max - 1) * N;

    // evaluate the function
    eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, x_d->point, x_d->time, ED_d, eval_func_d);

    // calculate minimum singular value of Jv
    min_svd_d(&min_sing_d, e_d->Jv, error_tol_d);

    // calculated norm of x_d
    norm_x = infNormVec_d(x_d->point);

    // calculate L
    L = pow(norm_x, d_max - 2);

    // calculate tol
    tol = K * L * M;
    if (tol == 0) // fail-safe
      tol = min_sing_d; 
    else
      tol = 2 * min_sing_d / tol;
  }
  else
  { // use MP
    int digits = prec_to_digits(x_prec);
    size_t size;
    char *str = NULL;
    mpf_t min_sing_mp, error_tol_mp;

    mpf_init2(min_sing_mp, x_prec);
    mpf_init2(error_tol_mp, x_prec);

    if (T->MPType == 2)
    { // setup precision
      initMP(x_prec);
      change_prec(ED_mp, x_prec);
      setprec_eval_struct_mp(*e_mp, x_prec);
    }

    // setup error_tol
    size = 1 + snprintf(NULL, 0, "1e-%d", digits - 4);
    str = (char *)bmalloc(size * sizeof(char));
    sprintf(str, "1e-%d", digits - 4);
    mpf_set_str(error_tol_mp, str, 10);

    // setup n
    n = x_mp->point->size;

    // setup N
    if (n <= 1)
      N = d_max;
    else
      N = combination(d_max + n - 1, n - 1);

    // setup M
    M = d_max * (d_max - 1) * N;

    // evaluate the function
    eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, x_mp->point, x_mp->time, ED_mp, eval_func_mp);

    // calculate the minium singular value of Jv
    min_svd_mp(min_sing_mp, e_mp->Jv, error_tol_mp);

    // calculated norm of x_mp
    norm_x = infNormVec_mp(x_mp->point);

    // calculate L
    L = pow(norm_x, d_max - 2);

    // calculate tol
    tol  = K * L * M;
    if (tol == 0) // fail-safe
      tol = mpf_get_d(min_sing_mp);
    else
    {
      tol = 2 / tol;
      tol *= mpf_get_d(min_sing_mp);
    }

    mpf_clear(min_sing_mp);
    mpf_clear(error_tol_mp);
    free(str);
  }

  // make sure that tol is between min_tol & max_tol
  if (tol > max_tol) // tol needs to be <= max_tol
    tol = max_tol;
  if (tol < min_tol) // tol needs to be >= min_tol
    tol = min_tol;

  return tol;
}

