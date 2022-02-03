// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

#define NUM_AGREEMENT 3 // number of cycle agreement for 'pre-Cauchy endgame'
#define FAIL_SAFE_MAX 250 // maximum number of loops before giving up
#define C_OVER_K_AGREE 0.75 // agreement for c/k
#define MAX_CAUCHY_RATIO 0.5 // ratio of min / max

int CauchyEG_preEG_d2(double min_closed_loop_tol, double max_closed_loop_tol, point_d preEG_approx, int *cycle_num, point_data_d *Final, comp_d finalTime, point_data_d **CauchySamples, int *curr_cycle_size, int samples_per_loop, point_data_d *Start, int num_agreement, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int preCauchy_loops_d(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d **CauchySamples, point_data_d *PSEGsamples, int *cycle_num, int num_samples, int *curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int CauchyEG_actual_d2(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int samples_per_loop, point_d preEG_approx, comp_d finalTime, point_data_d *CauchySamples, int curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int CauchyEG_preEG_mp2(double min_closed_loop_tol, double max_closed_loop_tol, point_mp preEG_approx, int *cycle_num, point_data_mp *Final, comp_mp finalTime, point_data_mp **CauchySamples, int *curr_cycle_size, int samples_per_loop, point_data_mp *Start, int num_agreement, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int preCauchy_loops_mp(double min_closed_loop_tol, double max_closed_loop_tol, point_data_mp **Samples, point_data_mp *PSEGsamples, int *cycle_num, int num_samples, int *curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int CauchyEG_actual_mp2(double min_closed_loop_tol, double max_closed_loop_tol, point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int samples_per_loop, point_mp preEG_approx, comp_mp finalTime, point_data_mp *CauchySamples, int curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int CauchyEG_preEG_amp2(double min_closed_loop_tol, double max_closed_loop_tol, point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, int *cycle_num, point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d **CauchySamples_d, point_data_mp **CauchySamples_mp, int *CauchySamples_prec, int *curr_cycle_size, int samples_per_loop, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, int num_agreement, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int preCauchy_loops_amp(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d **Samples_d, point_data_mp **Samples_mp, int *Samples_prec, point_data_d *PSEGsamples_d, point_data_mp *PSEGsamples_mp, int *PSEGsamples_prec, int *cycle_num, int num_samples, int *curr_cycle_size, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int CauchyEG_actual_amp2(double min_closed_loop_tol, double max_closed_loop_tol, int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int samples_per_loop, point_d preEG_approx_d, point_mp preEG_approx_mp, int preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *CauchySamples_d, point_data_mp *CauchySamples_mp, int CauchySamples_prec, int curr_cycle_size, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

////// Double Precision ///////

double compute_c_over_k_d(point_d randVec, point_data_d *Pts, int start, double sample_factor, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximates c/k using the sample factor and the 3     *
* sample points {start, start+1, start_2}                       *
\***************************************************************/
{
  int i, size;
  double c_over_k;
  comp_d tempComp, sum1, sum2;

  // find the size
  size = Pts[start].point->size;

  // find a random sum of the subtractions
  set_zero_d(sum1);
  set_zero_d(sum2);
  for (i = 0; i < size; i++)
  { // sum1 += rand * ([start+1] - [start])
    sub_d(tempComp, &Pts[start+1].point->coord[i], &Pts[start].point->coord[i]);
    sum_mul_d(sum1, &randVec->coord[i], tempComp);

    // sum2 += rand * ([start+2] - [start+1])
    sub_d(tempComp, &Pts[start+2].point->coord[i], &Pts[start+1].point->coord[i]);
    sum_mul_d(sum2, &randVec->coord[i], tempComp);
  }

  // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1
  if (d_abs_d(sum1) > tol && d_abs_d(sum2) > tol)
  { // sum2 / sum1
    div_d(tempComp, sum2, sum1);

    // find the norm of sum2 / sum1
    c_over_k = d_abs_d(tempComp);

    // approximate c/k
    c_over_k = log(sample_factor) / log(c_over_k);
  }
  else
  { // c/k = 1
    c_over_k = 1;
  }

  return c_over_k;
}

int check_agreement(int num_agreement, double *c_over_k, double min_agreement)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: see if c/k agrees well enough                          *
\***************************************************************/
{
  int i, retVal = 1;
  double a, b, divide;

  for (i = 1; i < num_agreement && retVal; i++)
  { // divide
    a = fabs(c_over_k[i-1]);
    b = fabs(c_over_k[i]);

    if (a < b)
      divide = a / b;
    else
      divide = b / a;

    // see if we have agreement 
    if (divide < min_agreement)
      retVal = 0;
  }

  return retVal;
}

int CauchyEG_main_d2(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in double precision            *
\***************************************************************/
{
  int i, retVal = 0, num_agreement = NUM_AGREEMENT;
  double min_closed_loop_tol = MAX(1e-10, T->final_tol_times_mult), max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  int curr_cycle_size = 1;
  point_data_d *CauchySamples = NULL;
  eval_struct_d e;

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // initialize
  *samples_per_loop = T->num_PSEG_sample_points + 1;
  init_eval_struct_d(e, 0, 0, 0);
  CauchySamples = (point_data_d *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));
  for (i = 0; i <= *samples_per_loop * curr_cycle_size; i++)
  {
    init_point_data_d(&CauchySamples[i], 0);
  }

  // perform the 'pre-Cauchy EG' - use cycle number approximations & min/max of Cauchy loops to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_d2(min_closed_loop_tol, max_closed_loop_tol, last_approx, cycle_num, Final, finalTime, &CauchySamples, &curr_cycle_size, *samples_per_loop, Start, num_agreement, T, OUT, &e, ED, eval_func);

  if (retVal == 0) 
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_d2(min_closed_loop_tol, max_closed_loop_tol, Final, last_approx, endLoopSamples, cycle_num, *samples_per_loop, last_approx, finalTime, CauchySamples, curr_cycle_size, T, OUT, &e, ED, eval_func, find_dehom);

    // Cauchy samples cleared during Cauchy endgame
  }
  else
  { // clear CauchySamples
    for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    {
      clear_point_data_d(&CauchySamples[i]);
    }
    free(CauchySamples);
  }
 
  clear_eval_struct_d(e);

  return retVal;
}

int CauchyEG_preEG_d2(double min_closed_loop_tol, double max_closed_loop_tol, point_d preEG_approx, int *cycle_num, point_data_d *Final, comp_d finalTime, point_data_d **CauchySamples, int *curr_cycle_size, int samples_per_loop, point_data_d *Start, int num_agreement, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  *
* can efficiently start using the actual Cauchy EG.             *
\***************************************************************/
{
  int i, last_sample = 0, retVal = 0, agreement = 0;
  double compute_tol = 1e-12, *c_over_k = (double *)bmalloc(num_agreement * sizeof(double));
  point_d randVec;
  point_data_d *samples_d = (point_data_d *)bmalloc(3 * sizeof(point_data_d)); // only 3 sample points are needed for approximating c/k
  comp_d endTime;

  // initialize samples_d
  for (i = 0; i < 3; i++)
    init_point_data_d(&samples_d[i], 0);

  // initialize
  *cycle_num = Final->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    c_over_k[i] = -1;

  init_vec_d(randVec, 0);
  make_vec_random_d(randVec, Start->point->size);

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

    // find the first approximation for c/k
    c_over_k[0] = compute_c_over_k_d(randVec, samples_d, 0, T->power_series_sample_factor, compute_tol);

    // loop to find the other cycle numbers
    i = 1;
    while (retVal == 0 && i < num_agreement && d_abs_d(samples_d[2].time) > T->cutoffCycleTime)
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

        // find the next approximation for c/k
        c_over_k[i] = compute_c_over_k_d(randVec, samples_d, 0, T->power_series_sample_factor, compute_tol);
      }

      // increment i;
      i++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && d_abs_d(samples_d[2].time) > T->cutoffCycleTime)
      { // update c/k
        for (i = 1; i < num_agreement; i++)
          c_over_k[i-1] = c_over_k[i];

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
          c_over_k[num_agreement - 1] = compute_c_over_k_d(randVec, samples_d, 0, T->power_series_sample_factor, compute_tol);

          // check for agreement
          agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);
        }
      }
    }
  }

  // check for errors
  if (retVal == 0)
  { // perform Cauchy loops
    do
    { // try to do preCauchy loops
      retVal = preCauchy_loops_d(min_closed_loop_tol, max_closed_loop_tol, CauchySamples, samples_d, cycle_num, samples_per_loop, curr_cycle_size, T, OUT, e_d, ED, eval_func);

      if (retVal == retVal_cycle_num_too_high)
      { // see if we can still continue
        if (d_abs_d(samples_d[2].time) > T->minTrackT)
        { // track to the next sample point

          // update samples_d
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
            // make sure that we do another loop
            retVal = retVal_cycle_num_too_high;
          }
        }
        else
        { // we need to break out of the loop
          break;
        }
      }
    } while (retVal == retVal_cycle_num_too_high);

    // check for errors
    if (retVal == 0)
    { // find the power series approximation described by the 3 sample points
      find_Cauchy_preEG_approx_d(preEG_approx, finalTime, samples_d, 3, *cycle_num, OUT, e_d, ED, eval_func);
    }
    else
    { // set cycle to 0 & setup preEG_approx
      *cycle_num = 0;

      // copy samples_d[last_sample] to preEG_approx
      point_cp_d(preEG_approx, samples_d[last_sample].point);
    }
  }
  else
  { // set cycle to 0 & setup preEG_approx
    *cycle_num = 0;

    // copy samples_d[last_sample] to preEG_approx
    point_cp_d(preEG_approx, samples_d[last_sample].point);
  }

  // update Final
  point_data_cp_d(Final, &samples_d[last_sample]);
  Final->cycle_num = *cycle_num;
  T->t_val_at_latest_sample_point = samples_d[last_sample].time->r;

  // clear memory
  clear_vec_d(randVec);
  for (i = 2; i >= 0; i--)
    clear_point_data_d(&samples_d[i]);
  free(samples_d);
  free(c_over_k);

  return retVal;
}

int compare_Cauchy_ratios_d(point_data_d **Samples, int num_samples, double good_ratio, double tol, double time_cutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - good ratio, 1 - bad ratio                  *
* NOTES: determine if the minimum / maximum ratios are good     *
* enough for convergence                                        *
\***************************************************************/
{
  int i, j, retVal = 0, size = (*Samples)[0].point->size;
  double max, min, norm;

  // see if we are past the time cuttoff
  norm = d_abs_d((*Samples)[0].time);
  if (norm < time_cutoff)
  { // we need to return a good ratio
    retVal = 0;
  }
  else
  { // find the minimum and maximum coordinate-wise
    for (i = 0; i < size && !retVal; i++)
    { // find the minimum and maximum difference
      min = 1e300;
      max = 0;
      for (j = 0; j <= num_samples; j++)
      { // find norm
        norm = norm_sqr_d(&(*Samples)[j].point->coord[i]);
        if (norm > max) 
          max = norm;
        if (norm < min)
          min = norm;
      }
      max = sqrt(max);
      min = sqrt(min);

      // find the ratio if they are large enough
      if (min > tol && max > tol)
      {
        norm = min / max;
        if (norm < good_ratio && max - min > tol)
          retVal = 1; // bad ratio and too far apart
      }
    }
  }

  return retVal;
}

int preCauchy_loops_d(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d **Samples, point_data_d *PSEGsamples, int *cycle_num, int num_samples, int *curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy loops until we have a good ratio of minimum*
* and maximum around the first loop and then compute samples    *
\***************************************************************/
{
  int i, tempInt, contLoop = 1, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  double closed_loop_tol;
  comp_d time;

  // copy Start to Samples[0]
  point_data_cp_d(&(*Samples)[0], &PSEGsamples[2]);

  do
  { // initialize
    *cycle_num = 0;
    set_d(time, (*Samples)[0].time);
    // find the tolerance to close the first loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &(*Samples)[0], NULL, 52, T, e_d, NULL, ED, NULL, eval_func, NULL, NULL);

    // perform the first loop
    retVal = circle_track_d(*Samples, &s_digits[0], &(*Samples)[0], num_samples, T, OUT, e_d, ED, eval_func);

    // increment cycle_num
    (*cycle_num)++;

    // check for error
    if (retVal)
    { // return error
      fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
      contLoop = 0;
    }
    else
    { // find the ratio of the maximum and minimum coordinate-wise for the loop
      contLoop = compare_Cauchy_ratios_d(Samples, num_samples, MAX_CAUCHY_RATIO, T->final_tolerance, T->cutoffRatioTime);

      if (!contLoop)
      { // we have a decent approximation - complete the loops
        do
        { // check to see if we have closed the loop
          if (check_closed_loop_d(closed_loop_tol, &(*Samples)[0], &s_digits[0], &(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time, e_d, T, OUT, ED, eval_func))
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

          // compute next loop
          retVal = circle_track_d(&(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], &(*Samples)[*cycle_num * num_samples], num_samples, T, OUT, e_d, ED, eval_func);

          // increment the cycle_num 
          (*cycle_num)++;

          // check for error
          if (retVal)
          { // return error
            fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
            break;
          }
        } while (1);

        if (retVal == retVal_cycle_num_too_high)
        { // see if we should continue to the next sample point
          if (d_abs_d((*Samples)[0].time) < T->minTrackT)
            contLoop = 0;
          else
            contLoop = 1;
        }
      }

      if (contLoop)
      { // find the time for the next sample point
        mul_rdouble_d(time, (*Samples)[0].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_d(&(*Samples)[0], &(*Samples)[0], time, T, OUT, ED, eval_func);

        // update PSEGsamples
        for (i = 1; i < 3; i++)
        { // update precision
          point_data_cp_d(&PSEGsamples[i-1], &PSEGsamples[i]);
        }
        point_data_cp_d(&PSEGsamples[2], &(*Samples)[0]);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_d_basic(&(*Samples)[0], T, OUT, e_d, ED, eval_func);
        }
        else
        { // quit loop
          contLoop = 0;
        }
      }
    }
  } while (contLoop);

  // release memory
  free(s_digits);

  return retVal;
}

int CauchyEG_actual_d2(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int samples_per_loop, point_d preEG_approx, comp_d finalTime, point_data_d *CauchySamples, int curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual Cauchy integral endgame in double prec *
\***************************************************************/
{
  int i, retVal = 0, dehom_prec = 0;
  double closed_loop_tol, approx_err, norm_last = 0, norm_new = 0;
  comp_d endTime;
  point_d new_approx, dehom;

  init_point_d(new_approx, 0);
  init_point_d(dehom, 0);

  // setup last_approx to be preEG_approx
  point_cp_d(last_approx, preEG_approx);
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom, NULL, &dehom_prec, last_approx, NULL, 52, ED, NULL);
    norm_last = infNormVec_d(dehom);
  }

  // find new_approx
  find_Cauchy_approx_d(new_approx, CauchySamples, *cycle_num, samples_per_loop);
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom, NULL, &dehom_prec, new_approx, NULL, 52, ED, NULL);
    norm_new = infNormVec_d(dehom);
  }

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

      T->t_val_at_latest_sample_point = CauchySamples[0].time->r;

      retVal = 0;
      break;
    }
    else if (d_abs_d(CauchySamples[0].time) < T->minTrackT)
    { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

      // setup Final
      point_cp_d(Final->point, new_approx);
      set_d(Final->time, CauchySamples[0].time);
      Final->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = CauchySamples[0].time->r;

      retVal = retVal_EG_failed_to_converge;
      break;
    }
    else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
    { // we are too large
      // break out of the loop to return error
      retVal = retVal_security_max;
      break;
    }
    else
    { // copy new_approx to last_approx and find the next approximation
      point_cp_d(last_approx, new_approx);
      norm_last = norm_new;

      // loop until we compute the next approximation
      do
      { // find end time for the next sample point
        mul_rdouble_d(endTime, CauchySamples[0].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_d(&CauchySamples[0], &CauchySamples[0], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_d_basic(&CauchySamples[0], T, OUT, e_d, ED, eval_func);

          // find the tolerance to close the next loop
          closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &CauchySamples[0], NULL, 52, T, e_d, NULL, ED, NULL, eval_func, NULL, NULL);

          // find the next set of samples
          retVal = find_Cauchy_samples_d(closed_loop_tol, &CauchySamples, cycle_num, samples_per_loop, &curr_cycle_size, &CauchySamples[0], T, OUT, e_d, ED, eval_func);

          // check for errors
          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_d(new_approx, CauchySamples, *cycle_num, samples_per_loop);
            if (T->securityLevel <= 0)
            { // find dehom & its norm
              find_dehom(dehom, NULL, &dehom_prec, new_approx, NULL, 52, ED, NULL);
              norm_new = infNormVec_d(dehom);
            }
          }
        }

        // make sure that we are not too small
        if (retVal == retVal_cycle_num_too_high && d_abs_d(CauchySamples[0].time) < T->minTrackT)
        { // break out of the while loop with the error
          break;
        }
      } while (retVal == retVal_cycle_num_too_high);
    }
  } while (!retVal);

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    point_cp_d(Final->point, CauchySamples[0].point);
    set_d(Final->time, CauchySamples[0].time);
    Final->cycle_num = *cycle_num;
    T->t_val_at_latest_sample_point = CauchySamples[0].time->r;
  }
  else
  { // allocate endLoopSamples
    *endLoopSamples = (point_data_d *)bmalloc(samples_per_loop * (*cycle_num) * sizeof(point_data_d));

    // setup endLoopSamples
    for (i = samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
    {
      init_point_data_d(&(*endLoopSamples)[i], CauchySamples[i].point->size);
      point_data_cp_d(&(*endLoopSamples)[i], &CauchySamples[i]);
    }
  }

  // release memory
  for (i = samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_d(&CauchySamples[i]);
  free(CauchySamples);
  clear_point_d(new_approx);
  clear_point_d(dehom);

  return retVal;
}

////// Multi Precision ///////

double compute_c_over_k_mp(point_mp randVec, point_data_mp *Pts, int start, double sample_factor, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximates c/k using the sample factor and the 3     *
* sample points {start, start+1, start_2}                       *
\***************************************************************/
{
  int i, size; 
  double c_over_k;
  comp_mp tempComp, sum1, sum2;

  init_mp(tempComp); init_mp(sum1); init_mp(sum2);

  // find the size
  size = Pts[start].point->size;

  // find a random sum of the subtractions
  set_zero_mp(sum1);
  set_zero_mp(sum2);
  for (i = 0; i < size; i++)
  { // sum1 += rand * ([start+1] - [start])
    sub_mp(tempComp, &Pts[start+1].point->coord[i], &Pts[start].point->coord[i]);
    sum_mul_mp(sum1, &randVec->coord[i], tempComp);

    // sum2 += rand * ([start+2] - [start+1])
    sub_mp(tempComp, &Pts[start+2].point->coord[i], &Pts[start+1].point->coord[i]);
    sum_mul_mp(sum2, &randVec->coord[i], tempComp);
  }

  // make sure sum1 > 0 & sum2 > 0 - otherwise the points are constant and so cycle 1
  if (d_abs_mp(sum1) > tol && d_abs_mp(sum2) > tol)
  { // sum2 / sum1
    div_mp(tempComp, sum2, sum1);
  
    // find the norm of sum2 / sum1
    c_over_k = d_abs_mp(tempComp);
  
    // approximate c/k
    c_over_k = log(sample_factor) / log(c_over_k);
  }
  else
  { // c/k = 1
    c_over_k = 1;
  }

  clear_mp(tempComp); clear_mp(sum1); clear_mp(sum2);
  
  return c_over_k;
}

int CauchyEG_main_mp2(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in multi precision             *
\***************************************************************/
{
  int i, retVal = 0, num_agreement = NUM_AGREEMENT;
  double min_closed_loop_tol = MAX(1e-150, pow(10, -(prec_to_digits(T->Precision) - 5))), max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  int curr_cycle_size = 1;
  point_data_mp *CauchySamples = NULL;
  eval_struct_mp e;

  // initialize
  *samples_per_loop = T->num_PSEG_sample_points + 1;
  init_eval_struct_mp(e, 0, 0, 0);
  CauchySamples = (point_data_mp *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = 0; i <= *samples_per_loop * curr_cycle_size; i++)
  {
    init_point_data_mp(&CauchySamples[i], 0);
  }

  // setup last_approx in case of failure
  point_cp_mp(last_approx, Start->point);

  // perform the 'pre-Cauchy EG' - use cycle number approximations & min/max of Cauchy loops to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_mp2(min_closed_loop_tol, max_closed_loop_tol, last_approx, cycle_num, Final, finalTime, &CauchySamples, &curr_cycle_size, *samples_per_loop, Start, num_agreement, T, OUT, &e, ED, eval_func);

  if (retVal == 0)
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_mp2(min_closed_loop_tol, max_closed_loop_tol, Final, last_approx, endLoopSamples, cycle_num, *samples_per_loop, last_approx, finalTime, CauchySamples, curr_cycle_size, T, OUT, &e, ED, eval_func, find_dehom);

    // Cauchy samples cleared during Cauchy endgame
  }
  else
  { // clear CauchySamples
    for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    {
      clear_point_data_mp(&CauchySamples[i]);
    }
    free(CauchySamples);
  }

  // clear
  clear_eval_struct_mp(e);

  return retVal;
}

int CauchyEG_preEG_mp2(double min_closed_loop_tol, double max_closed_loop_tol, point_mp preEG_approx, int *cycle_num, point_data_mp *Final, comp_mp finalTime, point_data_mp **CauchySamples, int *curr_cycle_size, int samples_per_loop, point_data_mp *Start, int num_agreement, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  * 
* can efficiently start using the actual Cauchy EG.             * 
\***************************************************************/
{
  int i, last_sample = 0, retVal = 0, agreement = 0, digits = prec_to_digits(T->Precision);
  double compute_tol, *c_over_k = (double *)bmalloc(num_agreement * sizeof(double));
  point_mp randVec;
  point_data_mp *samples = (point_data_mp *)bmalloc(3 * sizeof(point_data_mp)); // only 3 sample points are needed for a cycle number approximation
  comp_mp endTime;

  init_mp(endTime);
  init_point_mp(randVec, 0);
  for (i = 0; i < 3; i++)
    init_point_data_mp(&samples[i], 0);

  // initialize
  compute_tol = MAX(1e-150, pow(10, -digits + 5));
  *cycle_num = Final->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    c_over_k[i] = -1;

  make_vec_random_mp(randVec, Start->point->size);

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

    // find the first approximation for c/k
    c_over_k[0] = compute_c_over_k_mp(randVec, samples, 0, T->power_series_sample_factor, compute_tol);

    // loop to find the other cycle numbers
    i = 1;
    while (retVal == 0 && i < num_agreement && d_abs_mp(samples[2].time) > T->cutoffCycleTime)
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

        // find the next approximation for c/k
        c_over_k[i] = compute_c_over_k_mp(randVec, samples, 0, T->power_series_sample_factor, compute_tol);
      }

      // increment i;
      i++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && d_abs_mp(samples[2].time) > T->cutoffCycleTime)
      { // update c/k
        for (i = 1; i < num_agreement; i++)
          c_over_k[i-1] = c_over_k[i];

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

          // find the next approximation for c/k
          c_over_k[num_agreement - 1] = compute_c_over_k_mp(randVec, samples, 0, T->power_series_sample_factor, compute_tol);

          // check for agreement
          agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);
        }
      }
    }
  }

  // check for errors
  if (retVal == 0)
  { // perform Cauchy loops
    do 
    { // try to do preCauchy loops
      retVal = preCauchy_loops_mp(min_closed_loop_tol, max_closed_loop_tol, CauchySamples, samples, cycle_num, samples_per_loop, curr_cycle_size, T, OUT, e_mp, ED, eval_func);

      if (retVal == retVal_cycle_num_too_high)
      { // see if we can still continue
        if (d_abs_mp(samples[2].time) > T->minTrackT)
        { // track to the next sample point

          // update samples_mp
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
            // make sure that we do another loop
            retVal = retVal_cycle_num_too_high;
          }
        }
        else
        { // we need to break out of the loop
          break;
        }
      }
    } while (retVal == retVal_cycle_num_too_high);

    // check for errors
    if (retVal == 0)
    { // find the power series approximation described by the 3 sample points
      find_Cauchy_preEG_approx_mp(preEG_approx, finalTime, samples, 3, *cycle_num, OUT, e_mp, ED, eval_func);
    }
    else
    { // set cycle to 0 & setup preEG_approx
      *cycle_num = 0;

      // copy samples[last_sample] to preEG_approx
      point_cp_mp(preEG_approx, samples[last_sample].point);
    }
  }
  else
  { // set cycle to 0 & setup preEG_approx
    *cycle_num = 0;

    // copy samples[last_sample] to preEG_approx
    point_cp_mp(preEG_approx, samples[last_sample].point);
  }

  // update Final
  point_data_cp_mp(Final, &samples[last_sample]);
  Final->cycle_num = *cycle_num;
  T->t_val_at_latest_sample_point = mpf_get_d(samples[last_sample].time->r);

  // clear memory
  clear_vec_mp(randVec);
  clear_mp(endTime);
  for (i = 2; i >= 0; i--)
    clear_point_data_mp(&samples[i]);
  free(samples);
  free(c_over_k);

  return retVal;
}

int compare_Cauchy_ratios_mp(point_data_mp **Samples, int num_samples, double good_ratio, double tol, double time_cutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - good ratio, 1 - bad ratio                  *
* NOTES: determine if the minimum / maximum ratios are good     *
* enough for convergence                                        *
\***************************************************************/
{
  int i, j, retVal = 0, size = (*Samples)[0].point->size;
  double max, min, norm;
  mpf_t normSqr;

  mpf_init(normSqr);

  // see if we are past the time cutoff
  norm = d_abs_mp((*Samples)[0].time);
  if (norm < time_cutoff)
  { // we need to return a good ratio
    retVal = 0;
  }
  else
  { // find the minimum and maximum coordinate-wise
    for (i = 0; i < size && !retVal; i++)
    { // find the minimum and maximum difference
      min = 1e300;
      max = 0;
      for (j = 0; j <= num_samples; j++)
      { // find norm
        norm_sqr_mp(normSqr, &(*Samples)[j].point->coord[i]);
        norm = mpf_get_d(normSqr);
        if (norm > max)
          max = norm;
        if (norm < min)
          min = norm;
      }
      max = sqrt(max);
      min = sqrt(min);

      // find the ratio if they are large enough
      if (min > tol && max > tol)
      {
        norm = min / max;
        if (norm < good_ratio && max - min > tol)
          retVal = 1; // bad ratio and too far apart
      }
    }
  }

  mpf_clear(normSqr);

  return retVal;
}

int preCauchy_loops_mp(double min_closed_loop_tol, double max_closed_loop_tol, point_data_mp **Samples, point_data_mp *PSEGsamples, int *cycle_num, int num_samples, int *curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy loops until we have a good ratio of minimum*
* and maximum around the first loop and then compute samples    *
\***************************************************************/
{
  int i, tempInt, contLoop = 1, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  double closed_loop_tol;
  comp_mp time;

  init_mp(time);

  // copy Start to Samples[0]
  point_data_cp_mp(&(*Samples)[0], &PSEGsamples[2]);

  do
  { // initialize
    *cycle_num = 0;
    set_mp(time, (*Samples)[0].time);
    // find the tolerance to close the first loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, &(*Samples)[0], T->Precision, T, NULL, e_mp, NULL, ED, NULL, eval_func, NULL);

    // perform the first loop
    retVal = circle_track_mp(*Samples, &s_digits[0], &(*Samples)[0], num_samples, T, OUT, e_mp, ED, eval_func);

    // increment cycle_num
    (*cycle_num)++;

    // check for error
    if (retVal)
    { // return error
      fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
      contLoop = 0;
    }
    else
    { // find the ratio of the maximum and minimum coordinate-wise for the loop
      contLoop = compare_Cauchy_ratios_mp(Samples, num_samples, MAX_CAUCHY_RATIO, T->final_tolerance, T->cutoffRatioTime);

      if (!contLoop)
      { // we have a decent approximation - complete the loops
        do
        { // check to see if we have closed the loop
          if (check_closed_loop_mp(closed_loop_tol, &(*Samples)[0], &s_digits[0], &(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time, e_mp, T, OUT, ED, eval_func))
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
            tempInt = *curr_cycle_size * num_samples;
            *Samples = (point_data_mp *)bmalloc((2 * (*curr_cycle_size) * num_samples + 1) * sizeof(point_data_mp));
            for (i = 2 * (*curr_cycle_size) * num_samples; i >= 0; i--)
            {
              init_point_data_mp(&(*Samples)[i], 0);

              if (i <= tempInt)
              { // copy temp_ptr to Samples
                point_data_cp_mp(&(*Samples)[i], &temp_ptr[i]);
                clear_point_data_mp(&temp_ptr[i]);
              }
            }
            free(temp_ptr);
            temp_ptr = NULL;

            // update curr_cycle_size
            *curr_cycle_size = 2 * (*curr_cycle_size);
          }

          // compute next loop
          retVal = circle_track_mp(&(*Samples)[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], &(*Samples)[*cycle_num * num_samples], num_samples, T, OUT, e_mp, ED, eval_func);

          // increment the cycle_num
          (*cycle_num)++;

          // check for error
          if (retVal)
          { // return error
            fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
            break;
          }
        } while (1);

        if (retVal == retVal_cycle_num_too_high)
        { // see if we should continue to the next sample point
          if (d_abs_mp((*Samples)[0].time) < T->minTrackT)
            contLoop = 0;
          else
            contLoop = 1;
        }
      }

      if (contLoop)
      { // track to the next sample point
        mul_rdouble_mp(time, (*Samples)[0].time, T->power_series_sample_factor);
        retVal = track_mp(&(*Samples)[0], &(*Samples)[0], time, T, OUT, ED, eval_func);

        // update PSEGsamples
        for (i = 1; i < 3; i++)
        { // update precision
          point_data_cp_mp(&PSEGsamples[i-1], &PSEGsamples[i]);
        }
        point_data_cp_mp(&PSEGsamples[2], &(*Samples)[0]);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_mp_basic(&(*Samples)[0], T, OUT, e_mp, ED, eval_func);
        }
        else
        { // quit loop
          contLoop = 0;
        }
      }
    }
  } while (contLoop);

  // release memory
  free(s_digits);
  clear_mp(time);

  return retVal;
}

int CauchyEG_actual_mp2(double min_closed_loop_tol, double max_closed_loop_tol, point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int samples_per_loop, point_mp preEG_approx, comp_mp finalTime, point_data_mp *CauchySamples, int curr_cycle_size, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual Cauchy integral endgame in multi prec  *
\***************************************************************/
{
  int i, retVal = 0, dehom_prec = 0;
  double closed_loop_tol, norm_last = 0, norm_new = 0;
  mpf_t approx_err;
  comp_mp endTime;
  point_mp new_approx, dehom;
  
  // initialize MP
  mpf_init(approx_err);
  init_mp(endTime); 
  init_point_mp(new_approx, 0);
  init_point_mp(dehom, 0);

  // setup last_approx to be preEG_approx
  point_cp_mp(last_approx, preEG_approx);
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(NULL, dehom, &dehom_prec, NULL, last_approx, T->Precision, NULL, ED);
    norm_last = infNormVec_mp(dehom);
  }

  // find new_approx
  find_Cauchy_approx_mp(new_approx, CauchySamples, *cycle_num, samples_per_loop);
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(NULL, dehom, &dehom_prec, NULL, new_approx, T->Precision, NULL, ED);
    norm_last = infNormVec_mp(dehom);
  }

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

      T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples[0].time->r);

      retVal = 0;
      break;
    }
    else if (d_abs_mp(CauchySamples[0].time) < T->minTrackT)
    { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

      // setup Final
      point_cp_mp(Final->point, new_approx);
      set_mp(Final->time, CauchySamples[0].time);
      Final->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples[0].time->r);

      retVal = retVal_EG_failed_to_converge;
      break;
    }
    else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
    { // we are too large
      // break out of the loop to return error
      retVal = retVal_security_max;
      break;
    }
    else
    { // copy new_approx to last_approx
      point_cp_mp(last_approx, new_approx);
      norm_last = norm_new;

      // loop until we compute the next approximation
      do
      { // find end time for the next sample point
        mul_rdouble_mp(endTime, CauchySamples[0].time, T->power_series_sample_factor);

        // track to the next sample point
        retVal = track_mp(&CauchySamples[0], &CauchySamples[0], endTime, T, OUT, ED, eval_func);

        // check for errors
        if (!retVal)
        { // refine the sample point
          refine_mp_basic(&CauchySamples[0], T, OUT, e_mp, ED, eval_func);

          // find the tolerance to close the next loop
	  closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, &CauchySamples[0], T->Precision, T, NULL, e_mp, NULL, ED, NULL, eval_func, NULL);

          // find the next set of samples
          retVal = find_Cauchy_samples_mp(closed_loop_tol, &CauchySamples, cycle_num, samples_per_loop, &curr_cycle_size, &CauchySamples[0], T, OUT, e_mp, ED, eval_func);

          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_mp(new_approx, CauchySamples, *cycle_num, samples_per_loop);
            if (T->securityLevel <= 0)
            { // find dehom & its norm
              find_dehom(NULL, dehom, &dehom_prec, NULL, new_approx, T->Precision, NULL, ED);
              norm_new = infNormVec_mp(dehom);
            }
          }
        }

        // make sure that we are not too small
        if (retVal == retVal_cycle_num_too_high && d_abs_mp(CauchySamples[0].time) < T->minTrackT)
        { // break out of the while loop with the error
          break;
        }
      } while (retVal == retVal_cycle_num_too_high);
    }
  } while (!retVal);

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    point_cp_mp(Final->point, CauchySamples[0].point);
    set_mp(Final->time, CauchySamples[0].time);
    Final->cycle_num = *cycle_num;

    T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples[0].time->r);
  }
  else
  { // allocate endLoopSamples
    *endLoopSamples = (point_data_mp *)bmalloc(samples_per_loop * (*cycle_num) * sizeof(point_data_mp));

    // setup endLoopSamples
    for (i = samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
    {
      init_point_data_mp(&(*endLoopSamples)[i], CauchySamples[i].point->size);
      point_data_cp_mp(&(*endLoopSamples)[i], &CauchySamples[i]);
    }
  }

  // clear memory
  mpf_clear(approx_err);
  clear_mp(endTime); 
  clear_point_mp(new_approx); 
  clear_point_mp(dehom);
  for (i = samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_mp(&CauchySamples[i]);
  free(CauchySamples);

  return retVal;
}

////// Adaptive Multi Precision ///////

int CauchyEG_main_amp2(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, double *time_first_increase, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endLoopSamples - sample points at the end of   *
* loop for the successful time - number of them == 'cycle_num'  *
* NOTES: Cauchy integral endgame in double precision            *
\***************************************************************/
{
  int i, retVal = 0, num_agreement = NUM_AGREEMENT;
  double min_closed_loop_tol = 1e-150, max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  int CauchySamples_prec, curr_cycle_size = 1;
  point_data_d *CauchySamples_d = NULL;
  point_data_mp *CauchySamples_mp = NULL;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  // initialize
  *samples_per_loop = T->num_PSEG_sample_points + 1;
  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);
  CauchySamples_d = (point_data_d *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));
  CauchySamples_mp = (point_data_mp *)bmalloc((*samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = 0; i <= *samples_per_loop * curr_cycle_size; i++)
  {
    init_point_data_d(&CauchySamples_d[i], 0);
    init_point_data_mp(&CauchySamples_mp[i], 0);
  }

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

  // perform the 'pre-Cauchy EG' - use cycle number approximations & min/max of Cauchy loops to determine when we should actually start using Cauchy EG
  retVal = CauchyEG_preEG_amp2(min_closed_loop_tol, max_closed_loop_tol, last_approx_d, last_approx_mp, last_approx_prec, cycle_num, Final_d, Final_mp, prec_out, finalTime_d, finalTime_mp, &CauchySamples_d, &CauchySamples_mp, &CauchySamples_prec, &curr_cycle_size, *samples_per_loop, Start_d, Start_mp, prec_in, num_agreement, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  if (retVal == 0)
  { // we can now do the actual Cauchy endgame
    retVal = CauchyEG_actual_amp2(min_closed_loop_tol, max_closed_loop_tol, prec_out, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, endLoopSamples_d, endLoopSamples_mp, cycle_num, *samples_per_loop, last_approx_d, last_approx_mp, *last_approx_prec, finalTime_d, finalTime_mp, CauchySamples_d, CauchySamples_mp, CauchySamples_prec, curr_cycle_size, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    // Cauchy samples cleared during Cauchy endgame
  }
  else
  { // clear CauchySamples
    for (i = *samples_per_loop * curr_cycle_size; i >= 0; i--)
    {
      clear_point_data_d(&CauchySamples_d[i]);
      clear_point_data_mp(&CauchySamples_mp[i]);
    }
    free(CauchySamples_d);
    free(CauchySamples_mp);
  }
 
  // clear
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  return retVal;
}

int CauchyEG_preEG_amp2(double min_closed_loop_tol, double max_closed_loop_tol, point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, int *cycle_num, point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d **CauchySamples_d, point_data_mp **CauchySamples_mp, int *CauchySamples_prec, int *curr_cycle_size, int samples_per_loop, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, int num_agreement, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Uses cycle number approximations to determine when we  *
* can efficiently starting using the actual Cauchy EG.  Used    *
\***************************************************************/
{ 
  int i, size, curr_sample, last_sample = 0, retVal = 0, agreement = 0;
  double compute_tol, *c_over_k = (double *)bmalloc(num_agreement * sizeof(double));
  int *samples_prec = (int *)bmalloc(3 * sizeof(int));
  point_data_d *samples_d = (point_data_d *)bmalloc(3 * sizeof(point_data_d)); // only 3 sample points are needed for a cycle number approximation
  point_data_mp *samples_mp = (point_data_mp *)bmalloc(3 * sizeof(point_data_mp));
  point_d randVec_d;
  point_mp randVec_mp;
  mpq_t **randVec_rat;
  comp_d endTime_d;
  comp_mp endTime_mp;

  init_mp(endTime_mp);
  for (i = 0; i < 3; i++)
  {
    init_point_data_d(&samples_d[i], 0);
    init_point_data_mp(&samples_mp[i], 0);
  }

  size = prec_in < 64 ? Start_d->point->size : Start_mp->point->size;
  init_vec_d(randVec_d, size);
  init_vec_mp(randVec_mp, size);
  init_vec_rat(randVec_rat, size);

  make_vec_random_rat(randVec_d, randVec_mp, randVec_rat, size, T->Precision, T->AMP_max_prec, 0, 0);

  // initialize
  compute_tol = 1e-12;
  *cycle_num = Final_d->cycle_num = Final_mp->cycle_num = 0;
  for (i = 0; i < num_agreement; i++)
    c_over_k[i] = -1;

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
    setprec_point_data_mp(&samples_mp[0], samples_prec[0]);

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

      // find approximation for c/k using double precision
      compute_tol = 1e-12;
      c_over_k[0] = compute_c_over_k_d(randVec_d, samples_d, 0, T->power_series_sample_factor, compute_tol);
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
          setprec_point_data_mp(&samples_mp[i], samples_prec[2]);

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

      // find approximation for c/k using mulit precision
      size = prec_to_digits(samples_prec[2]);
      compute_tol = MAX(1e-150, pow(10, -size + 5));
      change_prec_point_mp_rat(randVec_mp, samples_prec[2], randVec_rat);
      c_over_k[0] = compute_c_over_k_mp(randVec_mp, samples_mp, 0, T->power_series_sample_factor, compute_tol);
    }

    // loop to find the other cycle numbers
    curr_sample = 1;
    while (retVal == 0 && curr_sample < num_agreement && (samples_prec[2] < 64 ? d_abs_d(samples_d[2].time) > T->cutoffCycleTime : d_abs_mp(samples_mp[2].time) > T->cutoffCycleTime))
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
          setprec_point_data_mp(&samples_mp[i-1], samples_prec[i]);

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

          // find approximation for c/k using double precision
          compute_tol = 1e-12;
          c_over_k[curr_sample] = compute_c_over_k_d(randVec_d, samples_d, 0, T->power_series_sample_factor, compute_tol);
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
              setprec_point_data_mp(&samples_mp[i], samples_prec[2]);

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

          // find approximation for c/k using mulit precision
          size = prec_to_digits(samples_prec[2]);
          compute_tol = MAX(1e-150, pow(10, -size + 5));
          change_prec_point_mp_rat(randVec_mp, samples_prec[2], randVec_rat);
          c_over_k[curr_sample] = compute_c_over_k_mp(randVec_mp, samples_mp, 0, T->power_series_sample_factor, compute_tol);
        }
      }

      // increment curr_sample
      curr_sample++;
    }

    // see if we should continue
    if (!retVal)
    { // check for agreement
      agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);

      // loop until we have an error, we have agreement OR we reach timeCutoff
      while (retVal == 0 && !agreement && (samples_prec[2] < 64 ? d_abs_d(samples_d[2].time) > T->cutoffCycleTime : d_abs_mp(samples_mp[2].time) > T->cutoffCycleTime))
      { // update c/k
        for (i = 1; i < num_agreement; i++)
          c_over_k[i-1] = c_over_k[i];

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
            setprec_point_data_mp(&samples_mp[i-1], samples_prec[i]);

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

            // find approximation for c/k using double precision
            compute_tol = 1e-12;
            c_over_k[num_agreement - 1] = compute_c_over_k_d(randVec_d, samples_d, 0, T->power_series_sample_factor, compute_tol);
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
                setprec_point_data_mp(&samples_mp[i], samples_prec[2]);

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

            // find approximation for c/k using mulit precision
            size = prec_to_digits(samples_prec[2]);
            compute_tol = MAX(1e-150, pow(10, -size + 5));
            change_prec_point_mp_rat(randVec_mp, samples_prec[2], randVec_rat);
            c_over_k[num_agreement - 1] = compute_c_over_k_mp(randVec_mp, samples_mp, 0, T->power_series_sample_factor, compute_tol);
          }

          // check for agreement
          agreement = check_agreement(num_agreement, c_over_k, C_OVER_K_AGREE);
        }
      }
    }
  }

  // check for errors
  if (retVal == 0)
  { // perform Cauchy loops
    do
    { // try to do preCauchy loops
      retVal = preCauchy_loops_amp(min_closed_loop_tol, max_closed_loop_tol, CauchySamples_d, CauchySamples_mp, CauchySamples_prec, samples_d, samples_mp, samples_prec, cycle_num, samples_per_loop, curr_cycle_size, time_first_increase, T, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

      if (retVal == retVal_cycle_num_too_high || retVal == retVal_refining_failed)
      { // see if we can still continue
        if ((samples_prec[2] < 64 && d_abs_d(samples_d[2].time) > T->minTrackT) || (samples_prec[2] >= 64 && d_abs_mp(samples_mp[2].time) > T->minTrackT))
        { // track to the next sample point

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
              setprec_point_data_mp(&samples_mp[i-1], samples_prec[i]);
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
                  setprec_point_data_mp(&samples_mp[i], samples_prec[2]);
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
            }
            // make sure that we do another loop
            retVal = retVal_cycle_num_too_high;
          }
        }
        else
        { // we need to break out of the loop
          break;
        }
      }
    } while (retVal == retVal_cycle_num_too_high);

    // check for errors
    if (retVal == 0)
    { // find the power series approximation described by the 3 sample points
      find_Cauchy_preEG_approx_amp(preEG_approx_d, preEG_approx_mp, preEG_prec, finalTime_d, finalTime_mp, samples_d, samples_mp, samples_prec, 3, *cycle_num, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }
    else
    { // set cycle to 0 & setup preEG_approx
      *cycle_num = 0;

      *preEG_prec = samples_prec[last_sample];
      if (*preEG_prec < 64)
      { // copy samples_d[last_sample] to preEG_approx_d
        point_cp_d(preEG_approx_d, samples_d[last_sample].point);
      }
      else
      { // copy samples_mp[last_sample] to preEG_appprox_mp
        setprec_point_mp(preEG_approx_mp, *preEG_prec);
        point_cp_mp(preEG_approx_mp, samples_mp[last_sample].point);
      }
    }
  }
  else
  { // set cycle to 0 & setup preEG_approx
    *cycle_num = 0;

    *preEG_prec = samples_prec[last_sample];
    if (*preEG_prec < 64)
    { // copy samples_d[last_sample] to preEG_approx_d
      point_cp_d(preEG_approx_d, samples_d[last_sample].point);
    }
    else
    { // copy samples_mp[last_sample] to preEG_appprox_mp
      setprec_point_mp(preEG_approx_mp, *preEG_prec);
      point_cp_mp(preEG_approx_mp, samples_mp[last_sample].point);
    }
  }

  // update Final
  *prec_out = samples_prec[last_sample];
  if (*prec_out < 64)
  { // copy _d
    point_data_cp_d(Final_d, &samples_d[last_sample]);
    Final_d->cycle_num = *cycle_num;
    T->t_val_at_latest_sample_point = samples_d[last_sample].time->r;
  }
  else
  { // copy _mp
    setprec_point_data_mp(Final_mp, *prec_out);

    point_data_cp_mp(Final_mp, &samples_mp[last_sample]);
    Final_mp->cycle_num = *cycle_num;
    T->t_val_at_latest_sample_point = mpf_get_d(samples_mp[last_sample].time->r);
  }

  // clear memory
  clear_vec(randVec_d, randVec_mp, randVec_rat, T->MPType);
  clear_mp(endTime_mp);
  free(c_over_k);
  for (i = 2; i >= 0; i--)
  {
    clear_point_data_d(&samples_d[i]);
    clear_point_data_mp(&samples_mp[i]);
  }
  free(samples_mp);
  free(samples_d);
  free(samples_prec);

  return retVal;
}

int preCauchy_loops_amp(double min_closed_loop_tol, double max_closed_loop_tol, point_data_d **Samples_d, point_data_mp **Samples_mp, int *Samples_prec, point_data_d *PSEGsamples_d, point_data_mp *PSEGsamples_mp, int *PSEGsamples_prec, int *cycle_num, int num_samples, int *curr_cycle_size, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int)) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy loops until we have a good ratio of minimum*
* and maximum around the first loop and then compute samples    *
\***************************************************************/
{
  int j, tempInt, contLoop = 1, retVal = 0, fail_safe_cycle_max = MAX(FAIL_SAFE_MAX, T->cycle_num_max);
  int *s_prec = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int)), *s_digits = (int *)bmalloc((*curr_cycle_size * num_samples + 1) * sizeof(int));
  double closed_loop_tol;
  comp_d time_d;
  comp_mp time_mp;

  init_mp(time_mp);

  // setup the starting point
  *Samples_prec = PSEGsamples_prec[2];
  if (*Samples_prec < 64)
  { // copy to Samples_d
    point_data_cp_d(&(*Samples_d)[0], &PSEGsamples_d[2]);
  }
  else
  { // copy to Samples_mp
    setprec_point_data_mp(&(*Samples_mp)[0], *Samples_prec);
    point_data_cp_mp(&(*Samples_mp)[0], &PSEGsamples_mp[2]);
  }

  do
  { // initialize
    *cycle_num = 0;
    if (*Samples_prec < 64)
    { // use _d
      set_d(time_d, (*Samples_d)[0].time);
      d_to_mp(time_mp, time_d);
    }
    else
    { // use _mp
      setprec_mp(time_mp, *Samples_prec);
      set_mp(time_mp, (*Samples_mp)[0].time);
      mp_to_d(time_d, time_mp);
    }
  
    // find the tolerance to close the first loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &(*Samples_d)[0], &(*Samples_mp)[0], *Samples_prec, T, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // refine Samples[0] to the correct number of digits
    s_digits[0] = ceil(-log10(closed_loop_tol));
    retVal = refine_digits_amp(T->outputLevel, s_digits[0], &T->latest_newton_residual_d, T->latest_newton_residual_mp, floor(-log10(T->currentNewtonTol)), &(*Samples_d)[0], &(*Samples_mp)[0], &s_prec[0], &(*Samples_d)[0], &(*Samples_mp)[0], *Samples_prec, time_d, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // make sure that refining worked
    if (!retVal)
    { // perform the first loop
      retVal = circle_track_amp(num_samples, *Samples_d, *Samples_mp, &tempInt, time_first_increase, &(*Samples_d)[0], &(*Samples_mp)[0], s_prec[0], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

      // increment cycle_num
      (*cycle_num)++;
    }

    // check for errors
    if (retVal)
    { // return error
      fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
      contLoop = 0;
    }
    else
    { // save the precision & number of digits correct for these samples
      for (j = 0; j <= num_samples; j++)
      {
        s_prec[j] = tempInt;
        s_digits[j] = floor(-log10(T->currentNewtonTol));
      }
   
      // find the ratio of the maximum and minimum coordinate-wise for the loop
      if (tempInt < 64)
        contLoop = compare_Cauchy_ratios_d(Samples_d, num_samples, MAX_CAUCHY_RATIO, T->final_tolerance, T->cutoffRatioTime);
      else 
        contLoop = compare_Cauchy_ratios_mp(Samples_mp, num_samples, MAX_CAUCHY_RATIO, T->final_tolerance, T->cutoffRatioTime);

      if (!contLoop)
      { // we have a decent approximation - complete the loops
        do
        { // check to see if we have closed the loop
          if (check_closed_loop_amp(closed_loop_tol, &(*Samples_d)[0], &(*Samples_mp)[0], &s_prec[0], &s_digits[0], &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], &s_prec[*cycle_num * num_samples], &s_digits[*cycle_num * num_samples], time_d, e_d, e_mp, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec))
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

          // compute next loop
          retVal = circle_track_amp(num_samples, &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], &tempInt, time_first_increase, &(*Samples_d)[*cycle_num * num_samples], &(*Samples_mp)[*cycle_num * num_samples], s_prec[*cycle_num * num_samples], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

          // save the precision & number of digits correct for these samples
          for (j = 0; j <= num_samples; j++)
          {
            s_prec[*cycle_num * num_samples + j] = tempInt;
            s_digits[*cycle_num * num_samples + j] = floor(-log10(T->currentNewtonTol));
          }

          // increment the cycle_num
          (*cycle_num)++;

          // check for error
          if (retVal)
          { // return error
            fprintf(OUT, "Note: The tracking failed while going around the origin (%d)!\n", retVal);
            break;
          }
        } while (1);

        if (retVal == retVal_cycle_num_too_high)
        { // see if we should continue to the next sample point
          if ((s_prec[0] < 64 && d_abs_d((*Samples_d)[0].time) < T->cutoffRatioTime) || (s_prec[0] >= 64 && d_abs_mp((*Samples_mp)[0].time) < T->cutoffRatioTime))
            contLoop = 0;
          else
            contLoop = 1;
        }
      }

      if (contLoop)
      { // track to the next sample point
        if (*Samples_prec < 64)
        { // setup time_d
          mul_rdouble_d(time_d, (*Samples_d)[0].time, T->power_series_sample_factor);
        }
        else
        { // setup time_mp
          setprec_mp(time_mp, *Samples_prec);
          mul_rdouble_mp(time_mp, (*Samples_mp)[0].time, T->power_series_sample_factor);
        }

        retVal = AMPtrack(&(*Samples_d)[0], &(*Samples_mp)[0], Samples_prec, time_first_increase, &(*Samples_d)[0], &(*Samples_mp)[0], *Samples_prec, time_d, time_mp, *Samples_prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // update PSEGsamples
        for (j = 1; j < 3; j++)
        { // update precision
          PSEGsamples_prec[j-1] = PSEGsamples_prec[j];
          // copy samples
          if (PSEGsamples_prec[j] < 64)
          { // copy _d
            point_data_cp_d(&PSEGsamples_d[j-1], &PSEGsamples_d[j]);
          }
          else
          { // copy _mp
            setprec_point_data_mp(&PSEGsamples_mp[j-1], PSEGsamples_prec[j]);
            point_data_cp_mp(&PSEGsamples_mp[j-1], &PSEGsamples_mp[j]);
          }
        }

        // setup PSEGsamples[2]
        PSEGsamples_prec[2] = *Samples_prec;
        if (*Samples_prec < 64)
        { // copy to PSEGsamples_d[2]
          point_data_cp_d(&PSEGsamples_d[2], &(*Samples_d)[0]);
        }
        else
        { // copy to PSEGsamples_mp[2]
          setprec_point_data_mp(&PSEGsamples_mp[2], PSEGsamples_prec[2]);
          point_data_cp_mp(&PSEGsamples_mp[2], &(*Samples_mp)[0]);
        }          

        // check for errors
        if (retVal)
        { // quit loop
          contLoop = 0;
        }
      }
    } 
  } while (contLoop);

  if (retVal == 0)
  { // convert all to the maximum precision used and then refine them
    *Samples_prec = 0;
    for (j = 0; j <= *cycle_num * num_samples; j++)
      if (s_prec[j] > *Samples_prec)
        *Samples_prec = s_prec[j];

    if (*Samples_prec < 64)
    { // all are in double precision
      for (j = 0; j <= *cycle_num * num_samples; j++)
      { // refine_d
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
          setprec_point_data_mp(&(*Samples_mp)[j], *Samples_prec);
          convert_point_data_d_to_mp(&(*Samples_mp)[j], &(*Samples_d)[j]);
        }
        else if (s_prec[j] < *Samples_prec)
        { // increase precision
          change_prec_point_data_mp(&(*Samples_mp)[j], *Samples_prec);
        }

        // refine _mp 
        refine_mp_basic(&(*Samples_mp)[j], T, OUT, e_mp, ED_mp, eval_func_mp);
      }
    }
  }

  // release memory
  free(s_prec);
  free(s_digits);
  clear_mp(time_mp);

  return retVal;  
}

int CauchyEG_actual_amp2(double min_closed_loop_tol, double max_closed_loop_tol, int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int samples_per_loop, point_d preEG_approx_d, point_mp preEG_approx_mp, int preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *CauchySamples_d, point_data_mp *CauchySamples_mp, int CauchySamples_prec, int curr_cycle_size, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                * 
* NOTES: does the actual Cauchy integral endgame using AMP      *
\***************************************************************/
{
  int i, retVal = 0, new_approx_prec = 52, dehom_prec = 52;
  double closed_loop_tol, norm_last = 0, norm_new = 0;
  mpf_t approx_err;
  comp_d endTime_d;
  comp_mp endTime_mp;
  point_d new_approx_d, dehom_d;
  point_mp new_approx_mp, dehom_mp;

  // initialize MP
  mpf_init(approx_err);
  init_mp(endTime_mp); 
  init_point_d(new_approx_d, 0); init_point_d(dehom_d, 0);
  init_point_mp(new_approx_mp, 0); init_point_mp(dehom_mp, 0);

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
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom_d, dehom_mp, &dehom_prec, last_approx_d, last_approx_mp, *last_approx_prec, ED_d, ED_mp);
    if (dehom_prec < 64)
      norm_last = infNormVec_d(dehom_d);  
    else
      norm_last = infNormVec_mp(dehom_mp);  
  }

  // find new_approx based on the Cauchy samples
  find_Cauchy_approx_amp(new_approx_d, new_approx_mp, &new_approx_prec, CauchySamples_d, CauchySamples_mp, CauchySamples_prec, *cycle_num, samples_per_loop);
  if (T->securityLevel <= 0)
  { // find dehom & its norm
    find_dehom(dehom_d, dehom_mp, &dehom_prec, new_approx_d, new_approx_mp, new_approx_prec, ED_d, ED_mp);
    if (dehom_prec < 64)
      norm_new = infNormVec_d(dehom_d);
    else
      norm_new = infNormVec_mp(dehom_mp);
  }

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

        T->t_val_at_latest_sample_point = CauchySamples_d[0].time->r;
      }
      else
      { // setup Final_mp
        setprec_point_data_mp(Final_mp, *prec_out);
        point_cp_mp(Final_mp->point, new_approx_mp);
        d_to_mp(Final_mp->time, finalTime_d);
        Final_mp->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples_mp[0].time->r);
      }
      retVal = 0;
      break;
    }
    else if ((CauchySamples_prec < 64 ? d_abs_d(CauchySamples_d[0].time) : d_abs_mp(CauchySamples_mp[0].time)) < T->minTrackT)
    { // we are too close to t = 0 but we do not have the correct tolerance - so we exit

      // setup Final
      *prec_out = new_approx_prec;
      if (new_approx_prec < 64)
      { // setup Final_d
        point_cp_d(Final_d->point, new_approx_d);
        set_d(Final_d->time, CauchySamples_d[0].time);
        Final_d->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = CauchySamples_d[0].time->r;
      }
      else
      { // setup Final_mp
        setprec_point_data_mp(Final_mp, *prec_out);
        point_cp_mp(Final_mp->point, new_approx_mp);
        set_mp(Final_mp->time, CauchySamples_mp[0].time);
        Final_mp->cycle_num = *cycle_num;

        T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples_mp[0].time->r);
      }
      retVal = retVal_EG_failed_to_converge;
      break;
    }
    else if (T->securityLevel <= 0 && norm_last > T->securityMaxNorm && norm_new > T->securityMaxNorm)
    { // we are too large
      // break out of the loop to return error
      retVal = retVal_security_max;
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
      norm_last = norm_new;

      // loop until we compute the next approximation
      do
      { // find the end time for the next sample
        if (CauchySamples_prec < 64)
        { // setup endTime_d
          mul_rdouble_d(endTime_d, CauchySamples_d[0].time, T->power_series_sample_factor);
        }
        else
        { // setup endTime_mp
          setprec_mp(endTime_mp, CauchySamples_prec);
          mul_rdouble_mp(endTime_mp, CauchySamples_mp[0].time, T->power_series_sample_factor);
        }

        // track to the next sample point
        retVal = AMPtrack(&CauchySamples_d[0], &CauchySamples_mp[0], &CauchySamples_prec, time_first_increase, &CauchySamples_d[0], &CauchySamples_mp[0], CauchySamples_prec, endTime_d, endTime_mp, CauchySamples_prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // check for errors
        if (!retVal)
        { // refine the sample point
          if (CauchySamples_prec < 64)
          { // refine using _d
            refine_d_basic(&CauchySamples_d[0], T, OUT, e_d, ED_d, eval_func_d);
          }
          else
          { // refine using _mp
            initMP(CauchySamples_prec);
            change_prec(ED_mp, CauchySamples_prec);
            setprec_eval_struct_mp(*e_mp, CauchySamples_prec);
  
            refine_mp_basic(&CauchySamples_mp[0], T, OUT, e_mp, ED_mp, eval_func_mp);
          }

          // find the tolerance to close the next loop
          closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &CauchySamples_d[0], &CauchySamples_mp[0], CauchySamples_prec, T, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

          // find the next set of samples
          retVal = find_Cauchy_samples_amp(closed_loop_tol, &CauchySamples_d, &CauchySamples_mp, &CauchySamples_prec, cycle_num, samples_per_loop, &curr_cycle_size, &CauchySamples_d[0], &CauchySamples_mp[0], CauchySamples_prec, time_first_increase, T, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

          if (!retVal)
          { // find the next approximation
            find_Cauchy_approx_amp(new_approx_d, new_approx_mp, &new_approx_prec, CauchySamples_d, CauchySamples_mp, CauchySamples_prec, *cycle_num, samples_per_loop);
            if (T->securityLevel <= 0)
            { // find dehom & its norm
              find_dehom(dehom_d, dehom_mp, &dehom_prec, new_approx_d, new_approx_mp, new_approx_prec, ED_d, ED_mp);
              if (dehom_prec < 64)
                norm_new = infNormVec_d(dehom_d);
              else
                norm_new = infNormVec_mp(dehom_mp);
            }
          }
        }

        // make sure that we are not too small
        if (retVal == retVal_cycle_num_too_high && ((CauchySamples_prec < 64 && d_abs_d(CauchySamples_d[0].time) < T->minTrackT) || (CauchySamples_prec >= 64 && d_abs_mp(CauchySamples_mp[0].time) < T->minTrackT)))
        { // break out of the while loop with the error
          break;
        }
      } while (retVal == retVal_cycle_num_too_high);
    }
  } while (!retVal);

  if (retVal != 0 && retVal != retVal_EG_failed_to_converge)
  { // we need to setup Final since we had an error
    *prec_out = CauchySamples_prec;
    if (CauchySamples_prec < 64)
    { // setup Final_d
      point_cp_d(Final_d->point, CauchySamples_d[0].point);
      set_d(Final_d->time, CauchySamples_d[0].time);
      Final_d->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = CauchySamples_d[0].time->r;
    }
    else
    { // setup Final_mp
      setprec_point_data_mp(Final_mp, *prec_out);
      point_data_cp_mp(Final_mp, &CauchySamples_mp[0]);
      Final_mp->cycle_num = *cycle_num;

      T->t_val_at_latest_sample_point = mpf_get_d(CauchySamples_mp[0].time->r);
    }
  }
  else
  { // setup endLoopSamples
    if (*prec_out < 64)
    { // allocate endLoopSamples_d
      *endLoopSamples_d = (point_data_d *)bmalloc((samples_per_loop * (*cycle_num) + 1) * sizeof(point_data_d));

      // setup endLoopSamples
      for (i = samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
      {
        init_point_data_d(&(*endLoopSamples_d)[i], CauchySamples_d[i].point->size);
        point_data_cp_d(&(*endLoopSamples_d)[i], &CauchySamples_d[i]);
      }

      // NULL out endLoopSamples_mp
      endLoopSamples_mp = NULL;
    }
    else
    { // allocate endLoopSamples_mp
      *endLoopSamples_mp = (point_data_mp *)bmalloc((samples_per_loop * (*cycle_num) + 1) * sizeof(point_data_mp));

      // setup endLoopSamples_mp
      for (i = samples_per_loop * (*cycle_num) - 1; i >= 0; i--)
      {
        init_point_data_mp(&(*endLoopSamples_mp)[i], CauchySamples_mp[i].point->size);
        point_data_cp_mp(&(*endLoopSamples_mp)[i], &CauchySamples_mp[i]);
      }

      // NULL out endLoopSamples_d
      endLoopSamples_d = NULL;
    }
  }

  // clear memory
  mpf_clear(approx_err);
  clear_mp(endTime_mp);
  clear_point_d(new_approx_d); clear_point_d(dehom_d);
  clear_point_mp(new_approx_mp); clear_point_mp(dehom_mp);
  for (i = samples_per_loop * curr_cycle_size; i >= 0; i--)
  {
    clear_point_data_d(&CauchySamples_d[i]);
    clear_point_data_mp(&CauchySamples_mp[i]);
  }
  free(CauchySamples_d);
  free(CauchySamples_mp);

  return retVal;
}


