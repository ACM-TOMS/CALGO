// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

////////// USE CAUCHY ENDGAME TO FIND EITHER (0) IF IT IS RANK DEFICIENT OR (1) ITS CORANK ////////////////

// rankType == 0 means only find if it is rank deficient or not,
// otherwise, find its actual corank
#define MIN_SV_CHECK 1e-6 // minimum singular value that is acceptable for the quick non-singular check

// _d functions
int Cauchy_rankDef_d(double *CN, point_d approx0, point_d approx1, comp_d finalTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int Cauchy_nonsing_check_d(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

// _mp functions
int Cauchy_rankDef_mp(double *CN, point_mp approx0, point_mp approx1, comp_mp finalTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int Cauchy_nonsing_check_mp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// _amp functions
int Cauchy_rankDef_amp(double *CN, mpf_t minSV0, int ratioTol, int prec0, mpf_t minSV1, int prec1, mpf_t mat_norm, FILE *OUT);
int Cauchy_nonsing_check_amp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, point_data_d *Final_d, point_data_mp *Final_mp, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

/////// Double Precision ////////

int CauchyEG_rank_d(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in double precision along *
* with rank determination                                       *
\***************************************************************/
{
  int i, retVal, cycle, samples_per_loop;
  comp_d finalTime;
  point_data_d *endSamples = NULL;

  // setup finalTime
  set_double_d(finalTime, T->targetT, 0);

  // setup last_approx in case of failure
  point_cp_d(last_approx, Start->point);

  // perform the Cauchy endgame - returns the samples so that we can use them for rank determination
  retVal = CauchyEG_d2(pathNum, Final, last_approx, &endSamples, &cycle, &samples_per_loop, finalTime, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_d(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &endSamples[0], cycle, samples_per_loop, T, OUT, ED, eval_func);
  
    // clear endSamples
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_d(&endSamples[i]);
    free(endSamples);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;
  }

  return retVal;
}

int Cauchy_rank_main_d(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, point_data_d *Final, point_d last_approx, comp_d finalTime, point_data_d *samplePt, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: uses Cauchy endgame to determine if rank deficient     *
\***************************************************************/
{ // error checking - make sure that we actually have some information
  if (Cauchy_retVal != 0 && Cauchy_retVal != retVal_EG_failed_to_converge)
    return Cauchy_retVal;

  int retVal = 0;
  eval_struct_d e_d;
  init_eval_struct_d(e_d, 0, 0, 0);

  // do a quick non-singular check
  if (Cauchy_nonsing_check_d(condNum, smallest_nonzero_SV, largest_zero_SV, Final, cycle_num, T, OUT, &e_d, ED, eval_func))
  { // definitely non-singular
    if (rankType == 0)
    { // print condition number and rankDef
      *rankDef = 0;
      fprintf(OUT, "CN: %e\nrankDef: %d\n", *condNum, *rankDef);
    }
    else
    { // print condition number and corank
      *corank = 0;
      fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
    }

    retVal = 0;
    clear_eval_struct_d(e_d);
    return retVal;
  }
  
  // perform the main rank code
  int i, useSecondOption = 0, curr_cycle_size = cycle_num;
  double closed_loop_tol, min_closed_loop_tol = 1e-12, max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  comp_d endTime;
  point_d approx_double;
  // initially allocate for Samples - want to use 2 * samples on loop
  point_data_d *Samples = (point_data_d *)bmalloc((2 * samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));

  // initialize
  init_point_d(approx_double, 0);
  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
    init_point_data_d(&Samples[i], 0);

  // setup the time of the next sample point
  mul_rdouble_d(endTime, samplePt->time, (1 + T->power_series_sample_factor) / 2.0);

  // track to the next sample point
  retVal = track_d(&Samples[0], samplePt, endTime, T, OUT, ED, eval_func);

  // check for errors
  if (!retVal)
  { // find the tolerance to close the loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &Samples[0], NULL, 52, T, &e_d, NULL, ED, NULL, eval_func, NULL, NULL);

    // find the set of samples
    retVal = find_Cauchy_samples_d(closed_loop_tol, &Samples, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, &Samples[0], T, OUT, &e_d, ED, eval_func);

    // check for errors
    if (!retVal)
    { // find the 'double' approximation
      find_Cauchy_approx_d(approx_double, Samples, cycle_num, 2 * samples_per_loop);

      // use the last approximation and the new 'double' approximation to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_d(condNum, last_approx, approx_double, finalTime, T, OUT, &e_d, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_d(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, approx_double, finalTime, T, OUT, &e_d, ED, eval_func); 
      }

      // see if we need to check for convergence
      if (Cauchy_retVal == 0)
      { // just update Final
        point_cp_d(Final->point, approx_double);
        set_d(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        retVal = 0;
      }
      else
      { // check for convergence
        vec_d tempVec;
        init_vec_d(tempVec, 0); 

        vec_sub_d(tempVec, approx_double, Final->point); // compare the last 2 approximations
        T->error_at_latest_sample_point = infNormVec_d(tempVec);
        T->t_val_at_latest_sample_point = samplePt[0].time->r;
  
        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        point_cp_d(Final->point, approx_double);
        set_d(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        clear_vec_d(tempVec);
      }
    }
    else
    { // need to try something different
      useSecondOption = 1;
    }
  }
  else
  { // need to try something different
    useSecondOption = 1;
  }

  if (useSecondOption)
  { // sample more points around the current loop

    // find the tolerance to close the next loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, samplePt, NULL, 52, T, &e_d, NULL, ED, NULL, eval_func, NULL, NULL);

    // find the set of samples
    retVal = find_Cauchy_samples_d(closed_loop_tol, &Samples, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, samplePt, T, OUT, &e_d, ED, eval_func);

    // check for errors
    if (!retVal)
    { // find the 'double' approximation
      find_Cauchy_approx_d(approx_double, Samples, cycle_num, 2 * samples_per_loop);

      // use the last approximation and the new 'double' approximation to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_d(condNum, last_approx, approx_double, finalTime, T, OUT, &e_d, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_d(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, approx_double, finalTime, T, OUT, &e_d, ED, eval_func);
      }

      // see if we need to check for convergence
      if (Cauchy_retVal == 0)
      { // just update Final
        point_cp_d(Final->point, approx_double);
        set_d(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        retVal = 0;
      }
      else
      { // check for convergence
        vec_d tempVec;
        init_vec_d(tempVec, 0);

        vec_sub_d(tempVec, approx_double, Final->point); // compare the last 2 approximations
        T->error_at_latest_sample_point = infNormVec_d(tempVec);
        T->t_val_at_latest_sample_point = samplePt[0].time->r;

        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        point_cp_d(Final->point, approx_double);
        set_d(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        clear_vec_d(tempVec);
      }
    }
    else
    { // use the last approximation and Final to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'Final'
        *rankDef = Cauchy_rankDef_d(condNum, last_approx, Final->point, finalTime, T, OUT, &e_d, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'Final'
        *corank = Cauchy_corank_d(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, Final->point, finalTime, T, OUT, &e_d, ED, eval_func);
      }

      // we report failure
      retVal = retVal_sharpening_failed;
    }
  }

  // clear memory
  clear_point_d(approx_double);
  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_d(&Samples[i]);
  free(Samples);
  clear_eval_struct_d(e_d);

  return retVal;
}

int Cauchy_rankDef_d(double *CN, point_d approx0, point_d approx1, comp_d finalTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if this endpoint*
* is rank deficient or not                                      *
\***************************************************************/
{
  int rankDef = 0;
  double mat_norm, max_CN = 1e13, max_SV_ratio = T->ratioTol, SV_tol = MAX(T->sing_val_zero_tol, 1e-15), minSV[2] = {0, 0};

  // find the minimum singular values of the Jacobian at the two points
  eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, approx0, finalTime, ED, eval_func);
  approx_min_svd_d(&minSV[0], e_d->Jv);

  eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, approx1, finalTime, ED, eval_func);
  approx_min_svd_d(&minSV[1], e_d->Jv);

  // find the norm of e_d->Jv
  mat_norm = infNormMat_d(e_d->Jv);

  // determine if we have rank deficiency or not
  rankDef = rankDef_d(CN, minSV[0], minSV[1], mat_norm, max_CN, max_SV_ratio, SV_tol);

  // print CN & rankDef
  fprintf(OUT, "CN: %e\nrankDef: %d\n", *CN, rankDef);

  return rankDef;
}

int Cauchy_corank_d(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_d approx0, point_d approx1, comp_d finalTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the endpoint                                                  *
\***************************************************************/
{
  int corank = 0;
  double max_CN = 1e13, max_SV_ratio = T->ratioTol, SV_tol = MAX(T->sing_val_zero_tol, 1e-15);
  mat_d tempMat;

  init_mat_d(tempMat, 0, 0);

  // find the Jacobian at the two points
  eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, tempMat, e_d->Jp, approx0, finalTime, ED, eval_func);
  eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, approx1, finalTime, ED, eval_func);

  // determine the corank
  corank = corank_d(CN, smallest_nonzero_SV, largest_zero_SV, tempMat, e_d->Jv, max_CN, max_SV_ratio, SV_tol);

  // print CN & corank
  fprintf(OUT, "CN: %e\ncorank: %d\n", *CN, corank);

  clear_mat_d(tempMat);

  return corank;
}

int Cauchy_nonsing_check_d(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - non-singular, 0 - otherwise                *
* NOTES: does a quick non-singular check                        *
\***************************************************************/
{ // check cycle_num == 1
  if (cycle_num != 1)
    return 0;

  // so we have cycle_num == 1 and so we need to check minimum singular value
  int retVal;
  double norm, minSV, minSV_acceptable = MIN_SV_CHECK;

  // evaluate function
  eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, Final->point, Final->time, ED, eval_func);

  // approximate minSV & norm
  approx_min_svd_d(&minSV, e_d->Jv);
  norm = infNormMat_d(e_d->Jv);

  // check to normalize and find CN
  if (norm > 1)
  { // normalize
    minSV /= norm;
    // find CN
    *CN = 1 / minSV;
  }
  else
  { // find CN
    *CN = norm / minSV;
  }

  // see if it is large enough
  if (minSV > minSV_acceptable)
    retVal = 1;
  else
    retVal = 0;

  // do a final refinement if this is non-singular
  if (retVal == 1)
  {
    int refine_retVal = 0;

    // refine to full precision
    refine_retVal = refine_d_basic(Final, T, OUT, e_d, ED, eval_func);

    // make sure refine was successful
    if (refine_retVal)
    { // failure - this should not happen!
      retVal = 0;
      *smallest_nonzero_SV = *largest_zero_SV = 0;
    }
    else
    { // setup smallest_nonzero_SV & largest_zero_SV
      *smallest_nonzero_SV = minSV;
      *largest_zero_SV = 0;
    }
  }

  return retVal;
}

///////// Multi Precision /////////

int CauchyEG_rank_mp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in fixed precision along  *
* with rank determination                                       *
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

  // perform the Cauchy endgame - returns the samples so that we can use them for rank determination
  retVal = CauchyEG_mp2(pathNum, Final, last_approx, &endSamples, &cycle, &samples_per_loop, finalTime, Start, T, OUT, midOUT, ED, eval_func, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    retVal = Cauchy_rank_main_mp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, Final, last_approx, finalTime, &endSamples[0], cycle, samples_per_loop, T, OUT, ED, eval_func);

    // clear endSamples
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_mp(&endSamples[i]);
    free(endSamples);
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;
  }

  clear_mp(finalTime);

  return retVal;
}

int Cauchy_rank_main_mp(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, point_data_mp *Final, point_mp last_approx, comp_mp finalTime, point_data_mp *samplePt, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: uses Cauchy endgame to determine if rank deficient     *
\***************************************************************/
{ // error checking - make sure that we actually have some information
  if (Cauchy_retVal != 0 && Cauchy_retVal != retVal_EG_failed_to_converge)
    return Cauchy_retVal;

  int retVal = 0;
  eval_struct_mp e_mp;
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // do a quick non-singular check
  if (Cauchy_nonsing_check_mp(condNum, smallest_nonzero_SV, largest_zero_SV, Final, cycle_num, T, OUT, &e_mp, ED, eval_func))
  { // definitely non-singular
    if (rankType == 0)
    { // print condition number and rankDef
      *rankDef = 0;
      fprintf(OUT, "CN: %e\nrankDef: %d\n", *condNum, *rankDef);
    }
    else
    { // print condition number and corank
      *corank = 0;
      fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
    }

    clear_eval_struct_mp(e_mp);
    retVal = 0;
    return retVal;
  }

  // perform the main rank code
  int i, useSecondOption = 0, curr_cycle_size = cycle_num;
  double closed_loop_tol, min_closed_loop_tol = MAX(1e-150, pow(10, -(prec_to_digits(T->Precision) - 5))), max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  comp_mp endTime;
  point_mp approx_double;
  point_data_mp *Samples = NULL;

  // initialize MP
  init_mp(endTime);
  init_point_mp(approx_double, 0);

  // initially allocate for Samples - want to use 2 * samples on next loop
  Samples = (point_data_mp *)bmalloc((2 * samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
    init_point_data_mp(&Samples[i], 0);

  // setup the time of the next sample point
  mul_rdouble_mp(endTime, samplePt->time, (1 + T->power_series_sample_factor) / 2.0);

  // track to the next sample point
  retVal = track_mp(&Samples[0], samplePt, endTime, T, OUT, ED, eval_func);

  // check for errors
  if (!retVal)
  { // find the tolerance to close the next loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, &Samples[0], T->Precision, T, NULL, &e_mp, NULL, ED, NULL, eval_func, NULL);

    // find the set of samples
    retVal = find_Cauchy_samples_mp(closed_loop_tol, &Samples, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, &Samples[0], T, OUT, &e_mp, ED, eval_func);

    // check for errors
    if (!retVal)
    { // find the 'double' approximation
      find_Cauchy_approx_mp(approx_double, Samples, cycle_num, 2 * samples_per_loop);

      // use the 2 approximations to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_mp(condNum, last_approx, approx_double, finalTime, T, OUT, &e_mp, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_mp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, approx_double, finalTime, T, OUT, &e_mp, ED, eval_func);
      }

      // see if we need to check for convergence
      if (Cauchy_retVal == 0)
      { // just update Final
        point_cp_mp(Final->point, approx_double);
        set_mp(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        retVal = 0;
      }
      else
      { // check for convergence
        mpf_t error;
        mpf_init(error);

        findDiff_point(error, NULL, approx_double, T->Precision, NULL, Final->point, T->Precision); // compare the last 2 approximations
        T->error_at_latest_sample_point = mpf_get_d(error);
        T->t_val_at_latest_sample_point = mpf_get_d(endTime->r);

        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        point_cp_mp(Final->point, approx_double);
        set_mp(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        mpf_clear(error);
      }
    }
    else
    { // need to try something different
      useSecondOption = 1;
    }
  }
  else
  { // need to try something different
    useSecondOption = 1;
  }

  if (useSecondOption)
  { // sample more points around the current loop

    // find the tolerance to close the next loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, NULL, samplePt, T->Precision, T, NULL, &e_mp, NULL, ED, NULL, eval_func, NULL);

    // find the set of samples
    retVal = find_Cauchy_samples_mp(closed_loop_tol, &Samples, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, samplePt, T, OUT, &e_mp, ED, eval_func);

    // check for errors
    if (!retVal)
    { // find the 'double' approximation
      find_Cauchy_approx_mp(approx_double, Samples, cycle_num, 2 * samples_per_loop);

      // use the 2 approximations to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_mp(condNum, last_approx, approx_double, finalTime, T, OUT, &e_mp, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_mp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, approx_double, finalTime, T, OUT, &e_mp, ED, eval_func);
      }

      // see if we need to check for convergence
      if (Cauchy_retVal == 0)
      { // just update Final
        point_cp_mp(Final->point, approx_double);
        set_mp(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        retVal = 0;
      }
      else
      { // check for convergence
        mpf_t error;
        mpf_init(error);

        findDiff_point(error, NULL, approx_double, T->Precision, NULL, Final->point, T->Precision); // compare the last 2 approximations
        T->error_at_latest_sample_point = mpf_get_d(error);
        T->t_val_at_latest_sample_point = mpf_get_d(endTime->r);

        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        point_cp_mp(Final->point, approx_double);
        set_mp(Final->time, finalTime);
        Final->cycle_num = cycle_num;

        mpf_clear(error);
      }
    }
    else
    { // use the last approximation and Final to determine either rank deficiency or corank
      if (rankType == 0)
      { // determine rank deficiency by looking at 'last_approx' & 'Final'
        *rankDef = Cauchy_rankDef_mp(condNum, last_approx, Final->point, finalTime, T, OUT, &e_mp, ED, eval_func);
      }
      else
      { // determine corank by looking at 'last_approx' & 'Final'
        *corank = Cauchy_corank_mp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx, Final->point, finalTime, T, OUT, &e_mp, ED, eval_func);
      }

      // we report failure
      retVal = retVal_sharpening_failed;
    }
  }

  // clear MP
  clear_mp(endTime);
  clear_point_mp(approx_double);
  clear_eval_struct_mp(e_mp);
  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
    clear_point_data_mp(&Samples[i]);
  free(Samples);

  return retVal;
}

int Cauchy_rankDef_mp(double *CN, point_mp approx0, point_mp approx1, comp_mp finalTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if this endpoint*
* is rank deficient or not                                      *
\***************************************************************/
{
  int rankDef = 0, curr_prec = T->Precision, prec_digits = prec_to_digits(T->Precision) - 4;
  double max_CN = MIN(1e150, pow(10, prec_digits)), max_SV_ratio = T->ratioTol, SV_tol = MAX(T->sing_val_zero_tol, pow(10, -prec_digits - 2));
  mpf_t mat_norm, minSV[2];

  mpf_init(mat_norm); mpf_init(minSV[0]); mpf_init(minSV[1]);

  // find the minimum singular values of the Jacobian at the two points
  eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, approx0, finalTime, ED, eval_func);
  approx_min_svd_mp_prec(minSV[0], e_mp->Jv, curr_prec);

  eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, approx1, finalTime, ED, eval_func);
  approx_min_svd_mp_prec(minSV[1], e_mp->Jv, curr_prec);

  // find the norm of e_mp->Jv
  infNormMat_mp2(mat_norm, e_mp->Jv);

  // determine if we have rank deficiency or not
  rankDef = rankDef_mp(CN, minSV[0], minSV[1], mat_norm, curr_prec, max_CN, max_SV_ratio, SV_tol);

  // print CN & rankDef
  fprintf(OUT, "CN: %e\nrankDef: %d\n", *CN, rankDef);

  mpf_clear(mat_norm); mpf_clear(minSV[0]); mpf_clear(minSV[1]);

  return rankDef;
}

int Cauchy_corank_mp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_mp approx0, point_mp approx1, comp_mp finalTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the endpoint                                                  *
\***************************************************************/
{
  int corank = 0, curr_prec = T->Precision, prec_digits = prec_to_digits(T->Precision) - 4;
  double max_CN = MIN(1e150, pow(10, prec_digits)), max_SV_ratio = T->ratioTol, SV_tol = MAX(T->sing_val_zero_tol, pow(10, -prec_digits - 2));
  mat_mp tempMat;

  init_mat_mp(tempMat, 0, 0);

  // find the Jacobian at the two points
  eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, tempMat, e_mp->Jp, approx0, finalTime, ED, eval_func);
  eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, approx1, finalTime, ED, eval_func);

  // determine the corank
  corank = corank_mp(CN, smallest_nonzero_SV, largest_zero_SV, tempMat, e_mp->Jv, curr_prec, max_CN, max_SV_ratio, SV_tol);

  // print CN & corank
  fprintf(OUT, "CN: %e\ncorank: %d\n", *CN, corank);

  clear_mat_mp(tempMat);

  return corank;
}

int Cauchy_nonsing_check_mp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - non-singular, 0 - otherwise                *
* NOTES: does a quick non-singular check                        *
\***************************************************************/
{ // check cycle_num == 1
  if (cycle_num != 1)
    return 0;

  // so we have cycle_num == 1 and so we need to check minimum singular value
  int retVal;
  double minSV_acceptable = MIN_SV_CHECK;
  mpf_t minSV, norm;

  mpf_init(minSV);
  mpf_init(norm);

  // evaluate function
  eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, Final->point, Final->time, ED, eval_func);

  // approximate minSV & norm
  approx_min_svd_mp_prec(minSV, e_mp->Jv, T->Precision);
  infNormMat_mp2(norm, e_mp->Jv);

  // check to normalize and find CN
  if (mpf_cmp_ui(norm, 1) > 0)
  { // normalize
    mpf_div(minSV, minSV, norm);
    // find CN
    mpf_ui_div(norm, 1, minSV);
    *CN = mpf_get_d(norm);
  }
  else
  { // find CN
    mpf_div(norm, norm, minSV);
    *CN = mpf_get_d(norm);
  }

  // see if it is large enough
  if (mpf_cmp_d(minSV, minSV_acceptable) > 0)
    retVal = 1;
  else
    retVal = 0;

  // do a final refinement if this is non-singular
  if (retVal == 1)
  {
    int refine_retVal = 0;

    // refine to full precision
    refine_retVal = refine_mp_basic(Final, T, OUT, e_mp, ED, eval_func);

    // make sure refine was successful
    if (refine_retVal)
    { // failure - this should not happen!
      retVal = 0;
      *smallest_nonzero_SV = *largest_zero_SV = 0;
    }
    else
    { // setup smallest_nonzero_SV & largest_zero_SV
      *smallest_nonzero_SV = mpf_get_d(minSV);
      *largest_zero_SV = 0;
    }
  }
  
  // clear
  mpf_clear(minSV);
  mpf_clear(norm);

  return retVal;
}

//////// Adaptive Multi Precision /////////

int CauchyEG_rank_amp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does Cauchy integral endgame in adaptive precision     *
* along with rank determination                                 *
\***************************************************************/
{
  int i, retVal, cycle, samples_per_loop, samples_prec;
  comp_d finalTime_d;
  comp_mp finalTime_mp;
  point_data_d *endSamples_d = NULL;
  point_data_mp *endSamples_mp = NULL;

  init_mp(finalTime_mp);

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
    point_cp_mp(last_approx_mp, Start_mp->point);
  }

  // perform the Cauchy endgame - returns the samples so that we can use them for rank determination
  retVal = CauchyEG_amp2(pathNum, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, &endSamples_d, &endSamples_mp, &cycle, &samples_per_loop, finalTime_d, finalTime_mp, Start_d, Start_mp, prec_in, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

  if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination

    // save the precision on which endSamples is setup with
    samples_prec = *prec;

    retVal = Cauchy_rank_main_amp(condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, retVal, prec, time_first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, *last_approx_prec, finalTime_d, finalTime_mp, &endSamples_d[0], &endSamples_mp[0], samples_prec, cycle, samples_per_loop, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // clear endSamples
    if (samples_prec < 64)
    { // clear endSamples_d
      for (i = samples_per_loop * cycle - 1; i >= 0; i--)
        clear_point_data_d(&endSamples_d[i]);
      free(endSamples_d);
    }
    else
    { // clear endSamples_mp
      for (i = samples_per_loop * cycle - 1; i >= 0; i--)
        clear_point_data_mp(&endSamples_mp[i]);
      free(endSamples_mp);
    }
  }
  else
  { // setup condNum since we have failure
    *condNum = -1;
    *smallest_nonzero_SV = *largest_zero_SV = 0;
  }

  clear_mp(finalTime_mp);

  return retVal;
}

int Cauchy_rank_main_amp(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *samplePt_d, point_data_mp *samplePt_mp, int prec_in, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: uses Cauchy endgame to determine if rank deficient     *
\***************************************************************/
{ // since this will leave the precision of Final alone, set prec
  *prec = prec_in;

  // error checking - make sure that we actually have some information
  if (Cauchy_retVal != 0 && Cauchy_retVal != retVal_EG_failed_to_converge)
    return Cauchy_retVal;

  int retVal = 0;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // do a quick non-singular check
  if (Cauchy_nonsing_check_amp(condNum, smallest_nonzero_SV, largest_zero_SV, prec, Final_d, Final_mp, cycle_num, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec))
  { // definitely non-singular
    if (rankType == 0)
    { // print condition number and rankDef
      *rankDef = 0;
      fprintf(OUT, "CN: %e\nrankDef: %d\n", *condNum, *rankDef);
    }
    else
    { // print condition number and corank
      *corank = 0;
      fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
    }

    clear_eval_struct_d(e_d);     
    clear_eval_struct_mp(e_mp);
    retVal = 0;
    return retVal;
  }
 
  // perform the main rank code 
  int i, useSecondOption = 0, curr_cycle_size = cycle_num, curr_prec = MAX(prec_in, 64), digits_in = prec_to_digits(prec_in), orig_sv_prec, new_digits;
  double closed_loop_tol, min_closed_loop_tol = 1e-150, max_closed_loop_tol = MAX(T->currentNewtonTol, T->final_tol_times_mult);
  mpf_t error, orig_min_sv, new_min_sv, new_mat_norm;
  comp_d endTime_d;
  comp_mp endTime_mp;
  point_d approx_double_d;
  point_mp approx_double_mp;
  point_data_d *Samples_d = NULL;
  point_data_mp *Samples_mp = NULL;
  int *Samples_prec = NULL;

  // initialize to curr_prec
  initMP(curr_prec);
  change_prec(ED_mp, curr_prec);

  // initialize
  mpf_init(error); mpf_init(orig_min_sv); mpf_init(new_min_sv); mpf_init(new_mat_norm);
  init_mp(endTime_mp);
  init_point_d(approx_double_d, 0);
  init_point_mp(approx_double_mp, 0);
  setprec_eval_struct_mp(e_mp, curr_prec);

  // initially allocate for Samples - want to use 2 * samples on next loop
  Samples_d = (point_data_d *)bmalloc((2 * samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_d));
  Samples_mp = (point_data_mp *)bmalloc((2 * samples_per_loop * curr_cycle_size + 1) * sizeof(point_data_mp));
  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
  {
    init_point_data_d(&Samples_d[i], 0);
    init_point_data_mp(&Samples_mp[i], 0);
  }

  // find the end time for the next sample
  if (prec_in < 64)
  { // setup endTime_d
    mul_rdouble_d(endTime_d, samplePt_d->time, (1 + T->power_series_sample_factor) / 2.0);
    d_to_mp(endTime_mp, endTime_d);
  }
  else
  { // setup endTime_mp
    mul_rdouble_mp(endTime_mp, samplePt_mp->time, (1 + T->power_series_sample_factor) / 2.0);
    mp_to_d(endTime_d, endTime_mp);
  }

  // track to the next sample point
  retVal = AMPtrack(&Samples_d[0], &Samples_mp[0], &curr_prec, time_first_increase, samplePt_d, samplePt_mp, prec_in, endTime_d, endTime_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // check for errors
  if (!retVal)
  { // refine the sample point
    if (curr_prec < 64)
    { // refine using _d
      refine_d_basic(&Samples_d[0], T, OUT, &e_d, ED_d, eval_func_d);
    }
    else
    { // refine using _mp
      initMP(curr_prec);
      change_prec(ED_mp, curr_prec);
      setprec_eval_struct_mp(e_mp, curr_prec);

      refine_mp_basic(&Samples_mp[0], T, OUT, &e_mp, ED_mp, eval_func_mp);
    }

    // find the tolerance to close the next loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, &Samples_d[0], &Samples_mp[0], curr_prec, T, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // find the set of samples
    retVal = find_Cauchy_samples_amp(closed_loop_tol, &Samples_d, &Samples_mp, &curr_prec, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, &Samples_d[0], &Samples_mp[0], curr_prec, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for errors
    if (!retVal)
    { // setup Samples_prec - do not care about the last sample since it is equal to the first
      Samples_prec = (int *)bmalloc((2 * samples_per_loop * cycle_num) * sizeof(int));
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 0; i--)
        Samples_prec[i] = curr_prec;

      // find the precision needed to see the minimum singular of the last approximation as non-zero
      new_digits = MAX(64, last_approx_prec);
      setprec_point_data_mp(&Samples_mp[2 * samples_per_loop * cycle_num], new_digits);

      d_to_mp(Samples_mp[2 * samples_per_loop * cycle_num].time, finalTime_d);
      Samples_mp[2 * samples_per_loop * cycle_num].cycle_num = cycle_num;
      if (last_approx_prec < 64)
      {
        point_d_to_mp(Samples_mp[2 * samples_per_loop * cycle_num].point, last_approx_d);
      }
      else
      {
        point_cp_mp(Samples_mp[2 * samples_per_loop * cycle_num].point, last_approx_mp);
      }

      // finds the original minium singular value and the precision to see it as non-zero
      orig_sv_prec = prec_needed(T->MPType, orig_min_sv, error, &Samples_mp[2 * samples_per_loop * cycle_num], new_digits, &e_mp, ED_mp, eval_func_mp, change_prec);

      // find the number of digits to use - add 10 digits on to the precision needed
      new_digits = prec_to_digits(orig_sv_prec) + 10;

      // refine the samples that are needed to create the approximation to this number of digits
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 0; i--)
      { // setup endTime_d
        if (Samples_prec[i] < 64)
        {
          set_d(endTime_d, Samples_d[i].time);
        }
        else
        {
          mp_to_d(endTime_d, Samples_mp[i].time);
        }
      
        refine_digits_amp(T->outputLevel, new_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, digits_in, &Samples_d[i], &Samples_mp[i], &Samples_prec[i], &Samples_d[i], &Samples_mp[i], Samples_prec[i], endTime_d, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      }

      // all should have the same precision - let's find it
      curr_prec = Samples_prec[0];
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 1; i--)
        if (curr_prec > Samples_prec[i])
          curr_prec = Samples_prec[i];

      // find the 'double' approximation
      find_Cauchy_approx_amp(approx_double_d, approx_double_mp, &curr_prec, Samples_d, Samples_mp, curr_prec, cycle_num, 2 * samples_per_loop);

      // use the 2 approximations to determine either rank deficiency or corank
      if (rankType == 0)
      { // find the minimum singular value for the 'double' approximation
        initMP(curr_prec);
        change_prec(ED_mp, curr_prec);
        setprec_eval_struct_mp(e_mp, curr_prec);
        mpf_set_prec(new_min_sv, curr_prec);
        mpf_set_prec(new_mat_norm, curr_prec);
        setprec_mp(finalTime_mp, curr_prec);
        d_to_mp(finalTime_mp, finalTime_d);

        // evaluate the function
        eval_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, approx_double_mp, finalTime_mp, ED_mp, eval_func_mp);
  
        // find the minimum singular value and matrix norm
        approx_min_svd_mp_prec(new_min_sv, e_mp.Jv, curr_prec);
        infNormMat_mp2(new_mat_norm, e_mp.Jv);

        // determine rank deficiency by looking at the minimum singular values of 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_amp(condNum, orig_min_sv, T->ratioTol, orig_sv_prec, new_min_sv, curr_prec, new_mat_norm, OUT);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_amp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx_d, last_approx_mp, last_approx_prec, approx_double_d, approx_double_mp, curr_prec, finalTime_d, finalTime_mp, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // error checking
        if (cycle_num > 1 && *corank == 0)
          *corank = 1;

        // print CN & corank
        fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
      }

      // see if we need to check for convergence & update Final - do not change the precision on Final since that was good enough before
      if (Cauchy_retVal == 0)
      { // just update Final
        if (prec_in < 64)
        { // update Final_d
          if (curr_prec < 64)
          { // use approx_double_d
            point_cp_d(Final_d->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_mp_to_d(Final_d->point, approx_double_mp);
          }
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = cycle_num; 
        }
        else
        { // update Final_mp
          if (curr_prec < 64)
          { // use approx_double_d - should not happen
            point_d_to_mp(Final_mp->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_cp_mp(Final_mp->point, approx_double_mp);
          }
          set_mp(Final_mp->time, finalTime_mp);
          Final_mp->cycle_num = cycle_num;
        }

        retVal = 0;
      }
      else
      { // check for convergence
        findDiff_point(error, Final_d->point, Final_mp->point, prec_in, approx_double_d, approx_double_mp, curr_prec); // compare the last 2 approximations

        T->error_at_latest_sample_point = mpf_get_d(error);

        if (prec_in < 64)
          T->t_val_at_latest_sample_point = endTime_d->r;
        else
          T->t_val_at_latest_sample_point = mpf_get_d(endTime_mp->r);

        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        if (prec_in < 64)
        { // update Final_d
          if (curr_prec < 64)
          { // use approx_double_d
            point_cp_d(Final_d->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_mp_to_d(Final_d->point, approx_double_mp);
          }
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = cycle_num;
        }
        else
        { // update Final_mp
          if (curr_prec < 64)
          { // use approx_double_d - should not happen
            point_d_to_mp(Final_mp->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_cp_mp(Final_mp->point, approx_double_mp);
          }
          set_mp(Final_mp->time, finalTime_mp);
          Final_mp->cycle_num = cycle_num;
        }
      }
    }
    else
    { // need to try something different
      useSecondOption = 1;
    }
  }
  else
  { // need to try something different
    useSecondOption = 1;
  }

  if (useSecondOption) // since using AMP, this should only happen in drastic cases!
  { // sample more points around the current loop

    // find the tolerance to close the next loop
    closed_loop_tol = find_closed_loop_tol(min_closed_loop_tol, max_closed_loop_tol, samplePt_d, samplePt_mp, prec_in, T, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // find the set of samples
    retVal = find_Cauchy_samples_amp(closed_loop_tol, &Samples_d, &Samples_mp, &curr_prec, &cycle_num, 2 * samples_per_loop, &curr_cycle_size, samplePt_d, samplePt_mp, prec_in, time_first_increase, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check for errors
    if (!retVal)
    { // setup Samples_prec - do not care about the last sample since it is equal to the first
      Samples_prec = (int *)bmalloc((2 * samples_per_loop * cycle_num) * sizeof(int));
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 0; i--)
        Samples_prec[i] = curr_prec;

      // find the precision needed to see the minimum singular value of the current approximation as non-zero
      new_digits = MAX(64, last_approx_prec);
      setprec_point_data_mp(&Samples_mp[2 * samples_per_loop * cycle_num], new_digits);

      d_to_mp(Samples_mp[2 * samples_per_loop * cycle_num].time, finalTime_d);
      Samples_mp[2 * samples_per_loop * cycle_num].cycle_num = cycle_num;
      if (last_approx_prec < 64)
      {
        point_d_to_mp(Samples_mp[2 * samples_per_loop * cycle_num].point, last_approx_d);
      }
      else
      {
        point_cp_mp(Samples_mp[2 * samples_per_loop * cycle_num].point, last_approx_mp);
      }

      // finds the original minium singular value and the precision to see it as non-zero
      orig_sv_prec = prec_needed(T->MPType, orig_min_sv, error, &Samples_mp[2 * samples_per_loop * cycle_num], new_digits, &e_mp, ED_mp, eval_func_mp, change_prec);

      // find the number of digits to use - add 10 digits on to the precision needed
      new_digits = prec_to_digits(orig_sv_prec) + 10;

      // refine the samples that are needed to create the approximation to this number of digits
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 0; i--)
      { // setup endTime_d
        if (Samples_prec[i] < 64)
        {
          set_d(endTime_d, Samples_d[i].time);
        }
        else
        {
          mp_to_d(endTime_d, Samples_mp[i].time);
        }

        refine_digits_amp(T->outputLevel, new_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, digits_in, &Samples_d[i], &Samples_mp[i], &Samples_prec[i], &Samples_d[i], &Samples_mp[i], Samples_prec[i], endTime_d, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      }

      // all should have the same precision - let's find it
      curr_prec = Samples_prec[0];
      for (i = 2 * samples_per_loop * cycle_num - 1; i >= 1; i--)
        if (curr_prec > Samples_prec[i])
          curr_prec = Samples_prec[i];

      // find the 'double' approximation
      find_Cauchy_approx_amp(approx_double_d, approx_double_mp, &curr_prec, Samples_d, Samples_mp, curr_prec, cycle_num, 2 * samples_per_loop);

      // use the 2 approximations to determine either rank deficiency or corank
      if (rankType == 0)
      { // find the minimum singular values for the 'last approximation' and 'Final'
        initMP(curr_prec);
        change_prec(ED_mp, curr_prec);
        setprec_eval_struct_mp(e_mp, curr_prec);
        mpf_set_prec(new_min_sv, curr_prec);
        mpf_set_prec(new_mat_norm, curr_prec);
        setprec_mp(finalTime_mp, curr_prec);
        d_to_mp(finalTime_mp, finalTime_d);

        // evaluate the function
        eval_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, approx_double_mp, finalTime_mp, ED_mp, eval_func_mp);

        // find the minimum singular value and matrix norm
        approx_min_svd_mp_prec(new_min_sv, e_mp.Jv, curr_prec);
        infNormMat_mp2(new_mat_norm, e_mp.Jv);

        // determine rank deficiency by looking at the minimum singular values of 'last_approx' & 'approx_double'
        *rankDef = Cauchy_rankDef_amp(condNum, orig_min_sv, T->ratioTol, orig_sv_prec, new_min_sv, curr_prec, new_mat_norm, OUT);
      }
      else
      { // determine corank by looking at 'last_approx' & 'approx_double'
        *corank = Cauchy_corank_amp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx_d, last_approx_mp, last_approx_prec, approx_double_d, approx_double_mp, curr_prec, finalTime_d, finalTime_mp, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // error checking
        if (cycle_num > 1 && *corank == 0)
          *corank = 1;

        // print CN & corank
        fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
      }

      // see if we need to check for convergence & update Final - do not change the precision on Final since that was good enough before
      if (Cauchy_retVal == 0)
      { // just update Final
        if (prec_in < 64)
        { // update Final_d
          if (curr_prec < 64)
          { // use approx_double_d
            point_cp_d(Final_d->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_mp_to_d(Final_d->point, approx_double_mp);
          }
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = cycle_num;
        }
        else
        { // update Final_mp
          if (curr_prec < 64)
          { // use approx_double_d - should not happen
            point_d_to_mp(Final_mp->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_cp_mp(Final_mp->point, approx_double_mp);
          }
          set_mp(Final_mp->time, finalTime_mp);
          Final_mp->cycle_num = cycle_num;
        }

        retVal = 0;
      }
      else
      { // check for convergence
        findDiff_point(error, Final_d->point, Final_mp->point, prec_in, approx_double_d, approx_double_mp, curr_prec); // compare the last 2 approximations

        T->error_at_latest_sample_point = mpf_get_d(error);

        if (prec_in < 64)
          T->t_val_at_latest_sample_point = endTime_d->r;
        else
          T->t_val_at_latest_sample_point = mpf_get_d(endTime_mp->r);

        if (T->error_at_latest_sample_point < T->final_tolerance)
        { // we have convergence
          retVal = 0;
        }
        else
        { // still no convergence
          retVal = retVal_EG_failed_to_converge;
        }

        // update Final
        if (prec_in < 64)
        { // update Final_d
          if (curr_prec < 64)
          { // use approx_double_d
            point_cp_d(Final_d->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_mp_to_d(Final_d->point, approx_double_mp);
          }
          set_d(Final_d->time, finalTime_d);
          Final_d->cycle_num = cycle_num;
        }
        else
        { // update Final_mp
          if (curr_prec < 64)
          { // use approx_double_d - should not happen
            point_d_to_mp(Final_mp->point, approx_double_d);
          }
          else
          { // use approx_double_mp
            point_cp_mp(Final_mp->point, approx_double_mp);
          }
          set_mp(Final_mp->time, finalTime_mp);
          Final_mp->cycle_num = cycle_num;
        }
      }
    }
    else
    { // use last_approx & Final to compute
      if (rankType == 0)
      {
        int corank = Cauchy_corank_amp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx_d, last_approx_mp, last_approx_prec, Final_d->point, Final_mp->point, prec_in, finalTime_d, finalTime_mp, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        *rankDef = corank > 0;

        // print CN & rankDef
        fprintf(OUT, "CN: %e\nrankDef: %d\n", *condNum, *rankDef);
      }
      else
      {
        *corank = Cauchy_corank_amp(condNum, smallest_nonzero_SV, largest_zero_SV, last_approx_d, last_approx_mp, last_approx_prec, Final_d->point, Final_mp->point, prec_in, finalTime_d, finalTime_mp, T, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // print CN & corank
        fprintf(OUT, "CN: %e\ncorank: %d\n", *condNum, *corank);
      }

      // we report failure
      retVal = retVal_sharpening_failed;
    }
  }

  // clear memory
  mpf_clear(error); mpf_clear(orig_min_sv); mpf_clear(new_min_sv); mpf_clear(new_mat_norm);
  clear_mp(endTime_mp);
  clear_point_d(approx_double_d);
  clear_point_mp(approx_double_mp);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  for (i = 2 * samples_per_loop * curr_cycle_size; i >= 0; i--)
  {
    clear_point_data_d(&Samples_d[i]);
    clear_point_data_mp(&Samples_mp[i]);
  }
  free(Samples_d);
  free(Samples_mp);
  free(Samples_prec);

  return retVal;
}

int Cauchy_rankDef_amp(double *CN, mpf_t minSV0, int ratioTol, int prec0, mpf_t minSV1, int prec1, mpf_t mat_norm, FILE *OUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if this endpoint*
* is rank deficient or not                                      *
\***************************************************************/
{
  int rankDef = 0;

  // determine if we have rank deficiency or not
  rankDef = rankDef_amp(CN, minSV0, prec0, minSV1, prec1, mat_norm, ratioTol);

  // print CN & rankDef
  fprintf(OUT, "CN: %e\nrankDef: %d\n", *CN, rankDef);

  return rankDef;
}

int Cauchy_corank_amp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_d approx0_d, point_mp approx0_mp, int prec0, point_d approx1_d, point_mp approx1_mp, int prec1, comp_d finalTime_d, comp_mp finalTime_mp, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the endpoint                                                  *
\***************************************************************/
{
  int corank = 0;
  double max_SV_ratio = T->ratioTol; 
  mat_d Jv0_d;
  mat_mp Jv0_mp;

  init_mat_d(Jv0_d, 0, 0);
  init_mat_mp(Jv0_mp, 0, 0);

  // find the first Jacobian
  if (prec0 < 64)
  { // evaluate in double precision
    eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, Jv0_d, e_d->Jp, approx0_d, finalTime_d, ED_d, eval_func_d);
  }
  else
  { // set the precision
    initMP(prec0);
    change_prec(ED_mp, prec0);
    setprec_eval_struct_mp(*e_mp, prec0);
    change_prec_mat_mp(Jv0_mp, prec0);

    // evaluate
    eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, Jv0_mp, e_mp->Jp, approx0_mp, finalTime_mp, ED_mp, eval_func_mp);
  }

  // find the second Jacobian
  if (prec1 < 64)
  { // evaluate in double precision
    eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, approx1_d, finalTime_d, ED_d, eval_func_d);
  }
  else
  { // set the precision
    initMP(prec1);
    change_prec(ED_mp, prec1);
    setprec_eval_struct_mp(*e_mp, prec1);

    // evaluate
    eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, approx1_mp, finalTime_mp, ED_mp, eval_func_mp);
  }

  // determine the corank
  corank = corank_amp(CN, smallest_nonzero_SV, largest_zero_SV, Jv0_d, Jv0_mp, prec0, e_d->Jv, e_mp->Jv, prec1, max_SV_ratio);

  // clear
  clear_mat_d(Jv0_d);
  clear_mat_mp(Jv0_mp);

  return corank;
}

int Cauchy_nonsing_check_amp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, point_data_d *Final_d, point_data_mp *Final_mp, int cycle_num, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - non-singular, 0 - otherwise                *
* NOTES: does a quick non-singular check                        *
\***************************************************************/
{ // check cycle_num == 1
  if (cycle_num != 1)
    return 0;

  // so we have cycle_num == 1 and so we need to check minimum singular value
  int retVal;
  double minSV_acceptable = MIN_SV_CHECK;

  if (*prec < 64)
  { // use _d
    double norm_d, minSV_d;

    eval_d(e_d->funcVals, e_d->parVals, e_d->parDer, e_d->Jv, e_d->Jp, Final_d->point, Final_d->time, ED_d, eval_func_d);

    // approximate the minSV & norm
    approx_min_svd_d(&minSV_d, e_d->Jv);
    norm_d = infNormMat_d(e_d->Jv);

    // check to normalize and find CN
    if (norm_d > 1)
    { // normalize
      minSV_d /= norm_d;
      // find CN
      *CN = 1 / minSV_d;
    }
    else
    { // find CN
      *CN = norm_d / minSV_d;
    }

    // see if it is large enough
    if (minSV_d > minSV_acceptable)
    {
      retVal = 1;
      *smallest_nonzero_SV = minSV_d;
      *largest_zero_SV = 0;
    }
    else
      retVal = 0;
  }
  else
  { // use _mp
    mpf_t norm_mp, minSV_mp;

    initMP(*prec);
    change_prec(ED_mp, *prec);
    setprec_eval_struct_mp(*e_mp, *prec);
    mpf_init(norm_mp);
    mpf_init(minSV_mp);

    eval_mp(e_mp->funcVals, e_mp->parVals, e_mp->parDer, e_mp->Jv, e_mp->Jp, Final_mp->point, Final_mp->time, ED_mp, eval_func_mp);

    // approximate the minSV & norm
    approx_min_svd_mp_prec(minSV_mp, e_mp->Jv, *prec);
    infNormMat_mp2(norm_mp, e_mp->Jv);

    // check to normalize and find CN
    if (mpf_cmp_ui(norm_mp, 1) > 0)
    { // normalize
      mpf_div(minSV_mp, minSV_mp, norm_mp);
      // find CN
      mpf_ui_div(norm_mp, 1, minSV_mp);
      *CN = mpf_get_d(norm_mp);
    }
    else
    { // find CN
      mpf_div(norm_mp, norm_mp, minSV_mp);
      *CN = mpf_get_d(norm_mp);
    }

    // see if it is large enough
    if (mpf_cmp_d(minSV_mp, minSV_acceptable) > 0)
    {
      retVal = 1;
      *smallest_nonzero_SV = mpf_get_d(minSV_mp);
      *largest_zero_SV = 0;
    }
    else
      retVal = 0;
 
    // clear
    mpf_clear(norm_mp);
    mpf_clear(minSV_mp);
  }

  // do a final refinement if this is non-singular
  if (retVal == 1)
  { // increase the number of digits
    int refine_retVal = 0, prec_out, digits_new = prec_to_digits(*prec) + 10, digits_in = prec_to_digits(*prec);
    comp_d time;
    point_data_d tempPD_d;
    point_data_mp tempPD_mp;

    init_point_data_d(&tempPD_d, 0);
    init_point_data_mp(&tempPD_mp, 0);

    if (*prec < 64)
    { // setup time
      set_d(time, Final_d->time);
    }
    else
    { // setup time
      mp_to_d(time, Final_mp->time);
    }

    // refine to 'digits'
    refine_retVal = refine_digits_amp(T->outputLevel, digits_new, &T->latest_newton_residual_d, T->latest_newton_residual_mp, digits_in, &tempPD_d, &tempPD_mp, &prec_out, Final_d, Final_mp, *prec, time, OUT, e_d, e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // make sure refine was successful
    if (refine_retVal == 0) 
    { // success - update Final
      if (*prec < 64)
      { // update Final_d
        if (prec_out < 64)
        { // use tempPD_d
          point_cp_d(Final_d->point, tempPD_d.point);
        }
        else
        { // use tempPD_mp
          point_mp_to_d(Final_d->point, tempPD_mp.point);
        }
      }
      else
      { // update Final_mp
        if (prec_out < 64)
        { // use tempPD_d
          point_d_to_mp(Final_mp->point, tempPD_d.point);
        }
        else
        { // use tempPD_mp
          point_cp_mp(Final_mp->point, tempPD_mp.point);
        }
      }
    }
    else
    { // failure - this should not happen!
      retVal = 0;
      *smallest_nonzero_SV = *largest_zero_SV = 0;
    }

    // clear
    clear_point_data_d(&tempPD_d);
    clear_point_data_mp(&tempPD_mp);
  }

  return retVal;
}


