// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

// decompose a pure-dimensional witness set into its irreducible components 

void pureDecomp_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode);
int traceTrack(comp_d value_d, comp_mp value_mp, int *value_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);
void componentDecomposition(int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int trace_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode);

int monodromy(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec, double *norms);

void mergeTempGroups(int new_gp_number, int old_gp_number, int num_paths, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec);

int checkTraceArray(int curr_size, int *curr_array, double tol, comp_d *temp_gp_trace_d, comp_mp *temp_gp_trace_mp, int trace_prec);

void copy_slice_moving_slice(membership_slice_moving_t *sliceMover, prog_t *fullRankProg, vec_d start_d, vec_mp start_mp, mpq_t **start_rat, vec_d target_d, vec_mp target_mp, mpq_t **target_rat, comp_d gamma_d, comp_mp gamma_mp, mpq_t *gamma_rat, int MPType);


// do the break-up into irreducible components
void pureDecomp(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decomposes all pure-dim witness sets into irred comp   *
\***************************************************************/
{
  if (num_processes > 1)
  {
#ifdef _HAVE_MPI
    pureDecomp_par(sliceMover, fullRankProgs, fullRankProgsInfo, endPts_d, endPts_mp, endPts_amp, W, T, pathMod, my_id, num_processes, headnode);
#endif
  }
  else
  {
    pureDecomp_seq(sliceMover, fullRankProgs, fullRankProgsInfo, endPts_d, endPts_mp, endPts_amp, W, T, pathMod, my_id, num_processes, headnode);
  }

  return;
}

void pureDecomp_seq(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decomposes all pure-dim witness sets into irred comp   *
* using sequential                                              *
\***************************************************************/
{
  int i, num_codim = W->num_codim;
  FILE *OUT = NULL;
  char outName[] = "output_decomp", midName[] = "midpath_data";

  // setup OUT
  OUT = fopen(outName, "w");

  // loop over the codimensions
  for (i = 0; i < num_codim; i++)
    if (W->codim[i].num_set > 0)
    { // decompose codim_index 'i' 
      if (T->MPType == 0)
        pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], endPts_d[i], NULL, NULL, W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
      else if (T->MPType == 1)
        pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], NULL, endPts_mp[i], NULL, W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
      else
        pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], NULL, NULL, endPts_amp[i], W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
    }

  // close OUT
  fclose(OUT);

  return;
}

void pureDecomp_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decomposes codim 'codim_index' into irred components   *
\***************************************************************/
{
  int i, gammaSetup = 1, amp_trace_prec, num_points = W->codim[codim_index].num_set;
  int *traces_prec = T->MPType == 2 ? (int *)bmalloc(num_points * sizeof(int)) : NULL, *isClassified = (int *)bmalloc(num_points * sizeof(int));
  comp_d s_d[2], *traces_d = (T->MPType == 0 || T->MPType == 2) ? (comp_d *)bmalloc(num_points * sizeof(comp_d)) : NULL;
  comp_mp s_mp[2], *traces_mp = (T->MPType == 1 || T->MPType == 2) ? (comp_mp *)bmalloc(num_points * sizeof(comp_mp)) : NULL;
  vec_d proj_d, sliceVec_d;
  vec_mp proj_mp, sliceVec_mp;
  mpq_t s_rat[2][2], **proj_rat = NULL, **sliceVec_rat = NULL;
  FILE *MIDOUT = fopen(midName, "w");

  // allocate for component_nums
  W->codim[codim_index].num_components = 0;
  W->codim[codim_index].component_nums = (int *)bmalloc(num_points * sizeof(int));

  // display messages
  printf("\nCalculating traces for codimension %d.\n", W->codim[codim_index].codim);
  fprintf(OUT, "\n*****************************************************\n");   
  fprintf(OUT, "Calculating traces for codimension %d.\n", W->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // setup projection vector, slice vector and s
  if (T->MPType == 0)
  { // setup proj_d
    init_vec_d(proj_d, sliceMover->orig_variables);
    make_vec_random_d(proj_d, sliceMover->orig_variables);
    // setup sliceVec_d
    init_vec_d(sliceVec_d, sliceMover->B_d->rows);
    make_vec_random_d(sliceVec_d, sliceMover->B_d->rows);
    // setup s_d
    get_comp_rand_d(s_d[0]);
    get_comp_rand_d(s_d[1]);
  }
  else if (T->MPType == 1)
  { // setup proj_mp
    init_vec_mp(proj_mp, sliceMover->orig_variables);
    make_vec_random_mp(proj_mp, sliceMover->orig_variables);
    // setup sliceVec_mp
    init_vec_mp(sliceVec_mp, sliceMover->B_mp->rows);
    make_vec_random_mp(sliceVec_mp, sliceMover->B_mp->rows);
    // setup s_mp
    init_mp(s_mp[0]); init_mp(s_mp[1]);
    get_comp_rand_mp(s_mp[0]);
    get_comp_rand_mp(s_mp[1]);
  }
  else
  { // setup proj_d, proj_mp & proj_rat
    init_vec_d(proj_d, sliceMover->orig_variables);
    init_vec_mp2(proj_mp, sliceMover->orig_variables, sliceMover->curr_precision);
    init_vec_rat(proj_rat, sliceMover->orig_variables);
    make_vec_random_rat(proj_d, proj_mp, proj_rat, sliceMover->orig_variables, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);

    // setup sliceVec_d, sliceVec_mp & sliceVec_rat
    init_vec_d(sliceVec_d, sliceMover->B_mp->rows);
    init_vec_mp2(sliceVec_mp, sliceMover->B_mp->rows, sliceMover->curr_precision);
    init_vec_rat(sliceVec_rat, sliceMover->B_mp->rows);
    make_vec_random_rat(sliceVec_d, sliceVec_mp, sliceVec_rat, sliceMover->B_mp->rows, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);

    // setup s_d, s_mp & s_rat
    get_comp_rand_rat(s_d[0], s_mp[0], s_rat[0], sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
    get_comp_rand_rat(s_d[1], s_mp[1], s_rat[1], sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // setup memory and calculate the traces
  for (i = 0; i < num_points; i++)
  { // setup memory
    isClassified[i] = 0; 
    W->codim[codim_index].component_nums[i] = -1;

    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Calculating %d of %d\n", i, num_points);

    // see if deflation worked
    if (fullRankProgsInfo[i] == -1)
    { // deflation did not work - set isClassified to 1
      isClassified[i] = 1;

      // setup traces so that everything else works
      if (T->MPType == 0)
      {
        set_zero_d(traces_d[i]);
      }
      else if (T->MPType == 1)
      {
        init_mp(traces_mp[i]);
        set_zero_mp(traces_mp[i]);
      }
      else
      {
        traces_prec[i] = 52;
        set_zero_d(traces_d[i]);
      }
    }
    else
    { // calculate the trace
      if (T->MPType == 0)
      { // calculate using double precision
        calculateTrace(traces_d[i], NULL, NULL, sliceMover, fullRankProgs[i], &endPts_d[i], NULL, NULL, W->codim[codim_index].codim, i, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gammaSetup);
      }
      else if (T->MPType == 1)
      { // calculate in fixed multiprecision
        init_mp(traces_mp[i]);
        calculateTrace(NULL, traces_mp[i], NULL, sliceMover, fullRankProgs[i], NULL, &endPts_mp[i], NULL, W->codim[codim_index].codim, i, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gammaSetup);
      }
      else
      { // calculate using AMP
        calculateTrace(traces_d[i], traces_mp[i], &traces_prec[i], sliceMover, fullRankProgs[i], NULL, NULL, &endPts_amp[i], W->codim[codim_index].codim, i, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gammaSetup);
      }

      // update gammaSetup - since it has been setup
      gammaSetup = 0;
    }
  }

  // close MIDOUT
  fclose(MIDOUT);

  if (T->MPType == 2)
  { // to simplify AMP, we find the maximum precision and set all traces to that precision
    amp_trace_prec = 52;
    for (i = 0; i < num_points; i++)
      if (traces_prec[i] > amp_trace_prec)
        amp_trace_prec = traces_prec[i];

    // if all are in double precision, we do not have to do anything
    // if some are in higher precision, change all of them to the highest precision
    if (amp_trace_prec >= 64)
    { // set all traces to this precision
      for (i = 0; i < num_points; i++)
      {
        if (traces_prec[i] < 64)
        { // initialize traces_mp[i] and setup
          init_mp2(traces_mp[i], amp_trace_prec);
          d_to_mp(traces_mp[i], traces_d[i]);
        }
        else
        { // change the precision on traces_mp[i]
          change_prec_mp2(traces_mp[i], amp_trace_prec);
        }
        // set the precision correctly
        traces_prec[i] = amp_trace_prec;
      }
    }
  }
  else
    amp_trace_prec = 0;

  // do the actual classification now that we have traces
  componentDecomposition(isClassified, traces_d, traces_mp, amp_trace_prec, sliceMover, fullRankProgs, endPts_d, endPts_mp, endPts_amp, W, codim_index, T, OUT, midName, pathMod, my_id, num_processes, headnode);

  // clear memory
  clear_vec(proj_d, proj_mp, proj_rat, T->MPType);
  clear_vec(sliceVec_d, sliceVec_mp, sliceVec_rat, T->MPType);
  clear_d_mp_rat(s_d[0], s_mp[0], s_rat[0], T->MPType);
  clear_d_mp_rat(s_d[1], s_mp[1], s_rat[1], T->MPType);
  free(isClassified);
  if (T->MPType == 0)
  { // clear traces_d
    free(traces_d);
  }
  else if (T->MPType == 1)
  { // clear traces_mp
    for (i = num_points - 1; i >= 0; i--)
    { 
      clear_mp(traces_mp[i]);
    }
    free(traces_mp);
  }
  else
  { // clear traces_d, traces_mp & traces_prec
    for (i = num_points - 1; i >= 0; i--)
      if (traces_prec[i] < 64)
      {
        clear_d(traces_d[i]);
      }
      else
      {
        clear_mp(traces_mp[i]);
      }

    free(traces_d);
    free(traces_mp);
    free(traces_prec);
  }

  return;
}

void calculateTrace(comp_d trace_d, comp_mp trace_mp, int *trace_prec, membership_slice_moving_t *sliceMover, prog_t *fullRankProg, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, int codim, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], int setupGamma)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the trace for 'codim_index' & pathNum       *
\***************************************************************/
{
  int i, j, its = 0, maxIts = 20, retVal[2], rows, size = sliceMover->orig_variables;

  // setup rows
  if (T->MPType == 0 || T->MPType == 2)
    rows = sliceMover->B_d->rows;
  else
    rows = sliceMover->B_mp->rows;

  // see if the point has dimension 0
  if (codim == sliceMover->orig_variables - 1)
  { // since this point is zero-dimensional, the trace is 0
    if (T->MPType == 0)
    {
      set_zero_d(trace_d);
    }
    else if (T->MPType == 1)
    { 
      set_zero_mp(trace_mp);
    }
    else
    {
      *trace_prec = endPt_amp->curr_prec;
      if (*trace_prec < 64)
      {
        set_zero_d(trace_d);
      }
      else
      {
        init_mp2(trace_mp, *trace_prec);
        set_zero_mp(trace_mp);
      }
    }
  }
  else
  { // we need to actual track

    // finish the setup for sliceMover
    initialize_slice_moving_sliceVec(sliceMover, rows, T->MPType);
    final_setup_slice_moving(sliceMover, fullRankProg, T->MPType, T->AMP_max_prec, setupGamma);

    if (T->MPType == 0)
    { // use double precision
      double lambda;
      point_data_d PD;
      comp_d tempComp, q[3], lambda_s[2];

      // set gamma to 1 - just want to move the slice linearly
      set_one_d(sliceMover->gamma_d);

      // setup lambda & lambda_s
      lambda = 0.8;
      set_d(lambda_s[0], s_d[0]);
      set_d(lambda_s[1], s_d[1]);

      // setup PD
      init_point_data_d(&PD, size);
      point_cp_d(PD.point, endPt_d->endPt);
      set_one_d(PD.time);

      // setup startSliceVec == 0 and setup targetSliceVec == 0 for refinement - does not matter what it is, as long as it is setup
      sliceMover->startSliceVec_d->size = rows;
      sliceMover->targetSliceVec_d->size = rows;
      for (i = 0; i < rows; i++)
      {
        set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
        set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
      }

      // refine PD
      refine_d(&PD, T, OUT, sliceMover, slice_mover_eval_d);

      // setup q[0] to be the original point under projection
      set_zero_d(q[0]);
      for (i = 0; i < size; i++)
      {
        sum_mul_d(q[0], &proj_d->coord[i], &PD.point->coord[i]);
      }  

      // calculate q[1] & q[2]
      for (i = 0; i < 2; i++)
      { // initialize
        its = 0;

        // loop until success or too many iterations
        do
        { // increment the number of iterations
          its++;

          // setup targetSliceVec properly
          vec_mulcomp_d(sliceMover->targetSliceVec_d, v_d, lambda_s[i]);

          // track the path and find q[i+1]
          retVal[i] = traceTrack(q[i+1], NULL, NULL, sliceMover, endPt_d->endPt, NULL, 52, pathNum, proj_d, proj_mp, proj_rat, T, OUT, MIDOUT);

          // check for errors
          if (its > maxIts)
          { // this should never happen!
            printf("\nERROR: Failed to calculate traces after %d attempts! Try tightening the tracking tolerances.\n\n", maxIts);
            bexit(ERROR_LOOP_FAILURE);
          }
          else if (retVal[i] != 0)
          { // adjust lambda_s
            mul_rdouble_d(lambda_s[i], lambda_s[i], lambda);
          }
 
        } while (its <= maxIts && retVal[i] != 0);
      }

      // successful track - find the trace = (q1 - q0) / s0 - (q2 - q0) / s1
      sub_d(tempComp, q[1], q[0]);
      div_d(tempComp, tempComp, lambda_s[0]);
      sub_d(trace_d, q[2], q[0]);
      div_d(trace_d, trace_d, lambda_s[1]);
      sub_d(trace_d, tempComp, trace_d);

      // clear memory
      clear_point_data_d(&PD);
    }
    else if (T->MPType == 1)
    { // use fixed multi precision
      mpf_t lambda;
      point_data_mp PD;
      comp_mp tempComp, q[3], lambda_s[2];

      // initialize memory
      mpf_init(lambda);
      init_point_data_mp(&PD, size);
      init_mp(tempComp); init_mp(q[0]); init_mp(q[1]); init_mp(q[2]); 
      init_mp(lambda_s[0]); init_mp(lambda_s[1]);
      
      // set gamma to 1 - just want to move the slice linearly
      set_one_mp(sliceMover->gamma_mp); 

      // setup lambda & lambda_s
      mpf_set_ui(lambda, 8);
      mpf_div_ui(lambda, lambda, 10); // 8 / 10
      set_mp(lambda_s[0], s_mp[0]);
      set_mp(lambda_s[1], s_mp[1]);
      
      // setup PD
      point_cp_mp(PD.point, endPt_mp->endPt);
      set_one_mp(PD.time);

      // setup startSliceVec == 0 and setup targetSliceVec == 0 for refinement - does not matter what it is, as long as it is setup
      sliceMover->startSliceVec_mp->size = rows;
      sliceMover->targetSliceVec_mp->size = rows;
      for (i = 0; i < rows; i++)
      {
        set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
        set_zero_mp(&sliceMover->targetSliceVec_mp->coord[i]);
      }

      // refine PD
      refine_mp(&PD, T, OUT, sliceMover, slice_mover_eval_mp);

      // setup q[0] to be the original point under projection
      set_zero_mp(q[0]);
      for (i = 0; i < size; i++)
      {
        sum_mul_mp(q[0], &proj_mp->coord[i], &PD.point->coord[i]);
      }

      // calculate q[1] & q[2]
      for (i = 0; i < 2; i++)
      { // initialize
        its = 0;

        // loop until success or too many iterations
        do
        { // increment the number of iterations
          its++;

          // setup targetSliceVec properly
          vec_mulcomp_mp(sliceMover->targetSliceVec_mp, v_mp, lambda_s[i]);

          // track the path and find q[i+1]
          retVal[i] = traceTrack(NULL, q[i+1], NULL, sliceMover, NULL, endPt_mp->endPt, T->Precision, pathNum, proj_d, proj_mp, proj_rat, T, OUT, MIDOUT);

          // check for errors
          if (its > maxIts)
          { // this should never happen!
            printf("\nERROR: Failed to calculate traces after %d attempts! Try tightening the tracking tolerances.\n\n", maxIts);
            bexit(ERROR_LOOP_FAILURE);
          }
          else if (retVal[i] != 0)
          { // adjust lambda_s
            mul_rmpf_mp(lambda_s[i], lambda_s[i], lambda);
          }

        } while (its <= maxIts && retVal[i] != 0);
      }

      // successful track - find the trace = (q1 - q0) / s0 - (q2 - q0) / s1
      sub_mp(tempComp, q[1], q[0]);
      div_mp(tempComp, tempComp, lambda_s[0]);
      sub_mp(trace_mp, q[2], q[0]);
      div_mp(trace_mp, trace_mp, lambda_s[1]);
      sub_mp(trace_mp, tempComp, trace_mp);

      // clear memory
      mpf_clear(lambda);
      clear_point_data_mp(&PD);
      clear_mp(tempComp); clear_mp(q[0]); clear_mp(q[1]); clear_mp(q[2]);
      clear_mp(lambda_s[0]); clear_mp(lambda_s[1]);
    }
    else 
    { // use AMP
      int q_prec[3], min_digits = ceil(-log10(T->final_tolerance) + 1.5);
      mpq_t lambda;
      point_data_d PD_d, in_d;
      point_data_mp PD_mp, in_mp;
      comp_d tempComp_d, q_d[3], s_temp_d;
      comp_mp tempComp_mp, q_mp[3], s_temp_mp;
      mpq_t lambda_s[2][2];
      eval_struct_d e_d; 
      eval_struct_mp e_mp;

      // initialize memory
      mpq_init(lambda);
      init_point_data_d(&PD_d, size); init_point_data_d(&in_d, size);
      init_point_data_mp(&PD_mp, size); init_point_data_mp(&in_mp, size);
      init_mp(tempComp_mp); init_mp(q_mp[0]); init_mp(q_mp[1]); init_mp(q_mp[2]); 
      init_mp(s_temp_mp); 
      init_rat(lambda_s[0]); init_rat(lambda_s[1]);
      init_eval_struct_d(e_d, 0, 0, 0);
      init_eval_struct_mp(e_mp, 0, 0, 0);

      // set gamma to 1 - just want to move the slice linearly
      set_one_d(sliceMover->gamma_d);
      set_one_mp(sliceMover->gamma_mp);
      set_one_rat(sliceMover->gamma_rat);

      // setup lambda & lambda_s
      mpq_set_ui(lambda, 8, 10); // 8 / 10
      for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
          mpq_set(lambda_s[i][j], s_rat[i][j]);

      // setup startSliceVec == 0 and setup targetSliceVec == 0 for refinement - does not matter what it is, as long as it is setup
      sliceMover->startSliceVec_d->size = sliceMover->startSliceVec_mp->size = rows;
      sliceMover->targetSliceVec_d->size = sliceMover->targetSliceVec_mp->size = rows;
      for (i = 0; i < rows; i++)
      {
        set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
        set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
        set_zero_rat(sliceMover->startSliceVec_rat[i]);

        set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
        set_zero_mp(&sliceMover->targetSliceVec_mp->coord[i]);
        set_zero_rat(sliceMover->targetSliceVec_rat[i]);
      }

      // refine the endpoint to enough precision
      set_one_d(tempComp_d);
      if (endPt_amp->curr_prec < 64)
      { // setup in_d
        point_cp_d(in_d.point, endPt_amp->endPt_d);
        set_one_d(in_d.time);
      }
      else
      { // setup in_mp
        change_prec_point_data_mp(&in_mp, endPt_amp->curr_prec);
        point_cp_mp(in_mp.point, endPt_amp->endPt_mp);
        set_one_mp(in_mp.time);
      }
      refine_digits_amp(T->outputLevel, min_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, floor(-log10(T->currentNewtonTol)), &PD_d, &PD_mp, &q_prec[0], &in_d, &in_mp, endPt_amp->curr_prec, tempComp_d, OUT, &e_d, &e_mp, sliceMover, sliceMover, slice_mover_eval_d, slice_mover_eval_mp, slice_mover_change_prec);

      // setup q[0]
      if (q_prec[0] < 64)
      { // setup q_d[0] to be the original point under projection
        set_zero_d(q_d[0]);
        for (i = 0; i < size; i++)
        {
          sum_mul_d(q_d[0], &proj_d->coord[i], &PD_d.point->coord[i]);
        }
      }
      else
      { // setup q_mp[0] to be the original point under projection
        set_zero_mp(q_mp[0]);
        for (i = 0; i < size; i++)
        {
          sum_mul_mp(q_mp[0], &proj_mp->coord[i], &PD_mp.point->coord[i]);
        }
      }

      // calculate q[1] & q[2]
      for (i = 0; i < 2; i++)
      { // initialize
        its = 0;

        // loop until success or too many iterations
        do
        { // increment the number of iterations
          its++;

          // setup targetSliceVec properly
          sliceMover->targetSliceVec_d->size = sliceMover->targetSliceVec_mp->size = rows;
          for (j = 0; j < rows; j++)
          {
            mul_rat(sliceMover->targetSliceVec_rat[j], v_rat[j], lambda_s[i]);

            mpf_set_q(sliceMover->targetSliceVec_mp->coord[j].r, sliceMover->targetSliceVec_rat[j][0]);
            mpf_set_q(sliceMover->targetSliceVec_mp->coord[j].i, sliceMover->targetSliceVec_rat[j][1]);
            sliceMover->targetSliceVec_d->coord[j].r = mpq_get_d(sliceMover->targetSliceVec_rat[j][0]);
            sliceMover->targetSliceVec_d->coord[j].i = mpq_get_d(sliceMover->targetSliceVec_rat[j][1]);
          }

          // track the path and find q[i+1]
          retVal[i] = traceTrack(q_d[i+1], q_mp[i+1], &q_prec[i+1], sliceMover, PD_d.point, PD_mp.point, q_prec[0], pathNum, proj_d, proj_mp, proj_rat, T, OUT, MIDOUT);

          // check for errors
          if (its > maxIts)
          { // this should never happen!
            printf("\nERROR: Failed to calculate traces after %d attempts! Try tightening the tracking tolerances.\n\n", maxIts);
            bexit(ERROR_LOOP_FAILURE);
          }
          else if (retVal[i] != 0)
          { // adjust lambda_s
            for (j = 0; j < 2; j++)
              mpq_mul(lambda_s[i][j], lambda_s[i][j], lambda);
          }

        } while (its <= maxIts && retVal[i] != 0);
      }

      // successful track - find the trace = (q1 - q0) / s0 - (q2 - q0) / s1 and exit loop
      *trace_prec = MAX(q_prec[0], q_prec[1]);
      *trace_prec = MAX(*trace_prec, q_prec[2]);

      if (*trace_prec < 64)
      { // all of the q's are in double precision
        sub_d(tempComp_d, q_d[1], q_d[0]);
        s_temp_d->r = mpq_get_d(lambda_s[0][0]);
        s_temp_d->i = mpq_get_d(lambda_s[0][1]);        
        div_d(tempComp_d, tempComp_d, s_temp_d);
        sub_d(trace_d, q_d[2], q_d[0]);
        s_temp_d->r = mpq_get_d(lambda_s[1][0]);
        s_temp_d->i = mpq_get_d(lambda_s[1][1]);
        div_d(trace_d, trace_d, s_temp_d);
        sub_d(trace_d, tempComp_d, trace_d);
      }
      else
      { // convert all of the q's to the maximum precision of the q's
        for (i = 0; i < 3; i++)
          if (q_prec[i] < 64)
          { // set precision and convert _d to _mp
            setprec_mp(q_mp[i], *trace_prec);
            d_to_mp(q_mp[i], q_d[i]);
          }
          else
          { // change the precision
            change_prec_mp2(q_mp[i], *trace_prec); 
          }

        // set precision on trace_mp
        init_mp2(trace_mp, *trace_prec);
        // set precision on s_temp_mp
        setprec_mp(s_temp_mp, *trace_prec);

        // calculate the trace
        sub_mp(tempComp_mp, q_mp[1], q_mp[0]);
        mpf_set_q(s_temp_mp->r, lambda_s[0][0]);
        mpf_set_q(s_temp_mp->i, lambda_s[0][1]);
        div_mp(tempComp_mp, tempComp_mp, s_temp_mp);
        sub_mp(trace_mp, q_mp[2], q_mp[0]);
        mpf_set_q(s_temp_mp->r, lambda_s[1][0]);
        mpf_set_q(s_temp_mp->i, lambda_s[1][1]);
        div_mp(trace_mp, trace_mp, s_temp_mp);
        sub_mp(trace_mp, tempComp_mp, trace_mp);
      }

      // clear memory
      mpq_clear(lambda);
      clear_point_data_d(&PD_d); clear_point_data_d(&in_d);
      clear_point_data_mp(&PD_mp); clear_point_data_mp(&in_mp);
      clear_mp(tempComp_mp); clear_mp(q_mp[0]); clear_mp(q_mp[1]); clear_mp(q_mp[2]);
      clear_mp(s_temp_mp); 
      clear_rat(lambda_s[0]); clear_rat(lambda_s[1]);
      clear_eval_struct_d(e_d);
      clear_eval_struct_mp(e_mp);
    }
  }

  return;
}

int traceTrack(comp_d value_d, comp_mp value_mp, int *value_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - successful tracking, otherwise track error *
* NOTES: tracks the path and find the value of the projection   *
\***************************************************************/
{
  int i, size, retVal = 0;
  endgame_data_t endPt;

  // initialize endPt
  init_endgame_data(&endPt, T->Precision);

  // move the slice
  retVal = slice_moving_track(&endPt, sliceMover, startPt_d, startPt_mp, startPt_prec, pathNum_startPt, 0, T, OUT, MIDOUT);

  // calculate the trace if the slice moving was successful
  if (!retVal)
  { // find the value of the projection
    if (T->MPType == 0)
    { // use _d
      size = proj_d->size;

      // refine end point
      set_zero_d(endPt.PD_d.time);
      refine_d(&endPt.PD_d, T, OUT, sliceMover, slice_mover_eval_d);

      set_zero_d(value_d);
      for (i = 0; i < size; i++)
      {
        sum_mul_d(value_d, &proj_d->coord[i], &endPt.PD_d.point->coord[i]);
      }
    }
    else if (T->MPType == 1)
    { // use _mp
      size = proj_mp->size;

      // refine end point
      set_zero_mp(endPt.PD_mp.time);
      refine_mp(&endPt.PD_mp, T, OUT, sliceMover, slice_mover_eval_mp);

      set_zero_mp(value_mp);
      for (i = 0; i < size; i++)
      {
        sum_mul_mp(value_mp, &proj_mp->coord[i], &endPt.PD_mp.point->coord[i]);
      }
    }
    else
    { // use _amp
      int out_prec, min_digits = ceil(-log10(T->final_tolerance) + 1.5);
      comp_d time;
      point_data_d PD_d;
      point_data_mp PD_mp;
      eval_struct_d e_d; 
      eval_struct_mp e_mp;

      // initialize
      init_point_data_d(&PD_d, 0);
      init_point_data_mp(&PD_mp, 0);
      init_eval_struct_d(e_d, 0, 0, 0); 
      init_eval_struct_mp(e_mp, 0, 0, 0);

      *value_prec = endPt.prec;
      size = proj_d->size;

      // refine the endpoint to enough precision
      set_zero_d(time);
      refine_digits_amp(T->outputLevel, min_digits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, floor(-log10(T->currentNewtonTol)), &PD_d, &PD_mp, &out_prec, &endPt.PD_d, &endPt.PD_mp, endPt.prec, time, OUT, &e_d, &e_mp, sliceMover, sliceMover, slice_mover_eval_d, slice_mover_eval_mp, slice_mover_change_prec);

      if (out_prec < 64)
      { // calculate value using double precision
        set_zero_d(value_d);
        for (i = 0; i < size; i++)
        {
          sum_mul_d(value_d, &proj_d->coord[i], &PD_d.point->coord[i]);
        }
      }
      else
      { // calculate value using multi precision
        setprec_mp(value_mp, out_prec);
        set_zero_mp(value_mp);

        // set the precision for proj_mp
        change_size_vec_mp(proj_mp, size);
        setprec_vec_mp(proj_mp, endPt.prec);
        proj_mp->size = size;

        // update proj_mp & calculate value_mp
        for (i = 0; i < size; i++)
        {
          mpf_set_q(proj_mp->coord[i].r, proj_rat[i][0]);
          mpf_set_q(proj_mp->coord[i].i, proj_rat[i][1]);

          sum_mul_mp(value_mp, &proj_mp->coord[i], &PD_mp.point->coord[i]);
        }
      }

      // clear
      clear_point_data_d(&PD_d);
      clear_point_data_mp(&PD_mp);
      clear_eval_struct_d(e_d);
      clear_eval_struct_mp(e_mp);
    }
  }

  // clear endPt
  clear_endgame_data(&endPt);

  return retVal;
}

void componentDecomposition(int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int trace_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: uses the traces and monodromy to do the decomposition  *
*  Assume that all traces are in the same precision             *
\***************************************************************/
{
  int i, its, num_no_connections, retVal, num_to_classify = 0, num_temp_gps = 0, num_points = W->codim[codim_index].num_set;
  int *temp_gp_nums = (int *)bmalloc(num_points * sizeof(int));
  comp_d *temp_gp_trace_d = NULL;
  comp_mp *temp_gp_trace_mp = NULL;

  // count the number that we need to classify and assign temporary groups
  for (i = 0; i < num_points; i++)
    if (!isClassified[i])
    { // assign a temporary group number
      temp_gp_nums[i] = num_to_classify;
      // increment num_to_classify
      num_to_classify++;
    }
    else
    { // set temporary group to -1 since they will be ignored
      temp_gp_nums[i] = -1;
    }
  num_temp_gps = num_to_classify;

  // setup the trace of the temporary groups
  if (T->MPType == 0 || (T->MPType == 2 && trace_prec < 64))
  { // use temp_gp_trace_d
    trace_prec = 52;
    temp_gp_trace_d = (comp_d *)bmalloc(num_temp_gps * sizeof(comp_d));

    // calculate the trace of each group
    for (i = 0; i < num_temp_gps; i++)
    {
      calculateGroupTrace_d(temp_gp_trace_d[i], num_points, temp_gp_nums, i, traces_d);
    }
  }
  else 
  { // use temp_gp_trace_mp
    if (T->MPType == 1) // use current fixed precision, otherwise trace_prec is already setup
      trace_prec = T->Precision;
    
    temp_gp_trace_mp = (comp_mp *)bmalloc(num_temp_gps * sizeof(comp_mp));
    for (i = 0; i < num_temp_gps; i++)
    { // initialize
      init_mp2(temp_gp_trace_mp[i], trace_prec);
      // calculate the trace of each group
      calculateGroupTrace_mp(temp_gp_trace_mp[i], num_points, temp_gp_nums, i, traces_mp);
    }
  }

  // check for completed components - updates isClassified, temp_gp_nums, temp_gp_trace, num_to_classify, num_temp_gps for each completed component found
  checkCompleteComponents(W, codim_index, T->final_tol_times_mult, &num_to_classify, &isClassified, &num_temp_gps, &temp_gp_nums, &temp_gp_trace_d, &temp_gp_trace_mp, trace_prec);

  // use monodromy if there are unclassified endpoints and enough temporary groups to justify it
  if (num_to_classify > 0 && num_temp_gps > T->max_num_pts_for_trace && T->max_num_mon_linears > 0 && T->max_num_bad_loops_in_mon > 0)
  { // display messages
    printf("\nUsing monodromy to decompose codimension %d.\n", W->codim[codim_index].codim);
    fprintf(OUT, "\n*****************************************************\n");
    fprintf(OUT, "Using monodromy to decompose codimension %d.\n", W->codim[codim_index].codim);
    fprintf(OUT, "*****************************************************\n");

    // find a double approximation to the norms - helps to classify them faster
    double *norms = (double *)bmalloc(num_points * sizeof(double));
    if (T->MPType == 0)
    { // use Pts_d
      for (i = 0; i < num_points; i++)
      {
        norms[i] = infNormVec_d(W->codim[codim_index].witnessPts_d[i].endPt);
      }
    }
    else if (T->MPType == 1)
    { // use Pts_mp
      for (i = 0; i < num_points; i++)
      {
        norms[i] = infNormVec_mp(W->codim[codim_index].witnessPts_mp[i].endPt);
      }
    }
    else
    { // use Pts_amp
      for (i = 0; i < num_points; i++)
        if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
        { // endPt_d
          norms[i] = infNormVec_d(W->codim[codim_index].witnessPts_amp[i].endPt_d);
        }
        else
        { // endPt_mp
          norms[i] = infNormVec_mp(W->codim[codim_index].witnessPts_amp[i].endPt_mp);
        }
    }

    // do monodromy to find connections
    its = num_no_connections = 0;
    while (num_to_classify > 0 && num_temp_gps > T->max_num_pts_for_trace && its < T->max_num_mon_linears && num_no_connections < T->max_num_bad_loops_in_mon)
    { // perform a monodromy loop and update the structures if any connections were made
      printf("Performing monodromy loops: %d point%s left to classify\n", num_to_classify, num_to_classify == 1 ? "" : "s");

      retVal = monodromy(sliceMover, fullRankProgs, endPts_d, endPts_mp, endPts_amp, W, codim_index, T, OUT, midName, T->final_tol_times_mult, &num_to_classify, &isClassified, &num_temp_gps, &temp_gp_nums, &temp_gp_trace_d, &temp_gp_trace_mp, trace_prec, norms);

      // check to see if any connections were made
      if (retVal > 0)
      { // made connections
        num_no_connections = 0;
      }
      else
      { // made no connections
        num_no_connections++;
      }

      // increment the number of iterations
      its++;   
    }

    // free norms
    free(norms);
  }

  // use trace test to finish the decomposition if there are unclassified endpoints
  if (num_to_classify > 0)
  { // display messages
    printf("\nUsing combinatorial trace test to decompose codimension %d.\n", W->codim[codim_index].codim);

    // complete breakup using traces
    traceDecomposition(W, codim_index, T->final_tol_times_mult, &num_to_classify, &isClassified, &num_temp_gps, &temp_gp_nums, &temp_gp_trace_d, &temp_gp_trace_mp, trace_prec);
  }

  // clear memory
  free(temp_gp_nums);
  if (num_temp_gps > 0)
  { // clear for extra temporary groups - hopefully there are none!
    free(temp_gp_trace_d);
    if (trace_prec >= 64)
    {
      for (i = num_temp_gps - 1; i >= 0; i--)
      { 
        clear_mp(temp_gp_trace_mp[i]);
      }
    }
    free(temp_gp_trace_mp);    
  }

  return;
}

void calculateGroupTrace_d(comp_d gp_trace, int num_points, int *gp_nums, int curr_gp, comp_d *traces)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the trace of the current group              *
\***************************************************************/
{
  int i;

  set_zero_d(gp_trace);
  for (i = 0; i < num_points; i++)
    if (gp_nums[i] == curr_gp)
      add_d(gp_trace, gp_trace, traces[i]);

  return;
}

void calculateGroupTrace_mp(comp_mp gp_trace, int num_points, int *gp_nums, int curr_gp, comp_mp *traces)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the trace of the current group              *
\***************************************************************/
{
  int i;

  set_zero_mp(gp_trace);
  for (i = 0; i < num_points; i++)
    if (gp_nums[i] == curr_gp)
      add_mp(gp_trace, gp_trace, traces[i]);

  return;
}

void checkCompleteComponents(witness_t *W, int codim_index, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: checks for a trace of zero - if found - updates all    *
\***************************************************************/
{
  int i, j, count, remove_gp, num_complete = 0, *gp_complete = NULL, num_points = W->codim[codim_index].num_set;

  // loop through to find the temporary groups that have a trace of zero
  if (trace_prec < 64)
  { // use _d
    for (i = 0; i < *num_temp_gps; i++)
      if (d_abs_d((*temp_gp_trace_d)[i]) < tol)
      { // this one has a trace of 0!!
        gp_complete = (int *)brealloc(gp_complete, (num_complete + 1) * sizeof(int));
        gp_complete[num_complete] = i;
        num_complete++;
      }
  }
  else
  { // use _mp
    for (i = 0; i < *num_temp_gps; i++)
      if (d_abs_mp((*temp_gp_trace_mp)[i]) < tol)
      { // this one has a trace of 0!!
        gp_complete = (int *)brealloc(gp_complete, (num_complete + 1) * sizeof(int));
        gp_complete[num_complete] = i;
        num_complete++;
      }
  }

  // see if we have any complete groups
  if (num_complete > 0)
  { // update everything that needs updated
    for (i = 0; i < num_complete; i++)
    { // find the points that are in this temporary group and assign them to a component
      for (j = 0; j < num_points; j++)
        if ((*temp_gp_nums)[j] == gp_complete[i])
        { // this one is in the current complete component
          W->codim[codim_index].component_nums[j] = W->codim[codim_index].num_components;
          // update the classified status
          (*isClassified)[j] = 1;
          (*num_to_classify)--;
        }
      // increment the number of components that we have found
      W->codim[codim_index].num_components++;
    }

    // update the temporary group structures
    if (trace_prec < 64)
    { // update _d
      comp_d *temp_traces_d = *temp_gp_trace_d;
      *temp_gp_trace_d = (comp_d *)bmalloc((*num_temp_gps - num_complete) * sizeof(comp_d));

      // remove the completed groups and keep the non-completed groups
      count = 0;
      for (i = 0; i < *num_temp_gps; i++)
      { // see if we should remove this group
        remove_gp = 0;
        for (j = 0; j < num_complete && !remove_gp; j++)
          if (i == gp_complete[j])
            remove_gp = 1;

        if (!remove_gp)
        { // this group is not complete so it will be kept
          set_d((*temp_gp_trace_d)[count], temp_traces_d[i]);
          // increment count
          count++;
        }
      }

      // update num_temp_gps
      *num_temp_gps = *num_temp_gps - num_complete;
      // free old memory
      free(temp_traces_d);
    }
    else
    { // update _mp
      comp_mp *temp_traces_mp = *temp_gp_trace_mp;
      *temp_gp_trace_mp = (comp_mp *)bmalloc((*num_temp_gps - num_complete) * sizeof(comp_mp));

      // remove the completed groups and keep the non-completed groups
      count = 0;
      for (i = 0; i < *num_temp_gps; i++)
      { // see if we should remove this group
        remove_gp = 0;
        for (j = 0; j < num_complete && !remove_gp; j++)
          if (i == gp_complete[j])
            remove_gp = 1;

        if (!remove_gp)
        { // this group is not complete so it will be kept
          init_mp2((*temp_gp_trace_mp)[count], trace_prec);
          set_mp((*temp_gp_trace_mp)[count], temp_traces_mp[i]);
          // increment count
          count++;
        }

        // clear old memory
        clear_mp(temp_traces_mp[i]);
      }

      // update num_temp_gps
      *num_temp_gps = *num_temp_gps - num_complete;
      // free old memory
      free(temp_traces_mp);
    }

    // update the group numbers
    for (i = 0; i < num_points; i++)
      if (!(*isClassified)[i])
      { // count how many groups were removed before this one
        count = 0;
        for (j = 0; j < num_complete; j++)
          if ((*temp_gp_nums)[i] > gp_complete[j])
            count++;

        // update the group number based on count
        (*temp_gp_nums)[i] = (*temp_gp_nums)[i] - count;
      }
      else
      { // set to -1 since it has been classified
        (*temp_gp_nums)[i] = -1;
      }
  }

  // clear memory
  free(gp_complete);

  return;
}

int monodromy(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec, double *norms)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of mergers created by this monodromy    *
* NOTES: performs one loop of monodromy on all of the           *
* unclassified points for this codim and update accordingly     *
\***************************************************************/
{
  int i, num_connections = 0, pathNum, size;
  int *loop_results = (int *)bmalloc(*num_to_classify * sizeof(int));
  FILE *MIDOUT = fopen(midName, "w");

  if (T->MPType == 0)
  { // monodromy using _d
    comp_d gamma_out_d, gamma_in_d;
    vec_d v_out_d, v_in_d;

    // setup gamma_out_d & gamma_in_d
    get_comp_rand_d(gamma_out_d);
    get_comp_rand_d(gamma_in_d);  

    // setup size
    size = sliceMover->B_d->rows;

    // setup v_out_d & v_in_d
    init_vec_d(v_out_d, size);
    init_vec_d(v_in_d, size);
    v_out_d->size = v_in_d->size = size;
    // out is random, in is zero
    make_vec_random_d(v_out_d, size);
    for (i = 0; i < size; i++)
    {
      set_zero_d(&v_in_d->coord[i]);
    }
    
    // do a monodromy loop for each of the points that still need to be classified
    pathNum = 0;
    for (i = 0; i < *num_to_classify; i++)
    { // find the next one that needs to be classified
      while ((*isClassified)[pathNum])
        pathNum++;

      // finish the setup for sliceMover
      initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
      final_setup_slice_moving(sliceMover, fullRankProgs[pathNum], T->MPType, T->AMP_max_prec, 0);

      // perform the monodromy loop and see where it ends up at
      loop_results[i] = monodromyLoop(sliceMover, fullRankProgs[pathNum], &endPts_d[pathNum], NULL, NULL, W, codim_index, pathNum, T, OUT, MIDOUT, tol, *isClassified, gamma_out_d, NULL, NULL, gamma_in_d, NULL, NULL, v_out_d, NULL, NULL, v_in_d, NULL, NULL, norms);

      // move past this path number
      pathNum++;
    }

    // clear v_out_d & v_in_d
    clear_vec_d(v_out_d);
    clear_vec_d(v_in_d);
  }
  else if (T->MPType == 1)
  { // monodromy using _mp
    comp_mp gamma_out_mp, gamma_in_mp;
    vec_mp v_out_mp, v_in_mp;

    // setup size
    size = sliceMover->B_mp->rows;

    // init gammas & vs
    init_mp(gamma_out_mp); init_mp(gamma_in_mp);
    init_vec_mp(v_out_mp, size); init_vec_mp(v_in_mp, size);

    // setup gammas
    get_comp_rand_mp(gamma_out_mp);
    get_comp_rand_mp(gamma_in_mp);

    // out is random, in is zero
    v_out_mp->size = v_in_mp->size = size;
    make_vec_random_mp(v_out_mp, size);
    for (i = 0; i < size; i++)
    {
      set_zero_mp(&v_in_mp->coord[i]);
    }
     
    // do a monodromy loop for each of the points that still need to be classified
    pathNum = 0;
    for (i = 0; i < *num_to_classify; i++)
    { // find the next one that needs to be classified
      while ((*isClassified)[pathNum])
        pathNum++;

      // finish the setup for sliceMover
      initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
      final_setup_slice_moving(sliceMover, fullRankProgs[pathNum], T->MPType, T->AMP_max_prec, 0);

      // perform the monodromy loop and see where it ends up at
      loop_results[i] = monodromyLoop(sliceMover, fullRankProgs[pathNum], NULL, &endPts_mp[pathNum], NULL, W, codim_index, pathNum, T, OUT, MIDOUT, tol, *isClassified, NULL, gamma_out_mp, NULL, NULL, gamma_in_mp, NULL, NULL, v_out_mp, NULL, NULL, v_in_mp, NULL, norms);

      // move past this path number
      pathNum++;
    }

    // clear gammas & vs
    clear_mp(gamma_out_mp); clear_mp(gamma_in_mp);
    clear_vec_mp(v_out_mp); clear_vec_mp(v_in_mp);
  }
  else
  { // monodromy using AMP
    comp_d gamma_out_d, gamma_in_d;
    comp_mp gamma_out_mp, gamma_in_mp;
    mpq_t gamma_out_rat[2], gamma_in_rat[2];
    vec_d v_out_d, v_in_d;
    vec_mp v_out_mp, v_in_mp;
    mpq_t **v_out_rat, **v_in_rat;

    // setup gammas
    get_comp_rand_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
    get_comp_rand_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, sliceMover->curr_precision, T->AMP_max_prec, 1, 1);

    // find size
    size = sliceMover->B_d->rows;
    // setup v
    init_vec_d(v_out_d, size);
    init_vec_mp2(v_out_mp, size, sliceMover->curr_precision);
    init_vec_rat(v_out_rat, size);
    init_vec_d(v_in_d, size);
    init_vec_mp2(v_in_mp, size, sliceMover->curr_precision);
    init_vec_rat(v_in_rat, size);
    // out is random, in is zero
    make_vec_random_rat(v_out_d, v_out_mp, v_out_rat, size, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);
    v_in_d->size = v_in_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&v_in_d->coord[i]);
      set_zero_mp(&v_in_mp->coord[i]);
      set_zero_rat(v_in_rat[i]);
    }

    // do a monodromy loop for each of the points that still need to be classified
    pathNum = 0;
    for (i = 0; i < *num_to_classify; i++)
    { // find the next one that needs to be classified
      while ((*isClassified)[pathNum])
        pathNum++;

      // finish the setup for sliceMover
      initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
      final_setup_slice_moving(sliceMover, fullRankProgs[pathNum], T->MPType, T->AMP_max_prec, 0);

      // perform the monodromy loop and see where it ends up at
      loop_results[i] = monodromyLoop(sliceMover, fullRankProgs[pathNum], NULL, NULL, &endPts_amp[pathNum], W, codim_index, pathNum, T, OUT, MIDOUT, tol, *isClassified, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, v_out_d, v_out_mp, v_out_rat, v_in_d, v_in_mp, v_in_rat, norms);

      // move past this path number
      pathNum++;
    }

    // clear gammas & vs
    clear_d_mp_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, T->MPType);
    clear_d_mp_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, T->MPType);
    clear_vec(v_out_d, v_out_mp, v_out_rat, T->MPType);
    clear_vec(v_in_d, v_in_mp, v_in_rat, T->MPType);
  }

  // now that we have the results, update the temporary groups and find the number of new connections made by this monodromy loop
  num_connections = updateTempGroupsFromLoops(W->codim[codim_index].num_set, *num_to_classify, *isClassified, loop_results, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);

  if (num_connections > 0)
  { // check for completed components and update accordingly
    checkCompleteComponents(W, codim_index, T->final_tol_times_mult, num_to_classify, isClassified, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);
  }

  // clear
  fclose(MIDOUT);
  free(loop_results);

  return num_connections;
}

int basicMonodromyLoop(endgame_data_t *monodromyPt, membership_slice_moving_t *sliceMover, prog_t *fullRankProgs, point_d Pt_d, point_mp Pt_mp, int Pt_prec, witness_t *W, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - tracking success, 1 - tracking failure     *
* NOTES: performs one loop of monodromy on the given point      *
\***************************************************************/
{
  int rV = 0;

  // setup to track 'out'
  copy_slice_moving_slice(sliceMover, fullRankProgs, vec_in_d, vec_in_mp, vec_in_rat, vec_out_d, vec_out_mp, vec_out_rat, gamma_out_d, gamma_out_mp, gamma_out_rat, T->MPType);

  // track 'out'
  rV = slice_moving_track(monodromyPt, sliceMover, Pt_d, Pt_mp, Pt_prec, pathNum, 0, T, OUT, MIDOUT);

  // check for tracking success
  if (!rV)
  { // setup to track 'in'
    int startPt_prec;
    point_d startPt_d;
    point_mp startPt_mp;

    init_point_d(startPt_d, 0);
    init_point_mp(startPt_mp, 0);

    copy_slice_moving_slice(sliceMover, fullRankProgs, vec_out_d, vec_out_mp, vec_out_rat, vec_in_d, vec_in_mp, vec_in_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, T->MPType);

    // setup startPt for tracking 'in'
    startPt_prec = monodromyPt->prec;
    if (startPt_prec < 64)
    { // use _d
      point_cp_d(startPt_d, monodromyPt->PD_d.point);
    }
    else
    { // use _mp
      setprec_point_mp(startPt_mp, startPt_prec);
      point_cp_mp(startPt_mp, monodromyPt->PD_mp.point);
    }

    // track 'in'
    rV = slice_moving_track(monodromyPt, sliceMover, startPt_d, startPt_mp, startPt_prec, pathNum, 0, T, OUT, MIDOUT);
  
    clear_point_d(startPt_d);
    clear_point_mp(startPt_mp);
  }

  return rV;
}

int monodromyLoop(membership_slice_moving_t *sliceMover, prog_t *fullRankProgs, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, witness_t *W, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, double tol, int *isClassified, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat, double *norms) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: path number for where this monodromy loop ended*
*  if there was a tracking error, it returns the input pathNum  *
* NOTES: performs one loop of monodromy on the given pathNum    *
\***************************************************************/
{
  int i, size, retVal, matchingPathNum, startPt_prec;
  double curr_norm;
  point_d startPt_d;
  point_mp startPt_mp;
  endgame_data_t endPt;

  init_point_d(startPt_d, 0);
  init_point_mp(startPt_mp, 0);
  init_endgame_data(&endPt, T->Precision);

  if (T->MPType == 0)
    retVal = basicMonodromyLoop(&endPt, sliceMover, fullRankProgs, endPt_d->endPt, NULL, 52, W, pathNum, T, OUT, MIDOUT, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, vec_out_d, vec_out_mp, vec_out_rat, vec_in_d, vec_in_mp, vec_in_rat);
  else if (T->MPType == 1)
    retVal = basicMonodromyLoop(&endPt, sliceMover, fullRankProgs, NULL, endPt_mp->endPt, T->Precision, W, pathNum, T, OUT, MIDOUT, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, vec_out_d, vec_out_mp, vec_out_rat, vec_in_d, vec_in_mp, vec_in_rat);
  else
    retVal = basicMonodromyLoop(&endPt, sliceMover, fullRankProgs, endPt_amp->endPt_d, endPt_amp->endPt_mp, endPt_amp->curr_prec, W, pathNum, T, OUT, MIDOUT, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, vec_out_d, vec_out_mp, vec_out_rat, vec_in_d, vec_in_mp, vec_in_rat);

  if (retVal)
  { // there was a tracking failure - return pathNum
    matchingPathNum = pathNum;
  }
  else
  { // see where the endpoint returned - check against those that are not classified
    mpf_t norm;
    mpf_init(norm);

    // default - return the starting path number
    matchingPathNum = pathNum;

    // project onto the original coordinates and find the norm
    startPt_prec = endPt.prec;
    if (endPt.prec < 64)
    {
      point_cp_d(startPt_d, endPt.PD_d.point);
      startPt_d->size = sliceMover->orig_variables;
      curr_norm = infNormVec_d(startPt_d);
    }
    else
    {
      setprec_point_mp(startPt_mp, endPt.prec);
      point_cp_mp(startPt_mp, endPt.PD_mp.point);
      startPt_mp->size = sliceMover->orig_variables;
      curr_norm = infNormVec_mp(startPt_mp);
    }

    // find the ones that are not classified and have a similar norm
    size = W->codim[codim_index].num_set;
    for (i = 0; i < size; i++)
      if (!isClassified[i] && fabs(curr_norm - norms[i]) < MAX(tol,1e-14))
      { // subtract the points to see if they are indeed equal
        if (T->MPType == 0)
        { // subtract using _d
          findDiff_point(norm, startPt_d, startPt_mp, startPt_prec, W->codim[codim_index].witnessPts_d[i].endPt, NULL, 52);
        }
        else if (T->MPType == 1)
        { // subtract using _mp
          findDiff_point(norm, startPt_d, startPt_mp, startPt_prec, NULL, W->codim[codim_index].witnessPts_mp[i].endPt, T->Precision);
        }
        else
        { // subtract using _amp
          findDiff_point(norm, startPt_d, startPt_mp, startPt_prec, W->codim[codim_index].witnessPts_amp[i].endPt_d, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);
        }

        // check to see if they are the same
        if (mpf_cmp_d(norm, tol) < 0)
        { // they are the same - store the path number and break out of the loop
          matchingPathNum = i;
          break;
        }
      }

    mpf_clear(norm);
  }

  clear_point_d(startPt_d);
  clear_point_mp(startPt_mp);
  clear_endgame_data(&endPt);
  
  return matchingPathNum;
}

int updateTempGroupsFromLoops(int num_paths, int num_to_classify, int *isClassified, int *loop_results, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of new connections made by loop_results *
* NOTES: finds the temporary group that should be together based*
*  on monodromy loops                                           *
\***************************************************************/
{
  int i, pathNum, new_connections = 0, old_gp_number, new_gp_number;

  // loop over the unclassified paths to see check for new connections
  pathNum = 0;
  for (i = 0; i < num_to_classify; i++)
  { // find the path number of the next one that needs to be classified
    while (isClassified[pathNum])
      pathNum++;

    // determine if pathNum & loop_results[i] have a different group number
    if ((*temp_gp_nums)[pathNum] != (*temp_gp_nums)[loop_results[i]])
    { // we have a new connection
      new_connections++;

      // find the new and old group number (MIN & MAX)
      new_gp_number = MIN((*temp_gp_nums)[pathNum], (*temp_gp_nums)[loop_results[i]]);
      old_gp_number = MAX((*temp_gp_nums)[pathNum], (*temp_gp_nums)[loop_results[i]]);

      // merge these temporary groups
      mergeTempGroups(new_gp_number, old_gp_number, num_paths, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);
    }

    // move past pathNum
    pathNum++;
  }

  return new_connections;
}

void mergeTempGroups(int new_gp_number, int old_gp_number, int num_paths, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: merges old_gp_number & new_gp_number                   *
\***************************************************************/
{
  int j;

  // combine these temporary groups and move all the ones past old_gp_number up by 1
  for (j = 0; j < num_paths; j++)
  {
    if ((*temp_gp_nums)[j] == old_gp_number)
      (*temp_gp_nums)[j] = new_gp_number;
    else if ((*temp_gp_nums)[j] > old_gp_number)
      (*temp_gp_nums)[j] = (*temp_gp_nums)[j] - 1;
  }

  // update the trace of the group
  if (trace_prec < 64)
  { // update _d
    comp_d *old_trace_d = *temp_gp_trace_d;
    *temp_gp_trace_d = (comp_d *)bmalloc((*num_temp_gps - 1) * sizeof(comp_d));

    for (j = 0; j < *num_temp_gps; j++)
      if (j < old_gp_number)
      { // simply copy over the trace
        set_d((*temp_gp_trace_d)[j], old_trace_d[j]);
      }
      else if (j == old_gp_number)
      { // add the traces of old and new groups to find the trace of the new group
        add_d((*temp_gp_trace_d)[new_gp_number], (*temp_gp_trace_d)[new_gp_number], old_trace_d[old_gp_number]);
      }
      else
      { // copy over the trace, reducing the index by 1 since we have eliminated 1 group
        set_d((*temp_gp_trace_d)[j-1], old_trace_d[j]);
      }

    // update num_temp_gps
    *num_temp_gps = *num_temp_gps - 1;
    // free old memory
    free(old_trace_d);
  }
  else
  { // update _mp
    comp_mp *old_trace_mp = *temp_gp_trace_mp;
    *temp_gp_trace_mp = (comp_mp *)bmalloc((*num_temp_gps - 1) * sizeof(comp_mp));
    for (j = 0; j < *num_temp_gps; j++)
    {
      if (j < old_gp_number)
      { // simply copy over the trace
        init_mp2((*temp_gp_trace_mp)[j], trace_prec);
        set_mp((*temp_gp_trace_mp)[j], old_trace_mp[j]);
      }
      else if (j == old_gp_number)
      { // add the traces of old and new groups to find the trace of the new group
        add_mp((*temp_gp_trace_mp)[new_gp_number], (*temp_gp_trace_mp)[new_gp_number], old_trace_mp[old_gp_number]);
      }
      else
      { // copy over the trace, reducing the index by 1 since we have eliminated 1 group
        init_mp2((*temp_gp_trace_mp)[j-1], trace_prec);
        set_mp((*temp_gp_trace_mp)[j-1], old_trace_mp[j]);
      }

      // clear old memory
      clear_mp(old_trace_mp[j]);
    }

    // update num_temp_gps
    *num_temp_gps = *num_temp_gps - 1;
    // free old memory
    free(old_trace_mp);
  }

  return;
}

void traceDecomposition(witness_t *W, int codim_index, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the decomposition using combinatorial trace test  *
\***************************************************************/
{
  int i, j, curr_size, size_complete, isComponent, num_paths = W->codim[codim_index].num_set, *curr_array = NULL;

  // check for group size of 1
  checkCompleteComponents(W, codim_index, tol, num_to_classify, isClassified, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);

  // so we can start at group size of 2
  curr_size = 2;
  // loop over all of the possible sizes 
  while (curr_size <= *num_temp_gps)
  { // setup curr_array to the proper size
    curr_array = (int *)bmalloc(curr_size * sizeof(int));
    // initialize curr_array
    for (i = 0; i < curr_size; i++)
      curr_array[i] = i;

    // loop over all possible combinations for this size - always have 0 <= curr_array[0] < .. < curr_array[curr_size - 1] <= *num_temp_gps - 1
    do
    { // check to see if curr_array describes a complete component
      isComponent = checkTraceArray(curr_size, curr_array, tol, *temp_gp_trace_d, *temp_gp_trace_mp, trace_prec);

      // determine if we have a component
      if (isComponent)
      { // curr_array describes a component - so we need to merge the temporary groups, clear the completed component and update curr_array

        // merge the temp groups curr_array[curr_size - 1], .., curr_array[1] into curr_array[0] - going in descending order so the group numbers do not change
        for (i = curr_size - 1; i > 0; i--)
          mergeTempGroups(curr_array[0], curr_array[i], num_paths, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);

        // clear the completed component
        checkCompleteComponents(W, codim_index, tol, num_to_classify, isClassified, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);

        // update curr_array
        if (curr_size > *num_temp_gps) // this should really only happen when *num_temp_gps == 0
        { // we are finished (completely finished classifying)
          size_complete = 1;
        }
        else if (curr_array[0] + curr_size > *num_temp_gps)
        { // we are finished with this size since there are no more groups of 'curr_size' that we have not already checked
          size_complete = 1;
        }
        else
        { // there are more possibilities that need to be checked
          size_complete = 0;

          // curr_array[i] = curr_array[0] + i, i = 1,..,curr_size - 1
          for (i = 1; i < curr_size; i++)
            curr_array[i] = curr_array[0] + i;
        }
      }
      else
      { // we do not have a complete component, so we update curr_array using the normal incrementing process

        // find the location in curr_array which we can increment, if possible
        for (i = curr_size - 1; i >= 0; i--)
          if (curr_array[i] < *num_temp_gps - curr_size + i)
          { // this location can be incremented - break out of for loop
            break;
          }

        // check to see if we can increment or we are finished
        if (i < 0)
        { // we have checked everything of this size
          size_complete = 1;
        }
        else
        { // we have more to check
          size_complete = 0;

          // increment curr_array[i]
          curr_array[i]++;
          // initialize curr_array[j] for j = i+1,..,curr_size-1
          for (j = i + 1; j < curr_size; j++)
            curr_array[j] = curr_array[i] + j - i;
        }
      }
    } while (!size_complete);

    // since we have tried every combination of this size, we increment the size and continue
    curr_size++;
  }

  // release memory
  free(curr_array);

  return;
}

int checkTraceArray(int curr_size, int *curr_array, double tol, comp_d *temp_gp_trace_d, comp_mp *temp_gp_trace_mp, int trace_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - not a complete component, 1 - otherwise    *
* NOTES: determines if curr_array describes a complete component*
\***************************************************************/
{
  int i, isComponent;

  if (trace_prec < 64)
  { // add using _d
    comp_d sum;
    set_zero_d(sum);

    for (i = 0; i < curr_size; i++)
    {
      add_d(sum, sum, temp_gp_trace_d[curr_array[i]]);
    }

    // determine if the sum is small enough
    if (d_abs_d(sum) < tol)
      isComponent = 1;
    else
      isComponent = 0;
  }
  else
  { // add using _mp
    comp_mp sum;
    init_mp2(sum, trace_prec);
    set_zero_mp(sum);

    for (i = 0; i < curr_size; i++)
    {
      add_mp(sum, sum, temp_gp_trace_mp[curr_array[i]]);
    }

    // determine if the sum is small enough
    if (d_abs_mp(sum) < tol)
      isComponent = 1;
    else
      isComponent = 0;

    clear_mp(sum);
  }

  return isComponent;
}

void copy_slice_moving_slice(membership_slice_moving_t *sliceMover, prog_t *fullRankProg, vec_d start_d, vec_mp start_mp, mpq_t **start_rat, vec_d target_d, vec_mp target_mp, mpq_t **target_rat, comp_d gamma_d, comp_mp gamma_mp, mpq_t *gamma_rat, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                *
* NOTES: setup sliceMover                                       *
\***************************************************************/
{
  int i, start_size, target_size;

  // setup sizes - should be equal!
  start_size = (MPType == 0 || MPType == 2 ? start_d->size : start_mp->size);
  target_size = (MPType == 0 || MPType == 2 ? target_d->size : target_mp->size);

  if (start_size != target_size)
  {
    printf("ERROR: The slice sizes are not the same!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // make sure that everything is initialized
  initialize_slice_moving_sliceVec(sliceMover, target_size, MPType);

  // setup prog
  sliceMover->Prog = fullRankProg;

  // setup slice vectors & gamma
  if (MPType == 0)
  { // copy _d
    vec_cp_d(sliceMover->startSliceVec_d, start_d);
    vec_cp_d(sliceMover->targetSliceVec_d, target_d);
    set_d(sliceMover->gamma_d, gamma_d);
  }
  else if (MPType == 1)
  { // setup _mp
    vec_cp_mp(sliceMover->startSliceVec_mp, start_mp);
    vec_cp_mp(sliceMover->targetSliceVec_mp, target_mp);
    set_mp(sliceMover->gamma_mp, gamma_mp);
  }
  else
  { // setup _d, _mp, _rat
    sliceMover->startSliceVec_d->size = sliceMover->startSliceVec_mp->size = start_size;
    sliceMover->targetSliceVec_d->size = sliceMover->targetSliceVec_mp->size = target_size;
    for (i = 0; i < target_size; i++)
    { // copy start
      set_rat(sliceMover->startSliceVec_rat[i], start_rat[i]);
      mpf_set_q(sliceMover->startSliceVec_mp->coord[i].r, start_rat[i][0]);
      mpf_set_q(sliceMover->startSliceVec_mp->coord[i].i, start_rat[i][1]);
      sliceMover->startSliceVec_d->coord[i].r = mpq_get_d(start_rat[i][0]);
      sliceMover->startSliceVec_d->coord[i].i = mpq_get_d(start_rat[i][1]);
      // copy target
      set_rat(sliceMover->targetSliceVec_rat[i], target_rat[i]);
      mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].r, target_rat[i][0]);
      mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].i, target_rat[i][1]);
      sliceMover->targetSliceVec_d->coord[i].r = mpq_get_d(target_rat[i][0]);
      sliceMover->targetSliceVec_d->coord[i].i = mpq_get_d(target_rat[i][1]);
    }
    // set gamma
    set_rat(sliceMover->gamma_rat, gamma_rat);
    mpf_set_q(sliceMover->gamma_mp->r, gamma_rat[0]);
    mpf_set_q(sliceMover->gamma_mp->i, gamma_rat[1]);
    sliceMover->gamma_d->r = mpq_get_d(gamma_rat[0]);
    sliceMover->gamma_d->i = mpq_get_d(gamma_rat[1]);
  }

  return;
}


