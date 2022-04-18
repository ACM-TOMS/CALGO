// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

int junkRemoval_membershipTest_seq(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT, char *midName);

int junkRemoval_membershipTest(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - the endpoint given by pathNum is a member  *
* of the pure_codim_index pure-dim witness set, 0 - otherwise   *
* NOTES: perform membership test                                *
\***************************************************************/
{
  int isMember = 0;

  isMember = junkRemoval_membershipTest_seq(W, pure_codim_index, pathNum, pathNum_codim_index, sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, T, OUT, midName);

  return isMember;
}

int junkRemoval_membershipTest_seq(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT, char *midName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - the endpoint given by pathNum is a member  *
* of the pure_codim_index pure-dim witness set, 0 - otherwise   *
* NOTES: perform membership test                                *
\***************************************************************/
{
  int i, retVal = 0, setupGamma = 1, badCount = 0, isMember = 0, num_paths = W->codim[pure_codim_index].num_set, maxBadCount = 10;
  FILE *MIDOUT = fopen(midName, "w");

  // setup the random slice in sliceMover for 'pathNum' endpoint 
  if (T->MPType == 0)
    setup_slice_moving_slice(sliceMover, W->codim[pathNum_codim_index].witnessPts_d[pathNum].endPt, NULL, 52, T->MPType, T->AMP_max_prec);
  else if (T->MPType == 1)
    setup_slice_moving_slice(sliceMover, NULL, W->codim[pathNum_codim_index].witnessPts_mp[pathNum].endPt, T->Precision, T->MPType, T->AMP_max_prec);
  else
    setup_slice_moving_slice(sliceMover, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_d, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_mp, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].curr_prec, T->MPType, T->AMP_max_prec);

  do
  { // we need to setup gamma
    setupGamma = 1;
    // initialize retVal
    retVal = 0;

    // loop through the witness points and move the slice to see if any of the endpoints hit the 'pathNum' endpoint
    // continue until we have checked all of them successfully OR we have hit the endpoint OR we have a tracking failure
    for (i = 0; i < num_paths && !isMember && !retVal; i++)
    { // verify that deflation worked
      if (fullRankProgInfo[i] != -1)
      { // finish the setup for sliceMover for this endpoint
        final_setup_slice_moving(sliceMover, fullRankProgs[i], T->MPType, T->AMP_max_prec, setupGamma);

        // we only setup gamma once for the whole set of witness points to avoid monodromy
        if (setupGamma)
          setupGamma = 0;

        // track the path and check for equality
        if (T->MPType == 0)
        { // compare using _d
          retVal = membership_slice_moving_track(&isMember, W->codim[pathNum_codim_index].witnessPts_d[pathNum].endPt, NULL, 52, sliceMover, endPts_d[i].endPt, NULL, 52, pathNum, T->final_tol_times_mult, T, OUT, MIDOUT);
        }
        else if (T->MPType == 1)
        { // compare using _mp
          retVal = membership_slice_moving_track(&isMember, NULL, W->codim[pathNum_codim_index].witnessPts_mp[pathNum].endPt, T->Precision, sliceMover, NULL, endPts_mp[i].endPt, T->Precision, pathNum, T->final_tol_times_mult, T, OUT, MIDOUT);
        }
        else
        { // compare using _amp
          retVal = membership_slice_moving_track(&isMember, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_d, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_mp, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].curr_prec, sliceMover, endPts_amp[i].endPt_d, endPts_amp[i].endPt_mp, endPts_amp[i].curr_prec, pathNum, T->final_tol_times_mult, T, OUT, MIDOUT);
        }
      }
    }

    if (isMember)
    { // we have hit the endpoint - exit loop
      break;
    }
    else if (retVal)
    { // we had bad random numbers - increment the number of bad loops and try again
      badCount++;

      // see if we had too many bad random numbers
      if (badCount >= maxBadCount)
      { // exit since we have failure - this should never happen
        printf("\nERROR: Failed to eliminate junk correctly (after trying %d different choices of random numbers).\n", maxBadCount);
        bexit(ERROR_LOOP_FAILURE);
      }
    }
  } while (retVal);

  // close MIDOUT
  fclose(MIDOUT);

  return isMember;
}

int membership_slice_moving_track(int *isMember, point_d testPt_d, point_mp testPt_mp, int testPt_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, double tol, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether tracking was a success or not          *
*  if success, checks to see if testPt and the endPt are same   *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;
  mpf_t error;
  endgame_data_t endPt;

  // initialize error
  mpf_init(error);

  // initialize endPt
  init_endgame_data(&endPt, T->Precision);
  // move the slice
  retVal = slice_moving_track(&endPt, sliceMover, startPt_d, startPt_mp, startPt_prec, pathNum_startPt, 0, T, OUT, MIDOUT);

  // check for membership if the slice moving was successful
  if (!retVal)
  { // find the difference between endpoint of the path and the given endPt - automatically uses the smallest size
    findDiff_point(error, endPt.PD_d.point, endPt.PD_mp.point, endPt.prec, testPt_d, testPt_mp, testPt_prec);
    if (mpf_get_d(error) < tol)
    { // they are the same!
      *isMember = 1;
    }
  }

  // clear error
  mpf_clear(error);
  // clear endPt
  clear_endgame_data(&endPt);

  return retVal;
}

int slice_moving_track(endgame_data_t *endPt, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, int checkSharpen, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether slice moving was a success or not      *
* NOTES: moves the slice                                        *
\***************************************************************/
{
  int retVal = 0;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  ptr_to_eval_d = &slice_mover_eval_d;
  ptr_to_eval_mp = &slice_mover_eval_mp;
  change_prec = &slice_mover_change_prec;
  find_dehom = &zero_dehom;

  // setup for tracking
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;
    init_point_data_d(&startPt, startPt_d->size);

    // setup for tracking the path
    point_cp_d(startPt.point, startPt_d);
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, sliceMover, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(pathNum_startPt, endPt, &startPt, OUT, MIDOUT, T, sliceMover, sliceMover, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // see if we should look to sharpen this endPt
    if (checkSharpen && endPt->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endPt, T, OUT, sliceMover, sliceMover, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer for the path and find retVal
    retVal = printMembershipFooter_d(&endPt->PD_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, endPt->retVal, T);

    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;
    init_point_data_mp(&startPt, startPt_mp->size);
  
    // setup for tracking the path
    point_cp_mp(startPt.point, startPt_mp);
    set_one_mp(startPt.time);
  
    // print the header for the path
    printPathHeader_mp(OUT, &startPt, T, pathNum_startPt, sliceMover, ptr_to_eval_mp);
  
    // track the path
    zero_dim_track_path_mp(pathNum_startPt, endPt, &startPt, OUT, MIDOUT, T, sliceMover, ptr_to_eval_mp, find_dehom);
    
    // see if we should look to sharpen this endPt
    if (checkSharpen && endPt->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endPt, T, OUT, sliceMover, sliceMover, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer in multi precision and find retVal
    retVal = printMembershipFooter_mp(&endPt->PD_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, endPt->retVal, T);

    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    point_data_d startPt;
    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    if (startPt_prec < 64)
    {
      point_cp_d(startPt.point, startPt_d);
    }
    else
    {
      point_mp_to_d(startPt.point, startPt_mp);
    }
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, sliceMover, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(pathNum_startPt, endPt, &startPt, OUT, MIDOUT, T, sliceMover, sliceMover, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // see if we should look to sharpen this endPt
    if (checkSharpen && endPt->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endPt, T, OUT, sliceMover, sliceMover, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer for the path and find retVal
    if (endPt->prec < 64)
    { // print the footer in double precision and find retVal
      retVal = printMembershipFooter_d(&endPt->PD_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, endPt->retVal, T);
    }
    else
    { // print the footer in multi precision and find retVal
      retVal = printMembershipFooter_mp(&endPt->PD_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, endPt->retVal, T);
    }

    clear_point_data_d(&startPt);
  }

  return retVal;
}

void initialize_slice_moving_sliceVec(membership_slice_moving_t *sliceMover, int size, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize startSliceVec & targetSliceVec in sliceMover*
\***************************************************************/
{ // setup startSliceVec
  if (sliceMover->startSliceVec_init)
  { // make sure that it is set to the correct size
    if (MPType == 0 && size != sliceMover->startSliceVec_d->size)
    { // reset the size
      change_size_vec_d(sliceMover->startSliceVec_d, size);
    }
    else if (MPType == 1 && size != sliceMover->startSliceVec_mp->size)
    { // reset the size
      change_size_vec_mp(sliceMover->startSliceVec_mp, size);
    }
    else if (MPType == 2 && size != sliceMover->startSliceVec_d->size)
    { // clear _rat
      clear_vec_rat(sliceMover->startSliceVec_rat, sliceMover->startSliceVec_d->size);

      // reset the size on _d, _mp & _rat
      change_size_vec_d(sliceMover->startSliceVec_d, size);
      change_size_vec_mp(sliceMover->startSliceVec_mp, size);
      init_vec_rat(sliceMover->startSliceVec_rat, size);
    }
  }
  else if (MPType == 0)
  { // initialize _d
    init_vec_d(sliceMover->startSliceVec_d, size);
  }
  else if (MPType == 1)
  { // initialize _mp
    init_vec_mp2(sliceMover->startSliceVec_mp, size, sliceMover->curr_precision);
  }
  else if (MPType == 2)
  { // initialize _d, _mp & _rat
    init_vec_d(sliceMover->startSliceVec_d, size);
    init_vec_mp2(sliceMover->startSliceVec_mp, size, sliceMover->curr_precision);
    init_vec_rat(sliceMover->startSliceVec_rat, size);
  }
  sliceMover->startSliceVec_init = 1;

  // setup targetSliceVec
  if (sliceMover->targetSliceVec_init)
  { // make sure that it is set to the correct size
    if (MPType == 0 && size != sliceMover->targetSliceVec_d->size)
    { // reset the size
      change_size_vec_d(sliceMover->targetSliceVec_d, size);
    }
    else if (MPType == 1 && size != sliceMover->targetSliceVec_mp->size)
    { // reset the size
      change_size_vec_mp(sliceMover->targetSliceVec_mp, size);
    }
    else if (MPType == 2 && size != sliceMover->targetSliceVec_d->size)
    { // clear _rat
      clear_vec_rat(sliceMover->targetSliceVec_rat, sliceMover->targetSliceVec_d->size);

      // reset the size on _d, _mp & _rat
      change_size_vec_d(sliceMover->targetSliceVec_d, size); 
      change_size_vec_mp(sliceMover->targetSliceVec_mp, size);
      init_vec_rat(sliceMover->targetSliceVec_rat, size);
    }
  }
  else if (MPType == 0)
  { // initialize _d
    init_vec_d(sliceMover->targetSliceVec_d, size);
  }
  else if (MPType == 1)
  { // initialize _mp
    init_vec_mp2(sliceMover->targetSliceVec_mp, size, sliceMover->curr_precision);
  }
  else if (MPType == 2)
  { // initialize _d, _mp & _rat
    init_vec_d(sliceMover->targetSliceVec_d, size);
    init_vec_mp2(sliceMover->targetSliceVec_mp, size, sliceMover->curr_precision);
    init_vec_rat(sliceMover->targetSliceVec_rat, size);
  }
  sliceMover->targetSliceVec_init = 1;

  return;
}

void setup_slice_moving_slice(membership_slice_moving_t *sliceMover, point_d endPt_d, point_mp endPt_mp, int endPt_prec, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the basic elements in sliceMover                 *
\***************************************************************/
{
  int i, j, size;

  if (MPType == 0)
  { // find the size
    size = sliceMover->B_d->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // setup startSliceVec_d == 0
    change_size_vec_d(sliceMover->startSliceVec_d, size);
    sliceMover->startSliceVec_d->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
    }

    // setup targetSliceVec_d
    change_size_vec_d(sliceMover->targetSliceVec_d, size);
    sliceMover->targetSliceVec_d->size = size;
    for (i = 0; i < size; i++)
    { // B * endPt
      set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
      for (j = 0; j < sliceMover->B_d->cols; j++)
      {
        sum_mul_d(&sliceMover->targetSliceVec_d->coord[i], &sliceMover->B_d->entry[i][j], &endPt_d->coord[j]);
      }      
    }
  }
  else if (MPType == 1)
  { // find the size
    size = sliceMover->B_mp->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // setup startSliceVec_mp == 0
    change_size_vec_mp(sliceMover->startSliceVec_mp, size);
    sliceMover->startSliceVec_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
    }

    // setup targetSliceVec_mp
    change_size_vec_mp(sliceMover->targetSliceVec_mp, size);
    sliceMover->targetSliceVec_mp->size = size;
    for (i = 0; i < size; i++)
    { // B * endPt
      set_zero_mp(&sliceMover->targetSliceVec_mp->coord[i]);
      for (j = 0; j < sliceMover->B_mp->cols; j++)
      {
        sum_mul_mp(&sliceMover->targetSliceVec_mp->coord[i], &sliceMover->B_mp->entry[i][j], &endPt_mp->coord[j]);
      }
    }
  }
  else // MPType == 2
  { // find the size
    size = sliceMover->B_d->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // setup startSliceVec_d, _mp, _rat == 0
    change_size_vec_d(sliceMover->startSliceVec_d, size);
    change_size_vec_mp(sliceMover->startSliceVec_mp, size);
    sliceMover->startSliceVec_d->size = sliceMover->startSliceVec_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
      set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
      set_zero_rat(sliceMover->startSliceVec_rat[i]);
    }

    // setup targetSliceVec_d, _mp, _rat
    change_size_vec_d(sliceMover->targetSliceVec_d, size);
    change_size_vec_mp(sliceMover->targetSliceVec_mp, size);
    sliceMover->targetSliceVec_d->size = sliceMover->targetSliceVec_mp->size = size;
    if (endPt_prec < 64)
    { // use endPt_d
      for (i = 0; i < size; i++)
      { // B * endPt
        set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
        for (j = 0; j < sliceMover->B_d->cols; j++)
        {
          sum_mul_d(&sliceMover->targetSliceVec_d->coord[i], &sliceMover->B_d->entry[i][j], &endPt_d->coord[j]);
        }

        // setup _d, _mp & _rat
        mpq_set_d(sliceMover->targetSliceVec_rat[i][0], sliceMover->targetSliceVec_d->coord[i].r);
        mpq_set_d(sliceMover->targetSliceVec_rat[i][1], sliceMover->targetSliceVec_d->coord[i].i);
        mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].r, sliceMover->targetSliceVec_rat[i][0]);
        mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].i, sliceMover->targetSliceVec_rat[i][1]);
        sliceMover->targetSliceVec_d->coord[i].r = mpq_get_d(sliceMover->targetSliceVec_rat[i][0]);
        sliceMover->targetSliceVec_d->coord[i].i = mpq_get_d(sliceMover->targetSliceVec_rat[i][1]);
      }
    }
    else
    {
      comp_mp tempComp, tempSum;
      init_mp2(tempComp, endPt_prec);
      init_mp2(tempSum, endPt_prec);

      for (i = 0; i < size; i++)
      { // B * endPt
        set_zero_mp(tempSum);
        for (j = 0; j < sliceMover->B_mp->cols; j++)
        {
          mpf_set_q(tempComp->r, sliceMover->B_rat[i][j][0]);
          mpf_set_q(tempComp->i, sliceMover->B_rat[i][j][1]);

          sum_mul_mp(tempSum, tempComp, &endPt_mp->coord[j]);
        }

        // setup _d, _mp & _rat
        mpf_t_to_rat(sliceMover->targetSliceVec_rat[i][0], tempSum->r);
        mpf_t_to_rat(sliceMover->targetSliceVec_rat[i][1], tempSum->i);
        mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].r, sliceMover->targetSliceVec_rat[i][0]);
        mpf_set_q(sliceMover->targetSliceVec_mp->coord[i].i, sliceMover->targetSliceVec_rat[i][1]);
        sliceMover->targetSliceVec_d->coord[i].r = mpq_get_d(sliceMover->targetSliceVec_rat[i][0]);
        sliceMover->targetSliceVec_d->coord[i].i = mpq_get_d(sliceMover->targetSliceVec_rat[i][1]);
      }

      clear_mp(tempComp);
      clear_mp(tempSum);
    }
  }

  return;
}

void basic_setup_slice_moving(membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the basic elements in sliceMover for codim_index *
\***************************************************************/
{
  int i, j;

  sliceMover->curr_codim = W->codim[codim_index].codim;
  sliceMover->orig_variables = W->orig_variables;
  sliceMover->curr_precision = W->curr_precision;

  // initialize A even though it will not be setup here
  if (MPType == 0 || MPType == 2)
  {
    init_mat_d(sliceMover->A_d, 0, 0);
  }
  if (MPType == 1 || MPType == 2)
  {
    init_mat_mp2(sliceMover->A_mp, 0, 0, sliceMover->curr_precision);
  }
  sliceMover->A_rat = NULL;
  sliceMover->A_d->rows = sliceMover->A_d->cols = sliceMover->A_mp->rows = sliceMover->A_mp->cols = 0;
  sliceMover->startSliceVec_init = 0;
  sliceMover->targetSliceVec_init = 0;

  // setup gamma
  if (MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(sliceMover->gamma_d);
  }
  else if (MPType == 1)
  { // setup gamma_mp
    init_mp(sliceMover->gamma_mp);
    get_comp_rand_mp(sliceMover->gamma_mp);
  }
  else
  { // setup gamma_d, gamma_mp, gamma_rat
    get_comp_rand_rat(sliceMover->gamma_d, sliceMover->gamma_mp, sliceMover->gamma_rat, sliceMover->curr_precision, max_prec, 1, 1);
  }

  // copy B
  if (MPType == 0 || MPType == 2)
  { // setup _d
    init_mat_d(sliceMover->B_d, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    mat_cp_d(sliceMover->B_d, W->codim[codim_index].B_d);
  }
  if (MPType == 1 || MPType == 2)
  { // setup _mp
    init_mat_mp2(sliceMover->B_mp, W->codim[codim_index].B_mp->rows, W->codim[codim_index].B_mp->cols, sliceMover->curr_precision);
    mat_cp_mp(sliceMover->B_mp, W->codim[codim_index].B_mp);
  }
  if (MPType == 2)
  { // setup _amp
    init_mat_rat(sliceMover->B_rat, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    mat_cp_rat(sliceMover->B_rat, W->codim[codim_index].B_rat, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
  }

  // setup p
  if (MPType == 0 || MPType == 2)
  { // setup _d
    init_vec_d(sliceMover->p_d, W->codim[codim_index].p_d->size);
    vec_cp_d(sliceMover->p_d, W->codim[codim_index].p_d);
  }
  if (MPType == 1 || MPType == 2)
  { // setup _mp
    init_vec_mp2(sliceMover->p_mp, W->codim[codim_index].p_mp->size, sliceMover->curr_precision);
    vec_cp_mp(sliceMover->p_mp, W->codim[codim_index].p_mp);
  }
  if (MPType == 2)
  { // setup _rat
    init_vec_rat(sliceMover->p_rat, W->codim[codim_index].p_d->size);
    vec_cp_rat(sliceMover->p_rat, W->codim[codim_index].p_rat, W->codim[codim_index].p_d->size);
  }

  // setup K_rat
  if (MPType == 0)
  { // setup K_rat
    double tol_pivot = 1e-14, tol_sign = 1e-18, largeChange = 1e13;
    mat_d tempMat, Q, R, P;

    init_mat_d(tempMat, 0, 0); init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

    // setup tempMat to be [[B];[p]];
    mat_cp_d(tempMat, sliceMover->B_d);

    // add a row to tempMat
    i = tempMat->rows + 1;
    increase_rows_mat_d(tempMat, i);
    tempMat->rows = i;

    // setup this bottom row
    i = sliceMover->B_d->rows;
    for (j = 0; j < tempMat->cols; j++)
    {
      set_d(&tempMat->entry[i][j], &sliceMover->p_d->coord[j]);
    }

    // transpose tempMat
    transpose_d(tempMat, tempMat);

    // perform a QR-decomposition on tempMat
    QR_d(Q, R, P, tempMat, tol_pivot, tol_sign, largeChange, 0);

    // setup K using Q
    sliceMover->K_rows = Q->rows;
    sliceMover->K_cols = Q->cols - tempMat->cols;
    init_mat_rat(sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols);

    for (i = 0; i < sliceMover->K_rows; i++)
      for (j = 0; j < sliceMover->K_cols; j++)
      {
        mpq_set_d(sliceMover->K_rat[i][j][0], Q->entry[i][tempMat->cols+j].r);
        mpq_set_d(sliceMover->K_rat[i][j][1], Q->entry[i][tempMat->cols+j].i);
      }
  
    // clear mat
    clear_mat_d(tempMat); clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
  }
  else if (MPType == 1)
  { // setup K_rat
    mat_mp tempMat, Q, R, P;
 
    init_mat_mp(tempMat, 0, 0); init_mat_mp(Q, 0, 0); init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0); 

    // setup tempMat to be [[B];[p]];
    mat_cp_mp(tempMat, sliceMover->B_mp);

    // add a row to tempMat
    i = tempMat->rows + 1;
    increase_rows_mat_mp(tempMat, i);
    tempMat->rows = i;

    // setup this bottom row 
    i = sliceMover->B_mp->rows;
    for (j = 0; j < tempMat->cols; j++)
    {
      set_mp(&tempMat->entry[i][j], &sliceMover->p_mp->coord[j]);
    }

    // transpose tempMat
    transpose_mp(tempMat, tempMat);

    // perform a QR-decomposition on tempMat
    QR_mp_prec(Q, R, P, tempMat, 0, sliceMover->curr_precision);

    // setup K using Q
    sliceMover->K_rows = Q->rows;
    sliceMover->K_cols = Q->cols - tempMat->cols;
    init_mat_rat(sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols);

    for (i = 0; i < sliceMover->K_rows; i++)
      for (j = 0; j < sliceMover->K_cols; j++)
      {
        mpf_t_to_rat(sliceMover->K_rat[i][j][0], Q->entry[i][j+tempMat->cols].r);
        mpf_t_to_rat(sliceMover->K_rat[i][j][1], Q->entry[i][j+tempMat->cols].i);
      }

    // clear mat
    clear_mat_mp(tempMat); clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);
  }
  else
  { // setup K_rat
    mat_mp tempMat, Q, R, P;

    init_mat_mp2(tempMat, 0, 0, max_prec); init_mat_mp2(Q, 0, 0, max_prec); init_mat_mp2(R, 0, 0, max_prec); init_mat_mp2(P, 0, 0, max_prec);

    // set precision to max_prec
    initMP(max_prec);

    change_size_mat_mp(tempMat, sliceMover->B_d->rows + 1, sliceMover->B_d->cols);
    tempMat->rows = sliceMover->B_d->rows + 1;
    tempMat->cols = sliceMover->B_d->cols;
    // setup tempMat to be [[B];[p]];
    for (j = 0; j < tempMat->cols; j++)
    {
      for (i = 0; i < sliceMover->B_d->rows; i++)
      {
        mpf_set_q(tempMat->entry[i][j].r, sliceMover->B_rat[i][j][0]);
        mpf_set_q(tempMat->entry[i][j].i, sliceMover->B_rat[i][j][1]);
      }
      mpf_set_q(tempMat->entry[i][j].r, sliceMover->p_rat[j][0]);
      mpf_set_q(tempMat->entry[i][j].i, sliceMover->p_rat[j][1]);
    }

    // transpose tempMat
    transpose_mp(tempMat, tempMat);

    // perform a QR-decomposition on tempMat
    QR_mp_prec(Q, R, P, tempMat, 0, sliceMover->curr_precision);

    // set precison back to curr_precision
    initMP(sliceMover->curr_precision);

    // setup K using Q
    sliceMover->K_rows = Q->rows;
    sliceMover->K_cols = Q->cols - tempMat->cols;
    init_mat_rat(sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols);

    for (i = 0; i < sliceMover->K_rows; i++)
      for (j = 0; j < sliceMover->K_cols; j++)
      {
        mpf_t_to_rat(sliceMover->K_rat[i][j][0], Q->entry[i][j+tempMat->cols].r);
        mpf_t_to_rat(sliceMover->K_rat[i][j][1], Q->entry[i][j+tempMat->cols].i);
      }

    // clear mat
    clear_mat_mp(tempMat); clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);
  }

  return;
}

void final_setup_slice_moving(membership_slice_moving_t *sliceMover, prog_t *fullRankProg, int MPType, int max_prec, int setupGamma)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup Prog & A                                         *
\***************************************************************/
{
  int A_rows, A_cols;

  sliceMover->Prog = fullRankProg;

  if (MPType == 0)
  { // setup A_d
    A_rows = sliceMover->Prog->numVars - sliceMover->B_d->rows - 1;
    A_cols = sliceMover->Prog->numFuncs - A_rows;

    if (A_rows < 0)
    {
      printf("ERROR: Incorrect number of variables and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }
    if (A_cols < 0)
    {
      printf("ERROR: Incorrect number of functions, variables, and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }

    make_matrix_random_d(sliceMover->A_d, A_rows, A_cols);

    if (setupGamma)
    { // setup gamma_d
      get_comp_rand_d(sliceMover->gamma_d);
    }
  }
  else if (MPType == 1)
  { // setup A_mp
    A_rows = sliceMover->Prog->numVars - sliceMover->B_mp->rows - 1;
    A_cols = sliceMover->Prog->numFuncs - A_rows;

    if (A_rows < 0)
    {
      printf("ERROR: Incorrect number of variables and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }
    if (A_cols < 0)
    {
      printf("ERROR: Incorrect number of functions, variables, and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }

    make_matrix_random_mp(sliceMover->A_mp, A_rows, A_cols, sliceMover->curr_precision);

    if (setupGamma)
    { // setup gamma_mp
      get_comp_rand_mp(sliceMover->gamma_mp);
    }
  }
  else
  { // setup A
    A_rows = sliceMover->Prog->numVars - sliceMover->B_d->rows - 1;
    A_cols = sliceMover->Prog->numFuncs - A_rows;

    if (A_rows < 0)
    {
      printf("ERROR: Incorrect number of variables and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }
    if (A_cols < 0)
    {
      printf("ERROR: Incorrect number of functions, variables, and slices!\n");
      bexit(ERROR_CONFIGURATION);
    }

    setprec_mat_mp(sliceMover->A_mp, sliceMover->curr_precision);
    if (sliceMover->A_rat != NULL)
    {
      clear_mat_rat(sliceMover->A_rat, sliceMover->A_d->rows, sliceMover->A_d->cols);
    }
    init_mat_rat(sliceMover->A_rat, A_rows, A_cols);

    make_matrix_random_rat(sliceMover->A_d, sliceMover->A_mp, sliceMover->A_rat, A_rows, A_cols, sliceMover->curr_precision, max_prec, 0, 0);

    if (setupGamma)
    { // setup gamma
      setprec_mp(sliceMover->gamma_mp, sliceMover->curr_precision);
      get_comp_rand_rat(sliceMover->gamma_d, sliceMover->gamma_mp, sliceMover->gamma_rat, sliceMover->curr_precision, max_prec, 0, 0);
    }
  }

  return;
}

int slice_mover_change_prec(void const *ED, int prec) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                *
* NOTES: change precision for slice moving                      *
\***************************************************************/
{
  int i, j;
  membership_slice_moving_t *SM = (membership_slice_moving_t *)ED;

  // set the SLP to the correct precision
  SM->Prog->precision = prec;

  // check to see if we need to change the precision
  if (prec != SM->curr_precision)
  { // change precision on gamma
    setprec_mp(SM->gamma_mp, prec);
    mpf_set_q(SM->gamma_mp->r, SM->gamma_rat[0]);
    mpf_set_q(SM->gamma_mp->i, SM->gamma_rat[1]);

    // update A
    setprec_mat_mp(SM->A_mp, prec);
    for (i = 0; i < SM->A_mp->rows; i++)
      for (j = 0; j < SM->A_mp->cols; j++)
      {
        mpf_set_q(SM->A_mp->entry[i][j].r, SM->A_rat[i][j][0]);
        mpf_set_q(SM->A_mp->entry[i][j].i, SM->A_rat[i][j][1]); 
      }

    // update B
    setprec_mat_mp(SM->B_mp, prec);
    for (i = 0; i < SM->B_mp->rows; i++)
      for (j = 0; j < SM->B_mp->cols; j++)
      {
        mpf_set_q(SM->B_mp->entry[i][j].r, SM->B_rat[i][j][0]);
        mpf_set_q(SM->B_mp->entry[i][j].i, SM->B_rat[i][j][1]);
      }

    // update p
    setprec_vec_mp(SM->p_mp, prec);
    for (i = 0; i < SM->p_mp->size; i++)
    {
      mpf_set_q(SM->p_mp->coord[i].r, SM->p_rat[i][0]);
      mpf_set_q(SM->p_mp->coord[i].i, SM->p_rat[i][1]);
    }

    // update startSliceVec, if needed
    if (SM->startSliceVec_init)
    {
      setprec_vec_mp(SM->startSliceVec_mp, prec);
      for (i = 0; i < SM->startSliceVec_mp->size; i++)
      {
        mpf_set_q(SM->startSliceVec_mp->coord[i].r, SM->startSliceVec_rat[i][0]);
        mpf_set_q(SM->startSliceVec_mp->coord[i].i, SM->startSliceVec_rat[i][1]);
      }
    }

    // update targetSliceVec, if needed
    if (SM->targetSliceVec_init)
    {
      setprec_vec_mp(SM->targetSliceVec_mp, prec);
      for (i = 0; i < SM->targetSliceVec_mp->size; i++)
      {
        mpf_set_q(SM->targetSliceVec_mp->coord[i].r, SM->targetSliceVec_rat[i][0]);
        mpf_set_q(SM->targetSliceVec_mp->coord[i].i, SM->targetSliceVec_rat[i][1]);
      }
    }

    // set the current precision
    SM->curr_precision = prec;
  }

  return 0;
}

void clear_slice_mover(membership_slice_moving_t *SM, int MPType) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                *
* NOTES: clear slice moving                                     *
\***************************************************************/
{
  // set the SLP to the correct precision
  SM->Prog = NULL;

  // clear structures
  if (MPType == 2)
  {
    clear_rat(SM->gamma_rat);
    clear_mat_rat(SM->A_rat, SM->A_d->rows, SM->A_d->cols);
    clear_mat_rat(SM->B_rat, SM->B_d->rows, SM->B_d->cols);
    clear_vec_rat(SM->p_rat, SM->p_d->size);
    if (SM->startSliceVec_init)
    {
      clear_vec_rat(SM->startSliceVec_rat, SM->startSliceVec_d->size);
    }
    if (SM->targetSliceVec_init)
    {
      clear_vec_rat(SM->targetSliceVec_rat, SM->targetSliceVec_d->size);
    }
  }
  if (MPType == 0 || MPType == 2)
  {
    clear_d(SM->gamma_d);
    clear_mat_d(SM->A_d);
    clear_mat_d(SM->B_d);
    clear_vec_d(SM->p_d);
    if (SM->startSliceVec_init)
    {
      clear_vec_d(SM->startSliceVec_d);
    }
    if (SM->targetSliceVec_init)
    {
      clear_vec_d(SM->targetSliceVec_d);
    }
  }
  if (MPType == 1 || MPType == 2)
  {
    clear_mp(SM->gamma_mp);
    clear_mat_mp(SM->A_mp);
    clear_mat_mp(SM->B_mp);
    clear_vec_mp(SM->p_mp);
    if (SM->startSliceVec_init)
    {
      clear_vec_mp(SM->startSliceVec_mp);
    }
    if (SM->targetSliceVec_init)
    {
      clear_vec_mp(SM->targetSliceVec_mp);
    }
  }

  clear_mat_rat(SM->K_rat, SM->K_rows, SM->K_cols);

  return;
}

