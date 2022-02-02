// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"
#include "parallel.h"

void regen_pos_dim_TrackLevel(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void regen_pos_dim_TrackPath(regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void regen_pos_dim_TrackLevel_trackBack(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void regen_pos_dim_SortLevel(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *WITSUPER, int codim_index, int num_codim, FILE *ENDPTS, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int regen_pos_dim_MoveNextCodim(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, char *nextStartFile, int codim_index, int num_codim, point_d *movePts_d, point_mp *movePts_mp, int *Pts_prec, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int regen_pos_dim_LinearTrackPath(int path_num, point_data_d *startPt_d, point_data_mp *startPt_mp, int new_degree, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int determineRPDSoln_d(double tol, double ratio, point_d point, point_d last_point, comp_d time, regen_pos_dim_t *RPD, int codim_index);
int determineRPDSoln_mp(double tol, double ratio, point_mp point, point_mp last_point, comp_mp time, int prec, regen_pos_dim_t *RPD, int codim_index);
void printRPDRawOut(regen_pos_dim_t *RPD, FILE *RAWOUT, endgame_data_t *endPt, int corank, int retVal, double smallest_nonzero_SV, double largest_zero_SV);
void regenPDFindOrigVarsDehom_d(point_d orig_vars_d, point_d dehom_d, point_d P_d, regen_pos_dim_t *RPD, int codim_index);
void regenPDFindOrigVarsDehom_mp(point_mp orig_vars_mp, point_mp dehom_mp, point_mp P_mp, regen_pos_dim_t *RPD, int codim_index);

int regen_pos_dim_main(witness_t *witnessSuperset, int startLevel, int maxCodim, int specificCodim, tracker_config_t *T, int pathMod, double midpoint_tol, double intrinsicCutoffMultiplier, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if specificCodim is top dimension              *
* NOTES: does regeneration to find witness superset             *
\***************************************************************/
{
  int startCodimIndex = 0, retVal = 0;
  regen_pos_dim_t RPD;
  trackingStats trackCount;
  char startName[] = "startRPD", witName[] = "witness_superset";

  // initialize trackCount
  init_trackingStats(&trackCount);

  // setup RPD
  regen_pos_dim_setup(startLevel, &maxCodim, specificCodim, T, &RPD, "preproc_data", "deg.out", startName, intrinsicCutoffMultiplier);

  // setup startCodimIndex
  startCodimIndex = startLevel <= 1 ? 0 : startLevel - 1;

  // setup the rest of T, if needed
  if (T->MPType == 1)
  { // initialize latest_newton_residual_mp
    mpf_init(T->latest_newton_residual_mp);
  }
  else if (T->MPType == 2)
  { // initialize latest_newton_residual_mp
    mpf_init2(T->latest_newton_residual_mp, T->AMP_max_prec);

    // setup eps, Phi & Psi
    T->AMP_eps = (double) T->numVars * T->numVars;
    T->AMP_Phi = T->AMP_bound_on_degree * (T->AMP_bound_on_degree - 1) * T->AMP_bound_on_abs_vals_of_coeffs;
    T->AMP_Psi = T->AMP_bound_on_degree * T->AMP_bound_on_abs_vals_of_coeffs;
  }

  // ready to do the tracking - go to either parallel or serial tracking
  if (num_processes > 1)
  { // do parallel tracking using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
    worker_info sendType;
    sendType.dataType = REGEN_POS_DIM;
    bcast_worker_info(&sendType, my_id, headnode);

    // do parallel tracking
    regen_pos_dim_par_track(startCodimIndex, maxCodim, &trackCount, pathMod, midpoint_tol, T, &RPD, startName, witName, my_id, num_processes, headnode);
#endif
  }
  else
  { // do sequential tracking
    regen_pos_dim_seq_track(startCodimIndex, maxCodim, &trackCount, pathMod, midpoint_tol, T, &RPD, startName, witName);
  }

  // print output chart
  regen_pos_dim_OutputChart(&RPD, stdout, maxCodim);

  // setup witnessSuperset & clear RPD
  retVal = regen_pos_dim_copyWitness_clear(witnessSuperset, &RPD, witName, T->MPType, T->AMP_max_prec, specificCodim);

  return retVal;
}

void regen_pos_dim_seq_track(int startCodimIndex, int maxCodim, trackingStats *trackCount, int pathMod, double midpoint_tol, tracker_config_t *T, regen_pos_dim_t *RPD, char *startName, char *witName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does regen pos dim tracking in serial                  *
\***************************************************************/
{
  int codim_index, num_paths, num_crossings = 0;
  size_t size;
  char *str = NULL;
  FILE *MIDOUT = NULL, *RAWOUT = NULL, *START = NULL, *FAIL = NULL, *OUT = NULL;
  char outName[] = "output", midName[] = "midpath_data", failName[] = "fail", rawTrackFile[] = "rawout_track", rawSortFile[] = "rawout_sort", rawPrepareFile[] = "rawout_prepare";

  // open OUT & FAIL
  OUT = fopen(outName, "w");
  FAIL = fopen(failName, "w");

  // find the witness supersets
  for (codim_index = startCodimIndex; codim_index < maxCodim; codim_index++)
  { // find the number of paths that are tracking on this codim
    num_paths = RPD->codim[codim_index].num_paths;

    // open MIDOUT
    MIDOUT = fopen(midName, "w");

    // setup RAWOUT to track this level
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, RPD->codim[codim_index].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, RPD->codim[codim_index].codim);
    RAWOUT = fopen(str, "w");

    // setup START
    size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[codim_index].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", startName, RPD->codim[codim_index].codim);
    START = fopen(str, "r");
    // make sure START exits
    if (START == NULL)
    {
      printf("ERROR: The file to contain the start points, '%s', does not exist!\n", str);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // read in the number of paths from START
    fscanf(START, "%d", &num_crossings);
    // make sure that we have agreement
    if (num_paths != num_crossings)
    {
      printf("ERROR: The number of paths (%d vs %d) described in '%s' is not correct!\n", num_paths, num_crossings, str);
      bexit(ERROR_INVALID_SIZE);
    }

    // track the codimension pointed to by codim_index
    if (T->endgameNumber == 3)
    { // use the trackBack endgame
      regen_pos_dim_TrackLevel_trackBack(pathMod, RPD, T, START, OUT, RAWOUT, MIDOUT, FAIL, codim_index, maxCodim, trackCount, regen_pos_dim_eval_d, regen_pos_dim_eval_mp, change_regen_pos_dim_prec, regen_pos_dim_dehom);
    }
    else
    { // use the standard endgame
      regen_pos_dim_TrackLevel(pathMod, RPD, T, START, OUT, RAWOUT, MIDOUT, FAIL, codim_index, maxCodim, trackCount, regen_pos_dim_eval_d, regen_pos_dim_eval_mp, change_regen_pos_dim_prec, regen_pos_dim_dehom);
    }

    // close START, RAWOUT & MIDOUT
    fclose(START);
    fclose(RAWOUT);
    fclose(MIDOUT);

    // check for path crossings
    num_crossings = 0;

    // check to see if using intrinsic slice
    if (RPD->codim[codim_index].useIntrinsicSlice)
      midpoint_checker(num_paths, RPD->codim[codim_index].codim, midpoint_tol, &num_crossings);
    else
      midpoint_checker(num_paths, RPD->new_variables, midpoint_tol, &num_crossings);

    // print message to screen about path crossing
    if (num_crossings > 0)
      printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

    // reopen RAWOUT for sorting - new start points
    size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, RPD->codim[codim_index].codim); 
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawSortFile, RPD->codim[codim_index].codim);
    RAWOUT = fopen(str, "w");

    // reopen MIDOUT for sorting - witness superset points
    size = 1 + snprintf(NULL, 0, "%s_%d", witName, RPD->codim[codim_index].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", witName, RPD->codim[codim_index].codim);
    MIDOUT = fopen(str, "w");

    // setup START to be the file containing the endpoints
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, RPD->codim[codim_index].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, RPD->codim[codim_index].codim);
    START = fopen(str, "r");
    // make sure START exits
    if (START == NULL)
    {
      printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // sort the endpoints
    regen_pos_dim_SortLevel(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, codim_index, maxCodim, START, regen_pos_dim_eval_d, regen_pos_dim_eval_mp, change_regen_pos_dim_prec);

    // close START, RAWOUT & MIDOUT
    fclose(START);
    fclose(RAWOUT);
    fclose(MIDOUT);

    // print the regen summary
    RAWOUT = fopen("regenSummary", "w");
    printRPDSummaryData(RPD, codim_index, RAWOUT);
    fclose(RAWOUT);

    // see if we need to prepare the next codim
    if (codim_index + 1 < maxCodim)
    { // reopen all of the files since we have another codim to track
      MIDOUT = fopen(midName, "w");

      // setup START to be the file containing the ones that need moved forward
      size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, RPD->codim[codim_index].codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawSortFile, RPD->codim[codim_index].codim);
      START = fopen(str, "r");
      // make sure START exits
      if (START == NULL)
      {
        printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // reopen RAWOUT for preparing
      size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, RPD->codim[codim_index + 1].codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawPrepareFile, RPD->codim[codim_index + 1].codim);
      RAWOUT = fopen(str, "w");

      // setup the name of the file to contain the next set of start points
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[codim_index + 1].codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, RPD->codim[codim_index + 1].codim);

      // create the next set of start points for the next codim
      num_paths = regen_pos_dim_PrepareNextCodim(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, FAIL, START, str, codim_index, maxCodim, trackCount, regen_pos_dim_moving_linear_eval_d, regen_pos_dim_moving_linear_eval_mp, change_regen_pos_dim_prec);

      // close RAWOUT, MIDOUT, START
      fclose(RAWOUT);
      fclose(MIDOUT);
      fclose(START);

      // check for path crossings
      num_crossings = 0;

      // check to see if using intrinsic slice
      if (RPD->codim[codim_index].useIntrinsicSlice)
        midpoint_checker(num_paths, RPD->codim[codim_index].codim + 1, midpoint_tol, &num_crossings);
      else
        midpoint_checker(num_paths, RPD->new_variables, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }
    else
    { // print the bottom codim file
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[codim_index].codim + 1);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, RPD->codim[codim_index].codim + 1);
      RAWOUT = fopen(str, "w");
      printRPDRelevantData(RPD, T->MPType, codim_index + 1, RAWOUT);
      fclose(RAWOUT);
    }
  }

  // close OUT & FAIL
  fclose(OUT);
  fclose(FAIL);

  return;
}

void regen_pos_dim_TrackPath_found(int indexJ, trackBack_samples_t *EGsample, int corank, double smallest, double largest, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the path                                         *
\***************************************************************/
{
  int j, size, old_pathNum = EGsample->endPt.pathNum;

  // setup pathNum
  EGsample->endPt.pathNum = path_num;

  // seutp for tracking the path
  RPD->curr_codim = codim_index;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // print the mid point
    fprintf(MIDOUT, "%d\n", path_num);
    size = EGsample->midPt_d[indexJ]->size;
    for (j = 0; j < size; j++)
      fprintf(MIDOUT, "%.15e %.15e\n", EGsample->midPt_d[indexJ]->coord[j].r, EGsample->midPt_d[indexJ]->coord[j].i);
    fprintf(MIDOUT, "\n");
  }
  else if (T->MPType == 1)
  { // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, RPD, ptr_to_eval_mp);

    // print the mid point
    fprintf(MIDOUT, "%d\n", path_num);
    size = EGsample->midPt_mp[indexJ]->size;
    for (j = 0; j < size; j++)
    {
      print_mp(MIDOUT, 0, &EGsample->midPt_mp[indexJ]->coord[j]);
      fprintf(MIDOUT, "\n");
    }
    fprintf(MIDOUT, "\n");
  }
  else
  { // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // print the mid point to MIDOUT
    fprintf(MIDOUT, "%d\n", path_num);
    if (EGsample->midPt_prec < 64)
    { // print midPt_d
      size = EGsample->midPt_d[indexJ]->size;
      for (j = 0; j < size; j++)
        fprintf(MIDOUT, "%.15e %.15e\n", EGsample->midPt_d[indexJ]->coord[j].r, EGsample->midPt_d[indexJ]->coord[j].i);
      fprintf(MIDOUT, "\n");
    }
    else
    { // print midPt_mp
      size = EGsample->midPt_mp[indexJ]->size;
      for (j = 0; j < size; j++)
      {
        print_mp(MIDOUT, 0, &EGsample->midPt_mp[indexJ]->coord[j]);
        fprintf(MIDOUT, "\n");
      }
      fprintf(MIDOUT, "\n");
    }
  }

  // print the footer to OUT and print all info to RAWOUT
  printRPDFooter(RPD, codim_index, &EGsample->endPt, corank, OUT, RAWOUT, FAIL, T, trackCount, smallest, largest);

  // reset pathNum
  EGsample->endPt.pathNum = old_pathNum;

  return;
}

void regen_pos_dim_TrackPath_trackBack(double trackBack_final_tol, trackBack_samples_t *EGsample, int *corank, double *smallest, double *largest, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path                                        *
\***************************************************************/
{
  int rankType = 1;

  // seutp for tracking the path
  RPD->curr_codim = codim_index;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision - find corank

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_rank_d(path_num, rankType, NULL, corank, smallest, largest, trackBack_final_tol, EGsample, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision - find corank

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, RPD, ptr_to_eval_mp);

    // track the path
    zero_dim_trackBack_path_rank_mp(path_num, rankType, NULL, corank, smallest, largest, trackBack_final_tol, EGsample, startPt_mp, OUT, MIDOUT, T, RPD, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP - find corank

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_rank_d(path_num, rankType, NULL, corank, smallest, largest, trackBack_final_tol, EGsample, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // check to see if it should be sharpened
  if (EGsample->endPt.retVal == 0 && T->sharpenDigits > 0 && *corank == 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&EGsample->endPt, T, OUT, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // print the footer to OUT and print all info to RAWOUT
  printRPDFooter(RPD, codim_index, &EGsample->endPt, *corank, OUT, RAWOUT, FAIL, T, trackCount, *smallest, *largest);

  return;
}

void regen_pos_dim_TrackLevel_trackBack(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this level                        *
\***************************************************************/
{
  int i, j, codim = RPD->codim[codim_index].codim, num_paths = RPD->codim[codim_index].num_paths;
  int num_input_vars;
  point_data_d *startPts_d = NULL;
  point_data_mp *startPts_mp = NULL;
  int retVal, indexI, indexJ, numSamples = 0, *corank = NULL;
  double trackBack_final_tol = 1e-2 * T->basicNewtonTol, *startPoint_norm_d = NULL, *smallest = NULL, *largest = NULL;
  mpf_t *startPoint_norm_mp = NULL;
  trackBack_samples_t *EGsamples = NULL;

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  // read in all of the start points
  if (T->MPType == 0 || T->MPType == 2)
  { // read in using double precision
    startPts_d = (point_data_d *)bmalloc(num_paths * sizeof(point_data_d));
    startPoint_norm_d = (double *)bmalloc(num_paths * sizeof(double));

    for (i = 0; i < num_paths; i++)
    { // setup startPts_d[i]
      init_point_data_d(&startPts_d[i], num_input_vars);
      startPts_d[i].point->size = num_input_vars;
      set_one_d(startPts_d[i].time);
      fscanf(START, "\n");
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(START, "%lf%lf", &startPts_d[i].point->coord[j].r, &startPts_d[i].point->coord[j].i);
        scanRestOfLine(START);
      }
      startPoint_norm_d[i] = infNormVec_d(startPts_d[i].point);
    }
  }
  else
  { // read in using fixed multi precision
    startPts_mp = (point_data_mp *)bmalloc(num_paths * sizeof(point_data_mp));
    startPoint_norm_mp = (mpf_t *)bmalloc(num_paths * sizeof(mpf_t));

    for (i = 0; i < num_paths; i++)
    { // setup startPts_mp[i]
      init_point_data_mp(&startPts_mp[i], num_input_vars);
      startPts_mp[i].point->size = num_input_vars;
      set_one_mp(startPts_mp[i].time);
      fscanf(START, "\n");
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(startPts_mp[i].point->coord[j].r, START, 10);
        mpf_inp_str(startPts_mp[i].point->coord[j].i, START, 10);
        scanRestOfLine(START);
      }
      mpf_init(startPoint_norm_mp[i]);
      infNormVec_mp2(startPoint_norm_mp[i], startPts_mp[i].point);
    }
  }

  // display messages
  printf("\nTracking regeneration codim %d of %d: %d path%s to track.\n", codim, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration codim %d.\n", codim);
  fprintf(OUT, "*****************************************************\n");

  // track each path for this level
  if (T->MPType == 0)
  { // loop over the paths
    for (i = 0; i < num_paths; i++)
    {     
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, startPoint_norm_d[i], NULL, startPts_d[i].point, NULL, 52);
      if (retVal)
      { // the point was found
        regen_pos_dim_TrackPath_found(indexJ, &EGsamples[indexI], corank[indexI], smallest[indexI], largest[indexI], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, &startPts_d[i], NULL, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);
        // increase size of other structures
        corank = (int *)brealloc(corank, (numSamples + 1) * sizeof(int));
        smallest = (double *)brealloc(smallest, (numSamples + 1) * sizeof(double));
        largest = (double *)brealloc(largest, (numSamples + 1) * sizeof(double));

        // track the path
        regen_pos_dim_TrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &corank[numSamples], &smallest[numSamples], &largest[numSamples], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, &startPts_d[i], NULL, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }
  }
  else if (T->MPType == 1)
  { // loop over the paths
    for (i = 0; i < num_paths; i++)
    { 
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, 0, startPoint_norm_mp[i], NULL, startPts_mp[i].point, T->Precision);
      if (retVal)
      { // the point was found
        regen_pos_dim_TrackPath_found(indexJ, &EGsamples[indexI], corank[indexI], smallest[indexI], largest[indexI], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, NULL, &startPts_mp[i], i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);
        // increase size of other structures
        corank = (int *)brealloc(corank, (numSamples + 1) * sizeof(int));
        smallest = (double *)brealloc(smallest, (numSamples + 1) * sizeof(double));
        largest = (double *)brealloc(largest, (numSamples + 1) * sizeof(double));

        // track the path
        regen_pos_dim_TrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &corank[numSamples], &smallest[numSamples], &largest[numSamples], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, NULL, &startPts_mp[i], i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }
  }
  else
  { // loop over the paths
    for (i = 0; i < num_paths; i++)
    { 
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, startPoint_norm_d[i], NULL, startPts_d[i].point, NULL, 52);
      if (retVal)
      { // the point was found
        regen_pos_dim_TrackPath_found(indexJ, &EGsamples[indexI], corank[indexI], smallest[indexI], largest[indexI], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, &startPts_d[i], NULL, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);
        // increase size of other structures
        corank = (int *)brealloc(corank, (numSamples + 1) * sizeof(int));
        smallest = (double *)brealloc(smallest, (numSamples + 1) * sizeof(double));
        largest = (double *)brealloc(largest, (numSamples + 1) * sizeof(double));

        // track the path
        regen_pos_dim_TrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &corank[numSamples], &smallest[numSamples], &largest[numSamples], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, &startPts_d[i], NULL, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }
  }

  // clear memory
  for (i = 0; i < numSamples; i++)
  {
    clear_trackBack_sample(&EGsamples[i]);
  }
  free(EGsamples);
  free(corank);
  free(smallest);
  free(largest);
  if (T->MPType == 0 || T->MPType == 2)
  { // double precision
    for (i = 0; i < num_paths; i++)
    { // clear startPts_d[i]
      clear_point_data_d(&startPts_d[i]);
    }
    free(startPts_d);
    free(startPoint_norm_d);
  }
  else
  { // fixed multi precision
    for (i = 0; i < num_paths; i++)
    { // clear startPts_mp[i]
      clear_point_data_mp(&startPts_mp[i]);
      mpf_clear(startPoint_norm_mp[i]);
    }
    free(startPts_mp);
    free(startPoint_norm_mp);
  }

  return;
}

void regen_pos_dim_TrackLevel(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this level                        *
\***************************************************************/
{
  int i, j, codim = RPD->codim[codim_index].codim, num_paths = RPD->codim[codim_index].num_paths;
  int num_input_vars;
  point_data_d startPts_d;
  point_data_mp startPts_mp;

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  init_point_data_d(&startPts_d, num_input_vars);
  init_point_data_mp(&startPts_mp, num_input_vars);
  startPts_d.point->size = startPts_mp.point->size = num_input_vars;

  // setup the current codim
  RPD->curr_codim = codim_index; 

  // display messages
  printf("\nTracking regeneration codim %d of %d: %d path%s to track.\n", codim, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration codim %d.\n", codim);
  fprintf(OUT, "*****************************************************\n");

  // track each path for this level
  for (i = 0; i < num_paths; i++)
  { 
    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // read in the start point and track the path
    if (T->MPType == 0 || T->MPType == 2)
    { // setup startPts_d
      for (j = 0; j < num_input_vars; j++)
      { // scan in the coord
        fscanf(START, "%lf%lf", &startPts_d.point->coord[j].r, &startPts_d.point->coord[j].i);
        scanRestOfLine(START);
      }
      set_one_d(startPts_d.time);

      // track the path
      regen_pos_dim_TrackPath(RPD, T, OUT, RAWOUT, MIDOUT, FAIL, &startPts_d, NULL, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
    }
    else
    { // setup startPts_mp
      for (j = 0; j < num_input_vars; j++)
      { // scan in the coord
        mpf_inp_str(startPts_mp.point->coord[j].r, START, 10);
        mpf_inp_str(startPts_mp.point->coord[j].i, START, 10);
        scanRestOfLine(START);
      }
      set_one_mp(startPts_mp.time);

      // track the path
      regen_pos_dim_TrackPath(RPD, T, OUT, RAWOUT, MIDOUT, FAIL, NULL, &startPts_mp, i, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
    }
  }

  clear_point_data_d(&startPts_d);
  clear_point_data_mp(&startPts_mp);

  return;
}

void regen_pos_dim_TrackPath(regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *)) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path                                        *
\***************************************************************/
{
  int rankType = 1, corank = 0;
  double smallest, largest;
  endgame_data_t endPt;

  init_endgame_data(&endPt, T->Precision);

  // seutp for tracking the path
  RPD->curr_codim = codim_index;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision - find corank
    endPt.last_approx_prec = 52;

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision - find corank
    endPt.last_approx_prec = T->Precision;

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, RPD, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_rank_mp(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, startPt_mp, OUT, MIDOUT, T, RPD, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP - find corank

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // check to see if it should be sharpened 
  if (endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // print the footer to OUT and print all info to RAWOUT
  printRPDFooter(RPD, codim_index, &endPt, corank, OUT, RAWOUT, FAIL, T, trackCount, smallest, largest);

  // clear memory
  clear_endgame_data(&endPt);

  return;
}

void regen_pos_dim_SortLevel(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *WITSUPER, int codim_index, int num_codim, FILE *ENDPTS, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts all paths for this level                         *
\***************************************************************/
{
  int i, j, num_input_vars, num_paths = RPD->codim[codim_index].num_paths, codim = RPD->codim[codim_index].codim;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_nonsoln = 0;
  int finite, soln, pathNum, prec, type;
  double tempD;
  
  // print header for level
  printf("\nSorting codimension %d of %d: %d path%s to sort.\n", codim, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting codimension %d.\n", codim);
  fprintf(OUT, "*****************************************************\n");

  // setup the number of variable that the points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = codim;
  else
    num_input_vars = RPD->new_variables;

  if (T->MPType == 0)
  { // sort in double precision
    endpoint_data_d endPt;
    init_endpoint_data_d(&endPt);

    // loop through the points
    for (i = 0; i < num_paths; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // pathNum & precision
      fscanf(ENDPTS, "%d\n%d\n", &pathNum, &prec);

      // coordinates
      change_size_point_d(endPt.endPt, num_input_vars);
      endPt.endPt->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(ENDPTS, "%lf%lf", &endPt.endPt->coord[j].r, &endPt.endPt->coord[j].i);
        scanRestOfLine(ENDPTS);
      }
      // time
      fscanf(ENDPTS, "%lf%lf", &endPt.finalT->r, &endPt.finalT->i);
      scanRestOfLine(ENDPTS);

      // last approx
      change_size_point_d(endPt.last_approx, num_input_vars);
      endPt.last_approx->size = num_input_vars;
      fscanf(ENDPTS, "%d\n", &j);
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(ENDPTS, "%lf%lf", &endPt.last_approx->coord[j].r, &endPt.last_approx->coord[j].i);
        scanRestOfLine(ENDPTS);
      }

      // other info
      fscanf(ENDPTS, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n", &tempD, &endPt.cond_num, &tempD, &tempD, &tempD, &tempD, &j);

      // corank, retVal, smallest non-zero & largest zero SV
      fscanf(ENDPTS, "%d\n%d\n%lf\n%lf\n", &endPt.corank, &endPt.retVal, &endPt.smallest_nonzero_SV, &endPt.largest_zero_SV);

      // determine if each path is (a) rank deficient (b) at infinity and (c) solution to the original functions
      finite = soln = 0;

      // see if this one was a success
      if (endPt.retVal == 0)
      { // determine finite & soln
        regen_pos_dim_SortEndpoint(&endPt.corank, &finite, &soln, endPt.cond_num, RPD, codim_index, pathNum, T, OUT, endPt.endPt, NULL, endPt.finalT, NULL, prec, endPt.last_approx, NULL, prec, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

        // classify the end point
        if (!finite)
        { // dehom point is infinite
          type = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (endPt.corank > 0)
          { // classify as singular
            if (soln)
              type = SOLUTION_AND_SING;
            else
              type = NONSOLUTION_AND_SING;
          }
          else
          { // classify as non-singular
            if (soln)
              type = SOLUTION_AND_NONSING;
            else 
              type = NONSOLUTION_AND_NONSING;
          }
        }
      }
      else
      { // this path was not a success - copy over error code
        type = endPt.retVal;
      }

      // check to see if this is good
      if (type == NONSOLUTION_AND_NONSING)
      { // we need to print to RAWOUT
        fprintf(RAWOUT, "%d\n%d\n", pathNum, prec);
        for (j = 0; j < num_input_vars; j++)
          fprintf(RAWOUT, "%.15e %.15e\n", endPt.endPt->coord[j].r, endPt.endPt->coord[j].i);
      }
      else if (type == SOLUTION_AND_SING || type == SOLUTION_AND_NONSING)
      { // we need to print to WITSUPER
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates first
          intrinsicToExtrinsic_d(endPt.endPt, endPt.endPt, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          intrinsicToExtrinsic_d(endPt.last_approx, endPt.last_approx, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
        }
        print_endpoint_data_d(WITSUPER, &endPt);
      }

      // add to count
      if (type == retVal_going_to_infinity || type == retVal_security_max)
        num_inf++;
      else if (type == NONSOLUTION_AND_NONSING)
        num_nonsoln++;
      else if (type == SOLUTION_AND_SING)
        num_sing++;
      else if (type == SOLUTION_AND_NONSING)
        num_nonsing++;
      else
        num_bad++;      
    }

    clear_endpoint_data_d(&endPt);
  }
  else if (T->MPType == 1)
  { // sort in fixed multi precision
    endpoint_data_mp endPt;
    init_endpoint_data_mp(&endPt);

    // loop through the points
    for (i = 0; i < num_paths; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // pathNum & precision
      fscanf(ENDPTS, "%d\n%d\n", &pathNum, &prec);

      // coordinates
      change_size_point_mp(endPt.endPt, num_input_vars);
      endPt.endPt->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(endPt.endPt->coord[j].r, ENDPTS, 10);
        mpf_inp_str(endPt.endPt->coord[j].i, ENDPTS, 10);
        scanRestOfLine(ENDPTS);
      }
      // time
      mpf_inp_str(endPt.finalT->r, ENDPTS, 10);
      mpf_inp_str(endPt.finalT->i, ENDPTS, 10);
      scanRestOfLine(ENDPTS);

      // last_approx
      change_size_point_mp(endPt.last_approx, num_input_vars);
      endPt.last_approx->size = num_input_vars;
      fscanf(ENDPTS, "%d\n", &j);
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(endPt.last_approx->coord[j].r, ENDPTS, 10);
        mpf_inp_str(endPt.last_approx->coord[j].i, ENDPTS, 10);
        scanRestOfLine(ENDPTS);
      }

      // other info
      fscanf(ENDPTS, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n", &tempD, &endPt.cond_num, &tempD, &tempD, &tempD, &tempD, &j);

      // corank, retVal, smallest non-zero & largest zero SV
      fscanf(ENDPTS, "%d\n%d\n%lf\n%lf\n", &endPt.corank, &endPt.retVal, &endPt.smallest_nonzero_SV, &endPt.largest_zero_SV);

      // determine if each path is (a) rank deficient (b) at infinity and (c) solution to the original functions
      finite = soln = 0;

      // see if this one was a success
      if (endPt.retVal == 0)
      { // determine finite & soln
        regen_pos_dim_SortEndpoint(&endPt.corank, &finite, &soln, endPt.cond_num, RPD, codim_index, pathNum, T, OUT, NULL, endPt.endPt, NULL, endPt.finalT, prec, NULL, endPt.last_approx, prec, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

        // classify the end point
        if (!finite)
        { // dehom point is infinite
          type = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (endPt.corank > 0)
          { // classify as singular
            if (soln)
              type = SOLUTION_AND_SING;
            else
              type = NONSOLUTION_AND_SING;
          }
          else
          { // classify as non-singular
            if (soln)
              type = SOLUTION_AND_NONSING;
            else
              type = NONSOLUTION_AND_NONSING;
          }
        }
      }
      else
      { // this path was not a success - copy over error code
        type = endPt.retVal;
      }

      // check to see if this is good
      if (type == NONSOLUTION_AND_NONSING)
      { // we need to print to RAWOUT
        fprintf(RAWOUT, "%d\n%d\n", pathNum, prec);
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_out_str(RAWOUT, 10, 0, endPt.endPt->coord[j].r);
          fprintf(RAWOUT, " ");
          mpf_out_str(RAWOUT, 10, 0, endPt.endPt->coord[j].i);
          fprintf(RAWOUT, "\n");
        }
      }
      else if (type == SOLUTION_AND_SING || type == SOLUTION_AND_NONSING)
      { // we need to print to WITSUPER
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates first
          intrinsicToExtrinsic_mp(endPt.endPt, endPt.endPt, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          intrinsicToExtrinsic_mp(endPt.last_approx, endPt.last_approx, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
        }
        print_endpoint_data_mp(WITSUPER, &endPt, prec);
      }

      // add to count
      if (type == retVal_going_to_infinity || type == retVal_security_max)
        num_inf++;
      else if (type == NONSOLUTION_AND_NONSING)
        num_nonsoln++;
      else if (type == SOLUTION_AND_SING)
        num_sing++;
      else if (type == SOLUTION_AND_NONSING)
        num_nonsing++;
      else
        num_bad++;
    }

    clear_endpoint_data_mp(&endPt);
  }
  else 
  { // sort using AMP
    endpoint_data_amp endPt;
    init_endpoint_data_amp(&endPt, 64, 64);
    change_size_point_d(endPt.endPt_d, num_input_vars);
    change_size_point_d(endPt.last_approx_d, num_input_vars);
    change_size_point_mp(endPt.endPt_mp, num_input_vars);
    change_size_point_mp(endPt.last_approx_mp, num_input_vars);
    endPt.endPt_d->size = endPt.last_approx_d->size = endPt.endPt_mp->size = endPt.last_approx_mp->size = num_input_vars;

    // loop through the points
    for (i = 0; i < num_paths; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // pathNum & precision
      fscanf(ENDPTS, "%d\n%d\n", &pathNum, &prec);

      // save the precision for endPt
      endPt.curr_prec = prec;

      // coordinates
      if (prec < 64)
      { // use _d
        change_size_point_d(endPt.endPt_d, num_input_vars);
        endPt.endPt_d->size = num_input_vars;
        for (j = 0; j < num_input_vars; j++)
        {
          fscanf(ENDPTS, "%lf%lf", &endPt.endPt_d->coord[j].r, &endPt.endPt_d->coord[j].i);
          scanRestOfLine(ENDPTS);
        }
        // time
        fscanf(ENDPTS, "%lf%lf", &endPt.finalT_d->r, &endPt.finalT_d->i);
      }
      else
      { // set precision & size correctly
        setprec_point_mp(endPt.endPt_mp, prec);
        change_size_point_mp(endPt.endPt_mp, num_input_vars);
        endPt.endPt_mp->size = num_input_vars;
        setprec_mp(endPt.finalT_mp, prec);

        // read in
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(endPt.endPt_mp->coord[j].r, ENDPTS, 10);
          mpf_inp_str(endPt.endPt_mp->coord[j].i, ENDPTS, 10);
          scanRestOfLine(ENDPTS);
        }
        mpf_inp_str(endPt.finalT_mp->r, ENDPTS, 10);
        mpf_inp_str(endPt.finalT_mp->i, ENDPTS, 10);
        scanRestOfLine(ENDPTS);
      }

      // last_approx
      fscanf(ENDPTS, "%d\n", &endPt.last_approx_prec);
      if (endPt.last_approx_prec < 64)
      { // use _d
        change_size_point_d(endPt.last_approx_d, num_input_vars);
        endPt.last_approx_d->size = num_input_vars;
        for (j = 0; j < num_input_vars; j++)
        {
          fscanf(ENDPTS, "%lf%lf", &endPt.last_approx_d->coord[j].r, &endPt.last_approx_d->coord[j].i);
          scanRestOfLine(ENDPTS);
        }
      }
      else
      { // set precision & size correctly
        setprec_point_mp(endPt.last_approx_mp, prec);
        change_size_point_mp(endPt.last_approx_mp, num_input_vars);
        endPt.last_approx_mp->size = num_input_vars;


        // read in
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(endPt.last_approx_mp->coord[j].r, ENDPTS, 10);
          mpf_inp_str(endPt.last_approx_mp->coord[j].i, ENDPTS, 10);
          scanRestOfLine(ENDPTS);
        }
      }

      // other info
      fscanf(ENDPTS, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n", &tempD, &endPt.cond_num, &tempD, &tempD, &tempD, &tempD, &j);

      // corank, retVal, smallest non-zero & largest zero SV
      fscanf(ENDPTS, "%d\n%d\n%lf\n%lf\n", &endPt.corank, &endPt.retVal, &endPt.smallest_nonzero_SV, &endPt.largest_zero_SV);

      // determine if each path is (a) rank deficient (b) at infinity and (c) solution to the original functions
      finite = soln = 0;

      // see if this one was a success
      if (endPt.retVal == 0)
      { // determine finite & soln
        regen_pos_dim_SortEndpoint(&endPt.corank, &finite, &soln, endPt.cond_num, RPD, codim_index, pathNum, T, OUT, endPt.endPt_d, endPt.endPt_mp, endPt.finalT_d, endPt.finalT_mp, endPt.curr_prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

        // classify the end point
        if (!finite)
        { // dehom point is infinite
          type = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (endPt.corank > 0)
          { // classify as singular
            if (soln)
              type = SOLUTION_AND_SING;
            else
              type = NONSOLUTION_AND_SING;
          }
          else
          { // classify as non-singular
            if (soln)
              type = SOLUTION_AND_NONSING;
            else
              type = NONSOLUTION_AND_NONSING;
          }
        }
      }
      else
      { // this path was not a success - copy over error code
        type = endPt.retVal;
      }

      // check to see if this is good
      if (type == NONSOLUTION_AND_NONSING)
      { // we need to print to RAWOUT
        fprintf(RAWOUT, "%d\n%d\n", pathNum, prec);
        if (prec < 64)
        { // print using _d
          for (j = 0; j < num_input_vars; j++)
            fprintf(RAWOUT, "%.15e %.15e\n", endPt.endPt_d->coord[j].r, endPt.endPt_d->coord[j].i);
        }
        else
        { // print using _mp
          for (j = 0; j < num_input_vars; j++)
          {
            mpf_out_str(RAWOUT, 10, 0, endPt.endPt_mp->coord[j].r);
            fprintf(RAWOUT, " ");
            mpf_out_str(RAWOUT, 10, 0, endPt.endPt_mp->coord[j].i);
            fprintf(RAWOUT, "\n");
          }
        }
      }
      else if (type == SOLUTION_AND_SING || type == SOLUTION_AND_NONSING)
      { // we need to print to WITSUPER
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates first
          if (endPt.curr_prec < 64)
          { // convert using _d
            intrinsicToExtrinsic_d(endPt.endPt_d, endPt.endPt_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          }
          else
          { // set the precision properly
            initMP(endPt.curr_prec);
            change_prec(RPD, endPt.curr_prec);
            intrinsicToExtrinsic_mp(endPt.endPt_mp, endPt.endPt_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          }

          if (endPt.last_approx_prec < 64)
          { // convert using _d
            intrinsicToExtrinsic_d(endPt.last_approx_d, endPt.last_approx_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          }
          else
          { // set the precision properly
            initMP(endPt.last_approx_prec);
            change_prec(RPD, endPt.last_approx_prec);
            intrinsicToExtrinsic_mp(endPt.last_approx_mp, endPt.last_approx_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          }
        }
        print_endpoint_data_amp(WITSUPER, &endPt);
      }

      // add to count
      if (type == retVal_going_to_infinity || type == retVal_security_max)
        num_inf++;
      else if (type == NONSOLUTION_AND_NONSING)
        num_nonsoln++;
      else if (type == SOLUTION_AND_SING)
        num_sing++;
      else if (type == SOLUTION_AND_NONSING)
        num_nonsing++;
      else
        num_bad++;
    }

    clear_endpoint_data_amp(&endPt);
  }

  // store the counts
  RPD->codim[codim_index].num_inf = num_inf;
  RPD->codim[codim_index].num_superset = num_sing + num_nonsing;
  RPD->codim[codim_index].num_sing = num_sing;
  RPD->codim[codim_index].num_nonsing = num_nonsing;
  RPD->codim[codim_index].num_nonsolns = num_nonsoln;
  RPD->codim[codim_index].num_bad = num_bad;

  return;
}

void regen_pos_dim_SortEndpoint(int *corank, int *finite, int *soln, double condNum, regen_pos_dim_t *RPD, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT, point_d endPt_d, point_mp endPt_mp, comp_d endTime_d, comp_mp endTime_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the point                                        *
\***************************************************************/
{ 
  int prec = MAX(64, Pt_prec);
  double norm = 0, residTol = 0;
  point_d orig_vars_d, orig_last_d, dehom_d;
  point_mp orig_vars_mp, orig_last_mp, dehom_mp;

  init_point_d(orig_vars_d, 0); init_point_d(orig_last_d, 0); init_point_d(dehom_d, 0);
  init_point_mp2(orig_vars_mp, 0, prec); init_point_mp2(orig_last_mp, 0, prec); init_point_mp2(dehom_mp, 0, prec);
  
  // initialize
  *finite = *soln = 1;

  if (Pt_prec >= 64 && T->MPType == 2)
  { // set the precision correctly
    initMP(Pt_prec);
    change_prec(RPD, Pt_prec);
  }

  // compute the original coordinates & the dehomogenized point
  if (Pt_prec < 64)
  { // compute orig_vars_d & dehom_d
    regenPDFindOrigVarsDehom_d(orig_vars_d, dehom_d, endPt_d, RPD, codim_index);
    norm = infNormVec_d(dehom_d);
  }
  else
  { // compute orig_vars_mp & dehom_mp
    regenPDFindOrigVarsDehom_mp(orig_vars_mp, dehom_mp, endPt_mp, RPD, codim_index);
    norm = infNormVec_mp(dehom_mp);
  }

  // check to see if finite
  if (T->finiteThreshold < norm)
    *finite = 0;

  // check to see if is a solution
  if (*finite)
  { // check to see if it is a solution
    if (Pt_prec < 64)
    { // setup orig_last_d
      if (last_approx_prec < 64)
      {  
        regenPDFindOrigVarsDehom_d(orig_last_d, dehom_d, last_approx_d, RPD, codim_index);
      }
      else
      { // compute _mp & then convert
        regenPDFindOrigVarsDehom_mp(orig_last_mp, dehom_mp, last_approx_mp, RPD, codim_index);
        point_mp_to_d(orig_last_d, orig_last_mp);
      }
      residTol = MAX(T->funcResTol, 1e-15);
      *soln = determineRPDSoln_d(residTol, T->ratioTol, orig_vars_d, orig_last_d, endTime_d, RPD, codim_index);
    }
    else
    { // setup orig_last_mp
      if (last_approx_prec < 64)
      { // compute _d & then convert
        regenPDFindOrigVarsDehom_d(orig_last_d, dehom_d, last_approx_d, RPD, codim_index);
        point_d_to_mp(orig_last_mp, orig_last_d);
      }
      else
      {
        regenPDFindOrigVarsDehom_mp(orig_last_mp, dehom_mp, last_approx_mp, RPD, codim_index);
      }
      residTol = pow(10, -prec_to_digits(Pt_prec) + 1);
      residTol = MAX(T->funcResTol, residTol);
      *soln = determineRPDSoln_mp(residTol, T->ratioTol, orig_vars_mp, orig_last_mp, endTime_mp, Pt_prec, RPD, codim_index);
    }
  }

  // we already know rankDef, so we can print the info
  fprintf(OUT, "Path number: %d Finite: %d Soln: %d Corank: %d CN: %e\n", pathNum, *finite, *soln, *corank, condNum);

  clear_point_d(dehom_d); clear_point_d(orig_last_d); clear_point_d(orig_vars_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last_mp); clear_point_mp(orig_vars_mp);

  return;
}

int regen_pos_dim_PrepareNextCodim(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *START, char *nextStartFile, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths actually tracked - for MIDOUT  *
* NOTES: prepares the next regen codim                          *
\***************************************************************/
{
  int i, j, retVal, num_input_vars, path_count = RPD->codim[codim_index].num_nonsolns;
  int *prec = (int *)bmalloc(path_count * sizeof(int));

  // we do not want to remove paths for having too big dehom when preparing next level
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);
  find_dehom = zero_dehom;

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  if (T->MPType == 0)
  { // prepare using double precision
    point_d *movePts = (point_d *)bmalloc(path_count * sizeof(point_d));
    
    // read in the points and convert to extrinsic coordinates
    if (RPD->codim[codim_index].useIntrinsicSlice)
    { // read in intrinsic coordinates & convert to extrinsic coordinates
      point_d tempPoint;
      init_point_d(tempPoint, num_input_vars);
      tempPoint->size = num_input_vars;

      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        fscanf(START, "%d\n%d\n", &j, &prec[i]);
        for (j = 0; j < num_input_vars; j++)
        {
          fscanf(START, "%lf%lf", &tempPoint->coord[j].r, &tempPoint->coord[j].i);
          scanRestOfLine(START);
        }

        // convert to extrinsic coordinates
        init_point_d(movePts[i], 0);
        intrinsicToExtrinsic_d(movePts[i], tempPoint, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
      }

      clear_point_d(tempPoint);
    }
    else
    { // read in extrinsic coordinates
      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        init_point_d(movePts[i], num_input_vars);
        movePts[i]->size = num_input_vars;

        fscanf(START, "%d\n%d\n", &j, &prec[i]);
        for (j = 0; j < num_input_vars; j++)
        {
          fscanf(START, "%lf%lf", &movePts[i]->coord[j].r, &movePts[i]->coord[j].i);
          scanRestOfLine(START);
        }
      }
    }

    // compute the endpoints for the next codim
    retVal = regen_pos_dim_MoveNextCodim(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, FAIL, nextStartFile, codim_index, num_codim, movePts, NULL, prec, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // clear memory
    for (i = 0; i < path_count; i++)
      clear_point_d(movePts[i]);
    free(movePts);
  }
  else if (T->MPType == 1)
  { // prepare using fixed multi precision
    point_mp *movePts = (point_mp *)bmalloc(path_count * sizeof(point_mp));

    // read in the points and convert to extrinsic coordinates
    if (RPD->codim[codim_index].useIntrinsicSlice)
    { // read in intrinsic coordinates & convert to extrinsic coordinates
      point_mp tempPoint;
      init_point_mp(tempPoint, num_input_vars);
      tempPoint->size = num_input_vars;

      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        fscanf(START, "%d\n%d\n", &j, &prec[i]);
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(tempPoint->coord[j].r, START, 10);
          mpf_inp_str(tempPoint->coord[j].i, START, 10);
          scanRestOfLine(START);
        }

        // convert to extrinsic coordinates
        init_point_mp(movePts[i], 0);
        intrinsicToExtrinsic_mp(movePts[i], tempPoint, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
      }

      clear_point_mp(tempPoint);
    }
    else
    { // read in extrinsic coordinates
      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        init_point_mp(movePts[i], num_input_vars);
        movePts[i]->size = num_input_vars;

        fscanf(START, "%d\n%d\n", &j, &prec[i]);
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(movePts[i]->coord[j].r, START, 10);
          mpf_inp_str(movePts[i]->coord[j].i, START, 10);
          scanRestOfLine(START);
        }
      }
    }

    // compute the endpoints for the next codim
    retVal = regen_pos_dim_MoveNextCodim(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, FAIL, nextStartFile, codim_index, num_codim, NULL, movePts, prec, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // clear memory
    for (i = 0; i < path_count; i++)
      clear_point_mp(movePts[i]);
    free(movePts);
  }
  else
  { // prepare using AMP
    point_d *movePts_d = (point_d *)bmalloc(path_count * sizeof(point_d));
    point_mp *movePts_mp = (point_mp *)bmalloc(path_count * sizeof(point_mp));

    // read in the points and convert to extrinsic coordinates
    if (RPD->codim[codim_index].useIntrinsicSlice)
    { // read in intrinsic coordinates & convert to extrinsic coordinates
      mat_mp B_mp;
      vec_mp p_mp;
      point_d tempPoint_d;
      point_mp tempPoint_mp;

      init_mat_mp2(B_mp, RPD->codim[codim_index].B_mp->rows, RPD->codim[codim_index].B_mp->cols, 64);
      B_mp->rows =  RPD->codim[codim_index].B_mp->rows;
      B_mp->cols =  RPD->codim[codim_index].B_mp->cols;
      change_prec_mat_mp_rat(B_mp, 96, RPD->codim[codim_index].B_rat);

      init_vec_mp2(p_mp, RPD->codim[codim_index].p_mp->size, 64);
      p_mp->size =  RPD->codim[codim_index].p_mp->size;
      change_prec_vec_mp_rat(p_mp, 96, RPD->codim[codim_index].p_rat);

      init_point_d(tempPoint_d, num_input_vars);
      init_point_mp(tempPoint_mp, num_input_vars);
      tempPoint_d->size = tempPoint_mp->size = num_input_vars;

      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        fscanf(START, "%d\n%d\n", &j, &prec[i]);

        if (prec[i] < 64)
        { // read in using double precision
          for (j = 0; j < num_input_vars; j++)
          {
            fscanf(START, "%lf%lf", &tempPoint_d->coord[j].r, &tempPoint_d->coord[j].i);
            scanRestOfLine(START);
          }
          // convert to extrinsic coordinates
          init_point_d(movePts_d[i], 0);
          intrinsicToExtrinsic_d(movePts_d[i], tempPoint_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
        }
        else
        { // read in using multi preicion
          change_prec_mat_mp_rat(B_mp, prec[i], RPD->codim[codim_index].B_rat);
          change_prec_vec_mp_rat(p_mp, prec[i], RPD->codim[codim_index].p_rat);
          setprec_point_mp(tempPoint_mp, prec[i]);
   
          for (j = 0; j < num_input_vars; j++)
          {
            mpf_inp_str(tempPoint_mp->coord[j].r, START, 10);
            mpf_inp_str(tempPoint_mp->coord[j].i, START, 10);
            scanRestOfLine(START);
          }

          // convert to extrinsic coordinates
          init_point_mp2(movePts_mp[i], 0, prec[i]);
          intrinsicToExtrinsic_mp(movePts_mp[i], tempPoint_mp, B_mp, p_mp);
        }
      }

      clear_mat_mp(B_mp);
      clear_vec_mp(p_mp);
      clear_point_d(tempPoint_d);
      clear_point_mp(tempPoint_mp);
    }
    else
    { // read in extrinsic coordinates
      for (i = 0; i < path_count; i++)
      { // read in the ith endpoint
        fscanf(START, "%d\n%d\n", &j, &prec[i]);

        if (prec[i] < 64)
        { // seutp _d
          init_point_d(movePts_d[i], num_input_vars);
          movePts_d[i]->size = num_input_vars;
          for (j = 0; j < num_input_vars; j++)
          {
            fscanf(START, "%lf%lf", &movePts_d[i]->coord[j].r, &movePts_d[i]->coord[j].i);
            scanRestOfLine(START);
          }
        }
        else
        { // setup _mp
          init_point_mp2(movePts_mp[i], num_input_vars, prec[i]);
          movePts_mp[i]->size = num_input_vars;
          for (j = 0; j < num_input_vars; j++)
          {
            mpf_inp_str(movePts_mp[i]->coord[j].r, START, 10);
            mpf_inp_str(movePts_mp[i]->coord[j].i, START, 10);
            scanRestOfLine(START);
          }
        }
      }
    }

    // compute the endpoints for the next codim
    retVal = regen_pos_dim_MoveNextCodim(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, FAIL, nextStartFile, codim_index, num_codim, movePts_d, movePts_mp, prec, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // clear memory
    for (i = 0; i < path_count; i++)
      if (prec[i] < 64)
      {
        clear_point_d(movePts_d[i]);
      }
      else
      {
        clear_point_mp(movePts_mp[i]);
      }
    free(movePts_d);
    free(movePts_mp);
  }

  // clear memory
  free(prec);

  return retVal;
}

int regen_pos_dim_MoveNextCodim(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, char *nextStartFile, int codim_index, int num_codim, point_d *movePts_d, point_mp *movePts_mp, int *Pts_prec, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths actually tracked - for MIDOUT  *
* NOTES: prepares the next regen codim                          *
\***************************************************************/
{
  int i, j, count, numNewPts, numOrigPts = RPD->codim[codim_index].num_nonsolns, next_codim_index = codim_index + 1;
  int num_bad = 0, pathCount = RPD->new_degrees[next_codim_index], *startPathNum = NULL, *rV = NULL, *new_degree = NULL, *badPaths = NULL;
  double basicTol, endgameTol, finalTol;
  FILE *NEXTSTARTPTS = fopen(nextStartFile, "w");
  mat_d B_transpose_d;
  mat_mp B_transpose_mp;
  point_data_d tempPoint_d;
  point_data_mp tempPoint_mp;

  init_mat_d(B_transpose_d, 0, 0);
  init_mat_mp(B_transpose_mp, 0, 0);
  init_point_data_d(&tempPoint_d, 0);
  init_point_data_mp(&tempPoint_mp, 0);

  // setup the number of paths that the next level 'should' have
  RPD->codim[next_codim_index].num_paths = numNewPts = numOrigPts * pathCount;

  // print the number of new start points
  fprintf(NEXTSTARTPTS, "%d\n\n", numNewPts);

  // setup a new gamma
  if (T->MPType == 0)
  { // setup _d
    get_comp_rand_d(RPD->gamma_d);
  }
  else if (T->MPType == 1)
  { // setup _mp
    get_comp_rand_mp(RPD->gamma_mp);
  }
  else
  { // setup _d, _mp, _rat
    get_comp_rand_rat(RPD->gamma_d, RPD->gamma_mp, RPD->gamma_rat, RPD->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup each of the paths to be tracked
  startPathNum = (int *)bmalloc((numNewPts - numOrigPts) * sizeof(int));
  rV = (int *)bmalloc((numNewPts - numOrigPts) * sizeof(int));
  new_degree = (int *)bmalloc((numNewPts - numOrigPts) * sizeof(int));

  // setup B_transpose, if needed
  if (RPD->codim[next_codim_index].useIntrinsicSlice)
  { // setup B_transpose
    if (T->MPType == 0 || T->MPType == 2)
    { // setup B_transpose_d
      transpose_d(B_transpose_d, RPD->codim[next_codim_index].B_d);
    }
    else
    { // setup B_transpose_mp
      transpose_mp(B_transpose_mp, RPD->codim[next_codim_index].B_mp);
    }
  }

  count = 0;
  for (i = 0; i < numOrigPts; i++)
  { // print this point to NEXTSTARTPTS
    fprintf(NEXTSTARTPTS, "%d\n", i);
    if (T->MPType == 0 || T->MPType == 2)
    { // print using double precision
      if (Pts_prec[i] < 64)
      { // see if we need to convert to intrinsic coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
          extrinsicToIntrinsic_d(movePts_d[i], movePts_d[i], B_transpose_d, RPD->codim[next_codim_index].p_d);
        }
        // print it
        for (j = 0; j < movePts_d[i]->size; j++)
        {
          fprintf(NEXTSTARTPTS, "%.15e %.15e\n", movePts_d[i]->coord[j].r, movePts_d[i]->coord[j].i);
        }
      }
      else
      { // convert to _d
        point_mp_to_d(tempPoint_d.point, movePts_mp[i]);
        // see if we need to convert to intrinsic coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
          extrinsicToIntrinsic_d(tempPoint_d.point, tempPoint_d.point, B_transpose_d, RPD->codim[next_codim_index].p_d);
        }
        // print it
        for (j = 0; j < tempPoint_d.point->size; j++)
        {
          fprintf(NEXTSTARTPTS, "%.15e %.15e\n", tempPoint_d.point->coord[j].r, tempPoint_d.point->coord[j].i);
        }
        // update movePts_mp
        point_d_to_mp(movePts_mp[i], tempPoint_d.point);
      }
    }
    else
    { // print using MP
      if (Pts_prec[i] < 64) // should NEVER happen!
      { // convert to _mp
        point_d_to_mp(tempPoint_mp.point, movePts_d[i]);
        // see if we need to convert to intrinsic coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
          extrinsicToIntrinsic_mp(tempPoint_mp.point, tempPoint_mp.point, B_transpose_mp, RPD->codim[next_codim_index].p_mp);
        }
        // print it
        for (j = 0; j < tempPoint_mp.point->size; j++)
        {
          mpf_out_str(NEXTSTARTPTS, 10, 0, tempPoint_mp.point->coord[j].r);
          fprintf(NEXTSTARTPTS, " ");
          mpf_out_str(NEXTSTARTPTS, 10, 0, tempPoint_mp.point->coord[j].i);
          fprintf(NEXTSTARTPTS, "\n");
        }
        // update movePts_d
        point_mp_to_d(movePts_d[i], tempPoint_mp.point);
      }
      else
      { // see if we need to convert to intrinsic coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
          extrinsicToIntrinsic_mp(movePts_mp[i], movePts_mp[i], B_transpose_mp, RPD->codim[next_codim_index].p_mp);
        }
        // print it
        for (j = 0; j < movePts_mp[i]->size; j++)
        {
          mpf_out_str(NEXTSTARTPTS, 10, 0, movePts_mp[i]->coord[j].r);
          fprintf(NEXTSTARTPTS, " ");
          mpf_out_str(NEXTSTARTPTS, 10, 0, movePts_mp[i]->coord[j].i);
          fprintf(NEXTSTARTPTS, "\n");
        }
      }
    }

    // setup the information for the paths that this one will generate on the next level
    for (j = 1; j < pathCount; j++)
    { // initialize values
      startPathNum[count] = i;
      rV[count] = 0;
      new_degree[count] = j;

      // increment count
      count++;
    }
  }

  // print header for level
  printf("\nPreparing regeneration codim %d of %d: %d witness point%s to move.\n", RPD->codim[next_codim_index].codim, num_codim, numNewPts - numOrigPts, (numNewPts - numOrigPts) == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Moving linears to regeneration codim %d.\n", RPD->codim[next_codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // setup the correct tracking tolerances
  basicTol = T->basicNewtonTol;
  endgameTol = T->endgameNewtonTol;
  finalTol = T->final_tolerance;
  T->basicNewtonTol = T->sliceBasicNewtonTol;
  T->endgameNewtonTol = T->sliceEndgameNewtonTol;
  T->final_tolerance = T->sliceFinalTol;

  // now that everything is setup, track the paths
  for (i = 0; i < numNewPts - numOrigPts; i++)
  {
    if (pathMod > 0 && !(i % pathMod))
      printf("Moving %d of %d\n", i, numNewPts - numOrigPts);

    // setup startPoint
    j = startPathNum[i];
    if (T->MPType == 0 || T->MPType == 2)
    { // setup tempPoint_d
      if (Pts_prec[j] < 64)
      { // use movePts_d
        point_cp_d(tempPoint_d.point, movePts_d[j]);
      }
      else
      { // use movePts_mp
        point_mp_to_d(tempPoint_d.point, movePts_mp[j]);
      }
      set_one_d(tempPoint_d.time);  
    }
    else
    { // setup tempPoint_mp
      if (Pts_prec[j] < 64)
      { // use movePts_d
        point_d_to_mp(tempPoint_mp.point, movePts_d[j]);
      }
      else
      { // use movePts_mp
        point_cp_mp(tempPoint_mp.point, movePts_mp[j]);
      }
      set_one_mp(tempPoint_mp.time);
    }

    // track the path
    rV[i] = regen_pos_dim_LinearTrackPath(numOrigPts + i, &tempPoint_d, &tempPoint_mp, new_degree[i], RPD, T, OUT, RAWOUT, MIDOUT, FAIL, NEXTSTARTPTS, codim_index, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // restore the correct tracking tolerances
  T->basicNewtonTol = basicTol;
  T->endgameNewtonTol = endgameTol;
  T->final_tolerance = finalTol;

  // check for errors
  num_bad = 0;
  for (i = 0; i < numNewPts - numOrigPts; i++)
  {
    if (rV[i] != 0 && rV[i] != retVal_reached_minTrackT)
    { // store this path to bad_paths
      printf("WARNING: When preparing path %d in codim %d, it had retVal %d.\n", i, RPD->codim[next_codim_index].codim, rV[i]);
      badPaths = (int *)brealloc(badPaths, (num_bad + 1) * sizeof(int));
      badPaths[num_bad] = startPathNum[i];
      num_bad++;
    }
  }

  // close the file
  fclose(NEXTSTARTPTS);
  // move the existing file containing the start points
  rename(nextStartFile, "regenTempFile");
  // open files
  NEXTSTARTPTS = fopen(nextStartFile, "w");
  FILE *INPUT = fopen("regenTempFile", "r");
  // setup the good ones and remove the bad ones
  if (RPD->codim[next_codim_index].useIntrinsicSlice)
  { 
    regen_pos_dim_RemoveBadPaths(NEXTSTARTPTS, INPUT, &RPD->codim[next_codim_index].num_paths, RPD->codim[next_codim_index].codim, T->MPType, numOrigPts, numNewPts - numOrigPts, startPathNum, num_bad, badPaths);
  }
  else
  {
    regen_pos_dim_RemoveBadPaths(NEXTSTARTPTS, INPUT, &RPD->codim[next_codim_index].num_paths, RPD->new_variables, T->MPType, numOrigPts, numNewPts - numOrigPts, startPathNum, num_bad, badPaths);
  }

  // print the relevant data to NEXTSTARTPTS so that it can be used to rerun this exact same problem
  printRPDRelevantData(RPD, T->MPType, next_codim_index, NEXTSTARTPTS);

  // close files
  fclose(INPUT);
  fclose(NEXTSTARTPTS);
  // remove the temporary file
  remove("regenTempFile");

  // clear memory
  free(startPathNum);
  free(rV);
  free(new_degree);
  free(badPaths);
  clear_mat_d(B_transpose_d);
  clear_mat_mp(B_transpose_mp);
  clear_point_data_d(&tempPoint_d);
  clear_point_data_mp(&tempPoint_mp);

  return (numNewPts - numOrigPts);
}

int regen_pos_dim_LinearTrackPath(int path_num, point_data_d *startPt_d, point_data_mp *startPt_mp, int new_degree, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int codim_index, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: retVal of the path                             *
* NOTES: moves the path to the next level                       *
\***************************************************************/
{
  int i, retVal = 0, next_codim_index = codim_index + 1, corank = 0; // moving full rank linears
  endgame_data_t endPt;

  init_endgame_data(&endPt, T->Precision);

  // setup for tracking the path
  RPD->curr_codim = codim_index;
  RPD->moving_degree = new_degree;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(path_num, &endPt, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }
  else if (T->MPType == 1)
  { // track the path in multi precision

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, RPD, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_mp(path_num, &endPt, startPt_mp, OUT, MIDOUT, T, RPD, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, RPD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(path_num, &endPt, startPt_d, OUT, MIDOUT, T, RPD, RPD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // print the footer
  printRPDFooter(RPD, next_codim_index, &endPt, corank, OUT, RAWOUT, FAIL, T, trackCount, 0, 0);

  // print to NEXTSTARTPTS
  fprintf(NEXTSTARTPTS, "%d\n", path_num);
  if (endPt.prec < 64)
  { 
    for (i = 0; i < endPt.PD_d.point->size; i++)
    {
      fprintf(NEXTSTARTPTS, "%.15e %.15e\n", endPt.PD_d.point->coord[i].r, endPt.PD_d.point->coord[i].i);
    }
  }
  else
  {
    for (i = 0; i < endPt.PD_mp.point->size; i++)
    {
      mpf_out_str(NEXTSTARTPTS, 10, 0, endPt.PD_mp.point->coord[i].r);
      fprintf(NEXTSTARTPTS, " ");
      mpf_out_str(NEXTSTARTPTS, 10, 0, endPt.PD_mp.point->coord[i].i);
      fprintf(NEXTSTARTPTS, "\n");
    }
  }

  // save the retVal
  retVal = endPt.retVal;

  // clear memory
  clear_endgame_data(&endPt);

  return retVal;
}

void regen_pos_dim_RemoveBadPaths(FILE *NEXTSTARTPTS, FILE *INPUT, int *num_new_paths, int num_vars, int MPType, int numOrigPts, int numMovePts, int *startPathNum, int num_bad, int *badPaths)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: removes all of the bad paths                           *
\***************************************************************/
{
  int i, j, num_curr_paths = numOrigPts + numMovePts;
  int *isGood = NULL;

  // read in the number of paths - to do error checking
  fscanf(INPUT, "%d\n\n", &i);
  if (i != num_curr_paths)
  {
    printf("ERROR: The number of paths are not equal!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup isGood
  isGood = (int *)bmalloc(num_curr_paths * sizeof(int));
  for (i = 0; i < num_curr_paths; i++)
    isGood[i] = 1;

  for (i = 0; i < numOrigPts; i++)
  { // see if this original one caused bad paths
    for (j = 0; j < num_bad && isGood[i]; j++)
      if (i == badPaths[j])
        isGood[i] = 0;
  }
  // find all the other paths that need to be removed
  for (i = 0; i < numOrigPts; i++)
    if (!isGood[i])
    { // find the ones that this one generated
      for (j = 0; j < numMovePts; j++)
        if (startPathNum[j] == i)
          isGood[j + numOrigPts] = 0;
    }

  // loop through to count the number of good paths
  *num_new_paths = 0;
  for (i = 0; i < num_curr_paths; i++)
    if (isGood[i])
      (*num_new_paths)++;

  // print the number of start points
  fprintf(NEXTSTARTPTS, "%d\n\n", *num_new_paths);

  // loop through to create the file of good start points
  if (MPType == 0 || MPType == 2)
  { // setup using double precision
    comp_d tempComp_d;

    for (i = 0; i < num_curr_paths; i++)
    { // read in the path number
      fscanf(INPUT, "%d\n", &j);
      // determine if this one is good
      if (isGood[j])
      { // print to NEXTSTARTPTS
        for (j = 0; j < num_vars; j++)
        { // read in coordinate j
          fscanf(INPUT, "%lf %lf\n", &tempComp_d->r, &tempComp_d->i);
          // print coordinate j
          fprintf(NEXTSTARTPTS, "%.15e %.15e;\n", tempComp_d->r, tempComp_d->i);
        }
        fscanf(INPUT, "\n");
        fprintf(NEXTSTARTPTS, "\n");
      }
      else
      { // move past this one
        for (j = 0; j < num_vars; j++)
        { // read in coordinate j
          fscanf(INPUT, "%lf %lf\n", &tempComp_d->r, &tempComp_d->i);
        }
        fscanf(INPUT, "\n");
      }
    }
  }
  else
  { // setup using fixed multi precision
    comp_mp tempComp_mp;
    init_mp(tempComp_mp);

    for (i = 0; i < num_curr_paths; i++)
    { // read in the path number
      fscanf(INPUT, "%d\n", &j);
      // determine if this one is good
      if (isGood[j])
      { // print to NEXTSTARTPTS
        for (j = 0; j < num_vars; j++)
        { // read in coordinate j
          mpf_inp_str(tempComp_mp->r, INPUT, 10);
          mpf_inp_str(tempComp_mp->i, INPUT, 10);
          // print coordinate j
          print_mp(NEXTSTARTPTS, 0, tempComp_mp);
          fprintf(NEXTSTARTPTS, ";\n");
        }
        fscanf(INPUT, "\n");
        fprintf(NEXTSTARTPTS, "\n");
      }
      else
      { // move past this one
        for (j = 0; j < num_vars; j++)
        { // read in coordinate j
          mpf_inp_str(tempComp_mp->r, INPUT, 10);
          mpf_inp_str(tempComp_mp->i, INPUT, 10);
        }
        fscanf(INPUT, "\n");
      }
    }

    clear_mp(tempComp_mp);
  }

  // clear
  free(isGood);

  return;
}



// EXTRA FUNCTIONS //

int determineRPDSoln_d(double tol, double ratio, point_d point, point_d last_point, comp_d time, regen_pos_dim_t *RPD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is a solution, 0 - nonsolution       *
* NOTES: determines if point satisfies the original system      *
\***************************************************************/
{
  int retVal = 1, codim = RPD->codim[codim_index].codim;

  if (codim != RPD->num_funcs)
  { // evaluate to test
    retVal = nonsolutions_check_d(RPD->num_funcs, codim, point, last_point, time, tol, ratio, RPD->Prog);
  }

  return retVal;
}

int determineRPDSoln_mp(double tol, double ratio, point_mp point, point_mp last_point, comp_mp time, int prec, regen_pos_dim_t *RPD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is a solution, 0 - nonsolution       *
* NOTES: determines if point satisfies the original system      *
\***************************************************************/
{
  int retVal = 1, codim = RPD->codim[codim_index].codim;

  if (codim != RPD->num_funcs)
  { // evaluate to test
    retVal = nonsolutions_check_mp(RPD->num_funcs, codim, point, last_point, time, tol, ratio, RPD->Prog);
  }

  return retVal;
}

void printRPDFooter(regen_pos_dim_t *RPD, int codim_index, endgame_data_t *endPt, int corank, FILE *OUT, FILE *RAWOUT, FILE *FAIL, tracker_config_t *T, trackingStats *trackCount, double smallest_nonzero_SV, double largest_zero_SV)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the footer for codim_index & path_num           *
\***************************************************************/
{
  int i, retVal, isNumber = 1, reachedMinTrackT = 0, codim = RPD->codim[codim_index].codim, prec = MAX(endPt->prec, 64);
  point_d dehom_d, orig_vars_d;
  point_mp dehom_mp, orig_vars_mp;

  init_point_d(dehom_d, 0); init_point_d(orig_vars_d, 0);
  init_point_mp2(dehom_mp, 0, prec); init_point_mp2(orig_vars_mp, 0, prec);

  // convert to original variables & dehomogenize
  if (endPt->prec < 64)
  { // compute orig_vars_d & dehom_d
    regenPDFindOrigVarsDehom_d(orig_vars_d, dehom_d, endPt->PD_d.point, RPD, codim_index);
  }
  else
  { // compute orig_vars_mp & dehom_mp
    regenPDFindOrigVarsDehom_mp(orig_vars_mp, dehom_mp, endPt->PD_mp.point, RPD, codim_index);
  }

  // print the footer

  if (endPt->retVal)
  { // we had some kind of error during tracking

    // check to see if we reached minTrackT
    if (endPt->prec < 64 && d_abs_d(endPt->PD_d.time) < T->minTrackT)
      reachedMinTrackT = 1;
    else if (endPt->prec >= 64 && d_abs_mp(endPt->PD_mp.time) < T->minTrackT)
      reachedMinTrackT = 1;

    // determine what the failure was
    if (reachedMinTrackT && endPt->retVal != retVal_Bertini_Junk)
    { // path that reached beyond minTrackT with an error code
      if (T->screenOut)
        printf("WARNING: Path %d of codim %d had retVal %d.\n", endPt->pathNum, codim, endPt->retVal);
      fprintf(OUT, "WARNING: Path %d of codim %d had retVal %d.\n", endPt->pathNum, codim, endPt->retVal);

      // consider it a success for now
      retVal = 0;
    }
    else
    { // the path had an error code

      // see if at infinity
      if (endPt->prec < 64 && infNormVec_d(dehom_d) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (endPt->prec >= 64 && infNormVec_mp(dehom_mp) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else
        retVal = endPt->retVal;

      // disply message
      if (retVal_going_to_infinity == retVal)
      {
        if (T->screenOut)
          printf("Dehom point was at infinity (inf-norm of dehom point exceeded %e).\n", T->finiteThreshold);
        fprintf(OUT, "Dehom point was at infinity (inf-norm of dehom point exceeded %e).\n", T->finiteThreshold);
      }
      else
        printResultOfPath(OUT, retVal, T);
    }
  }
  else // no error code
  {
    retVal = endPt->retVal;
    printResultOfPath(OUT, retVal, T);
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal)
  { // make sure that the output value is a number
    if (endPt->prec < 64)
    {
      for (i = 0; i < endPt->PD_d.point->size && isNumber; i++)
        if (isnan(endPt->PD_d.point->coord[i].r) || isnan(endPt->PD_d.point->coord[i].i) || isinf(endPt->PD_d.point->coord[i].r) || isinf(endPt->PD_d.point->coord[i].i))
          isNumber = 0;
    }
    else
    {
      for (i = 0; i < endPt->PD_mp.point->size && isNumber; i++)
        if (!(mpfr_number_p(endPt->PD_mp.point->coord[i].r) && mpfr_number_p(endPt->PD_mp.point->coord[i].i)))
          isNumber = 0;
    }

    if (!isNumber)
      retVal = retVal_NAN;
  }

  // print the footer to OUT
  if (endPt->prec < 64)
    printPathFooterOut_d(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, dehom_d, T, NULL, 1, 0);
  else
    printPathFooterOut_mp(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_mp, endPt->condition_number, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, endPt->first_increase, dehom_mp, T, NULL, 1, 0);

  // print the data to RAWOUT
  printRPDRawOut(RPD, RAWOUT, endPt, corank, retVal, smallest_nonzero_SV, largest_zero_SV);

  // print the data to FAIL, if needed
  if (retVal && retVal != retVal_going_to_infinity)
  { // print the path number, error message, time and point to FAIL

    // update trackCount
    trackCount->failures++;

    if (endPt->prec < 64)
      printFailureMsg_d(FAIL, &endPt->PD_d, dehom_d, endPt->pathNum, retVal, isNumber, 0, trackCount, T);
    else
      printFailureMsg_mp(FAIL, &endPt->PD_mp, dehom_mp, endPt->pathNum, retVal, isNumber, 0, trackCount, T);
  }
  else
  { // this one was had successful resolution
    trackCount->successes++;
  }

  clear_point_d(dehom_d); clear_point_d(orig_vars_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_vars_mp);

  return;
}

void printRPDRawOut(regen_pos_dim_t *RPD, FILE *RAWOUT, endgame_data_t *endPt, int corank, int retVal, double smallest_nonzero_SV, double largest_zero_SV)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the info to RAWOUT                              *
\***************************************************************/
{
  int i;

  // print path number & precision
  fprintf(RAWOUT, "%d\n%d\n", endPt->pathNum, endPt->prec);

  // print the point and time
  if (endPt->prec < 64)
  { // print using _d
    for (i = 0; i < endPt->PD_d.point->size; i++)
      fprintf(RAWOUT, "%.15e %.15e\n", endPt->PD_d.point->coord[i].r, endPt->PD_d.point->coord[i].i);
    fprintf(RAWOUT, "%.15e %.15e\n", endPt->PD_d.time->r, endPt->PD_d.time->i);
  }
  else
  { // print using _mp
    for (i = 0; i < endPt->PD_mp.point->size; i++)
    {
      mpf_out_str(RAWOUT, 10, 0, endPt->PD_mp.point->coord[i].r);
      fprintf(RAWOUT, " ");
      mpf_out_str(RAWOUT, 10, 0, endPt->PD_mp.point->coord[i].i);
      fprintf(RAWOUT, "\n");
    }
    mpf_out_str(RAWOUT, 10, 0, endPt->PD_mp.time->r);
    fprintf(RAWOUT, " ");
    mpf_out_str(RAWOUT, 10, 0, endPt->PD_mp.time->i);
    fprintf(RAWOUT, "\n");
  }

  // print last_approx
  fprintf(RAWOUT, "%d\n", endPt->last_approx_prec);
  if (endPt->last_approx_prec < 64)
  { // print using _d
    for (i = 0; i < endPt->last_approx_d->size; i++)
      fprintf(RAWOUT, "%.15e %.15e\n", endPt->last_approx_d->coord[i].r, endPt->last_approx_d->coord[i].i);
  }
  else
  { // print using _mp
    for (i = 0; i < endPt->last_approx_mp->size; i++)
    {
      mpf_out_str(RAWOUT, 10, 0, endPt->last_approx_mp->coord[i].r);
      fprintf(RAWOUT, " ");
      mpf_out_str(RAWOUT, 10, 0, endPt->last_approx_mp->coord[i].i);
      fprintf(RAWOUT, "\n");
    }
  }

  // print other data
  if (endPt->prec < 64)
  {
    fprintf(RAWOUT, "%.15e\n%.15e\n%.15e\n%.15e\n%.15e\n", endPt->function_residual_d, endPt->condition_number, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d);
    fprintf(RAWOUT, "%.15e\n%d\n%d\n%d\n", 0.0, endPt->PD_d.cycle_num, corank, retVal); // since no increase in precision, its first increase did not occur
  }
  else
  { // print using _mp
    mpf_out_str(RAWOUT, 10, 15, endPt->function_residual_mp);
    fprintf(RAWOUT, "\n%.15e\n", endPt->condition_number);
    mpf_out_str(RAWOUT, 10, 15, endPt->latest_newton_residual_mp);
    fprintf(RAWOUT, "\n");
    mpf_out_str(RAWOUT, 10, 15, endPt->t_val_at_latest_sample_point_mp);
    fprintf(RAWOUT, "\n");
    mpf_out_str(RAWOUT, 10, 15, endPt->error_at_latest_sample_point_mp);
    fprintf(RAWOUT, "\n%.15e\n%d\n%d\n%d\n", endPt->first_increase, endPt->PD_mp.cycle_num, corank, retVal);
  }
  fprintf(RAWOUT, "%.15e\n%.15e\n", smallest_nonzero_SV, largest_zero_SV);

  return;
}

void regenPDFindOrigVarsDehom_d(point_d orig_vars_d, point_d dehom_d, point_d P_d, regen_pos_dim_t *RPD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the original coordinates and dehom  *
\***************************************************************/
{
  int i;

  if (RPD->codim[codim_index].useIntrinsicSlice)
  { // first convert to extrinsic
    point_d ext_d;
    init_point_d(ext_d, 0);

    intrinsicToExtrinsic_d(ext_d, P_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);

    // change ext_d if we had an 'intrinsic' homogeneous coordinate - when the variable group was already homogenized
    if (RPD->PPD.num_hom_var_gp)
    { // find the 'intrinsic' dehom coordinate
      comp_d dehomCoord;
      set_d(dehomCoord, RPD->homVarConst_d);
      for (i = 0; i < ext_d->size; i++)
      {
        sum_mul_d(dehomCoord, &RPD->H_d->coord[i], &ext_d->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (isnan(dehomCoord->r) || isnan(dehomCoord->i) || isinf(dehomCoord->r) || isinf(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < ext_d->size; i++)
        {
          set_double_d(&ext_d->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (d_abs_d(dehomCoord) == 0)
        { // generate a random perturbation so that we can divide
          get_comp_rand_d(dehomCoord);
          mul_rdouble_d(dehomCoord, dehomCoord, 1e-16);
          recip_d(dehomCoord, dehomCoord);
        }
        else
        { // reciprocate dehomCoord
          recip_d(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < ext_d->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_d(&ext_d->coord[i], &ext_d->coord[i], dehomCoord);
        }
      }
    }

    // find orig_vars_d from ext_d
    if (RPD->new_variables != RPD->orig_variables)
    {
      mul_mat_vec_d(orig_vars_d, RPD->C_d, ext_d);
    }
    else
    {
      point_cp_d(orig_vars_d, ext_d);
    }

    clear_point_d(ext_d);
  }
  else
  { // find the actual 'new_variables'
    if (RPD->PPD.num_var_gp)
    { // the variable group was un-homogenized - find orig_vars as expected

      // convert to original coordinates
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_d(orig_vars_d, RPD->C_d, P_d);
      }
      else
      {
        point_cp_d(orig_vars_d, P_d);
      }
    }
    else
    { // the variable group was homogenized - remove the intrinsic dehom coordinate and then find orig_vars as expected
      comp_d dehomCoord;
      point_d tempPoint;

      init_point_d(tempPoint, P_d->size);
      tempPoint->size = P_d->size;
      // find the 'intrinsic' dehom coordinate
      set_d(dehomCoord, RPD->homVarConst_d);
      for (i = 0; i < P_d->size; i++)
      {
        sum_mul_d(dehomCoord, &RPD->H_d->coord[i], &P_d->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (isnan(dehomCoord->r) || isnan(dehomCoord->i) || isinf(dehomCoord->r) || isinf(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < tempPoint->size; i++)
        {
          set_double_d(&tempPoint->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (d_abs_d(dehomCoord) == 0)
        { // generate a random perturbation so that we can divide
          get_comp_rand_d(dehomCoord);
          mul_rdouble_d(dehomCoord, dehomCoord, 1e-16);
          recip_d(dehomCoord, dehomCoord);
        }
        else
        { // reciprocate dehomCoord
          recip_d(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < tempPoint->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_d(&tempPoint->coord[i], &P_d->coord[i], dehomCoord);
        }
      }

      // convert tempPoint to orig_vars_d
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_d(orig_vars_d, RPD->C_d, tempPoint);
      }
      else
      {
        point_cp_d(orig_vars_d, tempPoint);
      }

      clear_point_d(tempPoint);
    }
  }

  // find dehom_d using orig_vars_d
  getDehomPoint_d(dehom_d, orig_vars_d, orig_vars_d->size, &RPD->PPD);

  return;
}

void regenPDFindOrigVarsDehom_mp(point_mp orig_vars_mp, point_mp dehom_mp, point_mp P_mp, regen_pos_dim_t *RPD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the original coordinates and dehom  *
\***************************************************************/
{
  int i;

  if (RPD->codim[codim_index].useIntrinsicSlice)
  { // first convert to extrinsic
    point_mp ext_mp;
    init_point_mp(ext_mp, 0);

    intrinsicToExtrinsic_mp(ext_mp, P_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);

    // change ext_mp if we had an 'intrinsic' homogeneous coordinate - when the variable group was already homogenized
    if (RPD->PPD.num_hom_var_gp)
    { // find the 'intrinsic' dehom coordinate
      comp_mp dehomCoord;
      init_mp(dehomCoord);
      set_mp(dehomCoord, RPD->homVarConst_mp);
      for (i = 0; i < ext_mp->size; i++)
      {
        sum_mul_mp(dehomCoord, &RPD->H_mp->coord[i], &ext_mp->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (!mpfr_number_p(dehomCoord->r) || !mpfr_number_p(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < ext_mp->size; i++)
        {
          set_double_mp(&ext_mp->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (mpfr_zero_p(dehomCoord->r) && mpfr_zero_p(dehomCoord->i))
        { // generate a random perturbation so that we can divide
          mpf_t epsilon;
          mpf_init(epsilon);

          get_comp_rand_mp(dehomCoord);
          mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(RPD->curr_precision), __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(dehomCoord, dehomCoord, epsilon);
          recip_mp(dehomCoord, dehomCoord);

          mpf_clear(epsilon);
        }
        else
        { // reciprocate dehomCoord
          recip_mp(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < ext_mp->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_mp(&ext_mp->coord[i], &ext_mp->coord[i], dehomCoord);
        }
      }

      clear_mp(dehomCoord);
    }

    // find orig_vars_mp from ext_mp
    if (RPD->new_variables != RPD->orig_variables)
    {
      mul_mat_vec_mp(orig_vars_mp, RPD->C_mp, ext_mp);
    }
    else
    {
      point_cp_mp(orig_vars_mp, ext_mp);
    }

    clear_point_mp(ext_mp);
  }
  else
  { // find the actual 'new_variables'
    if (RPD->PPD.num_var_gp)
    { // the variable group was un-homogenized - find orig_vars as expected

      // convert to original coordinates
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_mp(orig_vars_mp, RPD->C_mp, P_mp);
      }
      else
      {
        point_cp_mp(orig_vars_mp, P_mp);
      }
    }
    else
    { // the variable group was homogenized - remove the intrinsic dehom coordinate and then find orig_vars as expected
      comp_mp dehomCoord;
      point_mp tempPoint;

      init_mp(dehomCoord);
      init_point_mp(tempPoint, P_mp->size);
      tempPoint->size = P_mp->size;
      // find the 'intrinsic' dehom coordinate
      set_mp(dehomCoord, RPD->homVarConst_mp);
      for (i = 0; i < P_mp->size; i++)
      {
        sum_mul_mp(dehomCoord, &RPD->H_mp->coord[i], &P_mp->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (!mpfr_number_p(dehomCoord->r) || !mpfr_number_p(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < tempPoint->size; i++)
        {
          set_double_mp(&tempPoint->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (mpfr_zero_p(dehomCoord->r) && mpfr_zero_p(dehomCoord->i))
        { // generate a random perturbation so that we can divide
          mpf_t epsilon;
          mpf_init(epsilon);

          get_comp_rand_mp(dehomCoord);
          mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(RPD->curr_precision), __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(dehomCoord, dehomCoord, epsilon);
          recip_mp(dehomCoord, dehomCoord);

          mpf_clear(epsilon);
        }
        else
        { // reciprocate dehomCoord
          recip_mp(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < tempPoint->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_mp(&tempPoint->coord[i], &P_mp->coord[i], dehomCoord);
        }
      }

      // convert tempPoint to orig_vars_mp
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_mp(orig_vars_mp, RPD->C_mp, tempPoint);
      }
      else
      {
        point_cp_mp(orig_vars_mp, tempPoint);
      }

      clear_mp(dehomCoord);
      clear_point_mp(tempPoint);
    }
  }

  // find dehom_mp using orig_vars_mp
  getDehomPoint_mp(dehom_mp, orig_vars_mp, orig_vars_mp->size, &RPD->PPD);

  return;
}

void regen_pos_dim_OutputChart(regen_pos_dim_t *RPD, FILE *fp, int maxCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int codim_index, total_paths = 0;

  fprintf(fp, "\n\n************ Regenerative Cascade Summary ************\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|codim|   paths   |witness superset| nonsingular | singular |nonsolutions| inf endpoints | other bad endpoints\n");
  fprintf(fp, "----------------------------------------------------------------------------------------------------------------\n");

  for (codim_index = 0; codim_index < maxCodim; codim_index++)
  {
    // add to the total paths tracked
    total_paths += RPD->codim[codim_index].num_paths;

    fprintf(fp, "| %-4d|   %-8d|   %-13d|  %-11d|  %-8d|  %-10d|   %-12d|  %d\n", RPD->codim[codim_index].codim, RPD->codim[codim_index].num_paths, RPD->codim[codim_index].num_superset, RPD->codim[codim_index].num_nonsing, RPD->codim[codim_index].num_sing, RPD->codim[codim_index].num_nonsolns, RPD->codim[codim_index].num_inf, RPD->codim[codim_index].num_bad);
  }
  fprintf(fp, "----------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|total|   %d\n\n", total_paths);

  fprintf(fp, "****************************************************\n\n");

  return;
}

void print_endpoint_data_d(FILE *fp, endpoint_data_d *endPt)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints endPt to fp                                     *
\***************************************************************/
{
  int i, prec = 52; 

  // print prec for endPt
  fprintf(fp, "%d\n", prec);
  // print endPt
  fprintf(fp, "%d\n", endPt->endPt->size);
  for (i = 0; i < endPt->endPt->size; i++)
    fprintf(fp, "%.15e %.15e\n", endPt->endPt->coord[i].r, endPt->endPt->coord[i].i);
  // print prec for last_approx
  fprintf(fp, "%d\n", prec);
  // print last_approx
  for (i = 0; i < endPt->last_approx->size; i++)
    fprintf(fp, "%.15e %.15e\n", endPt->last_approx->coord[i].r, endPt->last_approx->coord[i].i);
  // print finalT
  fprintf(fp, "%.15e %.15e\n", endPt->finalT->r, endPt->finalT->i);
  // print other data
  fprintf(fp, "%.15e\n%d\n%.15e\n%.15e\n%d\n", endPt->cond_num, endPt->corank, endPt->smallest_nonzero_SV, endPt->largest_zero_SV, endPt->retVal);

  return;
}

void setup_endpoint_data_d(endpoint_data_d *endPt, FILE *fp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup endPt from fp                                    *
\***************************************************************/
{
  int i, size;

  // read in prec & size
  fscanf(fp, "%d\n%d\n", &i, &size);
  // setup endPt & last_approx
  change_size_vec_d(endPt->endPt, size);
  change_size_vec_d(endPt->last_approx, size);
  endPt->endPt->size = endPt->last_approx->size = size;
  // read in endPt
  for (i = 0; i < size; i++)
    fscanf(fp, "%lf%lf\n", &endPt->endPt->coord[i].r, &endPt->endPt->coord[i].i);
  // read in prec
  fscanf(fp, "%d\n", &i);
  // read in last_approx
  for (i = 0; i < size; i++)
    fscanf(fp, "%lf%lf\n", &endPt->last_approx->coord[i].r, &endPt->last_approx->coord[i].i);
  // read in finalT
  fscanf(fp, "%lf%lf\n", &endPt->finalT->r, &endPt->finalT->i);
  // read in other data
  fscanf(fp, "%lf\n%d\n%lf\n%lf\n%d\n", &endPt->cond_num, &endPt->corank, &endPt->smallest_nonzero_SV, &endPt->largest_zero_SV, &endPt->retVal);

  return;
}

void print_endpoint_data_mp(FILE *fp, endpoint_data_mp *endPt, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints endPt to fp                                     *
\***************************************************************/
{
  int i;

  // print prec for endPt
  fprintf(fp, "%d\n", prec);
  // print endPt
  fprintf(fp, "%d\n", endPt->endPt->size);
  for (i = 0; i < endPt->endPt->size; i++)
  {
    mpf_out_str(fp, 10, 0, endPt->endPt->coord[i].r);
    fprintf(fp, " ");
    mpf_out_str(fp, 10, 0, endPt->endPt->coord[i].i);
    fprintf(fp, "\n");
  }
  // print prec for endPt
  fprintf(fp, "%d\n", prec);
  // print last_approx
  for (i = 0; i < endPt->last_approx->size; i++)
  {
    mpf_out_str(fp, 10, 0, endPt->last_approx->coord[i].r);
    fprintf(fp, " ");
    mpf_out_str(fp, 10, 0, endPt->last_approx->coord[i].i);
    fprintf(fp, "\n");
  }
  // print finalT
  mpf_out_str(fp, 10, 0, endPt->finalT->r);
  fprintf(fp, " ");
  mpf_out_str(fp, 10, 0, endPt->finalT->i);
  fprintf(fp, "\n");
  // print other data
  fprintf(fp, "%.15e\n%d\n%.15e\n%.15e\n%d\n", endPt->cond_num, endPt->corank, endPt->smallest_nonzero_SV, endPt->largest_zero_SV, endPt->retVal);

  return;
}

void setup_endpoint_data_mp(endpoint_data_mp *endPt, FILE *fp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup endPt from fp                                    *
\***************************************************************/
{
  int i, size;

  // read in prec & size
  fscanf(fp, "%d\n%d\n", &i, &size);
  // setup endPt & last_approx
  change_size_vec_mp(endPt->endPt, size);
  change_size_vec_mp(endPt->last_approx, size);
  endPt->endPt->size = endPt->last_approx->size = size;
  // read in endPt
  for (i = 0; i < size; i++)
  {
    mpf_inp_str(endPt->endPt->coord[i].r, fp, 10);
    mpf_inp_str(endPt->endPt->coord[i].i, fp, 10);
    scanRestOfLine(fp);
  }
  // read in prec
  fscanf(fp, "%d\n", &i);
  // read in last_approx
  for (i = 0; i < size; i++)
  {
    mpf_inp_str(endPt->last_approx->coord[i].r, fp, 10);
    mpf_inp_str(endPt->last_approx->coord[i].i, fp, 10);
    scanRestOfLine(fp);
  }
  // read in finalT
  mpf_inp_str(endPt->finalT->r, fp, 10);
  mpf_inp_str(endPt->finalT->i, fp, 10);
  scanRestOfLine(fp);
  // read in other data
  fscanf(fp, "%lf\n%d\n%lf\n%lf\n%d\n", &endPt->cond_num, &endPt->corank, &endPt->smallest_nonzero_SV, &endPt->largest_zero_SV, &endPt->retVal);

  return;
}

void print_endpoint_data_amp(FILE *fp, endpoint_data_amp *endPt)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints endPt to fp                                     *
\***************************************************************/
{
  int i;

  // print prec for endPt
  fprintf(fp, "%d\n", endPt->curr_prec);
  if (endPt->curr_prec < 64)
  { // print endPt_d
    fprintf(fp, "%d\n", endPt->endPt_d->size);
    for (i = 0; i < endPt->endPt_d->size; i++)
      fprintf(fp, "%.15e %.15e\n", endPt->endPt_d->coord[i].r, endPt->endPt_d->coord[i].i);
  }
  else
  { // print endPt_mp
    fprintf(fp, "%d\n", endPt->endPt_mp->size);
    for (i = 0; i < endPt->endPt_mp->size; i++)
    {
      mpf_out_str(fp, 10, 0, endPt->endPt_mp->coord[i].r);
      fprintf(fp, " ");
      mpf_out_str(fp, 10, 0, endPt->endPt_mp->coord[i].i);
      fprintf(fp, "\n");
    } 
  }
  // print prec for last_approx
  fprintf(fp, "%d\n", endPt->last_approx_prec);
  if (endPt->last_approx_prec < 64)
  { // print last_approx_d
    for (i = 0; i < endPt->last_approx_d->size; i++)
      fprintf(fp, "%.15e %.15e\n", endPt->last_approx_d->coord[i].r, endPt->last_approx_d->coord[i].i);
  }
  else
  { // print last_approx_mp
    for (i = 0; i < endPt->last_approx_mp->size; i++)
    {
      mpf_out_str(fp, 10, 0, endPt->last_approx_mp->coord[i].r);
      fprintf(fp, " ");
      mpf_out_str(fp, 10, 0, endPt->last_approx_mp->coord[i].i);
      fprintf(fp, "\n");
    }
  }
  
  if (endPt->curr_prec < 64)
  { // print finalT_d
    fprintf(fp, "%.15e %.15e\n", endPt->finalT_d->r, endPt->finalT_d->i);
  }
  else
  { // print finalT_mp
    mpf_out_str(fp, 10, 0, endPt->finalT_mp->r);
    fprintf(fp, " ");
    mpf_out_str(fp, 10, 0, endPt->finalT_mp->i);
    fprintf(fp, "\n");
  }
  // print other data
  fprintf(fp, "%.15e\n%d\n%.15e\n%.15e\n%d\n", endPt->cond_num, endPt->corank, endPt->smallest_nonzero_SV, endPt->largest_zero_SV, endPt->retVal);

  return;
}

void setup_endpoint_data_amp(endpoint_data_amp *endPt, FILE *fp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup endPt from fp                                    *
\***************************************************************/
{
  int i, size;

  // read in prec & size
  fscanf(fp, "%d\n%d\n", &endPt->curr_prec, &size);
  // setup endPt
  if (endPt->curr_prec < 64)
  { // setup endPt_d
    change_size_vec_d(endPt->endPt_d, size);
    endPt->endPt_d->size = size;
    // read in endPt
    for (i = 0; i < size; i++)
      fscanf(fp, "%lf%lf\n", &endPt->endPt_d->coord[i].r, &endPt->endPt_d->coord[i].i);
  }
  else
  { // setup endPt_mp
    setprec_vec_mp(endPt->endPt_mp, endPt->curr_prec);
    change_size_vec_mp(endPt->endPt_mp, size);
    endPt->endPt_mp->size = size;
    // read in endPt
    for (i = 0; i < size; i++)
    {
      mpf_inp_str(endPt->endPt_mp->coord[i].r, fp, 10);
      mpf_inp_str(endPt->endPt_mp->coord[i].i, fp, 10);
      scanRestOfLine(fp);
    }
  }
  // read in prec
  fscanf(fp, "%d\n", &endPt->last_approx_prec);
  // setup last_approx
  if (endPt->last_approx_prec < 64)
  { // setup last_approx_d
    change_size_vec_d(endPt->last_approx_d, size);
    endPt->last_approx_d->size = size;
    // read in last_approx
    for (i = 0; i < size; i++)
      fscanf(fp, "%lf%lf\n", &endPt->last_approx_d->coord[i].r, &endPt->last_approx_d->coord[i].i);
  }
  else
  { // setup last_approx_mp
    setprec_vec_mp(endPt->last_approx_mp, endPt->last_approx_prec);
    change_size_vec_mp(endPt->last_approx_mp, size);
    endPt->last_approx_mp->size = size;
    // read in last_approx
    for (i = 0; i < size; i++)
    {
      mpf_inp_str(endPt->last_approx_mp->coord[i].r, fp, 10);
      mpf_inp_str(endPt->last_approx_mp->coord[i].i, fp, 10);
      scanRestOfLine(fp);
    }
  }

  // read in finalT
  if (endPt->curr_prec < 64)
  { // setup finalT_d
    fscanf(fp, "%lf%lf\n", &endPt->finalT_d->r, &endPt->finalT_d->i);
  }
  else
  { // setup finalT_mp
    setprec_mp(endPt->finalT_mp, endPt->curr_prec);
    mpf_inp_str(endPt->finalT_mp->r, fp, 10);
    mpf_inp_str(endPt->finalT_mp->i, fp, 10);
    scanRestOfLine(fp);
  }
  // read in other data
  fscanf(fp, "%lf\n%d\n%lf\n%lf\n%d\n", &endPt->cond_num, &endPt->corank, &endPt->smallest_nonzero_SV, &endPt->largest_zero_SV, &endPt->retVal);

  return;
}

int regen_pos_dim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  regen_pos_dim_t *RPD = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    RPD = (regen_pos_dim_t *)ED_d;
    vec_d orig_d;
    init_vec_d(orig_d, 0);

    regenPDFindOrigVarsDehom_d(orig_d, out_d, in_d, RPD, RPD->curr_codim);

    clear_vec_d(orig_d);
  }
  else
  { // compute out_mp
    RPD = (regen_pos_dim_t *)ED_mp;
    vec_mp orig_mp;
    init_vec_mp2(orig_mp, 0, *out_prec);
    setprec_vec_mp(out_mp, *out_prec);

    regenPDFindOrigVarsDehom_mp(orig_mp, out_mp, in_mp, RPD, RPD->curr_codim);

    clear_vec_mp(orig_mp);
  }

  RPD = NULL;

  return 0;
}

int regen_pos_dim_setupNextStart(regen_pos_dim_t *RPD, int codim_index, tracker_config_t *T, FILE *START, FILE *NEXTSTARTPTS, FILE *PREPAREPTS, int **startPathNum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of points to move for the next codim    *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, count, num_input_vars, prec, path_count = RPD->codim[codim_index].num_nonsolns, next_codim_index = codim_index + 1;
  int new_deg = RPD->new_degrees[next_codim_index], retVal = (RPD->new_degrees[next_codim_index] - 1) * path_count;

  // allocate startPathNum
  *startPathNum = (int *)bmalloc(retVal * sizeof(int));

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  // print the number of next start points
  fprintf(NEXTSTARTPTS, "%d\n\n", new_deg * path_count);

  // read in the points and setup the files
  count = 0;
  if (T->MPType == 0)
  { // use _d
    point_d movePt_d;
    mat_d B_transpose_d;

    init_point_d(movePt_d, 0);
    init_mat_d(B_transpose_d, 0, 0);

    // setup B_tranpose_d, if needed
    if (RPD->codim[next_codim_index].useIntrinsicSlice)
    {
      transpose_d(B_transpose_d, RPD->codim[next_codim_index].B_d);
    }

    for (i = 0; i < path_count; i++)
    { // read in the point & put into extrinsic coordinates
      change_size_point_d(movePt_d, num_input_vars);
      movePt_d->size = num_input_vars;

      fscanf(START, "%d\n%d\n", &j, &prec);
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(START, "%lf%lf", &movePt_d->coord[j].r, &movePt_d->coord[j].i);
        scanRestOfLine(START);
      }
      if (RPD->codim[codim_index].useIntrinsicSlice)
      { // convert to extrinsic coordinates
        intrinsicToExtrinsic_d(movePt_d, movePt_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
      }

      // convert to the next set of coordinates
      if (RPD->codim[next_codim_index].useIntrinsicSlice)
      { // convert to intrinsic coords
        extrinsicToIntrinsic_d(movePt_d, movePt_d, B_transpose_d, RPD->codim[next_codim_index].p_d);
      }

      // print point to NEXTSTARTPTS
      fprintf(NEXTSTARTPTS, "%d\n", i);
      for (j = 0; j < movePt_d->size; j++)
        fprintf(NEXTSTARTPTS, "%.15e %.15e\n", movePt_d->coord[j].r, movePt_d->coord[j].i);

      // setup the information for the paths that this one will generate on the next level
      for (k = 1; k < new_deg; k++)
      { // initialize the value
        (*startPathNum)[count] = i;

        // print to PREPAREPTS
        for (j = 0; j < movePt_d->size; j++)
          fprintf(PREPAREPTS, "%.15e %.15e\n", movePt_d->coord[j].r, movePt_d->coord[j].i);
        fprintf(PREPAREPTS, "%d\n\n", k);

        // increment count
        count++;
      }
    }

    clear_point_d(movePt_d);
    clear_mat_d(B_transpose_d);
  }
  else if (T->MPType == 1)
  { // use _mp
    point_mp movePt_mp;
    mat_mp B_transpose_mp;

    init_point_mp(movePt_mp, 0);
    init_mat_mp(B_transpose_mp, 0, 0);

    // setup B_tranpose_d, if needed
    if (RPD->codim[next_codim_index].useIntrinsicSlice)
    {
      transpose_mp(B_transpose_mp, RPD->codim[next_codim_index].B_mp);
    }

    for (i = 0; i < path_count; i++)
    { // read in the point & put into extrinsic coordinates
      change_size_point_mp(movePt_mp, num_input_vars);
      movePt_mp->size = num_input_vars;

      fscanf(START, "%d\n%d\n", &j, &prec);
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(movePt_mp->coord[j].r, START, 10);
        mpf_inp_str(movePt_mp->coord[j].i, START, 10);
        scanRestOfLine(START);
      }
      if (RPD->codim[codim_index].useIntrinsicSlice)
      { // convert to extrinsic coordinates
        intrinsicToExtrinsic_mp(movePt_mp, movePt_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
      }

      // convert to the next set of coordinates
      if (RPD->codim[next_codim_index].useIntrinsicSlice)
      { // convert to intrinsic coords
        extrinsicToIntrinsic_mp(movePt_mp, movePt_mp, B_transpose_mp, RPD->codim[next_codim_index].p_mp);
      }

      // print point to NEXTSTARTPTS
      fprintf(NEXTSTARTPTS, "%d\n", i);
      for (j = 0; j < movePt_mp->size; j++)
      {
        mpf_out_str(NEXTSTARTPTS, 10, 0, movePt_mp->coord[j].r);
        fprintf(NEXTSTARTPTS, " ");
        mpf_out_str(NEXTSTARTPTS, 10, 0, movePt_mp->coord[j].i);
        fprintf(NEXTSTARTPTS, "\n");
      }

      // setup the information for the paths that this one will generate on the next level
      for (k = 1; k < new_deg; k++)
      { // initialize the value
        (*startPathNum)[count] = i;

        // print to PREPAREPTS
        for (j = 0; j < movePt_mp->size; j++)
        {
          mpf_out_str(PREPAREPTS, 10, 0, movePt_mp->coord[j].r);
          fprintf(PREPAREPTS, " ");
          mpf_out_str(PREPAREPTS, 10, 0, movePt_mp->coord[j].i);
          fprintf(PREPAREPTS, "\n");
        }
        fprintf(PREPAREPTS, "%d\n\n", k);

        // increment count
        count++;
      }
    }

    clear_point_mp(movePt_mp);
    clear_mat_mp(B_transpose_mp);
  }
  else 
  { // use _amp
    point_d movePt_d;
    point_mp movePt_mp;

    init_point_d(movePt_d, 0);
    init_point_mp(movePt_mp, 0);

    for (i = 0; i < path_count; i++)
    { // read in the point & put into extrinsic coordinates
      fscanf(START, "%d\n%d\n", &j, &prec);
      if (prec < 64)
      { // use _d   
        change_size_point_d(movePt_d, num_input_vars);
        movePt_d->size = num_input_vars;

        for (j = 0; j < num_input_vars; j++)
        {
          fscanf(START, "%lf%lf", &movePt_d->coord[j].r, &movePt_d->coord[j].i);
          scanRestOfLine(START);
        }
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_d(movePt_d, movePt_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
        }

        // convert to the next set of coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
	  mat_d B_transpose_d;
          init_mat_d(B_transpose_d, RPD->codim[next_codim_index].B_d->cols, RPD->codim[next_codim_index].B_d->rows);
          transpose_d(B_transpose_d, RPD->codim[next_codim_index].B_d);

          extrinsicToIntrinsic_d(movePt_d, movePt_d, B_transpose_d, RPD->codim[next_codim_index].p_d);

          clear_mat_d(B_transpose_d);
        }

        // print point to NEXTSTARTPTS
        fprintf(NEXTSTARTPTS, "%d\n", i);
        for (j = 0; j < movePt_d->size; j++)
          fprintf(NEXTSTARTPTS, "%.15e %.15e\n", movePt_d->coord[j].r, movePt_d->coord[j].i);

        // setup the information for the paths that this one will generate on the next level
        for (k = 1; k < new_deg; k++)
        { // initialize the value
          (*startPathNum)[count] = i;

          // print to PREPAREPTS
          for (j = 0; j < movePt_d->size; j++)
            fprintf(PREPAREPTS, "%.15e %.15e\n", movePt_d->coord[j].r, movePt_d->coord[j].i);
          fprintf(PREPAREPTS, "%d\n\n", k);
 
          // increment count
          count++;
        }
      }
      else
      { // use _mp
        setprec_point_mp(movePt_mp, prec);
        change_size_point_mp(movePt_mp, num_input_vars);
        movePt_mp->size = num_input_vars;

        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(movePt_mp->coord[j].r, START, 10);
          mpf_inp_str(movePt_mp->coord[j].i, START, 10);
          scanRestOfLine(START);
        }
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_mp(movePt_mp, movePt_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
        }

        // convert to the next set of coordinates
        if (RPD->codim[next_codim_index].useIntrinsicSlice)
        { // convert to intrinsic coords
          mat_mp B_transpose_mp;
          init_mat_mp2(B_transpose_mp, RPD->codim[next_codim_index].B_mp->cols, RPD->codim[next_codim_index].B_mp->rows, prec);
          transpose_mp(B_transpose_mp, RPD->codim[next_codim_index].B_mp);

          extrinsicToIntrinsic_mp(movePt_mp, movePt_mp, B_transpose_mp, RPD->codim[next_codim_index].p_mp);

          clear_mat_mp(B_transpose_mp);
        }

        // print point to NEXTSTARTPTS
        fprintf(NEXTSTARTPTS, "%d\n", i);
        for (j = 0; j < movePt_mp->size; j++)
        {
          mpf_out_str(NEXTSTARTPTS, 10, 0, movePt_mp->coord[j].r);
          fprintf(NEXTSTARTPTS, " ");
          mpf_out_str(NEXTSTARTPTS, 10, 0, movePt_mp->coord[j].i);
          fprintf(NEXTSTARTPTS, "\n");
        }

        // setup the information for the paths that this one will generate on the next level
        for (k = 1; k < new_deg; k++)
        { // initialize the value
          (*startPathNum)[count] = i;

          // print to PREPAREPTS
          for (j = 0; j < movePt_mp->size; j++)
          {
            mpf_out_str(PREPAREPTS, 10, 0, movePt_mp->coord[j].r);
            fprintf(PREPAREPTS, " ");
            mpf_out_str(PREPAREPTS, 10, 0, movePt_mp->coord[j].i);
            fprintf(PREPAREPTS, "\n");
          }
          fprintf(PREPAREPTS, "%d\n\n", k);
  
          // increment count
          count++;
        }
      }  
    }

    clear_point_d(movePt_d);
    clear_point_mp(movePt_mp);
  }

  return retVal;
}




