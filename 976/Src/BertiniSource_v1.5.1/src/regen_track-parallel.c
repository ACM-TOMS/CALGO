// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regeneration.h"
#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

void printRegenRawSort(regen_t *regen, int level_num, FILE *RAWSORT, endgame_data_t *endPt);
void printRegenRawOut(regen_t *regen, int level_num, FILE *RAWOUT, endgame_data_t *endPt, int rankDef, int retVal);

void regenTrackLevel_trackBack(int pathMod, regen_t *regen, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg);

void regen_track_seq(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, mat_d finalPatch_d, mat_mp finalPatch_mp, mpq_t ***finalPatch_rat, prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: assuming that regen is setup properly to do sequential *
* (possibly OpenMP) tracking, we want to do the regeneration    *
* process and return the non-singular, isolated endpoints in the*
* file FINALSOLN                                                *
\***************************************************************/
{
  int level, num_paths, num_crossings, num_levels = regen->num_levels;
  double midpoint_tol_prepare = T->sliceBasicNewtonTol;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);
  char *str = NULL, rawTrackFile[] = "rawout_track", rawSortFile[] = "rawout_sort", rawPrepareFile[] = "rawout_prepare";
  size_t size;
  FILE *MIDOUT = fopen(midFile, "w"), *RAWOUT = NULL, *RAWSORT = NULL, *START = NULL, *NONSOLN = NULL;

  if (regen->PPD.num_funcs > num_levels) 
  { // setup NONSOLN since system is overdetermined
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // setup the pointer to change the precision & find dehom point
  change_prec = &change_regen_prec;
  find_dehom = &regen_dehom_stack;

  for (level = startLevel; level < num_levels; level++)
  { // find the number of paths for this level
    num_paths = regen->level[level].num_paths;

    // setup RAWOUT to track this level
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, level);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, level);
    RAWOUT = fopen(str, "w");

    // setup RAWSORT to sort this level
    size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, level);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawSortFile, level);
    RAWSORT = fopen(str, "w");

    // setup START
    size = 1 + snprintf(NULL, 0, "%s_%d", startName, level);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", startName, level);
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

    // setup the evaluators for standard regeneration tracking
    ptr_to_eval_d = &regen_eval_d;
    ptr_to_eval_mp = &regen_eval_mp;

    // track each path for this level
    if (T->endgameNumber == 3)
    { // use the trackBack endgame
      regenTrackLevel_trackBack(pathMod, regen, T, START, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, level, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);
    }
    else
    { // use the standard endgame
      regenTrackLevel(pathMod, regen, T, START, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, level, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);
    }

    // close START, RAWOUT, MIDOUT & RAWSORT
    fclose(START);
    fclose(RAWOUT);
    fclose(MIDOUT);
    fclose(RAWSORT);

    // print the regeneration summary
    RAWOUT = fopen("regenSummary", "w");
    printRegenSummaryData(regen, level, RAWOUT);
    fclose(RAWOUT);

    // check for path crossings - if this is not the last level
    if (level + 1 < num_levels)
    { // make sure that no path crossing happened on this level
      num_crossings = 0;
      if (regen->level[level].useIntrinsicSlice)
        midpoint_checker(num_paths, regen->level[level].level + regen->level[level].depth, midpoint_tol, &num_crossings);
      else
        midpoint_checker(num_paths, regen->num_variables, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // setup START to be the file containing the nonsingular endpoints
    size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, level);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawSortFile, level);
    START = fopen(str, "r");
    // make sure START exits
    if (START == NULL)
    {
      printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // see if we need to prepare the next level OR create FINALSOLN
    if (level + 1 < num_levels)
    { // setup the evaluators for preparing the next level
      ptr_to_eval_d = &regen_moving_linear_eval_d;
      ptr_to_eval_mp = &regen_moving_linear_eval_mp;

      // setup the files for moving
      MIDOUT = fopen(midFile, "w"); // MIDOUT needs opened again

      size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, level + 1);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawPrepareFile, level + 1);
      RAWOUT = fopen(str, "w"); // RAWOUT needs opened again

      // setup the name of the file to contain the next set of start points
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, level + 1);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, level + 1);

      // create the next set of start points for the next level
      num_paths = regenPrepareNextLevel(pathMod, regen, T, OUT, RAWOUT, MIDOUT, FAIL, START, str, level, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

      // close RAWOUT & MIDOUT
      fclose(RAWOUT);
      fclose(MIDOUT);

      // check for path crossings
      num_crossings = 0;
      if (regen->level[level+1].useIntrinsicSlice)
        midpoint_checker(num_paths, regen->level[level+1].level + regen->level[level+1].depth, midpoint_tol_prepare, &num_crossings);
      else
        midpoint_checker(num_paths, regen->num_variables, midpoint_tol_prepare, &num_crossings);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

      // reopen MIDOUT
      MIDOUT = fopen(midFile, "w"); // MIDOUT needs opened again

      // release the memory on the previous level that is no longer needed
      clearRegenLevelStructures(1, regen, T->MPType, level);
    }
    else
    { // setup RAWOUT to be the raw track data file 
      size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, level);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawTrackFile, level);
      RAWOUT = fopen(str, "r");
      // make sure RAWOUT exits
      if (RAWOUT == NULL)
      {
        printf("ERROR: The file to contain the raw datas, '%s', does not exist!\n", str);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // setup FINALSOLN
      regenSetupFinalSoln(regen, T, START, RAWOUT, FINALSOLN, finalPatch_d, finalPatch_mp, finalPatch_rat);

      // close RAWOUT
      fclose(RAWOUT);
    }

    // close START
    fclose(START);

    // delete files that are not longer needed
    delete_omp_file(level + 1, rawTrackFile);
    delete_omp_file(level + 1, rawSortFile);
    delete_omp_file(level + 1, rawPrepareFile);
  }

  if (regen->PPD.num_funcs > num_levels) 
  { // complete NONSOLN since system was overdetermined
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // print the summary
  regenOutputChart(regen, stdout, 1);

  // set the number of paths tracked
  trackCount->numPoints = trackCount->successes + trackCount->failures;

  // release the memory on the last level
  clearRegenLevelStructures(1, regen, T->MPType, regen->num_levels - 1);

  // clear memory
  free(str);

  return;
}

void regenTrackLevel(int pathMod, regen_t *regen, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this level-use OpenMP if available*
\***************************************************************/
{
  int i, j, eachStartPt, startPointIndex, num_input_vars;
  int num_paths = regen->level[level_num].num_paths, num_levels = regen->num_levels;
  int linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  point_data_d *startPts_d = NULL;  
  point_data_mp *startPts_mp = NULL;
  int **curr_linear = NULL, **curr_linear_degree = NULL;

  // initialize the counts
  regen->level[level_num].num_inf = 0;
  regen->level[level_num].num_sing = 0;
  regen->level[level_num].num_nonsing = 0;
  regen->level[level_num].num_higher_dim = 0;
  regen->level[level_num].num_bad = 0;

  // setup the number of variables that the start points have
  if (regen->level[level_num].useIntrinsicSlice)
    num_input_vars = linearSize; // intrinsic size
  else
    num_input_vars = regen->num_variables; // extrinsic size

  // determine if we should read in each start point at the beginning
  if (num_paths < 1000) 
    eachStartPt = 0; // read in all of the start points at the beginning
  else
    eachStartPt = 1; // read in each start point when it is needed

  if (eachStartPt)
  { // only need 1 start point - it will be setup on each iteration of the loop
    curr_linear = (int **)bmalloc(1 * sizeof(int *));
    curr_linear[0] = (int *)bmalloc(linearSize * sizeof(int));
    curr_linear_degree = (int **)bmalloc(1 * sizeof(int *));
    curr_linear_degree[0] = (int *)bmalloc(linearSize * sizeof(int));

    if (T->MPType == 0 || T->MPType == 2)
    {
      startPts_d = (point_data_d *)bmalloc(1 * sizeof(point_data_d));
      init_point_data_d(&startPts_d[0], 0);
    }
    else
    {
      startPts_mp = (point_data_mp *)bmalloc(1 * sizeof(point_data_mp));
      init_point_data_mp(&startPts_mp[0], 0);
    }
  }
  else
  { // read in all of the start points
    curr_linear = (int **)bmalloc(num_paths * sizeof(int *));
    curr_linear_degree = (int **)bmalloc(num_paths * sizeof(int *));

    if (T->MPType == 0 || T->MPType == 2)
    { // read them in using double precision
      startPts_d = (point_data_d *)bmalloc(num_paths * sizeof(point_data_d));

      for (i = 0; i < num_paths; i++)
      { // setup startPts_d[i]
        init_point_data_d(&startPts_d[i], num_input_vars);
        startPts_d[i].point->size = num_input_vars;
        set_one_d(startPts_d[i].time);
        fscanf(START, "\n");
        for (j = 0; j < num_input_vars; j++)
          fscanf(START, "%lf %lf;\n", &startPts_d[i].point->coord[j].r, &startPts_d[i].point->coord[j].i);

        // setup curr_linear & curr_linear_degree
        curr_linear[i] = (int *)bmalloc(linearSize * sizeof(int));
        curr_linear_degree[i] = (int *)bmalloc(linearSize * sizeof(int));
        for (j = 0; j < linearSize; j++)
          fscanf(START, "%d %d\n", &curr_linear[i][j], &curr_linear_degree[i][j]); 
      }
    }
    else
    { // read them in using fixed multi precision
      startPts_mp = (point_data_mp *)bmalloc(num_paths * sizeof(point_data_mp));

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
          fscanf(START, ";\n");
        }

        // setup curr_linear & curr_linear_degree
        curr_linear[i] = (int *)bmalloc(linearSize * sizeof(int));
        curr_linear_degree[i] = (int *)bmalloc(linearSize * sizeof(int));
        for (j = 0; j < linearSize; j++)
          fscanf(START, "%d %d\n", &curr_linear[i][j], &curr_linear_degree[i][j]);
      }
    }
  }

  // display messages
  printf("\nTracking regeneration level %d of %d: %d path%s to track.\n", level_num, num_levels, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration level %d.\n", level_num);
  fprintf(OUT, "*****************************************************\n");

  // track each path for this level
  for (i = 0; i < num_paths; i++)
  { 
    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    if (eachStartPt)
    { // setup the next start point - THIS WILL ONLY BE DONE IN THE NON-OPENMP CASE!
      startPointIndex = 0;
      if (T->MPType == 0 || T->MPType == 2)
      { // setup startPts_d
        change_size_point_d(startPts_d[startPointIndex].point, num_input_vars);
        startPts_d[startPointIndex].point->size = num_input_vars;
        set_one_d(startPts_d[startPointIndex].time);
        fscanf(START, "\n");
        for (j = 0; j < num_input_vars; j++)
          fscanf(START, "%lf %lf;\n", &startPts_d[startPointIndex].point->coord[j].r, &startPts_d[startPointIndex].point->coord[j].i);

        // setup curr_linear & curr_linear_degree
        for (j = 0; j < linearSize; j++)
          fscanf(START, "%d %d\n", &curr_linear[startPointIndex][j], &curr_linear_degree[startPointIndex][j]);
      }
      else
      { // setup startPts_mp
        change_size_point_mp(startPts_mp[startPointIndex].point, num_input_vars);
        startPts_mp[startPointIndex].point->size = num_input_vars;
        set_one_mp(startPts_mp[startPointIndex].time);
        fscanf(START, "\n");
        for (j = 0; j < num_input_vars; j++)
        {
          mpf_inp_str(startPts_mp[startPointIndex].point->coord[j].r, START, 10);
          mpf_inp_str(startPts_mp[startPointIndex].point->coord[j].i, START, 10);
          fscanf(START, ";\n");
        }

        // setup curr_linear & curr_linear_degree
        for (j = 0; j < linearSize; j++)
          fscanf(START, "%d %d\n", &curr_linear[startPointIndex][j], &curr_linear_degree[startPointIndex][j]);
      }
    }
    else
    { // next start point is setup at index i
      startPointIndex = i;
    }

    // track the path
    if (T->MPType == 0 || T->MPType == 2)
    { // use the start point in double precision
      regenTrackPath(regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, &startPts_d[startPointIndex], NULL, i, level_num, curr_linear[startPointIndex], curr_linear_degree[startPointIndex], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);
    }
    else
    { // use the start point in fixed multi precision
      regenTrackPath(regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, NULL, &startPts_mp[startPointIndex], i, level_num, curr_linear[startPointIndex], curr_linear_degree[startPointIndex], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);
    }
  }

  // clear memory
  if (eachStartPt)
  { // only 1 start point
    free(curr_linear[0]);
    free(curr_linear_degree[0]);

    if (T->MPType == 0 || T->MPType == 2)
    {
      clear_point_data_d(&startPts_d[0]);
      free(startPts_d);
    }
    else
    {
      clear_point_data_mp(&startPts_mp[0]);
      free(startPts_mp);
    }
  }
  else
  { // all of the start points
    if (T->MPType == 0 || T->MPType == 2)
    { // double precision
      for (i = 0; i < num_paths; i++)
      { // clear startPts_d[i]
        clear_point_data_d(&startPts_d[i]);
        free(curr_linear[i]);
        free(curr_linear_degree[i]);
      }
      free(startPts_d);
    }
    else
    { // fixed multi precision
      for (i = 0; i < num_paths; i++)
      { // clear startPts_mp[i]
        clear_point_data_mp(&startPts_mp[i]);
        free(curr_linear[i]);
        free(curr_linear_degree[i]);
      }
      free(startPts_mp);
    }
  }
  free(curr_linear);
  free(curr_linear_degree);

  return;
}

void regenTrackPath(regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int level_num, int *curr_linear, int *curr_linear_degree, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the regeneration path                           *
\***************************************************************/
{
  int rankType = 0, rankDef = 0, finite = 1, higherDim = 0;
  double smallest, largest;
  endgame_data_t endPt;
  point_d dehom_d, orig_d, orig_last_d;
  point_mp dehom_mp, orig_mp, orig_last_mp;

  init_endgame_data(&endPt, T->Precision);
  init_point_d(dehom_d, 0); init_point_d(orig_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0); init_point_mp(orig_last_mp, 0);

  // seutp for tracking the path
  regen->curr_level_num = level_num;
  regen->curr_linear = curr_linear;
  regen->curr_linear_degree = curr_linear_degree;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision - find rankDef

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, &rankDef, NULL, &smallest, &largest, &endPt, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision - find rankDef

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, regen, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_rank_mp(path_num, rankType, &rankDef, NULL, &smallest, &largest, &endPt, startPt_mp, OUT, MIDOUT, T, regen, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(path_num, &endPt, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // check to see if it should be sharpened - - only endpoints at the top level need to be sharpened
  if (level_num + 1 == regen->num_levels && endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_last_d, orig_last_mp, &endPt.last_approx_prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, regen, regen);
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &endPt.prec, endPt.PD_d.point, endPt.PD_mp.point, endPt.prec, regen, regen);

  // compute the retVal
  endPt.retVal = computeRegenRetval(regen, level_num, &endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, level_num + 1 == regen->num_levels, origProg);

  if (endPt.retVal == 0)
  { // we need to check the endpoint
    regenSortEndpoint_basic(1, &rankDef, &finite, &higherDim, endPt.condition_number, regen, level_num, path_num, T, OUT, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, curr_linear, curr_linear_degree, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  } 

  // print the footer to OUT and print all info to RAWOUT & RAWSORT
  printRegenTrackFooter(endPt.retVal, regen, level_num, &endPt, rankDef, finite, higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, RAWSORT, FAIL, NONSOLN, T, trackCount, level_num + 1 == regen->num_levels);

  // clear memory
  clear_endgame_data(&endPt);
  clear_point_d(dehom_d); clear_point_d(orig_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_mp); clear_point_mp(orig_last_mp);

  return;
}

void regenTrackPath_trackBack_found(int indexJ, trackBack_samples_t *EGsample, int rankDef, int finite, int higherDim, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int level_num, int *curr_linear, int *curr_linear_degree, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the regeneration path                            *
\***************************************************************/
{
  int j, size;
  point_d dehom_d, orig_d;
  point_mp dehom_mp, orig_mp;

  init_point_d(dehom_d, 0); init_point_d(orig_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0);

  // seutp for tracking the path
  regen->curr_level_num = level_num;
  regen->curr_linear = curr_linear;
  regen->curr_linear_degree = curr_linear_degree;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // print the mid point
    fprintf(MIDOUT, "%d\n", path_num);
    size = EGsample->midPt_d[indexJ]->size;
    for (j = 0; j < size; j++)
      fprintf(MIDOUT, "%.15e %.15e\n", EGsample->midPt_d[indexJ]->coord[j].r, EGsample->midPt_d[indexJ]->coord[j].i);
    fprintf(MIDOUT, "\n");

  }
  else if (T->MPType == 1)
  { // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, regen, ptr_to_eval_mp);

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
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

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

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &EGsample->endPt.prec, EGsample->endPt.PD_d.point, EGsample->endPt.PD_mp.point, EGsample->endPt.prec, regen, regen);

  // print the footer to OUT and print all info to RAWOUT & RAWSORT
  printRegenTrackFooter(EGsample->endPt.retVal, regen, level_num, &EGsample->endPt, rankDef, finite, higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, RAWSORT, FAIL, NONSOLN, T, trackCount, level_num + 1 == regen->num_levels);

  // clear memory
  clear_point_d(dehom_d); clear_point_d(orig_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_mp);

  return;
}

void regenTrackPath_trackBack(double trackBack_final_tol, trackBack_samples_t *EGsample, int *rankDef, int *finite, int *higherDim, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int level_num, int *curr_linear, int *curr_linear_degree, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the regeneration path                           *
\***************************************************************/
{
  int rankType = 0;
  double smallest, largest;
  point_d dehom_d, orig_d, orig_last_d;
  point_mp dehom_mp, orig_mp, orig_last_mp;

  init_point_d(dehom_d, 0); init_point_d(orig_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0); init_point_mp(orig_last_mp, 0);

  // seutp for tracking the path
  regen->curr_level_num = level_num;
  regen->curr_linear = curr_linear;
  regen->curr_linear_degree = curr_linear_degree;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision - find rankDef

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_rank_d(path_num, rankType, rankDef, NULL, &smallest, &largest, trackBack_final_tol, EGsample, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision - find rankDef

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, regen, ptr_to_eval_mp);

    // track the path
    zero_dim_trackBack_path_rank_mp(path_num, rankType, rankDef, NULL, &smallest, &largest, trackBack_final_tol, EGsample, startPt_mp, OUT, MIDOUT, T, regen, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_d(path_num, trackBack_final_tol, EGsample, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // check to see if it should be sharpened - - only endpoints at the top level need to be sharpened
  if (level_num + 1 == regen->num_levels && EGsample->endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&EGsample->endPt, T, OUT, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_last_d, orig_last_mp, &EGsample->endPt.last_approx_prec, EGsample->endPt.last_approx_d, EGsample->endPt.last_approx_mp, EGsample->endPt.last_approx_prec, regen, regen);
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &EGsample->endPt.prec, EGsample->endPt.PD_d.point, EGsample->endPt.PD_mp.point, EGsample->endPt.prec, regen, regen);

  // compute the retVal
  EGsample->endPt.retVal = computeRegenRetval(regen, level_num, &EGsample->endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, level_num + 1 == regen->num_levels, origProg);

  if (EGsample->endPt.retVal == 0)
  { // we need to check the endpoint
    regenSortEndpoint_basic(1, rankDef, finite, higherDim, EGsample->endPt.condition_number, regen, level_num, path_num, T, OUT, &EGsample->endPt.PD_d, &EGsample->endPt.PD_mp, EGsample->endPt.prec, EGsample->endPt.last_approx_d, EGsample->endPt.last_approx_mp, EGsample->endPt.last_approx_prec, curr_linear, curr_linear_degree, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // print the footer to OUT and print all info to RAWOUT & RAWSORT
  printRegenTrackFooter(EGsample->endPt.retVal, regen, level_num, &EGsample->endPt, *rankDef, *finite, *higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, RAWSORT, FAIL, NONSOLN, T, trackCount, level_num + 1 == regen->num_levels);

  // clear memory
  clear_point_d(dehom_d); clear_point_d(orig_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_mp); clear_point_mp(orig_last_mp);

  return;
}

void regenTrackLevel_trackBack(int pathMod, regen_t *regen, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this level                        *
\***************************************************************/
{
  int i, j, num_input_vars;
  int num_paths = regen->level[level_num].num_paths, num_levels = regen->num_levels;
  int linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  point_data_d *startPts_d = NULL;
  point_data_mp *startPts_mp = NULL;
  int retVal, indexI, indexJ, numSamples = 0, **curr_linear = NULL, **curr_linear_degree = NULL;
  double trackBack_final_tol = 1e-2 * T->basicNewtonTol, *startPoint_norm_d = NULL;
  mpf_t *startPoint_norm_mp = NULL;
  trackBack_samples_t *EGsamples = NULL;
  int *rankDef = NULL, *finite = NULL, *higherDim = NULL;

  // initialize the counts
  regen->level[level_num].num_inf = 0;
  regen->level[level_num].num_sing = 0;
  regen->level[level_num].num_nonsing = 0;
  regen->level[level_num].num_higher_dim = 0;
  regen->level[level_num].num_bad = 0;

  // setup the number of variables that the start points have
  if (regen[0].level[level_num].useIntrinsicSlice)
    num_input_vars = linearSize; // intrinsic size
  else
    num_input_vars = regen->num_variables; // extrinsic size

  // read in all of the start points
  curr_linear = (int **)bmalloc(num_paths * sizeof(int *));
  curr_linear_degree = (int **)bmalloc(num_paths * sizeof(int *));
  rankDef = (int *)bmalloc(num_paths * sizeof(int));
  finite = (int *)bmalloc(num_paths * sizeof(int));
  higherDim = (int *)bmalloc(num_paths * sizeof(int));
  for (i = 0; i < num_paths; i++)
  {
    rankDef[i] = higherDim[i] = 0;
    finite[i] = 1;
  }

  if (T->MPType == 0 || T->MPType == 2)
  { // read them in using double precision
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
        fscanf(START, "%lf %lf;\n", &startPts_d[i].point->coord[j].r, &startPts_d[i].point->coord[j].i);
      }
      startPoint_norm_d[i] = infNormVec_d(startPts_d[i].point);

      // setup curr_linear & curr_linear_degree
      curr_linear[i] = (int *)bmalloc(linearSize * sizeof(int));
      curr_linear_degree[i] = (int *)bmalloc(linearSize * sizeof(int));
      for (j = 0; j < linearSize; j++)
        fscanf(START, "%d %d\n", &curr_linear[i][j], &curr_linear_degree[i][j]);
    }
  }
  else
  { // read them in using fixed multi precision
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
        fscanf(START, ";\n");
      }

      mpf_init(startPoint_norm_mp[i]);
      infNormVec_mp2(startPoint_norm_mp[i], startPts_mp[i].point);

      // setup curr_linear & curr_linear_degree
      curr_linear[i] = (int *)bmalloc(linearSize * sizeof(int));
      curr_linear_degree[i] = (int *)bmalloc(linearSize * sizeof(int));
      for (j = 0; j < linearSize; j++)
        fscanf(START, "%d %d\n", &curr_linear[i][j], &curr_linear_degree[i][j]);
    }
  }

  // display messages
  printf("\nTracking regeneration level %d of %d: %d path%s to track.\n", level_num, num_levels, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration level %d.\n", level_num);
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
        regenTrackPath_trackBack_found(indexJ, &EGsamples[indexI], rankDef[indexI], finite[indexI], higherDim[indexI], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, &startPts_d[i], NULL, i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, origProg);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);

        // track the path
        regenTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &rankDef[numSamples], &finite[numSamples], &higherDim[numSamples], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, &startPts_d[i], NULL, i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);

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
        regenTrackPath_trackBack_found(indexJ, &EGsamples[indexI], rankDef[indexI], finite[indexI], higherDim[indexI], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, NULL, &startPts_mp[i], i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, origProg);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);

        // track the path
        regenTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &rankDef[numSamples], &finite[numSamples], &higherDim[numSamples], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, NULL, &startPts_mp[i], i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);

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
        regenTrackPath_trackBack_found(indexJ, &EGsamples[indexI], rankDef[indexI], finite[indexI], higherDim[indexI], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, &startPts_d[i], NULL, i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, origProg);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);

        // track the path
        regenTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &rankDef[numSamples], &finite[numSamples], &higherDim[numSamples], regen, T, OUT, RAWOUT, RAWSORT, MIDOUT, FAIL, NONSOLN, &startPts_d[i], NULL, i, level_num, curr_linear[i], curr_linear_degree[i], trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom, origProg);

        // increment the size
        numSamples++;
      }
    }
  }

  // clear memory
  for (i = 0; i < numSamples; i++)
    clear_trackBack_sample(&EGsamples[i]);
  free(EGsamples);
  if (T->MPType == 0 || T->MPType == 2)
  { // double precision
    for (i = 0; i < num_paths; i++)
    { // clear startPts_d[i]
      clear_point_data_d(&startPts_d[i]);
      free(curr_linear[i]);
      free(curr_linear_degree[i]);
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
      free(curr_linear[i]);
      free(curr_linear_degree[i]);
    }
    free(startPts_mp);
    free(startPoint_norm_mp);
  }
  free(curr_linear);
  free(curr_linear_degree);
  free(rankDef);
  free(finite);
  free(higherDim);

  return;
}

int computeRegenRetval(regen_t *regen, int level_num, endgame_data_t *endPt, point_d orig_d, point_mp orig_mp, point_d orig_last_d, point_mp orig_last_mp, point_d dehom_d, point_mp dehom_mp, tracker_config_t *T, int junkCheck, prog_t *origProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computes the actual retVal                             *
\***************************************************************/
{
  int i, retVal, isNumber = 1, reachedMinTrackT, isSoln = 1;

  if (endPt->retVal)
  { // we had some kind of error during tracking

    // check to see if we reached minTrackT
    if (endPt->prec < 64 && d_abs_d(endPt->PD_d.time) < T->minTrackT)
      reachedMinTrackT = 1;
    else if (endPt->prec >= 64 && d_abs_mp(endPt->PD_mp.time) < T->minTrackT)
      reachedMinTrackT = 1;
    else
      reachedMinTrackT = 0;

    // determine what the failure was
    if (reachedMinTrackT && endPt->retVal != retVal_Bertini_Junk)
    { // consider it a success for now
      retVal = 0;
    }
    else
    { // the path had an error code

      // see if went to infinity
      if (endPt->prec < 64 && infNormVec_d(orig_d) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (endPt->prec >= 64 && infNormVec_mp(orig_mp) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (T->regen_remove_inf && endPt->prec < 64 && infNormVec_d(dehom_d) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (T->regen_remove_inf && endPt->prec >= 64 && infNormVec_mp(dehom_mp) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else
        retVal = endPt->retVal;
    }
  }
  else 
  { // no error code
    retVal = endPt->retVal;
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

  // run bertini junk checker if it is a number that looks like a successful path on the top level
  if (junkCheck && !retVal && endPt->prec < 64)
  { // check for junk
    if (endPt->last_approx_prec > 52)
    { // convert to _d
      point_mp_to_d(orig_last_d, orig_last_mp);
    }

    isSoln = nonsolutions_check_d(origProg->numFuncs, regen->num_funcs, orig_d, orig_last_d, endPt->PD_d.time, T->funcResTol, T->ratioTol, origProg);

    if (!isSoln)
    { // classify as junk, unless it is at infinity and then we can classify as infinite
      if (infNormVec_d(orig_d) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (T->regen_remove_inf && infNormVec_d(dehom_d) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else
        retVal = retVal_Bertini_Junk;
    }
  }
  else if (junkCheck && !retVal)
  { // check for junk
    if (endPt->last_approx_prec < 64)
    { // convert to _mp
      point_d_to_mp(orig_last_mp, orig_last_d);
    }

    isSoln = nonsolutions_check_mp(origProg->numFuncs, regen->num_funcs, orig_mp, orig_last_mp, endPt->PD_mp.time, T->funcResTol, T->ratioTol, origProg);

    if (!isSoln)
    { // classify as junk, unless it is at infinity and then we can classify as infinite
      if (infNormVec_mp(orig_mp) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else if (T->regen_remove_inf && infNormVec_mp(dehom_mp) > T->finiteThreshold)
        retVal = retVal_going_to_infinity;
      else
        retVal = retVal_Bertini_Junk;
    }
  }

  return retVal;
}

void printRegenTrackFooter(int retVal, regen_t *regen, int level_num, endgame_data_t *endPt, int rankDef, int finite, int higherDim, point_d orig_d, point_mp orig_mp, point_d dehom_d, point_mp dehom_mp, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *FAIL, FILE *NONSOLN, tracker_config_t *T, trackingStats *trackCount, int updateTrackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the footer                                      *
\***************************************************************/
{
  int i, moveToNext = 0, isNumber = 1, isJunk = 0;
  
  printResultOfPath(OUT, retVal, T);

  isNumber = (retVal != retVal_NAN);
  isJunk = (retVal == retVal_Bertini_Junk);

  // increment the counts
  if (retVal == retVal_going_to_infinity || retVal == retVal_security_max || !finite)
    regen->level[level_num].num_inf++;
  else if (retVal == retVal_higher_dim || higherDim)
    regen->level[level_num].num_higher_dim++;
  else if (rankDef == 1)
    regen->level[level_num].num_sing++;
  else if (retVal)
    regen->level[level_num].num_bad++;
  else
  { // this one gets moved to next level
    moveToNext = 1;
    regen->level[level_num].num_nonsing++;
  }

  // print the footer to OUT
  if (endPt->prec < 64)
    printPathFooterOut_d(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, dehom_d, T, NULL, 1, 0);
  else
    printPathFooterOut_mp(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_mp, endPt->condition_number, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, endPt->first_increase, dehom_mp, T, NULL, 1, 0);

  // print the data to RAWOUT
  printRegenRawOut(regen, level_num, RAWOUT, endPt, rankDef, retVal);

  // print the data to RAWSORT
  if (moveToNext)
  { 
    printRegenRawSort(regen, level_num, RAWSORT, endPt);
  }

  // print the data to FAIL, if needed
  if (retVal && retVal != retVal_going_to_infinity)
  { // print the path number, error message, time and point to FAIL
    if (updateTrackCount)
    { // update trackCount
      trackCount->failures++;

      if (endPt->prec < 64)
      {
        printFailureMsg_d(FAIL, &endPt->PD_d, dehom_d, endPt->pathNum, retVal, isNumber, isJunk, trackCount, T);

        if (isJunk)
        { // print to NONSOLN
          for (i = 0; i < dehom_d->size; i++)
          {
            print_d(NONSOLN, 0, &dehom_d->coord[i]);
            fprintf(NONSOLN, "\n");
          } 
          fprintf(NONSOLN, "\n");
        }
      }
      else
      {
        printFailureMsg_mp(FAIL, &endPt->PD_mp, dehom_mp, endPt->pathNum, retVal, isNumber, isJunk, trackCount, T);

        if (isJunk)
        { // print to NONSOLN
          for (i = 0; i < dehom_mp->size; i++)
          {
            print_mp(NONSOLN, 0, &dehom_mp->coord[i]);
            fprintf(NONSOLN, "\n");
          }
          fprintf(NONSOLN, "\n");
        }
      }
    }
    else
    { // update a fake trackCount
      trackingStats tC;

      if (endPt->prec < 64)
        printFailureMsg_d(FAIL, &endPt->PD_d, dehom_d, endPt->pathNum, retVal, isNumber, 0, &tC, T);
      else
        printFailureMsg_mp(FAIL, &endPt->PD_mp, dehom_mp, endPt->pathNum, retVal, isNumber, 0, &tC, T);
    }
  }
  else if (updateTrackCount)
  { // this one was had successful resolution
    trackCount->successes++;
  }

  return;
}

int printRegenLinearTrackFooter(int retVal, regen_t *regen, int level_num, int path_num, endgame_data_t *endPt, int rankDef, int finite, int higherDim, point_d orig_d, point_mp orig_mp, point_d dehom_d, point_mp dehom_mp, FILE *OUT, FILE *RAWOUT, FILE *NEXTSTARTPTS, FILE *FAIL, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - good point, 1 - not good point             *
* NOTES: prints the footer                                      *
\***************************************************************/
{
  int notGoodPoint = 0, isNumber = 1, isJunk = 0;
  int linearSize = regen->level[level_num].level + regen->level[level_num].depth;

  printResultOfPath(OUT, retVal, T);

  isNumber = (retVal != retVal_NAN);
  isJunk = (retVal == retVal_Bertini_Junk);

  // setup retVal
  if (retVal == 0 && rankDef == 0 && higherDim == 0 && finite)
    notGoodPoint = 0; // good point
  else
    notGoodPoint = 1; // not good point

  // print the footer to OUT
  if (endPt->prec < 64)
    printPathFooterOut_d(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, dehom_d, T, NULL, 1, 0);
  else
    printPathFooterOut_mp(OUT, RAWOUT, 0, endPt->pathNum, &endPt->PD_mp, endPt->condition_number, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, endPt->first_increase, dehom_mp, T, NULL, 1, 0);

  // print the data to RAWOUT
  printRegenRawOut(regen, level_num, RAWOUT, endPt, rankDef, retVal);

  // print to NEXTSTARTPTS
  fprintf(NEXTSTARTPTS, "%d\n", path_num);
  printPointLinearDegree(NEXTSTARTPTS, endPt->PD_d.point, endPt->PD_mp.point, endPt->prec, regen->curr_linear, regen->curr_linear_degree, linearSize);

  // print the data to FAIL, if needed
  if (retVal && retVal != retVal_going_to_infinity)
  { // print the path number, error message, time and point to FAIL
    // update a fake trackCount
    trackingStats tC;

    if (endPt->prec < 64)
      printFailureMsg_d(FAIL, &endPt->PD_d, dehom_d, endPt->pathNum, retVal, isNumber, 0, &tC, T);
    else
      printFailureMsg_mp(FAIL, &endPt->PD_mp, dehom_mp, endPt->pathNum, retVal, isNumber, 0, &tC, T);
  }

  return notGoodPoint;
}


void printRegenRawSort(regen_t *regen, int level_num, FILE *RAWSORT, endgame_data_t *endPt)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the info to RAWSORT                             *
\***************************************************************/
{
  int j, linearSize = regen->level[level_num].level + regen->level[level_num].depth;

  fprintf(RAWSORT, "%d\n%d\n", endPt->pathNum, endPt->prec);
  if (endPt->prec < 64)
  { // print using _d
    for (j = 0; j < endPt->PD_d.point->size; j++)
      fprintf(RAWSORT, "%.15e %.15e\n", endPt->PD_d.point->coord[j].r, endPt->PD_d.point->coord[j].i);
  }
  else
  { // print using _mp
    for (j = 0; j < endPt->PD_mp.point->size; j++)
    {
      mpf_out_str(RAWSORT, 10, 0, endPt->PD_mp.point->coord[j].r);
      fprintf(RAWSORT, " ");
      mpf_out_str(RAWSORT, 10, 0, endPt->PD_mp.point->coord[j].i);
      fprintf(RAWSORT, "\n");
    }
  }
  
  for (j = 0; j < linearSize; j++)
    fprintf(RAWSORT, "%d %d\n", regen->curr_linear[j], regen->curr_linear_degree[j]);

  return;
}

void printRegenRawOut(regen_t *regen, int level_num, FILE *RAWOUT, endgame_data_t *endPt, int rankDef, int retVal)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the info to RAWOUT                              *
\***************************************************************/
{
  int i, linearSize = regen->level[level_num].level + regen->level[level_num].depth;

  // print path number & precision
  fprintf(RAWOUT, "%d\n%d\n", endPt->pathNum, endPt->prec);

  // print the point and other data
  if (endPt->prec < 64)
  { // print using _d
    for (i = 0; i < endPt->PD_d.point->size; i++)
      fprintf(RAWOUT, "%.15e %.15e\n", endPt->PD_d.point->coord[i].r, endPt->PD_d.point->coord[i].i);

    fprintf(RAWOUT, "%.15e\n%.15e\n%.15e\n%.15e\n%.15e\n", endPt->function_residual_d, endPt->condition_number, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d);
    fprintf(RAWOUT, "%.15e\n%d\n%d\n%d\n", 0.0, endPt->PD_d.cycle_num, rankDef, retVal); // since no increase in precision, its first increase did not occur
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

    mpf_out_str(RAWOUT, 10, 15, endPt->function_residual_mp);
    fprintf(RAWOUT, "\n%.15e\n", endPt->condition_number);
    mpf_out_str(RAWOUT, 10, 15, endPt->latest_newton_residual_mp);
    fprintf(RAWOUT, "\n");
    mpf_out_str(RAWOUT, 10, 15, endPt->t_val_at_latest_sample_point_mp);
    fprintf(RAWOUT, "\n");
    mpf_out_str(RAWOUT, 10, 15, endPt->error_at_latest_sample_point_mp);
    fprintf(RAWOUT, "\n%.15e\n%d\n%d\n%d\n", endPt->first_increase, endPt->PD_mp.cycle_num, rankDef, retVal);
  }

  // print curr_linear & curr_linear_degree
  for (i = 0; i < linearSize; i++)
    fprintf(RAWOUT, "%d %d\n", regen->curr_linear[i], regen->curr_linear_degree[i]);

  return;
}

void regenSortEndpoint_basic(int checkRankDef, int *rankDef, int *finite, int *higherDim, double condNum, regen_t *regen, int level_num, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int *curr_linear, int *curr_linear_degree, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines rankDef, finite & higherDim                 *
*   1) finite 2) higherDim 3) rankDef - if any fail, we quit!   *
\***************************************************************/
{ // initialize higherDim & finite (rankDef may be known!)
  *higherDim = 0;
  *finite = 1;

  // setup for sorting the endpoint
  regen->curr_level_num = level_num;
  regen->curr_linear = curr_linear;
  regen->curr_linear_degree = curr_linear_degree;

  // check to see if finite
  *finite = determineRegenFinite(T->regen_remove_inf, T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, regen, level_num);

  if (T->MPType == 2)
  { // determine if we need to do the rank deficient test (which also refines the endpoint)
    if (checkRankDef && *finite)
    { // determine if it is rank deficient
      *rankDef = determineRankDef(&condNum, T->final_tolerance, NULL, 52, endPt_d, endPt_mp, endPt_prec, T, OUT, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

      // determine if it is finite or infinite
      *finite = determineRegenFinite(T->regen_remove_inf, T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, regen, level_num);
    }
     
    // if we are finite and smooth, check to see if higher-dimensional
    if (*rankDef == 0 && *finite && T->regen_higher_dim_check)
    { // determine if it is on higher dimensional component
      *higherDim = determineRegenHigherDim(T->funcResTol, T->ratioTol, endPt_d, endPt_mp, endPt_prec, last_approx_d, last_approx_mp, last_approx_prec, regen, level_num);
    }
  }
  else if (*finite && T->regen_higher_dim_check)
  { // check to see if on a higher-dimensional set
    *higherDim = determineRegenHigherDim(T->funcResTol, T->ratioTol, endPt_d, endPt_mp, endPt_prec, last_approx_d, last_approx_mp, last_approx_prec, regen, level_num);
  }

  // print the info
  if (T->regen_higher_dim_check)
  {
    fprintf(OUT, "Higher Dim'l: %d", *higherDim);
    if (checkRankDef)
      fprintf(OUT, " RankDef: %d", *rankDef);
    fprintf(OUT, " Finite: %d", *finite);
  }
  else 
  {
    if (checkRankDef)
    {
      fprintf(OUT, "RankDef: %d", *rankDef);
      fprintf(OUT, " Finite: %d", *finite);
    }
    else
      fprintf(OUT, "Finite: %d", *finite);
  } 
  fprintf(OUT, " CN: %e\n", condNum);

  return;
}

int determineRegenFinite(int remove_inf, double maxNorm, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, regen_t *regen, int level_num)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether it is finite or not                    *
* NOTES:                                                        *
\***************************************************************/
{
  int finite = 0, out_prec = Pt_prec;

  if (remove_inf)
  { // do the normal test
    if (Pt_prec < 64)
    { // check using _d
      vec_d orig_d, dehom_d;
      init_vec_d(orig_d, 0); init_vec_d(dehom_d, 0);

      regen_dehom2(dehom_d, NULL, orig_d, NULL, &out_prec, Pt_d->point, NULL, Pt_prec, regen, regen);

      if (infNormVec_d(dehom_d) < maxNorm && infNormVec_d(orig_d) < maxNorm)
        finite = 1;
      else
        finite = 0;

      clear_vec_d(orig_d); clear_vec_d(dehom_d);
    }
    else
    { // check using _mp
      vec_mp orig_mp, dehom_mp;
      init_point_mp2(orig_mp, 0, Pt_prec); init_point_mp2(dehom_mp, 0, Pt_prec);

      regen_dehom2(NULL, dehom_mp, NULL, orig_mp, &out_prec, NULL, Pt_mp->point, Pt_prec, regen, regen);

      if (infNormVec_mp(dehom_mp) < maxNorm && infNormVec_mp(orig_mp) < maxNorm)
        finite = 1;
      else
        finite = 0;

      clear_point_mp(orig_mp); clear_point_mp(dehom_mp);
    }
  }
  else
  { // only check that orig coordinates are finite
    if (Pt_prec < 64)
    { // check using _d
      vec_d orig_d;
      init_vec_d(orig_d, 0);

      regen_dehom_hom(orig_d, NULL, &out_prec, Pt_d->point, NULL, Pt_prec, regen, regen);

      if (infNormVec_d(orig_d) < maxNorm)
        finite = 1;
      else
        finite = 0;

      clear_vec_d(orig_d);
    }
    else
    { // check using _mp
      vec_mp orig_mp;
      init_point_mp2(orig_mp, 0, Pt_prec); 

      regen_dehom_hom(NULL, orig_mp, &out_prec, NULL, Pt_mp->point, Pt_prec, regen, regen);

      if (infNormVec_mp(orig_mp) < maxNorm)
        finite = 1;
      else
        finite = 0;

      clear_point_mp(orig_mp);
    }
  }

  return finite;
}

int determineRegenHigherDim(double tol, double ratio, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, regen_t *regen, int level_num)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether it satisfies later functions or not    *
* NOTES:                                                        *
\***************************************************************/
{
  int higherDim = 0, endFunc = regen->level[level_num].level + regen->level[level_num].depth;

  // see if there is a possibility of being on a higher dim component
  if (endFunc == regen->num_funcs)
  { // no 'extra' functions
    higherDim = 0;
  }
  else
  { // there are 'extra' functions to test
    int i, *isZero = (int *)bmalloc(regen->num_funcs * sizeof(int));

    if (Pt_prec < 64)
    { // check using _d
      vec_d orig_d, orig_last_d;
      eval_struct_d e;
      init_vec_d(orig_d, 0);
      init_vec_d(orig_last_d, 0);
      init_eval_struct_d(e, 0, 0, 0);

      // setup orig
      if (regen->level[level_num].useIntrinsicSlice)
      { // convert to original coordinates first
        intrinsicToExtrinsic_d(orig_d, Pt_d->point, regen->level[level_num].B_d, regen->level[level_num].p_d);
      }
      else
      { // copy 
        point_cp_d(orig_d, Pt_d->point);
      }
      // remove extra coordinate
      orig_d->size--;

      // evaluate
      if (regen->noChanges)
        eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_d, Pt_d->time, regen->square_d, regen->square_eval_d);
      else
        regen_square_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_d, Pt_d->time, regen);

      // setup orig_last
      if (last_approx_prec < 64)
      { // compute orig_last_d
        if (regen->level[level_num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_d(orig_last_d, last_approx_d, regen->level[level_num].B_d, regen->level[level_num].p_d);
        }
        else
        { // copy
          point_cp_d(orig_last_d, last_approx_d);
        }
        // remove extra coordinate
        orig_last_d->size--;
      }        
      else
      { // convert to _d
        point_mp_to_d(orig_last_d, last_approx_mp);

        if (regen->level[level_num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_d(orig_last_d, orig_last_d, regen->level[level_num].B_d, regen->level[level_num].p_d);
        }
        // remove extra coordinate
        orig_last_d->size--;
      } 

      // evaluate
      if (regen->noChanges)
        eval_d(orig_d, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_d, Pt_d->time, regen->square_d, regen->square_eval_d);
      else
        regen_square_eval_d(orig_d, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_d, Pt_d->time, regen);

      // compare
      nonsolutions_check_compare_d(isZero, e.funcVals, orig_d, endFunc, regen->num_funcs, tol, ratio);

      // check the functions below
      for (i = endFunc; i < regen->num_funcs; i++)
        if (isZero[i])
          higherDim = 1;

      // clear
      clear_vec_d(orig_d);
      clear_vec_d(orig_last_d);
      clear_eval_struct_d(e);
    }
    else
    { // check using _mp
      vec_mp orig_mp, orig_last_mp;
      eval_struct_mp e;
      init_vec_mp2(orig_mp, 0, Pt_prec);
      init_vec_mp2(orig_last_mp, 0, Pt_prec);
      init_eval_struct_mp2(e, 0, 0, 0, Pt_prec);

      // setup orig
      if (regen->level[level_num].useIntrinsicSlice)
      { // convert to original coordinates first
        intrinsicToExtrinsic_mp(orig_mp, Pt_mp->point, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
      }
      else
      { // copy
        point_cp_mp(orig_mp, Pt_mp->point);
      }
      // remove extra coordinate
      orig_mp->size--;

      // evaluate
      if (regen->noChanges)
        eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_mp, Pt_mp->time, regen->square_mp, regen->square_eval_mp);
      else
        regen_square_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_mp, Pt_mp->time, regen);

      // setup orig_last
      if (last_approx_prec < 64)
      { // convert to _mp 
        point_d_to_mp(orig_last_mp, last_approx_d);

        if (regen->level[level_num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_mp(orig_last_mp, orig_last_mp, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
        }
        // remove extra coordinate
        orig_last_mp->size--;
      }
      else
      { // compute orig_last 
        if (regen->level[level_num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_mp(orig_last_mp, last_approx_mp, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
        }
        else
        { // copy
          point_cp_mp(orig_last_mp, last_approx_mp);
        }
        // remove extra coordinate
        orig_last_mp->size--;
      }

      // evaluate
      if (regen->noChanges)
        eval_mp(orig_mp, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_mp, Pt_mp->time, regen->square_mp, regen->square_eval_mp);
      else
        regen_square_eval_mp(orig_mp, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_mp, Pt_mp->time, regen);

      // compare
      nonsolutions_check_compare_mp(isZero, e.funcVals, orig_mp, endFunc, regen->num_funcs, tol, ratio);

      // check the functions below
      for (i = endFunc; i < regen->num_funcs; i++)
        if (isZero[i])
          higherDim = 1;

      clear_vec_mp(orig_mp);
      clear_vec_mp(orig_last_mp);
      clear_eval_struct_mp(e);
    }

    free(isZero);
  }

  return higherDim;
}

int regenMoveOrigPt(int *pathNum, int totalPaths, int pathMod, regen_t *regen, tracker_config_t *T, int **new_linear, int **new_linear_degree, int new_count, mat_d B_transpose_d, mat_mp B_transpose_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *START, FILE *NEXTSTARTPTS, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths moved to NEXTSTARTPTS          *
* NOTES: tracks the paths to prepare the next level             *
\***************************************************************/
{
  int i, j, prec, num_input_vars, num_output_vars, next_level_num = level_num + 1;
  int linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  int new_linearSize = regen->level[next_level_num].level + regen->level[next_level_num].depth;
  int input_IntrinsicSlice = regen->level[level_num].useIntrinsicSlice;
  int curr_linear, curr_linear_degree; // ignore them!
  int num_bad, num_moved = 0;
  int initial_pathNum = *pathNum;
  int *temp_linear = (int *)bmalloc(new_linearSize * sizeof(int)), *temp_linear_degree = (int *)bmalloc(new_linearSize * sizeof(int));
  FILE *TEMPSTARTPTS = fopen("regenTempFile", "w+");
  int *rV = (int *)bmalloc(new_count * sizeof(int));
  for (i = 0; i < new_count; i++)
    rV[i] = 0;

  // top of temp_linear & temp_linear_degree are all 0s
  for (i = 0; i < linearSize; i++)
    temp_linear[i] = temp_linear_degree[i] = 0;

  // setup the number of variables that the start points have
  if (input_IntrinsicSlice)
    num_input_vars = linearSize;
  else
    num_input_vars = regen->num_variables;

  if (regen->level[next_level_num].useIntrinsicSlice)
    num_output_vars = new_linearSize;
  else
    num_output_vars = regen->num_variables;

  if (T->MPType == 0)
  { // setup using double precision
    comp_d tempComp;
    point_data_d origPt;
    init_point_data_d(&origPt, num_input_vars);
    origPt.point->size = num_input_vars;
    set_one_d(origPt.time);

    // read in original point
    fscanf(START, "%d\n%d\n", &i, &prec);
    for (i = 0; i < num_input_vars; i++)
      fscanf(START, "%lf %lf\n", &origPt.point->coord[i].r, &origPt.point->coord[i].i);
    for (i = 0; i < linearSize; i++)
      fscanf(START, "%d %d\n", &curr_linear, &curr_linear_degree);

    if (input_IntrinsicSlice)
    { // convert to extrinsic coordinates
      intrinsicToExtrinsic_d(origPt.point, origPt.point, regen->level[level_num].B_d, regen->level[level_num].p_d);
    }

    // move it to the new patch
    move_to_patch_vec_d(origPt.point, origPt.point, regen->main_homVar_d);

    // see if we need to convert to intrinsic coordinates for next level
    if (regen->level[next_level_num].useIntrinsicSlice)
    { // setup the point using intrinsic coordinates
      extrinsicToIntrinsic_d(origPt.point, origPt.point, B_transpose_d, regen->level[next_level_num].p_d);
    }

    // loop over the new points
    for (i = 0; i < new_count; i++)
    { // setup temp linear & degree
      for (j = linearSize; j < new_linearSize; j++)
      {
        temp_linear[j] = new_linear[i][j - linearSize];
        temp_linear_degree[j] = new_linear_degree[i][j - linearSize];
      }

      // see if we need to track this path
      if (i > 0 || regen->num_var_gps > 1)
      { // we need to track this path
        if (pathMod > 0 && !(*pathNum % pathMod))
          printf("Moving %d of %d\n", *pathNum, totalPaths);

        rV[i] = regenLinearTrackPath(*pathNum, &origPt, NULL, temp_linear, temp_linear_degree, regen, T, OUT, RAWOUT, MIDOUT, FAIL, TEMPSTARTPTS, next_level_num, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment path number
        (*pathNum)++;
      }
    }

    if (regen->num_var_gps == 1)
    { // there should be no errors on the paths actually tracked!
      num_bad = 0;
      for (i = 1; i < new_count; i++)
      {
        if (rV[i] != 0 && rV[i] != retVal_reached_minTrackT)
        { // store this path to bad_paths
          printf("WARNING: When preparing path %d in level %d, it had retVal %d.\n", initial_pathNum + i, next_level_num, rV[i]);
          num_bad++;
        }
      }
      if (num_bad == 0)
      { // print points to NEXTSTARTPTS
        rewind(TEMPSTARTPTS);
        // setup first one
        for (j = linearSize; j < new_linearSize; j++)
        {
          temp_linear[j] = new_linear[0][j - linearSize];
          temp_linear_degree[j] = new_linear_degree[0][j - linearSize]; // should be 0
        }
        printPointLinearDegree(NEXTSTARTPTS, origPt.point, NULL, prec, temp_linear, temp_linear_degree, new_linearSize);

        // copy the rest
        for (i = 1; i < new_count; i++)
        { // read in path number
          fscanf(TEMPSTARTPTS, "%d\n", &j);
          // copy the coordinates
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp->r, &tempComp->i);
            // print coordinate j
            fprintf(NEXTSTARTPTS, "%.15e %.15e;\n", tempComp->r, tempComp->i);
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            // print linear & degree for function j
            fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
          } 
          fprintf(NEXTSTARTPTS, "\n");
        }
        num_moved = new_count;
      }
      else
      { // do not move any
        num_moved = 0;
      }
    }
    else
    { // move only those that have a good rV
      num_moved = 0;
      rewind(TEMPSTARTPTS);
      for (i = 0; i < new_count; i++)
      { // read in the path number
        fscanf(TEMPSTARTPTS, "%d\n", &j);

        if (rV[j - initial_pathNum] == 0)
        { // copy this path to NEXTSTARTPTS
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp->r, &tempComp->i);
            // print coordinate j
            print_d(NEXTSTARTPTS, 0, tempComp);
            fprintf(NEXTSTARTPTS, ";\n");
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            // print linear & degree for function j
            fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
          }
          fprintf(NEXTSTARTPTS, "\n");
          num_moved++;
        }
        else
        { // move past this one
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp->r, &tempComp->i);
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
          }
        }
      }
    }

    // clear memory 
    clear_point_data_d(&origPt);
  }
  else if (T->MPType == 1)
  { // setup using multi precision
    comp_mp tempComp;
    point_data_mp origPt;
    init_mp(tempComp);
    init_point_data_mp(&origPt, num_input_vars);
    origPt.point->size = num_input_vars;
    set_one_mp(origPt.time);

    // read in original point
    fscanf(START, "%d\n%d\n", &i, &prec);
    for (i = 0; i < num_input_vars; i++)
    {
      mpf_inp_str(origPt.point->coord[i].r, START, 10);
      mpf_inp_str(origPt.point->coord[i].i, START, 10);
    }
    for (i = 0; i < linearSize; i++)
      fscanf(START, "%d %d\n", &curr_linear, &curr_linear_degree);

    if (input_IntrinsicSlice)
    { // convert to extrinsic coordinates
      intrinsicToExtrinsic_mp(origPt.point, origPt.point, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
    }

    // move it to the new patch
    move_to_patch_vec_mp(origPt.point, origPt.point, regen->main_homVar_mp);

    // see if we need to convert to intrinsic coordinates for next level
    if (regen->level[next_level_num].useIntrinsicSlice)
    { // setup the point using intrinsic coordinates
      extrinsicToIntrinsic_mp(origPt.point, origPt.point, B_transpose_mp, regen->level[next_level_num].p_mp);
    }

    // loop over the new points
    for (i = 0; i < new_count; i++)
    { // setup temp linear & degree
      for (j = linearSize; j < new_linearSize; j++)
      {
        temp_linear[j] = new_linear[i][j - linearSize];
        temp_linear_degree[j] = new_linear_degree[i][j - linearSize];
      }
  
      // see if we need to track this path
      if (i > 0 || regen->num_var_gps > 1)
      { // we need to track this path
        if (pathMod > 0 && !(*pathNum % pathMod))
          printf("Moving %d of %d\n", *pathNum, totalPaths);
  
        rV[i] = regenLinearTrackPath(*pathNum, NULL, &origPt, temp_linear, temp_linear_degree, regen, T, OUT, RAWOUT, MIDOUT, FAIL, TEMPSTARTPTS, next_level_num, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  
        // increment path number
        (*pathNum)++;
      }
    }

    if (regen->num_var_gps == 1)
    { // there should be no errors on the paths actually tracked!
      num_bad = 0;
      for (i = 1; i < new_count; i++)
      {
        if (rV[i] != 0 && rV[i] != retVal_reached_minTrackT)
        { // store this path to bad_paths
          printf("WARNING: When preparing path %d in level %d, it had retVal %d.\n", initial_pathNum + i, next_level_num, rV[i]);
          num_bad++;
        }
      }
      if (num_bad == 0)
      { // print points to NEXTSTARTPTS
        rewind(TEMPSTARTPTS);
        // setup first one
        for (j = linearSize; j < new_linearSize; j++)
        {
          temp_linear[j] = new_linear[0][j - linearSize];
          temp_linear_degree[j] = new_linear_degree[0][j - linearSize]; // should be 0
        }
        printPointLinearDegree(NEXTSTARTPTS, NULL, origPt.point, prec, temp_linear, temp_linear_degree, new_linearSize);

        // copy the rest
        for (i = 1; i < new_count; i++)
        { // read in path number
          fscanf(TEMPSTARTPTS, "%d\n", &j);
          // copy the coordinates
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            mpf_inp_str(tempComp->r, TEMPSTARTPTS, 10);
            mpf_inp_str(tempComp->i, TEMPSTARTPTS, 10);
            // print coordinate j
            print_mp(NEXTSTARTPTS, 0, tempComp);
            fprintf(NEXTSTARTPTS, ";\n");
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            // print linear & degree for function j
            fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
          }
          fprintf(NEXTSTARTPTS, "\n");
        }
        num_moved = new_count;
      }
      else
      { // do not move any
        num_moved = 0;
      }
    }
    else
    { // move only those that have a good rV
      num_moved = 0;
      rewind(TEMPSTARTPTS);
      for (i = 0; i < new_count; i++)
      { // read in the path number
        fscanf(TEMPSTARTPTS, "%d\n", &j);

        if (rV[j - initial_pathNum] == 0)
        { // copy this path to NEXTSTARTPTS
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            mpf_inp_str(tempComp->r, TEMPSTARTPTS, 10);
            mpf_inp_str(tempComp->i, TEMPSTARTPTS, 10);
            // print coordinate j
            print_mp(NEXTSTARTPTS, 0, tempComp);
            fprintf(NEXTSTARTPTS, ";\n");
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            // print linear & degree for function j
            fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
          }
          fprintf(NEXTSTARTPTS, "\n");
          num_moved++;
        }
        else
        { // move past this one
          for (j = 0; j < num_output_vars; j++)
          { // read in jth coordinate
            mpf_inp_str(tempComp->r, TEMPSTARTPTS, 10);
            mpf_inp_str(tempComp->i, TEMPSTARTPTS, 10);
          }
          for (j = 0; j < new_linearSize; j++)
          { // read in linear & degree for function j
            fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
          }
        }
      }
    }

    // clear memory 
    clear_mp(tempComp);
    clear_point_data_mp(&origPt);
  }
  else
  { // setup using adaptive precision
    comp_d tempComp_d;
    comp_mp tempComp_mp;
    point_data_d origPt_d;
    point_data_mp origPt_mp;
    init_mp(tempComp_mp);
    init_point_data_d(&origPt_d, num_input_vars);
    init_point_data_mp(&origPt_mp, num_input_vars);
    origPt_d.point->size = origPt_mp.point->size = num_input_vars;

    // read in precision
    fscanf(START, "%d\n%d\n", &i, &prec);
    // read in point
    if (prec < 64)
    { // read in using double precision
      for (i = 0; i < num_input_vars; i++)
        fscanf(START, "%lf %lf\n", &origPt_d.point->coord[i].r, &origPt_d.point->coord[i].i);
    }
    else
    { // read in using multi precision
      change_prec(regen, prec);
      change_prec_mp(tempComp_mp, prec);
      change_prec_point_data_mp(&origPt_mp, prec);

      for (i = 0; i < num_input_vars; i++)
      {
        mpf_inp_str(origPt_mp.point->coord[i].r, START, 10);
        mpf_inp_str(origPt_mp.point->coord[i].i, START, 10);
      }
    }
    set_one_d(origPt_d.time);
    set_one_mp(origPt_mp.time);

    // move past slices
    for (i = 0; i < linearSize; i++)
      fscanf(START, "%d %d\n", &curr_linear, &curr_linear_degree);

    if (input_IntrinsicSlice)
    { // convert to extrinsic coordinates
      if (prec < 64)
        intrinsicToExtrinsic_d(origPt_d.point, origPt_d.point, regen->level[level_num].B_d, regen->level[level_num].p_d);
      else
        intrinsicToExtrinsic_mp(origPt_mp.point, origPt_mp.point, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
    }

    // move it to the new patch
    if (prec < 64)
      move_to_patch_vec_d(origPt_d.point, origPt_d.point, regen->main_homVar_d);
    else
      move_to_patch_vec_mp(origPt_mp.point, origPt_mp.point, regen->main_homVar_mp);

    // see if we need to convert to intrinsic coordinates for next level
    if (regen->level[next_level_num].useIntrinsicSlice)
    { // setup the point using intrinsic coordinates
      if (prec < 64)
        extrinsicToIntrinsic_d(origPt_d.point, origPt_d.point, B_transpose_d, regen->level[next_level_num].p_d);
      else
      { // setup B_transpose_mp properly
        change_prec_mat_mp(B_transpose_mp, prec);
        transpose_mp(B_transpose_mp, regen->level[next_level_num].B_mp);

        extrinsicToIntrinsic_mp(origPt_mp.point, origPt_mp.point, B_transpose_mp, regen->level[next_level_num].p_mp);
      }
    }

    // setup origPt_d if prec >= 64
    if (prec >= 64)
      convert_point_data_mp_to_d(&origPt_d, &origPt_mp);

    // loop over the new points
    for (i = 0; i < new_count; i++)
    { // setup temp linear & degree
      for (j = linearSize; j < new_linearSize; j++)
      {
        temp_linear[j] = new_linear[i][j - linearSize];
        temp_linear_degree[j] = new_linear_degree[i][j - linearSize];
      }

      // see if we need to track this path
      if (i > 0 || regen->num_var_gps > 1)
      { // we need to track this path
        if (pathMod > 0 && !(*pathNum % pathMod))
          printf("Moving %d of %d\n", *pathNum, totalPaths);

        rV[i] = regenLinearTrackPath(*pathNum, &origPt_d, NULL, temp_linear, temp_linear_degree, regen, T, OUT, RAWOUT, MIDOUT, FAIL, TEMPSTARTPTS, next_level_num, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment path number
        (*pathNum)++;
      }
    }

    if (regen->num_var_gps == 1)
    { // there should be no errors on the paths actually tracked!
      num_bad = 0;
      for (i = 1; i < new_count; i++)
      {
        if (rV[i] != 0 && rV[i] != retVal_reached_minTrackT)
        { // store this path to bad_paths
          printf("WARNING: When preparing path %d in level %d, it had retVal %d.\n", initial_pathNum + i, next_level_num, rV[i]);
          num_bad++;
        }
      }
      if (num_bad == 0)
      { // print points to NEXTSTARTPTS
        rewind(TEMPSTARTPTS);
        // setup first one
        for (j = linearSize; j < new_linearSize; j++)
        {
          temp_linear[j] = new_linear[0][j - linearSize];
          temp_linear_degree[j] = new_linear_degree[0][j - linearSize]; // should be 0
        }
        printPointLinearDegree(NEXTSTARTPTS, origPt_d.point, origPt_mp.point, prec, temp_linear, temp_linear_degree, new_linearSize);

        // copy the rest
        if (prec < 64)
        { // copy using double precision
          for (i = 1; i < new_count; i++)
          { // read in path number
            fscanf(TEMPSTARTPTS, "%d\n", &j);
            // copy the coordinates
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
              // print coordinate j
              fprintf(NEXTSTARTPTS, "%.15e %.15e;\n", tempComp_d->r, tempComp_d->i);
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
              // print linear & degree for function j
              fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
            }
            fprintf(NEXTSTARTPTS, "\n");
          }
          num_moved = new_count;
        }
        else 
        { // copy using mulitprecision
          for (i = 1; i < new_count; i++)
          { // read in path number
            fscanf(TEMPSTARTPTS, "%d\n", &j);
            // copy the coordinates
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              mpf_inp_str(tempComp_mp->r, TEMPSTARTPTS, 10);
              mpf_inp_str(tempComp_mp->i, TEMPSTARTPTS, 10);
              // print coordinate j
              print_mp(NEXTSTARTPTS, 0, tempComp_mp);
              fprintf(NEXTSTARTPTS, ";\n");
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
              // print linear & degree for function j
              fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
            }
            fprintf(NEXTSTARTPTS, "\n");
          }
          num_moved = new_count;
        }
      }
      else
      { // do not move any
        num_moved = 0;
      }
    }
    else
    { // move only those that have a good rV
      num_moved = 0;
      rewind(TEMPSTARTPTS);
      if (prec < 64)
      { // move using double precision
        for (i = 0; i < new_count; i++)
        { // read in the path number
          fscanf(TEMPSTARTPTS, "%d\n", &j);

          if (rV[j - initial_pathNum] == 0)
          { // copy this path to NEXTSTARTPTS
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
              // print coordinate j
              print_d(NEXTSTARTPTS, 0, tempComp_d);
              fprintf(NEXTSTARTPTS, ";\n");
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
              // print linear & degree for function j
              fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
            }
            fprintf(NEXTSTARTPTS, "\n");
            num_moved++;
          }
          else
          { // move past this one
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            }
          }
        }
      }
      else
      { // move using multiprecision
        for (i = 0; i < new_count; i++)
        { // read in the path number
          fscanf(TEMPSTARTPTS, "%d\n", &j);

          if (rV[j - initial_pathNum] == 0)
          { // copy this path to NEXTSTARTPTS
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              mpf_inp_str(tempComp_mp->r, TEMPSTARTPTS, 10);
              mpf_inp_str(tempComp_mp->i, TEMPSTARTPTS, 10);
              // print coordinate j
              print_mp(NEXTSTARTPTS, 0, tempComp_mp);
              fprintf(NEXTSTARTPTS, ";\n");
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
              // print linear & degree for function j
              fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
            }
            fprintf(NEXTSTARTPTS, "\n");
            num_moved++;
          }
          else
          { // move past this one
            for (j = 0; j < num_output_vars; j++)
            { // read in jth coordinate
              mpf_inp_str(tempComp_mp->r, TEMPSTARTPTS, 10);
              mpf_inp_str(tempComp_mp->i, TEMPSTARTPTS, 10);
            }
            for (j = 0; j < new_linearSize; j++)
            { // read in linear & degree for function j
              fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
            }
          }
        }
      }
    }

    clear_mp(tempComp_mp);
    clear_point_data_d(&origPt_d);
    clear_point_data_mp(&origPt_mp);
  }

  // close & delete file
  fclose(TEMPSTARTPTS);
  remove("regenTempFile");

  free(rV);
  free(temp_linear);
  free(temp_linear_degree);

  return num_moved;
}

int regenPrepareNextLevel(int pathMod, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *START, char *nextStartFile, int curr_level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths actually tracked - for MIDOUT  *
* NOTES: prepares the next regen level                          *
\***************************************************************/
{
  int i, j, k, l, count, var_gp, numNewPts = 0, numOrigPts = regen->level[curr_level_num].num_nonsing, next_level_num = curr_level_num + 1;
  int startFunc = regen->level[next_level_num].level, depth = regen->level[next_level_num].depth;
  int pathNum, pathsTracked = 0, endFunc = startFunc + depth;
  int **new_linear = NULL, **new_linear_degree = NULL;
  double basicTol, endgameTol, finalTol;
  mat_d B_transpose_d;
  mat_mp B_transpose_mp;
  FILE *NEXTSTARTPTS = fopen(nextStartFile, "w");
  int good_size = 0, **goodLoc = NULL;
  int (*find_dehom_temp)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);

  // initialize
  init_mat_d(B_transpose_d, 0, 0);
  init_mat_mp(B_transpose_mp, 0, 0);

  // setup a new patch
  if (T->MPType == 0)
  { // create new patch coeff in _d
    make_vec_random_d(regen->main_homVar_d, regen->num_variables);
  }
  else if (T->MPType == 1)
  { // create new patch coeff in _mp
    make_vec_random_mp(regen->main_homVar_mp, regen->num_variables);
  }
  else
  { // create new patch coeff in _d, _mp, _rat
    make_vec_random_rat(regen->main_homVar_d, regen->main_homVar_mp, regen->main_homVar_rat, regen->num_variables, regen->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup the structures for the next level
  setupRegenLevelStructures(regen, T->MPType, next_level_num, T->AMP_max_prec);

  // setup B_transpose, if needed
  if (regen->level[next_level_num].useIntrinsicSlice)
  {
    if (T->MPType == 0 || T->MPType == 2)
    { // setup B_transpose_d
      transpose_d(B_transpose_d, regen->level[next_level_num].B_d);
    }
    if (T->MPType == 1 || T->MPType == 2)
    { // setup B_transpose_mp
      transpose_mp(B_transpose_mp, regen->level[next_level_num].B_mp);
    }
  }

  // since we are using slices across different variable groups, the linears mean nothing!
  // setup the possible decompositions for the next move
  createDecomp(&goodLoc, &good_size, regen, startFunc, endFunc);

  // compute number of paths per original point
  count = 0;
  for (i = 0; i < good_size; i++)
  {
    numNewPts = 1;
    for (j = 0; j < depth; j++)
    {
      var_gp = goodLoc[i][j];
      numNewPts *= regen->degrees[startFunc + j][var_gp];
    }

    count += numNewPts;
  }
  // NOTE: count * numOrigPts bounds the number of paths for the next level
  // Need to track all of them if m-hom and (count - 1) * numOrigPts if 1-hom

  if (regen->num_var_gps == 1)
    pathsTracked = (count - 1) * numOrigPts;
  else
    pathsTracked = count * numOrigPts;

  // setup the linear/degree structures
  new_linear = (int **)bmalloc(count * sizeof(int *));
  new_linear_degree = (int **)bmalloc(count * sizeof(int *));
  for (i = 0; i < count; i++)
  {
    new_linear[i] = (int *)bmalloc(depth * sizeof(int));
    new_linear_degree[i] = (int *)bmalloc(depth * sizeof(int));
  }

  l = 0;
  for (i = 0; i < good_size; i++)
  {
    numNewPts = 1;
    for (j = 0; j < depth; j++)
    {
      var_gp = goodLoc[i][j];
      numNewPts *= regen->degrees[startFunc + j][var_gp];
    }

    // setup the first new_linear & new_linear_degree for this one
    for (k = 0; k < depth; k++)
    {
      new_linear[l][k] = goodLoc[i][k];
      new_linear_degree[l][k] = 0;
    }
    // increment l
    l++;

    // setup other ones by increment new_linear_degree
    for (j = 1; j < numNewPts; j++)
    { // copy new_linear & new_linear_degree from previous
      for (k = 0; k < depth; k++)
      { 
        new_linear[l][k] = new_linear[l-1][k];
        new_linear_degree[l][k] = new_linear_degree[l-1][k];
      }
      // move degree to next one
      for (k = depth - 1; k >= 0; k--)
      {
        var_gp = new_linear[l][k];
        if (new_linear_degree[l][k] + 1 == regen->degrees[startFunc + k][var_gp])
        { // set back to 0
          new_linear_degree[l][k] = 0;
        }
        else
        { // increment and exit
          new_linear_degree[l][k]++;
          k = -10;
        }
      }
      // increment l
      l++;
    }
  }

  if (regen->num_var_gps == 1)
  { // all slices are the same general setup, so every path should track correctly!
    find_dehom_temp = &zero_dehom;
  }
  else
  { // points could go to infinity
    find_dehom_temp = find_dehom;
  }

  fprintf(NEXTSTARTPTS, "                                                   \n\n");

  // print header for level
  printf("\nPreparing regeneration level %d of %d: %d witness point%s to move.\n", next_level_num, regen->num_levels, pathsTracked, pathsTracked == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Moving linears to regeneration level %d.\n", next_level_num);
  fprintf(OUT, "*****************************************************\n");

  // setup the correct tracking tolerances
  basicTol = T->basicNewtonTol;
  endgameTol = T->endgameNewtonTol;
  finalTol = T->final_tolerance;
  T->basicNewtonTol = T->sliceBasicNewtonTol;
  T->endgameNewtonTol = T->sliceEndgameNewtonTol;
  T->final_tolerance = T->sliceFinalTol;

  // loop over the original paths
  pathNum = numNewPts = 0;
  for (i = 0; i < numOrigPts; i++)
  { // prepare the set of points associated with this original point
    numNewPts += regenMoveOrigPt(&pathNum, pathsTracked, pathMod, regen, T, new_linear, new_linear_degree, count, B_transpose_d, B_transpose_mp, OUT, RAWOUT, MIDOUT, FAIL, START, NEXTSTARTPTS, curr_level_num, trackCount, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom_temp);
  }

  // store the number of new paths for next level
  regen->level[next_level_num].num_paths = numNewPts;

  // restore the correct tracking tolerances
  T->basicNewtonTol = basicTol;
  T->endgameNewtonTol = endgameTol;
  T->final_tolerance = finalTol;

  // print the relevant data to NEXTSTARTPTS so that it can be used to rerun this exact same problem
  printRegenRelevantData(regen, T->MPType, next_level_num, NEXTSTARTPTS);

  // rewind to print number of new points
  rewind(NEXTSTARTPTS);
  fprintf(NEXTSTARTPTS, "%d", numNewPts);

  // close files
  fclose(NEXTSTARTPTS);

  // clear memory
  clear_mat_d(B_transpose_d);
  clear_mat_mp(B_transpose_mp);
  for (i = good_size - 1; i >= 0; i--)
    free(goodLoc[i]);
  free(goodLoc);
  for (i = count - 1; i >= 0; i--)
  {
    free(new_linear[i]);
    free(new_linear_degree[i]);
  }
  free(new_linear);
  free(new_linear_degree);

  return pathsTracked;
}

int regenLinearTrackPath(int path_num, point_data_d *startPt_d, point_data_mp *startPt_mp, int *curr_linear, int *curr_linear_degree, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: retVal of the path                             *
* NOTES: moves the path to the next level                       *
\***************************************************************/
{
  int retVal, rankType = 0, rankDef = 0, finite = 1, higherDim = 0; 
  double smallest, largest;
  endgame_data_t endPt;
  point_d dehom_d, orig_d, orig_last_d;
  point_mp dehom_mp, orig_mp, orig_last_mp;

  init_endgame_data(&endPt, T->Precision);
  init_point_d(dehom_d, 0); init_point_d(orig_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0); init_point_mp(orig_last_mp, 0);

  // seutp for tracking the path
  regen->curr_level_num = level_num;
  regen->curr_linear = curr_linear;
  regen->curr_linear_degree = curr_linear_degree;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision  - find rankDef

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, &rankDef, NULL, &smallest, &largest, &endPt, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision - find rankDef

    // print the path header for the point
    printPathHeader_mp(OUT, startPt_mp, T, path_num, regen, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_rank_mp(path_num, rankType, &rankDef, NULL, &smallest, &largest, &endPt, startPt_mp, OUT, MIDOUT, T, regen, ptr_to_eval_mp, find_dehom);
  }
  else
  { // track the path using AMP

    // print the header for the point
    printPathHeader_d(OUT, startPt_d, T, path_num, regen, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(path_num, &endPt, startPt_d, OUT, MIDOUT, T, regen, regen, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_last_d, orig_last_mp, &endPt.last_approx_prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, regen, regen);
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &endPt.prec, endPt.PD_d.point, endPt.PD_mp.point, endPt.prec, regen, regen);

  // compute the retVal
  endPt.retVal = computeRegenRetval(regen, level_num, &endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, 0, NULL);

  // only need to check when using multi homogeneous
  if (endPt.retVal == 0 && regen->num_var_gps > 1)
  { // we need to check the endpoint
    regenSortEndpoint_basic(1, &rankDef, &finite, &higherDim, endPt.condition_number, regen, level_num, path_num, T, OUT, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, curr_linear, curr_linear_degree, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // print the footer to OUT and print all info to RAWOUT & RAWSORT
  retVal = printRegenLinearTrackFooter(endPt.retVal, regen, level_num, path_num, &endPt, rankDef, finite, higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, NEXTSTARTPTS, FAIL, T);

  // clear memory
  clear_endgame_data(&endPt);
  clear_point_d(dehom_d); clear_point_d(orig_d); clear_point_d(orig_last_d);
  clear_point_d(dehom_mp); clear_point_mp(orig_mp); clear_point_mp(orig_last_mp);

  return retVal;
}

void regenRemoveBadPaths(FILE *NEXTSTARTPTS, FILE *INPUT, int *num_new_paths, int num_funcs, int num_vars, int MPType, int numOrigPts, int numMovePts, int *startPathNum, int num_bad, int *badPaths)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: removes all of the bad paths                           *
\***************************************************************/
{
  int i, j, linear, degree, num_curr_paths = numOrigPts + numMovePts;
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
          fscanf(INPUT, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
          // print coordinate j
          fprintf(NEXTSTARTPTS, "%.15e %.15e;\n", tempComp_d->r, tempComp_d->i);
        }
        for (j = 0; j < num_funcs; j++)
        { // read in linear & degree for function j
          fscanf(INPUT, "%d %d\n", &linear, &degree);
          // print linear & degree for function j
          fprintf(NEXTSTARTPTS, "%d %d\n", linear, degree);
        }
        fscanf(INPUT, "\n");
        fprintf(NEXTSTARTPTS, "\n");
      }
      else
      { // move past this one
        for (j = 0; j < num_vars; j++)
        { // read in coordinate j
          fscanf(INPUT, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
        }
        for (j = 0; j < num_funcs; j++)
        { // read in linear & degree for function j
          fscanf(INPUT, "%d %d\n", &linear, &degree);
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
        for (j = 0; j < num_funcs; j++)
        { // read in linear & degree for function j
          fscanf(INPUT, "%d %d\n", &linear, &degree);
          // print linear & degree for function j
          fprintf(NEXTSTARTPTS, "%d %d\n", linear, degree);
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
        for (j = 0; j < num_funcs; j++)
        { // read in linear & degree for function j
          fscanf(INPUT, "%d %d\n", &linear, &degree);
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

void regenOutputChart(regen_t *regen, FILE *fp, int infRemoved)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int level, discarded = 0, total_paths = 0, num_levels = regen->num_levels;

  fprintf(fp, "\n\n*************** Regeneration summary ***************\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|level|   paths   | nonsing endpoints | total discarded | sing endpoints |");
  if (infRemoved)
    fprintf(fp, " inf endpoints |");
  fprintf(fp, " higher dim'l | other bad endpoints\n");
  fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");

  for (level = 0; level < num_levels; level++)
  { // add to the total paths tracked
    total_paths += regen->level[level].num_paths;

    // find the number discarded at this level
    discarded = regen->level[level].num_sing + regen->level[level].num_higher_dim + regen->level[level].num_bad;
    if (infRemoved)
      discarded += regen->level[level].num_inf;

    fprintf(fp, "| %-4d|   %-8d| %-18d| %-16d| %-15d|", regen->level[level].level, regen->level[level].num_paths, regen->level[level].num_nonsing, discarded, regen->level[level].num_sing);
    if (infRemoved)
      fprintf(fp, " %-14d|", regen->level[level].num_inf);
    fprintf(fp, "  %-12d| %d", regen->level[level].num_higher_dim, regen->level[level].num_bad);

    if (level + 1 == num_levels && regen->level[level].num_bad > 0)
      fprintf(fp, " ** see failure summary\n");
    else
      fprintf(fp, "\n");
  }
  fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|total|   %d\n\n", total_paths);

  fprintf(fp, "****************************************************\n\n");

  return;
}

void regenSetupFinalSoln(regen_t *regen, tracker_config_t *T, FILE *NONSINGIN, FILE *RAWIN, FILE *FINALSOLN, mat_d patch_d, mat_mp patch_mp, mpq_t ***patch_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: turns the raw output file from regen into true raw out *
\***************************************************************/
{
  int i, j, pathNum, tempInt, max_prec = 52, num_input_vars, top_level = regen->num_levels - 1, num_funcs = regen->num_funcs, num_output_vars = 0;
  int num_points = regen->level[top_level].num_paths, num_good_points = regen->level[top_level].num_nonsing, useIntrinsic = regen->level[top_level].useIntrinsicSlice;
  comp_d tempComp;
  point_d tempPoint_d;
  point_mp tempPoint_mp;
  double *funcRes = (double *)bmalloc(num_points * sizeof(double)), *condNum = (double *)bmalloc(num_points * sizeof(double)), *newtonRes = (double *)bmalloc(num_points * sizeof(double));
  double *tVal = (double *)bmalloc(num_points * sizeof(double)), *error = (double *)bmalloc(num_points * sizeof(double)), *firstInc = (double *)bmalloc(num_points * sizeof(double));
  int *cycleNum = (int *)bmalloc(num_points * sizeof(int));

  // initialize
  init_point_d(tempPoint_d, 0);
  init_point_mp(tempPoint_mp, 0);

  // setup the number of variables that the points have
  if (useIntrinsic)
    num_input_vars = regen->level[top_level].level + regen->level[top_level].depth;
  else
    num_input_vars = regen->num_variables;

  // setup the number of variables we will print
  for (i = 0; i < regen->PPD.num_hom_var_gp + regen->PPD.num_var_gp; i++)
    num_output_vars += regen->PPD.type[i] + regen->PPD.size[i];

  // read in the data from RAWIN
  for (i = 0; i < num_points; i++)
  { // read in the path number
    pathNum = -1;
    fscanf(RAWIN, "%d\n", &pathNum);
    // verify that pathNum is valid
    if (pathNum < 0 || pathNum >= num_points)
    {
      printf("ERROR: The path number is invalid!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // see if the precision is larger than max_prec
    fscanf(RAWIN, "%d\n", &j);
    if (j > max_prec)
      max_prec = j;

    // move past the coordinates
    for (j = 0; j < num_input_vars; j++)
      fscanf(RAWIN, "%lf%lf\n", &tempComp->r, &tempComp->i);

    // read in the data
    fscanf(RAWIN, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n%d\n%d\n", &funcRes[pathNum], &condNum[pathNum], &newtonRes[pathNum], &tVal[pathNum], &error[pathNum], &firstInc[pathNum], &cycleNum[pathNum], &tempInt, &j);

    // move past the linear & degree info
    for (j = 0; j < num_funcs; j++)
      fscanf(RAWIN, "%d%d\n", &tempInt, &tempInt);
  }

  // print number of variables and dimension to FINALSOLN
  fprintf(FINALSOLN, "%d\n%d\n", num_output_vars, 0);

  // setup patch_mp if using AMP and we have high precision solutions
  if (T->MPType == 2 && max_prec >= 64)
  { // set the precision to max_prec
    change_prec_mat_mp_rat(patch_mp, max_prec, patch_rat);
  }

  // read in the data from NONSINGIN and print to FINALSOLN
  for (i = 0; i < num_good_points; i++)
  { // read in the path number & the precision
    pathNum = -1; 
    fscanf(NONSINGIN, "%d\n%d\n", &pathNum, &tempInt);
    // verify that pathNum is valid
    if (pathNum < 0 || pathNum >= num_points)
    { 
      printf("ERROR: The path number is invalid!\n");
      bexit(ERROR_CONFIGURATION);
    } 

    // print the path number & precision to FINALSOLN
    fprintf(FINALSOLN, "%d\n%d\n", pathNum, tempInt);

    // read in the point and print correctly to FINALSOLN
    if (tempInt < 64)
    { // read in tempPoint_d
      increase_size_point_d(tempPoint_d, num_input_vars);
      tempPoint_d->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
        fscanf(NONSINGIN, "%lf%lf\n", &tempPoint_d->coord[j].r, &tempPoint_d->coord[j].i);

      // convert to original coordinates
      regen->curr_level_num = top_level;
      regen_dehom_hom(tempPoint_d, NULL, &tempInt, tempPoint_d, NULL, tempInt, regen, regen);
   
      // print to FINALSOLN
      for (j = 0; j < tempPoint_d->size; j++)
        fprintf(FINALSOLN, "%.15e %.15e\n", tempPoint_d->coord[j].r, tempPoint_d->coord[j].i);
    }
    else
    { // set the precision on tempPoint_mp
      setprec_point_mp(tempPoint_mp, tempInt);
      increase_size_point_mp(tempPoint_mp, num_input_vars);
      // read in tempPoint_mp
      tempPoint_mp->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(tempPoint_mp->coord[j].r, NONSINGIN, 10);
        mpf_inp_str(tempPoint_mp->coord[j].i, NONSINGIN, 10);
        fscanf(NONSINGIN, "\n"); 
      }

      // convert to original coordinates
      regen->curr_level_num = top_level;
      regen_dehom_hom(NULL, tempPoint_mp, &tempInt, NULL, tempPoint_mp, tempInt, regen, regen);

      // print to FINALSOLN
      for (j = 0; j < tempPoint_mp->size; j++)
      {
        mpf_out_str(FINALSOLN, 10, 0, tempPoint_mp->coord[j].r);
        fprintf(FINALSOLN, " ");
        mpf_out_str(FINALSOLN, 10, 0, tempPoint_mp->coord[j].i);
        fprintf(FINALSOLN, "\n");
      }
    }

    // print the data to FINALSOLN - all of these are successes and so the last number is '1'
    fprintf(FINALSOLN, "%.15e\n%.15e\n%.15e\n%.15e\n%.15e\n%.15e\n%d\n%d\n", funcRes[pathNum], condNum[pathNum], newtonRes[pathNum], tVal[pathNum], error[pathNum], firstInc[pathNum], cycleNum[pathNum], 1);

    // move past the linear & degree info
    for (j = 0; j < num_funcs; j++)
      fscanf(NONSINGIN, "%d %d\n", &tempInt, &tempInt);
  }

  // finish the output to FINALSOLN
  fprintf(FINALSOLN, "%d\n\n", -1);  // bottom of FINALSOLN

  // clear memory
  free(funcRes); 
  free(condNum);
  free(newtonRes);
  free(tVal);
  free(error); 
  free(firstInc);
  free(cycleNum);

  // clear tempPoint
  clear_point_d(tempPoint_d);
  clear_point_mp(tempPoint_mp);

  return;
}

int regen_dehom_hom(point_d hom_d, point_mp hom_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the original hom variables                     *
\***************************************************************/
{
  int j, coord, level_num;
  regen_t *regen = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute hom_d
    comp_d recip;

    regen = (regen_t *)ED_d;
    level_num = regen->curr_level_num;

    if (regen->level[level_num].useIntrinsicSlice)
    { // convert to original coordinates first, normalize by last variable (new dehom) and then compute dehom
      intrinsicToExtrinsic_d(hom_d, in_d, regen->level[level_num].B_d, regen->level[level_num].p_d);
    }
    else
    { // simply copy
      vec_cp_d(hom_d, in_d);
    }

    // adjust size - remove extra variable
    hom_d->size--;
    coord = hom_d->size;

    // normalize by last variable
    if (hom_d->coord[coord].r == 0 && hom_d->coord[coord].i == 0)
    { // generate a random perturbation 
      get_comp_rand_d(recip);
      mul_rdouble_d(recip, recip, 1e-16);
      recip_d(recip, recip);
    }
    else
    { // reciprocate
      recip_d(recip, &hom_d->coord[coord]);
    }
    // dehomogenize by this extra variable
    for (j = 0; j < coord; j++)
      mul_d(&hom_d->coord[j], &hom_d->coord[j], recip);
  }
  else
  { // compute hom_mp
    mpf_t epsilon;
    comp_mp recip;
    mpf_init2(epsilon, *out_prec);
    init_mp2(recip, *out_prec);

    regen = (regen_t *)ED_mp;
    level_num = regen->curr_level_num;

    // set prec on hom_mp
    setprec_point_mp(hom_mp, *out_prec);

    if (regen->level[level_num].useIntrinsicSlice)
    { // convert to original coordinates first and then compute dehom
      intrinsicToExtrinsic_mp(hom_mp, in_mp, regen->level[level_num].B_mp, regen->level[level_num].p_mp);
    }
    else
    { // simply copy
      vec_cp_mp(hom_mp, in_mp);
    }

    // adjust size - remove extra variable
    hom_mp->size--;
    coord = hom_mp->size;

    // normalize by last variable
    if (mpfr_zero_p(hom_mp->coord[coord].r) && mpfr_zero_p(hom_mp->coord[coord].i))
    { // generate a random perturbation 
      get_comp_rand_mp(recip);
      mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(*out_prec), __gmp_default_rounding_mode);
      mpf_ui_div(epsilon, 1, epsilon);
      mul_rmpf_mp(recip, recip, epsilon);
      recip_mp(recip, recip);
    }
    else
    { // reciprocate
      recip_mp(recip, &hom_mp->coord[coord]);
    }
    // dehomogenize by this extra variable
    for (j = 0; j < coord; j++)
      mul_mp(&hom_mp->coord[j], &hom_mp->coord[j], recip);

    mpf_clear(epsilon);
    clear_mp(recip);
  }

  regen = NULL;

  return 0;
}

int regen_dehom2(point_d nonhom_d, point_mp nonhom_mp, point_d hom_d, point_mp hom_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the 2 dehom points for regen                   *
\***************************************************************/
{
  regen_t *regen = NULL;

  // compute hom
  regen_dehom_hom(hom_d, hom_mp, out_prec, in_d, in_mp, in_prec, ED_d, ED_mp);

  if (in_prec < 64)
  { // compute nonhom_d
    regen = (regen_t *)ED_d;

    // dehomogenize normally
    getDehomPoint_d(nonhom_d, hom_d, hom_d->size, &regen->PPD);
  }
  else
  { // compute nonhom_mp
    regen = (regen_t *)ED_mp;

    // dehomogenize normally
    getDehomPoint_mp(nonhom_mp, hom_mp, hom_mp->size, &regen->PPD);
  }

  regen = NULL;

  return 0;
}

int regen_dehom_stack(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stack the orig and dehom points into one vector        *
*  USED FOR PATH TRUNCATION                                     *
\***************************************************************/
{
  int j, coord, size1, size2;

  if (in_prec < 64)
  { // compute out_d
    vec_d tempVec1, tempVec2;
    init_vec_d(tempVec1, 0);
    init_vec_d(tempVec2, 0);

    regen_dehom2(tempVec1, NULL, tempVec2, NULL, out_prec, in_d, in_mp, in_prec, ED_d, ED_mp);

    size1 = tempVec1->size;
    size2 = tempVec2->size;

    // setup out_d
    change_size_vec_d(out_d, size1 + size2);
    out_d->size = size1 + size2;

    for (j = 0; j < size1; j++)
      set_d(&out_d->coord[j], &tempVec1->coord[j]);
    for (j = 0, coord = size1; j < size2; j++, coord++)
      set_d(&out_d->coord[coord], &tempVec2->coord[j]);

    clear_vec_d(tempVec1);
    clear_vec_d(tempVec2);
  }
  else
  { // compute out_mp
    vec_mp tempVec1, tempVec2;
    init_vec_mp(tempVec1, 0);
    init_vec_mp(tempVec2, 0);

    regen_dehom2(NULL, tempVec1, NULL, tempVec2, out_prec, in_d, in_mp, in_prec, ED_d, ED_mp);

    size1 = tempVec1->size;
    size2 = tempVec2->size;

    // setup out_mp
    change_size_vec_mp(out_mp, size1 + size2);
    out_mp->size = size1 + size2;

    for (j = 0; j < size1; j++)
      set_mp(&out_mp->coord[j], &tempVec1->coord[j]);
    for (j = 0, coord = size1; j < size2; j++, coord++)
      set_mp(&out_mp->coord[coord], &tempVec2->coord[j]);

    clear_vec_mp(tempVec1);
    clear_vec_mp(tempVec2);
  }

  return 0;
}

int regen_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  if (in_prec < 64)
  { // compute out_d
    vec_d tempVec;
    init_vec_d(tempVec, 0);

    regen_dehom2(out_d, out_mp, tempVec, NULL, out_prec, in_d, in_mp, in_prec, ED_d, ED_mp);

    clear_vec_d(tempVec);
  }
  else
  { // compute out_mp
    vec_mp tempVec;
    init_vec_mp(tempVec, 0);

    regen_dehom2(out_d, out_mp, NULL, tempVec, out_prec, in_d, in_mp, in_prec, ED_d, ED_mp);

    clear_vec_mp(tempVec);
  }

  return 0;
}


