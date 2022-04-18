// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "regeneration.h"
#include "parallel.h"

#ifdef _HAVE_MPI

void head_regen_track_zero_dim(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does regen tracking for zero dimensional - 'headnode'  *
\***************************************************************/
{
  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  if (T->MPType == 0)
  { // send ED_d - only double
    bcast_basic_eval_data_d(ED_d, T->MPType, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // send ED_mp - only fixed MP
    bcast_basic_eval_data_mp(ED_mp, T->MPType, my_id, headnode);
  }
  else
  { // send ED_d to the workers - using AMP
    bcast_basic_eval_data_amp(ED_d, my_id, headnode);
  }

  // send regen to the workers
  bcast_regen_t(regen, T->MPType, my_id, headnode);

  // now that everything is sent, do the main regeneration tracking
  if (T->MPType == 0)
    head_regen_track(startLevel, regen, T, pathMod, midpoint_tol, startName, OUT, midFile, FAIL, FINALSOLN, trackCount, ED_d->patch.patchCoeff, NULL, NULL, my_id, num_processes, headnode);
  else if (T->MPType == 1)
    head_regen_track(startLevel, regen, T, pathMod, midpoint_tol, startName, OUT, midFile, FAIL, FINALSOLN, trackCount, NULL, ED_mp->patch.patchCoeff, ED_mp->patch.patchCoeff_rat, my_id, num_processes, headnode);
  else
    head_regen_track(startLevel, regen, T, pathMod, midpoint_tol, startName, OUT, midFile, FAIL, FINALSOLN, trackCount, ED_d->patch.patchCoeff, ED_mp->patch.patchCoeff, ED_mp->patch.patchCoeff_rat, my_id, num_processes, headnode);

  return;
}

void head_regen_track(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, mat_d finalPatch_d, mat_mp finalPatch_mp, mpq_t ***finalPatch_rat, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does standard regen tracking - 'headnode'              *
* Assume that everything is setup/sent so that we can track     *
\***************************************************************/
{
  int level, num_crossings = 0, size, num_moved, num_levels = regen->num_levels, minPacketSize = 1, maxPacketSize = 20;
  double midpoint_tol_prepare = T->sliceBasicNewtonTol;
  char *str = NULL, rawTrackFile[] = "rawout_track", rawSortFile[] = "rawout_sort", rawPrepareFile[] = "rawout_prepare", midTrackFile[] = "midout_track", midPrepareFile[] = "midout_prepare";
  FILE *RAWOUT = NULL, *RAWSORT = NULL, *START = NULL, *NONSOLN = NULL;

  if (regen->PPD.num_funcs > num_levels) 
  { // setup NONSOLN since system is overdetermined
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // initialize trackCount
  init_trackingStats(trackCount);

  // send 'startLevel'
  MPI_Bcast(&startLevel, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // send the info regarding the starting level - other levels sent when moving linears
  bcast_regenLevel_t(&regen->level[startLevel], T->MPType, regen->curr_precision, my_id, headnode);

  // main tracking loop
  for (level = startLevel; level < num_levels; level++)
  { // setup RAWOUT for this level
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, level);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, level);
    RAWOUT = fopen(str, "w");

    // setup RAWSORT for this level
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

    // read in the number of start points
    fscanf(START, "%d\n", &num_crossings);

    // error checking
    if (num_crossings != regen->level[level].num_paths)
    {
      printf("ERROR: The number of paths (%d vs %d) described in '%s' is not correct!\n", regen->level[level].num_paths, num_crossings, str);
      bexit(ERROR_INVALID_SIZE);
    }

    // track this level
    head_regenTrackLevel(minPacketSize, maxPacketSize, trackCount, level, pathMod, T, regen, START, OUT, RAWOUT, RAWSORT, FAIL, NONSOLN, my_id, headnode, num_processes);

    // close START, RAWOUT & RAWSORT
    fclose(START);
    fclose(RAWOUT);
    fclose(RAWSORT);

    // print the regeneration summary
    RAWOUT = fopen("regenSummary", "w");
    printRegenSummaryData(regen, level, RAWOUT);
    fclose(RAWOUT);

    // wait until the files are closed
    MPI_Barrier(MPI_COMM_WORLD);

    // check for path crossings - if this is not the last level - do not delete!
    if (level + 1 < num_levels)
    { // make sure that no path crossing happened on this level
      if (regen->level[level].useIntrinsicSlice)
        num_crossings = parallel_midpoint_checking(midFile, midTrackFile, 0, regen->level[level].num_paths, regen->level[level].level + regen->level[level].depth, midpoint_tol, my_id, headnode, num_processes);
      else
        num_crossings = parallel_midpoint_checking(midFile, midTrackFile, 0, regen->level[level].num_paths, regen->num_variables, midpoint_tol, my_id, headnode, num_processes);

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
    { // setup RAWOUT for moving
      size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, level + 1);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawPrepareFile, level + 1);
      RAWOUT = fopen(str, "w"); 

      // setup the name of the file to contain the next set of start points
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, level + 1);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, level + 1);

      // prepare the next set of start points
      num_moved = head_regenPrepareNextLevel(minPacketSize, maxPacketSize, trackCount, level, pathMod, T, regen, START, str, OUT, RAWOUT, FAIL, my_id, headnode, num_processes);

      // close RAWOUT
      fclose(RAWOUT);

      // wait until the files are closed
      MPI_Barrier(MPI_COMM_WORLD);

      // make sure that no path crossing happened on this level - do not delete!
      if (regen->level[level+1].useIntrinsicSlice)
        num_crossings = parallel_midpoint_checking(midFile, midPrepareFile, 0, num_moved, regen->level[level+1].level + regen->level[level+1].depth, midpoint_tol_prepare, my_id, headnode, num_processes);
      else
        num_crossings = parallel_midpoint_checking(midFile, midPrepareFile, 0, num_moved, regen->num_variables, midpoint_tol_prepare, my_id, headnode, num_processes);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

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

  // release the memory on the previous level that is no longer needed
  clearRegenLevelStructures(1, regen, T->MPType, num_levels - 1);

  // set the number of paths tracked
  trackCount->numPoints = trackCount->successes + trackCount->failures;

  // wait until the files are closed
  MPI_Barrier(MPI_COMM_WORLD);

  // combine the midpath_data files - okay to delete them
  combine_midpath_data(midFile, midTrackFile, 1, headnode, num_processes);

  // clear memory
  free(str);

  return;
}

void head_regenTrackLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int curr_level, int pathMod, tracker_config_t *T, regen_t *regen, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this level in parallel     *
\***************************************************************/
{
  int (*change_prec)(void const *, int) = NULL;
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  change_prec = &change_regen_prec;
  create_send_packet = &regen_create_send_packet_track;
  recv_store_packet = &regen_recv_store_packet_track;

  // setup the current level
  regen->curr_level_num = curr_level;

  // initialize the counts
  regen->level[curr_level].num_inf = 0;
  regen->level[curr_level].num_sing = 0;
  regen->level[curr_level].num_nonsing = 0;
  regen->level[curr_level].num_higher_dim = 0;
  regen->level[curr_level].num_bad = 0;

  // display messages
  printf("\nTracking regeneration level %d of %d: %d path%s to track.\n", curr_level, regen->num_levels, regen->level[curr_level].num_paths, regen->level[curr_level].num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration level %d.\n", curr_level);
  fprintf(OUT, "*****************************************************\n");

  // do the actual parallel tracking
  head_trackPaths2("Tracking", 0, regen->level[curr_level].num_paths, minPacketSize, maxPacketSize, trackCount, pathMod, T, regen, regen, change_prec, START, OUT, RAWOUT, FAIL, RAWSORT, NONSOLN, NULL, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

int head_regenPrepareNextLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int curr_level, int pathMod, tracker_config_t *T, regen_t *regen, FILE *START, char *nextStartFile, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths actually tracked               *
* NOTES: prepares the next set of start points in parallel      *
\***************************************************************/
{
  int i, j, k, l, count, var_gp, numNewPts = 0, numOrigPts = regen->level[curr_level].num_nonsing, next_level_num = curr_level + 1;
  int startFunc = regen->level[next_level_num].level, depth = regen->level[next_level_num].depth;
  int pathsTracked = 0, endFunc = startFunc + depth;
  int **new_linear = NULL, **new_linear_degree = NULL;
  mat_d B_transpose_d;
  mat_mp B_transpose_mp;
  FILE *NEXTSTARTPTS = fopen(nextStartFile, "w");
  int good_size = 0, **goodLoc = NULL;
  int (*change_prec)(void const *, int) = &change_regen_prec;

  // initialize
  init_mat_d(B_transpose_d, 0, 0);
  init_mat_mp(B_transpose_mp, 0, 0);

  // setup a new patch & send it to the workers
  if (T->MPType == 0)
  { // create new patch coeff in _d
    make_vec_random_d(regen->main_homVar_d, regen->num_variables);
    bcast_vec_d(regen->main_homVar_d, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // create new patch coeff in _mp
    make_vec_random_mp(regen->main_homVar_mp, regen->num_variables);
    bcast_vec_mp(regen->main_homVar_mp, my_id, headnode);
  } 
  else
  { // create new patch coeff in _d, _mp, _rat 
    make_vec_random_rat(regen->main_homVar_d, regen->main_homVar_mp, regen->main_homVar_rat, regen->num_variables, regen->curr_precision, T->AMP_max_prec, 0, 0);
    bcast_vec_rat(&regen->main_homVar_rat, regen->num_variables, my_id, headnode);
  }

  // setup the structures for the next level
  setupRegenLevelStructures(regen, T->MPType, next_level_num, T->AMP_max_prec);

  // send the structures for the next level to the workers
  bcast_regenLevel_t(&regen->level[next_level_num], T->MPType, regen->curr_precision, my_id, headnode);

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

  // track the paths
  numNewPts = head_regenMoveOrigPts(minPacketSize, maxPacketSize, numOrigPts, pathsTracked, pathMod, regen, T, new_linear, new_linear_degree, count, B_transpose_d, B_transpose_mp, OUT, RAWOUT, FAIL, START, NEXTSTARTPTS, curr_level, trackCount, change_prec, my_id, headnode, num_processes);

  // sotre the number of new paths for next level
  regen->level[next_level_num].num_paths = numNewPts;

  // print the relevant data to NEXTSTARTPTS so that it can be used to rerun this exact same problem
  printRegenRelevantData(regen, T->MPType, next_level_num, NEXTSTARTPTS);

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

int head_regenMoveOrigPts(int minPacketSize, int maxPacketSize, int numOrigPts, int pathsTracked, int pathMod, regen_t *regen, tracker_config_t *T, int **new_linear, int **new_linear_degree, int new_count, mat_d B_transpose_d, mat_mp B_transpose_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *START, FILE *NEXTSTARTPTS, int level_num, trackingStats *trackCount, int (*change_prec)(void const *, int), int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths moved to NEXTSTARTPTS          *
* NOTES: tracks the paths to prepare the next level             *
\***************************************************************/
{
  int i, j, k, num_input_vars, num_output_vars, linearSize, new_linearSize, allGood, numNewPts = 0, num_good = 0, next_level_num = level_num + 1;
  int input_IntrinsicSlice = regen->level[level_num].useIntrinsicSlice;
  int output_IntrinsicSlice = regen->level[next_level_num].useIntrinsicSlice;
  int prec, curr_linear, curr_linear_degree, currPathNum = 0, num_var_gps = regen->num_var_gps;
  int *temp_linear = NULL, *temp_linear_degree = NULL, *rV = NULL;
  FILE *LINEARSTART = fopen("start_temp", "w+"), *TEMPSTARTPTS = fopen("regenTempFile", "w+");
  comp_d tempComp_d;
  comp_mp tempComp_mp;
  point_d tempPoint_d;
  point_mp tempPoint_mp;

  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  create_send_packet = &regen_create_send_packet_move;
  recv_store_packet = &regen_recv_store_packet_move;

  init_mp(tempComp_mp);
  init_point_d(tempPoint_d, 0);
  init_point_mp(tempPoint_mp, 0);

  linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  if (input_IntrinsicSlice)
    num_input_vars = linearSize;
  else
    num_input_vars = regen->num_variables;

  new_linearSize = regen->level[next_level_num].level + regen->level[next_level_num].depth;
  if (regen->level[next_level_num].useIntrinsicSlice)
    num_output_vars = new_linearSize;
  else
    num_output_vars = regen->num_variables;

  temp_linear = (int *)bmalloc(new_linearSize * sizeof(int));
  temp_linear_degree = (int *)bmalloc(new_linearSize * sizeof(int));
  for (i = 0; i < new_linearSize; i++)
    temp_linear[i] = temp_linear_degree[i] = 0;

  numNewPts = numOrigPts * new_count;
  rV = (int *)bmalloc(numNewPts * sizeof(int));
  for (i = 0; i < numNewPts; i++)
    rV[i] = 0;

  // loop over the points, printing the data to LINEARSTART and to TEMPSTARTPTS (if needed)
  for (i = 0; i < numOrigPts; i++)
  { // read in the point for file
    fscanf(START, "%d\n%d\n", &j, &prec);
    if (prec < 64)
    { // _d
      increase_size_point_d(tempPoint_d, num_input_vars);
      tempPoint_d->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
        fscanf(START, "%lf %lf\n", &tempPoint_d->coord[j].r, &tempPoint_d->coord[j].i);

      // convert to extrinsic coordinates
      if (input_IntrinsicSlice)
        intrinsicToExtrinsic_d(tempPoint_d, tempPoint_d, regen->level[level_num].B_d, regen->level[level_num].p_d);

      // move it to the new patch
      move_to_patch_vec_d(tempPoint_d, tempPoint_d, regen->main_homVar_d);

      // convert to coordinates for next level
      if (output_IntrinsicSlice)
        extrinsicToIntrinsic_d(tempPoint_d, tempPoint_d, B_transpose_d, regen->level[next_level_num].p_d);
    }
    else
    { // _mp
      change_prec(regen, prec);
      change_prec_point_mp(tempPoint_mp, prec);
      increase_size_point_mp(tempPoint_mp, num_input_vars);
      tempPoint_mp->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(tempPoint_mp->coord[j].r, START, 10);
        mpf_inp_str(tempPoint_mp->coord[j].i, START, 10);
      }

      // convert to extrinsic coordinates
      if (input_IntrinsicSlice)
        intrinsicToExtrinsic_mp(tempPoint_mp, tempPoint_mp, regen->level[level_num].B_mp, regen->level[level_num].p_mp);

      // move it to the new patch
      move_to_patch_vec_mp(tempPoint_mp, tempPoint_mp, regen->main_homVar_mp);

      // convert to coordinates for next level
      if (output_IntrinsicSlice)
      { // setup B_transpose_mp properly
        change_prec_mat_mp(B_transpose_mp, prec);
        transpose_mp(B_transpose_mp, regen->level[next_level_num].B_mp);
        extrinsicToIntrinsic_mp(tempPoint_mp, tempPoint_mp, B_transpose_mp, regen->level[next_level_num].p_mp);
      }
    } 
    // move past linear information
    for (j = 0; j < linearSize; j++)
      fscanf(START, "%d %d\n", &curr_linear, &curr_linear_degree);

    // print the points to either LINEARSTART or TEMPSTARTPTS
    for (j = 0; j < new_count; j++)
    { // setup linear information
      for (k = linearSize; k < new_linearSize; k++)
      {
        temp_linear[k] = new_linear[j][k - linearSize];
        temp_linear_degree[k] = new_linear_degree[j][k - linearSize];
      }

      // print data
      if (j == 0 && num_var_gps == 1)
      { // print this one to TEMPSTARTPTS
        fprintf(TEMPSTARTPTS, "%d\n", currPathNum);
        printPointLinearDegree(TEMPSTARTPTS, tempPoint_d, tempPoint_mp, prec, temp_linear, temp_linear_degree, new_linearSize);      
      }
      else
      { // print this one to LINEARSTART
        fprintf(LINEARSTART, "%d\n", currPathNum);
        printPointLinearDegree(LINEARSTART, tempPoint_d, tempPoint_mp, prec, temp_linear, temp_linear_degree, new_linearSize); 
      }
      currPathNum++;
    }
  } 

  // rewind so that we can use this to track the points
  rewind(LINEARSTART);

  // update to the next level
  regen->curr_level_num = next_level_num;

  // print header for level
  printf("\nPreparing regeneration level %d of %d: %d witness point%s to move.\n", next_level_num, regen->num_levels, pathsTracked, pathsTracked == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Moving linears to regeneration level %d.\n", next_level_num);
  fprintf(OUT, "*****************************************************\n");

  // do the actual parallel tracking
  head_trackPaths2("Preparing", 0, pathsTracked, minPacketSize, maxPacketSize, trackCount, pathMod, T, regen, regen, change_prec, LINEARSTART, OUT, RAWOUT, FAIL, TEMPSTARTPTS, NULL, rV, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  // close & delete LINEARSTART
  fclose(LINEARSTART);
  remove("start_temp");

  // rewind TEMPSTARTPTS
  rewind(TEMPSTARTPTS);

  // loop over the points setting up rV properly
  if (num_var_gps == 1)
  { // all should be good!
    currPathNum = 0;
    for (i = 0; i < numOrigPts; i++)
    {
      allGood = 1;
      for (j = 0; j < new_count; j++)
      {
        if (rV[j + currPathNum] != 0)
        { // store this path to bad_paths
          printf("WARNING: When preparing path %d in level %d, it had retVal %d.\n", i, next_level_num, rV[j + currPathNum]);
          allGood = 0;
        }
      }

      if (!allGood)
      { // set each as bad
        for (j = 0; j < new_count; j++)
          rV[j + currPathNum] = 1;
      }

      // increment currPathNum
      currPathNum += new_count;
    }
  }

  // count the number of good points
  num_good = 0;
  for (i = 0; i < numNewPts; i++)
    if (rV[i] == 0)
      num_good++;

  fprintf(NEXTSTARTPTS, "%d\n\n", num_good);

  // loop over the points printing to NEXTSTARTPOINTS if needed
  for (i = 0; i < numNewPts; i++)
  { // read in the path number
    fscanf(TEMPSTARTPTS, "%d\n", &j);
    // see if we should move it
    if (rV[j] == 0)
    { // move it to NEXTSTARTPTS
      if (T->MPType == 0 || T->MPType == 2)
      { // move double precision
        for (j = 0; j < num_output_vars; j++)
        { // read in jth coordinate
          fscanf(TEMPSTARTPTS, "%lf %lf;\n", &tempComp_d->r, &tempComp_d->i);
          // print coordinate j
          print_d(NEXTSTARTPTS, 0, tempComp_d);
          fprintf(NEXTSTARTPTS, ";\n");
        }
      }
      else
      { // move mulit precision
        for (j = 0; j < num_output_vars; j++)
        { // read in jth coordinate
          mpf_inp_str(tempComp_mp->r, TEMPSTARTPTS, 10);
          mpf_inp_str(tempComp_mp->i, TEMPSTARTPTS, 10);
          // print coordinate j
          print_mp(NEXTSTARTPTS, 0, tempComp_mp);
          fprintf(NEXTSTARTPTS, ";\n");
        }
      }
      for (j = 0; j < new_linearSize; j++)
      { // read in linear & degree for function j
        fscanf(TEMPSTARTPTS, "%d %d\n", &curr_linear, &curr_linear_degree);
        // print linear & degree for function j
        fprintf(NEXTSTARTPTS, "%d %d\n", curr_linear, curr_linear_degree);
      } 
      fprintf(NEXTSTARTPTS, "\n");
    }
    else
    { // move past this
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

  fclose(TEMPSTARTPTS);
  remove("regenTempFile");

  clear_mp(tempComp_mp);
  clear_point_d(tempPoint_d);
  clear_point_mp(tempPoint_mp);

  free(temp_linear);
  free(temp_linear_degree);
  free(rV);

  return num_good;
}

///////////////// WORKER FUNCTIONS ///////////////////

void worker_regen_zero_dim(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does regen tracking for zero dimensional - 'worker'    *
\***************************************************************/
{
  trackingStats trackCount;
  tracker_config_t T;
  regen_t regen;
  basic_eval_data_d BED_d;
  basic_eval_data_mp BED_mp;

  // setup T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision - does not hurt if doing only double precision
  initMP(T.Precision);

  if (T.MPType == 0)
  { // send ED_d - only double
    bcast_basic_eval_data_d(&BED_d, T.MPType, my_id, headnode);
  }
  else if (T.MPType == 1)
  { // send ED_mp - only fixed MP
    bcast_basic_eval_data_mp(&BED_mp, T.MPType, my_id, headnode);
  }
  else
  { // send ED_d to the workers - using AMP
    bcast_basic_eval_data_amp(&BED_d, my_id, headnode);
  }

  // recv regen
  bcast_regen_t(&regen, T.MPType, my_id, headnode);

  // setup square 
  if (T.MPType == 0)
  { // setup square_d
    regen.square_d = &BED_d.squareSystem;
    regen.square_mp = NULL;
  }
  else if (T.MPType == 1)
  { // setup square_mp
    regen.square_mp = &BED_mp.squareSystem;
    regen.square_d = NULL;
  }
  else
  { // setup square_d & square_mp
    regen.square_d = &BED_d.squareSystem;
    regen.square_mp = &BED_d.BED_mp->squareSystem;
  }

  // setup the evaluators inside regen
  regen.square_eval_d = square_system_eval_d;
  regen.square_eval_mp = square_system_eval_mp;
  regen.change_square_prec = change_square_prec;

  // now that everything is sent, do the main regeneration tracking
  worker_regen_track(&regen, &T, &trackCount, my_id, num_processes, headnode);

  // clear memory
  clearRegenRandom_zero_dim(1, &regen, T.MPType);
  if (T.MPType == 0 || T.MPType == 2)
  {
    basic_eval_clear_d(&BED_d, 0, T.MPType);
  }
  else
  {
    basic_eval_clear_mp(&BED_mp, 0, 1);
  }
  tracker_config_clear(&T);
  clearMP();

  return;
}

void worker_regen_track(regen_t *regen, tracker_config_t *T, trackingStats *trackCount, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does standard regen tracking - 'worker'                *
* Assume that everything is setup/sent so that we can track     *
\***************************************************************/
{
  int startLevel, level;
  size_t size;
  char *str = NULL, midTrackFile[] = "midout_track", midPrepareFile[] = "midout_prepare";
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // allocate levels
  regen->level = (regenLevel_t *)bmalloc(regen->num_levels * sizeof(regenLevel_t));

  // recv 'startLevel'
  MPI_Bcast(&startLevel, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // setup the starting level
  bcast_regenLevel_t(&regen->level[startLevel], T->MPType, regen->curr_precision, my_id, headnode);

  // setup the local files of OUT & FAIL & RAWOUT
  size = 1 + snprintf(NULL, 0, "output_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "output_%d", my_id);
  OUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  FAIL = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  RAWOUT = fopen(str, "w");

  // main tracking loop
  for (level = startLevel; level < regen->num_levels; level++)
  { // setup MIDOUT for this level
    size = 1 + snprintf(NULL, 0, "%s_%d", midTrackFile, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", midTrackFile, my_id);
    MIDOUT = fopen(str, "w");

    // track this level
    worker_regenTrackLevel(trackCount, level, T, regen, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    // close MIDOUT
    fclose(MIDOUT);

    // wait until everybody has closed MIDOUT
    MPI_Barrier(MPI_COMM_WORLD);

    // consider doing midpoint checking in parallel

    // see if we need to prepare the next level
    if (level + 1 < regen->num_levels)
    { // reopen MIDOUT for preparing next level
      size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", midPrepareFile, my_id);
      MIDOUT = fopen(str, "w");

      // prepare the next set of start 
      worker_regenPrepareLevel(trackCount, level, T, regen, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

      // close MIDOUT
      fclose(MIDOUT);

      // wait until everybody has closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel

      // release the memory on the previous level that is no longer needed
      clearRegenLevelStructures(1, regen, T->MPType, level);
    }
  }

  // close the files
  fclose(OUT);
  fclose(FAIL);
  fclose(RAWOUT);

  // delete RAWOUT
  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  remove(str);

  // delete MIDOUT for preparing (tracking is deleted by head)
  size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "%s_%d", midPrepareFile, my_id);
  remove(str);

  // delete FAIL
  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  remove(str);

  // wait until everybody has closed their files
  MPI_Barrier(MPI_COMM_WORLD);

  // clear memory
  free(str);

  return;
}

void worker_regenTrackLevel(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_t *regen, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this level in parallel     *
\***************************************************************/
{
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*useSharpener)(int, int, void const *, void const *);
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

  // setup the evaluators
  ptr_to_eval_d = &regen_eval_d;
  ptr_to_eval_mp = &regen_eval_mp;
  change_prec = &change_regen_prec;
  useSharpener = &regen_useSharpener;
  recv_track_send_packet = &regen_recv_track_send_packet;

  // setup the current codimension
  regen->curr_level_num = curr_level;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking level %d.\n", curr_level);
  fprintf(OUT, "*****************************************************\n");

  // do the actual tracking
  worker_trackPaths2(trackCount, T, regen, regen, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  return;
}

void worker_regenPrepareLevel(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_t *regen, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prepares the start point for this level in parallel    *
\***************************************************************/
{
  int i, next_level_num = curr_level + 1;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*useSharpener)(int, int, void const *, void const *);
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

  // store the old tolerances
  double basicTol = T->basicNewtonTol, endgameTol = T->endgameNewtonTol, finalTol = T->final_tolerance;
  // setup the new tolerances
  T->basicNewtonTol = T->sliceBasicNewtonTol;
  T->endgameNewtonTol = T->sliceEndgameNewtonTol;
  T->final_tolerance = T->sliceFinalTol;

  // recv the new patch
  if (T->MPType == 0)
  { // recv in _d
    bcast_vec_d(regen->main_homVar_d, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // recv in _mp
    bcast_vec_mp(regen->main_homVar_mp, my_id, headnode);
  }
  else
  { // recv in _rat and copy to _d & _mp 
    bcast_vec_rat(&regen->main_homVar_rat, regen->num_variables, my_id, headnode);
    for (i = 0; i < regen->num_variables; i++)
    { 
      rat_to_d(&regen->main_homVar_d->coord[i], regen->main_homVar_rat[i]);
      rat_to_mp(&regen->main_homVar_mp->coord[i], regen->main_homVar_rat[i]);
    }
  }

  // recv the structures for the next level
  bcast_regenLevel_t(&regen->level[next_level_num], T->MPType, regen->curr_precision, my_id, headnode);

  // setup the evaluators
  ptr_to_eval_d = &regen_moving_linear_eval_d;
  ptr_to_eval_mp = &regen_moving_linear_eval_mp;
  change_prec = &change_regen_prec;
  useSharpener = &regen_useSharpener_move;
  recv_track_send_packet = &regen_recv_linear_track_send_packet;

  // setup the current codimension
  regen->curr_level_num = next_level_num;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Moving linears to regeneration level %d.\n", next_level_num);
  fprintf(OUT, "*****************************************************\n");

  // do the actual tracking
  worker_trackPaths2(trackCount, T, regen, regen, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  // set the tolerances back
  T->basicNewtonTol = basicTol;
  T->endgameNewtonTol = endgameTol;
  T->final_tolerance = finalTol;

  return;
}

////////////////////// OTHER FUNCTIONS //////////////////////////

void create_regen_int(MPI_Datatype *mpi_regen_int)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the MPI datatype mpi_regen_int                 *
\***************************************************************/
{
  MPI_Datatype mpi_preproc_int;

  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  regen_t_int tempRegen;

  // create datatypes needed
  create_preproc_data_int(&mpi_preproc_int);

  // arrays for length, displacement and datatypes in mpi_regen_int
  int regen_int_length[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype regen_int_datatypes[13] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, mpi_preproc_int, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint regen_int_displacements[13];

  // calculate displacements
  regen_int_displacements[0] = 0;
  MPI_Address(&tempRegen.num_levels, &startaddress);
  MPI_Address(&tempRegen.num_funcs, &address);
  regen_int_displacements[1] = address - startaddress;
  MPI_Address(&tempRegen.num_variables, &address);
  regen_int_displacements[2] = address - startaddress;
  MPI_Address(&tempRegen.num_var_gps, &address);
  regen_int_displacements[3] = address - startaddress;
  MPI_Address(&tempRegen.PPD_int, &address);
  regen_int_displacements[4] = address - startaddress;
  MPI_Address(&tempRegen.patch_rows, &address);
  regen_int_displacements[5] = address - startaddress;
  MPI_Address(&tempRegen.patch_cols, &address);
  regen_int_displacements[6] = address - startaddress;
  MPI_Address(&tempRegen.curr_precision, &address);
  regen_int_displacements[7] = address - startaddress;
  MPI_Address(&tempRegen.num_comp_d, &address);
  regen_int_displacements[8] = address - startaddress;
  MPI_Address(&tempRegen.totalLength, &address);
  regen_int_displacements[9] = address - startaddress;
  MPI_Address(&tempRegen.noChanges, &address);
  regen_int_displacements[10] = address - startaddress;
  MPI_Address(&tempRegen.numSubFuncs, &address);
  regen_int_displacements[11] = address - startaddress;
  MPI_Address(&tempRegen.numInts, &address);
  regen_int_displacements[12] = address - startaddress;

  // build the mpi datatype mpi_regen_int
  MPI_Type_struct(13, regen_int_length, regen_int_displacements, regen_int_datatypes, mpi_regen_int);
  MPI_Type_commit(mpi_regen_int);

  return;
}

void create_regenLevel_int(MPI_Datatype *mpi_regenLevel_int)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the MPI datatype mpi_regenLevel_int            *
\***************************************************************/
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  regenLevel_t_int tempRegen;

  // arrays for length, displacement and datatypes in mpi_regenLevel_int
  int regen_int_length[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype regen_int_datatypes[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint regen_int_displacements[8];

  // calculate displacements
  regen_int_displacements[0] = 0;
  MPI_Address(&tempRegen.level, &startaddress);
  MPI_Address(&tempRegen.depth, &address);
  regen_int_displacements[1] = address - startaddress;
  MPI_Address(&tempRegen.useIntrinsicSlice, &address);
  regen_int_displacements[2] = address - startaddress;
  MPI_Address(&tempRegen.B_rows, &address);
  regen_int_displacements[3] = address - startaddress;
  MPI_Address(&tempRegen.B_cols, &address);
  regen_int_displacements[4] = address - startaddress;
  MPI_Address(&tempRegen.p_size, &address);
  regen_int_displacements[5] = address - startaddress;
  MPI_Address(&tempRegen.num_comp_d, &address);
  regen_int_displacements[6] = address - startaddress;
  MPI_Address(&tempRegen.totalLength, &address);
  regen_int_displacements[7] = address - startaddress;

  // build the mpi datatype mpi_regenLevel_int
  MPI_Type_struct(8, regen_int_length, regen_int_displacements, regen_int_datatypes, mpi_regenLevel_int);
  MPI_Type_commit(mpi_regenLevel_int);

  return;
}

void cp_regen_int(void *regen_out, void *regen_in, int MPType, int **Degrees, int **PPD_type, int **PPD_size, char **regenStr, comp_d **coeff_d, int freeStr, int **instCount, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - regen_out is _int, regen_in is _t      *
*                otherwise - regen_out is _t, regen_in is _int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of regen_in to regen_out                 *
\***************************************************************/
{
  int i, j, k, l, count, base = 10;

  if (inType == 0)
  { // _t to _int so that it can be sent
    regen_t *regen_right = (regen_t *)regen_in;
    regen_t_int *regen_left = (regen_t_int *)regen_out;

    // setup basic info
    regen_left->num_levels = regen_right->num_levels;
    regen_left->num_funcs = regen_right->num_funcs;
    regen_left->num_variables = regen_right->num_variables;
    regen_left->num_var_gps = regen_right->num_var_gps;
    regen_left->curr_precision = regen_right->curr_precision;

    // setup PPD_int
    cp_preproc_data_int(&regen_left->PPD_int, &regen_right->PPD, PPD_type, PPD_size, inType);

    // setup Degrees
    *Degrees = (int *)bmalloc(regen_right->num_funcs * regen_right->num_var_gps * sizeof(int));
    count = 0;
    for (i = 0; i < regen_right->num_funcs; i++)
      for (j = 0; j < regen_right->num_var_gps; j++)
      {
        (*Degrees)[count] = regen_right->degrees[i][j];
        count++;
      }

    // setup instCount, if needed
    regen_left->noChanges = regen_right->noChanges;
    regen_left->numSubFuncs = regen_right->numSubFuncs;
    if (regen_left->noChanges)
    { // count the number of things to send
      regen_left->numInts = 4 * (regen_right->numSubFuncs + regen_right->num_funcs) + regen_right->numSubFuncs * regen_right->num_funcs;
      *instCount = (int *)bmalloc(regen_left->numInts * sizeof(int));

      // setup instCount
      for (i = 0; i < regen_right->numSubFuncs; i++)
      { // copy the things related to subfunctions
        (*instCount)[i] = regen_right->startSub[i];
        (*instCount)[i + regen_right->numSubFuncs] = regen_right->endSub[i];
        (*instCount)[i + 2 * regen_right->numSubFuncs] = regen_right->startJvsub[i];
        (*instCount)[i + 3 * regen_right->numSubFuncs] = regen_right->endJvsub[i];
      }
      for (i = 4 * regen_right->numSubFuncs; i < 4 * regen_right->numSubFuncs + regen_right->num_funcs; i++)
      { // copy the things releted to functions
        (*instCount)[i] = regen_right->startFunc[i - 4 * regen_right->numSubFuncs];
        (*instCount)[i + regen_right->num_funcs] = regen_right->endFunc[i - 4 * regen_right->numSubFuncs];
        (*instCount)[i + 2 * regen_right->num_funcs] = regen_right->startJv[i - 4 * regen_right->numSubFuncs];
        (*instCount)[i + 3 * regen_right->num_funcs] = regen_right->endJv[i - 4 * regen_right->numSubFuncs];
      }
      for (i = 0; i < regen_right->num_funcs; i++)
        for (j = 0; j < regen_right->numSubFuncs; j++)
          (*instCount)[4 * (regen_right->num_funcs + regen_right->numSubFuncs) + i * regen_right->numSubFuncs + j] = regen_right->subFuncsBelow[i][j];
    }
    else
    { // set to NULL
      regen_left->numInts = 0;
      *instCount = NULL;
    }

    // copy gamma, patchCoeff, coeff, mainCoeff & main_homVar
    if (MPType == 0)
    { // use _d
      regen_left->patch_rows = regen_right->patchCoeff_d->rows;
      regen_left->patch_cols = regen_right->patchCoeff_d->cols;

      // count the total number of comp_d that are to be sent
      regen_left->num_comp_d = 1 + regen_left->patch_rows * regen_left->patch_cols + regen_left->num_funcs * regen_left->num_variables + regen_left->num_variables;
      for (i = 0; i < regen_right->num_funcs; i++)
        for (j = 0; j < regen_right->num_var_gps; j++)
          regen_left->num_comp_d += regen_right->degrees[i][j] * regen_right->num_variables;

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(regen_left->num_comp_d * sizeof(comp_d));
      count = 0;
      // copy gamma
      set_d((*coeff_d)[count], regen_right->gamma_d);
      count++;
      // copy patchCoeff
      for (i = 0; i < regen_left->patch_rows; i++)
        for (j = 0; j < regen_left->patch_cols; j++)
        {
          set_d((*coeff_d)[count], &regen_right->patchCoeff_d->entry[i][j]);
          count++;
        }
      // copy coeff
      for (i = 0; i < regen_right->num_funcs; i++)
        for (j = 0; j < regen_right->num_var_gps; j++)
          for (k = 0; k < regen_right->degrees[i][j]; k++)
            for (l = 0; l < regen_right->num_variables; l++)
            {
              set_d((*coeff_d)[count], regen_right->coeff_d[i][j][k][l]);
              count++;
            }
      // copy mainCoeff
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        {
          set_d((*coeff_d)[count], &regen_right->mainCoeff_d->entry[i][j]);
          count++;
        }
      // copy main_homVar
      for (i = 0; i < regen_left->num_variables; i++)
      {
        set_d((*coeff_d)[count], &regen_right->main_homVar_d->coord[i]);
        count++;
      }

      // setup totalLength
      regen_left->totalLength = 0;     
    }
    else if (MPType == 1)
    { // use _mp
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      regen_left->patch_rows = regen_right->patchCoeff_mp->rows;
      regen_left->patch_cols = regen_right->patchCoeff_mp->cols;
      regen_left->num_comp_d = 0; // no comp_d setup

      // setup regenStr
      *regenStr = NULL;
      count = 0;
      // copy gamma_mp to regenStr
      strR = mpf_to_str(regen_right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      strI = mpf_to_str(regen_right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
      strcpy(&(*regenStr)[count], strR);
      strcpy(&(*regenStr)[count + sizeR], strI);
      count += sizeR + sizeI;
      free(strR);
      free(strI);
      // copy patchCoeff_mp to regenStr
      for (i = 0; i < regen_left->patch_rows; i++)
        for (j = 0; j < regen_left->patch_cols; j++)
        { // real part
          strR = mpf_to_str(regen_right->patchCoeff_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(regen_right->patchCoeff_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update regenStr
          *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*regenStr)[count], strR);
          strcpy(&(*regenStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy coeff_mp to regenStr
      for (i = 0; i < regen_right->num_funcs; i++)
        for (j = 0; j < regen_right->num_var_gps; j++)
          for (k = 0; k < regen_right->degrees[i][j]; k++)
            for (l = 0; l < regen_right->num_variables; l++)
            { // real part
              strR = mpf_to_str(regen_right->coeff_mp[i][j][k][l]->r, base);
              sizeR = strlen(strR) + 1; // +1 for '\0'
              // imag part
              strI = mpf_to_str(regen_right->coeff_mp[i][j][k][l]->i, base);
              sizeI = strlen(strI) + 1; // +1 for '\0'

              // update regenStr
              *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
              strcpy(&(*regenStr)[count], strR);
              strcpy(&(*regenStr)[count + sizeR], strI);
              // update count
              count += sizeR + sizeI;

              // free strR, strI
              free(strR);
              free(strI);
            }
      // copy mainCoeff_mp to regenStr
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        { // real part
          strR = mpf_to_str(regen_right->mainCoeff_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(regen_right->mainCoeff_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update regenStr
          *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*regenStr)[count], strR);
          strcpy(&(*regenStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy main_homVar_mp to regenStr
      for (i = 0; i < regen_left->num_variables; i++)
      { // real part
        strR = mpf_to_str(regen_right->main_homVar_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(regen_right->main_homVar_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update regenStr
        *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
        strcpy(&(*regenStr)[count], strR);
        strcpy(&(*regenStr)[count + sizeR], strI);
        // update count
        count += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      regen_left->totalLength = count;
      regen_left->num_comp_d = 0;
    }
    else
    { // use _rat
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      regen_left->patch_rows = regen_right->patchCoeff_d->rows;
      regen_left->patch_cols = regen_right->patchCoeff_d->cols;
      regen_left->num_comp_d = 0; // no comp_d setup

      // setup regenStr
      *regenStr = NULL;
      count = 0;
      // copy gamma_rat to regenStr
      strR = mpq_get_str(NULL, base, regen_right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      strI = mpq_get_str(NULL, base, regen_right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
      strcpy(&(*regenStr)[count], strR);
      strcpy(&(*regenStr)[count + sizeR], strI);
      count += sizeR + sizeI;
      free(strR);
      free(strI);
      // copy patchCoeff_rat to regenStr
      for (i = 0; i < regen_left->patch_rows; i++)
        for (j = 0; j < regen_left->patch_cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, regen_right->patchCoeff_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, regen_right->patchCoeff_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update regenStr
          *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*regenStr)[count], strR);
          strcpy(&(*regenStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy coeff_rat to regenStr
      for (i = 0; i < regen_right->num_funcs; i++)
        for (j = 0; j < regen_right->num_var_gps; j++)
          for (k = 0; k < regen_right->degrees[i][j]; k++)
            for (l = 0; l < regen_right->num_variables; l++)
            { // real part
              strR = mpq_get_str(NULL, base, regen_right->coeff_rat[i][j][k][l][0]);
              sizeR = strlen(strR) + 1; // +1 for '\0'
              // imag part
              strI = mpq_get_str(NULL, base, regen_right->coeff_rat[i][j][k][l][1]);
              sizeI = strlen(strI) + 1; // +1 for '\0'

              // update regenStr
              *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
              strcpy(&(*regenStr)[count], strR);
              strcpy(&(*regenStr)[count + sizeR], strI);
              // update count
              count += sizeR + sizeI;

              // free strR, strI
              free(strR);
              free(strI);
            }
      // copy mainCoeff_rat to regenStr
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        { // real part
          strR = mpq_get_str(NULL, base, regen_right->mainCoeff_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, regen_right->mainCoeff_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update regenStr
          *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*regenStr)[count], strR);
          strcpy(&(*regenStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy main_homVar_rat to regenStr
      for (i = 0; i < regen_left->num_variables; i++)
      { // real part
        strR = mpq_get_str(NULL, base, regen_right->main_homVar_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, regen_right->main_homVar_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update regenStr
        *regenStr = (char *)brealloc(*regenStr, (count + sizeR + sizeI) * sizeof(char));
        strcpy(&(*regenStr)[count], strR);
        strcpy(&(*regenStr)[count + sizeR], strI);
        // update count
        count += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      regen_left->totalLength = count;
      regen_left->num_comp_d = 0;
    }
  }
  else
  { // _int to _t so that is can be used
    regen_t_int *regen_right = (regen_t_int *)regen_in;
    regen_t *regen_left = (regen_t *)regen_out;

    // setup basic info
    regen_left->num_levels = regen_right->num_levels;
    regen_left->num_funcs = regen_right->num_funcs;
    regen_left->num_variables = regen_right->num_variables;
    regen_left->num_var_gps = regen_right->num_var_gps;
    regen_left->curr_precision = regen_right->curr_precision;

    // setup PPD
    cp_preproc_data_int(&regen_left->PPD, &regen_right->PPD_int, PPD_type, PPD_size, inType);

    // setup degrees
    regen_left->degrees = (int **)bmalloc(regen_right->num_funcs * sizeof(int *));
    count = 0;
    for (i = 0; i < regen_right->num_funcs; i++)
    {
      regen_left->degrees[i] = (int *)bmalloc(regen_right->num_var_gps * sizeof(int));
      for (j = 0; j < regen_right->num_var_gps; j++)
      {
        regen_left->degrees[i][j] = (*Degrees)[count];
        count++;
      }
    }
    // clear degrees
    free(*Degrees);

    // setup instCount
    regen_left->noChanges = regen_right->noChanges;
    regen_left->numSubFuncs = regen_right->numSubFuncs;
    if (regen_left->noChanges)
    { // allocate
      regen_left->startSub = (int *)bmalloc(regen_right->numSubFuncs * sizeof(int));
      regen_left->endSub = (int *)bmalloc(regen_right->numSubFuncs * sizeof(int));
      regen_left->startFunc = (int *)bmalloc(regen_right->num_funcs * sizeof(int));
      regen_left->endFunc = (int *)bmalloc(regen_right->num_funcs * sizeof(int));
      regen_left->startJvsub = (int *)bmalloc(regen_right->numSubFuncs * sizeof(int));
      regen_left->endJvsub = (int *)bmalloc(regen_right->numSubFuncs * sizeof(int));
      regen_left->startJv = (int *)bmalloc(regen_right->num_funcs * sizeof(int));
      regen_left->endJv = (int *)bmalloc(regen_right->num_funcs * sizeof(int));

      // setup
      for (i = 0; i < regen_right->numSubFuncs; i++)
      { // copy the things related to subfunctions
        regen_left->startSub[i] = (*instCount)[i];
        regen_left->endSub[i] = (*instCount)[i + regen_right->numSubFuncs];
        regen_left->startJvsub[i] = (*instCount)[i + 2 * regen_right->numSubFuncs];
        regen_left->endJvsub[i] = (*instCount)[i + 3 * regen_right->numSubFuncs];
      }
      for (i = 0; i < regen_right->num_funcs; i++)
      { // copy the things releted to functions
        regen_left->startFunc[i] = (*instCount)[i + 4 * regen_right->numSubFuncs];
        regen_left->endFunc[i] = (*instCount)[i + regen_right->num_funcs + 4 * regen_right->numSubFuncs];
        regen_left->startJv[i] = (*instCount)[i + 2 * regen_right->num_funcs + 4 * regen_right->numSubFuncs];
        regen_left->endJv[i] = (*instCount)[i + 3 * regen_right->num_funcs + 4 * regen_right->numSubFuncs];
      }
      if (regen_right->numSubFuncs > 0)
      { // setup subFuncsBelow
        regen_left->subFuncsBelow = (int **)bmalloc(regen_right->num_funcs * sizeof(int *));
        for (i = 0; i < regen_right->num_funcs; i++)
          regen_left->subFuncsBelow[i] = (int *)bmalloc(regen_right->numSubFuncs * sizeof(int));

        for (i = 0; i < regen_right->num_funcs; i++)
          for (j = 0; j < regen_right->numSubFuncs; j++)
            regen_left->subFuncsBelow[i][j] = (*instCount)[4 * (regen_right->num_funcs + regen_right->numSubFuncs) + i * regen_right->numSubFuncs + j];
      }
      else
        regen_left->subFuncsBelow = NULL;

      // clear instCount
      free(*instCount);
    }
    else
    { // set to NULL
      regen_left->startSub = regen_left->endSub = regen_left->startFunc = regen_left->endFunc = regen_left->startJvsub = regen_left->endJvsub = regen_left->startJv = regen_left->endJv = NULL;
      regen_left->subFuncsBelow = NULL;
    }

    // setup gamma, patchCoeff, coeff, mainCoeff & main_homVar
    if (MPType == 0)
    { // setup _d
      init_mat_d(regen_left->patchCoeff_d, regen_right->patch_rows, regen_right->patch_cols);
      init_mat_d(regen_left->mainCoeff_d, regen_left->num_funcs, regen_left->num_variables);
      init_vec_d(regen_left->main_homVar_d, regen_left->num_variables);
      regen_left->patchCoeff_d->rows = regen_right->patch_rows;
      regen_left->patchCoeff_d->cols = regen_right->patch_cols;
      regen_left->mainCoeff_d->rows = regen_left->num_funcs;
      regen_left->mainCoeff_d->cols = regen_left->num_variables;
      regen_left->main_homVar_d->size = regen_left->num_variables;

      // setup gamma_d
      count = 0;
      set_d(regen_left->gamma_d, (*coeff_d)[count]);
      count++;
      // setup patchCoeff_d
      for (i = 0; i < regen_right->patch_rows; i++)
        for (j = 0; j < regen_right->patch_cols; j++)
        {
          set_d(&regen_left->patchCoeff_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup coeff_d
      regen_left->coeff_d = (comp_d ****)bmalloc(regen_left->num_funcs * sizeof(comp_d ***));
      for (i = 0; i < regen_left->num_funcs; i++)
      {
        regen_left->coeff_d[i] = (comp_d ***)bmalloc(regen_left->num_var_gps * sizeof(comp_d **));
        for (j = 0; j < regen_left->num_var_gps; j++)
        {
          regen_left->coeff_d[i][j] = (comp_d **)bmalloc(regen_left->degrees[i][j] * sizeof(comp_d *));
          for (k = 0; k < regen_left->degrees[i][j]; k++)
          {
            regen_left->coeff_d[i][j][k] = (comp_d *)bmalloc(regen_left->num_variables * sizeof(comp_d));
            for (l = 0; l < regen_left->num_variables; l++)
            {
              set_d(regen_left->coeff_d[i][j][k][l], (*coeff_d)[count]);
              count++;
            }
          }
        }
      }
      // setup mainCoeff
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        {
          set_d(&regen_left->mainCoeff_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup main_homVar
      for (i = 0; i < regen_left->num_variables; i++)
      {
        set_d(&regen_left->main_homVar_d->coord[i], (*coeff_d)[count]);
        count++;
      }

      // free coeff_d
      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup _mp
      init_mp(regen_left->gamma_mp);
      init_mat_mp(regen_left->patchCoeff_mp, regen_right->patch_rows, regen_right->patch_cols);
      init_mat_mp(regen_left->mainCoeff_mp, regen_left->num_funcs, regen_left->num_variables);
      init_vec_mp(regen_left->main_homVar_mp, regen_left->num_variables);
      regen_left->patchCoeff_mp->rows = regen_right->patch_rows;
      regen_left->patchCoeff_mp->cols = regen_right->patch_cols;
      regen_left->mainCoeff_mp->rows = regen_left->num_funcs;
      regen_left->mainCoeff_mp->cols = regen_left->num_variables;
      regen_left->main_homVar_mp->size = regen_left->num_variables;

      // setup gamma_mp
      count = 0;
      mpf_set_str(regen_left->gamma_mp->r, &(*regenStr)[count], base);
      count += 1 + strlen(&(*regenStr)[count]);
      mpf_set_str(regen_left->gamma_mp->i, &(*regenStr)[count], base);
      count += 1 + strlen(&(*regenStr)[count]);
      // setup patchCoeff_mp
      for (i = 0; i < regen_right->patch_rows; i++)
        for (j = 0; j < regen_right->patch_cols; j++)
        { // setup real part
          mpf_set_str(regen_left->patchCoeff_mp->entry[i][j].r, &(*regenStr)[count], base);
          count += 1 + strlen(&(*regenStr)[count]);
          // setup imag part
          mpf_set_str(regen_left->patchCoeff_mp->entry[i][j].i, &(*regenStr)[count], base);
          count += 1 + strlen(&(*regenStr)[count]);
        }
      // setup coeff_mp
      regen_left->coeff_mp = (comp_mp ****)bmalloc(regen_left->num_funcs * sizeof(comp_mp ***));
      for (i = 0; i < regen_left->num_funcs; i++)
      {
        regen_left->coeff_mp[i] = (comp_mp ***)bmalloc(regen_left->num_var_gps * sizeof(comp_mp **));
        for (j = 0; j < regen_left->num_var_gps; j++)
        {
          regen_left->coeff_mp[i][j] = (comp_mp **)bmalloc(regen_left->degrees[i][j] * sizeof(comp_mp *));
          for (k = 0; k < regen_left->degrees[i][j]; k++)
          {
            regen_left->coeff_mp[i][j][k] = (comp_mp *)bmalloc(regen_left->num_variables * sizeof(comp_mp));
            for (l = 0; l < regen_left->num_variables; l++)
            {
              init_mp(regen_left->coeff_mp[i][j][k][l]);
              // setup real part
              mpf_set_str(regen_left->coeff_mp[i][j][k][l]->r, &(*regenStr)[count], base);
              count += 1 + strlen(&(*regenStr)[count]);
              // setup imag part
              mpf_set_str(regen_left->coeff_mp[i][j][k][l]->i, &(*regenStr)[count], base);
              count += 1 + strlen(&(*regenStr)[count]);
            }
          }
        }
      }
      // setup mainCoeff_mp
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        { // setup real part
          mpf_set_str(regen_left->mainCoeff_mp->entry[i][j].r, &(*regenStr)[count], base);
          count += 1 + strlen(&(*regenStr)[count]);
          // setup imag part
          mpf_set_str(regen_left->mainCoeff_mp->entry[i][j].i, &(*regenStr)[count], base);
          count += 1 + strlen(&(*regenStr)[count]);
        }
      // setup main_homVar_mp
      for (i = 0; i < regen_left->num_variables; i++)
      { // setup real part
        mpf_set_str(regen_left->main_homVar_mp->coord[i].r, &(*regenStr)[count], base);
        count += 1 + strlen(&(*regenStr)[count]);
        // setup imag part
        mpf_set_str(regen_left->main_homVar_mp->coord[i].i, &(*regenStr)[count], base);
        count += 1 + strlen(&(*regenStr)[count]);
      }

      if (freeStr)
        free(*regenStr);
    }
    else
    { // setup _d, _mp, _rat
      init_d(regen_left->gamma_d);
      init_mp2(regen_left->gamma_mp, regen_left->curr_precision);
      mpq_init(regen_left->gamma_rat[0]); mpq_init(regen_left->gamma_rat[1]);
      init_mat_d(regen_left->patchCoeff_d, regen_right->patch_rows, regen_right->patch_cols);
      init_mat_mp2(regen_left->patchCoeff_mp, regen_right->patch_rows, regen_right->patch_cols, regen_left->curr_precision);
      init_mat_rat(regen_left->patchCoeff_rat, regen_right->patch_rows, regen_right->patch_cols);
      init_mat_d(regen_left->mainCoeff_d, regen_left->num_funcs, regen_left->num_variables);
      init_mat_mp2(regen_left->mainCoeff_mp, regen_left->num_funcs, regen_left->num_variables, regen_left->curr_precision);
      init_mat_rat(regen_left->mainCoeff_rat, regen_left->num_funcs, regen_left->num_variables);
      init_vec_d(regen_left->main_homVar_d, regen_left->num_variables);
      init_vec_mp(regen_left->main_homVar_mp, regen_left->num_variables);
      init_vec_rat(regen_left->main_homVar_rat, regen_left->num_variables);

      regen_left->patchCoeff_d->rows = regen_left->patchCoeff_mp->rows = regen_right->patch_rows;
      regen_left->patchCoeff_d->cols = regen_left->patchCoeff_mp->cols = regen_right->patch_cols;
      regen_left->mainCoeff_d->rows = regen_left->mainCoeff_mp->rows = regen_left->num_funcs;
      regen_left->mainCoeff_d->cols = regen_left->mainCoeff_mp->cols = regen_left->num_variables;
      regen_left->main_homVar_d->size = regen_left->main_homVar_mp->size = regen_left->num_variables;

      // setup gamma
      count = 0;
      mpq_set_str(regen_left->gamma_rat[0], &(*regenStr)[count], base);
      mpq_canonicalize(regen_left->gamma_rat[0]);
      count += 1 + strlen(&(*regenStr)[count]);
      mpq_set_str(regen_left->gamma_rat[1], &(*regenStr)[count], base);
      mpq_canonicalize(regen_left->gamma_rat[1]);
      count += 1 + strlen(&(*regenStr)[count]);

      // setup gamma_d & gamma_mp
      rat_to_d(regen_left->gamma_d, regen_left->gamma_rat);
      rat_to_mp(regen_left->gamma_mp, regen_left->gamma_rat);

      // setup patchCoeff
      for (i = 0; i < regen_right->patch_rows; i++)
        for (j = 0; j < regen_right->patch_cols; j++)
        { // setup real part
          mpq_set_str(regen_left->patchCoeff_rat[i][j][0], &(*regenStr)[count], base);
          mpq_canonicalize(regen_left->patchCoeff_rat[i][j][0]);
          count += 1 + strlen(&(*regenStr)[count]);
          // setup imag part
          mpq_set_str(regen_left->patchCoeff_rat[i][j][1], &(*regenStr)[count], base);
          mpq_canonicalize(regen_left->patchCoeff_rat[i][j][1]);
          count += 1 + strlen(&(*regenStr)[count]);

          // setup _d & _mp
          rat_to_d(&regen_left->patchCoeff_d->entry[i][j], regen_left->patchCoeff_rat[i][j]);
          rat_to_mp(&regen_left->patchCoeff_mp->entry[i][j], regen_left->patchCoeff_rat[i][j]);
        }
      // setup coeff
      regen_left->coeff_d = (comp_d ****)bmalloc(regen_left->num_funcs * sizeof(comp_d ***));
      regen_left->coeff_mp = (comp_mp ****)bmalloc(regen_left->num_funcs * sizeof(comp_mp ***));
      regen_left->coeff_rat = (mpq_t *****)bmalloc(regen_left->num_funcs * sizeof(mpq_t ****));
      for (i = 0; i < regen_left->num_funcs; i++)
      {
        regen_left->coeff_d[i] = (comp_d ***)bmalloc(regen_left->num_var_gps * sizeof(comp_d **));
        regen_left->coeff_mp[i] = (comp_mp ***)bmalloc(regen_left->num_var_gps * sizeof(comp_mp **));
        regen_left->coeff_rat[i] = (mpq_t ****)bmalloc(regen_left->num_var_gps * sizeof(mpq_t ***));
        for (j = 0; j < regen_left->num_var_gps; j++)
        {
          regen_left->coeff_d[i][j] = (comp_d **)bmalloc(regen_left->degrees[i][j] * sizeof(comp_d *));
          regen_left->coeff_mp[i][j] = (comp_mp **)bmalloc(regen_left->degrees[i][j] * sizeof(comp_mp *));
          regen_left->coeff_rat[i][j] = (mpq_t ***)bmalloc(regen_left->degrees[i][j] * sizeof(mpq_t **));
          for (k = 0; k < regen_left->degrees[i][j]; k++)
          {
            regen_left->coeff_d[i][j][k] = (comp_d *)bmalloc(regen_left->num_variables * sizeof(comp_d));
            regen_left->coeff_mp[i][j][k] = (comp_mp *)bmalloc(regen_left->num_variables * sizeof(comp_mp));
            regen_left->coeff_rat[i][j][k] = (mpq_t **)bmalloc(regen_left->num_variables * sizeof(mpq_t *));
            for (l = 0; l < regen_left->num_variables; l++)
            {
              regen_left->coeff_rat[i][j][k][l] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
              mpq_init(regen_left->coeff_rat[i][j][k][l][0]); mpq_init(regen_left->coeff_rat[i][j][k][l][1]);

              // setup real part
              mpq_set_str(regen_left->coeff_rat[i][j][k][l][0], &(*regenStr)[count], base);
              mpq_canonicalize(regen_left->coeff_rat[i][j][k][l][0]);
              count += 1 + strlen(&(*regenStr)[count]);
              // setup imag part
              mpq_set_str(regen_left->coeff_rat[i][j][k][l][1], &(*regenStr)[count], base);
              mpq_canonicalize(regen_left->coeff_rat[i][j][k][l][1]);
              count += 1 + strlen(&(*regenStr)[count]);

              // setup _d & _mp
              rat_to_d(regen_left->coeff_d[i][j][k][l], regen_left->coeff_rat[i][j][k][l]);
              init_mp2(regen_left->coeff_mp[i][j][k][l], regen_left->curr_precision);
              rat_to_mp(regen_left->coeff_mp[i][j][k][l], regen_left->coeff_rat[i][j][k][l]);
            }
          }
        }
      }
      // setup mainCoeff
      for (i = 0; i < regen_left->num_funcs; i++)
        for (j = 0; j < regen_left->num_variables; j++)
        { // setup real part
          mpq_set_str(regen_left->mainCoeff_rat[i][j][0], &(*regenStr)[count], base);
          mpq_canonicalize(regen_left->mainCoeff_rat[i][j][0]);
          count += 1 + strlen(&(*regenStr)[count]);
          // setup imag part
          mpq_set_str(regen_left->mainCoeff_rat[i][j][1], &(*regenStr)[count], base);
          mpq_canonicalize(regen_left->mainCoeff_rat[i][j][1]);
          count += 1 + strlen(&(*regenStr)[count]);

          // setup _d & _mp
          rat_to_d(&regen_left->mainCoeff_d->entry[i][j], regen_left->mainCoeff_rat[i][j]);
          rat_to_mp(&regen_left->mainCoeff_mp->entry[i][j], regen_left->mainCoeff_rat[i][j]);
        }
      // setup mainCoeff
      for (i = 0; i < regen_left->num_variables; i++)
      { // setup real part
        mpq_set_str(regen_left->main_homVar_rat[i][0], &(*regenStr)[count], base);
        mpq_canonicalize(regen_left->main_homVar_rat[i][0]);
        count += 1 + strlen(&(*regenStr)[count]);
        // setup imag part
        mpq_set_str(regen_left->main_homVar_rat[i][1], &(*regenStr)[count], base);
        mpq_canonicalize(regen_left->main_homVar_rat[i][1]);
        count += 1 + strlen(&(*regenStr)[count]);

        // setup _d & _mp
        rat_to_d(&regen_left->main_homVar_d->coord[i], regen_left->main_homVar_rat[i]);
        rat_to_mp(&regen_left->main_homVar_mp->coord[i], regen_left->main_homVar_rat[i]);
      }

      if (freeStr)
        free(*regenStr);
    }
  }

  return;
}

void cp_regenLevel_int(void *rL_out, void *rL_in, int MPType, int curr_prec, char **rlStr, comp_d **coeff_d, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - rL_out is _int, rL_in is _t            *
*                otherwise - rL_out is _t, rL_in is _int        *
* RETURN VALUES:                                                *
* NOTES: stores a copy of rL_in to rL_out                       *
\***************************************************************/
{
  int i, j, count, base = 10;

  if (inType == 0)
  { // _t to _int so that it can be sent
    regenLevel_t *rL_right = (regenLevel_t *)rL_in;
    regenLevel_t_int *rL_left = (regenLevel_t_int *)rL_out;

    // setup the basic info
    rL_left->level = rL_right->level;
    rL_left->depth = rL_right->depth;
    rL_left->useIntrinsicSlice = rL_right->useIntrinsicSlice;

    if (rL_right->useIntrinsicSlice)
    { // setup the rest
      if (MPType == 0)
      { // use _d
        rL_left->B_rows = rL_right->B_d->rows;
        rL_left->B_cols = rL_right->B_d->cols;
        rL_left->p_size = rL_right->p_d->size;

        // count the total number of comp_d that are to be sent
        rL_left->num_comp_d = rL_left->B_rows * rL_left->B_cols + rL_left->p_size;
        // setup coeff_d
        *coeff_d = (comp_d *)bmalloc(rL_left->num_comp_d * sizeof(comp_d));
        count = 0;
        for (i = 0; i < rL_left->B_rows; i++)
          for (j = 0; j < rL_left->B_cols; j++)
          {
            set_d((*coeff_d)[count], &rL_right->B_d->entry[i][j]);
            count++;
          }
        for (i = 0; i < rL_left->p_size; i++)
        {
          set_d((*coeff_d)[count], &rL_right->p_d->coord[i]);
          count++;
        }
        // setup totalLength
        rL_left->totalLength = 0;
      }
      else if (MPType == 1)
      { // use _mp
        int sizeR, sizeI;
        char *strR = NULL, *strI = NULL;
  
        rL_left->B_rows = rL_right->B_mp->rows;
        rL_left->B_cols = rL_right->B_mp->cols;
        rL_left->p_size = rL_right->p_mp->size;
        rL_left->num_comp_d = 0; // no comp_d setup

        // setup rlStr
        *rlStr = NULL;
        count = 0;
        // copy B_mp to rlStr
        for (i = 0; i < rL_left->B_rows; i++)
          for (j = 0; j < rL_left->B_cols; j++)
          { // real part
            strR = mpf_to_str(rL_right->B_mp->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(rL_right->B_mp->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rlStr
            *rlStr = (char *)brealloc(*rlStr, (count + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rlStr)[count], strR);
            strcpy(&(*rlStr)[count + sizeR], strI);
            // update count
            count += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // copy p_mp to rlStr
        for (i = 0; i < rL_left->p_size; i++)
        { // real part
          strR = mpf_to_str(rL_right->p_mp->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(rL_right->p_mp->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update rlStr
          *rlStr = (char *)brealloc(*rlStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*rlStr)[count], strR);
          strcpy(&(*rlStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
        // setup totalLength
        rL_left->totalLength = count;
      }
      else
      { // use _rat
        int sizeR, sizeI;
        char *strR = NULL, *strI = NULL;

        rL_left->B_rows = rL_right->B_d->rows;
        rL_left->B_cols = rL_right->B_d->cols;
        rL_left->p_size = rL_right->p_d->size;
        rL_left->num_comp_d = 0; // no comp_d setup

        // setup rlStr
        *rlStr = NULL;
        count = 0;
        // copy B_rat to rlStr
        for (i = 0; i < rL_left->B_rows; i++)
          for (j = 0; j < rL_left->B_cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, rL_right->B_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, rL_right->B_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rlStr
            *rlStr = (char *)brealloc(*rlStr, (count + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rlStr)[count], strR);
            strcpy(&(*rlStr)[count + sizeR], strI);
            // update count
            count += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // copy p_rat to rlStr
        for (i = 0; i < rL_left->p_size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, rL_right->p_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, rL_right->p_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update rlStr
          *rlStr = (char *)brealloc(*rlStr, (count + sizeR + sizeI) * sizeof(char));
          strcpy(&(*rlStr)[count], strR);
          strcpy(&(*rlStr)[count + sizeR], strI);
          // update count
          count += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
        // setup totalLength
        rL_left->totalLength = count;
      }
    }
    else
    { // set all to 0
      rL_left->B_rows = rL_left->B_cols = rL_left->p_size = rL_left->num_comp_d = rL_left->totalLength = 0;
    }
  }
  else
  { // _int to _t so that is can be used
    regenLevel_t_int *rL_right = (regenLevel_t_int *)rL_in;
    regenLevel_t *rL_left = (regenLevel_t *)rL_out;

    // setup basic info
    rL_left->level = rL_right->level;
    rL_left->depth = rL_right->depth;
    rL_left->useIntrinsicSlice = rL_right->useIntrinsicSlice;

    // zero out other values
    rL_left->num_paths = rL_left->num_sing = rL_left->num_nonsing = rL_left->num_inf = rL_left->num_higher_dim = rL_left->num_bad = 0;

    if (rL_right->useIntrinsicSlice)
    { // setup the rest
      if (MPType == 0)
      { // setup _d
        init_mat_d(rL_left->B_d, rL_right->B_rows, rL_right->B_cols);
        init_vec_d(rL_left->p_d, rL_right->p_size);
        rL_left->B_d->rows = rL_right->B_rows;
        rL_left->B_d->cols = rL_right->B_cols;
        rL_left->p_d->size = rL_right->p_size;

        // setup B_d
        count = 0;
        for (i = 0; i < rL_right->B_rows; i++)
          for (j = 0; j < rL_right->B_cols; j++)
          {
            set_d(&rL_left->B_d->entry[i][j], (*coeff_d)[count]);
            count++;
          }
        // setup p_d
        for (i = 0; i < rL_right->p_size; i++)
        {
          set_d(&rL_left->p_d->coord[i], (*coeff_d)[count]);
          count++;
        }
        // free coeff_d
        free(*coeff_d);
      }
      else if (MPType == 1)
      { // setup _mp
        init_mat_mp(rL_left->B_mp, rL_right->B_rows, rL_right->B_cols);
        init_vec_mp(rL_left->p_mp, rL_right->p_size);

        rL_left->B_mp->rows = rL_right->B_rows;
        rL_left->B_mp->cols = rL_right->B_cols;
        rL_left->p_mp->size = rL_right->p_size;

        // setup B_mp
        count = 0;
        for (i = 0; i < rL_right->B_rows; i++)
          for (j = 0; j < rL_right->B_cols; j++)
          { // setup real part
            mpf_set_str(rL_left->B_mp->entry[i][j].r, &(*rlStr)[count], base);
            count += 1 + strlen(&(*rlStr)[count]);
            // setup imag part
            mpf_set_str(rL_left->B_mp->entry[i][j].i, &(*rlStr)[count], base);
            count += 1 + strlen(&(*rlStr)[count]);
          }
        // setup p_mp
        for (i = 0; i < rL_right->p_size; i++)
        { // setup real part
          mpf_set_str(rL_left->p_mp->coord[i].r, &(*rlStr)[count], base);
          count += 1 + strlen(&(*rlStr)[count]);
          // setup imag part
          mpf_set_str(rL_left->p_mp->coord[i].i, &(*rlStr)[count], base);
          count += 1 + strlen(&(*rlStr)[count]);
        }
        if (freeStr)
          free(*rlStr);
      }
      else
      { // setup _d, _mp & _rat
        init_mat_d(rL_left->B_d, rL_right->B_rows, rL_right->B_cols);
        init_mat_mp2(rL_left->B_mp, rL_right->B_rows, rL_right->B_cols, curr_prec);
        init_mat_rat(rL_left->B_rat, rL_right->B_rows, rL_right->B_cols);
        init_vec_d(rL_left->p_d, rL_right->p_size);
        init_vec_mp2(rL_left->p_mp, rL_right->p_size, curr_prec);
        init_vec_rat(rL_left->p_rat, rL_right->p_size);

        rL_left->B_d->rows = rL_left->B_mp->rows = rL_right->B_rows;
        rL_left->B_d->cols = rL_left->B_mp->cols = rL_right->B_cols;
        rL_left->p_d->size = rL_left->p_mp->size = rL_right->p_size;

        // setup B
        count = 0;
        for (i = 0; i < rL_right->B_rows; i++)
          for (j = 0; j < rL_right->B_cols; j++)
          { // setup real part
            mpq_set_str(rL_left->B_rat[i][j][0], &(*rlStr)[count], base);
            mpq_canonicalize(rL_left->B_rat[i][j][0]);
            count += 1 + strlen(&(*rlStr)[count]);
            // setup imag part
            mpq_set_str(rL_left->B_rat[i][j][1], &(*rlStr)[count], base);
            mpq_canonicalize(rL_left->B_rat[i][j][1]);
            count += 1 + strlen(&(*rlStr)[count]);

            // setup _d & _mp
            mpf_set_q(rL_left->B_mp->entry[i][j].r, rL_left->B_rat[i][j][0]);
            mpf_set_q(rL_left->B_mp->entry[i][j].i, rL_left->B_rat[i][j][1]);
            rL_left->B_d->entry[i][j].r = mpq_get_d(rL_left->B_rat[i][j][0]);
            rL_left->B_d->entry[i][j].i = mpq_get_d(rL_left->B_rat[i][j][1]);
          }
        // setup p
        for (i = 0; i < rL_right->p_size; i++)
        { // setup real part
          mpq_set_str(rL_left->p_rat[i][0], &(*rlStr)[count], base);
          mpq_canonicalize(rL_left->p_rat[i][0]);
          count += 1 + strlen(&(*rlStr)[count]);
          // setup imag part
          mpq_set_str(rL_left->p_rat[i][1], &(*rlStr)[count], base);
          mpq_canonicalize(rL_left->p_rat[i][1]);
          count += 1 + strlen(&(*rlStr)[count]);

          // setup _d & _mp
          mpf_set_q(rL_left->p_mp->coord[i].r, rL_left->p_rat[i][0]);
          mpf_set_q(rL_left->p_mp->coord[i].i, rL_left->p_rat[i][1]);
          rL_left->p_d->coord[i].r = mpq_get_d(rL_left->p_rat[i][0]);
          rL_left->p_d->coord[i].i = mpq_get_d(rL_left->p_rat[i][1]);
        }
        if (freeStr)
          free(*rlStr);
      }
    }
    else
    { // 'zero' out everything
      rL_left->B_d->rows = rL_left->B_mp->rows = rL_left->B_d->cols = rL_left->B_mp->cols = rL_left->p_d->size = rL_left->p_mp->size = 0;
      rL_left->B_rat = NULL;
      rL_left->p_rat = NULL;
    }
  }

  return;
}

void bcast_regen_t(regen_t *regen, int MPType, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts regen                                       *
\***************************************************************/
{
  MPI_Datatype mpi_comp_d, mpi_regen_int;
  regen_t_int regen_int;
  int *regenDegrees = NULL, *PPD_type = NULL, *PPD_size = NULL, *instCount = NULL;
  char *regenStr = NULL;
  comp_d *regenComp = NULL;

  // create the datatype mpi_comp_d & mpi_regen_int
  create_comp_d(&mpi_comp_d);
  create_regen_int(&mpi_regen_int);

  if (my_id == headnode)
  { // setup regen_int & the other structures
    cp_regen_int(&regen_int, regen, MPType, &regenDegrees, &PPD_type, &PPD_size, &regenStr, &regenComp, 1, &instCount, 0);

    // broadcast everything
    MPI_Bcast(&regen_int, 1, mpi_regen_int, headnode, MPI_COMM_WORLD);
    MPI_Bcast(regenDegrees, regen_int.num_funcs * regen_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(PPD_type, regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(PPD_size, regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(regenStr, regen_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    MPI_Bcast(regenComp, regen_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    MPI_Bcast(instCount, regen_int.numInts, MPI_INT, headnode, MPI_COMM_WORLD);

    // clear memory
    free(regenDegrees);
    free(PPD_type);
    free(PPD_size);
    free(regenStr);
    free(regenComp);
  }
  else
  { // recv regen_int
    MPI_Bcast(&regen_int, 1, mpi_regen_int, headnode, MPI_COMM_WORLD);
    // recv regenDegrees
    regenDegrees = (int *)bmalloc(regen_int.num_funcs * regen_int.num_var_gps * sizeof(int));
    MPI_Bcast(regenDegrees, regen_int.num_funcs * regen_int.num_var_gps, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv PPD_type & PPD_size
    PPD_type = (int *)bmalloc((regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp) * sizeof(int));
    PPD_size = (int *)bmalloc((regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp) * sizeof(int));
    MPI_Bcast(PPD_type, regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(PPD_size, regen_int.PPD_int.num_hom_var_gp + regen_int.PPD_int.num_var_gp, MPI_INT, headnode, MPI_COMM_WORLD);
    // recv regenStr
    regenStr = (char *)bmalloc(regen_int.totalLength * sizeof(char));
    MPI_Bcast(regenStr, regen_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv regenComp
    regenComp = (comp_d *)bmalloc(regen_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(regenComp, regen_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);
    // recv instCount
    instCount = (int *)bmalloc(regen_int.numInts * sizeof(int));
    MPI_Bcast(instCount, regen_int.numInts, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup regen - clears all memory!
    cp_regen_int(regen, &regen_int, MPType, &regenDegrees, &PPD_type, &PPD_size, &regenStr, &regenComp, 1, &instCount, 1);
  }

  // free mpi_comp_d & mpi_regen_int
  MPI_Type_free(&mpi_comp_d);
  MPI_Type_free(&mpi_regen_int);

  return;
}

void bcast_regenLevel_t(regenLevel_t *rL, int MPType, int curr_prec, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcasts regen                                       *
\***************************************************************/
{
  MPI_Datatype mpi_comp_d, mpi_regenLevel_int;
  regenLevel_t_int rL_int;
  char *rlStr = NULL;
  comp_d *rlComp = NULL;

  // create the datatype mpi_comp_d & mpi_regen_int
  create_comp_d(&mpi_comp_d);
  create_regenLevel_int(&mpi_regenLevel_int);

  if (my_id == headnode)
  { // setup rL_int & the other structures
    cp_regenLevel_int(&rL_int, rL, MPType, curr_prec, &rlStr, &rlComp, 1, 0);

    // broadcast everything
    MPI_Bcast(&rL_int, 1, mpi_regenLevel_int, headnode, MPI_COMM_WORLD);
    MPI_Bcast(rlStr, rL_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    MPI_Bcast(rlComp, rL_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // clear memory
    free(rlStr);
    free(rlComp);
  }
  else
  { // recv rLn_int
    MPI_Bcast(&rL_int, 1, mpi_regenLevel_int, headnode, MPI_COMM_WORLD);
    // recv rlStr
    rlStr = (char *)bmalloc(rL_int.totalLength * sizeof(char));
    MPI_Bcast(rlStr, rL_int.totalLength, MPI_CHAR, headnode, MPI_COMM_WORLD);
    // recv rlComp
    rlComp = (comp_d *)bmalloc(rL_int.num_comp_d * sizeof(comp_d));
    MPI_Bcast(rlComp, rL_int.num_comp_d, mpi_comp_d, headnode, MPI_COMM_WORLD);

    // setup rL - clears all memory!
    cp_regenLevel_int(rL, &rL_int, MPType, curr_prec, &rlStr, &rlComp, 1, 1);
  }

  return;
}

/////////////////// CREATE & SEND, RECV & STORE //////////////////

int regen_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i, j, num_input_vars, linearSize, level_num, intCount, *sendInts = NULL;
  regen_t *regen = NULL;

  // setup regen
  if (MPType == 0 || MPType == 2)
    regen = (regen_t *)ED_d;
  else
    regen = (regen_t *)ED_mp;

  // setup the current level
  level_num = regen->curr_level_num;

  // setup the number of variables that the start points have
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  if (regen->level[level_num].useIntrinsicSlice)
    num_input_vars = linearSize;
  else
    num_input_vars = regen->num_variables;

  // count the number of integers (curr_linear, curr_linear_degree) that need to be sent
  intCount = 2 * size * linearSize;
  sendInts = (int *)bmalloc(intCount * sizeof(int));

  // create the packet
  intCount = 0;
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    // find the path number
    sendPts[i].pathNum = pathNum[startNum + i]; // should be == startNum + i

    // setup sendPts[i]
    if (MPType == 0 || MPType == 2)
    { // setup sendPts[i].PD_d
      sendPts[i].prec = 52; // start in double precision
      change_size_point_d(sendPts[i].PD_d.point, num_input_vars);
      sendPts[i].PD_d.point->size = num_input_vars;
      set_one_d(sendPts[i].PD_d.time);
      fscanf(START, "\n");
      for (j = 0; j < num_input_vars; j++)
        fscanf(START, "%lf %lf;\n", &sendPts[i].PD_d.point->coord[j].r, &sendPts[i].PD_d.point->coord[j].i);
    }
    else
    { // setup sendPts[i].PD_mp
      sendPts[i].prec = (int) mpf_get_default_prec();
      change_prec_point_data_mp(&sendPts[i].PD_mp, sendPts[i].prec);
      change_size_point_mp(sendPts[i].PD_mp.point, num_input_vars);
      sendPts[i].PD_mp.point->size = num_input_vars;
      set_one_mp(sendPts[i].PD_mp.time);
      fscanf(START, "\n");
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].r, START, 10);
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].i, START, 10);
        fscanf(START, ";\n");
      }
    }

    // read in curr_linear & curr_linear_degree
    for (j = 0; j < linearSize; j++, intCount += 2)
      fscanf(START, "%d %d\n", &sendInts[intCount], &sendInts[intCount + 1]);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);
  // send sendInts to 'sendProc'
  MPI_Send(sendInts, intCount, MPI_INT, sendProc, TAG_OTHER_DATA, MPI_COMM_WORLD);

  // clear memory
  free(sendInts);
  regen = NULL;

  return size;
}

int regen_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, j, retVal, linearSize, level_num, recvProc, intCount, rankDef, finite, higherDim;
  int *recvInts = NULL, *curr_linear = NULL, *curr_linear_degree = NULL;
  MPI_Status status;
  point_d orig_d, dehom_d;
  point_mp orig_mp, dehom_mp;
  regen_t *regen = NULL;

  // initialize
  init_point_d(orig_d, 0);
  init_point_d(dehom_d, 0);
  init_point_mp(orig_mp, 0);
  init_point_mp(dehom_mp, 0);

  if (T->MPType == 0 || T->MPType == 2)
    regen = (regen_t *)ED_d;
  else 
    regen = (regen_t *)ED_mp;

  // setup the current level
  level_num = regen->curr_level_num;
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;
 
  // allocate curr_linear & curr_linear_degree
  curr_linear = (int *)bmalloc(linearSize * sizeof(int));
  curr_linear_degree = (int *)bmalloc(linearSize * sizeof(int));

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // setup recvInts
  intCount = 2 * (*numRecvPts) * linearSize;
  recvInts = (int *)bmalloc(intCount * sizeof(int));
  MPI_Recv(recvInts, intCount, MPI_INT, recvProc, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

  // store the data
  intCount = 0;
  for (i = 0; i < *numRecvPts; i++)
  { // setup curr_linear & curr_linear_degree
    for (j = 0; j < linearSize; j++, intCount += 2)
    {
      curr_linear[j] = recvInts[intCount];
      curr_linear_degree[j] = recvInts[intCount + 1];
    }
    // point to curr_linear & curr_linear_degree
    regen->curr_linear = curr_linear;
    regen->curr_linear_degree = curr_linear_degree;

    // print the footer to OUT and print all info to RAWOUT
    fprintf(OUT, "Path number: %d (ID: %d)\n", (*recvPts)[i].pathNum, recvProc);

    // setup rankDef, finite & higherDim 
    retVal = (*recvPts)[i].retVal;
    rankDef = 0; finite = 1; higherDim = 0;
    if (retVal == retVal_going_to_infinity || retVal == retVal_security_max)
      finite = 0;
    if (retVal == retVal_higher_dim) 
    {
      retVal = 0;
      higherDim = 1;
    }
    if (retVal == DO_NOT_MOVE_TO_NEXT)
    {
      retVal = 0;
      rankDef = 1;
    }
    if (retVal == MOVE_TO_NEXT)
      retVal = 0;

    // compute dehom & orig
    regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &(*recvPts)[i].prec, (*recvPts)[i].PD_d.point, (*recvPts)[i].PD_mp.point, (*recvPts)[i].prec, regen, regen);

    printRegenTrackFooter(retVal, regen, level_num, &(*recvPts)[i], rankDef, finite, higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, OTHER, FAIL, OTHER2, T, trackCount, level_num + 1 == regen->num_levels);
  }

  // clear memory
  free(recvInts);
  free(curr_linear);
  free(curr_linear_degree);
  regen = NULL;
  clear_point_d(orig_d);
  clear_point_d(dehom_d);
  clear_point_mp(orig_mp);
  clear_point_mp(dehom_mp);

  return recvProc;
}

int regen_create_send_packet_move(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen slice moving    *
\***************************************************************/
{
  int i, j, num_input_vars, linearSize, level_num, intCount, *sendInts = NULL;
  regen_t *regen = NULL;

  // setup regen
  if (MPType == 0 || MPType == 2)
    regen = (regen_t *)ED_d;
  else
    regen = (regen_t *)ED_mp;

  // setup the current level
  level_num = regen->curr_level_num;
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  // setup the number of variables that the start points have
  if (regen->level[level_num].useIntrinsicSlice)
    num_input_vars = linearSize;
  else
    num_input_vars = regen->num_variables;

  // count the number of integers (curr_linear, curr_linear_degree) that need to be sent
  intCount = 2 * size * linearSize;
  sendInts = (int *)bmalloc(intCount * sizeof(int));

  // create the packet
  intCount = 0;
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    // read in the path number
    fscanf(START, "%d\n", &sendPts[i].pathNum);

    // setup sendPts[i]
    if (MPType == 0 || MPType == 2)
    { // setup sendPts[i].PD_d
      sendPts[i].prec = 52; // start in double precision
      change_size_point_d(sendPts[i].PD_d.point, num_input_vars);
      sendPts[i].PD_d.point->size = num_input_vars;
      set_one_d(sendPts[i].PD_d.time);
      for (j = 0; j < num_input_vars; j++)
        fscanf(START, "%lf %lf;\n", &sendPts[i].PD_d.point->coord[j].r, &sendPts[i].PD_d.point->coord[j].i);
    }
    else
    { // setup sendPts[i].PD_mp
      sendPts[i].prec = (int) mpf_get_default_prec();
      change_prec_point_data_mp(&sendPts[i].PD_mp, sendPts[i].prec);
      change_size_point_mp(sendPts[i].PD_mp.point, num_input_vars);
      sendPts[i].PD_mp.point->size = num_input_vars;
      set_one_mp(sendPts[i].PD_mp.time);
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].r, START, 10);
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].i, START, 10);
        fscanf(START, ";\n");
      }
    }

    // read in curr_linear & curr_linear_degree
    for (j = 0; j < linearSize; j++, intCount += 2)
      fscanf(START, "%d %d\n", &sendInts[intCount], &sendInts[intCount + 1]);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);
  // send sendInts to 'sendProc'
  MPI_Send(sendInts, intCount, MPI_INT, sendProc, TAG_OTHER_DATA, MPI_COMM_WORLD);

  // clear memory
  free(sendInts);
  regen = NULL;

  return size;
}

int regen_recv_store_packet_move(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, j, retVal, level_num, linearSize, recvProc, intCount, rankDef, finite, higherDim;
  int *recvInts = NULL, *curr_linear = NULL, *curr_linear_degree = NULL;
  MPI_Status status;
  point_d orig_d, dehom_d;
  point_mp orig_mp, dehom_mp;
  regen_t *regen = NULL;

  // initialize
  init_point_d(orig_d, 0);
  init_point_d(dehom_d, 0);
  init_point_mp(orig_mp, 0);
  init_point_mp(dehom_mp, 0);

  if (T->MPType == 0 || T->MPType == 2)
  { // setup regen 
    regen = (regen_t *)ED_d;
  }
  else
  { // setup regen
    regen = (regen_t *)ED_mp;
  }
  // setup the current level
  level_num = regen->curr_level_num;
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;
  
  // allocate curr_linear & curr_linear_degree
  curr_linear = (int *)bmalloc(linearSize * sizeof(int));
  curr_linear_degree = (int *)bmalloc(linearSize * sizeof(int));

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // setup recvInts
  intCount = 2 * (*numRecvPts) * linearSize;
  recvInts = (int *)bmalloc(intCount * sizeof(int));
  MPI_Recv(recvInts, intCount, MPI_INT, recvProc, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

  // store the data
  intCount = 0;
  for (i = 0; i < *numRecvPts; i++)
  { // setup curr_linear & curr_linear_degree
    for (j = 0; j < linearSize; j++, intCount += 2)
    {
      curr_linear[j] = recvInts[intCount];
      curr_linear_degree[j] = recvInts[intCount + 1];
    }
    // point to curr_linear & curr_linear_degree
    regen->curr_linear = curr_linear;
    regen->curr_linear_degree = curr_linear_degree;

    // print the footer to OUT and print all info to RAWOUT
    fprintf(OUT, "Path number: %d (ID: %d)\n", (*recvPts)[i].pathNum, recvProc);

    // setup rankDef, finite & higherDim 
    retVal = (*recvPts)[i].retVal;
    rankDef = 0; finite = 1; higherDim = 0;
    if (retVal == retVal_going_to_infinity || retVal == retVal_security_max)
      finite = 0;
    if (retVal == retVal_higher_dim)
    {
      retVal = 0;
      higherDim = 1;
    }
    if (retVal == DO_NOT_MOVE_TO_NEXT)
    {
      retVal = 0;
      rankDef = 1;
    }
    if (retVal == MOVE_TO_NEXT)
      retVal = 0;

    // compute dehom & orig
    regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &(*recvPts)[i].prec, (*recvPts)[i].PD_d.point, (*recvPts)[i].PD_mp.point, (*recvPts)[i].prec, regen, regen);

    // print the required data
    retVal = printRegenLinearTrackFooter(retVal, regen, level_num, (*recvPts)[i].pathNum, &(*recvPts)[i], rankDef, finite, higherDim, orig_d, orig_mp, dehom_d, dehom_mp, OUT, RAWOUT, OTHER, FAIL, T);

    // save the retVal
    rV[(*recvPts)[i].pathNum] = retVal;
  }

  // clear memory
  free(recvInts);
  free(curr_linear);
  free(curr_linear_degree);
  regen = NULL;
  clear_point_d(orig_d);
  clear_point_d(dehom_d);
  clear_point_mp(orig_mp);
  clear_point_mp(dehom_mp);

  return recvProc;
}

int regen_useSharpener(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether to sharpen or not                      *
* NOTES:                                                        *
\***************************************************************/
{
  int rV = 0;
  regen_t *regen = NULL;

  if (ED_d != NULL)
  { // setup regen 
    regen = (regen_t *)ED_d;
  }
  else
  { // setup regen
    regen = (regen_t *)ED_mp;
  }

  if (regen->curr_level_num + 1 == regen->num_levels && retVal == 0 && sharpenDigits > 0)
    rV = 1;
  else
    rV = 0;

  // clear
  regen = NULL;

  return rV;
}

int regen_useSharpener_move(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether to sharpen or not when moving slices   *
* NOTES:                                                        *
\***************************************************************/
{
  return 0;
}

int regen_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, j, level_num, linearSize, intCount, *recvInts = NULL, *curr_linear = NULL, *curr_linear_degree = NULL;
  MPI_Status status;
  regen_t *regen = NULL;

  // setup regen
  if (T->MPType == 0 || T->MPType == 2)
    regen = (regen_t *)ED_d;
  else
    regen = (regen_t *)ED_mp;

  level_num = regen->curr_level_num;
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;

  // recv next packet of start points
  send_recv_endgame_data_t(startPts, numStartPts, T->MPType, headnode, 0);

  // setup recvInts
  intCount = 2 * (*numStartPts) * linearSize;
  recvInts = (int *)bmalloc(intCount * sizeof(int));
  MPI_Recv(recvInts, intCount, MPI_INT, headnode, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

  // setup endPts
  if (*numEndPts != *numStartPts)
  { // clear endPts
    for (i = *numEndPts - 1; i >= 0; i--)
    {
      clear_endgame_data(&(*endPts)[i]);
    }

    // set the number to reallocate
    *numEndPts = *numStartPts;

    // reallocate
    *endPts = (endgame_data_t *)brealloc(*endPts, *numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < *numEndPts; i++)
    {
      init_endgame_data(&(*endPts)[i], T->Precision);
    }
  }

  // make sure that we have paths to track
  if (*numStartPts > 0) 
  { // allocate curr_linear & curr_linear_degree
    curr_linear = (int *)bmalloc(linearSize * sizeof(int));
    curr_linear_degree = (int *)bmalloc(linearSize * sizeof(int));

    // track the paths
    intCount = 0;
    for (i = 0; i < *numStartPts; i++)
    { // setup curr_linear & curr_linear_degree
      for (j = 0; j < linearSize; j++, intCount += 2)
      {
        curr_linear[j] = recvInts[intCount];
        curr_linear_degree[j] = recvInts[intCount + 1];
      }
      // point to curr_linear & curr_linear_degree
      regen->curr_linear = curr_linear;
      regen->curr_linear_degree = curr_linear_degree;

      // track the ith path
      regen_worker_track_path(&(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener);
    }

    // send the packet back
    send_recv_endgame_data_t(endPts, numEndPts, T->MPType, headnode, 1);

    // send recvInts back
    MPI_Send(recvInts, intCount, MPI_INT, headnode, TAG_OTHER_DATA, MPI_COMM_WORLD);
  }

  // clear memory
  free(recvInts);
  free(curr_linear);
  free(curr_linear_degree);
  regen = NULL;

  return 0;
}

int regen_recv_linear_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, j, level_num, linearSize, intCount, *recvInts = NULL, *curr_linear = NULL, *curr_linear_degree = NULL;
  MPI_Status status;
  regen_t *regen = NULL;

  // setup regen
  if (T->MPType == 0 || T->MPType == 2)
    regen = (regen_t *)ED_d;
  else
    regen = (regen_t *)ED_mp;

  level_num = regen->curr_level_num;
  linearSize = regen->level[level_num].level + regen->level[level_num].depth;

  // recv next packet of start points
  send_recv_endgame_data_t(startPts, numStartPts, T->MPType, headnode, 0);

  // setup recvInts
  intCount = 2 * (*numStartPts) * linearSize;
  recvInts = (int *)bmalloc(intCount * sizeof(int));
  MPI_Recv(recvInts, intCount, MPI_INT, headnode, TAG_OTHER_DATA, MPI_COMM_WORLD, &status);

  // setup endPts
  if (*numEndPts != *numStartPts)
  { // clear endPts
    for (i = *numEndPts - 1; i >= 0; i--)
    {
      clear_endgame_data(&(*endPts)[i]);
    }

    // set the number to reallocate
    *numEndPts = *numStartPts;

    // reallocate
    *endPts = (endgame_data_t *)brealloc(*endPts, *numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < *numEndPts; i++)
    {
      init_endgame_data(&(*endPts)[i], T->Precision);
    }
  }

  // make sure that we have paths to track
  if (*numStartPts > 0)
  { // allocate curr_linear & curr_linear_degree
    curr_linear = (int *)bmalloc(linearSize * sizeof(int));
    curr_linear_degree = (int *)bmalloc(linearSize * sizeof(int));

    // track the paths
    intCount = 0;
    for (i = 0; i < *numStartPts; i++)
    { // setup curr_linear & curr_linear_degree
      for (j = 0; j < linearSize; j++, intCount += 2)
      {
        curr_linear[j] = recvInts[intCount];
        curr_linear_degree[j] = recvInts[intCount + 1];
      }
      // point to curr_linear & curr_linear_degree
      regen->curr_linear = curr_linear;
      regen->curr_linear_degree = curr_linear_degree;

      // track the ith path
      regen_worker_linear_track_path(&(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener);
    }

    // send the packet back
    send_recv_endgame_data_t(endPts, numEndPts, T->MPType, headnode, 1);
    // send recvInts back
    MPI_Send(recvInts, intCount, MPI_INT, headnode, TAG_OTHER_DATA, MPI_COMM_WORLD);
  }

  // clear memory
  free(recvInts);
  free(curr_linear);
  free(curr_linear_degree);
  regen = NULL;

  return 0;
}

int regen_worker_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: worker process to track the path using regeneration    *
\***************************************************************/
{
  int rankType = 0, rankDef = 0, finite = 1, higherDim = 0;
  double smallest, largest;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = &regen_dehom;
  regen_t *regen = NULL;

  if (T->MPType == 0 || T->MPType == 2)
    regen = (regen_t *)CD_d;
  else
    regen = (regen_t *)CD_mp;

  point_d dehom_d, orig_d, orig_last_d;
  point_mp dehom_mp, orig_mp, orig_last_mp;
  init_point_d(dehom_d, 0); init_point_d(orig_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0); init_point_mp(orig_last_mp, 0);

  if (T->MPType == 2)
  { // track the path as expected
    worker_track_path(startPt, endPt, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
  }
  else
  { // track the path and determine if it is rank deficient
    worker_track_path_rank(rankType, &rankDef, NULL, &smallest, &largest, startPt, endPt, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
  }

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_last_d, orig_last_mp, &endPt->last_approx_prec, endPt->last_approx_d, endPt->last_approx_mp, endPt->last_approx_prec, CD_d, CD_mp);
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &endPt->prec, endPt->PD_d.point, endPt->PD_mp.point, endPt->prec, CD_d, CD_mp);

  // compute the retVal
  if (T->MPType == 0 || T->MPType == 2)
    endPt->retVal = computeRegenRetval(regen, regen->curr_level_num, endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, regen->curr_level_num + 1 == regen->num_levels, ((square_system_eval_data_d *)regen->square_d)->Prog);
  else
    endPt->retVal = computeRegenRetval(regen, regen->curr_level_num, endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, regen->curr_level_num + 1 == regen->num_levels, ((square_system_eval_data_mp *)regen->square_mp)->Prog);

  if (endPt->retVal == 0)
  { // we need to check the endpoint
    regenSortEndpoint_basic(1, &rankDef, &finite, &higherDim, endPt->condition_number, regen, regen->curr_level_num, endPt->pathNum, T, OUT, &endPt->PD_d, &endPt->PD_mp, endPt->prec, endPt->last_approx_d, endPt->last_approx_mp, endPt->last_approx_prec, regen->curr_linear, regen->curr_linear_degree, eval_d, eval_mp, change_prec);
 
    // determine what to tell the head node about this path
    if (!finite)
      endPt->retVal = retVal_going_to_infinity;
    else if (higherDim)
      endPt->retVal = retVal_higher_dim;
    else if (rankDef == 1)
      endPt->retVal = DO_NOT_MOVE_TO_NEXT;
    else
      endPt->retVal = MOVE_TO_NEXT;
  }

  clear_point_d(dehom_d); clear_point_d(orig_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_mp); clear_point_mp(orig_last_mp);

  regen = NULL;

  return 0;
}

int regen_worker_linear_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: worker process to track the path using regeneration    *
\***************************************************************/
{
  int rankType = 0, rankDef = 0, finite = 1, higherDim = 0;
  double smallest, largest;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
  regen_t *regen = NULL;

  if (T->MPType == 0 || T->MPType == 2)
    regen = (regen_t *)CD_d;
  else
    regen = (regen_t *)CD_mp;

  point_d dehom_d, orig_d, orig_last_d;
  point_mp dehom_mp, orig_mp, orig_last_mp;
  init_point_d(dehom_d, 0); init_point_d(orig_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_mp, 0); init_point_mp(orig_last_mp, 0);

  if (regen->num_var_gps == 1)
  { // all slices are the same general setup, so every path should track correctly!
    find_dehom = &zero_dehom;
  }
  else
  { // points could go to infinity
    find_dehom = &regen_dehom;
  }

  if (T->MPType == 2)
  { // track the path as expected
    worker_track_path(startPt, endPt, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
  }
  else
  { // track the path and determine if it is rank deficient
    worker_track_path_rank(rankType, &rankDef, NULL, &smallest, &largest, startPt, endPt, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
  }

  // compute orig & dehom
  regen_dehom2(dehom_d, dehom_mp, orig_last_d, orig_last_mp, &endPt->last_approx_prec, endPt->last_approx_d, endPt->last_approx_mp, endPt->last_approx_prec, CD_d, CD_mp);
  regen_dehom2(dehom_d, dehom_mp, orig_d, orig_mp, &endPt->prec, endPt->PD_d.point, endPt->PD_mp.point, endPt->prec, CD_d, CD_mp);

  // compute the retVal
  if (T->MPType == 0 || T->MPType == 2)
    endPt->retVal = computeRegenRetval(regen, regen->curr_level_num, endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, 0, NULL);
  else
    endPt->retVal = computeRegenRetval(regen, regen->curr_level_num, endPt, orig_d, orig_mp, orig_last_d, orig_last_mp, dehom_d, dehom_mp, T, 0, NULL);

  if (endPt->retVal == 0)
  {
    if (regen->num_var_gps > 1)
    { // we need to check the endpoint
      regenSortEndpoint_basic(1, &rankDef, &finite, &higherDim, endPt->condition_number, regen, regen->curr_level_num, endPt->pathNum, T, OUT, &endPt->PD_d, &endPt->PD_mp, endPt->prec, endPt->last_approx_d, endPt->last_approx_mp, endPt->last_approx_prec, regen->curr_linear, regen->curr_linear_degree, eval_d, eval_mp, change_prec);
    }

    // determine what to tell the head node about this path
    if (!finite)
      endPt->retVal = retVal_going_to_infinity;
    else if (higherDim)
      endPt->retVal = retVal_higher_dim;
    else if (rankDef)
      endPt->retVal = DO_NOT_MOVE_TO_NEXT;
    else
      endPt->retVal = MOVE_TO_NEXT;
  }

  clear_point_d(dehom_d); clear_point_d(orig_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_mp); clear_point_mp(orig_last_mp);

  regen = NULL;

  return 0;
}

#endif

