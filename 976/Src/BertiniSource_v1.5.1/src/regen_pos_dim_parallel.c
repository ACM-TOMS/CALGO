// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"
#include "parallel.h"

void head_regen_pos_dim_MoveNextCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, int newPathCount, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int *rV, int my_id, int headnode, int num_processes);

int regen_pos_dim_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int regen_pos_dim_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int regen_pos_dim_create_send_packet_sort(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int regen_pos_dim_recv_store_packet_sort(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int regen_pos_dim_create_send_packet_prepare(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int regen_pos_dim_recv_store_packet_prepare(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));

int regen_pos_dim_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));
int regen_pos_dim_recv_sort_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));
int regen_pos_dim_recv_prepare_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

#ifdef _HAVE_MPI

void regen_pos_dim_par_track(int startCodimIndex, int maxCodim, trackingStats *trackCount, int pathMod, double midpoint_tol, tracker_config_t *T, regen_pos_dim_t *RPD, char *startName, char *witName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does regen pos dim track in parallel  - 'headnode'     *
\***************************************************************/
{
  int codim, codim_index, num_paths, num_crossings = 0, num_codim = RPD->num_codim, minPacketSize = 1, maxPacketSize = 20;
  size_t size;
  char *str = NULL;
  FILE *RAWOUT = NULL, *START = NULL, *FAIL = NULL, *OUT = NULL, *WITSUPER = NULL;
  char outName[] = "output", midName[] = "midpath_data", failName[] = "failed_paths", rawTrackFile[] = "rawout_track", rawSortFile[] = "rawout_sort", rawPrepareFile[] = "rawout_prepare", midTrackFile[] = "midout_track", midPrepareFile[] = "midout_prepare";

  // open OUT & FAILE
  OUT = fopen(outName, "w");
  FAIL = fopen(failName, "w");

  // send maxCodim to the workers
  MPI_Bcast(&maxCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  // send RPD to the workers
  bcast_regen_pos_dim_t(RPD, T->MPType, my_id, headnode);

  // send codimension information to the workers
  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { // send the codimension data to the workers
    bcast_regenCodim_t(&RPD->codim[codim_index], RPD->curr_precision, T->MPType, my_id, headnode);
  }

  // send 'startCodimIndex'
  MPI_Bcast(&startCodimIndex, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // loop over the codimensions to find the witness supersets
  for (codim_index = startCodimIndex; codim_index < maxCodim; codim_index++)
  { // track the paths for this codim
    codim = RPD->codim[codim_index].codim;
    num_paths = RPD->codim[codim_index].num_paths;

    // setup START
    size = 1 + snprintf(NULL, 0, "%s_%d", startName, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", startName, codim);
    START = fopen(str, "r");
    // make sure START exits
    if (START == NULL)
    {
      printf("ERROR: The file to contain the start points, '%s', does not exist!\n", str);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // setup RAWOUT
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, codim);
    RAWOUT = fopen(str, "w");

    // read in the number of paths from START
    fscanf(START, "%d", &num_crossings);
    // make sure that we have agreement
    if (num_paths != num_crossings)
    {
      printf("ERROR: The number of paths (%d vs %d) described in '%s' is not correct!\n", num_paths, num_crossings, str);
      bexit(ERROR_INVALID_SIZE);
    }

    // track the paths for this codimension
    head_regen_pos_dim_TrackCodim(minPacketSize, maxPacketSize, trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, my_id, headnode, num_processes);

    // close START & RAWOUT
    fclose(START);
    fclose(RAWOUT);

    // wait until the files are closed
    MPI_Barrier(MPI_COMM_WORLD);

    // check for path crossings for this codimension - do not delete!
    if (RPD->codim[codim_index].useIntrinsicSlice)
      num_crossings = parallel_midpoint_checking(midName, midTrackFile, 0, num_paths, codim, midpoint_tol, my_id, headnode, num_processes);
    else
      num_crossings = parallel_midpoint_checking(midName, midTrackFile, 0, num_paths, RPD->new_variables, midpoint_tol, my_id, headnode, num_processes);

    if (num_crossings > 0)
      printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

    // open WITSUPER for sorting - witness superset points
    size = 1 + snprintf(NULL, 0, "%s_%d", witName, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", witName, codim);
    WITSUPER = fopen(str, "w");

    // setup START to be the file containing the endpoints
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, codim);
    START = fopen(str, "r");
    // make sure START exits
    if (START == NULL)
    {
      printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // reopen RAWOUT for sorting
    size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawSortFile, codim);
    RAWOUT = fopen(str, "w");
    
    // sort the paths for this codimension
    head_regen_pos_dim_SortCodim(minPacketSize, maxPacketSize, trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, WITSUPER, my_id, headnode, num_processes);

    // close START, RAWOUT & WITSUPER
    fclose(START);
    fclose(RAWOUT);
    fclose(WITSUPER);

    // print the regen summary
    RAWOUT = fopen("regenSummary", "w");
    printRPDSummaryData(RPD, codim_index, RAWOUT);
    fclose(RAWOUT);

    // see if we need to prepare the next codim
    if (codim_index + 1 < maxCodim)
    { // setup START to be the file containing the ones that need moved forward
      size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawSortFile, codim);
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
      num_paths = head_regen_pos_dim_PrepareNextCodim(minPacketSize, maxPacketSize, trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, str, my_id, headnode, num_processes); 

      // close RAWOUT & START
      fclose(RAWOUT);
      fclose(START);

      // wait until the files are closed
      MPI_Barrier(MPI_COMM_WORLD);

      // check for path crossings for this codimension - do not delete!
      if (RPD->codim[codim_index].useIntrinsicSlice)
        num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, codim + 1, midpoint_tol, my_id, headnode, num_processes);
      else
        num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, RPD->new_variables, midpoint_tol, my_id, headnode, num_processes);

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

    // delete RAWOUT files
    size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawTrackFile, codim);
    remove(str);

    size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawSortFile, codim);
    remove(str);

    size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawPrepareFile, codim);
    remove(str);
  }

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  // clear memory
  free(str);

  return;
}

void head_regen_pos_dim_TrackCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this codim in parallel     *
\***************************************************************/
{
  int (*change_prec)(void const *, int) = NULL;
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  change_prec = &change_regen_pos_dim_prec;
  create_send_packet = &regen_pos_dim_create_send_packet_track;
  recv_store_packet = &regen_pos_dim_recv_store_packet_track;

  // setup the current level
  RPD->curr_codim = codim_index;

  // display messages
  printf("\nTracking regeneration codim %d of %d: %d path%s to track.\n", RPD->codim[codim_index].codim, num_codim, RPD->codim[codim_index].num_paths, RPD->codim[codim_index].num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking regeneration codim %d.\n", RPD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual parallel tracking
  head_trackPaths2("Tracking", 0, RPD->codim[codim_index].num_paths, minPacketSize, maxPacketSize, trackCount, pathMod, T, RPD, RPD, change_prec, START, OUT, RAWOUT, FAIL, NULL, NULL, NULL, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

void head_regen_pos_dim_SortCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *WITSUPER, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts all of the paths for this codim in parallel      *
\***************************************************************/
{
  int (*change_prec)(void const *, int) = NULL;
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  change_prec = &change_regen_pos_dim_prec;
  create_send_packet = &regen_pos_dim_create_send_packet_sort;
  recv_store_packet = &regen_pos_dim_recv_store_packet_sort;

  // setup the current codim
  RPD->curr_codim = codim_index;

  // initialize the counts
  RPD->codim[codim_index].num_inf = RPD->codim[codim_index].num_superset = RPD->codim[codim_index].num_sing = RPD->codim[codim_index].num_nonsing 
    = RPD->codim[codim_index].num_nonsolns = RPD->codim[codim_index].num_bad = 0;

  // display messages
  printf("\nSorting regeneration codim %d of %d: %d path%s to sort.\n", RPD->codim[codim_index].codim, num_codim, RPD->codim[codim_index].num_paths, RPD->codim[codim_index].num_paths == 1 ? "" : "s");

  // do the actual parallel sorting
  head_trackPaths2("Sorting", 0, RPD->codim[codim_index].num_paths, minPacketSize, maxPacketSize, trackCount, pathMod, T, RPD, RPD, change_prec, START, OUT, RAWOUT, FAIL, WITSUPER, NULL, NULL, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

int head_regen_pos_dim_PrepareNextCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *nextStartFile, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of paths actually tracked - for MIDOUT  *
* NOTES: prepares the next regen codim                          *
\***************************************************************/
{
  int i, next_codim_index = codim_index + 1, trackPts = 0, num_bad = 0, *startPathNum = NULL, *rV = NULL, *badPaths = NULL;
  FILE *NEXTSTARTPTS = NULL, *PREPAREPTS = NULL;
  size_t size = 1 + snprintf(NULL, 0, "%s_temp", nextStartFile);
  char *str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "%s_temp", nextStartFile);

  PREPAREPTS = fopen(str, "w+");
  NEXTSTARTPTS = fopen(nextStartFile, "w+");

  // setup a new gamma & send it to the workers
  if (T->MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(RPD->gamma_d);
    bcast_comp_d(RPD->gamma_d, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    get_comp_rand_mp(RPD->gamma_mp);
    bcast_comp_mp(RPD->gamma_mp, my_id, headnode);
  } 
  else
  { // setup gamma_d, _mp & _rat
    get_comp_rand_rat(RPD->gamma_d, RPD->gamma_mp, RPD->gamma_rat, RPD->curr_precision, T->AMP_max_prec, 0, 0);
    bcast_comp_rat(RPD->gamma_rat, my_id, headnode);
  }

  // create a file containing the next start points & the new degrees
  trackPts = regen_pos_dim_setupNextStart(RPD, codim_index, T, START, NEXTSTARTPTS, PREPAREPTS, &startPathNum);

  // setup rV
  rV = (int *)bmalloc(trackPts * sizeof(int));
  for (i = 0; i < trackPts; i++)
    rV[i] = 0;

  // rewind PREPAREPTS
  rewind(PREPAREPTS);

  // prepare the next points
  head_regen_pos_dim_MoveNextCodim(minPacketSize, maxPacketSize, trackCount, codim_index, num_codim, pathMod, T, RPD, trackPts, PREPAREPTS, OUT, RAWOUT, FAIL, NEXTSTARTPTS, rV, my_id, headnode, num_processes);

  // check for errors
  num_bad = 0;
  for (i = 0; i < trackPts; i++)
  {
    if (rV[i] != 0 && rV[i] != retVal_reached_minTrackT)
    { // store this path to bad_paths
      printf("WARNING: When preparing path %d in codim %d, it had retVal %d.\n", i, RPD->codim[next_codim_index].codim, rV[i]);
      badPaths = (int *)brealloc(badPaths, (num_bad + 1) * sizeof(int));
      badPaths[num_bad] = startPathNum[i];
      num_bad++;
    }
  }

  // close the files
  fclose(PREPAREPTS);
  fclose(NEXTSTARTPTS);

  // delete old file
  remove(str);

  // move the existing file containing the start points
  rename(nextStartFile, "RPDTempFile");
  // open files
  NEXTSTARTPTS = fopen(nextStartFile, "w");
  PREPAREPTS = fopen("RPDTempFile", "r");

  // setup the good ones & remove the bad ones
  if (RPD->codim[next_codim_index].useIntrinsicSlice)
    regen_pos_dim_RemoveBadPaths(NEXTSTARTPTS, PREPAREPTS, &RPD->codim[next_codim_index].num_paths, RPD->codim[next_codim_index].codim, T->MPType, RPD->codim[codim_index].num_nonsolns, trackPts, startPathNum, num_bad, badPaths);
  else
    regen_pos_dim_RemoveBadPaths(NEXTSTARTPTS, PREPAREPTS, &RPD->codim[next_codim_index].num_paths, RPD->new_variables, T->MPType, RPD->codim[codim_index].num_nonsolns, trackPts, startPathNum, num_bad, badPaths);

  // print the relevant data to NEXTSTARTPTS so that it can be used to rerun this exact same problem
  printRPDRelevantData(RPD, T->MPType, next_codim_index, NEXTSTARTPTS);

  // close the files
  fclose(PREPAREPTS);
  fclose(NEXTSTARTPTS);
  // remove temporary file
  remove("RPDTempFile");

  // clear memory 
  free(startPathNum);
  free(rV);
  free(badPaths);
  free(str);

  return trackPts;
}

void head_regen_pos_dim_MoveNextCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, int newPathCount, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int *rV, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prepares the next regen codim                          *
\***************************************************************/
{
  int (*change_prec)(void const *, int) = NULL;
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  change_prec = &change_regen_pos_dim_prec;
  create_send_packet = &regen_pos_dim_create_send_packet_prepare;
  recv_store_packet = &regen_pos_dim_recv_store_packet_prepare;

  // setup the current level
  RPD->curr_codim = codim_index;

  // display messages
  printf("\nPreparing regeneration codim %d of %d: %d path%s to track.\n", RPD->codim[codim_index + 1].codim, num_codim, newPathCount, newPathCount == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Preparing regeneration codim %d.\n", RPD->codim[codim_index + 1].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual parallel tracking
  head_trackPaths2("Preparing", 0, newPathCount, minPacketSize, maxPacketSize, trackCount, pathMod, T, RPD, RPD, change_prec, START, OUT, RAWOUT, FAIL, NEXTSTARTPTS, NULL, rV, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

//////////////////// WORKER FUNCTIONS //////////////////////////

void worker_regen_pos_dim(int my_id, int num_processes, int headnode, int dataType) 
/***************************************************************\
* USAGE:                                                        * 
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: does regen pos dim tracking - 'worker' process         * 
\***************************************************************/
{
  int size, codim, codim_index, num_codim, maxCodim, startCodimIndex = 0;
  trackingStats trackCount;   
  tracker_config_t T;   
  regen_pos_dim_t RPD;
  char *str = NULL, midTrackFile[] = "midout_track", midPrepareFile[] = "midout_prepare";
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // recv the maxCodim
  MPI_Bcast(&maxCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // recv T
  bcast_tracker_config_t(&T, my_id, headnode); 

  // now that we know the precision, set the default to that precision - does not hurt if doing only double precision
  initMP(T.Precision);

  // recv RPD
  bcast_regen_pos_dim_t(&RPD, T.MPType, my_id, headnode);

  // recv codimension information
  num_codim = RPD.num_codim;
  RPD.codim = (regenCodim_t *)bmalloc(num_codim * sizeof(regenCodim_t));
  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { // recv codimension data from the headnode
    bcast_regenCodim_t(&RPD.codim[codim_index], RPD.curr_precision, T.MPType, my_id, headnode);
  }

  // recv 'startCodim'
  MPI_Bcast(&startCodimIndex, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // setup the local files - OUT & FAIL & RAWOUT
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
  RAWOUT = fopen(str, "w+");

  // loop over each codimension
  for (codim_index = startCodimIndex; codim_index < maxCodim; codim_index++)
  { // find the current codim
    codim = RPD.codim[codim_index].codim;

    // setup MIDOUT for this codim
    size = 1 + snprintf(NULL, 0, "%s_%d", midTrackFile, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", midTrackFile, my_id);
    MIDOUT = fopen(str, "w");

    // track this codim
    worker_regen_pos_dim_TrackCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    // close MIDOUT
    fclose(MIDOUT);

    // wait until everybody has closed MIDOUT
    MPI_Barrier(MPI_COMM_WORLD);

    // consider doing midpoint checking in parallel

    // setup MIDOUT for sorting
    MIDOUT = NULL;

    // sort this codim
    worker_regen_pos_dim_SortCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    // see if we need to prepare the next codim
    if (codim_index + 1 < maxCodim)        
    { // reopen MIDOUT for preparing next level
      size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", midPrepareFile, my_id);
      MIDOUT = fopen(str, "w");

      // prepare for the next codim
      worker_regen_pos_dim_PrepareNextCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

      // close MIDOUT
      fclose(MIDOUT);

      // wait until everybody has closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel
    }
  }

  // close OUT, RAWOUT & FAIL
  fclose(OUT);
  fclose(RAWOUT);
  fclose(FAIL);

  // delete RAWOUT
  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  remove(str);

  // delete MIDOUT
  size = 1 + snprintf(NULL, 0, "%s_%d", midTrackFile, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "%s_%d", midTrackFile, my_id);
  remove(str);

  size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "%s_%d", midPrepareFile, my_id);
  remove(str);

  // delete FAIL
  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  remove(str);

  // clear memory
  regen_pos_dim_clear(&RPD, T.MPType);
  tracker_config_clear(&T);
  clearMP();
  free(str);

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

void worker_regen_pos_dim_TrackCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this codim in parallel     *
\***************************************************************/
{
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*useSharpener)(int, int, void const *, void const *);
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

  // setup the evaluators
  ptr_to_eval_d = &regen_pos_dim_eval_d;
  ptr_to_eval_mp = &regen_pos_dim_eval_mp;
  change_prec = &change_regen_pos_dim_prec;
  useSharpener = &worker_useSharpener;
  recv_track_send_packet = &regen_pos_dim_recv_track_send_packet;

  // setup the current codimension
  RPD->curr_codim = codim_index;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking codimension %d.\n", RPD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual tracking
  worker_trackPaths2(trackCount, T, RPD, RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  return;
}

void worker_regen_pos_dim_SortCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts all of the paths for this codim in parallel      *
\***************************************************************/
{
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*useSharpener)(int, int, void const *, void const *);
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

  // setup the evaluators
  ptr_to_eval_d = &regen_pos_dim_eval_d;
  ptr_to_eval_mp = &regen_pos_dim_eval_mp;
  change_prec = &change_regen_pos_dim_prec;
  useSharpener = &worker_useSharpener;
  recv_track_send_packet = &regen_pos_dim_recv_sort_send_packet;

  // setup the current codimension
  RPD->curr_codim = codim_index;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting codimension %d.\n", RPD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual sorting
  worker_trackPaths2(trackCount, T, RPD, RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  return;
}

void worker_regen_pos_dim_PrepareNextCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prepare the start points for the next codim in parallel*
\***************************************************************/
{
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

  // setup the evaluators
  ptr_to_eval_d = &regen_pos_dim_moving_linear_eval_d;
  ptr_to_eval_mp = &regen_pos_dim_moving_linear_eval_mp;
  change_prec = &change_regen_pos_dim_prec;
  useSharpener = &worker_useSharpener;
  recv_track_send_packet = &regen_pos_dim_recv_prepare_send_packet;

  // setup the current codimension
  RPD->curr_codim = codim_index;

  // recv a new gamma
  if (T->MPType == 0)
  { // setup gamma_d
    bcast_comp_d(RPD->gamma_d, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    bcast_comp_mp(RPD->gamma_mp, my_id, headnode);
  }
  else
  { // setup gamma_d, _mp & _rat
    bcast_comp_rat(RPD->gamma_rat, my_id, headnode);
    // setup gamma_d & gamma_mp
    rat_to_d(RPD->gamma_d, RPD->gamma_rat);
    rat_to_mp(RPD->gamma_mp, RPD->gamma_rat);
  }

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Moving linears to regeneration codim %d.\n", RPD->codim[codim_index + 1].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual sorting
  worker_trackPaths2(trackCount, T, RPD, RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  // set the tolerances back
  T->basicNewtonTol = basicTol;
  T->endgameNewtonTol = endgameTol;
  T->final_tolerance = finalTol;

  return;
}

/////////////////// CREATE & SEND, RECV & STORE //////////////////

int regen_pos_dim_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i, j, num_input_vars, codim_index;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (MPType == 0 || MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  // create the packet
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
      {
        fscanf(START, "%lf%lf", &sendPts[i].PD_d.point->coord[j].r, &sendPts[i].PD_d.point->coord[j].i);
        scanRestOfLine(START);
      }
      sendPts[i].last_approx_prec = 52;
      sendPts[i].last_approx_d->size = 0;
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
        scanRestOfLine(START);
      }
      sendPts[i].last_approx_prec = sendPts[i].prec;
      sendPts[i].last_approx_mp->size = 0;
    }
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // clear memory
  RPD = NULL;

  return size;
}

int regen_pos_dim_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, codim_index, recvProc, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (T->MPType == 0 || T->MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else if (T->MPType == 1)
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // recv the last approximations
  corank = (int *)bmalloc(*numRecvPts * sizeof(int));
  sm = (double *)bmalloc(*numRecvPts * sizeof(double));
  lg = (double *)bmalloc(*numRecvPts * sizeof(double));

  send_recv_corank_data(corank, sm, lg, *numRecvPts, recvProc, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  { // print the footer 
    fprintf(OUT, "Path number: %d (ID: %d)\n", (*recvPts)[i].pathNum, recvProc);
    printRPDFooter(RPD, codim_index, &(*recvPts)[i], corank[i], OUT, RAWOUT, FAIL, T, trackCount, sm[i], lg[i]);
  }

  free(sm);
  free(lg);
  free(corank);
  RPD = NULL;

  return recvProc;
}

int regen_pos_dim_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, rankType = 1, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = regen_pos_dim_dehom;

  // recv next packet of start points
  send_recv_endgame_data_t(startPts, numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (*numEndPts != *numStartPts)
  { // clear endPts
    for (i = *numEndPts - 1; i >= 0; i--)
      clear_endgame_data(&(*endPts)[i]);

    // set the number to reallocate
    *numEndPts = *numStartPts;

    // reallocate
    *endPts = (endgame_data_t *)brealloc(*endPts, *numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < *numEndPts; i++)
      init_endgame_data(&(*endPts)[i], T->Precision);
  }

  // make sure that we have paths to track
  if (*numStartPts > 0)
  { // allocate memory
    corank = (int *)bmalloc(*numStartPts * sizeof(int));
    sm = (double *)bmalloc(*numStartPts * sizeof(double));
    lg = (double *)bmalloc(*numStartPts * sizeof(double));

    // track the paths
    for (i = 0; i < *numStartPts; i++)
    { // track the ith path
      if (T->MPType == 0)
      { // track in D
        worker_track_path_rank(rankType, NULL, &corank[i], &sm[i], &lg[i], &(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
      }
      else if (T->MPType == 1)
      { // track in MP
        worker_track_path_rank(rankType, NULL, &corank[i], &sm[i], &lg[i], &(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
      }
      else
      { // track using AMP
        worker_track_path_rank(rankType, NULL, &corank[i], &sm[i], &lg[i], &(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
      }
    }

    // send the packet back
    send_recv_endgame_data_t(endPts, numEndPts, T->MPType, headnode, 1);

    // send the last approximations back
    send_recv_corank_data(corank, sm, lg, *numEndPts, headnode, 1);

    // clear memory
    free(corank);
    free(sm);
    free(lg);
  }

  return 0;
}

int regen_pos_dim_create_send_packet_sort(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen sorting         *
\***************************************************************/
{
  int i, j, num_input_vars, codim_index, *corank = NULL;
  double tempD, *sm = NULL, *lg = NULL;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (MPType == 0 || MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;

  // setup the number of variables that the points have
  if (RPD->codim[codim_index].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index].codim;
  else
    num_input_vars = RPD->new_variables;

  // create the packet
  corank = (int *)bmalloc(size * sizeof(int));
  sm = (double *)bmalloc(size * sizeof(double));
  lg = (double *)bmalloc(size * sizeof(double));

  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    // read in the path number & precision
    fscanf(START, "%d\n%d\n", &sendPts[i].pathNum, &sendPts[i].prec);

    if (sendPts[i].prec < 64)
    { // setup sendPts[i].PD_d
      change_size_point_d(sendPts[i].PD_d.point, num_input_vars);
      sendPts[i].PD_d.point->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(START, "%lf%lf", &sendPts[i].PD_d.point->coord[j].r, &sendPts[i].PD_d.point->coord[j].i);
        scanRestOfLine(START);
      }
      // time
      fscanf(START, "%lf%lf", &sendPts[i].PD_d.time->r, &sendPts[i].PD_d.time->i);
    }
    else
    { // setup sendPts[i].PD_mp
      setprec_point_data_mp(&sendPts[i].PD_mp, sendPts[i].prec);
      change_size_point_mp(sendPts[i].PD_mp.point, num_input_vars);
      sendPts[i].PD_mp.point->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].r, START, 10);
        mpf_inp_str(sendPts[i].PD_mp.point->coord[j].i, START, 10);
        scanRestOfLine(START);
      }
      mpf_inp_str(sendPts[i].PD_mp.time->r, START, 10);
      mpf_inp_str(sendPts[i].PD_mp.time->i, START, 10);
      scanRestOfLine(START);
    }

    // last_approx
    fscanf(START, "%d\n", &sendPts[i].last_approx_prec);
    if (sendPts[i].last_approx_prec < 64)
    { // use _d
      change_size_point_d(sendPts[i].last_approx_d, num_input_vars);
      sendPts[i].last_approx_d->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        fscanf(START, "%lf%lf", &sendPts[i].last_approx_d->coord[j].r, &sendPts[i].last_approx_d->coord[j].i);
        scanRestOfLine(START);
      }
    }
    else
    { // set precision & size correctly
      setprec_point_mp(sendPts[i].last_approx_mp, sendPts[i].last_approx_prec);
      change_size_point_mp(sendPts[i].last_approx_mp, num_input_vars);
      sendPts[i].last_approx_mp->size = num_input_vars;
      for (j = 0; j < num_input_vars; j++)
      {
        mpf_inp_str(sendPts[i].last_approx_mp->coord[j].r, START, 10);
        mpf_inp_str(sendPts[i].last_approx_mp->coord[j].i, START, 10);
        scanRestOfLine(START);
      }
    }

    // other info
    fscanf(START, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n", &tempD, &sendPts[i].condition_number, &tempD, &tempD, &tempD, &tempD, &j);
    // read in corank, retVal, smallest non-zero & largest zero SV
    fscanf(START, "%d\n%d\n%lf\n%lf\n", &corank[i], &sendPts[i].retVal, &sm[i], &lg[i]);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // send the last approximations to 'sendProc'
  send_recv_corank_data(corank, sm, lg, size, sendProc, 1);

  // clear memory
  free(corank);
  free(sm);
  free(lg);
  RPD = NULL;

  return size;
}

int regen_pos_dim_recv_store_packet_sort(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, j, size, type, codim_index, recvProc, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  endpoint_data_d endPt_d;
  endpoint_data_mp endPt_mp;
  endpoint_data_amp endPt_amp;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (T->MPType == 0 || T->MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // recv the last approximations
  corank = (int *)bmalloc(*numRecvPts * sizeof(int));
  sm = (double *)bmalloc(*numRecvPts * sizeof(double));
  lg = (double *)bmalloc(*numRecvPts * sizeof(double));

  send_recv_corank_data(corank, sm, lg, *numRecvPts, recvProc, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  { // see what type it is
    type = (*recvPts)[i].retVal;
    if (type == NONSOLUTION_AND_NONSING)
    { // print to RAWOUT
      fprintf(RAWOUT, "%d\n%d\n", (*recvPts)[i].pathNum, (*recvPts)[i].prec);
      if ((*recvPts)[i].prec < 64)
      { // print using _d
        size = (*recvPts)[i].PD_d.point->size;
        for (j = 0; j < size; j++)
          fprintf(RAWOUT, "%.15e %.15e\n", (*recvPts)[i].PD_d.point->coord[j].r, (*recvPts)[i].PD_d.point->coord[j].i);
      }
      else
      { // print using _mp
        size = (*recvPts)[i].PD_mp.point->size;
        for (j = 0; j < size; j++)
        {
          mpf_out_str(RAWOUT, 10, 0, (*recvPts)[i].PD_mp.point->coord[j].r);
          fprintf(RAWOUT, " ");
          mpf_out_str(RAWOUT, 10, 0, (*recvPts)[i].PD_mp.point->coord[j].i);
          fprintf(RAWOUT, "\n");
        }
      }
    }
    else if (type == SOLUTION_AND_SING || type == SOLUTION_AND_NONSING)
    { // we need to print to WITSUPER
      if (T->MPType == 0) 
      { // setup endPt_d & print it
        init_endpoint_data_d(&endPt_d);

        // setup endPt & last_approx
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates first
          intrinsicToExtrinsic_d(endPt_d.endPt, (*recvPts)[i].PD_d.point, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          intrinsicToExtrinsic_d(endPt_d.last_approx, (*recvPts)[i].last_approx_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
        }  
        else
        { // copy
          point_cp_d(endPt_d.endPt, (*recvPts)[i].PD_d.point);
          point_cp_d(endPt_d.last_approx, (*recvPts)[i].last_approx_d);
        }
        // seutp other data
        set_d(endPt_d.finalT, (*recvPts)[i].PD_d.time);
        endPt_d.cond_num = (*recvPts)[i].condition_number;
        endPt_d.corank = corank[i];
        endPt_d.smallest_nonzero_SV = sm[i];
        endPt_d.largest_zero_SV = lg[i];
        endPt_d.retVal = (*recvPts)[i].retVal;

        // print the data
        print_endpoint_data_d(OTHER, &endPt_d);

        // clear      
        clear_endpoint_data_d(&endPt_d);
      }
      else if (T->MPType == 1)
      { // setup endPt_mp & print it
        init_endpoint_data_mp(&endPt_mp);

        // setup endPt & last_approx
        if (RPD->codim[codim_index].useIntrinsicSlice)
        { // convert to extrinsic coordinates first
          intrinsicToExtrinsic_mp(endPt_mp.endPt, (*recvPts)[i].PD_mp.point, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          intrinsicToExtrinsic_mp(endPt_mp.last_approx, (*recvPts)[i].last_approx_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
        }
        else
        { // copy
          point_cp_mp(endPt_mp.endPt, (*recvPts)[i].PD_mp.point);
          point_cp_mp(endPt_mp.last_approx, (*recvPts)[i].last_approx_mp);
        }
        // seutp other data
        set_mp(endPt_mp.finalT, (*recvPts)[i].PD_mp.time);
        endPt_mp.cond_num = (*recvPts)[i].condition_number;
        endPt_mp.corank = corank[i];
        endPt_mp.smallest_nonzero_SV = sm[i];
        endPt_mp.largest_zero_SV = lg[i];
        endPt_mp.retVal = (*recvPts)[i].retVal;

        // print the data
        print_endpoint_data_mp(OTHER, &endPt_mp, T->Precision);

        // clear
        clear_endpoint_data_mp(&endPt_mp);
      }
      else 
      { // setup endPt_amp & print it
        init_endpoint_data_amp(&endPt_amp, (*recvPts)[i].prec, (*recvPts)[i].last_approx_prec);
        endPt_amp.curr_prec = (*recvPts)[i].prec;
        endPt_amp.last_approx_prec = (*recvPts)[i].last_approx_prec;

        // setup endPt
        if (endPt_amp.curr_prec < 64)
        { // setup endPt_d & finalT_d
          set_d(endPt_amp.finalT_d, (*recvPts)[i].PD_d.time);
          if (RPD->codim[codim_index].useIntrinsicSlice)
          { // convert to extrinsic coordinates first
            intrinsicToExtrinsic_d(endPt_amp.endPt_d, (*recvPts)[i].PD_d.point, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          }
          else
          { // copy
            point_cp_d(endPt_amp.endPt_d, (*recvPts)[i].PD_d.point);
          }
        }
        else
        { // setup endPt_mp & finalT_mp
          set_mp(endPt_amp.finalT_mp, (*recvPts)[i].PD_mp.time);
          if (RPD->codim[codim_index].useIntrinsicSlice)
          { // convert to extrinsic coordinates first
            intrinsicToExtrinsic_mp(endPt_amp.endPt_mp, (*recvPts)[i].PD_mp.point, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          }
          else
          { // copy
            point_cp_mp(endPt_amp.endPt_mp, (*recvPts)[i].PD_mp.point);
          }
        }

        // setup last_approx
        if (endPt_amp.last_approx_prec < 64)
        { // setup last_approx_d
          if (RPD->codim[codim_index].useIntrinsicSlice)
          { // convert to extrinsic coordinates first
            intrinsicToExtrinsic_d(endPt_amp.last_approx_d, (*recvPts)[i].last_approx_d, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);
          }
          else
          { // copy
            point_cp_d(endPt_amp.last_approx_d, (*recvPts)[i].last_approx_d);
          }
        }
        else
        { // setup last_approx_mp
          if (RPD->codim[codim_index].useIntrinsicSlice)
          { // convert to extrinsic coordinates first
            intrinsicToExtrinsic_mp(endPt_amp.last_approx_mp, (*recvPts)[i].last_approx_mp, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);
          }
          else
          { // copy
            point_cp_mp(endPt_amp.last_approx_mp, (*recvPts)[i].last_approx_mp);
          }
        }
        // seutp other data
        endPt_amp.cond_num = (*recvPts)[i].condition_number;
        endPt_amp.corank = corank[i];
        endPt_amp.smallest_nonzero_SV = sm[i];
        endPt_amp.largest_zero_SV = lg[i];
        endPt_amp.retVal = (*recvPts)[i].retVal;

        // print the data
        print_endpoint_data_amp(OTHER, &endPt_amp);

        // clear
        clear_endpoint_data_amp(&endPt_amp);
      }
    }

    // add to the counts
    if (type == retVal_going_to_infinity || type == retVal_security_max)
      RPD->codim[codim_index].num_inf++;
    else if (type == NONSOLUTION_AND_NONSING)
      RPD->codim[codim_index].num_nonsolns++;
    else if (type == SOLUTION_AND_SING)
    {
      RPD->codim[codim_index].num_sing++;
      RPD->codim[codim_index].num_superset++;
    }
    else if (type == SOLUTION_AND_NONSING)
    {
      RPD->codim[codim_index].num_nonsing++;
      RPD->codim[codim_index].num_superset++;
    }
    else
      RPD->codim[codim_index].num_bad++;
  }

  // clear memory
  free(sm);
  free(lg);
  free(corank);
  RPD = NULL;

  return recvProc;
}

int regen_pos_dim_recv_sort_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, finite = 1, soln = 1, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  regen_pos_dim_t *RPD = NULL;

  // setup CD
  if (T->MPType == 0 || T->MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // recv next packet of start points
  send_recv_endgame_data_t(startPts, numStartPts, T->MPType, headnode, 0);

  // recv the last approximations
  corank = (int *)bmalloc(*numStartPts * sizeof(int));
  sm = (double *)bmalloc(*numStartPts * sizeof(double));
  lg = (double *)bmalloc(*numStartPts * sizeof(double));

  send_recv_corank_data(corank, sm, lg, *numStartPts, headnode, 0);

  // setup endPts
  if (*numEndPts != *numStartPts)
  { // clear endPts
    for (i = *numEndPts - 1; i >= 0; i--)
      clear_endgame_data(&(*endPts)[i]);

    // set the number to reallocate
    *numEndPts = *numStartPts;

    // reallocate
    *endPts = (endgame_data_t *)brealloc(*endPts, *numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < *numEndPts; i++)
      init_endgame_data(&(*endPts)[i], T->Precision);
  }

  // make sure that we have paths to sort
  if (*numStartPts > 0)
  { // loop over each point
    for (i = 0; i < *numStartPts; i++)
    { // initialize rankDef & endPts
      (*endPts)[i].prec = (*startPts)[i].prec;
      (*endPts)[i].pathNum = (*startPts)[i].pathNum;
      (*endPts)[i].condition_number = (*startPts)[i].condition_number;
      if ((*endPts)[i].prec < 64)
      {
        point_data_cp_d(&(*endPts)[i].PD_d, &(*startPts)[i].PD_d);
      }
      else
      {
        setprec_point_data_mp(&(*endPts)[i].PD_mp, (*endPts)[i].prec);
        point_data_cp_mp(&(*endPts)[i].PD_mp, &(*startPts)[i].PD_mp);
      }

      (*endPts)[i].last_approx_prec = (*startPts)[i].last_approx_prec;
      if ((*endPts)[i].last_approx_prec < 64)
      {
        point_cp_d((*endPts)[i].last_approx_d, (*startPts)[i].last_approx_d);
      }
      else
      {
        setprec_point_mp((*endPts)[i].last_approx_mp, (*endPts)[i].last_approx_prec);
        point_cp_mp((*endPts)[i].last_approx_mp, (*startPts)[i].last_approx_mp);
      }

      // do the sorting
      if ((*startPts)[i].retVal == 0)
      { // classify
        regen_pos_dim_SortEndpoint(&corank[i], &finite, &soln, (*endPts)[i].condition_number, RPD, RPD->curr_codim, (*endPts)[i].pathNum, T, OUT, (*endPts)[i].PD_d.point, (*endPts)[i].PD_mp.point, (*endPts)[i].PD_d.time, (*endPts)[i].PD_mp.time, (*endPts)[i].prec, (*endPts)[i].last_approx_d, (*endPts)[i].last_approx_mp, (*endPts)[i].last_approx_prec, eval_d, eval_mp, change_prec);

        // classify the end point
        if (!finite)
        { // dehom point is infinite
          (*endPts)[i].retVal = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (corank[i] > 0)
          { // classify as singular
            if (soln)
              (*endPts)[i].retVal = SOLUTION_AND_SING;
            else
              (*endPts)[i].retVal = NONSOLUTION_AND_SING;
          }
          else
          { // classify as non-singular
            if (soln)
              (*endPts)[i].retVal = SOLUTION_AND_NONSING;
            else
              (*endPts)[i].retVal = NONSOLUTION_AND_NONSING;
          }
        }
      }
      else
      { // path was not a success - copy over error code
        (*endPts)[i].retVal = (*startPts)[i].retVal;
      }
    }

    // send the packet back
    send_recv_endgame_data_t(endPts, numEndPts, T->MPType, headnode, 1);

    // send the last approximations back
    send_recv_corank_data(corank, sm, lg, *numEndPts, headnode, 1);
  }

  // clear
  free(corank);
  free(sm);
  free(lg);

  RPD = NULL;

  return 0;
}

int regen_pos_dim_recv_prepare_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = zero_dehom;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (T->MPType == 0 || T->MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // recv next packet of start points
  send_recv_endgame_data_t(startPts, numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (*numEndPts != *numStartPts)
  { // clear endPts
    for (i = *numEndPts - 1; i >= 0; i--)
      clear_endgame_data(&(*endPts)[i]);

    // set the number to reallocate
    *numEndPts = *numStartPts;

    // reallocate
    *endPts = (endgame_data_t *)brealloc(*endPts, *numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < *numEndPts; i++)
      init_endgame_data(&(*endPts)[i], T->Precision);
  }

  // make sure that we have paths to track
  if (*numStartPts > 0)
  {// track the paths
    for (i = 0; i < *numStartPts; i++)
    { // setup moving_degree
      RPD->moving_degree = (*startPts)[i].codim;

      // track the ith path
      worker_track_path(&(*startPts)[i], &(*endPts)[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener, find_dehom);
    }

    // send the packet back
    send_recv_endgame_data_t(endPts, numEndPts, T->MPType, headnode, 1);
  }

  RPD = NULL;

  return 0;
}

int regen_pos_dim_create_send_packet_prepare(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i, j, num_input_vars, codim_index;
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (MPType == 0 || MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;

  // setup the number of variables that the start points have
  if (RPD->codim[codim_index + 1].useIntrinsicSlice)
    num_input_vars = RPD->codim[codim_index + 1].codim;
  else
    num_input_vars = RPD->new_variables;

  // create the packet
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
      {
        fscanf(START, "%lf%lf", &sendPts[i].PD_d.point->coord[j].r, &sendPts[i].PD_d.point->coord[j].i);
        scanRestOfLine(START);
      }

      sendPts[i].last_approx_prec = 52;
      sendPts[i].last_approx_d->size = 0;
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
        scanRestOfLine(START);
      }

      sendPts[i].last_approx_prec = sendPts[i].prec;
      sendPts[i].last_approx_mp->size = 0;
    }
    // scan in the next degree (store in codim)
    fscanf(START, "%d", &sendPts[i].codim);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // clear memory
  RPD = NULL;

  return size;
}

int regen_pos_dim_recv_store_packet_prepare(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, j, codim_index, next_codim_index, recvProc, corank = 0; // moving full rank linears
  regen_pos_dim_t *RPD = NULL;

  // setup RPD
  if (T->MPType == 0 || T->MPType == 2)
    RPD = (regen_pos_dim_t *)ED_d;
  else if (T->MPType == 1)
    RPD = (regen_pos_dim_t *)ED_mp;

  // setup the current codim
  codim_index = RPD->curr_codim;
  next_codim_index = codim_index + 1;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  { // print the footer
    fprintf(OUT, "Path number: %d (ID: %d)\n", (*recvPts)[i].pathNum, recvProc);
    printRPDFooter(RPD, next_codim_index, &(*recvPts)[i], corank, OUT, RAWOUT, FAIL, T, trackCount, 0, 0);
    // save the retVal
    rV[(*recvPts)[i].pathNum] = (*recvPts)[i].retVal;
    // print to OTHER
    fprintf(OTHER, "%d\n", (*recvPts)[i].pathNum);
    if ((*recvPts)[i].prec < 64)
    {
      for (j = 0; j < (*recvPts)[i].PD_d.point->size; j++)
        fprintf(OTHER, "%.15e %.15e\n", (*recvPts)[i].PD_d.point->coord[j].r, (*recvPts)[i].PD_d.point->coord[j].i);
    }
    else
    {
      for (j = 0; j < (*recvPts)[i].PD_mp.point->size; j++)
      {
        mpf_out_str(OTHER, 10, 0, (*recvPts)[i].PD_mp.point->coord[j].r);
        fprintf(OTHER, " ");
        mpf_out_str(OTHER, 10, 0, (*recvPts)[i].PD_mp.point->coord[j].i);
        fprintf(OTHER, "\n");
      }
    }
  }

  RPD = NULL;

  return recvProc;
}


#endif


