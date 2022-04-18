// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "eqbyeq.h"
#include "parallel.h"

#ifdef _HAVE_MPI

void head_eqbyeqWitnessTrack_d(trackingStats *trackCount, int curr_sub, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes);

void head_eqbyeqSortWitnessEndpoints_d(basic_eval_data_d *BED, tracker_config_t *T, int curr_sub, double final_tol, int pathMod, int my_id, int headnode, int num_processes);

void head_eqbyeqStageTrack_d(trackingStats *trackCount, int stage, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes);

void head_eqbyeqSortStageEndpoints_d(basic_eval_data_d *BED, tracker_config_t *T, int stage, double final_tol, int pathMod, int my_id, int headnode, int num_processes);

void worker_eqbyeqWitnessTrack_d(trackingStats *trackCount, int curr_sub, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

void worker_eqbyeqSortWitnessEndpoints_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, int curr_sub, int my_id, int headnode, int num_processes);

void worker_eqbyeqStageTrack_d(trackingStats *trackCount, int curr_stage, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

void worker_eqbyeqSortStageEndpoints_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, int curr_stage, int my_id, int headnode, int num_processes);

void head_eqbyeqWitnessTrack_mp(trackingStats *trackCount, int curr_sub, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes);

void head_eqbyeqSortWitnessEndpoints_mp(basic_eval_data_mp *BED, tracker_config_t *T, int curr_sub, double final_tol, int pathMod, int my_id, int headnode, int num_processes);

void head_eqbyeqStageTrack_mp(trackingStats *trackCount, int stage, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes);

void head_eqbyeqSortStageEndpoints_mp(basic_eval_data_mp *BED, tracker_config_t *T, int stage, double final_tol, int pathMod, int my_id, int headnode, int num_processes);

void worker_eqbyeqWitnessTrack_mp(trackingStats *trackCount, int curr_sub, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

void worker_eqbyeqSortWitnessEndpoints_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, int curr_sub, int my_id, int headnode, int num_processes);

void worker_eqbyeqStageTrack_mp(trackingStats *trackCount, int curr_stage, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

void worker_eqbyeqSortStageEndpoints_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, int curr_stage, int my_id, int headnode, int num_processes);

void head_eqbyeq_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does eq-by-eq tracking - 'headnode'                    *
\***************************************************************/
{
  int stage, num_crossings = 0, size, num_paths, num_vars, num_subs = ED_d->EqD->num_subsystems;
  char *str = NULL, nonName[] = "nonsolutions";
  FILE *NONSOLN = NULL;

  if (!ED_d->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen(nonName, "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // initialize trackCount
  init_trackingStats(trackCount);

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  if (T->MPType == 2)
  { // send ED_d to the workers - using AMP
    bcast_basic_eval_data_amp(ED_d, my_id, headnode);
  }
  else
  { // send ED_d to the workers - only using double precision
    bcast_basic_eval_data_d(ED_d, T->MPType, my_id, headnode);
  }

  // send eqData_t to the workers
  bcast_eqData_t(ED_d->EqD, T->MPType, my_id, headnode);

  // send witness information to the workers
  for (stage = 0; stage < num_subs; stage++)
  { // send the 'stage'th witness data to the workers
    bcast_witnessData(ED_d->EqD, stage, T->MPType, my_id, headnode);
  }

  // loop to generate the witness sets for each of the subsystems
  for (stage = 0; stage < num_subs; stage++)
  { // find the number of paths and the number of variables used in this subsystem
    num_paths = ED_d->EqD->witnessData_d[stage].num_paths;
    num_vars = ED_d->EqD->witnessData_d[stage].depth;

    // find the witness points for this subsystem
    head_eqbyeqWitnessTrack_d(trackCount, stage, pathMod, T, ED_d, ED_mp, OUT, RAWOUT, FAIL, NONSOLN, my_id, headnode, num_processes);

    // check for path crossings if this is not the only subsystem
    if (num_subs > 1)
    { // wait unitl the files are close
      MPI_Barrier(MPI_COMM_WORLD);

      // setup str
      size = 1 + snprintf(NULL, 0, "midout_w%d", stage);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_w%d", stage);

      num_crossings = parallel_midpoint_checking(midFile, str, 1, num_paths, num_vars, midpoint_tol, my_id, headnode, num_processes);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // sort the endpoints
    head_eqbyeqSortWitnessEndpoints_d(ED_d, T, stage, target_tol, pathMod, my_id, headnode, num_processes);
  }

  // setup the first stage
  setupEqbyEqFirstStage_d(ED_d->EqD, T->MPType);

  // clear the first witness data
  clearEqbyEqFirstWitnessData_d(ED_d->EqD, T->MPType);

  // loop to track each of the stages
  for (stage = 1; stage < num_subs; stage++)
  { // setup the next stage
    setupEqbyEqNextStage_d(ED_d, stage, T->MPType, T->AMP_max_prec);

    // send the stage information to the workers
    bcast_stageData(ED_d->EqD, stage, T->MPType, my_id, headnode);

    // clear the stage data from stage 'stage - 1' since it is no longer needed
    clearEqbyEqStageData_d(ED_d->EqD, stage - 1, T->MPType);

    // clear the witness data from subsystem 'stage' since it is no longer needed
    clearEqbyEqWitnessData_d(ED_d->EqD, stage, T->MPType);

    // find the number of paths for this next stage
    num_paths = ED_d->EqD->stageData_d[stage].num_paths;
    if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
      num_vars = ED_d->EqD->stageData_d[stage].depth_x + ED_d->EqD->stageData_d[stage].depth_y;
    else
      num_vars = 2 * ED_d->EqD->num_vars;

    // track each path for this stage
    head_eqbyeqStageTrack_d(trackCount, stage, pathMod, T, ED_d, ED_mp, OUT, RAWOUT, FAIL, NONSOLN, my_id, headnode, num_processes);

    // check for path crossings if this is not the last subsystem
    if (stage + 1 < num_subs)
    { // wait unitl the files are close
      MPI_Barrier(MPI_COMM_WORLD);

      // setup str
      size = 1 + snprintf(NULL, 0, "midout_s%d", stage);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_s%d", stage);

      num_crossings = parallel_midpoint_checking(midFile, str, 1, num_paths, num_vars, midpoint_tol, my_id, headnode, num_processes);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // sort the endpoints
    head_eqbyeqSortStageEndpoints_d(ED_d, T, stage, target_tol, pathMod, my_id, headnode, num_processes);
  }

  // set the number of points that were tracked on the last level
  trackCount->numPoints = ED_d->EqD->stageData_d[num_subs - 1].num_paths;

  // wait until the files are closed
  MPI_Barrier(MPI_COMM_WORLD);

  // setup str
  if (num_subs == 1)
    size = 1 + snprintf(NULL, 0, "midout_w%d", num_subs - 1);
  else
    size = 1 + snprintf(NULL, 0, "midout_s%d", num_subs - 1);
  str = (char *)brealloc(str, size * sizeof(char));
  if (num_subs == 1)
    sprintf(str, "midout_w%d", num_subs - 1);
  else
    sprintf(str, "midout_s%d", num_subs - 1);

  // combine the midpath_data files
  combine_midpath_data(midFile, str, 1, headnode, num_processes);

  // delete the rawout & fail worker files
  delete_parallel_files("rawout", headnode, num_processes);
  delete_parallel_files("fail", headnode, num_processes);

  if (!ED_d->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // clear memory
  free(str);

  return;
}

void head_eqbyeqWitnessTrack_d(trackingStats *trackCount, int curr_sub, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the witness points for the current subsystem     *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize;
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  int path_num, num_subs = ED_d->EqD->num_subsystems, num_paths = ED_d->EqD->witnessData_d[curr_sub].num_paths;
  int (*change_prec)(void const *, int);

  point_d dehom_d, orig_last_d;
  point_mp dehom_mp, orig_last_mp, tempPoint_mp;
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  init_point_d(dehom_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last_mp, 0);
  init_point_mp(tempPoint_mp, 0);

  // setup the subsystem
  ED_d->EqD->curr_stage_num = curr_sub;

  // setup the evaluators for printing the footer for witness tracking
  change_prec = &change_eqbyeq_eval_prec;

  // print header for level
  printf("\nFinding witness points for subsystem %d of %d: %d path%s to track.\n", curr_sub, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Finding witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum - randomize the list of start points
  pathNum = (int *)brealloc(pathNum, num_paths * sizeof(int));
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to i
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
    init_endgame_data(&sendPts[i], 64);

  // initialize count
  count = 0;

  // send out the initial set of packets for this level
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = 52;  // start point in double precision
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num];
        point_cp_d(sendPts[j].PD_d.point, ED_d->EqD->witnessData_d[curr_sub].startPts[path_num]);
        set_double_d(sendPts[j].PD_d.time, 1, 0);
        sendPts[j].last_approx_prec = 52;
        sendPts[j].last_approx_d->size = 0;
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Tracking path %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = 52;  // start point in double precision
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num];
      point_cp_d(sendPts[j].PD_d.point, ED_d->EqD->witnessData_d[curr_sub].startPts[path_num]);
      set_double_d(sendPts[j].PD_d.time, 1, 0);
      sendPts[j].last_approx_prec = 52;
      sendPts[j].last_approx_d->size = 0;
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED_d->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED_d->EqD->witnessData_d[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED_d->EqD->witnessData_d[curr_sub].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, ED_d, ED_mp, curr_sub, 0);

      // setup orig_last
      if (recvPts[i].last_approx_prec < 64)
      { // use _d
        intrinsicToExtrinsic_d(orig_last_d, recvPts[i].last_approx_d, ED_d->EqD->witnessData_d[curr_sub].B, ED_d->EqD->witnessData_d[curr_sub].p);

        if (recvPts[i].prec > 52)
        { // convert to _mp
          setprec_point_mp(orig_last_mp, recvPts[i].prec);
          point_d_to_mp(orig_last_mp, orig_last_d);
        }
      }
      else
      { // use _mp
        setprec_point_mp(orig_last_mp, recvPts[i].last_approx_prec);

        intrinsicToExtrinsic_mp(orig_last_mp, recvPts[i].last_approx_mp, ED_mp->EqD->witnessData_mp[curr_sub].B, ED_mp->EqD->witnessData_mp[curr_sub].p);

        if (recvPts[i].prec < 64)
        { // convert _d
          point_mp_to_d(orig_last_d, orig_last_mp);
        }
      }

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        point_cp_d(ED_d->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_d.point);
        set_d(ED_d->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_d.time);

        // convert to the original coordinates
        intrinsicToExtrinsic_d(ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], recvPts[i].PD_d.point, ED_d->EqD->witnessData_d[curr_sub].B, ED_d->EqD->witnessData_d[curr_sub].p);

        // find dehom_d
        getDehomPoint_d(dehom_d, ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], ED_d->EqD->witnessData_d[curr_sub].endPts[path_num]->size, &ED_d->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num] = printWitnessFooter_d(ED_d, curr_sub, path_num, &recvPts[i].PD_d, ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], orig_last_d, dehom_d, recvPts[i].condition_number, recvPts[i].function_residual_d, recvPts[i].latest_newton_residual_d, recvPts[i].t_val_at_latest_sample_point_d, recvPts[i].error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
      else
      { // set the precision
        T->Precision = recvPts[i].prec;
        initMP(recvPts[i].prec);
        change_prec(ED_mp, recvPts[i].prec);

        // setup dehom_mp & tempPoint_mp
        setprec_point_mp(dehom_mp, recvPts[i].prec);
        setprec_point_mp(tempPoint_mp, recvPts[i].prec);

        // copy over to the appropriate spot - converting to double precision
        point_mp_to_d(ED_d->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
        mp_to_d(ED_d->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

        // covert to the original coordinates
        intrinsicToExtrinsic_mp(tempPoint_mp, recvPts[i].PD_mp.point, ED_mp->EqD->witnessData_mp[curr_sub].B, ED_mp->EqD->witnessData_mp[curr_sub].p);

        // copy tempPoint_mp to endPts
        point_mp_to_d(ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], tempPoint_mp);

        // find dehom_mp
        getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num] = printWitnessFooter_mp(ED_mp, curr_sub, path_num, &recvPts[i].PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED_d->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED_d->EqD->witnessData_d[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED_d->EqD->witnessData_d[curr_sub].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, ED_d, ED_mp, curr_sub, 0);

      // setup orig_last
      if (recvPts[i].last_approx_prec < 64)
      { // use _d
        intrinsicToExtrinsic_d(orig_last_d, recvPts[i].last_approx_d, ED_d->EqD->witnessData_d[curr_sub].B, ED_d->EqD->witnessData_d[curr_sub].p);

        if (recvPts[i].prec > 52)
        { // convert to _mp
          setprec_point_mp(orig_last_mp, recvPts[i].prec);
          point_d_to_mp(orig_last_mp, orig_last_d);
        }
      }
      else
      { // use _mp
        setprec_point_mp(orig_last_mp, recvPts[i].last_approx_prec);

        intrinsicToExtrinsic_mp(orig_last_mp, recvPts[i].last_approx_mp, ED_mp->EqD->witnessData_mp[curr_sub].B, ED_mp->EqD->witnessData_mp[curr_sub].p);

        if (recvPts[i].prec < 64)
        { // convert _d
          point_mp_to_d(orig_last_d, orig_last_mp);
        }
      }

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        point_cp_d(ED_d->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_d.point);
        set_d(ED_d->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_d.time);

        // convert to the original coordinates
        intrinsicToExtrinsic_d(ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], recvPts[i].PD_d.point, ED_d->EqD->witnessData_d[curr_sub].B, ED_d->EqD->witnessData_d[curr_sub].p);

        // find dehom_d
        getDehomPoint_d(dehom_d, ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], ED_d->EqD->witnessData_d[curr_sub].endPts[path_num]->size, &ED_d->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num] = printWitnessFooter_d(ED_d, curr_sub, path_num, &recvPts[i].PD_d, ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], orig_last_d, dehom_d, recvPts[i].condition_number, recvPts[i].function_residual_d, recvPts[i].latest_newton_residual_d, recvPts[i].t_val_at_latest_sample_point_d, recvPts[i].error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
      else
      { // set the precision
        T->Precision = recvPts[i].prec;
        initMP(recvPts[i].prec);
        change_prec(ED_mp, recvPts[i].prec);

        // setup dehom_mp & tempPoint_mp
        setprec_point_mp(dehom_mp, recvPts[i].prec);
        setprec_point_mp(tempPoint_mp, recvPts[i].prec);

        // copy over to the appropriate spot - converting to double precision
        point_mp_to_d(ED_d->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
        mp_to_d(ED_d->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

        // covert to the original coordinates
        intrinsicToExtrinsic_mp(tempPoint_mp, recvPts[i].PD_mp.point, ED_mp->EqD->witnessData_mp[curr_sub].B, ED_mp->EqD->witnessData_mp[curr_sub].p);

        // copy tempPoint_mp to endPts
        point_mp_to_d(ED_d->EqD->witnessData_d[curr_sub].endPts[path_num], tempPoint_mp);

        // find dehom_mp
        getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->witnessData_d[curr_sub].endPt_retVals[path_num] = printWitnessFooter_mp(ED_mp, curr_sub, path_num, &recvPts[i].PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  free(recvPts);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);

  clear_point_d(dehom_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last_mp);
  clear_point_mp(tempPoint_mp);

  return;
}

void head_eqbyeqSortWitnessEndpoints_d(basic_eval_data_d *BED, tracker_config_t *T, int curr_sub, double final_tol, int pathMod, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at for the current subsystem *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, indexI = 0, indexJ = 0, count, path_num, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize, num_paths = BED->EqD->witnessData_d[curr_sub].num_paths;
  int  num_subs = BED->EqD->num_subsystems, num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int)), *higherDim = NULL;
  vec_d tempVec;
  sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  init_vec_d(tempVec, 0);
  
  // print header for level
  printf("\nSorting witness points for subsystem %d of %d: %d path%s to sort.\n", curr_sub, num_subs, num_paths, num_paths == 1 ? "" : "s");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // randomize the list
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to j
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
    init_endgame_data(&sendPts[i], 64);
  higherDim = (int *)bmalloc(maxSize * sizeof(int));
  for (i = 0; i < maxSize; i++)
    higherDim[i] = 0;
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Sorting %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = 52; // point in double precision
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = BED->EqD->witnessData_d[curr_sub].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
        point_cp_d(sendPts[j].PD_d.point, BED->EqD->witnessData_d[curr_sub].endPts_in[path_num]);
        set_d(sendPts[j].PD_d.time, BED->EqD->witnessData_d[curr_sub].finalTs[path_num]);
        sendPts[j].last_approx_prec = 52;
        sendPts[j].last_approx_d->size = 0;
        higherDim[j] = BED->EqD->witnessData_d[curr_sub].higherDim[path_num];
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      // send higherDim
      MPI_Send(higherDim, lastSize[i], MPI_INT, i, TAG_INT, MPI_COMM_WORLD);      

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be sorted
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Sorting %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = 52;  // start point in double precision
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = BED->EqD->witnessData_d[curr_sub].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
      point_cp_d(sendPts[j].PD_d.point, BED->EqD->witnessData_d[curr_sub].endPts_in[path_num]);
      set_d(sendPts[j].PD_d.time, BED->EqD->witnessData_d[curr_sub].finalTs[path_num]);
      sendPts[j].last_approx_prec = 52;
      sendPts[j].last_approx_d->size = 0;
      higherDim[j] = BED->EqD->witnessData_d[curr_sub].higherDim[path_num];
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    // send higherDim
    MPI_Send(higherDim, lastSize[recvProc], MPI_INT, recvProc, TAG_INT, MPI_COMM_WORLD);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->witnessData_d[curr_sub].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->witnessData_d[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        point_cp_d(BED->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_d.point);
        set_d(BED->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_d.time);
      }
      else
      { // copy over to the appropriate spot - converting to double precision
        point_mp_to_d(BED->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
        mp_to_d(BED->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);
      }

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      sortPts[path_num].norm = infNormVec_d(BED->EqD->witnessData_d[curr_sub].endPts[path_num]);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->witnessData_d[curr_sub].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->witnessData_d[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        point_cp_d(BED->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_d.point);
        set_d(BED->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_d.time);
      }
      else
      { // copy over to the appropriate spot - converting to double precision
        point_mp_to_d(BED->EqD->witnessData_d[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
        mp_to_d(BED->EqD->witnessData_d[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);
      }

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      sortPts[path_num].norm = infNormVec_d(BED->EqD->witnessData_d[curr_sub].endPts[path_num]);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // sort sortPts
  qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

  // do the final classification
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      j = i + 1;
      while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
      {
        indexJ = sortPts[j].path_num;
        if (BED->EqD->witnessData_d[curr_sub].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          vec_sub_d(tempVec, BED->EqD->witnessData_d[curr_sub].endPts[indexI], BED->EqD->witnessData_d[curr_sub].endPts[indexJ]);
          if (infNormVec_d(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] = BED->EqD->witnessData_d[curr_sub].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
        j++;
      }
    }

    // add to count
    if (BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == retVal_going_to_infinity || BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (BED->EqD->witnessData_d[curr_sub].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  BED->EqD->witnessData_d[curr_sub].num_sing = num_sing;
  BED->EqD->witnessData_d[curr_sub].num_nonsing = num_nonsing;
  BED->EqD->witnessData_d[curr_sub].num_higher_dim = num_higher_dim;
  BED->EqD->witnessData_d[curr_sub].num_bad = num_bad;
  BED->EqD->witnessData_d[curr_sub].num_inf = num_inf;

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  free(recvPts);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  free(sendPts); // only sent in double precision
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  free(sortPts);
  free(higherDim);

  clear_vec_d(tempVec);

  return;
}

void head_eqbyeqStageTrack_d(trackingStats *trackCount, int stage, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the current stage of equation-by-equation       *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize;
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  int path_num, num_stages = ED_d->EqD->num_subsystems, num_paths = ED_d->EqD->stageData_d[stage].num_paths;
  int (*change_prec)(void const *, int);

  point_d dehom_d, orig_last_d;
  point_mp dehom_mp, orig_last_mp, tempPoint_mp;
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  init_point_d(dehom_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last_mp, 0);
  init_point_mp(tempPoint_mp, 0);

  // setup the function for changing precision
  change_prec = &change_eqbyeq_eval_prec;

  // setup the subsystem
  ED_d->EqD->curr_stage_num = stage;

  // print header for level
  printf("\nTracking points for stage %d of %d: %d path%s to track.\n", stage, num_stages, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking points for stage %d.\n", stage);
  fprintf(OUT, "*****************************************************\n");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum - randomize the list of start points
  pathNum = (int *)brealloc(pathNum, num_paths * sizeof(int));
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to i
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
    init_endgame_data(&sendPts[i], 64);

  // initialize count
  count = 0;

  // send out the initial set of packets for this level
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = 52;  // start point in double precision
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = ED_d->EqD->stageData_d[stage].endPt_retVals[path_num];
        point_cp_d(sendPts[j].PD_d.point, ED_d->EqD->stageData_d[stage].startPts[path_num]);
        set_double_d(sendPts[j].PD_d.time, 1, 0);
        sendPts[j].last_approx_prec = 52;
        sendPts[j].last_approx_d->size = 0;
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Tracking path %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = 52;  // start point in double precision
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = ED_d->EqD->stageData_d[stage].endPt_retVals[path_num];
      point_cp_d(sendPts[j].PD_d.point, ED_d->EqD->stageData_d[stage].startPts[path_num]);
      set_double_d(sendPts[j].PD_d.time, 1, 0);
      sendPts[j].last_approx_prec = 52;
      sendPts[j].last_approx_d->size = 0;
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED_d->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED_d->EqD->stageData_d[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED_d->EqD->stageData_d[stage].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, ED_d, ED_mp, stage, 1);

      // find orig_last
      if (recvPts[i].last_approx_prec < 64)
      {
        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_d(orig_last_d, recvPts[i].last_approx_d, ED_d->EqD->stageData_d[stage].B0, ED_d->EqD->stageData_d[stage].p0);
        }
        else
        { // copy over
          point_cp_d(orig_last_d, recvPts[i].last_approx_d);
        }
        orig_last_d->size /= 2; // remove the bottom half of the coordinates

        if (recvPts[i].prec > 52)
        { // convert to _mp
          setprec_point_mp(orig_last_mp, recvPts[i].prec);
          point_d_to_mp(orig_last_mp, orig_last_d);
        }
      }
      else
      {
        setprec_point_mp(orig_last_mp, recvPts[i].last_approx_prec);
        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_mp(orig_last_mp, recvPts[i].last_approx_mp, ED_d->EqD->stageData_mp[stage].B0, ED_d->EqD->stageData_mp[stage].p0);
        }
        else
        { // copy over
          point_cp_mp(orig_last_mp, recvPts[i].last_approx_mp);
        }
        orig_last_d->size /= 2; // remove the bottom half of the coordinates

        if (recvPts[i].prec < 64)
        { // convert to _d
          point_mp_to_d(orig_last_d, orig_last_mp);
        }
      }

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        set_d(ED_d->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_d.time);
        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store the intrinsic endpoint and its corresponding point in the original variables
          intrinsicToExtrinsic_d(dehom_d, recvPts[i].PD_d.point, ED_d->EqD->stageData_d[stage].B0, ED_d->EqD->stageData_d[stage].p0);
          dehom_d->size /= 2; // remove the bottom half of the coordinates

          point_cp_d(ED_d->EqD->stageData_d[stage].endPts[path_num], dehom_d);
 
          // convert to intrinsic coordinates
          mat_d B_transpose;
          init_mat_d(B_transpose, ED_d->EqD->stageData_d[stage].B->cols, ED_d->EqD->stageData_d[stage].B->rows);
          transpose_d(B_transpose, ED_d->EqD->stageData_d[stage].B);
 
          extrinsicToIntrinsic_d(ED_d->EqD->stageData_d[stage].endPts_in[path_num], dehom_d, B_transpose, ED_d->EqD->stageData_d[stage].p);

          clear_mat_d(B_transpose);
        }
        else
        { // use the top set of coordinates so that we can store the original coordinates
          point_cp_d(ED_d->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_d.point);
          ED_d->EqD->stageData_d[stage].endPts[path_num]->size /= 2;
        }

        // find dehom_d
        getDehomPoint_d(dehom_d, ED_d->EqD->stageData_d[stage].endPts[path_num], ED_d->EqD->stageData_d[stage].endPts[path_num]->size, &ED_d->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->stageData_d[stage].endPt_retVals[path_num] = printStageFooter_d(ED_d, stage, path_num, &recvPts[i].PD_d, ED_d->EqD->stageData_d[stage].endPts[path_num], orig_last_d, dehom_d, recvPts[i].condition_number, recvPts[i].function_residual_d, recvPts[i].latest_newton_residual_d, recvPts[i].t_val_at_latest_sample_point_d, recvPts[i].error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
      else
      { // set the precision
        T->Precision = recvPts[i].prec;
        initMP(recvPts[i].prec);
        change_prec(ED_mp, recvPts[i].prec);

        // change dehom_mp & tempPoint_mp
        setprec_point_mp(dehom_mp, recvPts[i].prec);
        setprec_point_mp(tempPoint_mp, recvPts[i].prec);

        // copy over to the appropriate spot - converting to double precision
        mp_to_d(ED_d->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_mp.time);

        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store the intrinsic endpoint and its corresponding point in the original variables
          intrinsicToExtrinsic_mp(tempPoint_mp, recvPts[i].PD_mp.point, ED_d->EqD->stageData_mp[stage].B0, ED_d->EqD->stageData_mp[stage].p0);
          tempPoint_mp->size /= 2; // remove the bottom half of the coordinates

          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts[path_num], tempPoint_mp);

          // convert to the original coordinates
          mat_mp B_transpose;
          init_mat_mp2(B_transpose, ED_d->EqD->stageData_mp[stage].B->cols, ED_d->EqD->stageData_mp[stage].B->rows, recvPts[i].prec);
          transpose_mp(B_transpose, ED_d->EqD->stageData_mp[stage].B);

          extrinsicToIntrinsic_mp(dehom_mp, tempPoint_mp, B_transpose, ED_d->EqD->stageData_mp[stage].p);

          // copy to endPts_in
          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts_in[path_num], dehom_mp);

          clear_mat_mp(B_transpose);
        }
        else
        { // use the top set of coordinates so that we can store the original coordinates
          point_cp_mp(tempPoint_mp, recvPts[i].PD_mp.point);
          tempPoint_mp->size /= 2;

          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts[path_num], tempPoint_mp);
        }

        // find dehom_mp
        getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->stageData_d[stage].endPt_retVals[path_num] = printStageFooter_mp(ED_mp, stage, path_num, &recvPts[i].PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED_d->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED_d->EqD->stageData_d[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED_d->EqD->stageData_d[stage].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, ED_d, ED_mp, stage, 1);

      // find orig_last
      if (recvPts[i].last_approx_prec < 64)
      {
        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_d(orig_last_d, recvPts[i].last_approx_d, ED_d->EqD->stageData_d[stage].B0, ED_d->EqD->stageData_d[stage].p0);
        }
        else
        { // copy over
          point_cp_d(orig_last_d, recvPts[i].last_approx_d);
        }
        orig_last_d->size /= 2; // remove the bottom half of the coordinates

        if (recvPts[i].prec > 52)
        { // convert to _mp
          setprec_point_mp(orig_last_mp, recvPts[i].prec);
          point_d_to_mp(orig_last_mp, orig_last_d);
        }
      }
      else
      {
        setprec_point_mp(orig_last_mp, recvPts[i].last_approx_prec);
        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // convert to extrinsic coordinates
          intrinsicToExtrinsic_mp(orig_last_mp, recvPts[i].last_approx_mp, ED_d->EqD->stageData_mp[stage].B0, ED_d->EqD->stageData_mp[stage].p0);
        }
        else
        { // copy over
          point_cp_mp(orig_last_mp, recvPts[i].last_approx_mp);
        }
        orig_last_d->size /= 2; // remove the bottom half of the coordinates

        if (recvPts[i].prec < 64)
        { // convert to _d
          point_mp_to_d(orig_last_d, orig_last_mp);
        }
      }

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        set_d(ED_d->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_d.time);

        if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store the intrinsic endpoint and its corresponding point in the original variables
          intrinsicToExtrinsic_d(dehom_d, recvPts[i].PD_d.point, ED_d->EqD->stageData_d[stage].B0, ED_d->EqD->stageData_d[stage].p0);
          dehom_d->size /= 2; // remove the bottom half of the coordinates

          point_cp_d(ED_d->EqD->stageData_d[stage].endPts[path_num], dehom_d);

          // convert to intrinsic coordinates
          mat_d B_transpose;
          init_mat_d(B_transpose, ED_d->EqD->stageData_d[stage].B->cols, ED_d->EqD->stageData_d[stage].B->rows);
          transpose_d(B_transpose, ED_d->EqD->stageData_d[stage].B);

          extrinsicToIntrinsic_d(ED_d->EqD->stageData_d[stage].endPts_in[path_num], dehom_d, B_transpose, ED_d->EqD->stageData_d[stage].p);

          clear_mat_d(B_transpose);
        }
        else
        { // use the top set of coordinates so that we can store the original coordinates
          point_cp_d(ED_d->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_d.point);
          ED_d->EqD->stageData_d[stage].endPts[path_num]->size /= 2;
        }

        // find dehom_d
        getDehomPoint_d(dehom_d, ED_d->EqD->stageData_d[stage].endPts[path_num], ED_d->EqD->stageData_d[stage].endPts[path_num]->size, &ED_d->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->stageData_d[stage].endPt_retVals[path_num] = printStageFooter_d(ED_d, stage, path_num, &recvPts[i].PD_d, ED_d->EqD->stageData_d[stage].endPts[path_num], orig_last_d, dehom_d, recvPts[i].condition_number, recvPts[i].function_residual_d, recvPts[i].latest_newton_residual_d, recvPts[i].t_val_at_latest_sample_point_d, recvPts[i].error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
      else
      { // set the precision
        T->Precision = recvPts[i].prec;
        initMP(recvPts[i].prec);
        change_prec(ED_mp, recvPts[i].prec);

        // setup dehom_mp & tempPoint_mp
        setprec_point_mp(dehom_mp, recvPts[i].prec);
        setprec_point_mp(tempPoint_mp, recvPts[i].prec);

        // copy over to the appropriate spot - converting to double precision
        mp_to_d(ED_d->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_mp.time);

        if (ED_mp->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store the intrinsic endpoint and its corresponding point in the original variables
          intrinsicToExtrinsic_mp(tempPoint_mp, recvPts[i].PD_mp.point, ED_d->EqD->stageData_mp[stage].B0, ED_d->EqD->stageData_mp[stage].p0);
          tempPoint_mp->size /= 2; // remove the bottom half of the coordinates

          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts[path_num], tempPoint_mp);

          // convert to the original coordinates
          mat_mp B_transpose;
          init_mat_mp2(B_transpose, ED_d->EqD->stageData_mp[stage].B->cols, ED_d->EqD->stageData_mp[stage].B->rows, recvPts[i].prec);
          transpose_mp(B_transpose, ED_d->EqD->stageData_mp[stage].B);

          extrinsicToIntrinsic_mp(dehom_mp, tempPoint_mp, B_transpose, ED_d->EqD->stageData_mp[stage].p);

          // copy to endPts_in
          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts_in[path_num], dehom_mp);

          clear_mat_mp(B_transpose);
        }
        else
        { // use the top set of coordinates so that we can store the original coordinates
          point_cp_mp(tempPoint_mp, recvPts[i].PD_mp.point);
          tempPoint_mp->size /= 2;

          point_mp_to_d(ED_d->EqD->stageData_d[stage].endPts[path_num], tempPoint_mp);
        }

        // find dehom_mp
        getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

        // print the footer to OUT for the point
        ED_d->EqD->stageData_d[stage].endPt_retVals[path_num] = printStageFooter_mp(ED_mp, stage, path_num, &recvPts[i].PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
      }
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  free(recvPts);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);

  clear_point_d(dehom_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last_mp);
  clear_point_mp(tempPoint_mp);

  return;
}

void head_eqbyeqSortStageEndpoints_d(basic_eval_data_d *BED, tracker_config_t *T, int stage, double final_tol, int pathMod, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at for the current stage     *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, indexI = 0, indexJ = 0, count, path_num, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize, num_paths = BED->EqD->stageData_d[stage].num_paths;
  int num_stages = BED->EqD->num_subsystems, num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int)), *higherDim = NULL;
  vec_d tempVec;
  sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  init_vec_d(tempVec, 0);

  // print header for level
  printf("\nSorting points for stage %d of %d: %d path%s to sort.\n", stage, num_stages, num_paths, num_paths == 1 ? "" : "s");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // randomize the list
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to j
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
    init_endgame_data(&sendPts[i], 64);
  higherDim = (int *)bmalloc(maxSize * sizeof(int));
  for (i = 0; i < maxSize; i++)
    higherDim[i] = 0;
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Sorting %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = 52; // point in double precision
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = BED->EqD->stageData_d[stage].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
        if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
        { // sort using the intrinsic coordinates
          point_cp_d(sendPts[j].PD_d.point, BED->EqD->stageData_d[stage].endPts_in[path_num]);
        }
        else
        { // sort using original coordinates
          point_cp_d(sendPts[j].PD_d.point, BED->EqD->stageData_d[stage].endPts[path_num]);
        }
        set_d(sendPts[j].PD_d.time, BED->EqD->stageData_d[stage].finalTs[path_num]);
        sendPts[j].last_approx_prec = 52;
        sendPts[j].last_approx_d->size = 0;
        higherDim[j] = BED->EqD->stageData_d[stage].higherDim[path_num];
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      // send higherDim
      MPI_Send(higherDim, lastSize[i], MPI_INT, i, TAG_INT, MPI_COMM_WORLD);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be sorted
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Sorting %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = 52;  // start point in double precision
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = BED->EqD->stageData_d[stage].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
      if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
      { // sort using the intrinsic coordinates
        point_cp_d(sendPts[j].PD_d.point, BED->EqD->stageData_d[stage].endPts_in[path_num]);
      }
      else
      { // sort using original coordinates
        point_cp_d(sendPts[j].PD_d.point, BED->EqD->stageData_d[stage].endPts[path_num]);
      }
      set_d(sendPts[j].PD_d.time, BED->EqD->stageData_d[stage].finalTs[path_num]);
      sendPts[j].last_approx_prec = 52;
      sendPts[j].last_approx_d->size = 0;
      higherDim[j] = BED->EqD->stageData_d[stage].higherDim[path_num];
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    // send higherDim
    MPI_Send(higherDim, lastSize[recvProc], MPI_INT, recvProc, TAG_INT, MPI_COMM_WORLD);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->stageData_d[stage].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->stageData_d[stage].condition_nums[path_num] = recvPts[i].condition_number;

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store to the intrinsic coordinates
          point_cp_d(BED->EqD->stageData_d[stage].endPts_in[path_num], recvPts[i].PD_d.point);
        }
        else
        { // store using original coordinates
          point_cp_d(BED->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_d.point);
        }
        set_d(BED->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_d.time);
      }
      else
      { // copy over to the appropriate spot - converting to double precision
        if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store to the intrinsic coordinates
          point_mp_to_d(BED->EqD->stageData_d[stage].endPts_in[path_num], recvPts[i].PD_mp.point);
        }
        else
        { // store using original coordinates
          point_mp_to_d(BED->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_mp.point);
        }
        mp_to_d(BED->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_mp.time);
      }

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      sortPts[path_num].norm = infNormVec_d(BED->EqD->stageData_d[stage].endPts[path_num]);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->stageData_d[stage].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->stageData_d[stage].condition_nums[path_num] = recvPts[i].condition_number;

      if (recvPts[i].prec < 64)
      { // copy over to the appropriate spot
        if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store to the intrinsic coordinates
          point_cp_d(BED->EqD->stageData_d[stage].endPts_in[path_num], recvPts[i].PD_d.point);
        }
        else
        { // store using original coordinates
          point_cp_d(BED->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_d.point);
        }
        set_d(BED->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_d.time);
      }
      else
      { // copy over to the appropriate spot - converting to double precision
        if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
        { // store to the intrinsic coordinates
          point_mp_to_d(BED->EqD->stageData_d[stage].endPts_in[path_num], recvPts[i].PD_mp.point);
        }
        else
        { // store using original coordinates
          point_mp_to_d(BED->EqD->stageData_d[stage].endPts[path_num], recvPts[i].PD_mp.point);
        }
        mp_to_d(BED->EqD->stageData_d[stage].finalTs[path_num], recvPts[i].PD_mp.time);
      }

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      sortPts[path_num].norm = infNormVec_d(BED->EqD->stageData_d[stage].endPts[path_num]);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // sort sortPts
  qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

  // do the final classification
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (BED->EqD->stageData_d[stage].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      j = i + 1;
      while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
      {
        indexJ = sortPts[j].path_num;
        if (BED->EqD->stageData_d[stage].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          if (BED->EqD->stageData_d[stage].useIntrinsicSlice)
          { // subtract using intrinsic coordinates
            vec_sub_d(tempVec, BED->EqD->stageData_d[stage].endPts_in[indexI], BED->EqD->stageData_d[stage].endPts_in[indexJ]);
          }
          else
          { // subtract using extrinsic coordinates
            vec_sub_d(tempVec, BED->EqD->stageData_d[stage].endPts[indexI], BED->EqD->stageData_d[stage].endPts[indexJ]);
          }

          if (infNormVec_d(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            BED->EqD->stageData_d[stage].endPt_types[indexI] = BED->EqD->stageData_d[stage].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
        j++;
      }
    }

    // add to count
    if (BED->EqD->stageData_d[stage].endPt_types[indexI] == retVal_going_to_infinity || BED->EqD->stageData_d[stage].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (BED->EqD->stageData_d[stage].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (BED->EqD->stageData_d[stage].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (BED->EqD->stageData_d[stage].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  BED->EqD->stageData_d[stage].num_sing = num_sing;
  BED->EqD->stageData_d[stage].num_nonsing = num_nonsing;
  BED->EqD->stageData_d[stage].num_higher_dim = num_higher_dim;
  BED->EqD->stageData_d[stage].num_bad = num_bad;
  BED->EqD->stageData_d[stage].num_inf = num_inf;

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  free(recvPts);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  free(sortPts);

  clear_vec_d(tempVec);

  return;
}

/////////// WORKER FUNCTIONS ////////////////

void worker_eqbyeq_d(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does equation-by-equation - 'worker' process           *
\***************************************************************/
{
  int i, size, num_subs;
  trackingStats trackCount;
  tracker_config_t T;
  basic_eval_data_d BED;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // recv T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision - does not hurt if doing only double precision
  initMP(T.Precision); 

  if (T.MPType == 2)
  { // recv BED - using AMP
    bcast_basic_eval_data_amp(&BED, my_id, headnode);
  }
  else
  { // recv ED - only double
    bcast_basic_eval_data_d(&BED, T.MPType, my_id, headnode);
  }

  // allocate for eq-by-eq
  BED.EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));

  // recv EqD
  bcast_eqData_t(BED.EqD, T.MPType, my_id, headnode);

  // store the number of subsystems
  num_subs = BED.EqD->num_subsystems;

  // allocate memory for subsystems
  BED.EqD->witnessData_d = (eqWitnessData_d *)bmalloc(num_subs * sizeof(eqWitnessData_d));
  BED.EqD->stageData_d = (eqStageData_d *)bmalloc(num_subs * sizeof(eqStageData_d));
  if (T.MPType == 2)
  { // allocate levels im MP
    BED.EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(num_subs * sizeof(eqWitnessData_mp));
    BED.EqD->stageData_mp = (eqStageData_mp *)bmalloc(num_subs * sizeof(eqStageData_mp));

    // point to EqD from inside BED->BED_mp
    BED.BED_mp->EqD = BED.EqD;
  }

  // recv witness information for all of the subsystems
  for (i = 0; i < num_subs; i++)
  { // recv the ith witness data from the head
    bcast_witnessData(BED.EqD, i, T.MPType, my_id, headnode);
  }
  // setup to increase precision on witness data
  BED.EqD->increase_witness_prec = 1;

  // setup the local files - OUT, MIDOUT, RAWOUT & FAIL
  size = 1 + snprintf(NULL, 0, "output_%d", my_id);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "output_%d", my_id);
  OUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "midout_w%d_%d", 0, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "midout_w%d_%d", 0, my_id);
  MIDOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  RAWOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  FAIL = fopen(str, "w");

  // loop to find the witness points for each subsystem
  for (i = 0; i < num_subs; i++)
  { // track this level
    worker_eqbyeqWitnessTrack_d(&trackCount, i, &T, &BED, BED.BED_mp, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    if (num_subs > 1)
    { // close MIDOUT
      fclose(MIDOUT);
      MIDOUT = NULL;
      // wait until all workers have closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel
    }

    // sort the endpoints
    worker_eqbyeqSortWitnessEndpoints_d(&BED, BED.BED_mp, &T, OUT, i, my_id, headnode, num_processes);

    // if this is not the last subsystem, repoen MIDOUT
    if (i + 1 < num_subs)
    { // reopen MIDOUT
      size = 1 + snprintf(NULL, 0, "midout_w%d_%d", i + 1, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_w%d_%d", i + 1, my_id);
      MIDOUT = fopen(str, "w");
    }
  }

  // clear the first witness data
  clearEqbyEqWitnessData_d(BED.EqD, 0, T.MPType);

  // setup to increase precision on stage data
  BED.EqD->increase_witness_prec = 0;

  // loop to track the stages
  for (i = 1; i < num_subs; i++)
  { // recv the ith stage from the head
    bcast_stageData(BED.EqD, i, T.MPType, my_id, headnode);

    if (i > 1)
    { // clear the stage data from stage 'i - 1' since it is no longer needed
      clearEqbyEqStageData_d(BED.EqD, i - 1, T.MPType);
    }

    // clear the witness data from subsystem 'stage' since it is no longer needed
    clearEqbyEqWitnessData_d(BED.EqD, i, T.MPType);

    // open MIDOUT
    size = 1 + snprintf(NULL, 0, "midout_s%d_%d", i, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "midout_s%d_%d", i, my_id);
    MIDOUT = fopen(str, "w");

    // track this stage
    worker_eqbyeqStageTrack_d(&trackCount, i, &T, &BED, BED.BED_mp, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    if (i + 1 < num_subs)
    { // close MIDOUT
      fclose(MIDOUT);
      MIDOUT = NULL;
      // wait until all workers have closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel
    }

    // sort the endpoints
    worker_eqbyeqSortStageEndpoints_d(&BED, BED.BED_mp, &T, OUT, i, my_id, headnode, num_processes);
  }

  // close the files
  fclose(OUT);
  fclose(MIDOUT);
  fclose(RAWOUT);
  fclose(FAIL);

  // clear other allocated memory
  basic_eval_clear_d(&BED, -59, T.MPType); // -59 since this does use equation-by-equation
  tracker_config_clear(&T);
  clearMP();
  free(str);

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

void worker_eqbyeqWitnessTrack_d(trackingStats *trackCount, int curr_sub, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this subsystem in parallel *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0, num_subs = ED_d->EqD->num_subsystems;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  // setup the function evaluators
  eval_func_d = &witness_eqbyeq_eval_d;
  eval_func_mp = &witness_eqbyeq_eval_mp;
  change_prec = &change_eqbyeq_eval_prec;
  find_dehom = &eqbyeq_witness_dehom;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Finding witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], 64);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED_d->EqD->curr_stage_num = curr_sub;
    if (T->MPType == 2)
      ED_mp->EqD->curr_stage_num = curr_sub;

    for (i = 0; i < numStartPts; i++)
    { // setup for tracking the path
      ED_d->EqD->curr_path_num = startPts[i].pathNum;
      if (T->MPType == 2)
        ED_mp->EqD->curr_path_num = startPts[i].pathNum;

      T->first_step_of_path = 1;
      T->endgameOnly = 0;

      // print the header of the path to OUT
      printPathHeader_d(OUT, &startPts[i].PD_d, T, startPts[i].pathNum, ED_d, eval_func_d);

      // track the path
      zero_dim_track_path_d(startPts[i].pathNum, &endPts[i], &startPts[i].PD_d, OUT, MIDOUT, T, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

      // check to see if it should be sharpened - only when this is the only subsystems
      if (num_subs == 1 && endPts[i].retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      }

      printResultOfPath(OUT, endPts[i].retVal, T);

      if (endPts[i].prec < 64)
      { // print footer in double precision
        printBasicFooter_d(OUT, &endPts[i].PD_d, T, endPts[i].function_residual_d);
      }
      else
      { // print footer in multi precision
        printBasicFooter_mp(OUT, &endPts[i].PD_mp, T, endPts[i].function_residual_mp);
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], 64);
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqSortWitnessEndpoints_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, int curr_sub, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points for the current subsystem     *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, rankDef, finite, numStartPts = 0, numEndPts = 0, *higherDim = NULL;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  MPI_Status status;

  // setup the function evaluators
  eval_func_d = &witness_eqbyeq_eval_d;
  eval_func_mp = &witness_eqbyeq_eval_mp;
  change_prec = &change_eqbyeq_eval_prec;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
  // recv higherDim
  higherDim = (int *)bmalloc(numStartPts * sizeof(int));
  MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], 64);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED_d->EqD->curr_stage_num = curr_sub;
    if (T->MPType == 2)
      ED_mp->EqD->curr_stage_num = curr_sub;

    for (i = 0; i < numStartPts; i++)
    { // setup for evaluation
      ED_d->EqD->curr_path_num = startPts[i].pathNum;
      if (T->MPType == 2)
        ED_mp->EqD->curr_path_num = startPts[i].pathNum;

      // setup endPts
      endPts[i].prec = 52;
      endPts[i].pathNum = startPts[i].pathNum;
      point_cp_d(endPts[i].PD_d.point, startPts[i].PD_d.point);
      set_d(endPts[i].PD_d.time, startPts[i].PD_d.time);
      endPts[i].last_approx_prec = 52;
      endPts[i].last_approx_d->size = 0;

      if (startPts[i].retVal == 0)
      { // set time to 0
        set_zero_d(endPts[i].PD_d.time);

        // determine if it is rank deficient, finite, higherDim
        eqbyeqWitnessSortEndpoint_d(&rankDef, &finite, higherDim[i], ED_d, ED_mp, curr_sub, startPts[i].pathNum, &endPts[i].condition_number, T, OUT, &endPts[i].PD_d, &endPts[i].PD_mp, endPts[i].prec, eval_func_d, eval_func_mp, change_prec); 

        // check for success on the jacobian analyzer
        if (T->regen_remove_inf && !finite)
        { // dehom point is infinite
          endPts[i].retVal = retVal_going_to_infinity;
        }
        else if (T->regen_higher_dim_check && higherDim[i])
        { // it lies on a higher dimensional component
          endPts[i].retVal = retVal_higher_dim;
        }
        else if (!rankDef)
        { // classify as non-singular
          endPts[i].retVal = MOVE_TO_NEXT;
        }
        else
        { // classify as singular
          endPts[i].retVal = DO_NOT_MOVE_TO_NEXT;
        }
      }
      else
      { // path was not a success - copy over error code
        endPts[i].retVal = startPts[i].retVal;
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
    // recv higherDim
    higherDim = (int *)brealloc(higherDim, numStartPts * sizeof(int));
    if (numStartPts > 0)
      MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], 64);
    }
  }

  // startPts & endPts & higherDim are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqStageTrack_d(trackingStats *trackCount, int curr_stage, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this stage in parallel     *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0, num_stages = ED_d->EqD->num_subsystems;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  // setup the function evaluators
  eval_func_d = &standard_eqbyeq_eval_d;
  eval_func_mp = &standard_eqbyeq_eval_mp;
  change_prec = &change_eqbyeq_eval_prec;
  find_dehom = &eqbyeq_stage_dehom;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking points for stage %d.\n", curr_stage);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], 64);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED_d->EqD->curr_stage_num = curr_stage;
    if (T->MPType == 2)
      ED_mp->EqD->curr_stage_num = curr_stage;

    for (i = 0; i < numStartPts; i++)
    { // setup for tracking the path
      ED_d->EqD->curr_path_num = startPts[i].pathNum;
      if (T->MPType == 2)
        ED_mp->EqD->curr_path_num = startPts[i].pathNum;

      T->first_step_of_path = 1;
      T->endgameOnly = 0;

      // print the header of the path to OUT
      printPathHeader_d(OUT, &startPts[i].PD_d, T, startPts[i].pathNum, ED_d, eval_func_d);

      // track the path
      zero_dim_track_path_d(startPts[i].pathNum, &endPts[i], &startPts[i].PD_d, OUT, MIDOUT, T, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);


      // check to see if it should be sharpened - only when this is the only subsystems
      if (curr_stage + 1 == num_stages && endPts[i].retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      }

      printResultOfPath(OUT, endPts[i].retVal, T);

      if (endPts[i].prec < 64)
      { // print footer in double precision
        printBasicFooter_d(OUT, &endPts[i].PD_d, T, endPts[i].function_residual_d);
      }
      else
      { // print footer in multi precision
        printBasicFooter_mp(OUT, &endPts[i].PD_mp, T, endPts[i].function_residual_mp);
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], 64);
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqSortStageEndpoints_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, int curr_stage, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the points for the current stage                 *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, rankDef, finite, numStartPts = 0, numEndPts = 0, *higherDim = NULL;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  MPI_Status status;

  // setup the function evaluators
  eval_func_d = &stage_sort_eqbyeq_eval_d;
  eval_func_mp = &stage_sort_eqbyeq_eval_mp;
  change_prec = &change_eqbyeq_eval_prec;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting points for stage %d.\n", curr_stage);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
  // recv higherDim
  higherDim = (int *)bmalloc(numStartPts * sizeof(int));
  MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], 64);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED_d->EqD->curr_stage_num = curr_stage;
    if (T->MPType == 2)
      ED_mp->EqD->curr_stage_num = curr_stage;

    for (i = 0; i < numStartPts; i++)
    { // setup for evaluation
      ED_d->EqD->curr_path_num = startPts[i].pathNum;
      if (T->MPType == 2)
        ED_mp->EqD->curr_path_num = startPts[i].pathNum;

      // setup endPts
      endPts[i].prec = 52;
      endPts[i].pathNum = startPts[i].pathNum;
      point_cp_d(endPts[i].PD_d.point, startPts[i].PD_d.point);
      set_d(endPts[i].PD_d.time, startPts[i].PD_d.time);
      endPts[i].last_approx_prec = 52;
      endPts[i].last_approx_d->size = 0;

      if (startPts[i].retVal == 0)
      { // set time to 0
        set_zero_d(endPts[i].PD_d.time);

        // determine if it is rank deficient, finite, higherDim
        eqbyeqStageSortEndpoint_d(&rankDef, &finite, higherDim[i], ED_d, ED_mp, curr_stage, startPts[i].pathNum, &endPts[i].condition_number, T, OUT, &endPts[i].PD_d, &endPts[i].PD_mp, endPts[i].prec, eval_func_d, eval_func_mp, change_prec);

        // check for success
        if (T->regen_remove_inf && !finite)
        { // dehom point is infinite
          endPts[i].retVal = retVal_going_to_infinity;
        }
        else if (T->regen_higher_dim_check && higherDim[i])
        { // it lies on a higher dimensional component
          endPts[i].retVal = retVal_higher_dim;
        }
        else if (!rankDef)
        { // classify as non-singular
          endPts[i].retVal = MOVE_TO_NEXT;
        }
        else
        { // classify as singular
          endPts[i].retVal = DO_NOT_MOVE_TO_NEXT;
        }
      }
      else
      { // path was not a success - copy over error code
        endPts[i].retVal = startPts[i].retVal;
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
    // recv higherDim
    higherDim = (int *)brealloc(higherDim, numStartPts * sizeof(int));
    if (numStartPts > 0)
      MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], 64);
    }
  }

  // startPts & endPts & higherDim are cleared since numStartPts, numEndPts <=0

  return;
}

///////////////////////// MP Parallel Functions //////////////////////////////

void head_eqbyeq_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_mp *ED, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does eq-by-eq tracking - 'headnode'                    *
\***************************************************************/
{
  int stage, num_crossings = 0, size, num_paths, num_vars, num_subs = ED->EqD->num_subsystems;
  char *str = NULL, nonName[] = "nonsolutions";
  FILE *NONSOLN = NULL;

  if (!ED->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen(nonName, "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // initialize trackCount
  init_trackingStats(trackCount);

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  // send ED to the workers - using MP
  bcast_basic_eval_data_mp(ED, 1, my_id, headnode);

  // send eqData_t to the workers
  bcast_eqData_t(ED->EqD, T->MPType, my_id, headnode);
  // send witness information to the workers
  for (stage = 0; stage < num_subs; stage++)
  { // send the 'stage'th witness data to the workers
    bcast_witnessData(ED->EqD, stage, T->MPType, my_id, headnode);
  }

  // loop to generate the witness sets for each of the subsystems
  for (stage = 0; stage < num_subs; stage++)
  { // find the number of paths and the number of variables used in this subsystem
    num_paths = ED->EqD->witnessData_mp[stage].num_paths;
    num_vars = ED->EqD->witnessData_mp[stage].depth;

    // find the witness points for this subsystem
    head_eqbyeqWitnessTrack_mp(trackCount, stage, pathMod, T, ED, OUT, RAWOUT, FAIL, NONSOLN, my_id, headnode, num_processes);

    // check for path crossings if this is not the only subsystem
    if (num_subs > 1)
    { // wait unitl the files are close
      MPI_Barrier(MPI_COMM_WORLD);

      // setup str
      size = 1 + snprintf(NULL, 0, "midout_w%d", stage);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_w%d", stage);

      num_crossings = parallel_midpoint_checking(midFile, str, 1, num_paths, num_vars, midpoint_tol, my_id, headnode, num_processes);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // sort the endpoints
    head_eqbyeqSortWitnessEndpoints_mp(ED, T, stage, target_tol, pathMod, my_id, headnode, num_processes);
  }

  // setup the first stage
  setupEqbyEqFirstStage_mp(ED->EqD);

  // clear the first witness data
  clearEqbyEqFirstWitnessData_mp(ED->EqD);

  // loop to track each of the stages
  for (stage = 1; stage < num_subs; stage++)
  { // setup the next stage
    setupEqbyEqNextStage_mp(ED, stage);

    // send the stage information to the workers
    bcast_stageData(ED->EqD, stage, T->MPType, my_id, headnode);

    // clear the stage data from stage 'stage - 1' since it is no longer needed
    clearEqbyEqStageData_mp(ED->EqD, stage - 1);

    // clear the witness data from subsystem 'stage' since it is no longer needed
    clearEqbyEqWitnessData_mp(ED->EqD, stage);

    // find the number of paths for this next stage
    num_paths = ED->EqD->stageData_mp[stage].num_paths;
    if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
      num_vars = ED->EqD->stageData_mp[stage].depth_x + ED->EqD->stageData_mp[stage].depth_y;
    else
      num_vars = 2 * ED->EqD->num_vars;

    // track each path for this stage
    head_eqbyeqStageTrack_mp(trackCount, stage, pathMod, T, ED, OUT, RAWOUT, FAIL, NONSOLN, my_id, headnode, num_processes);

    // check for path crossings if this is not the last subsystem
    if (stage + 1 < num_subs)
    { // wait unitl the files are close
      MPI_Barrier(MPI_COMM_WORLD);

      // setup str
      size = 1 + snprintf(NULL, 0, "midout_s%d", stage);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_s%d", stage);

      num_crossings = parallel_midpoint_checking(midFile, str, 1, num_paths, num_vars, midpoint_tol, my_id, headnode, num_processes);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // sort the endpoints
    head_eqbyeqSortStageEndpoints_mp(ED, T, stage, target_tol, pathMod, my_id, headnode, num_processes);
  }

  // set the number of points that were tracked on the last level
  trackCount->numPoints = ED->EqD->stageData_mp[num_subs - 1].num_paths;

  // wait until the files are closed
  MPI_Barrier(MPI_COMM_WORLD);

  // setup str
  if (num_subs == 1)
    size = 1 + snprintf(NULL, 0, "midout_w%d", num_subs - 1);
  else
    size = 1 + snprintf(NULL, 0, "midout_s%d", num_subs - 1);
  str = (char *)brealloc(str, size * sizeof(char));
  if (num_subs == 1)
    sprintf(str, "midout_w%d", num_subs - 1);
  else
    sprintf(str, "midout_s%d", num_subs - 1);

  // combine the midpath_data files
  combine_midpath_data(midFile, str, 1, headnode, num_processes);

  // delete the rawout & fail worker files
  delete_parallel_files("rawout", headnode, num_processes);
  delete_parallel_files("fail", headnode, num_processes);

  if (!ED->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // clear memory
  free(str);

  return;
}

void head_eqbyeqWitnessTrack_mp(trackingStats *trackCount, int curr_sub, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the witness points for the current subsystem     *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize;
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  int path_num, num_subs = ED->EqD->num_subsystems, num_paths = ED->EqD->witnessData_mp[curr_sub].num_paths;

  point_mp dehom_mp, orig_last;
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // initialize dehom_mp
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last, 0);

  // setup the subsystem
  ED->EqD->curr_stage_num = curr_sub;

  // print header for level
  printf("\nFinding witness points for subsystem %d of %d: %d path%s to track.\n", curr_sub, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Finding witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum - randomize the list of start points
  pathNum = (int *)brealloc(pathNum, num_paths * sizeof(int));
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to i
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  }

  // initialize count
  count = 0;

  // send out the initial set of packets for this level
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = T->Precision;
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = ED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num];
        point_cp_mp(sendPts[j].PD_mp.point, ED->EqD->witnessData_mp[curr_sub].startPts[path_num]);
        set_one_mp(sendPts[j].PD_mp.time);
        sendPts[j].last_approx_prec = T->Precision;
        sendPts[j].last_approx_mp->size = 0;
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Tracking path %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = T->Precision;
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = ED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num];
      point_cp_mp(sendPts[j].PD_mp.point, ED->EqD->witnessData_mp[curr_sub].startPts[path_num]);
      set_one_mp(sendPts[j].PD_mp.time);
      sendPts[j].last_approx_prec = T->Precision;
      sendPts[j].last_approx_mp->size = 0;
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED->EqD->witnessData_mp[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED->EqD->witnessData_mp[curr_sub].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, NULL, ED, curr_sub, 0);

      // copy over to the appropriate spot
      point_cp_mp(ED->EqD->witnessData_mp[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
      set_mp(ED->EqD->witnessData_mp[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

      // covert to the original coordinates
      intrinsicToExtrinsic_mp(orig_last, recvPts[i].last_approx_mp, ED->EqD->witnessData_mp[curr_sub].B, ED->EqD->witnessData_mp[curr_sub].p);
      intrinsicToExtrinsic_mp(ED->EqD->witnessData_mp[curr_sub].endPts[path_num], recvPts[i].PD_mp.point, ED->EqD->witnessData_mp[curr_sub].B, ED->EqD->witnessData_mp[curr_sub].p);

      // find dehom_mp
      getDehomPoint_mp(dehom_mp, ED->EqD->witnessData_mp[curr_sub].endPts[path_num], ED->EqD->witnessData_mp[curr_sub].endPts[path_num]->size, &ED->preProcData);

      // print the footer to OUT for the point
      ED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num] = printWitnessFooter_mp(ED, curr_sub, path_num, &recvPts[i].PD_mp, ED->EqD->witnessData_mp[curr_sub].endPts[path_num], orig_last, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED->EqD->witnessData_mp[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED->EqD->witnessData_mp[curr_sub].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, NULL, ED, curr_sub, 0);

      // copy over to the appropriate spot
      point_cp_mp(ED->EqD->witnessData_mp[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
      set_mp(ED->EqD->witnessData_mp[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

      // covert to the original coordinates
      intrinsicToExtrinsic_mp(orig_last, recvPts[i].last_approx_mp, ED->EqD->witnessData_mp[curr_sub].B, ED->EqD->witnessData_mp[curr_sub].p);
      intrinsicToExtrinsic_mp(ED->EqD->witnessData_mp[curr_sub].endPts[path_num], recvPts[i].PD_mp.point, ED->EqD->witnessData_mp[curr_sub].B, ED->EqD->witnessData_mp[curr_sub].p);

      // find dehom_mp
      getDehomPoint_mp(dehom_mp, ED->EqD->witnessData_mp[curr_sub].endPts[path_num], ED->EqD->witnessData_mp[curr_sub].endPts[path_num]->size, &ED->preProcData);

      // print the footer to OUT for the point
      ED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num] = printWitnessFooter_mp(ED, curr_sub, path_num, &recvPts[i].PD_mp, ED->EqD->witnessData_mp[curr_sub].endPts[path_num], orig_last, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last);
  free(recvPts);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  
  return;
}

void head_eqbyeqSortWitnessEndpoints_mp(basic_eval_data_mp *BED, tracker_config_t *T, int curr_sub, double final_tol, int pathMod, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at for the current subsystem *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, indexI = 0, indexJ = 0, count, cont, path_num, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize, num_paths = BED->EqD->witnessData_mp[curr_sub].num_paths;
  int  num_subs = BED->EqD->num_subsystems, num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int)), *higherDim = NULL;
  mpf_t norm_diff;
  vec_mp tempVec;
  sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  mpf_init(norm_diff);
  init_vec_mp(tempVec, 0);

  // print header for level
  printf("\nSorting witness points for subsystem %d of %d: %d path%s to sort.\n", curr_sub, num_subs, num_paths, num_paths == 1 ? "" : "s");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // randomize the list
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to j
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  }
  higherDim = (int *)bmalloc(maxSize * sizeof(int));
  for (i = 0; i < maxSize; i++)
    higherDim[i] = 0;
 
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Sorting %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = T->Precision;
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = BED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
        point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->witnessData_mp[curr_sub].endPts_in[path_num]);
        set_zero_mp(sendPts[j].PD_mp.time);
        sendPts[j].last_approx_prec = T->Precision;
        sendPts[j].last_approx_mp->size = 0;
        higherDim[j] = BED->EqD->witnessData_mp[curr_sub].higherDim[path_num];
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      // send higherDim
      MPI_Send(higherDim, lastSize[i], MPI_INT, i, TAG_INT, MPI_COMM_WORLD);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be sorted
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Sorting %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = T->Precision;
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = BED->EqD->witnessData_mp[curr_sub].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
      point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->witnessData_mp[curr_sub].endPts_in[path_num]);
      set_zero_mp(sendPts[j].PD_mp.time);
      sendPts[j].last_approx_prec = T->Precision;
      sendPts[j].last_approx_mp->size = 0;
      higherDim[j] = BED->EqD->witnessData_mp[curr_sub].higherDim[path_num];
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    // send higherDim
    MPI_Send(higherDim, lastSize[recvProc], MPI_INT, recvProc, TAG_INT, MPI_COMM_WORLD);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->witnessData_mp[curr_sub].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->witnessData_mp[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // copy over to the appropriate spot
      point_cp_mp(BED->EqD->witnessData_mp[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
      set_mp(BED->EqD->witnessData_mp[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      mpf_init(sortPts[path_num].norm);
      infNormVec_mp2(sortPts[path_num].norm, BED->EqD->witnessData_mp[curr_sub].endPts[path_num]);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->witnessData_mp[curr_sub].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->witnessData_mp[curr_sub].condition_nums[path_num] = recvPts[i].condition_number;

      // copy over to the appropriate spot
      point_cp_mp(BED->EqD->witnessData_mp[curr_sub].endPts_in[path_num], recvPts[i].PD_mp.point);
      set_mp(BED->EqD->witnessData_mp[curr_sub].finalTs[path_num], recvPts[i].PD_mp.time);

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      mpf_init(sortPts[path_num].norm);
      infNormVec_mp2(sortPts[path_num].norm, BED->EqD->witnessData_mp[curr_sub].endPts[path_num]);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // sort sortPts
  qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

  // do the final classification
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      cont = 1;
      j = i;
      do
      { // increment the counter (start at i + 1)
        j++;

        if (j < num_paths)
        { // check norm diff
          indexJ = sortPts[j].path_num;
          mpf_sub(norm_diff, sortPts[j].norm, sortPts[i].norm);
          if (mpf_get_d(norm_diff) > final_tol)
            cont = 0;
        }
        else
          cont = 0;

        // check to see if we can continue
        if (cont && BED->EqD->witnessData_mp[curr_sub].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          vec_sub_mp(tempVec, BED->EqD->witnessData_mp[curr_sub].endPts_in[indexI], BED->EqD->witnessData_mp[curr_sub].endPts_in[indexJ]);
          if (infNormVec_mp(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] = BED->EqD->witnessData_mp[curr_sub].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
      } while (cont);
    }

    // add to count
    if (BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == retVal_going_to_infinity || BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (BED->EqD->witnessData_mp[curr_sub].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  BED->EqD->witnessData_mp[curr_sub].num_sing = num_sing;
  BED->EqD->witnessData_mp[curr_sub].num_nonsing = num_nonsing;
  BED->EqD->witnessData_mp[curr_sub].num_higher_dim = num_higher_dim;
  BED->EqD->witnessData_mp[curr_sub].num_bad = num_bad;
  BED->EqD->witnessData_mp[curr_sub].num_inf = num_inf;

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  for (i = num_paths - 1; i >= 0; i--)
    mpf_clear(sortPts[i].norm);
  free(recvPts);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  free(sortPts);
  free(higherDim);

  mpf_clear(norm_diff);
  clear_vec_mp(tempVec);

  return;
}

void head_eqbyeqStageTrack_mp(trackingStats *trackCount, int stage, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the current stage of equation-by-equation       *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize;
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  int path_num, num_stages = ED->EqD->num_subsystems, num_paths = ED->EqD->stageData_mp[stage].num_paths;

  point_mp dehom_mp, orig_last;
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // initialize dehom_mp
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last, 0);

  // setup the subsystem
  ED->EqD->curr_stage_num = stage;

  // print header for level
  printf("\nTracking points for stage %d of %d: %d path%s to track.\n", stage, num_stages, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking points for stage %d.\n", stage);
  fprintf(OUT, "*****************************************************\n");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum - randomize the list of start points
  pathNum = (int *)brealloc(pathNum, num_paths * sizeof(int));
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to i
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  }

  // initialize count
  count = 0;

  // send out the initial set of packets for this level
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = T->Precision;
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = ED->EqD->stageData_mp[stage].endPt_retVals[path_num];
        point_cp_mp(sendPts[j].PD_mp.point, ED->EqD->stageData_mp[stage].startPts[path_num]);
        set_one_mp(sendPts[j].PD_mp.time);
        sendPts[j].last_approx_prec = T->Precision;
        sendPts[j].last_approx_mp->size = 0;
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Tracking path %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = T->Precision;
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = ED->EqD->stageData_mp[stage].endPt_retVals[path_num];
      point_cp_mp(sendPts[j].PD_mp.point, ED->EqD->stageData_mp[stage].startPts[path_num]);
      set_one_mp(sendPts[j].PD_mp.time);
      sendPts[j].last_approx_prec = T->Precision;
      sendPts[j].last_approx_mp->size = 0;
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED->EqD->stageData_mp[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED->EqD->stageData_mp[stage].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, NULL, ED, stage, 1);

      // copy over to the appropriate spot
      set_mp(ED->EqD->stageData_mp[stage].finalTs[path_num], recvPts[i].PD_mp.time);

      if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
      { // store the intrinsic endpoint and its cooresponding point in the original variables
        intrinsicToExtrinsic_mp(orig_last, recvPts[i].last_approx_mp, ED->EqD->stageData_mp[stage].B0, ED->EqD->stageData_mp[stage].p0);
        orig_last->size /= 2;
        intrinsicToExtrinsic_mp(dehom_mp, recvPts[i].PD_mp.point, ED->EqD->stageData_mp[stage].B0, ED->EqD->stageData_mp[stage].p0);
        dehom_mp->size /= 2; // remove the bottom half of the coordinates

        point_cp_mp(ED->EqD->stageData_mp[stage].endPts[path_num], dehom_mp);

        // convert to intrinsic coordinates
        mat_mp B_transpose;
        init_mat_mp(B_transpose, ED->EqD->stageData_mp[stage].B->cols, ED->EqD->stageData_mp[stage].B->rows);
        transpose_mp(B_transpose, ED->EqD->stageData_mp[stage].B);

        extrinsicToIntrinsic_mp(ED->EqD->stageData_mp[stage].endPts_in[path_num], dehom_mp, B_transpose, ED->EqD->stageData_mp[stage].p);

        clear_mat_mp(B_transpose);
      }
      else
      { // use the top set of coordinates so that we can store the original coordinates
        point_cp_mp(orig_last, recvPts[i].last_approx_mp);
        orig_last->size /= 2;      
        point_cp_mp(ED->EqD->stageData_mp[stage].endPts[path_num], recvPts[i].PD_mp.point);
        ED->EqD->stageData_mp[stage].endPts[path_num]->size /= 2;
      }

      // find dehom_mp
      getDehomPoint_mp(dehom_mp, ED->EqD->stageData_mp[stage].endPts[path_num], ED->EqD->stageData_mp[stage].endPts[path_num]->size, &ED->preProcData);

      // print the footer to OUT for the point
      ED->EqD->stageData_mp[stage].endPt_retVals[path_num] = printStageFooter_mp(ED, stage, path_num, &recvPts[i].PD_mp, ED->EqD->stageData_mp[stage].endPts[path_num], orig_last, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      ED->EqD->curr_path_num = path_num = recvPts[i].pathNum;
      fprintf(OUT, "Path number: %d (ID: %d)\n", path_num, recvProc);

      // store the condition number
      ED->EqD->stageData_mp[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // store if higher dimenaional
      ED->EqD->stageData_mp[stage].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &recvPts[i].PD_d, &recvPts[i].PD_mp, recvPts[i].prec, recvPts[i].last_approx_d, recvPts[i].last_approx_mp, recvPts[i].last_approx_prec, NULL, ED, stage, 1);

      // copy over to the appropriate spot
      set_mp(ED->EqD->stageData_mp[stage].finalTs[path_num], recvPts[i].PD_mp.time);

      if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
      { // store the intrinsic endpoint and its cooresponding point in the original variables
        intrinsicToExtrinsic_mp(orig_last, recvPts[i].last_approx_mp, ED->EqD->stageData_mp[stage].B0, ED->EqD->stageData_mp[stage].p0);
        orig_last->size /= 2;
        intrinsicToExtrinsic_mp(dehom_mp, recvPts[i].PD_mp.point, ED->EqD->stageData_mp[stage].B0, ED->EqD->stageData_mp[stage].p0);
        dehom_mp->size /= 2; // remove the bottom half of the coordinates

        point_cp_mp(ED->EqD->stageData_mp[stage].endPts[path_num], dehom_mp);

        // convert to intrinsic coordinates
        mat_mp B_transpose;
        init_mat_mp(B_transpose, ED->EqD->stageData_mp[stage].B->cols, ED->EqD->stageData_mp[stage].B->rows);
        transpose_mp(B_transpose, ED->EqD->stageData_mp[stage].B);

        extrinsicToIntrinsic_mp(ED->EqD->stageData_mp[stage].endPts_in[path_num], dehom_mp, B_transpose, ED->EqD->stageData_mp[stage].p);

        clear_mat_mp(B_transpose);
      }
      else
      { // use the top set of coordinates so that we can store the original coordinates
        point_cp_mp(orig_last, recvPts[i].last_approx_mp);
        orig_last->size /= 2;      
        point_cp_mp(ED->EqD->stageData_mp[stage].endPts[path_num], recvPts[i].PD_mp.point);
        ED->EqD->stageData_mp[stage].endPts[path_num]->size /= 2;
      }

      // find dehom_mp
      getDehomPoint_mp(dehom_mp, ED->EqD->stageData_mp[stage].endPts[path_num], ED->EqD->stageData_mp[stage].endPts[path_num]->size, &ED->preProcData);

      // print the footer to OUT for the point
      ED->EqD->stageData_mp[stage].endPt_retVals[path_num] = printStageFooter_mp(ED, stage, path_num, &recvPts[i].PD_mp, ED->EqD->stageData_mp[stage].endPts[path_num], orig_last, dehom_mp, recvPts[i].condition_number, recvPts[i].first_increase, recvPts[i].function_residual_mp, recvPts[i].latest_newton_residual_mp, recvPts[i].t_val_at_latest_sample_point_mp, recvPts[i].error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, recvPts[i].retVal, T, trackCount);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last);
  free(recvPts);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);

  return;
}

void head_eqbyeqSortStageEndpoints_mp(basic_eval_data_mp *BED, tracker_config_t *T, int stage, double final_tol, int pathMod, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at for the current stage     *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, indexI = 0, indexJ = 0, count, cont, path_num, recvProc, numRecvPts = 0, minPacketSize = 1, maxPacketSize = 50, maxSize, num_paths = BED->EqD->stageData_mp[stage].num_paths;
  int num_stages = BED->EqD->num_subsystems, num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int)), *higherDim = NULL;
  mpf_t norm_diff;
  vec_mp tempVec;
  sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  mpf_init(norm_diff);
  init_vec_mp(tempVec, 0);

  // print header for level
  printf("\nSorting points for stage %d of %d: %d path%s to sort.\n", stage, num_stages, num_paths, num_paths == 1 ? "" : "s");

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // randomize the list
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to j
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)   
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  } 
  higherDim = (int *)bmalloc(maxSize * sizeof(int));
  for (i = 0; i < maxSize; i++)
    higherDim[i] = 0;
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("Sorting %d of %d\n", count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        sendPts[j].prec = T->Precision;
        sendPts[j].pathNum = path_num; // setup path_num
        sendPts[j].retVal = BED->EqD->stageData_mp[stage].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
        if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
        { // sort using the intrinsic coordinates
          point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->stageData_mp[stage].endPts_in[path_num]);
        }
        else
        { // sort using original coordinates
          point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->stageData_mp[stage].endPts[path_num]);
        }
        set_zero_mp(sendPts[j].PD_mp.time);
        sendPts[j].last_approx_prec = T->Precision;
        sendPts[j].last_approx_mp->size = 0;
        higherDim[j] = BED->EqD->stageData_mp[stage].higherDim[path_num];
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      MPI_Send(higherDim, lastSize[i], MPI_INT, i, TAG_INT, MPI_COMM_WORLD);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be sorted
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0 && !((count + j) % pathMod))
        printf("Sorting %d of %d\n", count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      sendPts[j].prec = T->Precision;
      sendPts[j].pathNum = path_num; // setup path_num
      sendPts[j].retVal = BED->EqD->stageData_mp[stage].endPt_retVals[path_num]; // send retVal - it will return type (except for comparison)
      if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
      { // sort using the intrinsic coordinates
        point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->stageData_mp[stage].endPts_in[path_num]);
      }
      else
      { // sort using original coordinates
        point_cp_mp(sendPts[j].PD_mp.point, BED->EqD->stageData_mp[stage].endPts[path_num]);
      }
      set_zero_mp(sendPts[j].PD_mp.time);
      sendPts[j].last_approx_prec = T->Precision;
      sendPts[j].last_approx_mp->size = 0;
      higherDim[j] = BED->EqD->stageData_mp[stage].higherDim[path_num];
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    MPI_Send(higherDim, lastSize[recvProc], MPI_INT, recvProc, TAG_INT, MPI_COMM_WORLD);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->stageData_mp[stage].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->stageData_mp[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // copy over to the appropriate spot
      if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
      { // store to the intrinsic coordinates
        point_cp_mp(BED->EqD->stageData_mp[stage].endPts_in[path_num], recvPts[i].PD_mp.point);
      }
      else
      { // store using original coordinates
        point_cp_mp(BED->EqD->stageData_mp[stage].endPts[path_num], recvPts[i].PD_mp.point);
      }
      set_mp(BED->EqD->stageData_mp[stage].finalTs[path_num], recvPts[i].PD_mp.time);

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      mpf_init(sortPts[path_num].norm);
      if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
        infNormVec_mp2(sortPts[path_num].norm, BED->EqD->stageData_mp[stage].endPts_in[path_num]);
      else
        infNormVec_mp2(sortPts[path_num].norm, BED->EqD->stageData_mp[stage].endPts[path_num]);
    }
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    { // find the path number
      path_num = recvPts[i].pathNum;
      BED->EqD->stageData_mp[stage].endPt_types[path_num] = recvPts[i].retVal;
      BED->EqD->stageData_mp[stage].condition_nums[path_num] = recvPts[i].condition_number;

      // copy over to the appropriate spot
      if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
      { // store to the intrinsic coordinates
        point_cp_mp(BED->EqD->stageData_mp[stage].endPts_in[path_num], recvPts[i].PD_mp.point);
      }
      else
      { // store using original coordinates
        point_cp_mp(BED->EqD->stageData_mp[stage].endPts[path_num], recvPts[i].PD_mp.point);
      }
      set_mp(BED->EqD->stageData_mp[stage].finalTs[path_num], recvPts[i].PD_mp.time);

      // setup sortPts[path_num]
      sortPts[path_num].path_num = path_num;
      mpf_init(sortPts[path_num].norm);
      if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
        infNormVec_mp2(sortPts[path_num].norm, BED->EqD->stageData_mp[stage].endPts_in[path_num]);
      else
        infNormVec_mp2(sortPts[path_num].norm, BED->EqD->stageData_mp[stage].endPts[path_num]);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // sort sortPts
  qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

  // do the final classification
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (BED->EqD->stageData_mp[stage].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      cont = 1;
      j = i;
      do
      { // increment the counter (start at i+1)
        j++;

        if (j < num_paths)
        { // check norm diff
          indexJ = sortPts[j].path_num;
          mpf_sub(norm_diff, sortPts[j].norm, sortPts[i].norm);
          if (mpf_get_d(norm_diff) > final_tol)
            cont = 0;
        }
        else
          cont = 0;

        // check to see if we can continue
        if (cont && BED->EqD->stageData_mp[stage].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          if (BED->EqD->stageData_mp[stage].useIntrinsicSlice)
          { // subtract using intrinsic coordinates
            vec_sub_mp(tempVec, BED->EqD->stageData_mp[stage].endPts[indexI], BED->EqD->stageData_mp[stage].endPts[indexJ]);
          }
          else
          { // subtract using original coordinates
            vec_sub_mp(tempVec, BED->EqD->stageData_mp[stage].endPts_in[indexI], BED->EqD->stageData_mp[stage].endPts_in[indexJ]);
          }

          if (infNormVec_mp(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            BED->EqD->stageData_mp[stage].endPt_types[indexI] = BED->EqD->stageData_mp[stage].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
      } while (cont);
    }

    // add to count
    if (BED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_going_to_infinity || BED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (BED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (BED->EqD->stageData_mp[stage].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (BED->EqD->stageData_mp[stage].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  BED->EqD->stageData_mp[stage].num_sing = num_sing;
  BED->EqD->stageData_mp[stage].num_nonsing = num_nonsing;
  BED->EqD->stageData_mp[stage].num_higher_dim = num_higher_dim;
  BED->EqD->stageData_mp[stage].num_bad = num_bad;
  BED->EqD->stageData_mp[stage].num_inf = num_inf;

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
    clear_endgame_data(&recvPts[i]);
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  for (i = num_paths - 1; i >= 0; i--)
    mpf_clear(sortPts[i].norm);
  free(recvPts);
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  free(sortPts);
  free(higherDim);

  mpf_clear(norm_diff);
  clear_vec_mp(tempVec);

  return;
}

////////////////////////// WORKER FUNCTIONS /////////////////////////////

void worker_eqbyeq_mp(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does equation-by-equation - 'worker' process           *
\***************************************************************/
{
  int i, size, num_subs;
  trackingStats trackCount;
  tracker_config_t T;
  basic_eval_data_mp BED;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // recv T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision
  initMP(T.Precision);

  // recv BED
  bcast_basic_eval_data_mp(&BED, 1, my_id, headnode);

  // allocate for eq-by-eq
  BED.EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));

  // recv EqD
  bcast_eqData_t(BED.EqD, T.MPType, my_id, headnode);

  // store the number of subsystems
  num_subs = BED.EqD->num_subsystems;

 // allocate memory for subsystems
  BED.EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(num_subs * sizeof(eqWitnessData_mp));
  BED.EqD->stageData_mp = (eqStageData_mp *)bmalloc(num_subs * sizeof(eqStageData_mp));

  // recv witness information for all of the subsystems
  for (i = 0; i < num_subs; i++)
  { // recv the ith witness data from the head
    bcast_witnessData(BED.EqD, i, T.MPType, my_id, headnode);
  }

  // setup the local files - OUT, MIDOUT, RAWOUT & FAIL
  size = 1 + snprintf(NULL, 0, "output_%d", my_id);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "output_%d", my_id);
  OUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "midout_w%d_%d", 0, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "midout_w%d_%d", 0, my_id);
  MIDOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  RAWOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  FAIL = fopen(str, "w");

  // loop to find the witness points for each subsystem
  for (i = 0; i < num_subs; i++)
  { // track this level
    worker_eqbyeqWitnessTrack_mp(&trackCount, i, &T, &BED, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    if (num_subs > 1)
    { // close MIDOUT
      fclose(MIDOUT);
      MIDOUT = NULL;
      // wait until all workers have closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel
    }

    // sort the endpoints
    worker_eqbyeqSortWitnessEndpoints_mp(&BED, &T, OUT, i, my_id, headnode, num_processes);

    // if this is not the last subsystem, repoen MIDOUT
    if (i + 1 < num_subs)
    { // reopen MIDOUT
      size = 1 + snprintf(NULL, 0, "midout_w%d_%d", i + 1, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_w%d_%d", i + 1, my_id);
      MIDOUT = fopen(str, "w");
    }
  }

  // clear the first witness data
  clearEqbyEqWitnessData_mp(BED.EqD, 0);

  // loop to track the stages
  for (i = 1; i < num_subs; i++)
  { // recv the ith stage from the head
    bcast_stageData(BED.EqD, i, T.MPType, my_id, headnode);

    if (i > 1)
    { // clear stage data from stage 'i - 1' since it is no longer needed
      clearEqbyEqStageData_mp(BED.EqD, i - 1);
    }
  
    // clear the witness data from subsystem 'i' since it is no longer needed
    clearEqbyEqWitnessData_mp(BED.EqD, i);

    // open MIDOUT
    size = 1 + snprintf(NULL, 0, "midout_s%d_%d", i, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "midout_s%d_%d", i, my_id);
    MIDOUT = fopen(str, "w");

    // track this stage
    worker_eqbyeqStageTrack_mp(&trackCount, i, &T, &BED, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    if (i + 1 < num_subs)
    { // close MIDOUT
      fclose(MIDOUT);
      MIDOUT = NULL;
      // wait until all workers have closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel
    }

    // sort the endpoints
    worker_eqbyeqSortStageEndpoints_mp(&BED, &T, OUT, i, my_id, headnode, num_processes);
  }

  // close the files
  fclose(OUT);
  fclose(MIDOUT);
  fclose(RAWOUT);
  fclose(FAIL);

  // clear other allocated memory
  basic_eval_clear_mp(&BED, -59, 1); // -59 since this does use equation-by-equation
  tracker_config_clear(&T);
  clearMP();
  free(str);

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

void worker_eqbyeqWitnessTrack_mp(trackingStats *trackCount, int curr_sub, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this subsystem in parallel *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0, num_subs = ED->EqD->num_subsystems;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  // setup the function evaluators
  eval_func = &witness_eqbyeq_eval_mp;
  find_dehom = eqbyeq_witness_dehom;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Finding witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    // initialize for MP tracking
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], T->Precision);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED->EqD->curr_stage_num = curr_sub;

    for (i = 0; i < numStartPts; i++)
    { // setup for tracking the path
      ED->EqD->curr_path_num = startPts[i].pathNum;

      T->first_step_of_path = 1;
      T->endgameOnly = 0;

      // print the header of the path to OUT
      printPathHeader_mp(OUT, &startPts[i].PD_mp, T, startPts[i].pathNum, ED, eval_func);

      // track the path
      zero_dim_track_path_mp(startPts[i].pathNum, &endPts[i], &startPts[i].PD_mp, OUT, MIDOUT, T, ED, eval_func, find_dehom);

      // check to see if it should be sharpened - only when this is the only subsystems
      if (num_subs == 1 && endPts[i].retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], T, OUT, NULL, ED, NULL, eval_func, NULL);
      }

      // print footer
      printResultOfPath(OUT, endPts[i].retVal, T);
      printBasicFooter_mp(OUT, &endPts[i].PD_mp, T, endPts[i].function_residual_mp);
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      // initialize for MP tracking
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqSortWitnessEndpoints_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, int curr_sub, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points for the current subsystem     *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, rankDef, finite, numStartPts = 0, numEndPts = 0, *higherDim = NULL;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  MPI_Status status;

  // setup the function evaluators
  eval_func = &witness_eqbyeq_eval_mp;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting witness points for subsystem %d.\n", curr_sub);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
  // recv higherDim
  higherDim = (int *)bmalloc(numStartPts * sizeof(int));
  MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    // initialize for MP tracking
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], T->Precision);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED->EqD->curr_stage_num = curr_sub;

    for (i = 0; i < numStartPts; i++)
    { // setup for evaluation
      ED->EqD->curr_path_num = startPts[i].pathNum;

      // setup endPts
      endPts[i].prec = T->Precision;
      endPts[i].pathNum = startPts[i].pathNum;
      point_cp_mp(endPts[i].PD_mp.point, startPts[i].PD_mp.point);
      set_mp(endPts[i].PD_mp.time, startPts[i].PD_mp.time);
      endPts[i].last_approx_prec = T->Precision;
      endPts[i].last_approx_mp->size = 0;

      if (startPts[i].retVal == 0)
      { // set time to 0
        set_zero_mp(endPts[i].PD_mp.time);

        // determine if it is rank deficient, finite, higherDim
        eqbyeqWitnessSortEndpoint_mp(&rankDef, &finite, higherDim[i], ED, curr_sub, startPts[i].pathNum, &endPts[i].condition_number, T, OUT, &endPts[i].PD_mp, eval_func);

        // check for success on the jacobian analyzer
        if (T->regen_remove_inf && !finite)
        { // dehom point is infinite
          endPts[i].retVal = retVal_going_to_infinity;
        }
        else if (T->regen_higher_dim_check && higherDim[i])
        { // it lies on a higher dimensional component
          endPts[i].retVal = retVal_higher_dim;
        }
        else if (!rankDef)
        { // classify as non-singular
          endPts[i].retVal = MOVE_TO_NEXT;
        }
        else
        { // classify as singular
          endPts[i].retVal = DO_NOT_MOVE_TO_NEXT;
        }
      }
      else
      { // path was not a success - copy over error code
        endPts[i].retVal = startPts[i].retVal;
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
    // recv higherDim
    higherDim = (int *)brealloc(higherDim, numStartPts * sizeof(int));
    if (numStartPts > 0)
      MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      // initialize for MP tracking
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // startPts & endPts & higherDim are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqStageTrack_mp(trackingStats *trackCount, int curr_stage, tracker_config_t *T, basic_eval_data_mp *ED, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this stage in parallel     *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0, num_stages = ED->EqD->num_subsystems;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  // setup the function evaluators
  eval_func = &standard_eqbyeq_eval_mp;
  find_dehom = &eqbyeq_stage_dehom;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking points for stage %d.\n", curr_stage);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    // initialize for MP tracking
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], T->Precision);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED->EqD->curr_stage_num = curr_stage;

    for (i = 0; i < numStartPts; i++)
    { // setup for tracking the path
      ED->EqD->curr_path_num = startPts[i].pathNum;

      T->first_step_of_path = 1;
      T->endgameOnly = 0;

      // print the header of the path to OUT
      printPathHeader_mp(OUT, &startPts[i].PD_mp, T, startPts[i].pathNum, ED, eval_func);

      // track the path
      zero_dim_track_path_mp(startPts[i].pathNum, &endPts[i], &startPts[i].PD_mp, OUT, MIDOUT, T, ED, eval_func, find_dehom);

      // check to see if it should be sharpened - only when this is the only subsystems
      if (curr_stage + 1 == num_stages && endPts[i].retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], T, OUT, NULL, ED, NULL, eval_func, NULL);
      }

      // print footer
      printResultOfPath(OUT, endPts[i].retVal, T);
      printBasicFooter_mp(OUT, &endPts[i].PD_mp, T, endPts[i].function_residual_mp);
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      // initialize for MP tracking
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_eqbyeqSortStageEndpoints_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, int curr_stage, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the points for the current stage                 *
*  Uses parallel to refine and check for singular endpoints     *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, rankDef, finite, numStartPts = 0, numEndPts = 0, *higherDim = NULL;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  MPI_Status status;

  // setup the function evaluators
  eval_func = &stage_sort_eqbyeq_eval_mp;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting points for stage %d.\n", curr_stage);
  fprintf(OUT, "*****************************************************\n");

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
  // recv higherDim
  higherDim = (int *)bmalloc(numStartPts * sizeof(int));
  MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    // initialize for MP tracking
    for (i = 0; i < numEndPts; i++)
      init_endgame_data(&endPts[i], T->Precision);
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    ED->EqD->curr_stage_num = curr_stage;

    for (i = 0; i < numStartPts; i++)
    { // setup for evaluation
      ED->EqD->curr_path_num = startPts[i].pathNum;

      // setup endPts
      endPts[i].prec = T->Precision;
      endPts[i].pathNum = startPts[i].pathNum;
      point_cp_mp(endPts[i].PD_mp.point, startPts[i].PD_mp.point);
      set_mp(endPts[i].PD_mp.time, startPts[i].PD_mp.time);
      endPts[i].last_approx_prec = T->Precision;
      endPts[i].last_approx_mp->size = 0;

      if (startPts[i].retVal == 0)
      { // set time to 0
        set_zero_mp(endPts[i].PD_mp.time);

        // determine if it is rank deficient, finite, higherDim
        eqbyeqStageSortEndpoint_mp(&rankDef, &finite, higherDim[i], ED, curr_stage, startPts[i].pathNum, &endPts[i].condition_number, T, OUT, &endPts[i].PD_mp, eval_func);

        // check for success on the jacobian analyzer
        if (T->regen_remove_inf && !finite)
        { // dehom point is infinite
          endPts[i].retVal = retVal_going_to_infinity;
        }
        else if (T->regen_higher_dim_check && higherDim[i])
        { // it lies on a higher dimensional component
          endPts[i].retVal = retVal_higher_dim;
        }
        else if (!rankDef)
        { // classify as non-singular
          endPts[i].retVal = MOVE_TO_NEXT;
        }
        else
        { // classify as singular
          endPts[i].retVal = DO_NOT_MOVE_TO_NEXT;
        }
      }
      else
      { // path was not a success - copy over error code
        endPts[i].retVal = startPts[i].retVal;
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);
    // recv higherDim
    higherDim = (int *)brealloc(higherDim, numStartPts * sizeof(int));
    if (numStartPts > 0)
      MPI_Recv(higherDim, numStartPts, MPI_INT, headnode, TAG_INT, MPI_COMM_WORLD, &status);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
        clear_endgame_data(&endPts[i]);

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      // initialize for MP tracking
      for (i = 0; i < numEndPts; i++)
        init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

#endif

