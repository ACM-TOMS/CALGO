// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

// parallel endgame (either using the track back feature or not)

#ifdef _HAVE_MPI

int create_send_packet(trackingStats *trackCount, int pathMod, int *currLoc, int totalPaths, int proc, int packetSize, endgame_data_t *sendPts, tracker_config_t *T, int trackBack_count, trackBack_samples_t *trackBack_points, point_d *startPts_d, double *startPoint_norm_d, point_mp *startPts_mp, mpf_t *startPoint_norm_mp, double trackBack_final_tol, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, void const *ED_d, void const *ED_mp, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: create and send a packet to proc                       *
\***************************************************************/
{
  int i, j, retVal, indexI, indexJ, size, old_path_num, count = 0;
  int *goodI = NULL, *goodJ = NULL, *goodPathNum = NULL, goodCount = 0;

  for (i = 0; i < packetSize; i++)
  { // setup the next path
    do
    { // print the path number - if needed
      if (pathMod > 0)
        if (!(*currLoc % pathMod))
          printf("Tracking path %d of %d\n", *currLoc, totalPaths);

      // check the next start point
      if (T->MPType == 0 || T->MPType == 2)
      { // check in double precision
        retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &trackBack_points, trackBack_count, startPoint_norm_d[*currLoc], NULL, startPts_d[*currLoc], NULL, 52);
      }
      else
      { // check using fixed precision
        retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &trackBack_points, trackBack_count, 0, startPoint_norm_mp[*currLoc], NULL, startPts_mp[*currLoc], T->Precision);
      }

      if (retVal == 0)
      { // this is a new start point - setup sendPts
        if (T->MPType == 0 || T->MPType == 2)
        { // send out in double precision
          sendPts[count].prec = 52; // start point is in double precision
          sendPts[count].pathNum = *currLoc;
          point_cp_d(sendPts[count].PD_d.point, startPts_d[*currLoc]);
          set_one_d(sendPts[count].PD_d.time);

          sendPts[count].last_approx_prec = 52;
          sendPts[count].last_approx_d->size = 0;
        }
        else
        { // send out in MP
          sendPts[count].prec = T->Precision; // start point is in the fixed precision
          sendPts[count].pathNum = *currLoc;
          point_cp_mp(sendPts[count].PD_mp.point, startPts_mp[*currLoc]);
          set_one_mp(sendPts[count].PD_mp.time);

          sendPts[count].last_approx_prec = T->Precision;
          sendPts[count].last_approx_mp->size = 0;
        }

        // increment count
        count++;
      }
      else // this point is already known
      { // save indexI & indexJ
        goodI = (int *)brealloc(goodI, (goodCount + 1) * sizeof(int));
        goodJ = (int *)brealloc(goodJ, (goodCount + 1) * sizeof(int));
        goodPathNum = (int *)brealloc(goodPathNum, (goodCount + 1) * sizeof(int));
        goodI[goodCount] = indexI;
        goodJ[goodCount] = indexJ;
        goodPathNum[goodCount] = *currLoc;
        goodCount++;

        // say that this sample has been used
        trackBack_points[indexI].samplePts_notUsed[indexJ] = 0;
      }

      // increment currLoc
      (*currLoc)++;

      // make sure we are not at the end
      if (*currLoc >= totalPaths)
      { // break out of the loops
        i = packetSize;
        break;
      }
      else if (retVal == 0)
      { // this one has been copied so we can break the first loop
        break;
      }
    } while (1);
  }

  // send sendPts to proc
  send_recv_endgame_data_t(&sendPts, &count, T->MPType, proc, 1);

  // print the points that were found
  for (i = 0; i < goodCount; i++)
  { // print the data to OUT & RAWOUT
    fprintf(OUT, "Path number: %d (ID: %d)\n", goodPathNum[i], headnode);

    // save and update the path number so that it prints properly
    old_path_num = trackBack_points[goodI[i]].endPt.pathNum;
    trackBack_points[goodI[i]].endPt.pathNum = goodPathNum[i];
    if (T->MPType == 0 || (T->MPType == 2 && trackBack_points[goodI[i]].endPt.prec < 64))
    { // print footer in double precision
      printPathFooter_d(trackCount, &trackBack_points[goodI[i]].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
    }
    else
    { // print footer in multi precision
      printPathFooter_mp(trackCount, &trackBack_points[goodI[i]].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
    }
    // set the path number back
    trackBack_points[goodI[i]].endPt.pathNum = old_path_num;

    // print to MIDOUT
    fprintf(MIDOUT, "%d\n", goodPathNum[i]); 
    if (T->MPType == 0 || (T->MPType == 2 && trackBack_points[goodI[i]].midPt_prec < 64))
    { // print using midPt_d
      size = trackBack_points[goodI[i]].midPt_d[goodJ[i]]->size;
      for (j = 0; j < size; j++)
        fprintf(MIDOUT, "%.15e %.15e\n", trackBack_points[goodI[i]].midPt_d[goodJ[i]]->coord[j].r, trackBack_points[goodI[i]].midPt_d[goodJ[i]]->coord[j].i);
      fprintf(MIDOUT, "\n");
    }
    else
    { // print using midPt_mp
      size = trackBack_points[goodI[i]].midPt_mp[goodJ[i]]->size;
      for (j = 0; j < size; j++)
      {
        print_mp(MIDOUT, 0, &trackBack_points[goodI[i]].midPt_mp[goodJ[i]]->coord[j]);
        fprintf(MIDOUT, "\n");
      }
      fprintf(MIDOUT, "\n");
    }
  }

  // clear memory
  free(goodI);
  free(goodJ);
  free(goodPathNum);

  // return count
  return count;
}

void head_zero_dim_endgame(int useTrackBack, int minCycleTrackBack, trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'headnode'            *
\***************************************************************/
{
  int i, count, recvProc, maxSize, numPoints, numRecvPts = 0, maxPacketSize = 20, minPacketSize = 1;
  int *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  double trackBack_final_tol = 0, minDiff_d = 0, *startPts_norm_d = NULL;
  mpf_t minDiff_mp, *startPts_norm_mp = NULL;
  point_d *startPts_d = NULL;
  point_mp *startPts_mp = NULL;
  endgame_data_t *sendPts = NULL, *recvPts = NULL;
  char ch, *str = NULL;
  FILE *tempFILE = NULL, *NONSOLN = NULL;

  // setup NONSOLN
  if ((T->MPType == 1 && !ED_mp->squareSystem.noChanges) || ((T->MPType == 0 || T->MPType == 2) && !ED_d->squareSystem.noChanges))
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  int trackBack_count = 0;
  trackBack_samples_t *recvPts_trackBack = NULL, *trackBack_points = NULL;

  // initialize trackCount
  init_trackingStats(trackCount);

  // initialize minDiff_mp
  mpf_init(minDiff_mp);  

  // setup data about start points
  if (useTrackBack)
    find_minimum_separation(&numPoints, &minDiff_d, &startPts_norm_d, &startPts_d, minDiff_mp, &startPts_norm_mp, &startPts_mp, T, START);

  // record the number of start points
  trackCount->numPoints = numPoints;

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  if (T->MPType == 0)
  { // send ED_d to the workers - only using double precision
    bcast_basic_eval_data_d(ED_d, T->MPType, my_id, headnode);
  }
  else if (T->MPType == 1)
  { // send ED_mp to the workers
    bcast_basic_eval_data_mp(ED_mp, T->MPType, my_id, headnode);
  }
  else 
  { // send ED_d to the workers - using AMP
    bcast_basic_eval_data_amp(ED_d, my_id, headnode);
  }

  // send whether to use track back
  MPI_Bcast(&useTrackBack, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (useTrackBack)
  { // setup trackBack_final_tol
    if (T->MPType == 1)
      trackBack_final_tol = mpf_get_d(minDiff_mp) * 1e-2;
    else
      trackBack_final_tol = minDiff_d * 1e-2;

    // send the tolerance
    MPI_Bcast(&trackBack_final_tol, 1, MPI_DOUBLE, headnode, MPI_COMM_WORLD);
    MPI_Bcast(&minCycleTrackBack, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  }

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // setup packetSizes & determine maxSize  
  packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup sendPts based on maxSize
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
    init_endgame_data(&sendPts[i], 64);

  // send out the initial set of packets
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create & send packet to proc i
      lastSize[i] = create_send_packet(trackCount, pathMod, &count, trackCount->numPoints, i, packetSizes[i], sendPts, T, trackBack_count, trackBack_points, startPts_d, startPts_norm_d, startPts_mp, startPts_norm_mp, trackBack_final_tol, OUT, RAWOUT, MIDOUT, FAIL, NONSOLN, ED_d, ED_mp, headnode);
    }
  
  // loop until all of the paths have been sent out to the workers
  while (count < trackCount->numPoints)
  { // recv a packet back & who sent it
    if (useTrackBack)
    { // recv the track back points
      recvProc = send_recv_trackBack_samples_t(&recvPts_trackBack, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

      // update trackBack_points with the ones that have sample points
      for (i = 0; i < numRecvPts; i++)
        if (recvPts_trackBack[i].numSamples > 0)
        {
          trackBack_points = (trackBack_samples_t *)brealloc(trackBack_points, (trackBack_count + 1) * sizeof(trackBack_samples_t));
          cp_trackBack_samples_t(&trackBack_points[trackBack_count], &recvPts_trackBack[i], 1);
          trackBack_count++;
        }
    }  
    else
    { // recv the points back
      recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);
    }

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create & send the next packet to recvProc
    lastSize[recvProc] = create_send_packet(trackCount, pathMod, &count, trackCount->numPoints, recvProc, packetSizes[recvProc], sendPts, T, trackBack_count, trackBack_points, startPts_d, startPts_norm_d, startPts_mp, startPts_norm_mp, trackBack_final_tol, OUT, RAWOUT, MIDOUT, FAIL, NONSOLN, ED_d, ED_mp, headnode);

    // print the footers for each of the paths recvd
    if (useTrackBack)
    { 
      // print the footers 
      for (i = 0; i < numRecvPts; i++)
      {
        fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts_trackBack[i].endPt.pathNum, recvProc);
        if (T->MPType == 0 || (T->MPType == 2 && recvPts_trackBack[i].endPt.prec < 64))
        { // print footer in double precision
          printPathFooter_d(trackCount, &recvPts_trackBack[i].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
        }
        else
        { // print footer in multi precision
          printPathFooter_mp(trackCount, &recvPts_trackBack[i].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
        }
      }
    }
    else
    {
      for (i = 0; i < numRecvPts; i++)
      {
        fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
        if (T->MPType == 0 || (T->MPType == 2 && recvPts[i].prec < 64))
        { // print footer in double precision
          printPathFooter_d(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
        }
        else
        { // print footer in multi precision
          printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
        }
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
    if (useTrackBack)
    { // recv the track back points
      recvProc = send_recv_trackBack_samples_t(&recvPts_trackBack, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

      // print the footers for each of the paths
      for (i = 0; i < numRecvPts; i++)
      {
        fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts_trackBack[i].endPt.pathNum, recvProc);
        if (T->MPType == 0 || (T->MPType == 2 && recvPts_trackBack[i].endPt.prec < 64))
        { // print footer in double precision
          printPathFooter_d(trackCount, &recvPts_trackBack[i].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
        }
        else
        { // print footer in multi precision
          printPathFooter_mp(trackCount, &recvPts_trackBack[i].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
        }
      }
    }
    else
    { // recv the points back
      recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

      // print the footers for each of the paths
      for (i = 0; i < numRecvPts; i++)
      {
        fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
        if (T->MPType == 0 || (T->MPType == 2 && recvPts[i].prec < 64))
        { // print footer in double precision
          printPathFooter_d(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
        }
        else
        { // print footer in multi precision
          printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
        }
      }
    }

    // tell the worker that tracking is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // since all of the workers have sent back their data, we just need to wait until the files are closed
  MPI_Barrier(MPI_COMM_WORLD);

  // combine files
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // open midout_'i'
      count = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "midout_%d", i);
      tempFILE = fopen(str, "r");

      // make sure the file exists
      if (tempFILE == NULL)
      {
        printf("ERROR: %s file does not exist!\n", str);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // copy it over
      ch = fgetc(tempFILE);
      while (ch != EOF)
      {
        fprintf(MIDOUT, "%c", ch);
        ch = fgetc(tempFILE);
      }

      // close the file and delete it
      fclose(tempFILE);
      remove(str);

      // delete fail_'i' & rawout_'i' & nonsolutions_'i' - the headnode made a global version of these during printPathFooter
      count = 1 + snprintf(NULL, 0, "rawout_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "rawout_%d", i);
      remove(str);

      count = 1 + snprintf(NULL, 0, "fail_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "fail_%d", i);
      remove(str);

      count = 1 + snprintf(NULL, 0, "nonsolutions_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "nonsolutions_%d", i);
      remove(str);
    }

  if ((T->MPType == 1 && !ED_mp->squareSystem.noChanges) || ((T->MPType == 0 || T->MPType == 2) && !ED_d->squareSystem.noChanges))
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // clear memory
  free(str);
  mpf_clear(minDiff_mp);  
  for (i = maxSize - 1; i >= 0; i--)
    clear_endgame_data(&sendPts[i]);
  free(sendPts);
  if (useTrackBack)
  {
    for (i = numRecvPts - 1; i >= 0; i--)
      clear_trackBack_sample(&recvPts_trackBack[i]);
    free(recvPts_trackBack);
  }
  else
  {
    for (i = numRecvPts - 1; i >= 0; i--)
      clear_endgame_data(&recvPts[i]);
    free(recvPts);
  }
  free(packetSizes);
  free(lastSize);
  for (i = trackBack_count - 1; i >= 0; i--)
    clear_trackBack_sample(&trackBack_points[i]);
  free(trackBack_points);
  if (T->MPType == 0 || T->MPType == 2)
  {
    for (i = numPoints - 1; i >= 0; i--)
      clear_point_d(startPts_d[i]);
    free(startPts_d);
    free(startPts_norm_d);
  }
  else
  {
    for (i = numPoints - 1; i >= 0; i--)
    {
      clear_point_mp(startPts_mp[i]);
      mpf_clear(startPts_norm_mp[i]);
    }
    free(startPts_mp);
    free(startPts_norm_mp);
  }

  return;
}

void worker_zero_dim_endgame(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'worker' process      *
\***************************************************************/
{
  int i, j, rV, indexI, indexJ, size, useTrackBack = 0, numStartPts = 0, numEndPts = 0, minCycleTrackBack = 0;
  double norm_d, trackBack_final_tol = 0;
  mpf_t norm_mp;
  trackingStats trackCount;
  tracker_config_t T;
  basic_eval_data_d BED_d;
  basic_eval_data_mp BED_mp;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  trackBack_samples_t *endPts_trackBack = NULL;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL, *NONSOLN = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  mpf_init(norm_mp);

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // setup eval_func_d & eval_func_mp
  change_prec = &change_basic_eval_prec;
  find_dehom = &zero_dim_dehom;
  if (dataType == USER_HOM_D || dataType == USER_HOM_MP)
  { // user defined homotopy
    eval_func_d = &userHom_eval_d;
    eval_func_mp = &userHom_eval_mp;
  }
  else if (dataType == PARAM_HOM_D || dataType == PARAM_HOM_MP)
  { // parameter homotopy
    eval_func_d = &paramHom_eval_d;
    eval_func_mp = &paramHom_eval_mp;
  }
  else
  { // standard zero dim
    eval_func_d = &basic_eval_d;
    eval_func_mp = &basic_eval_mp;
  }

  // setup T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision - does not hurt if doing only double precision
  initMP(T.Precision);

  if (T.MPType == 0)
  { // setup BED only for double precision
    bcast_basic_eval_data_d(&BED_d, T.MPType, my_id, headnode);
  }
  else if (T.MPType == 1)
  { // setup BED only for fixed multi precision
    bcast_basic_eval_data_mp(&BED_mp, T.MPType, my_id, headnode);
  }
  else
  { // setup BED to use AMP
    bcast_basic_eval_data_amp(&BED_d, my_id, headnode);
  }

  // recv whether to track back
  MPI_Bcast(&useTrackBack, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (useTrackBack)
  { // recv the tolerance
    MPI_Bcast(&trackBack_final_tol, 1, MPI_DOUBLE, headnode, MPI_COMM_WORLD);
    // recv the minimum cycle for track back
    MPI_Bcast(&minCycleTrackBack, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  }

  // setup the local files - OUT, MIDOUT, RAWOUT & FAIL
  size = 1 + snprintf(NULL, 0, "output_%d", my_id);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "output_%d", my_id);
  OUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "midout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "midout_%d", my_id);
  MIDOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  RAWOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  FAIL = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "nonsolutions_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "nonsolutions_%d", my_id);
  NONSOLN = fopen(str, "w");

  free(str);
  str = NULL;

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T.MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    endPts_trackBack = (trackBack_samples_t *)bmalloc(numEndPts * sizeof(trackBack_samples_t));
    for (i = 0; i < numEndPts; i++)
    { // initialize
      init_endgame_data(&endPts[i], T.Precision);
      init_trackBack_sample(&endPts_trackBack[i], T.Precision);
    } 
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd
    for (i = 0; i < numStartPts; i++)
    { 
      if (T.MPType == 1)
      { // print the header of the path to OUT
        printPathHeader_mp(OUT, &startPts[i].PD_mp, &T, startPts[i].pathNum, &BED_mp, eval_func_mp);

        // track the path
        if (useTrackBack)
        { // use the track back tracker

          // check to see if we have computed this one locally
          infNormVec_mp2(norm_mp, startPts[i].PD_mp.point);
          rV = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &endPts_trackBack, i, 0, norm_mp, NULL, startPts[i].PD_mp.point, T.Precision);

          if (rV == 0)
          { // this is not known locally - track the path
            zero_dim_trackBack_path_mp(startPts[i].pathNum, trackBack_final_tol, &endPts_trackBack[i], &startPts[i].PD_mp, OUT, MIDOUT, &T, &BED_mp, eval_func_mp, find_dehom);

            // check to see if it should be sharpened
            if (endPts_trackBack[i].endPt.retVal == 0 && T.sharpenDigits > 0)
            { // use the sharpener for after an endgame
              sharpen_endpoint_endgame(&endPts_trackBack[i].endPt, &T, OUT, NULL, &BED_mp, NULL, eval_func_mp, NULL);
            }
          }
          else
          { // this point is already known locally
            // mark that this one has been used
            endPts_trackBack[indexI].samplePts_notUsed[indexJ] = 0;

            // print using midPt_mp
            fprintf(MIDOUT, "%d\n", startPts[i].pathNum);
            size = endPts_trackBack[indexI].midPt_mp[indexJ]->size;
            for (j = 0; j < size; j++)
            {
              print_mp(MIDOUT, 0, &endPts_trackBack[indexI].midPt_mp[indexJ]->coord[j]);
              fprintf(MIDOUT, "\n");
            }
            fprintf(MIDOUT, "\n");

            // copy the endpoint
            cp_endgame_data_t(&endPts_trackBack[i].endPt, &endPts_trackBack[indexI].endPt, 0);
            // null out other structures
            endPts_trackBack[i].numSamples = endPts_trackBack[i].samplePts_prec = endPts_trackBack[i].midPt_prec = 0;
            endPts_trackBack[i].normSamples_d = NULL;
            endPts_trackBack[i].normSamples_mp = NULL;
            endPts_trackBack[i].samplePts_d = endPts_trackBack[i].midPt_d = NULL;
            endPts_trackBack[i].samplePts_mp = endPts_trackBack[i].midPt_mp = NULL;
            endPts_trackBack[i].samplePts_notUsed = NULL;
          }

          // print footer
          printPathFooter_mp(&trackCount, &endPts_trackBack[i].endPt, &T, OUT, RAWOUT, FAIL, NONSOLN, &BED_mp);
        }
        else 
        { // use the standard tracker
          zero_dim_track_path_mp(startPts[i].pathNum, &endPts[i], &startPts[i].PD_mp, OUT, MIDOUT, &T, &BED_mp, eval_func_mp, find_dehom);

          // check to see if it should be sharpened
          if (endPts[i].retVal == 0 && T.sharpenDigits > 0)
          { // use the sharpener for after an endgame
            sharpen_endpoint_endgame(&endPts[i], &T, OUT, NULL, &BED_mp, NULL, eval_func_mp, NULL);
          }
 
          // print footer
          printPathFooter_mp(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, &BED_mp);
        }
      }
      else
      { // print the header of the path to OUT
        printPathHeader_d(OUT, &startPts[i].PD_d, &T, startPts[i].pathNum, &BED_d, eval_func_d);

        // track the path
        if (useTrackBack)
        { // use the track back tracker

          // check to see if we have computed this one locally
          norm_d = infNormVec_d(startPts[i].PD_d.point);
          rV = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &endPts_trackBack, i, norm_d, NULL, startPts[i].PD_d.point, NULL, 52);  
          if (rV == 0)
          { // this is not known locally - track the path
            zero_dim_trackBack_path_d(startPts[i].pathNum, trackBack_final_tol, &endPts_trackBack[i], &startPts[i].PD_d, OUT, MIDOUT, &T, &BED_d, BED_d.BED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

            // check to see if it should be sharpened
            if (endPts_trackBack[i].endPt.retVal == 0 && T.sharpenDigits > 0)
            { // use the sharpener for after an endgame
              sharpen_endpoint_endgame(&endPts_trackBack[i].endPt, &T, OUT, &BED_d, BED_d.BED_mp, eval_func_d, eval_func_mp, change_prec);
            }
          }
          else
          { // this point is already known locally 
            // mark that this one has been used
            endPts_trackBack[indexI].samplePts_notUsed[indexJ] = 0;

            // print midpoint
            fprintf(MIDOUT, "%d\n", startPts[i].pathNum);
            if (T.MPType == 0 || (T.MPType == 2 && endPts_trackBack[indexI].midPt_prec < 64))
            { // print using midPt_d
              size = endPts_trackBack[indexI].midPt_d[indexJ]->size;
              for (j = 0; j < size; j++)
                fprintf(MIDOUT, "%.15e %.15e\n", endPts_trackBack[indexI].midPt_d[indexJ]->coord[j].r, endPts_trackBack[indexI].midPt_d[indexJ]->coord[j].i);
              fprintf(MIDOUT, "\n");
            }
            else
            { // print using midPt_mp
              size = endPts_trackBack[indexI].midPt_mp[indexJ]->size;
              for (j = 0; j < size; j++)
              {
                print_mp(MIDOUT, 0, &endPts_trackBack[indexI].midPt_mp[indexJ]->coord[j]);
                fprintf(MIDOUT, "\n");
              }
              fprintf(MIDOUT, "\n");
            }

            // copy the endpoint
            cp_endgame_data_t(&endPts_trackBack[i].endPt, &endPts_trackBack[indexI].endPt, 0);
            // null out other structures
            endPts_trackBack[i].numSamples = endPts_trackBack[i].samplePts_prec = endPts_trackBack[i].midPt_prec = 0;
            endPts_trackBack[i].normSamples_d = NULL;
            endPts_trackBack[i].normSamples_mp = NULL;
            endPts_trackBack[i].samplePts_d = endPts_trackBack[i].midPt_d = NULL;
            endPts_trackBack[i].samplePts_mp = endPts_trackBack[i].midPt_mp = NULL;
            endPts_trackBack[i].samplePts_notUsed = NULL;
          }

          // print the footer
          if (endPts_trackBack[i].endPt.prec < 64)
          { // print footer in double precision
            printPathFooter_d(&trackCount, &endPts_trackBack[i].endPt, &T, OUT, RAWOUT, FAIL, NONSOLN, &BED_d);
          }
          else
          { // print footer in multi precision
            printPathFooter_mp(&trackCount, &endPts_trackBack[i].endPt, &T, OUT, RAWOUT, FAIL, NONSOLN, BED_d.BED_mp);
          }
        }
        else
        { // use the standard tracker
          zero_dim_track_path_d(startPts[i].pathNum, &endPts[i], &startPts[i].PD_d, OUT, MIDOUT, &T, &BED_d, BED_d.BED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

          // check to see if it should be sharpened
          if (endPts[i].retVal == 0 && T.sharpenDigits > 0)
          { // use the sharpener for after an endgame
            sharpen_endpoint_endgame(&endPts[i], &T, OUT, &BED_d, BED_d.BED_mp, eval_func_d, eval_func_mp, change_prec);
          }

          if (endPts[i].prec < 64)
          { // print footer in double precision
            printPathFooter_d(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, &BED_d);
          }
          else
          { // print footer in multi precision
            printPathFooter_mp(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, BED_d.BED_mp);
          }
        }
      }
    }

    if (useTrackBack)
    { // send track back data to head
      send_recv_trackBack_samples_t(&endPts_trackBack, &numEndPts, T.MPType, headnode, 1);
    }
    else
    { // send endPts back to head
      send_recv_endgame_data_t(&endPts, &numEndPts, T.MPType, headnode, 1);
    }

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T.MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
      {
        clear_endgame_data(&endPts[i]);
        clear_trackBack_sample(&endPts_trackBack[i]);
      }

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      endPts_trackBack = (trackBack_samples_t *)brealloc(endPts_trackBack, numEndPts * sizeof(trackBack_samples_t));
      for (i = 0; i < numEndPts; i++)
      {
        init_endgame_data(&endPts[i], T.Precision);
        init_trackBack_sample(&endPts_trackBack[i], T.Precision);
      }
    }
  }

  // close the files
  fclose(OUT);
  fclose(MIDOUT);
  fclose(RAWOUT);
  fclose(FAIL);
  fclose(NONSOLN);

  // the headnode needs to wait to copy over the files until all have been closed properly
  MPI_Barrier(MPI_COMM_WORLD);

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  // clear other allocated memory
  mpf_clear(norm_mp);
  if (T.MPType == 1)
    basic_eval_clear_mp(&BED_mp, 0, 1); // 0 since not regeneration, 1 for clear the Prog
  else
    basic_eval_clear_d(&BED_d, 0, T.MPType); // 0 since this does not use regeneration
  tracker_config_clear(&T);
  clearMP();

  return;
}



#endif

