// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"
#include "regen_pos_dim.h"

// This file contains the main parallel functions

void worker_zero_dim_track_d(int my_id, int num_processes, int headnode, int dataType);
void worker_zero_dim_track_mp(int my_id, int num_processes, int headnode, int dataType);

#ifdef _HAVE_MPI

void head_zero_dim_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'headnode'            *
\***************************************************************/
{
  int i, j, recvProc, count = 0, eachStartPt, numRecvPts = 0, maxSize = 0, minPacketSize = 1, maxPacketSize = 20;
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *startPts = NULL, *sendPts = NULL, *recvPts = NULL;
  char ch, *str = NULL;
  FILE *tempFILE = NULL, *NONSOLN = NULL;

  // setup NONSOLN
  if (!ED_d->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // initialize trackCount
  init_trackingStats(trackCount);

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

  // not using trackBack endgame
  count = 0;
  MPI_Bcast(&count, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
  // Find the number of start points
  fscanf(START, "%d\n", &trackCount->numPoints);

  if (trackCount->numPoints > num_processes * maxPacketSize)
  { // since there is a large number of points, read them in one by one
    eachStartPt = 1;
  }
  else
  { // read them all in at once
    eachStartPt = 0;
  }

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup startPts - if we do not need to setup startPts, they will be read directly from file to sendPts
  if (!eachStartPt)
  { // setup pathNum - randomize the list of start points
    pathNum = (int *)bmalloc(trackCount->numPoints * sizeof(int));
    for (i = 0; i < trackCount->numPoints; i++)
    { // find random integer from 0 to i
      j = rand() % (i + 1);
      // swap i & j
      pathNum[i] = pathNum[j];
      pathNum[j] = i;
    }

    // setup startPts
    startPts = (endgame_data_t *)bmalloc(trackCount->numPoints * sizeof(endgame_data_t));
    for (i = 0; i < trackCount->numPoints; i++)
    { // initialize pathNum[i]
      init_endgame_data(&startPts[pathNum[i]], 64);

      // read in startPts[pathNum[i]]
      startPts[pathNum[i]].prec = 52; // start point is in double precision
      startPts[pathNum[i]].pathNum = i; // path number
      setupStart_d(T, &startPts[pathNum[i]].PD_d, START);
    }
  }
  else
  { // setup sendPts based on maxSize
    sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
    for (i = 0; i < maxSize; i++)
      init_endgame_data(&sendPts[i], 64);
  }

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0)
          if (!((count + j)% pathMod))
            printf("Tracking path %d of %d\n", count + j, trackCount->numPoints);

        if (eachStartPt)
        { // read in from START
          sendPts[j].prec = 52; // start point is in double precision
          sendPts[j].pathNum = count + j; // path number is count + j
          setupStart_d(T, &sendPts[j].PD_d, START);
        }
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      if (!eachStartPt)
      {
        sendPts = &startPts[count];
      }
      // send sendPts
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      if (!eachStartPt)
      { 
        sendPts = NULL;
      }

      // update count
      count += lastSize[i];
    }

  // loop until all of the paths have been sent out to the workers to be tracked
  while (count < trackCount->numPoints)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0)
        if (!((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, trackCount->numPoints);

      if (eachStartPt)
      { // read in from START
        sendPts[j].prec = 52; // start point is in double precision
        sendPts[j].pathNum = count + j; // path number is count + j
        setupStart_d(T, &sendPts[j].PD_d, START);
      }
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to recvProc
    if (!eachStartPt)
    {
      sendPts = &startPts[count];
    }
    // send sendPts
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    if (!eachStartPt)
    {
      sendPts = NULL;
    }

    // update count
    count += lastSize[recvProc];

    // print the footers for each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    {
      fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
      if (recvPts[i].prec < 64)
      { // print footer in double precision
        printPathFooter_d(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
      }
      else
      { // print footer in multi precision
        printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
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

    // print the footers for each of the paths
    for (i = 0; i < numRecvPts; i++)
    {
      fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
      if (recvPts[i].prec < 64)
      { // print footer in double precision
        printPathFooter_d(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
      }
      else
      { // print footer in multi precision
        printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
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

  // the only thing that needs combined is all of the MIDOUT files so that the midpoint checker will work
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // open midout_'i'
      count = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "midout_%d", i);
      tempFILE = fopen(str, "r");

      // make sure this file exists
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

      // delete fail_'i' & rawout_'i' & 'nonsolutions_i' - the headnode made a global version of these during printPathFooter
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

  if (!ED_d->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // free the allocate memory
  free(str);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  if (!eachStartPt)
  { 
    for (i = trackCount->numPoints - 1; i >= 0; i--)
      clear_endgame_data(&startPts[i]);
  }
  else
  {
    for (i = maxSize - 1; i >= 0; i--)
      clear_endgame_data(&sendPts[i]);
  }
  free(startPts);
  free(sendPts);
  for (i = numRecvPts - 1; i >= 0; i--)
  {
    clear_endgame_data(&recvPts[i]);
  }
  free(recvPts);

  return;
}

void head_zero_dim_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'headnode'            *
\***************************************************************/
{
  int i, j, recvProc, count = 0, eachStartPt, numRecvPts = 0, maxSize = 0, minPacketSize = 1, maxPacketSize = 20; 
  int *pathNum = NULL, *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *startPts = NULL, *sendPts = NULL, *recvPts = NULL;
  char ch, *str = NULL;
  FILE *tempFILE = NULL, *NONSOLN = NULL;

  // setup NONSOLN
  if (!ED->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // initialize trackCount
  init_trackingStats(trackCount);

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  // send ED to the workers
  bcast_basic_eval_data_mp(ED, 1, my_id, headnode);

  // not using trackBack endgame
  count = 0;
  MPI_Bcast(&count, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
  // Find the number of start points
  fscanf(START, "%d\n", &trackCount->numPoints);

  if (trackCount->numPoints > num_processes * maxPacketSize)
  { // since there is a large number of points, read them in one by one
    eachStartPt = 1;
  }
  else
  { // read them all in at once
    eachStartPt = 0;
  }

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup startPts - if we do not need to setup startPts, they will be read directly from file to sendPts
  if (!eachStartPt)
  { // setup pathNum - randomize the list of start points
    pathNum = (int *)bmalloc(trackCount->numPoints * sizeof(int));
    for (i = 0; i < trackCount->numPoints; i++)
    { // find random integer from 0 to i
      j = rand() % (i + 1);
      // swap i & j
      pathNum[i] = pathNum[j];
      pathNum[j] = i;
    }

    // setup startPts
    startPts = (endgame_data_t *)bmalloc(trackCount->numPoints * sizeof(endgame_data_t));
    for (i = 0; i < trackCount->numPoints; i++)
    { // initialize startPts[pathNum[i]]
      init_endgame_data(&startPts[pathNum[i]], T->Precision);
      startPts[pathNum[i]].prec = T->Precision; // start point in fixed precision
      startPts[pathNum[i]].pathNum = i; // path number
      // read in startPts[pathNum[i]].PD_mp
      setupStart_mp(T, &startPts[pathNum[i]].PD_mp, START);
    }
  }
  else
  { // setup sendPts based on maxSize
    sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
    for (i = 0; i < maxSize; i++)
    { // initialize 
      init_endgame_data(&sendPts[i], T->Precision);
    }
  }

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0)
          if (!((count + j) % pathMod))
            printf("Tracking path %d of %d\n", count + j, trackCount->numPoints);

        if (eachStartPt)
        { // read in from START
          sendPts[j].prec = T->Precision; // start point is in the fixed multi precision
          sendPts[j].pathNum = count + j; // path number is count + j
          setupStart_mp(T, &sendPts[j].PD_mp, START);
        }
      }
      // store the size of this packet = j
      lastSize[i] = j;
      // send this packet to proc i
      if (!eachStartPt)
      { 
        sendPts = &startPts[count];
      }
      // send sendPts
      send_recv_endgame_data_t(&sendPts, &lastSize[i], T->MPType, i, 1);
      if (!eachStartPt)
      {
        sendPts = NULL;
      }

      // update count
      count += lastSize[i];
    }

  // loop until all of the paths have been sent out to the workers to be tracked
  while (count < trackCount->numPoints)
  { // recv a packet back & who sent it
    recvProc = send_recv_endgame_data_t(&recvPts, &numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, trackCount->numPoints - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create the next packet to send back to recvProc
    for (j = 0; j < packetSizes[recvProc]; j++)
    { // print the path number if needed
      if (pathMod > 0)
        if (!((count + j) % pathMod))
          printf("Tracking path %d of %d\n", count + j, trackCount->numPoints);

      if (eachStartPt)
      { // read in from START
        sendPts[j].prec = T->Precision; // start point is in the fixed multi precision
        sendPts[j].pathNum = count + j; // path number is count + j
        setupStart_mp(T, &sendPts[j].PD_mp, START);
      }
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to recvProc
    if (!eachStartPt)
    {
      sendPts = &startPts[count];
    }
 
    // send sendPts
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);
    if (!eachStartPt)
    {
      sendPts = NULL;
    }

    // update count
    count += lastSize[recvProc];

    // print the footers for each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    {
      fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
      // print footer in multi precision
      printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED);
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

    // tell the worker that tracking is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // print the footers for each of the paths
    for (i = 0; i < numRecvPts; i++)
    {
      fprintf(OUT, "Path number: %d (ID: %d)\n", recvPts[i].pathNum, recvProc);
      // print footer in multi precision
      printPathFooter_mp(trackCount, &recvPts[i], T, OUT, RAWOUT, FAIL, NONSOLN, ED);
    }

    // decrement count
    count--;
  }
  // since all of the workers have sent back their data, we just need to wait until the files are closed
  MPI_Barrier(MPI_COMM_WORLD);

  // the only thing that needs combined is all of the MIDOUT files so that the midpoint checker will work
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // open midout_'i'
      count = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, count * sizeof(char));
      sprintf(str, "midout_%d", i);
      tempFILE = fopen(str, "r");

      // make sure this file exists
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

  if (!ED->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // free the allocate memory
  free(str);
  free(pathNum);
  free(packetSizes);
  free(lastSize);
  if (!eachStartPt)
  {
    for (i = trackCount->numPoints - 1; i >= 0; i--)
    {
      clear_endgame_data(&startPts[i]);
    }
  }
  else
  {
    for (i = maxSize - 1; i >= 0; i--)
    {
      clear_endgame_data(&sendPts[i]);
    }
  }
  free(startPts);
  free(sendPts);
  for (i = numRecvPts - 1; i >= 0; i--)
  {
    clear_endgame_data(&recvPts[i]);
  }
  free(recvPts);

  return;
}

void worker_process_main(int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: controls work done by 'worker' processes               *
\***************************************************************/
{
  worker_info recvType;

  // recv the worker info
  bcast_worker_info(&recvType, my_id, headnode);

  if (recvType.dataType == ZERO_DIM_D || recvType.dataType == USER_HOM_D || recvType.dataType == PARAM_HOM_D)
  { // do zero dimensional tracking in either double or adaptive precision
    worker_zero_dim_endgame(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == ZERO_DIM_MP || recvType.dataType == USER_HOM_MP || recvType.dataType == PARAM_HOM_MP)
  { // do zero dimensional tracking in fixed multi precision
    worker_zero_dim_endgame(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == EQBYEQ_D)
  { // do equation-by-equation in either double or adaptive precision
    worker_eqbyeq_d(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == EQBYEQ_MP)
  { // do equation-by-equation in fixed multi precision
    worker_eqbyeq_mp(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == ZERO_DIM_REGEN)
  { // do zero dimensional regeneration
    worker_regen_zero_dim(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == IRREDDECOMP)
  { // do numerical irreducible decomposition
    worker_numericalIrredDecomp(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == REGEN_EXTEND)
  { // use regeneration extension
    worker_regenExtend(my_id, num_processes, headnode, recvType.dataType);
  }

  return;
}

void worker_zero_dim_track_d(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'worker' process      *
\***************************************************************/
{
  int i, size, numStartPts = 0, numEndPts = 0;
  trackingStats trackCount;
  tracker_config_t T;
  basic_eval_data_d BED;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL, *NONSOLN = NULL;
  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // setup eval_func_d & eval_func_mp
  change_prec = &change_basic_eval_prec;
  if (dataType == USER_HOM_D)
  { // user defined homotopy
    eval_func_d = &userHom_eval_d;
    eval_func_mp = &userHom_eval_mp;
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

  if (T.MPType == 2)
  { // setup BED to use AMP
    bcast_basic_eval_data_amp(&BED, my_id, headnode);
  }
  else 
  { // setup BED only for double precision
    bcast_basic_eval_data_d(&BED, T.MPType, my_id, headnode);
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
    for (i = 0; i < numEndPts; i++)
    {
      init_endgame_data(&endPts[i], 64);
    }
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd

    for (i = 0; i < numStartPts; i++)
    { // print the header of the path to OUT
      printPathHeader_d(OUT, &startPts[i].PD_d, &T, startPts[i].pathNum, &BED, eval_func_d);

      // track the path
      zero_dim_track_path_d(startPts[i].pathNum, &endPts[i], &startPts[i].PD_d, OUT, MIDOUT, &T, &BED, BED.BED_mp, eval_func_d, eval_func_mp, change_prec, zero_dim_dehom);

      // check to see if it should be sharpened
      if (endPts[i].retVal == 0 && T.sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], &T, OUT, &BED, BED.BED_mp, eval_func_d, eval_func_mp, change_prec);
      }
     
      if (endPts[i].prec < 64)
      { // print footer in double precision
        printPathFooter_d(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, &BED);
      }
      else
      { // print footer in multi precision
        printPathFooter_mp(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, BED.BED_mp);
      }
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T.MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T.MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
      {
        clear_endgame_data(&endPts[i]);
      }

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
      {
        init_endgame_data(&endPts[i], 64);
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
  basic_eval_clear_d(&BED, 0, T.MPType); // 0 since this does not use regeneration
  tracker_config_clear(&T);
  clearMP();

  return;
}

void worker_zero_dim_track_mp(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking - 'worker' process      *
\***************************************************************/
{ 
  int i, size, numStartPts = 0, numEndPts = 0;
  trackingStats trackCount;
  tracker_config_t T;
  basic_eval_data_mp BED;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL, *NONSOLN = NULL;
  int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // setup eval_func
  if (dataType == USER_HOM_MP)
  { // user defined homotopy
    eval_func = &userHom_eval_mp;
  }
  else
  { // standard zero dim
    eval_func = &basic_eval_mp;
  }

  // setup T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision
  initMP(T.Precision);

  // setup BED only for fixed multi precision
  bcast_basic_eval_data_mp(&BED, 1, my_id, headnode);

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
    for (i = 0; i < numEndPts; i++)
    { // initialize
      init_endgame_data(&endPts[i], T.Precision);
    }
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd

    for (i = 0; i < numStartPts; i++)
    { // print the header of the path to OUT
      printPathHeader_mp(OUT, &startPts[i].PD_mp, &T, startPts[i].pathNum, &BED, eval_func);

      // track the path
      zero_dim_track_path_mp(startPts[i].pathNum, &endPts[i], &startPts[i].PD_mp, OUT, MIDOUT, &T, &BED, eval_func, zero_dim_dehom);

      // check to see if it should be sharpened
      if (endPts[i].retVal == 0 && T.sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&endPts[i], &T, OUT, NULL, &BED, NULL, eval_func, NULL);
      }

      // print footer
      printPathFooter_mp(&trackCount, &endPts[i], &T, OUT, RAWOUT, FAIL, NONSOLN, &BED);
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T.MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T.MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
      {
        clear_endgame_data(&endPts[i]);
      }

      // set the number to reallocate
      numEndPts = numStartPts;

      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      // initialize 
      for (i = 0; i < numEndPts; i++)
      {
        init_endgame_data(&endPts[i], T.Precision);
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

  // startPts & endPts are cleared since numStartPts <= 0

  // clear other allocated memory
  basic_eval_clear_mp(&BED, 0, 1); // 0 since not regeneration, 1 for clear the Prog
  tracker_config_clear(&T);
  clearMP();

  return;
}

void packetSize_maker(int *packetSizes, int num_procs, int headnode, int totalToDistribute, int currProc, int minSize, int maxSize, int *max)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: currProc == headnode - initially setup packetSizes *
*      this will also setup max                                 *
*   otherwise, setup only packetSizes[currProc]                 *
* RETURN VALUES:                                                *
* NOTES: sets packetSizes using some form of guided allocation  *
*  if maxSize <= 0, uses whatever size it wants otherwise limits*
\***************************************************************/
{
  int i, currTotal = totalToDistribute;
  double denominator; // want as a double to divide correctly

  // error checking
  if (num_procs < 2)
  {
    printf("ERROR: The number of processes needs to be >= 2 to use packetSize_maker!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (currProc < 0 || currProc >= num_procs) // acceptable numbers are 0,1,..,num_procs - 1
  {
    printf("ERROR: The current process needs to be between %d and %d!\n", 0, num_procs - 1);
    bexit(ERROR_CONFIGURATION);
  }

  // set minimum size to be 1 if it is <= 0
  if (minSize < 1)
    minSize = 1;
  // if min > max > 0, set max == min
  if (maxSize > 0 && minSize > maxSize)
    maxSize = minSize;
  
  if (currProc == headnode)
  { // initially setup each of them
    denominator = 1.6 * (num_procs - 1); // large packages at first
    *max = 0;
    for (i = 0; i < num_procs; i++)
      if (i != headnode)
      { // make sure the total left is atleast minSize
        if (currTotal < minSize)
        { // the number left is smaller than the minimum packet size, so everything goes into the same packet
          packetSizes[i] = currTotal;
        }
        else // can create a 'full' packet
        { // find the suggested packet size
          packetSizes[i] = ceil(currTotal / denominator);

          // make sure this suggested size is smaller than the maximum
          if (packetSizes[i] > maxSize)
          { // take the size to be the maximum
            packetSizes[i] = maxSize;
          }
          // make sure this suggested size is larger than the minimum
          if (packetSizes[i] < minSize)
          { // take the size to be the minimum
            packetSizes[i] = minSize;
          }
        }
        // update currTotal
        currTotal -= packetSizes[i];
        // check for max
        if (packetSizes[i] > *max)
          *max = packetSizes[i];
      }
  } 
  else
  { // only update currProc
    if (currTotal < minSize)
    { // the number left is smaller than the minimum packet size, so everything goes into the same packet
      packetSizes[currProc] = currTotal;
    }
    else // can create a 'full' packet
    { // find the suggested packet size
      if (currTotal > 5 * maxSize)
        denominator = 1.9 * (num_procs - 1); // smaller packages afterwards
      else
        denominator = 2.5 * (num_procs - 1); // smaller packages afterwards

      packetSizes[currProc] = ceil(currTotal / denominator);

      // make sure this suggested size is smaller than the maximum
      if (packetSizes[currProc] > maxSize)
      { // take the size to be the maximum
        packetSizes[currProc] = maxSize;
      }
      // make sure this suggested size is larger than the minimum
      if (packetSizes[currProc] < minSize)
      { // take the size to be the minimum
        packetSizes[currProc] = minSize;
      }
    }
  }        

  return;
}

///////////// GENERAL PARALLEL FUNCTIONS ////////////////////////

void head_trackPaths(char *strJob, int num_paths, int minPacketSize, int maxPacketSize, trackingStats *trackCount, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*setup_startPoint)(endgame_data_t *, int, int, void const *, void const *), int (*store_endPoint)(endgame_data_t *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks paths in parallel                               *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, maxSize, path_num;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum - randomize the list of start points
  for (i = 0; i < num_paths; i++)
  { // find random integer from 0 to i
    j = rand() % (i + 1);
    // swap i & j
    pathNum[i] = pathNum[j];
    pathNum[j] = i;
  }

  // allocate & intialize sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  }
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create the packet
      for (j = 0; j < packetSizes[i]; j++)
      { // print the path number if needed
        if (pathMod > 0 && !((count + j) % pathMod))
          printf("%s path %d of %d\n", strJob, count + j, num_paths);

        // find the path number
        path_num = pathNum[count + j];

        // setup sendPts[j]
        setup_startPoint(&sendPts[j], path_num, T->MPType, ED_d, ED_mp);
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
        printf("%s path %d of %d\n", strJob, count + j, num_paths);

      // find the path number
      path_num = pathNum[count + j];

      // setup sendPts[j]
      setup_startPoint(&sendPts[j], path_num, T->MPType, ED_d, ED_mp);
    }
    // store the size of this packet = j
    lastSize[recvProc] = j;
    // send this packet to proc recvProc
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // update count
    count += lastSize[recvProc];

    // store each of the paths recvd
    for (i = 0; i < numRecvPts; i++)
    {
      store_endPoint(&recvPts[i], trackCount, T, OUT, RAWOUT, FAIL, ED_d, ED_mp);
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
    {
      store_endPoint(&recvPts[i], trackCount, T, OUT, RAWOUT, FAIL, ED_d, ED_mp);
    }

    // tell the worker that this level is complete
    lastSize[recvProc] = 0;
    send_recv_endgame_data_t(&sendPts, &lastSize[recvProc], T->MPType, recvProc, 1);

    // decrement count
    count--;
  }

  // clear the memory
  for (i = numRecvPts - 1; i >= 0; i--)
  {
    clear_endgame_data(&recvPts[i]);
  }
  free(recvPts);
  for (i = maxSize - 1; i >= 0; i--)
  {
    clear_endgame_data(&sendPts[i]);
  }
  free(sendPts);
  free(pathNum);
  free(packetSizes);
  free(lastSize);

  return;
}

void worker_trackPaths(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*worker_track)(endgame_data_t *, endgame_data_t *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *)))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks paths in parallel                               *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0;
  endgame_data_t *startPts = NULL, *endPts = NULL;

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
    {
      init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd

    for (i = 0; i < numStartPts; i++)
    { // track the path
      worker_track(&startPts[i], &endPts[i], OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener);
    }

    // send endPts
    send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);

    // recv next set of start points
    send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
      {
        clear_endgame_data(&endPts[i]);
      }

      // set the number to reallocate
      numEndPts = numStartPts;

      // reallocate
      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
      {
       init_endgame_data(&endPts[i], T->Precision);
      }
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_sortPaths(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*setup_worker_track)(int, int, void const *, void const *), int (*classifyEndpoint)(endgame_data_t *, tracker_config_t *, FILE *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int)))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sort paths in parallel                                 *
\***************************************************************/
{
  int i, numStartPts = 0, numEndPts = 0;
  endgame_data_t *startPts = NULL, *endPts = NULL;

  // recv the first set of start points
  send_recv_endgame_data_t(&startPts, &numStartPts, T->MPType, headnode, 0);

  // setup endPts
  if (numStartPts > 0)
  { // allocate endPts
    numEndPts = numStartPts;
    endPts = (endgame_data_t *)bmalloc(numEndPts * sizeof(endgame_data_t));
    for (i = 0; i < numEndPts; i++)
    {
      init_endgame_data(&endPts[i], T->Precision);
    }
  }

  // main loop
  while (numStartPts > 0)
  { // track the paths that were recvd

    for (i = 0; i < numStartPts; i++)
    { // setup for sorting the path
      setup_worker_track(startPts[i].pathNum, T->MPType, ED_d, ED_mp);

      // setup endPts
      endPts[i].prec = startPts[i].prec;
      endPts[i].pathNum = startPts[i].pathNum;
      if (endPts[i].prec < 64)
      {
        point_cp_d(endPts[i].PD_d.point, startPts[i].PD_d.point);
        set_d(endPts[i].PD_d.time, startPts[i].PD_d.time);
      }
      else
      {
        setprec_point_mp(endPts[i].PD_mp.point, endPts[i].prec);
        setprec_mp(endPts[i].PD_mp.time, endPts[i].prec);

        point_cp_mp(endPts[i].PD_mp.point, startPts[i].PD_mp.point);
        set_mp(endPts[i].PD_mp.time, startPts[i].PD_mp.time);
      }

      if (startPts[i].retVal == 0)
      { // determine the outcome
        classifyEndpoint(&endPts[i], T, OUT, ED_d, ED_mp, eval_d, eval_mp, change_prec);
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

    // setup endPts
    if (numEndPts != numStartPts)
    { // clear endPts
      for (i = numEndPts - 1; i >= 0; i--)
      {
        clear_endgame_data(&endPts[i]);
      }

      // set the number to reallocate
      numEndPts = numStartPts;

      // reallocate
      endPts = (endgame_data_t *)brealloc(endPts, numEndPts * sizeof(endgame_data_t));
      for (i = 0; i < numEndPts; i++)
      {
        init_endgame_data(&endPts[i], T->Precision);
      }
    }
  }

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

// general path tracker

int worker_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: worker process to track the path                       *
\***************************************************************/
{
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  // print the header
  if (startPt->prec < 64)
  { // print the header of the path to OUT using _d
    printPathHeader_d(OUT, &startPt->PD_d, T, startPt->pathNum, CD_d, eval_d);
  }
  else
  { // print the header of the path to OUT using _mp
    printPathHeader_mp(OUT, &startPt->PD_mp, T, startPt->pathNum, CD_mp, eval_mp);
  }

  // track the path
  if (T->MPType == 0 || T->MPType == 2)
  { // track with _d or _amp
    zero_dim_track_path_d(startPt->pathNum, endPt, &startPt->PD_d, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, find_dehom);
  }
  else
  { // track with _mp
    zero_dim_track_path_mp(startPt->pathNum, endPt, &startPt->PD_mp, OUT, MIDOUT, T, CD_mp, eval_mp, find_dehom);
  }

  if (useSharpener(endPt->retVal, T->sharpenDigits, CD_d, CD_mp))
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(endPt, T, OUT, CD_d, CD_mp, eval_d, eval_mp, change_prec);
  }

  printResultOfPath(OUT, endPt->retVal, T);

  if (endPt->prec < 64)
  { // print footer in double precision
    printBasicFooter_d(OUT, &endPt->PD_d, T, endPt->function_residual_d);
  }
  else
  { // print footer in multi precision
    printBasicFooter_mp(OUT, &endPt->PD_mp, T, endPt->function_residual_mp);
  }

  return 0;
}

int worker_track_path_rank(int rankType, int *rankDef, int *corank, double *smallest_nonzero_sv, double *largest_zero_sv, endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: worker process to track the path                       *
\***************************************************************/
{
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  // print the header
  if (startPt->prec < 64)
  { // print the header of the path to OUT using _d
    printPathHeader_d(OUT, &startPt->PD_d, T, startPt->pathNum, CD_d, eval_d);
  }
  else
  { // print the header of the path to OUT using _mp
    printPathHeader_mp(OUT, &startPt->PD_mp, T, startPt->pathNum, CD_mp, eval_mp);
  }

  // track the path
  if (T->MPType == 0 || T->MPType == 2)
  { // track with _d or _amp

    // track the path
    zero_dim_track_path_rank_d(startPt->pathNum, rankType, rankDef, corank, smallest_nonzero_sv, largest_zero_sv, endPt, &startPt->PD_d, OUT, MIDOUT, T, CD_d, CD_mp, eval_d, eval_mp, change_prec, find_dehom);
  }
  else
  { // track with _mp   

    // track the path
    zero_dim_track_path_rank_mp(startPt->pathNum, rankType, rankDef, corank, smallest_nonzero_sv, largest_zero_sv, endPt, &startPt->PD_mp, OUT, MIDOUT, T, CD_mp, eval_mp, find_dehom);
  }

  if (useSharpener(endPt->retVal, T->sharpenDigits, CD_d, CD_mp) && (rankType ? *corank == 0 : *rankDef == 0))
  { // use the sharpener for after an endgame, if possible
    sharpen_endpoint_endgame(endPt, T, OUT, CD_d, CD_mp, eval_d, eval_mp, change_prec);
  }

  printResultOfPath(OUT, endPt->retVal, T);

  if (endPt->prec < 64)
  { // print footer in double precision
    printBasicFooter_d(OUT, &endPt->PD_d, T, endPt->function_residual_d);
  }
  else
  { // print footer in multi precision
    printBasicFooter_mp(OUT, &endPt->PD_mp, T, endPt->function_residual_mp);
  }

  return 0;
}

// general sharpening condition

int worker_useSharpener(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine whether to use sharpener or not              *
\***************************************************************/
{
  int rV = 0;

  if (retVal == 0 && sharpenDigits > 0)
    rV = 1; // sharpen the endpoint
  else
    rV = 0; // do not sharpen the endpoint

  return rV;
}

// general setup 

int worker_setup(int pathNum, int MPType, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup to work with the path                            *
\***************************************************************/
{ // nothing is needed to be done
  return 0;
}

void delete_parallel_files(char *fileName, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: deletes all of the worker files of the form            *
*   "%s_%d", fileName, id                                       *
\***************************************************************/
{
  int i, size;
  char *str = NULL;

  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    {
      size = 1 + snprintf(NULL, 0, "%s_%d", fileName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", fileName, i);

      remove(str);
    }

  free(str);

  return;
}

void combine_midpath_data(char *midFile, char *parFile, int delete_files, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: combines all of the midpath_data files from workers    *
\***************************************************************/
{
  int i, size;
  char ch, *str = NULL;
  FILE *tempFILE = NULL, *midOUT = fopen(midFile, "w");
  if (midOUT == NULL)
  {
    printf("ERROR: '%s' is not a valid name for a file!\n", midFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // open the parallel files
      size = 1 + snprintf(NULL, 0, "%s_%d", parFile, i);
      str = (char *)brealloc(str, size * sizeof(int));
      sprintf(str, "%s_%d", parFile, i);
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
        fprintf(midOUT, "%c", ch);
        ch = fgetc(tempFILE);
      }

      // close file and delete it
      fclose(tempFILE);
      if (delete_files)
        remove(str);
    }

  // close midOUT
  fclose(midOUT);
  // free str
  free(str);

  return;
}

int parallel_midpoint_checking(char *midFile, char *parFile, int delete_files, int num_paths, int numVars, double midpoint_tol, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of path crossings                       *
* NOTES: checks for path crossing                               *
\***************************************************************/
{
  int retVal = 0;

  // combine the midpath_data files and then delete them
  combine_midpath_data(midFile, parFile, delete_files, headnode, num_processes);

  // find the number of crossings
  retVal = 0;
  midpoint_checker(num_paths, numVars, midpoint_tol, &retVal);

  // return the number of crossings
  return retVal;
}

void head_trackPaths2(char *strJob, int randomize_paths, int num_paths, int minPacketSize, int maxPacketSize, trackingStats *trackCount, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int), FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, int my_id, int headnode, int num_processes, int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int), int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks paths in parallel                               *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, maxSize;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum
  if (randomize_paths)
  { // randomize the list of start points
    for (i = 0; i < num_paths; i++)
    { // find random integer from 0 to i
      j = rand() % (i + 1);
      // swap i & j
      pathNum[i] = pathNum[j];
      pathNum[j] = i;
    }
  }
  else
  { // identity
    for (i = 0; i < num_paths; i++)
      pathNum[i] = i;
  }
  
  // allocate & intialize sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], T->Precision);
  }
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create & send the packet - returns the size of the packet
      lastSize[i] = create_send_packet(count, packetSizes[i], START, sendPts, pathNum, T->MPType, pathMod, ED_d, ED_mp, strJob, num_paths, i);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = recv_store_packet(&recvPts, &numRecvPts, trackCount, T, OUT, RAWOUT, FAIL, OTHER, OTHER2, rV, ED_d, ED_mp, change_prec);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create & send the next packet - returns the size of the packet
    lastSize[recvProc] = create_send_packet(count, packetSizes[recvProc], START, sendPts, pathNum, T->MPType, pathMod, ED_d, ED_mp, strJob, num_paths, recvProc);

    // update count
    count += lastSize[recvProc];
  }

  // now that all of the paths have been sent out, we need to loop to recv all the packets back
  // count the number of packets still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  while (count > 0)
  { // recv a packet back & who sent it
    recvProc = recv_store_packet(&recvPts, &numRecvPts, trackCount, T, OUT, RAWOUT, FAIL, OTHER, OTHER2, rV, ED_d, ED_mp, change_prec);

    // tell the worker that this level is complete
    packetSizes[recvProc] = 0;
    lastSize[recvProc] = create_send_packet(count, packetSizes[recvProc], START, sendPts, pathNum, T->MPType, pathMod, ED_d, ED_mp, strJob, num_paths, recvProc);

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

  return;
}

void worker_trackPaths2(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *)))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks paths in parallel                               *
\***************************************************************/
{
  int numStartPts = 0, numEndPts = 0;
  endgame_data_t *startPts = NULL, *endPts = NULL;

  do
  { // recv a packet, track the paths and send it back
    recv_track_send_packet(headnode, &startPts, &numStartPts, &endPts, &numEndPts, OUT, MIDOUT, T, ED_d, ED_mp, eval_d, eval_mp, change_prec, useSharpener);

  } while (numStartPts > 0);

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

#endif

