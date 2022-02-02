// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

#ifdef _HAVE_MPI

void head_cascadeTrackCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int pathMod, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void head_cascadeSortCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int pathMod, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, double final_tol, int my_id, int headnode, int num_processes);

int setup_cascade_startPoint(endgame_data_t *startPt, int path_num, int MPType, void const *CD_d, void const *CD_mp, int setupSort);

void worker_cascadeTrackCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

int cascade_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int cascade_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int cascade_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

void worker_cascadeSortCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes);

int cascade_create_send_packet_sort(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int cascade_recv_store_packet_sort(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int cascade_recv_sort_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

void cascade_par_track(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, double midpoint_tol, tracker_config_t *T, cascade_t *CD, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does cascade tracking - 'headnode'                     *
\***************************************************************/
{
  int size, codim_index, codim, num_paths, num_crossings = 0, num_codim = CD->num_codim, minPacketSize = 1, maxPacketSize = 50;
  char *str = NULL;

  // send T to the workers
  bcast_tracker_config_t(T, my_id, headnode);

  // send CD to the workers
  bcast_cascade_t(CD, T->MPType, my_id, headnode);

  // send the first codimension data to the workers
  bcast_cascadeCodim_t(&CD->codim[0], my_id, headnode);

  // loop over the codimensions to find the witness supersets
  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { // track the paths for this codimension
    codim = CD->codim[codim_index].codim;
    num_paths = CD->codim[codim_index].num_paths;

    head_cascadeTrackCodim(minPacketSize, maxPacketSize, trackCount, codim_index, pathMod, T, CD, OUT, RAWOUT, FAIL, my_id, headnode, num_processes);

    // check for path crossings for this codimension

    // wait until the files are closed
    MPI_Barrier(MPI_COMM_WORLD);

    // setup str
    size = 1 + snprintf(NULL, 0, "midout_%d", codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "midout_%d", codim);

    num_crossings = parallel_midpoint_checking(midFile, str, 1, num_paths, CD->new_variables, midpoint_tol, my_id, headnode, num_processes);

    if (num_crossings > 0)
      printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

    // allow the workers to open MIDOUT again
    MPI_Barrier(MPI_COMM_WORLD);

    // sort the endpoints for this codimension
    head_cascadeSortCodim(minPacketSize, maxPacketSize, trackCount, codim_index, pathMod, T, CD, OUT, RAWOUT, FAIL, T->final_tol_times_mult, my_id, headnode, num_processes);

    if (codim_index + 1 < num_codim)
    { // setup the next codim to track
      cascadePrepareNextCodim(CD, codim_index, T->MPType);

      // send the next codim to the workers
      bcast_cascadeCodim_t(&CD->codim[codim_index + 1], my_id, headnode);
    }

    // clear the start points for this codim since they are no longer needed
    cascade_clear_start_points(CD, codim_index, T->MPType);
  }

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  for (codim_index = 0; codim_index < num_processes; codim_index++)
    if (codim_index != headnode)
    { // delete fail_'i' & rawout_'i' - the headnode made a global version of these during printPathFooter
      size = 1 + snprintf(NULL, 0, "rawout_%d", codim_index);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", codim_index);
      remove(str);

      size = 1 + snprintf(NULL, 0, "fail_%d", codim_index);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", codim_index);
      remove(str);
    }

  free(str);

  return;
}

int setup_cascade_startPoint(endgame_data_t *startPt, int path_num, int MPType, void const *CD_d, void const *CD_mp, int setupSort)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup startPt                                          *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)CD_d;
  int codim_index = CD->curr_codim_index;

  // setup pathNum & retVal
  startPt->pathNum = path_num;
  startPt->retVal = 0;

  // setup startPt
  if (MPType == 0)
  { // setup startPt.PD_d
    startPt->prec = startPt->last_approx_prec = 52;
    if (setupSort)
    { // setup for sorting
      point_cp_d(startPt->PD_d.point, CD->codim[codim_index].endPts_d[path_num].endPt);
      set_zero_d(startPt->PD_d.time);

      // setup last_approx
      point_cp_d(startPt->last_approx_d, CD->codim[codim_index].endPts_d[path_num].last_approx);

      // copy retVal/corank/CN
      startPt->retVal = CD->codim[codim_index].endPts_d[path_num].retVal;
      startPt->codim = CD->codim[codim_index].endPts_d[path_num].corank;
      startPt->condition_number = CD->codim[codim_index].endPts_d[path_num].cond_num;
    }
    else
    { // setup for tracking
      point_cp_d(startPt->PD_d.point, CD->codim[codim_index].startPts_d[path_num]);
      set_one_d(startPt->PD_d.time);

      startPt->last_approx_d->size = 0;
    }
  }
  else if (MPType == 1)
  { // setup startPt.PD_mp
    startPt->prec = startPt->last_approx_prec = CD->curr_precision;
    setprec_point_mp(startPt->PD_mp.point, CD->curr_precision);
    setprec_mp(startPt->PD_mp.time, CD->curr_precision);
    if (setupSort)
    { // setup for sorting
      point_cp_mp(startPt->PD_mp.point, CD->codim[codim_index].endPts_mp[path_num].endPt);
      set_zero_mp(startPt->PD_mp.time);

      // setup last_approx
      point_cp_mp(startPt->last_approx_mp, CD->codim[codim_index].endPts_mp[path_num].last_approx);

      // copy retVal/corank/CN
      startPt->retVal = CD->codim[codim_index].endPts_mp[path_num].retVal;
      startPt->codim = CD->codim[codim_index].endPts_mp[path_num].corank;
      startPt->condition_number = CD->codim[codim_index].endPts_mp[path_num].cond_num;
    }
    else
    { // setup for tracking
      point_cp_mp(startPt->PD_mp.point, CD->codim[codim_index].startPts_mp[path_num]);
      set_one_mp(startPt->PD_mp.time);

      startPt->last_approx_mp->size = 0;
    }
  }
  else
  { // setup based on precision
    if (setupSort)
    { // setup for sorting
      startPt->prec = CD->codim[codim_index].endPts_amp[path_num].curr_prec;

      if (startPt->prec < 64)
      { // setup PD_d
        point_cp_d(startPt->PD_d.point, CD->codim[codim_index].endPts_amp[path_num].endPt_d);
        set_zero_d(startPt->PD_d.time);
      }
      else
      { // setup PD_mp
        setprec_point_data_mp(&startPt->PD_mp, startPt->prec);

        point_cp_mp(startPt->PD_mp.point, CD->codim[codim_index].endPts_amp[path_num].endPt_mp);
        set_zero_mp(startPt->PD_mp.time);
      }

      startPt->last_approx_prec = CD->codim[codim_index].endPts_amp[path_num].last_approx_prec;
      
      if (startPt->last_approx_prec < 64)
      { // setup _d
        point_cp_d(startPt->last_approx_d, CD->codim[codim_index].endPts_amp[path_num].last_approx_d);
      }
      else
      { // setup _mp
        setprec_point_mp(startPt->last_approx_mp, startPt->last_approx_prec);
        point_cp_mp(startPt->last_approx_mp, CD->codim[codim_index].endPts_amp[path_num].last_approx_mp);
      }

      // copy retVal/corank/CN
      startPt->retVal = CD->codim[codim_index].endPts_amp[path_num].retVal;
      startPt->codim = CD->codim[codim_index].endPts_amp[path_num].corank;
      startPt->condition_number = CD->codim[codim_index].endPts_amp[path_num].cond_num;
    }
    else
    { // setup for tracking
      startPt->prec = startPt->last_approx_prec = 52;
      point_cp_d(startPt->PD_d.point, CD->codim[codim_index].startPts_d[path_num]);
      set_one_d(startPt->PD_d.time);

      startPt->last_approx_d->size = 0;
    }
  }

  return 0;
}

int store_cascade_endPoint(endgame_data_t *endPt, int corank, double smallest, double largest, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, void const *CD_d, void const *CD_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: store endPt                                            *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)CD_d;
  int codim_index = CD->curr_codim_index, path_num = endPt->pathNum;

  if (T->MPType == 0)
  { // store to _d
    point_d dehom_d, orig_vars;
    init_point_d(dehom_d, 0); init_point_d(orig_vars, 0);

    // store the condition number
    CD->codim[codim_index].endPts_d[path_num].cond_num = endPt->condition_number;

    // copy over to the appropriate spot
    set_d(CD->codim[codim_index].endPts_d[path_num].finalT, endPt->PD_d.time);
    point_cp_d(CD->codim[codim_index].endPts_d[path_num].endPt, endPt->PD_d.point);

    // find the point in the original coordinates and dehomogenized point
    cascadeFindOrigVarsDehom_d(orig_vars, dehom_d, endPt->PD_d.point, CD);

    // print the footer for the point
    CD->codim[codim_index].endPts_d[path_num].retVal = printCascadeFooter_d(CD, codim_index, path_num, &endPt->PD_d, orig_vars, dehom_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

    // store the other info
    point_cp_d(CD->codim[codim_index].endPts_d[path_num].last_approx, endPt->last_approx_d);
    CD->codim[codim_index].endPts_d[path_num].corank = corank;
    CD->codim[codim_index].endPts_d[path_num].smallest_nonzero_SV = smallest;
    CD->codim[codim_index].endPts_d[path_num].largest_zero_SV = largest;

    clear_point_d(dehom_d); clear_point_d(orig_vars);
  }
  else if (T->MPType == 1)
  { // store to _mp
    point_mp dehom_mp, orig_vars;
    init_point_mp(dehom_mp, 0); init_point_mp(orig_vars, 0);

    // store the condition number
    CD->codim[codim_index].endPts_mp[path_num].cond_num = endPt->condition_number;

    // copy over to the appropriate spot
    set_mp(CD->codim[codim_index].endPts_mp[path_num].finalT, endPt->PD_mp.time);
    point_cp_mp(CD->codim[codim_index].endPts_mp[path_num].endPt, endPt->PD_mp.point);

    // find the point in the original coordinates and dehomogenized point
    cascadeFindOrigVarsDehom_mp(orig_vars, dehom_mp, endPt->PD_mp.point, CD);

    // print the footer for the point
    CD->codim[codim_index].endPts_mp[path_num].retVal = printCascadeFooter_mp(CD, codim_index, path_num, &endPt->PD_mp, orig_vars, dehom_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

    // store the other info
    point_cp_mp(CD->codim[codim_index].endPts_mp[path_num].last_approx, endPt->last_approx_mp);
    CD->codim[codim_index].endPts_mp[path_num].corank = corank;
    CD->codim[codim_index].endPts_mp[path_num].smallest_nonzero_SV = smallest;
    CD->codim[codim_index].endPts_mp[path_num].largest_zero_SV = largest;

    clear_point_mp(dehom_mp); clear_point_mp(orig_vars);
  }
  else
  { // store to _amp

    // store the condition number
    CD->codim[codim_index].endPts_amp[path_num].cond_num = endPt->condition_number;
    // store the precision
    CD->codim[codim_index].endPts_amp[path_num].curr_prec = endPt->prec;

    if (endPt->prec < 64)
    {
      point_d dehom_d, orig_vars;
      init_point_d(dehom_d, 0); init_point_d(orig_vars, 0);

      // copy over to the appropriate spot
      set_d(CD->codim[codim_index].endPts_amp[path_num].finalT_d, endPt->PD_d.time);
      point_cp_d(CD->codim[codim_index].endPts_amp[path_num].endPt_d, endPt->PD_d.point);

      // find the point in the original coordinates and dehomogenized point
      cascadeFindOrigVarsDehom_d(orig_vars, dehom_d, endPt->PD_d.point, CD);

      // print the footer for the point
      CD->codim[codim_index].endPts_amp[path_num].retVal = printCascadeFooter_d(CD, codim_index, path_num, &endPt->PD_d, orig_vars, dehom_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

      clear_point_d(dehom_d); clear_point_d(orig_vars);
    }
    else
    {
      point_mp dehom_mp, orig_vars;
      init_point_mp2(dehom_mp, 0, endPt->prec); init_point_mp2(orig_vars, 0, endPt->prec);

      // initialize memory in endPts_amp
      setprec_point_mp(CD->codim[codim_index].endPts_amp[path_num].endPt_mp, endPt->prec);
      setprec_mp(CD->codim[codim_index].endPts_amp[path_num].finalT_mp, endPt->prec);

      // copy over to the appropriate spot - converting to double precision
      set_mp(CD->codim[codim_index].endPts_amp[path_num].finalT_mp, endPt->PD_mp.time);
      point_cp_mp(CD->codim[codim_index].endPts_amp[path_num].endPt_mp, endPt->PD_mp.point);

      // find the point in the original coordinates and dehomogenized point
      cascadeFindOrigVarsDehom_mp(orig_vars, dehom_mp, endPt->PD_mp.point, CD);

      // print the footer for the point
      CD->codim[codim_index].endPts_amp[path_num].retVal = printCascadeFooter_mp(CD, codim_index, path_num, &endPt->PD_mp, orig_vars, dehom_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

      clear_point_mp(dehom_mp); clear_point_mp(orig_vars);
    }

    // store the other info
    CD->codim[codim_index].endPts_amp[path_num].corank = corank;
    CD->codim[codim_index].endPts_amp[path_num].smallest_nonzero_SV = smallest;
    CD->codim[codim_index].endPts_amp[path_num].largest_zero_SV = largest;

    CD->codim[codim_index].endPts_amp[path_num].last_approx_prec = endPt->last_approx_prec;
    if (endPt->last_approx_prec < 64)
    { // copy _d
      point_cp_d(CD->codim[codim_index].endPts_amp[path_num].last_approx_d, endPt->last_approx_d);
    }
    else
    { // copy _mp
      change_prec_point_mp(CD->codim[codim_index].endPts_amp[path_num].last_approx_mp, endPt->last_approx_prec);
      point_cp_mp(CD->codim[codim_index].endPts_amp[path_num].last_approx_mp, endPt->last_approx_mp);
    }
  }

  return 0;
}

int store_cascade_sortPoint(endgame_data_t *endPt, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, void const *CD_d, void const *CD_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: store endPt                                            *
\***************************************************************/
{
  cascade_t *CD = NULL;
  int codim_index, path_num = endPt->pathNum;

  if (T->MPType == 0 || T->MPType == 2)
    CD = (cascade_t *)CD_d; 
  else
    CD = (cascade_t *)CD_mp;

  codim_index = CD->curr_codim_index;

  // store the endpoint type
  CD->codim[codim_index].endPt_types[path_num] = endPt->retVal;

  if (T->MPType == 0)
  { // copy back to endPts_d
    point_cp_d(CD->codim[codim_index].endPts_d[path_num].endPt, endPt->PD_d.point);
    point_cp_d(CD->codim[codim_index].endPts_d[path_num].last_approx, endPt->last_approx_d);
  }
  else if (T->MPType == 1)
  { // copy back to endPts_mp
    point_cp_mp(CD->codim[codim_index].endPts_mp[path_num].endPt, endPt->PD_mp.point);
    point_cp_mp(CD->codim[codim_index].endPts_mp[path_num].last_approx, endPt->last_approx_mp);
  }
  else
  { // copy back to endPts_amp
    CD->codim[codim_index].endPts_amp[path_num].curr_prec = endPt->prec;
    if (endPt->prec < 64)
    { // copy to endPt_d
      point_cp_d(CD->codim[codim_index].endPts_amp[path_num].endPt_d, endPt->PD_d.point);
    }
    else
    { // copy to endPt_mp  
      setprec_point_mp(CD->codim[codim_index].endPts_amp[path_num].endPt_mp, endPt->prec);
      point_cp_mp(CD->codim[codim_index].endPts_amp[path_num].endPt_mp, endPt->PD_mp.point);
    }

    CD->codim[codim_index].endPts_amp[path_num].last_approx_prec = endPt->last_approx_prec;
    if (endPt->last_approx_prec < 64)
    { // copy to endPt_d
      point_cp_d(CD->codim[codim_index].endPts_amp[path_num].last_approx_d, endPt->last_approx_d);
    }
    else
    { // copy to endPt_mp  
      setprec_point_mp(CD->codim[codim_index].endPts_amp[path_num].last_approx_mp, endPt->last_approx_prec);
      point_cp_mp(CD->codim[codim_index].endPts_amp[path_num].last_approx_mp, endPt->last_approx_mp);
    }
  }

  CD = NULL;

  return 0;
}

void head_cascadeTrackCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int pathMod, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all of the paths for this codim in parallel     *
\***************************************************************/
{
  int (*change_prec)(void const *, int);
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the evaluators
  change_prec = &change_cascade_prec;
  create_send_packet = &cascade_create_send_packet_track;
  recv_store_packet = &cascade_recv_store_packet_track;

  // setup the current codimension
  CD->curr_codim_index = codim_index;

  // display messages
  printf("\nFinding witness superset for codimension %d of %d: %d path%s to track.\n", CD->codim[codim_index].codim, CD->num_codim, CD->codim[codim_index].num_paths, CD->codim[codim_index].num_paths == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Finding witness superset for codimension %d.\n", CD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual parallel tracking
  head_trackPaths2("Tracking", 0, CD->codim[codim_index].num_paths, minPacketSize, maxPacketSize, trackCount, pathMod, T, CD, CD, change_prec, NULL, OUT, RAWOUT, FAIL, NULL, NULL, NULL, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

void head_cascadeSortCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int pathMod, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, double final_tol, int my_id, int headnode, int num_processes) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at the current codim         *
\***************************************************************/
{
  int (*change_prec)(void const *, int);
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the evaluators
  change_prec = &change_cascade_prec;
  create_send_packet = &cascade_create_send_packet_sort;
  recv_store_packet = &cascade_recv_store_packet_sort;

  // setup the current codimension
  CD->curr_codim_index = codim_index;

  // display messages
  printf("\nSorting codimension %d of %d: %d path%s to sort.\n", CD->codim[codim_index].codim, CD->num_codim, CD->codim[codim_index].num_paths, CD->codim[codim_index].num_paths == 1 ? "" : "s");

  // do the actual parallel sorting
  head_trackPaths2("Sorting", 0, CD->codim[codim_index].num_paths, minPacketSize, maxPacketSize, trackCount, pathMod, T, CD, CD, change_prec, NULL, OUT, RAWOUT, FAIL, NULL, NULL, NULL, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  // do the final classificiation
  cascade_classifyCodim(CD, codim_index, final_tol, T->MPType);

  return;
}

void worker_cascade(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: MPType                                         *
* NOTES: does cascade tracking - 'worker' process               *
\***************************************************************/
{
  int size, codim_index, num_codim;
  trackingStats trackCount;
  tracker_config_t T;
  cascade_t CD;
  char *str = NULL;
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // initialize trackCount - even though it is not really used
  init_trackingStats(&trackCount);

  // recv T
  bcast_tracker_config_t(&T, my_id, headnode);

  // now that we know the precision, set the default to that precision - does not hurt if doing only double precision
  initMP(T.Precision);

  // recv CD
  bcast_cascade_t(&CD, T.MPType, my_id, headnode);

  // recv first codimension information
  num_codim = CD.num_codim;
  CD.codim = (cascadeCodim_t *)bmalloc(num_codim * sizeof(cascadeCodim_t));
  // recv first codimension data from the headnode
  bcast_cascadeCodim_t(&CD.codim[0], my_id, headnode);

  // setup the local files - OUT, MIDOUT, RAWOUT & FAIL
  size = 1 + snprintf(NULL, 0, "output_%d", my_id);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "output_%d", my_id);
  OUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "midout_%d_%d", CD.codim[0].codim, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "midout_%d_%d", CD.codim[0].codim, my_id);
  MIDOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "rawout_%d", my_id);
  RAWOUT = fopen(str, "w");

  size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "fail_%d", my_id);
  FAIL = fopen(str, "w");

  // main loop - over each codim
  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { // track this codimension
    worker_cascadeTrackCodim(&trackCount, codim_index, &T, &CD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    // close MIDOUT
    fclose(MIDOUT);
    MIDOUT = NULL;
    // wait until all workers have closed MIDOUT
    MPI_Barrier(MPI_COMM_WORLD);

    // consider doing midpoint checking in parallel

    // wait until the headnode has used MIDOUT
    MPI_Barrier(MPI_COMM_WORLD);

    // sort the endpoints for this codimension
    worker_cascadeSortCodim(&trackCount, codim_index, &T, &CD, OUT, RAWOUT, FAIL, my_id, headnode, num_processes);

    if (codim_index + 1 < num_codim)
    { // recv the next codim
      bcast_cascadeCodim_t(&CD.codim[codim_index + 1], my_id, headnode); 

      // reopen MIDOUT
      size = 1 + snprintf(NULL, 0, "midout_%d_%d", CD.codim[codim_index + 1].codim, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d_%d", CD.codim[codim_index + 1].codim, my_id);
      MIDOUT = fopen(str, "w");
    }
  }

  // close the files
  fclose(OUT);
  fclose(RAWOUT);
  fclose(FAIL);

  // clear other allocated memory
  cascade_clear(&CD, T.MPType);
  tracker_config_clear(&T);
  clearMP();
  free(str);

  // wait until all workers have closed the files
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

void worker_cascadeTrackCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
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
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));
  int (*useSharpener)(int, int, void const *, void const *);

  // setup the evaluators
  if (codim_index == 0)
  { // initial evaluation
    ptr_to_eval_d = &initial_codim_cascade_eval_d;
    ptr_to_eval_mp = &initial_codim_cascade_eval_mp;
  }
  else
  { // standard evaluation
    ptr_to_eval_d = &standard_codim_cascade_eval_d;
    ptr_to_eval_mp = &standard_codim_cascade_eval_mp;
  }
  change_prec = &change_cascade_prec;
  recv_track_send_packet = &cascade_recv_track_send_packet;
  useSharpener = &worker_useSharpener;

  // setup the current codimension
  CD->curr_codim_index = codim_index;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Tracking codimension %d.\n", CD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual tracking
  worker_trackPaths2(trackCount, T, CD, CD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  return;
}

void worker_cascadeSortCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, cascade_t *CD, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes)
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
  int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *)) = NULL;
  int (*useSharpener)(int, int, void const *, void const *) = NULL;

  // setup the evaluators
  if (codim_index == 0)
  { // initial evaluation
    ptr_to_eval_d = &initial_codim_cascade_eval_d;
    ptr_to_eval_mp = &initial_codim_cascade_eval_mp;
  }
  else
  { // standard evaluation
    ptr_to_eval_d = &standard_codim_cascade_eval_d;
    ptr_to_eval_mp = &standard_codim_cascade_eval_mp;
  }
  change_prec = &change_cascade_prec;
  recv_track_send_packet = &cascade_recv_sort_send_packet;
  useSharpener = &worker_useSharpener;

  // setup the current codimension
  CD->curr_codim_index = codim_index;

  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Sorting codimension %d.\n", CD->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // do the actual sorting
  worker_trackPaths2(trackCount, T, CD, CD, OUT, RAWOUT, NULL, FAIL, my_id, headnode, num_processes, ptr_to_eval_d, ptr_to_eval_mp, change_prec, useSharpener, recv_track_send_packet);

  return;
}

int cascade_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for cascade tracking      *
\***************************************************************/
{
  int i;
  cascade_t *CD = NULL;

  // setup CD
  if (MPType == 0 || MPType == 2)
    CD = (cascade_t *)ED_d;
  else
    CD = (cascade_t *)ED_mp;

  // create the packet
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    setup_cascade_startPoint(&sendPts[i], pathNum[startNum + i], MPType, CD, CD, 0);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // clear memory
  CD = NULL;

  return size;
}

int cascade_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, recvProc, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  cascade_t *CD = NULL;

  // setup CD
  if (T->MPType == 0 || T->MPType == 2)
    CD = (cascade_t *)ED_d;
  else
    CD = (cascade_t *)ED_mp;

  // recv the packet
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // recv the last approximations
  corank = (int *)bmalloc(*numRecvPts * sizeof(int));
  sm = (double *)bmalloc(*numRecvPts * sizeof(double));
  lg = (double *)bmalloc(*numRecvPts * sizeof(double));

  send_recv_corank_data(corank, sm, lg, *numRecvPts, recvProc, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  { // print the path number
    fprintf(OUT, "Path number: %d (ID: %d)\n", (*recvPts)[i].pathNum, recvProc);

    if (T->MPType == 0)
    { // store the endpoint using _d
      store_cascade_endPoint(&(*recvPts)[i], corank[i], sm[i], lg[i], trackCount, T, OUT, RAWOUT, FAIL, CD, CD);
    }
    else if (T->MPType == 1)
    { // store the endpoint using _mp
      store_cascade_endPoint(&(*recvPts)[i], corank[i], sm[i], lg[i], trackCount, T, OUT, RAWOUT, FAIL, CD, CD);
    }
    else
    { // store the endpoint using AMP
      store_cascade_endPoint(&(*recvPts)[i], corank[i], sm[i], lg[i], trackCount, T, OUT, RAWOUT, FAIL, CD, CD);
    }
  }

  free(lg);
  free(sm);
  free(corank);

  return recvProc;
}

int cascade_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, rankType = 1, *corank = NULL;
  double *sm = NULL, *lg = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = &cascade_dehom;

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

int cascade_create_send_packet_sort(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i;
  cascade_t *CD = NULL;

  // setup CD
  if (MPType == 0 || MPType == 2)
    CD = (cascade_t *)ED_d;
  else
    CD = (cascade_t *)ED_mp;

  // create the packet
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    setup_cascade_startPoint(&sendPts[i], pathNum[startNum + i], MPType, CD, CD, 1);
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // clear memory
  CD = NULL;

  return size;
}

int cascade_recv_store_packet_sort(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, recvProc;

  // recv the packet
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
    store_cascade_sortPoint(&(*recvPts)[i], trackCount, T, OUT, RAWOUT, FAIL, ED_d, ED_mp);

  return recvProc;
}

int cascade_recv_sort_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recvs packet, tracks paths and send packet back        *
\***************************************************************/
{
  int i, rankDef = 0, finite = 1, soln = 1;
  cascade_t *CD = NULL;

  // setup CD
  if (T->MPType == 0 || T->MPType == 2)
    CD = (cascade_t *)ED_d;
  else
    CD = (cascade_t *)ED_mp;

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

  // make sure that we have paths to sort
  if (*numStartPts > 0)
  { // loop over each point
    for (i = 0; i < *numStartPts; i++)
    { // initialize rankDef & endPts
      (*endPts)[i].prec = (*startPts)[i].prec;
      rankDef = (*startPts)[i].codim > 0; // contains corank
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
        cascadeSortEndpoint(&rankDef, &finite, &soln, (*endPts)[i].condition_number, CD, CD->curr_codim_index, (*endPts)[i].pathNum, T, OUT, &(*endPts)[i].PD_d, &(*endPts)[i].PD_mp, (*endPts)[i].prec, (*endPts)[i].last_approx_d, (*endPts)[i].last_approx_mp, (*endPts)[i].last_approx_prec, change_cascade_prec);

        if (!finite)
        { // dehom point is infinite
          (*endPts)[i].retVal = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (rankDef)
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
  }

  // clear
  CD = NULL;

  return 0;
}

#endif
