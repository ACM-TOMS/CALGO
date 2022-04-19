// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

#ifdef _HAVE_MPI

void head_junkRemoval_codim(int minPacketSize, int maxPacketSize, membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int pathMod, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);

void head_junkRemoval_codim_loop(int minPacketSize, int maxPacketSize, int reducedOnly, int *isJunk, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int pathMod, int my_id, int num_processes, int headnode);

void head_pureDecomp_codim(int minPacketSize, int maxPacketSize, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode);

void head_calculateTraces(int minPacketSize, int maxPacketSize, int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProg, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int pathMod, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int my_id, int num_processes, int headnode);

void head_computeTraces(comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, int *isClassified, int num_paths, int minPacketSize, int maxPacketSize, int pathMod, int MPType, int Precision, int my_id, int headnode, int num_processes);

void head_componentDecomposition(int minPacketSize, int maxPacketSize, int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int trace_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode);

int head_monodromy(int minPacketSize, int maxPacketSize, witness_t *W, int codim_index, tracker_config_t *T, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec, int my_id, int num_processes, int headnode);

void head_sendMonodromyStructures(membership_slice_moving_t *sliceMover, int MPType, int maxPrec, int my_id, int headnode);

void worker_junkRemoval_codim(witness_t *W, int codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);

void worker_decomposition(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int my_id, int num_processes, int headnode);

void worker_pureDecomp_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int my_id, int num_processes, int headnode);

void worker_calculateTraces_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int my_id, int num_processes, int headnode);

int junkRemoval_create_send_packet(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int junkRemoval_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));

int trace_send_packet(int startNum, int size, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, int totalPaths, int sendProc);
int trace_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, int MPType);

void head_computeMonodromy(int minPacketSize, int maxPacketSize, int *loop_results, int *isClassified, int num_paths, int MPType, int Precision, int my_id, int headnode, int num_processes);
int monodromy_send_packet(int startNum, int size, endgame_data_t *sendPts, int *pathNum, int MPType, int sendProc);
int monodromy_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, int *loop_results, int MPType);

void worker_junkRemoval_check(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int my_id, int num_processes, int headnode);

void worker_componentDecomposition(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int my_id, int num_processes, int headnode) ;

void worker_monodromyLoop(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat, double *norms, int my_id, int num_processes, int headnode);

/////////////// JUNK REMOVAL ////////////////

void junkRemoval_par(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup fullRankProgs and endPts and then remove junk    *
* from W to make it a true witness set in parallel              *
\***************************************************************/
{
  int i, num_codim = W->num_codim, minPacketSize = 1, maxPacketSize = 10;
  FILE *OUT = NULL;
  char outName[] = "output_junkRemoval";

  // send T
  bcast_tracker_config_t(T, my_id, headnode);

  // send W
  bcast_witness_t(W, T->MPType, my_id, headnode);

  // send specificCodim & topDimension
  MPI_Bcast(&specificCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  MPI_Bcast(&topDimension, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // make sure that we should continue
  if (num_codim < 1)
  { // no codimensions for junk removal!
    return;
  }

  // seutp OUT
  OUT = fopen(outName, "w");

  // allocate for each codimension
  *sliceMover = (membership_slice_moving_t *)bmalloc(num_codim * sizeof(membership_slice_moving_t));
  *fullRankProgs = (prog_t ***)bmalloc(num_codim * sizeof(prog_t **));
  *fullRankProgInfo = (int **)bmalloc(num_codim * sizeof(int *));
  if (T->MPType == 0)
    *endPts_d = (endpoint_data_d **)bmalloc(num_codim * sizeof(endpoint_data_d *));
  else if (T->MPType == 1)
    *endPts_mp = (endpoint_data_mp **)bmalloc(num_codim * sizeof(endpoint_data_mp *));
  else
    *endPts_amp = (endpoint_data_amp **)bmalloc(num_codim * sizeof(endpoint_data_amp *));

  // setup sliceMover for each codimension and send them
  for (i = 0; i < num_codim; i++)
  {
    basic_setup_slice_moving(&(*sliceMover)[i], W, i, T->MPType, T->AMP_max_prec);
    bcast_membership_slice_moving_t(&(*sliceMover)[i], 0, T->MPType, my_id, headnode);
  }

  // loop through each codimension removing the junk points
  for (i = 0; i < num_codim; i++)
  { // loop through the singular points determining which ones are junk
    head_junkRemoval_codim(minPacketSize, maxPacketSize, sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, i, T, OUT, pathMod, reducedOnly, specificCodim, topDimension, my_id, num_processes, headnode);
  }

  // print the witness set chart
  witnessSetOutputChart(W, stdout, T->MPType, reducedOnly);

  // close OUT
  fclose(OUT);

  return;
}

void head_junkRemoval_codim(int minPacketSize, int maxPacketSize, membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int pathMod, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: parallel junk removal for the given codimension        *
\***************************************************************/
{
  int i, j, num_sing = W->codim[codim_index].num_sing;
  int *isJunk = (int *)bmalloc(num_sing * sizeof(int));

  // setup isJunk - non-singular endpoints are always not junk so we only need to check the ones that are singular
  for (i = 0; i < num_sing; i++)
  { // initialize
    if (reducedOnly)
      isJunk[i] = 1;
    else
      isJunk[i] = 0;
  }

  if (codim_index > 0 || (codim_index == 0 && specificCodim > 0 && topDimension == 0)) // only need to remove junk if this is not the top codim
  { // loop through the singular endpoints and check if they are junk
    head_junkRemoval_codim_loop(minPacketSize, maxPacketSize, reducedOnly, isJunk, W, codim_index, T, OUT, pathMod, my_id, num_processes, headnode);
  }

  if (specificCodim > 0)
  { // count the number that failed
    j = 0;
    for (i = 0; i < num_sing; i++)
      if (isJunk[i] < 0)
      { // increment and say it is junk
        j++;
        isJunk[i] = 1;
      }

      if (j > 0)
      {
        printf("\nNOTE: The local dimension test was terminated before it could determine if %d point%s %s isolated.\n", j, j == 1 ? "" : "s", j == 1 ? "was" : "were");
        printf("      Consider either increasing MaxLDTDepth or using MaxCodimension instead of SpecificCodimension.\n");
      }
  }

  // remove the junk points and resend the dimension to the workers
  remove_junk_points(W, codim_index, T->MPType, isJunk);
  bcast_witnessCodim_t(&W->codim[codim_index], T->MPType, W->curr_precision, my_id, headnode);

  // setup everything at this codim are witness points
  W->codim[codim_index].deflations_needed = (int *)bmalloc(W->codim[codim_index].num_set * sizeof(int));
  if (T->MPType == 0)
    deflate_for_junkRemoval(&(*fullRankProgs)[codim_index], &(*fullRankProgInfo)[codim_index], &(*endPts_d)[codim_index], NULL, NULL, &(*sliceMover)[codim_index], W, codim_index, -1, T, OUT);
  else if (T->MPType == 1)
    deflate_for_junkRemoval(&(*fullRankProgs)[codim_index], &(*fullRankProgInfo)[codim_index], NULL, &(*endPts_mp)[codim_index], NULL, &(*sliceMover)[codim_index], W, codim_index, -1, T, OUT);
  else
    deflate_for_junkRemoval(&(*fullRankProgs)[codim_index], &(*fullRankProgInfo)[codim_index], NULL, NULL, &(*endPts_amp)[codim_index], &(*sliceMover)[codim_index], W, codim_index, -1, T, OUT);

  // look to sharpen the witness points, if needed
  if (T->sharpenDigits > 0)
  { // sharpen the witness points
    if (T->MPType == 0)
      sharpen_deflated(W, codim_index, (*fullRankProgs)[codim_index], (*fullRankProgInfo)[codim_index], &(*sliceMover)[codim_index], (*endPts_d)[codim_index], NULL, NULL, T, OUT);
    else if (T->MPType == 1)
      sharpen_deflated(W, codim_index, (*fullRankProgs)[codim_index], (*fullRankProgInfo)[codim_index], &(*sliceMover)[codim_index], NULL, (*endPts_mp)[codim_index], NULL, T, OUT);
    else
      sharpen_deflated(W, codim_index, (*fullRankProgs)[codim_index], (*fullRankProgInfo)[codim_index], &(*sliceMover)[codim_index], NULL, NULL, (*endPts_amp)[codim_index], T, OUT);
  }

  // send this new information to the workers
  bcast_witness_structures(&(*fullRankProgs)[codim_index], &(*fullRankProgInfo)[codim_index], &(*endPts_d)[codim_index], &(*endPts_mp)[codim_index], &(*endPts_amp)[codim_index], T->MPType, W, codim_index, my_id, headnode);

  // free isJunk
  free(isJunk);

  return;
}

void head_junkRemoval_codim_loop(int minPacketSize, int maxPacketSize, int reducedOnly, int *isJunk, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: parallel junk removal for the given codimension        *
\***************************************************************/
{
  trackingStats trackCount;
  int (*change_prec)(void const *, int) = NULL;
  int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int) = NULL;
  int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)) = NULL;

  // setup the pointers
  change_prec = NULL;
  create_send_packet = &junkRemoval_create_send_packet;
  recv_store_packet = &junkRemoval_recv_store_packet;

  // setup the current codimension
  W->curr_codim_index = codim_index;

  // display messages
  printf("\nRemoving junk points from codimension %d: %d endpoint%s to check.\n", W->codim[codim_index].codim, W->codim[codim_index].num_sing, W->codim[codim_index].num_sing == 1 ? "" : "s");
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Removing junk points from codimension %d.\n", W->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  // initialize trackCount - not used
  init_trackingStats(&trackCount);

  // do the actual parallel junk removal
  if (reducedOnly) // all singular points are already set to be removed
    head_trackPaths2("Checking", 0, 0, minPacketSize, maxPacketSize, &trackCount, pathMod, T, W, W, change_prec, NULL, OUT, NULL, NULL, NULL, NULL, isJunk, my_id, headnode, num_processes, create_send_packet, recv_store_packet);
  else // do normal junk removal
    head_trackPaths2("Checking", 0, W->codim[codim_index].num_sing, minPacketSize, maxPacketSize, &trackCount, pathMod, T, W, W, change_prec, NULL, OUT, NULL, NULL, NULL, NULL, isJunk, my_id, headnode, num_processes, create_send_packet, recv_store_packet);

  return;
}

/////////// IRREDUCIBLE DECOMPOSITION ///////////////

void pureDecomp_par(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decomposes all pure-dim witness sets into irred comp   *
* using parallel                                                *
\***************************************************************/
{
  int i, num_codim = W->num_codim, minPacketSize = 1, maxPacketSize = 10;
  FILE *OUT = NULL;
  char outName[] = "output_decomp", midName[] = "midpath_data";

  // setup OUT
  OUT = fopen(outName, "w");

  // loop over the codimensions
  for (i = 0; i < num_codim; i++)
    if (W->codim[i].num_set > 0)
    { // decompose codim_index 'i'
      if (T->MPType == 0)
        head_pureDecomp_codim(minPacketSize, maxPacketSize, &sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], endPts_d[i], NULL, NULL, W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
      else if (T->MPType == 1)
        head_pureDecomp_codim(minPacketSize, maxPacketSize, &sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], NULL, endPts_mp[i], NULL, W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
      else
        head_pureDecomp_codim(minPacketSize, maxPacketSize, &sliceMover[i], fullRankProgs[i], fullRankProgsInfo[i], NULL, NULL, endPts_amp[i], W, i, T, OUT, midName, pathMod, my_id, num_processes, headnode);
    }

  // close OUT
  fclose(OUT);

  return;
}

void head_pureDecomp_codim(int minPacketSize, int maxPacketSize, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decomposes codim 'codim_index' into irred components   *
\***************************************************************/
{
  int i, amp_trace_prec, num_points = W->codim[codim_index].num_set;
  int *traces_prec = T->MPType == 2 ? (int *)bmalloc(num_points * sizeof(int)) : NULL, *isClassified = (int *)bmalloc(num_points * sizeof(int));
  comp_d gamma_d, s_d[2], *traces_d = (T->MPType == 0 || T->MPType == 2) ? (comp_d *)bmalloc(num_points * sizeof(comp_d)) : NULL;
  comp_mp gamma_mp, s_mp[2], *traces_mp = (T->MPType == 1 || T->MPType == 2) ? (comp_mp *)bmalloc(num_points * sizeof(comp_mp)) : NULL;
  vec_d proj_d, sliceVec_d;
  vec_mp proj_mp, sliceVec_mp;
  mpq_t gamma_rat[2], s_rat[2][2], **proj_rat = NULL, **sliceVec_rat = NULL;
  FILE *MIDOUT = fopen(midName, "w");

  // allocate for component_nums
  W->codim[codim_index].num_components = 0;
  W->codim[codim_index].component_nums = (int *)bmalloc(num_points * sizeof(int));

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
    // setup gamma_d
    get_comp_rand_d(gamma_d);
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
    // setup gamma_mp
    init_mp(gamma_mp);
    get_comp_rand_mp(gamma_mp);
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

    // setup gamma_d, gamma_mp & gamma_rat
    get_comp_rand_rat(gamma_d, gamma_mp, gamma_rat, sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // compute traces for this codimension
  if (T->MPType == 0)
  { // compute traces in double precision
    head_calculateTraces(minPacketSize, maxPacketSize, isClassified, traces_d, NULL, NULL, sliceMover, fullRankProgs, fullRankProgsInfo, endPts_d, NULL, NULL, W, codim_index, T, OUT, MIDOUT, pathMod, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, my_id, num_processes, headnode);
  }
  else if (T->MPType == 1)
  { // compute traces in fixed multi precision
    head_calculateTraces(minPacketSize, maxPacketSize, isClassified, NULL, traces_mp, NULL, sliceMover, fullRankProgs, fullRankProgsInfo, NULL, endPts_mp, NULL, W, codim_index, T, OUT, MIDOUT, pathMod, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, my_id, num_processes, headnode);
  }
  else
  { // compute traces using adaptive precision
    head_calculateTraces(minPacketSize, maxPacketSize, isClassified, traces_d, traces_mp, traces_prec, sliceMover, fullRankProgs, fullRankProgsInfo, NULL, NULL, endPts_amp, W, codim_index, T, OUT, MIDOUT, pathMod, proj_d, proj_mp, proj_rat, sliceVec_d, sliceVec_mp, sliceVec_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, my_id, num_processes, headnode);
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

  // do the actual component classification
  head_componentDecomposition(minPacketSize, maxPacketSize, isClassified, traces_d, traces_mp, amp_trace_prec, sliceMover, fullRankProgs, endPts_d, endPts_mp, endPts_amp, W, codim_index, T, OUT, midName, pathMod, my_id, num_processes, headnode);

  // clear memory
  clear_vec(proj_d, proj_mp, proj_rat, T->MPType);
  clear_vec(sliceVec_d, sliceVec_mp, sliceVec_rat, T->MPType);
  clear_d_mp_rat(gamma_d, gamma_mp, gamma_rat, T->MPType);
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

void head_calculateTraces(int minPacketSize, int maxPacketSize, int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProg, int *fullRankProgsInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int pathMod, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the traces for 'codim_index'                   *
\***************************************************************/
{
  int i, num_traces = 0, num_set = W->codim[codim_index].num_set;

  // calculate the number of traces that need to be computed and setup memory
  for (i = 0; i < num_set; i++)
  { // setup memory
    isClassified[i] = 0;
    W->codim[codim_index].component_nums[i] = -1;

    // initialize traces
    if (T->MPType == 0)
    {
      set_zero_d(traces_d[i]);
    }
    if (T->MPType == 1)
    {
      init_mp(traces_mp[i]);
      set_zero_mp(traces_mp[i]);
    }
    else if (T->MPType == 2)
    {
      traces_prec[i] = 52;
      set_zero_d(traces_d[i]);
    }

    // see if deflation worked
    if (fullRankProgsInfo[i] == -1)
    { // deflation did not work - set isClassified to 1
      isClassified[i] = 1;
    }
    else
    { // increment the number of traces that need to be found
      num_traces++;
    }
  }

  // send the data
  bcast_trace_structures(proj_d, proj_mp, &proj_rat, v_d, v_mp, &v_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, T->MPType, sliceMover->curr_precision, my_id, headnode);

  // setup the current codimension
  W->curr_codim_index = codim_index;

  // display messages
  printf("\nCalculating traces for codimension %d.\n", W->codim[codim_index].codim);
  fprintf(OUT, "\n*****************************************************\n");
  fprintf(OUT, "Calculating traces for codimension %d.\n", W->codim[codim_index].codim);
  fprintf(OUT, "*****************************************************\n");

  head_computeTraces(traces_d, traces_mp, traces_prec, isClassified, num_traces, minPacketSize, maxPacketSize, pathMod, T->MPType, T->Precision, my_id, headnode, num_processes);

  return;
}

void head_computeTraces(comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, int *isClassified, int num_paths, int minPacketSize, int maxPacketSize, int pathMod, int MPType, int Precision, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computing the traces                                   *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, maxSize;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum
  j = 0;
  for (i = 0; i < num_paths; i++)
  { // find the next one that is not classified
    while (isClassified[j])
      j++;

    // save this one
    pathNum[i] = j;

    // move past it
    j++;
  }

  // allocate & intialize sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], Precision);
  }
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create & send the packet - returns the size of the packet
      lastSize[i] = trace_send_packet(count, packetSizes[i], sendPts, pathNum, MPType, pathMod, num_paths, i);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = trace_recv_store_packet(&recvPts, &numRecvPts, traces_d, traces_mp, traces_prec, MPType);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create & send the next packet - returns the size of the packet
    lastSize[recvProc] = trace_send_packet(count, packetSizes[recvProc], sendPts, pathNum, MPType, pathMod, num_paths, recvProc);

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
    recvProc = trace_recv_store_packet(&recvPts, &numRecvPts, traces_d, traces_mp, traces_prec, MPType);

    // tell the worker that this level is complete
    packetSizes[recvProc] = 0;
    lastSize[recvProc] = trace_send_packet(count, packetSizes[recvProc], sendPts, pathNum, MPType, pathMod, num_paths, recvProc);

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

void head_componentDecomposition(int minPacketSize, int maxPacketSize, int *isClassified, comp_d *traces_d, comp_mp *traces_mp, int trace_prec, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, char *midName, int pathMod, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: break up into irreducible components                   *
*   Assume that all traces are in the same precision            *
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
    printf("\nUsing mondoromy to decompose codimension %d.\n", W->codim[codim_index].codim);
    fprintf(OUT, "\n*****************************************************\n");
    fprintf(OUT, "Using monodromy to decompose codimension %d.\n", W->codim[codim_index].codim);
    fprintf(OUT, "*****************************************************\n");

    // do monodromy to find connections
    its = num_no_connections = 0;
    while (num_to_classify > 0 && num_temp_gps > T->max_num_pts_for_trace && its < T->max_num_mon_linears && num_no_connections < T->max_num_bad_loops_in_mon)
    { // setup and send the monodromy structures 
      head_sendMonodromyStructures(sliceMover, T->MPType, T->AMP_max_prec, my_id, headnode);

      // perform a monodromy loop and update the structures if any connections were made
      retVal = head_monodromy(minPacketSize, maxPacketSize, W, codim_index, T, &num_to_classify, &isClassified, &num_temp_gps, &temp_gp_nums, &temp_gp_trace_d, &temp_gp_trace_mp, trace_prec, my_id, num_processes, headnode);

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
  }

  // tell workers that monodromy is over with
  i = 0;
  bcast_monodromy_structures(&i, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, T->MPType, T->Precision, my_id, headnode);

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

int head_monodromy(int minPacketSize, int maxPacketSize, witness_t *W, int codim_index, tracker_config_t *T, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of new connections                      *
* NOTES: do monodromy to decompose a witness set                *
\***************************************************************/
{
  int num_connections = 0, *loop_results = (int *)bmalloc(*num_to_classify * sizeof(int));
   
  // setup the current codimension
  W->curr_codim_index = codim_index;
 
  // display messages
  printf("Performing monodromy loops: %d point%s left to classify\n", *num_to_classify, *num_to_classify == 1 ? "" : "s");

  // do the actual parallel monodromy looops
  head_computeMonodromy(minPacketSize, maxPacketSize, loop_results, *isClassified, *num_to_classify, T->MPType, T->Precision, my_id, headnode, num_processes);

  // update the temporary groups and find the number of new connections
  num_connections = updateTempGroupsFromLoops(W->codim[codim_index].num_set, *num_to_classify, *isClassified, loop_results, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);

  if (num_connections > 0)
  { // check for completed components and update
    checkCompleteComponents(W, codim_index, T->final_tol_times_mult, num_to_classify, isClassified, num_temp_gps, temp_gp_nums, temp_gp_trace_d, temp_gp_trace_mp, trace_prec);
  }
 
  // clear memory
  free(loop_results);

  return num_connections;
}

void head_sendMonodromyStructures(membership_slice_moving_t *sliceMover, int MPType, int maxPrec, int my_id, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup and send the structures for monodromy            *
\***************************************************************/
{
  int i, cont = 1, size;

  if (MPType == 0)
  { // send _d
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

    bcast_monodromy_structures(&cont, v_out_d, NULL, NULL, v_in_d, NULL, NULL, gamma_out_d, NULL, NULL, gamma_in_d, NULL, NULL, MPType, 52, my_id, headnode);

    // clear memory
    clear_vec_d(v_out_d);
    clear_vec_d(v_in_d);
  } 
  else if (MPType == 1)
  { // send _mp
    comp_mp gamma_out_mp, gamma_in_mp;
    vec_mp v_out_mp, v_in_mp;

    // setup size
    size = sliceMover->B_mp->rows;

    // initialize
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

    bcast_monodromy_structures(&cont, NULL, v_out_mp, NULL, NULL, v_in_mp, NULL, NULL, gamma_out_mp, NULL, NULL, gamma_in_mp, NULL, MPType, sliceMover->curr_precision, my_id, headnode);

    // clear memory
    clear_mp(gamma_out_mp);
    clear_mp(gamma_in_mp);
    clear_vec_mp(v_out_mp);
    clear_vec_mp(v_in_mp);
  }
  else
  { // send _rat
    comp_d gamma_out_d, gamma_in_d;
    comp_mp gamma_out_mp, gamma_in_mp;
    mpq_t gamma_out_rat[2], gamma_in_rat[2];
    vec_d v_out_d, v_in_d;
    vec_mp v_out_mp, v_in_mp;
    mpq_t **v_out_rat, **v_in_rat = NULL;

    // setup gammas
    get_comp_rand_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, sliceMover->curr_precision, maxPrec, 1, 1);
    get_comp_rand_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, sliceMover->curr_precision, maxPrec, 1, 1);
    // setup size
    size = sliceMover->B_d->rows;

    // setup v's
    init_vec_d(v_out_d, size); init_vec_d(v_in_d, size);
    init_vec_mp2(v_out_mp, size, sliceMover->curr_precision); init_vec_mp2(v_in_mp, size, sliceMover->curr_precision);
    init_vec_rat(v_out_rat, size); init_vec_rat(v_in_rat, size);

    // out is random, in is zero
    make_vec_random_rat(v_out_d, v_out_mp, v_out_rat, size, sliceMover->curr_precision, maxPrec, 0, 0);
    v_in_d->size = v_in_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&v_in_d->coord[i]);
      set_zero_mp(&v_in_mp->coord[i]);
      set_zero_rat(v_in_rat[i]);
    }

    bcast_monodromy_structures(&cont, v_out_d, v_out_mp, &v_out_rat, v_in_d, v_in_mp, &v_in_rat, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, MPType, sliceMover->curr_precision, my_id, headnode);

    // clear memory
    clear_d_mp_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, MPType);
    clear_d_mp_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, MPType);
    clear_vec(v_out_d, v_out_mp, v_out_rat, MPType);
    clear_vec(v_in_d, v_in_mp, v_in_rat, MPType);
  }

  return;
}

/////////// WORKER FUNCTIONS ///////////////

void worker_witness_superset_decomposition(int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: worker function for decomposing a witness superset     *
\***************************************************************/
{
  witness_t witnessSet;
  tracker_config_t T;
  int specificCodim, topDimension, **fullRankProgInfo = NULL;
  prog_t ***fullRankProgs = NULL;
  membership_slice_moving_t *sliceMover = NULL;
  endpoint_data_d **endPts_d = NULL;
  endpoint_data_mp **endPts_mp = NULL;
  endpoint_data_amp **endPts_amp = NULL;

  // recv T
  bcast_tracker_config_t(&T, my_id, headnode);
  initMP(T.Precision);

  // recv witnessSet
  bcast_witness_t(&witnessSet, T.MPType, my_id, headnode);

  // recv specificCodim & topDimension
  MPI_Bcast(&specificCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  MPI_Bcast(&topDimension, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // do junk removal in parallel
  worker_junkRemoval(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, &T, specificCodim, topDimension, my_id, num_processes, headnode);

  // decompose each witness set into irreducible pieces
  worker_decomposition(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, &witnessSet, &T, my_id, num_processes, headnode);

  // clear
  clear_sliceMover_fullRankProgs(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, T.MPType);
  witness_clear(&witnessSet, T.MPType);
  tracker_config_clear(&T);
  clearMP();

  return;
}

void worker_junkRemoval(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do junk removal in parallel                            *
\***************************************************************/
{
  int i, size, num_codim = W->num_codim;
  char *str = NULL, *midStr = NULL, outName[] = "output_junkRemoval", midName[] = "midpath_data";
  FILE *OUT = NULL;

  // make sure that we should continue
  if (num_codim < 1)
  { // no codimensions for junk removal!
    return;
  }

  // setup OUT
  size = 1 + snprintf(NULL, 0, "%s_%d", outName, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "%s_%d", outName, my_id);
  OUT = fopen(str, "w");

  // setup midStr
  size = 1 + snprintf(NULL, 0, "%s_%d", midName, my_id);
  midStr = (char *)brealloc(midStr, size * sizeof(char));
  sprintf(midStr, "%s_%d", midName, my_id);

  // setup sliceMover
  *sliceMover = (membership_slice_moving_t *)bmalloc(num_codim * sizeof(membership_slice_moving_t));
  for (i = 0; i < num_codim; i++)
  {
    bcast_membership_slice_moving_t(&(*sliceMover)[i], 0, T->MPType, my_id, headnode);
  }

  // allocate other structures
  *fullRankProgs = (prog_t ***)bmalloc(num_codim * sizeof(prog_t **));
  *fullRankProgInfo = (int **)bmalloc(num_codim * sizeof(int *));
  if (T->MPType == 0)
    *endPts_d = (endpoint_data_d **)bmalloc(num_codim * sizeof(endpoint_data_d *));
  else if (T->MPType == 1)
    *endPts_mp = (endpoint_data_mp **)bmalloc(num_codim * sizeof(endpoint_data_mp *));
  else
    *endPts_amp = (endpoint_data_amp **)bmalloc(num_codim * sizeof(endpoint_data_amp *));

  // loop through the each codimensions removing the junk points
  for (i = 0; i < num_codim; i++)
  { // remove the junk points for this codimension
    worker_junkRemoval_codim(W, i, *sliceMover, *fullRankProgs, *fullRankProgInfo, *endPts_d, *endPts_mp, *endPts_amp, T, OUT, midStr, specificCodim, topDimension, my_id, num_processes, headnode);
  }

  // close file
  fclose(OUT);

  // clear memory
  free(str);
  free(midStr);

  return;
}

void worker_junkRemoval_codim(witness_t *W, int codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: parallel junk removal for the given codimension        *
\***************************************************************/
{ 
  // do the junk removal for this codimension
  W->curr_codim_index = codim_index;

  if (codim_index > 0 || (codim_index == 0 && specificCodim > 0 && topDimension == 0)) // only need to remove junk if this is not the top codim
  { // remove the junk points
    worker_junkRemoval_check(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, OUT, midName, specificCodim, my_id, num_processes, headnode);
  }

  // clear the codimension and recv the codimension without the junk points
  witness_clear_codim(&W->codim[codim_index], T->MPType);
  bcast_witnessCodim_t(&W->codim[codim_index], T->MPType, W->curr_precision, my_id, headnode);

  // setup the other structures with the new data
  bcast_witness_structures(&fullRankProgs[codim_index], &fullRankProgInfo[codim_index], &endPts_d[codim_index], &endPts_mp[codim_index], &endPts_amp[codim_index], T->MPType, W, codim_index, my_id, headnode);

  return;
}

void worker_decomposition(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do irreducible decomposition in parallel               *
\***************************************************************/
{
  int i, size, num_codim = W->num_codim;
  char *str = NULL, outName[] = "output_decomp";
  FILE *OUT = NULL;

  // setup OUT
  size = 1 + snprintf(NULL, 0, "%s_%d", outName, my_id);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "%s_%d", outName, my_id);
  OUT = fopen(str, "w");

  // loop through each codimension to decompose it into irreducible pieces
  for (i = 0; i < num_codim; i++)
    if (W->codim[i].num_set > 0)
    { // decompose codim_index 'i'
      if (T->MPType == 0)
        worker_pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgInfo[i], endPts_d[i], NULL, NULL, W, i, T, OUT, my_id, num_processes, headnode);
      else if (T->MPType == 1)
        worker_pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgInfo[i], NULL, endPts_mp[i], NULL, W, i, T, OUT, my_id, num_processes, headnode);
      else
        worker_pureDecomp_codim(&sliceMover[i], fullRankProgs[i], fullRankProgInfo[i], NULL, NULL, endPts_amp[i], W, i, T, OUT, my_id, num_processes, headnode);
    }

  // close file
  fclose(OUT);

  // clear memory
  free(str);

  return;
}

void worker_pureDecomp_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, int codim_index, tracker_config_t *T, FILE *OUT, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do irreducible decomposition in parallel               *
\***************************************************************/
{
  int size;
  comp_d gamma_d, s_d[2];
  comp_mp gamma_mp, s_mp[2];
  vec_d proj_d, v_d;
  vec_mp proj_mp, v_mp;
  mpq_t gamma_rat[2], s_rat[2][2], **proj_rat = NULL, **v_rat = NULL;
  char *midStr = NULL, midName[] = "midpath_data";
  FILE *MIDOUT = NULL;

  // setup codim
  W->curr_codim_index = codim_index;

  // setup midStr & open MIDOUT
  size = 1 + snprintf(NULL, 0, "%s_%d", midName, my_id);
  midStr = (char *)brealloc(midStr, size * sizeof(char));
  sprintf(midStr, "%s_%d", midName, my_id);
  MIDOUT = fopen(midName, "w");

  // recv the structures
  bcast_trace_structures(proj_d, proj_mp, &proj_rat, v_d, v_mp, &v_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, T->MPType, sliceMover->curr_precision, my_id, headnode);

  // compute the traces
  worker_calculateTraces_codim(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, v_d, v_mp, v_rat, s_d, s_mp, s_rat, gamma_d, gamma_mp, gamma_rat, my_id, num_processes, headnode);

  // close file
  fclose(MIDOUT);

  // reopen MIDOUT
  MIDOUT = fopen(midName, "w");

  // do the actual component classification
  worker_componentDecomposition(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, OUT, MIDOUT, my_id, num_processes, headnode);
  
  // close file
  fclose(MIDOUT);

  // clear memory
  clear_vec(proj_d, proj_mp, proj_rat, T->MPType);
  clear_vec(v_d, v_mp, v_rat, T->MPType);
  clear_d_mp_rat(gamma_d, gamma_mp, gamma_rat, T->MPType);
  clear_d_mp_rat(s_d[0], s_mp[0], s_rat[0], T->MPType);
  clear_d_mp_rat(s_d[1], s_mp[1], s_rat[1], T->MPType);
  free(midStr);

  return;
}

/////////////////// CREATE & SEND, RECV & STORE //////////////////

int junkRemoval_create_send_packet(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i, j, codim_index, *singEndpoints = (int *)bmalloc(size * sizeof(int));
  witness_t *witSet = NULL;

  // setup witSet
  if (MPType == 0 || MPType == 2)
    witSet = (witness_t *)ED_d;
  else
    witSet = (witness_t *)ED_mp;

  // setup the current codim index
  codim_index = witSet->curr_codim_index;

  // move past the other path numbers
  j = 0;
  for (i = 0; i < startNum; i++)
  { // find the index of the jth singular endpoint
    while (witSet->codim[codim_index].witnessPt_types[j] != SINGULAR)
      j++;

    // increment j
    j++;
  }
  // find the path numbers
  for (i = 0; i < size; i++)
  { // find the index of the jth singular endpoint
    while (witSet->codim[codim_index].witnessPt_types[j] != SINGULAR)
      j++;

    singEndpoints[i] = j;

    // increment j
    j++;
  }

  // create the packet
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("%s path %d of %d\n", jobName, startNum + i, totalPaths);

    // find the path number
    sendPts[i].pathNum = singEndpoints[i]; 
    // store the 'singular number'
    sendPts[i].codim = startNum + i;
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  // clear memory
  witSet = NULL;
  free(singEndpoints);

  return size;
}

int junkRemoval_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, recvProc;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, T->MPType, MPI_ANY_SOURCE, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  { // store the value
    rV[(*recvPts)[i].codim] = (*recvPts)[i].retVal;
    // print the footer to OUT 
    fprintf(OUT, "Path number: %d isJunk: %d\n", (*recvPts)[i].pathNum, rV[(*recvPts)[i].codim]);
  }

  return recvProc;
}

int trace_send_packet(int startNum, int size, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, int totalPaths, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i;

  // create the packet
  for (i = 0; i < size; i++)
  { // print the path number if needed
    if (pathMod > 0 && !((startNum + i) % pathMod))
      printf("Calculating %d of %d\n", startNum + i, totalPaths);

    // find the path number
    sendPts[i].pathNum = pathNum[startNum + i];
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  return size;
}

int trace_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, comp_d *traces_d, comp_mp *traces_mp, int *traces_prec, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, recvProc, pathNum;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, MPType, MPI_ANY_SOURCE, 0);

  // store the data
  if (MPType == 0)
  { // store to _d
    for (i = 0; i < *numRecvPts; i++)
    { // store the values
      pathNum = (*recvPts)[i].pathNum;

      traces_d[pathNum]->r = (*recvPts)[i].function_residual_d;
      traces_d[pathNum]->i = (*recvPts)[i].latest_newton_residual_d;
    }
  }
  else if (MPType == 1)
  { // store to _mp
    for (i = 0; i < *numRecvPts; i++)
    { // store the values
      pathNum = (*recvPts)[i].pathNum;

      mpf_set(traces_mp[pathNum]->r, (*recvPts)[i].function_residual_mp);
      mpf_set(traces_mp[pathNum]->i, (*recvPts)[i].latest_newton_residual_mp);
    }
  }
  else
  { // store prec & setup either _d or _mp
    for (i = 0; i < *numRecvPts; i++)
    { // store the values
      pathNum = (*recvPts)[i].pathNum;

      traces_prec[pathNum] = (*recvPts)[i].prec;

      if (traces_prec[pathNum] < 64)
      { // setup traces_d
        traces_d[pathNum]->r = (*recvPts)[i].function_residual_d;
        traces_d[pathNum]->i = (*recvPts)[i].latest_newton_residual_d;
      }
      else
      { // setup traces_mp
        init_mp2(traces_mp[pathNum], traces_prec[pathNum]);
       
        mpf_set(traces_mp[pathNum]->r, (*recvPts)[i].function_residual_mp);
        mpf_set(traces_mp[pathNum]->i, (*recvPts)[i].latest_newton_residual_mp);
      }
    }
  }

  return recvProc;
}

int monodromy_send_packet(int startNum, int size, endgame_data_t *sendPts, int *pathNum, int MPType, int sendProc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of packet actually sent                   *
* NOTES: creates and sends the packet for regen tracking        *
\***************************************************************/
{
  int i;

  // create the packet
  for (i = 0; i < size; i++)
  { // find the path number
    sendPts[i].pathNum = pathNum[startNum + i];
    // store the location
    sendPts[i].codim = startNum + i;
  }

  // send sendPts to 'sendProc'
  send_recv_endgame_data_t(&sendPts, &size, MPType, sendProc, 1);

  return size;
}

int monodromy_recv_store_packet(endgame_data_t **recvPts, int *numRecvPts, int *loop_results, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: processor that it recv'd the packet from       *
* NOTES: recv's packet and stores it correctly                  *
\***************************************************************/
{
  int i, recvProc, pathNum;

  // recv a packet back & who sent it
  recvProc = send_recv_endgame_data_t(recvPts, numRecvPts, MPType, MPI_ANY_SOURCE, 0);

  // store the data
  for (i = 0; i < *numRecvPts; i++)
  {
    pathNum = (*recvPts)[i].pathNum;
    loop_results[(*recvPts)[i].codim] = (*recvPts)[i].retVal;
  }

  return recvProc;
}

void head_computeMonodromy(int minPacketSize, int maxPacketSize, int *loop_results, int *isClassified, int num_paths, int MPType, int Precision, int my_id, int headnode, int num_processes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computing the monodromy                                *
\***************************************************************/
{
  int i, j, count, recvProc, numRecvPts = 0, maxSize;
  int *pathNum = (int *)bmalloc(num_paths * sizeof(int)), *packetSizes = (int *)bmalloc(num_processes * sizeof(int)), *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  endgame_data_t *sendPts = NULL, *recvPts = NULL;

  // setup packetSizes and determine maxSize
  packetSize_maker(packetSizes, num_processes, headnode, num_paths, headnode, minPacketSize, maxPacketSize, &maxSize);

  // setup pathNum
  j = 0;
  for (i = 0; i < num_paths; i++)
  { // find the next one that is not classified
    while (isClassified[j])
      j++;

    // save this one
    pathNum[i] = j;

    // move past it
    j++;
  }

  // allocate & intialize sendPts
  sendPts = (endgame_data_t *)bmalloc(maxSize * sizeof(endgame_data_t));
  for (i = 0; i < maxSize; i++)
  { // initialize
    init_endgame_data(&sendPts[i], Precision);
  }
  // initialize count
  count = 0;

  // send out the initial set of packets
  for (i = 0; i < num_processes; i++)
    if (i != headnode)
    { // create & send the packet - returns the size of the packet
      lastSize[i] = monodromy_send_packet(count, packetSizes[i], sendPts, pathNum, MPType, i);

      // update count
      count += lastSize[i];
    }

  // loop until all the paths have been sent out to the workers to be tracked
  while (count < num_paths)
  { // recv a packet back & who sent it
    recvProc = monodromy_recv_store_packet(&recvPts, &numRecvPts, loop_results, MPType);

    // find the size of the next packet
    packetSize_maker(packetSizes, num_processes, headnode, num_paths - count, recvProc, minPacketSize, maxPacketSize, &maxSize);

    // create & send the next packet - returns the size of the packet
    lastSize[recvProc] = monodromy_send_packet(count, packetSizes[recvProc], sendPts, pathNum, MPType, recvProc);

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
    recvProc = monodromy_recv_store_packet(&recvPts, &numRecvPts, loop_results, MPType);

    // tell the worker that this level is complete
    packetSizes[recvProc] = 0;
    lastSize[recvProc] = monodromy_send_packet(count, packetSizes[recvProc], sendPts, pathNum, MPType, recvProc);

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

void worker_junkRemoval_check(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: junk removal in parallel                               *
\***************************************************************/
{
  int i, codim_index = W->curr_codim_index, numStartPts = 0, numEndPts = 0;
  endgame_data_t *startPts = NULL, *endPts = NULL;

  do
  { // recv a packet, check the paths and send it back
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

    // make sure that we have paths to check
    if (numStartPts > 0)
    {
      for (i = 0; i < numStartPts; i++)
      { // setup the ith end point
        endPts[i].pathNum = startPts[i].pathNum;
        endPts[i].codim = startPts[i].codim;

        // check the ith path
        if (T->junkRemovalTest == 0)
        { // use the standard membership test for junk removal
          endPts[i].retVal = junkRemoval_mem(W, startPts[i].pathNum, codim_index, sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, T, OUT, midName, my_id, num_processes, headnode);
        }
        else // (T->junkRemovalTest == 1
        { // use the local dimension test for junk removal
          endPts[i].retVal = junkRemoval_ldt(W, startPts[i].pathNum, codim_index, sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, T, OUT, midName, specificCodim, my_id, num_processes, headnode);
        }
      }
  
      // send the packet back
      send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);
    } 
  } while (numStartPts > 0);

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_calculateTraces_codim(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int my_id, int num_processes, int headnode) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute traces in parallel                             *
\***************************************************************/
{
  int i, pathNum, codim_index = W->curr_codim_index, numStartPts = 0, numEndPts = 0;
  endgame_data_t *startPts = NULL, *endPts = NULL;
  int trace_prec = 52;
  comp_d trace_d;
  comp_mp trace_mp;

  if (T->MPType == 1)
    init_mp(trace_mp);

  // setup gamma
  if (T->MPType == 0 || T->MPType == 2)
    set_d(sliceMover->gamma_d, gamma_d);
  if (T->MPType == 1 || T->MPType == 2)
    set_mp(sliceMover->gamma_mp, gamma_mp);
  if (T->MPType == 2)
    set_rat(sliceMover->gamma_rat, gamma_rat);
  
  do
  { // recv a packet, check the paths and send it back
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

    // make sure that we have paths to check
    if (numStartPts > 0)
    {
      for (i = 0; i < numStartPts; i++)
      { // setup the ith end point
        pathNum = endPts[i].pathNum = startPts[i].pathNum;
        // compute its trace
        if (T->MPType == 0)
        {
          calculateTrace(trace_d, trace_mp, &trace_prec, sliceMover, fullRankProgs[pathNum], &endPts_d[pathNum], NULL, NULL, W->codim[codim_index].codim, pathNum, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, v_d, v_mp, v_rat, s_d, s_mp, s_rat, 0);

          // setup function_residual_d, latest_newton_residual_d with trace_d
          endPts[i].prec = 52;
          endPts[i].function_residual_d = trace_d->r;
          endPts[i].latest_newton_residual_d = trace_d->i;
        }
        else if (T->MPType == 1)
        {
          calculateTrace(trace_d, trace_mp, &trace_prec, sliceMover, fullRankProgs[pathNum], NULL, &endPts_mp[pathNum], NULL, W->codim[codim_index].codim, pathNum, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, v_d, v_mp, v_rat, s_d, s_mp, s_rat, 0);

          // setup function_residual_mp, latest_newton_residual_mp with trace_mp
          endPts[i].prec = T->Precision;
          mpf_set(endPts[i].function_residual_mp, trace_mp->r);
          mpf_set(endPts[i].latest_newton_residual_mp, trace_mp->i);
        }
        else
        {
          calculateTrace(trace_d, trace_mp, &trace_prec, sliceMover, fullRankProgs[pathNum], NULL, NULL, &endPts_amp[pathNum], W->codim[codim_index].codim, pathNum, T, OUT, MIDOUT, proj_d, proj_mp, proj_rat, v_d, v_mp, v_rat, s_d, s_mp, s_rat, 0);

          if (trace_prec < 64)
          { // setup function_residual_d, latest_newton_residual_d with trace_d
            endPts[i].prec = 52;
            endPts[i].function_residual_d = trace_d->r;
            endPts[i].latest_newton_residual_d = trace_d->i;
          }
          else
          { // setup function_residual_mp, latest_newton_residual_mp with trace_mp
            endPts[i].prec = trace_prec;
            mpf_set_prec(endPts[i].function_residual_mp, trace_prec);
            mpf_set_prec(endPts[i].latest_newton_residual_mp, trace_prec);
            mpf_set(endPts[i].function_residual_mp, trace_mp->r);
            mpf_set(endPts[i].latest_newton_residual_mp, trace_mp->i);

            clear_mp(trace_mp);
          }  
        }
      }

      // send the packet back
      send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);
    }

  } while (numStartPts > 0);

  if (T->MPType == 1)
    clear_mp(trace_mp);

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  return;
}

void worker_componentDecomposition(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int my_id, int num_processes, int headnode) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform monodromy in parallel                          *
\***************************************************************/
{
  int i, contMonodromy = 0, codim_index = W->curr_codim_index;
  int num_points = W->codim[codim_index].num_set;
  double *norms = (double *)bmalloc(num_points * sizeof(double));
  comp_d gamma_out_d, gamma_in_d;
  comp_mp gamma_out_mp, gamma_in_mp;
  mpq_t gamma_out_rat[2], gamma_in_rat[2];
  vec_d v_out_d, v_in_d;
  vec_mp v_out_mp, v_in_mp;
  mpq_t **v_out_rat = NULL, **v_in_rat = NULL;

  // setup norms - a double approximation to the norms - helps to classify them faster
  if (T->MPType == 0)
  { // use Pts_d
    for (i = 0; i < num_points; i++)
      norms[i] = infNormVec_d(W->codim[codim_index].witnessPts_d[i].endPt);
  }
  else if (T->MPType == 1)
  { // use Pts_mp
    for (i = 0; i < num_points; i++)
      norms[i] = infNormVec_mp(W->codim[codim_index].witnessPts_mp[i].endPt);
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

  do
  { // recv the monodromy structures
    bcast_monodromy_structures(&contMonodromy, v_out_d, v_out_mp, &v_out_rat, v_in_d, v_in_mp, &v_in_rat, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, T->MPType, sliceMover->curr_precision, my_id, headnode);

    // see if we need to do another monodromy loop
    if (contMonodromy)
    { // perform another loop
      worker_monodromyLoop(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, OUT, MIDOUT, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, v_out_d, v_out_mp, v_out_rat, v_in_d, v_in_mp, v_in_rat, norms, my_id, num_processes, headnode);
    
      // clear the structures
      clear_d_mp_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, T->MPType);
      clear_d_mp_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, T->MPType);
      clear_vec(v_out_d, v_out_mp, v_out_rat, T->MPType);
      clear_vec(v_in_d, v_in_mp, v_in_rat, T->MPType);
    }

  } while (contMonodromy);

  // clear memory
  free(norms);

  return;
}

void worker_monodromyLoop(membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, witness_t *W, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat, double *norms, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform a monodromy loop in parallel                   *
\***************************************************************/
{
  int i, size, pathNum, codim_index = W->curr_codim_index, numStartPts = 0, numEndPts = 0;
  int numSet = W->codim[codim_index].num_set, *isClassified = NULL;
  endgame_data_t *startPts = NULL, *endPts = NULL;

  // setup size
  size = T->MPType == 1 ? sliceMover->B_mp->rows : sliceMover->B_d->rows;

  // setup isClassified
  isClassified = (int *)bmalloc(numSet * sizeof(int));
  for (i = 0; i < numSet; i++)
    isClassified[i] = 0;

  do
  { // recv a packet, check the paths and send it back
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

    // make sure that we have paths to check
    if (numStartPts > 0)
    {
      for (i = 0; i < numStartPts; i++)
      { // setup the ith end point
        pathNum = endPts[i].pathNum = startPts[i].pathNum;
        endPts[i].codim = startPts[i].codim;

        // setup the sliceMover
        initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
        final_setup_slice_moving(sliceMover, fullRankProgs[pathNum], T->MPType, T->AMP_max_prec, 0);

        // perform the monodromy loop on this point
        if (T->MPType == 0)
          endPts[i].retVal = monodromyLoop(sliceMover, fullRankProgs[pathNum], &endPts_d[pathNum], NULL, NULL, W, codim_index, pathNum, T, OUT, MIDOUT, T->final_tol_times_mult, isClassified, gamma_out_d, NULL, NULL, gamma_in_d, NULL, NULL, vec_out_d, NULL, NULL, vec_in_d, NULL, NULL, norms);
        else if (T->MPType == 1)
          endPts[i].retVal = monodromyLoop(sliceMover, fullRankProgs[pathNum], NULL, &endPts_mp[pathNum], NULL, W, codim_index, pathNum, T, OUT, MIDOUT, T->final_tol_times_mult, isClassified, NULL, gamma_out_mp, NULL, NULL, gamma_in_mp, NULL, NULL, vec_out_mp, NULL, NULL, vec_in_mp, NULL, norms);
        else
          endPts[i].retVal = monodromyLoop(sliceMover, fullRankProgs[pathNum], NULL, NULL, &endPts_amp[pathNum], W, codim_index, pathNum, T, OUT, MIDOUT, T->final_tol_times_mult, isClassified, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, vec_out_d, vec_out_mp, vec_out_rat, vec_in_d, vec_in_mp, vec_in_rat, norms);
      }

      // send the packet back
      send_recv_endgame_data_t(&endPts, &numEndPts, T->MPType, headnode, 1);
    }
  } while (numStartPts > 0);

  // startPts & endPts are cleared since numStartPts, numEndPts <=0

  // clear memory
  free(isClassified);

  return;
}

#endif

