// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "dimbydim.h"
#include "parallel.h"

// provides the functions to find the witness supersets using dimension-by-dimension

void dimbydim_seq_track(trackingStats *trackCount, FILE *OUT, char *outName, FILE *RAWOUT, char *rawName, FILE *MIDOUT, char *midName, FILE *FAIL, char *failName, int pathMod, double midpoint_tol, tracker_config_t *T, codim_t *CD);

void dimbydimTrackCodim(int max, int pathMod, codim_t CD[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, int codim_index, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void codimTrackPath(codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void dimbydimTrackCodim_trackBack(int max, int pathMod, codim_t CD[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, int codim_index, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void codimTrackPath_trackBack_found(int indexJ, trackBack_samples_t *EGsample, codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void codimTrackPath_trackBack(double trackBack_final_tol, trackBack_samples_t *EGsample, codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void dimbydimSortCodim(int max, int pathMod, codim_t CD_copy[], codim_t *CD, tracker_config_t T[], FILE **OUT, int codim_index, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));


void dimbydim_main(witness_t *witnessSuperset, int maxCodim, int specificCodim, tracker_config_t *T, int pathMod, double midpoint_tol, double intrinsicCutoffMultiplier, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does dim-by-dim tracking to find witness superset      *
\***************************************************************/
{
  codim_t CD;
  trackingStats trackCount;
  FILE *OUT, *RAWOUT, *MIDOUT, *FAIL;
  char outName[] = "output", rawName[] = "raw_data", midName[] = "midpath_data", failName[] = "failed_paths";

  // initialize trackCount
  init_trackingStats(&trackCount);

  // setup CD
  dimbydim_setup(&OUT, outName, &RAWOUT, rawName, &MIDOUT, midName, &FAIL, failName, T, &CD, "preproc_data", "deg.out", intrinsicCutoffMultiplier, maxCodim, specificCodim);

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
    fclose(MIDOUT); // not used
#ifdef _HAVE_MPI
    worker_info sendType;
    sendType.dataType = DIMBYDIM;
    bcast_worker_info(&sendType, my_id, headnode);

    dimbydim_par_track(&trackCount, OUT, RAWOUT, FAIL, midName, pathMod, midpoint_tol, T, &CD, my_id, num_processes, headnode);
#endif
  }
  else 
  { // do sequential tracking (including using OpenMP)
    dimbydim_seq_track(&trackCount, OUT, outName, RAWOUT, rawName, MIDOUT, midName, FAIL, failName, pathMod, midpoint_tol, T, &CD);
  }

  // print output chart
  dimbydimOutputChart(&CD, stdout, T->MPType);

  // copy the data in CD to witnessSuperset and clear CD
  dimbydim_copyWitness_clear(witnessSuperset, &CD, T->MPType, T->AMP_max_prec);

  return;
}

void dimbydim_seq_track(trackingStats *trackCount, FILE *OUT, char *outName, FILE *RAWOUT, char *rawName, FILE *MIDOUT, char *midName, FILE *FAIL, char *failName, int pathMod, double midpoint_tol, tracker_config_t *T, codim_t *CD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does dim-by-dim tracking in serial (or OpenMP)         *
\***************************************************************/
{
  int codim_index, num_paths, num_crossings = 0, max = max_threads();
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);

  // pointers for OpenMP tracking
  tracker_config_t *T_copy = NULL;
  codim_t *CD_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL;

  // setup the OpenMP structures
  setup_dimbydim_omp(max, &trackCount_copy, trackCount, &OUT_copy, OUT, outName, &RAWOUT_copy, RAWOUT, rawName, &MIDOUT_copy, MIDOUT, midName, &FAIL_copy, FAIL, failName, &T_copy, T, &CD_copy, CD);

  // setup the evaluators
  ptr_to_eval_d = &standard_dimbydim_eval_d;
  ptr_to_eval_mp = &standard_dimbydim_eval_mp;
  change_prec = &change_dimbydim_prec;
  find_dehom = &dimbydim_dehom;

  // find the witness supersets
  for (codim_index = 0; codim_index < CD->num_codim; codim_index++)
  { // find the number of paths that are tracking on this codim
    num_paths = CD->codim[codim_index].num_paths;

    // track the codimension pointed to by codim_index
    if (T->endgameNumber == 3)
    { // use the trackBack endgame
      dimbydimTrackCodim_trackBack(max, pathMod, CD_copy, T_copy, OUT_copy, RAWOUT_copy, MIDOUT_copy, FAIL_copy, codim_index, trackCount_copy, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
    }
    else
    { // use the standard endgame
      dimbydimTrackCodim(max, pathMod, CD_copy, T_copy, OUT_copy, RAWOUT_copy, MIDOUT_copy, FAIL_copy, codim_index, trackCount_copy, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
    }

    // combine all of the files
    combine_omp_file(OUT, &OUT_copy, max);
    combine_omp_file(MIDOUT, &MIDOUT_copy, max);
    combine_omp_file(RAWOUT, &RAWOUT_copy, max);
    combine_omp_file(FAIL, &FAIL_copy, max);

    // check for path crossings
    fclose(MIDOUT);
    num_crossings = 0;

    // check to see if using intrinsic slice
    if (CD->codim[codim_index].useIntrinsicSlice)
      midpoint_checker(num_paths, CD->codim[codim_index].codim, midpoint_tol, &num_crossings);
    else
      midpoint_checker(num_paths, CD->new_variables, midpoint_tol, &num_crossings);

    // setup OUT_copy for classifying
    setup_omp_file(&OUT_copy, OUT, outName, max);

    // sort the endpoints
    dimbydimSortCodim(max, pathMod, CD_copy, CD, T_copy, OUT_copy, codim_index, T->final_tol_times_mult, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

    // combine OUT
    combine_omp_file(OUT, &OUT_copy, max);

    if (codim_index + 1 < CD->num_codim)
    { // reopen all of the files since we have another codim to track
      MIDOUT = fopen(midName, "w");
      setup_omp_file(&MIDOUT_copy, MIDOUT, midName, max);
      setup_omp_file(&OUT_copy, OUT, outName, max);
      setup_omp_file(&RAWOUT_copy, RAWOUT, rawName, max);
      setup_omp_file(&FAIL_copy, FAIL, failName, max);
    }
    else
    { // close files
      fclose(OUT);
      fclose(RAWOUT);
      fclose(FAIL);
    }

    // clear the start points for this codim since they are no longer needed
    dimbydim_clear_start_points(CD, codim_index, T->MPType);
  }

  // clear the OpenMP structures
  clear_dimbydim_omp(max, &trackCount_copy, trackCount, outName, rawName, midName, failName, &T_copy, &CD_copy);

  return;
}

void setup_dimbydim_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, tracker_config_t **T_copy, tracker_config_t *T, codim_t **CD_copy, codim_t *CD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup everything needed to do dim-by-dim tracking      *
*  using OpenMP                                                 *
\***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int i;

  // error checking
  if (max_threads <= 0)
  {
    printf("\n\nERROR: The number of threads (%d) needs to be positive when setting up for tracking!\n", max_threads);
    bexit(ERROR_CONFIGURATION);
  }

  // setup the files
  setup_omp_file(OUT_copy, OUT, outName, max_threads);
  setup_omp_file(RAWOUT_copy, RAWOUT, rawName, max_threads);
  setup_omp_file(MIDOUT_copy, MIDOUT, midName, max_threads);
  setup_omp_file(FAIL_copy, FAIL, failName, max_threads);

  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *CD_copy = CD;
  }
  else // max_threads > 1
  { // allocate memory
    *trackCount_copy = (trackingStats *)bmalloc(max_threads * sizeof(trackingStats));
    *T_copy = (tracker_config_t *)bmalloc(max_threads * sizeof(tracker_config_t));
    *CD_copy = (codim_t *)bmalloc(max_threads * sizeof(codim_t));

    // copy T, CD & trackCount
    for (i = 0; i < max_threads; i++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[i], T);
      // copy CD
      setup_omp_codim_t(&(*CD_copy)[i], CD, T->MPType);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[i]);
      (*trackCount_copy)[i].numPoints = trackCount->numPoints;
    }
  }

  return;
}

void clear_dimbydim_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, tracker_config_t **T_copy, codim_t **CD_copy)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear everything needed to do dim-by-dim tracking      *
*  using OpenMP                                                 *
\***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int i;

  if (max_threads == 1)
  { // set the pointers to NULL since they just pointed to the actual values
    *trackCount_copy = NULL;
    *T_copy = NULL;
    *CD_copy = NULL;
  }
  else if (max_threads > 1)
  { // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);

    // clear the copies of T & ED_d
    for (i = max_threads - 1; i >= 0; i--)
    { // clear CD
      clear_omp_codim_t(&(*CD_copy)[i], (*T_copy)[i].MPType);
      // clear T_copy
      tracker_config_clear(&(*T_copy)[i]);
    }

    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*CD_copy);

    // delete the temporary files
    char *str = NULL;
    int size;
    for (i = 0; i < max_threads; i++)
    { // delete output
      size = 1 + snprintf(NULL, 0, "%s_%d", outName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", outName, i);
      remove(str);

      // delete rawout
      size = 1 + snprintf(NULL, 0, "%s_%d", rawName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawName, i);
      remove(str);

      // delete midout
      size = 1 + snprintf(NULL, 0, "%s_%d", midName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", midName, i);
      remove(str);

      // delete fail
      size = 1 + snprintf(NULL, 0, "%s_%d", failName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", failName, i);
      remove(str);
    }
  }

  return;
}

void dimbydimTrackCodim(int max, int pathMod, codim_t CD[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, int codim_index, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this codim-use OpenMP if available*
\***************************************************************/
{
  int i, oid, num_paths = CD[0].codim[codim_index].num_paths, codim = CD[0].codim[codim_index].codim, num_codim = CD[0].num_codim;

  // display messages
  printf("\nFinding a witness superset for codimension %d (%d of %d): %d path%s to track.\n", codim, codim_index, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Finding a witness superset for codimension %d.\n", codim);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this level
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid) shared(num_paths) schedule(runtime) if (num_paths > 0)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // track the path
    codimTrackPath(&CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  return;
}

void dimbydimTrackCodim_trackBack(int max, int pathMod, codim_t CD[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, int codim_index, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths for this codim                        *
\***************************************************************/
{
  int i, oid, num_paths = CD[0].codim[codim_index].num_paths, codim = CD[0].codim[codim_index].codim, num_codim = CD[0].num_codim;
  int retVal, indexI, indexJ, numSamples = 0;
  double trackBack_final_tol = 1e-2 * T->basicNewtonTol, *startPoint_norm_d = NULL;
  mpf_t *startPoint_norm_mp = NULL;
  trackBack_samples_t *EGsamples = NULL;

  // display messages
  printf("\nFinding a witness superset for codimension %d (%d of %d): %d path%s to track.\n", codim, codim_index, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Finding a witness superset for codimension %d.\n", codim);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this level
  if (T->MPType == 0)
  { // setup the norms
    startPoint_norm_d = (double *)bmalloc(num_paths * sizeof(double));
    for (i = 0; i < num_paths; i++)
      startPoint_norm_d[i] = infNormVec_d(CD->codim[codim_index].startPts_d[i]);

    // loop over the paths
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, startPoint_norm_d[i], NULL, CD->codim[codim_index].startPts_d[i], NULL, 52);
      if (retVal)
      { // the point was found
        codimTrackPath_trackBack_found(indexJ, &EGsamples[indexI], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], 52);

        // track the path
        codimTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }

    // clear memory
    free(startPoint_norm_d);
    for (i = 0; i < numSamples; i++)
      clear_trackBack_sample(&EGsamples[i]);
    free(EGsamples);
  }
  else if (T->MPType == 1)
  { // setup the norms
    startPoint_norm_mp = (mpf_t *)bmalloc(num_paths * sizeof(mpf_t));
    for (i = 0; i < num_paths; i++)
    {
      mpf_init(startPoint_norm_mp[i]);
      infNormVec_mp2(startPoint_norm_mp[i], CD->codim[codim_index].startPts_mp[i]);
    }

    // loop over the paths
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, 0, startPoint_norm_mp[i], NULL, CD->codim[codim_index].startPts_mp[i], T->Precision);
      if (retVal)
      { // the point was found
        codimTrackPath_trackBack_found(indexJ, &EGsamples[indexI], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], T->Precision);
         // track the path
        codimTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }

    // clear memory
    for (i = 0; i < num_paths; i++)
      mpf_clear(startPoint_norm_mp[i]);
    free(startPoint_norm_mp);
    for (i = 0; i < numSamples; i++)
      clear_trackBack_sample(&EGsamples[i]);
    free(EGsamples);
  }
  else
  { // use AMP
    startPoint_norm_d = (double *)bmalloc(num_paths * sizeof(double));
    for (i = 0; i < num_paths; i++)
      startPoint_norm_d[i] = infNormVec_d(CD->codim[codim_index].startPts_d[i]);

    // loop over the paths
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, num_paths);

      // determine if the start point is already known & find the indices if it is found
      retVal = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, startPoint_norm_d[i], NULL, CD->codim[codim_index].startPts_d[i], NULL, 52);
      if (retVal)
      { // the point was found
        codimTrackPath_trackBack_found(indexJ, &EGsamples[indexI], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec);
      }
      else
      { // increase the size of EGsamples
        EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));
        init_trackBack_sample(&EGsamples[numSamples], T->Precision);

        // track the path
        codimTrackPath_trackBack(trackBack_final_tol, &EGsamples[numSamples], &CD[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], codim_index, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

        // increment the size
        numSamples++;
      }
    }

    // clear memory
    free(startPoint_norm_d);
    for (i = 0; i < numSamples; i++)
      clear_trackBack_sample(&EGsamples[i]);
    free(EGsamples);
  }

  return;
}

void codimTrackPath_trackBack_found(int indexJ, trackBack_samples_t *EGsample, codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the path in codim_index & path_num of CD         *
\***************************************************************/
{
  int j, size, same_pathNum = EGsample->endPt.pathNum;
  EGsample->endPt.pathNum = path_num;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;

    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

    // print the mid point to MIDOUT
    size = EGsample->midPt_d[indexJ]->size;
    fprintf(MIDOUT, "%d\n", path_num);
    for (j = 0; j < size; j++)
      fprintf(MIDOUT, "%.15e %.15e\n", EGsample->midPt_d[indexJ]->coord[j].r, EGsample->midPt_d[indexJ]->coord[j].i);
    fprintf(MIDOUT, "\n");

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, CD->codim[codim_index].endPts_d[same_pathNum].corank, CD->codim[codim_index].endPts_d[same_pathNum].smallest_nonzero_SV, CD->codim[codim_index].endPts_d[same_pathNum].largest_zero_SV, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear
    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;

    init_point_data_mp(&startPt, 0);

    // setup for tracking the path
    point_cp_mp(startPt.point, CD->codim[codim_index].startPts_mp[path_num]);
    set_one_mp(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the path header for the point
    printPathHeader_mp(OUT, &startPt, T, path_num, CD, ptr_to_eval_mp);

    // print the mid point
    size = EGsample->midPt_mp[indexJ]->size;
    fprintf(MIDOUT, "%d\n", path_num);
    for (j = 0; j < size; j++)
    {
      print_mp(MIDOUT, 0, &EGsample->midPt_mp[indexJ]->coord[j]);
      fprintf(MIDOUT, "\n");
    }
    fprintf(MIDOUT, "\n");

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, CD->codim[codim_index].endPts_mp[same_pathNum].corank, CD->codim[codim_index].endPts_mp[same_pathNum].smallest_nonzero_SV, CD->codim[codim_index].endPts_mp[same_pathNum].largest_zero_SV, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    point_data_d startPt;

    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

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

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, CD->codim[codim_index].endPts_amp[same_pathNum].corank, CD->codim[codim_index].endPts_amp[same_pathNum].smallest_nonzero_SV, CD->codim[codim_index].endPts_amp[same_pathNum].largest_zero_SV, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_d(&startPt);
  }

  EGsample->endPt.pathNum = same_pathNum;

  return;
}

void codimTrackPath_trackBack(double trackBack_final_tol, trackBack_samples_t *EGsample, codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in codim_index & path_num of CD        *
\***************************************************************/
{
  int rankType = 1, corank;
  double smallest, largest;

  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;

    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, trackBack_final_tol, EGsample, &startPt, OUT, MIDOUT, T, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (EGsample->endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&EGsample->endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear
    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;

    // init MP
    init_point_data_mp(&startPt, 0);

    // setup for tracking the path
    point_cp_mp(startPt.point, CD->codim[codim_index].startPts_mp[path_num]);
    set_one_mp(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the path header for the point
    printPathHeader_mp(OUT, &startPt, T, path_num, CD, ptr_to_eval_mp);

    // track the path
    zero_dim_trackBack_path_rank_mp(path_num, rankType, NULL, &corank, &smallest, &largest, trackBack_final_tol, EGsample, &startPt, OUT, MIDOUT, T, CD, ptr_to_eval_mp, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (EGsample->endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&EGsample->endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    point_data_d startPt;

    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

    // track the path
    zero_dim_trackBack_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, trackBack_final_tol, EGsample, &startPt, OUT, MIDOUT, T, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (EGsample->endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&EGsample->endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&EGsample->endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_d(&startPt);
  }

  return;
}

void codimTrackPath(codim_t *CD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int codim_index, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in codim_index & path_num of CD        *
\***************************************************************/
{
  int rankType = 1, corank;
  double smallest, largest;
  endgame_data_t endPt;

  init_endgame_data(&endPt, T->Precision);

  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;
    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, &startPt, OUT, MIDOUT, T, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear
    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;

    // init MP
    init_point_data_mp(&startPt, 0);

    // setup for tracking the path
    point_cp_mp(startPt.point, CD->codim[codim_index].startPts_mp[path_num]);
    set_one_mp(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the path header for the point
    printPathHeader_mp(OUT, &startPt, T, path_num, CD, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_rank_mp(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, &startPt, OUT, MIDOUT, T, CD, ptr_to_eval_mp, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    point_data_d startPt;

    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, CD->codim[codim_index].startPts_d[path_num]);
    set_one_d(startPt.time);
    CD->curr_codim_index = codim_index;

    // print the header for the point
    printPathHeader_d(OUT, &startPt, T, path_num, CD, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(path_num, rankType, NULL, &corank, &smallest, &largest, &endPt, &startPt, OUT, MIDOUT, T, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // check to see if it should be sharpened (and can be)
    if (endPt.retVal == 0 && T->sharpenDigits > 0 && corank == 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&endPt, T, OUT, CD, CD, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // store the end point
    store_dimbydim_endPoint(&endPt, corank, smallest, largest, trackCount, T, OUT, RAWOUT, FAIL, CD, CD);

    // clear MP
    clear_point_data_d(&startPt);
  }

  clear_endgame_data(&endPt);

  return;
}

void dimbydimSortCodim(int max, int pathMod, codim_t CD_copy[], codim_t *CD, tracker_config_t T[], FILE **OUT, int codim_index, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoints found at the current codim         *
\***************************************************************/
/*
UNCLASSIFIED 0
NON_SINGULAR 10
SINGULAR 15
< 0 -> bad path (code = retVal returned from path)
*/
{
  int i, j, rankDef, soln, finite, oid, cont, indexI = 0, indexJ = 0, num_paths = CD->codim[codim_index].num_paths, codim = CD->codim[codim_index].codim, num_codim = CD->num_codim;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_nonsoln = 0;
 
  // print header for level
  printf("\nSorting codimension %d (%d of %d): %d path%s to sort.\n", codim, codim_index, num_codim, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Sorting codimension %d.\n", codim);
  fprintf(OUT[0], "*****************************************************\n");

  if (T[0].MPType == 0)
  { // sort in double precision
    point_data_d *PD_d = (point_data_d *)bmalloc(max * sizeof(point_data_d)); 
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
    vec_d tempVec;

    // initialize
    init_vec_d(tempVec, 0);
    for (i = 0; i < max; i++)
      init_point_data_d(&PD_d[i], 0);

    // determine if each path is (a) rank deficient (b) at infinity and (c) solution
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, soln, finite, oid) shared(num_paths) schedule(runtime) if (num_paths > 0)
#endif
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      rankDef = soln = finite = 0;

      // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // only check the ones that were successful
      if (CD_copy[oid].codim[codim_index].endPts_d[i].retVal == 0)
      { // setup for evaluation
        CD_copy[oid].curr_codim_index = codim_index;

        // setup PD_d[oid]
        point_cp_d(PD_d[oid].point, CD_copy[oid].codim[codim_index].endPts_d[i].endPt);
        set_zero_d(PD_d[oid].time);
  
        // determine if it is rank deficient - i.e. the corank is positive
        rankDef = CD_copy[oid].codim[codim_index].endPts_d[i].corank > 0;

        // do the classification
        dimbydimSortEndpoint(&rankDef, &finite, &soln, CD_copy[oid].codim[codim_index].endPts_d[i].cond_num, &CD_copy[oid], codim_index, i, &T[oid], OUT[oid], &PD_d[oid], NULL, 52, CD_copy[oid].codim[codim_index].endPts_d[i].last_approx, NULL, 52, change_dimbydim_prec);
 
        if (!finite)
        { // dehom point is infinite
          CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (!soln)
          { // classify as non-solution (junk)
            CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_Bertini_Junk;
          }
          else if (rankDef == 0)
          { // classify as non-singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = NON_SINGULAR;
          }
          else
          { // classify as singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = SINGULAR;
          }
        }
      }
      else
      { // path was not a success - copy over error code
        CD_copy[oid].codim[codim_index].endPt_types[i] = CD_copy[oid].codim[codim_index].endPts_d[i].retVal;
      }

      // setup sortPts
      sortPts[i].path_num = i;
      sortPts[i].norm = infNormVec_d(CD_copy[oid].codim[codim_index].endPts_d[i].endPt);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

    // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        j = i + 1;
        while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
        {
          indexJ = sortPts[j].path_num;
          if (CD->codim[codim_index].endPts_d[indexJ].retVal == 0)
          { // find difference if the jth path is successful
            vec_sub_d(tempVec, CD->codim[codim_index].endPts_d[indexI].endPt, CD->codim[codim_index].endPts_d[indexJ].endPt);
  
            if (infNormVec_d(tempVec) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
          j++;
        }
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // store the counts
    CD->codim[codim_index].num_superset = num_sing + num_nonsing;
    CD->codim[codim_index].num_sing = num_sing;
    CD->codim[codim_index].num_nonsing = num_nonsing;
    CD->codim[codim_index].num_nonsolns = num_nonsoln;
    CD->codim[codim_index].num_bad = num_bad;
    CD->codim[codim_index].num_inf = num_inf;

    // clear the memory
    clear_vec_d(tempVec);
    free(sortPts);
    for (i = max - 1; i >= 0; i--)
      clear_point_data_d(&PD_d[i]);
    free(PD_d);
  }
  else if (T[0].MPType == 1)
  { // sort using fixed multi precision
    point_data_mp *PD_mp = (point_data_mp *)bmalloc(max * sizeof(point_data_mp));
    sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));
    vec_mp tempVec;
    mpf_t norm_diff;

    // initialize tempVec, norm_diff & PD_mp
    init_vec_mp(tempVec, 0);
    mpf_init(norm_diff);
    for (i = 0; i < max; i++)
      init_point_data_mp(&PD_mp[i], 0);

    // determine if each path is (a) rank deficient (b) at infinity and (c) solution
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, soln, finite, oid) shared(num_paths) schedule(runtime) if (num_paths > 0)
#endif
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      rankDef = soln = finite = 0;

      // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // only check the ones that were successful
      if (CD_copy[oid].codim[codim_index].endPts_mp[i].retVal == 0)
      { // setup for evaluation
        CD_copy[oid].curr_codim_index = codim_index;

        // setup PD_mp[oid]
        point_cp_mp(PD_mp[oid].point, CD_copy[oid].codim[codim_index].endPts_mp[i].endPt);
        set_zero_mp(PD_mp[oid].time);

        // determine if it is rank deficient - i.e. the corank is positive
        rankDef = CD_copy[oid].codim[codim_index].endPts_mp[i].corank > 0;

        // do the classification
        dimbydimSortEndpoint(&rankDef, &finite, &soln, CD_copy[oid].codim[codim_index].endPts_mp[i].cond_num, &CD_copy[oid], codim_index, i, &T[oid], OUT[oid], NULL, &PD_mp[oid], T->Precision, NULL, CD_copy[oid].codim[codim_index].endPts_mp[i].last_approx, T->Precision, change_dimbydim_prec);

        if (!finite)
        { // dehom point is infinite
          CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (!soln)
          { // classify as non-solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_Bertini_Junk;
          }
          else if (rankDef == 0)
          { // classify as non-singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = NON_SINGULAR;
          }
          else
          { // classify as singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = SINGULAR;
          }
        }
      }
      else
      { // path was not a success - copy over error code
        CD_copy[oid].codim[codim_index].endPt_types[i] = CD_copy[oid].codim[codim_index].endPts_mp[i].retVal;
      }

      // setup sortPts
      mpf_init(sortPts[i].norm);
      sortPts[i].path_num = i;
      infNormVec_mp2(sortPts[i].norm, CD_copy[oid].codim[codim_index].endPts_mp[i].endPt);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

    // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        cont = 1;
        j = i;
        do
        { // increment the counter - start at i + 1
          j++;

          if (j < num_paths)
          { // check norm_diff
            indexJ = sortPts[j].path_num;
            mpf_sub(norm_diff, sortPts[j].norm, sortPts[i].norm);
            if (mpf_get_d(norm_diff) > final_tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we can continue
          if (cont && CD->codim[codim_index].endPts_mp[indexJ].retVal == 0)
          { // find difference if the jth path is successful
            vec_sub_mp(tempVec, CD->codim[codim_index].endPts_mp[indexI].endPt, CD->codim[codim_index].endPts_mp[indexJ].endPt);

            if (infNormVec_mp(tempVec) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
        } while (cont);
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // store the counts
    CD->codim[codim_index].num_superset = num_sing + num_nonsing;
    CD->codim[codim_index].num_sing = num_sing;
    CD->codim[codim_index].num_nonsing = num_nonsing;
    CD->codim[codim_index].num_nonsolns = num_nonsoln;
    CD->codim[codim_index].num_bad = num_bad;
    CD->codim[codim_index].num_inf = num_inf;

    // clear the memory
    for (i = num_paths - 1; i >= 0; i--)
      mpf_clear(sortPts[i].norm);
    free(sortPts);

    for (i = max - 1; i >= 0; i--)
      clear_point_data_mp(&PD_mp[i]);
    free(PD_mp);

    mpf_clear(norm_diff);
    clear_vec_mp(tempVec);
  }
  else
  { // sort using AMP
    point_data_d *PD_d = (point_data_d *)bmalloc(max * sizeof(point_data_d));
    point_data_mp *PD_mp = (point_data_mp *)bmalloc(max * sizeof(point_data_mp));
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
    mpf_t norm_diff_mp;

    // initialize MP
    mpf_init2(norm_diff_mp, 64);
    for (i = 0; i < max; i++)
    {
      init_point_data_d(&PD_d[i], 0);
      init_point_data_mp2(&PD_mp[i], 0, 64);
    }

    // determine if each path is (a) rank deficient (b) at infinity and (c) solution
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, soln, finite, oid) shared(num_paths) schedule(runtime) if (num_paths > 0)
#endif
    for (i = 0; i < num_paths; i++)
    { // get the current thread number
      oid = thread_num();

      rankDef = finite = soln = 0;

      // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Sorting %d of %d\n", i, num_paths);

      // only check the ones that were successful
      if (CD_copy[oid].codim[codim_index].endPts_amp[i].retVal == 0)
      { // setup for evaluation
        CD_copy[oid].curr_codim_index = codim_index;

        if (CD_copy[oid].codim[codim_index].endPts_amp[i].curr_prec < 64)
        { // setup PD_d[oid]
          point_cp_d(PD_d[oid].point, CD_copy[oid].codim[codim_index].endPts_amp[i].endPt_d);
          set_zero_d(PD_d[oid].time);
        }
        else
        { // setup PD_mp[oid]
          setprec_point_mp(PD_mp[oid].point, CD_copy[oid].codim[codim_index].endPts_amp[i].curr_prec);
          setprec_mp(PD_mp[oid].time, CD_copy[oid].codim[codim_index].endPts_amp[i].curr_prec);

          point_cp_mp(PD_mp[oid].point, CD_copy[oid].codim[codim_index].endPts_amp[i].endPt_mp);
          set_zero_mp(PD_mp[oid].time);
        }

        // determine if it is rank deficient - i.e. the corank is positive
        rankDef = CD_copy[oid].codim[codim_index].endPts_amp[i].corank > 0;

        // do the classification
        dimbydimSortEndpoint(&rankDef, &finite, &soln, CD_copy[oid].codim[codim_index].endPts_amp[i].cond_num, &CD_copy[oid], codim_index, i, &T[oid], OUT[oid], &PD_d[oid], &PD_mp[oid], CD_copy[oid].codim[codim_index].endPts_amp[i].curr_prec, CD_copy[oid].codim[codim_index].endPts_amp[i].last_approx_d, CD_copy[oid].codim[codim_index].endPts_amp[i].last_approx_mp, CD_copy[oid].codim[codim_index].endPts_amp[i].last_approx_prec, change_dimbydim_prec);

        if (!finite)
        { // dehom point is infinite
          CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_going_to_infinity;
        }
        else // dehom point is finite
        { // determine which category it belongs
          if (!soln)
          { // classify as non-solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = retVal_Bertini_Junk;
          }
          else if (rankDef == 0)
          { // classify as non-singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = NON_SINGULAR;
          }
          else
          { // classify as singular solution
            CD_copy[oid].codim[codim_index].endPt_types[i] = SINGULAR;
          }
        }
      }
      else
      { // path was not a success - copy over error code
        CD_copy[oid].codim[codim_index].endPt_types[i] = CD_copy[oid].codim[codim_index].endPts_amp[i].retVal;
      }

      // setup sortPts
      sortPts[i].path_num = i;
      if (CD_copy[oid].codim[codim_index].endPts_amp[i].curr_prec < 64)
        sortPts[i].norm = infNormVec_d(CD_copy[oid].codim[codim_index].endPts_amp[i].endPt_d);
      else
        sortPts[i].norm = infNormVec_mp(CD_copy[oid].codim[codim_index].endPts_amp[i].endPt_mp);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

    // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        cont = 1;
        j = i;
        do
        { // increment the counter - start at i + 1
          j++;

          if (j < num_paths)
          { // chech norm_diff
            indexJ = sortPts[j].path_num;
            if (sortPts[i].norm - sortPts[j].norm > final_tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we can continue
          if (cont && CD->codim[codim_index].endPts_amp[indexJ].retVal == 0)
          { // find the difference if the jth path is successful
            findDiff_endpoint_data(norm_diff_mp, CD->codim[codim_index].endPts_amp[indexI], CD->codim[codim_index].endPts_amp[indexJ]);

            if (mpf_get_d(norm_diff_mp) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
        } while (cont);
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // store the counts
    CD->codim[codim_index].num_superset = num_sing + num_nonsing;
    CD->codim[codim_index].num_sing = num_sing;
    CD->codim[codim_index].num_nonsing = num_nonsing;
    CD->codim[codim_index].num_nonsolns = num_nonsoln;
    CD->codim[codim_index].num_bad = num_bad;
    CD->codim[codim_index].num_inf = num_inf;

    // clear the memory
    free(sortPts);

    mpf_clear(norm_diff_mp);
    for (i = max - 1; i >= 0; i--)
    {
      clear_point_data_d(&PD_d[i]);
      clear_point_data_mp(&PD_mp[i]);
    }
    free(PD_mp);
    free(PD_d);
  }

  return;
}

int determineDimbyDimFinite_d(double max_norm, point_d point, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is finite, 0 - point is infinite     *
* NOTES: determines if point is finite or infinite              *
\***************************************************************/
{
  int retVal = 1;
  point_d orig_vars, dehom;
  init_point_d(orig_vars, 0);
  init_point_d(dehom, 0);

  dimbydimFindOrigVarsDehom_d(orig_vars, dehom, point, CD, codim_index);

  if (max_norm > infNormVec_d(dehom))
    retVal = 1;
  else
    retVal = 0;

  clear_point_d(orig_vars);
  clear_point_d(dehom);

  return retVal;
}

int determineDimbyDimFinite_mp(double max_norm, point_mp point, int prec, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is finite, 0 - point is infinite     *
* NOTES: determines if point is finite or infinite              *
\***************************************************************/
{
  int retVal = 1;
  point_mp orig_vars, dehom;

  init_point_mp2(orig_vars, 0, prec);
  init_point_mp2(dehom, 0, prec);

  dimbydimFindOrigVarsDehom_mp(orig_vars, dehom, point, CD, codim_index);

  if (max_norm > infNormVec_mp(dehom))
    retVal = 1;
  else
    retVal = 0;

  clear_point_mp(orig_vars);
  clear_point_mp(dehom);

  return retVal;
}

int determineDimbyDimSoln_d(double tol, double ratio, point_d point, point_d last_point, comp_d time, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is a solution, 0 - nonsolution       *
* NOTES: determines if point satisfies the linear slices        *
\***************************************************************/
{
  int retVal, codim = CD->codim[codim_index].codim;
  point_d orig_vars, last_approx, dehom;

  init_point_d(orig_vars, 0);
  init_point_d(last_approx, 0);
  init_point_d(dehom, 0);

  dimbydimFindOrigVarsDehom_d(last_approx, dehom, last_point, CD, codim_index);
  dimbydimFindOrigVarsDehom_d(orig_vars, dehom, point, CD, codim_index);

  // do the junk check
  retVal = nonsolutions_check_d(CD->num_funcs, codim, orig_vars, last_approx, time, tol, ratio, CD->Prog);

  clear_point_d(orig_vars);
  clear_point_d(last_approx);
  clear_point_d(dehom);

  return retVal;
}

int determineDimbyDimSoln_mp(double tol, double ratio, point_mp point, point_mp last_point, comp_mp time, int prec, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - point is a solution, 0 - nonsolution       *
* NOTES: determines if point satisfies the linear slices        *
\***************************************************************/
{
  int retVal, codim = CD->codim[codim_index].codim;
  point_mp orig_vars, last_approx, dehom;

  init_point_mp2(orig_vars, 0, prec);
  init_point_mp2(last_approx, 0, prec);
  init_point_mp2(dehom, 0, prec);

  dimbydimFindOrigVarsDehom_mp(last_approx, dehom, last_point, CD, codim_index); 
  dimbydimFindOrigVarsDehom_mp(orig_vars, dehom, point, CD, codim_index); 
 
  // do the junk check
  retVal = nonsolutions_check_mp(CD->num_funcs, codim, orig_vars, last_approx, time, tol, ratio, CD->Prog);

  clear_point_mp(orig_vars);
  clear_point_mp(last_approx);
  clear_point_mp(dehom);

  return retVal;
}

int printDimbyDimFooter_d(codim_t *CD, int codim_index, int path_num, point_data_d *endPoint, point_d orig_vars, point_d dehomP, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for codim_index & path_num           *
\***************************************************************/
{
  int i, retVal_out = 0, isNumber = 1, codim = CD->codim[codim_index].codim;

  if (retVal_in && d_abs_d(endPoint->time) < T->minTrackT)
  { // display a warning message
    if (T->screenOut)
      printf("WARNING: Path %d for codim %d reached the minimum value of T (%e < %e) with retVal %d.\n", path_num, codim, endPoint->time->r, T->minTrackT, retVal_in);
    fprintf(OUT, "WARNING: Path %d for codim %d reached the minimum value of T (%e < %e) with retVal %d.\n", path_num, codim, endPoint->time->r, T->minTrackT, retVal_in);

    // consider it a success
    retVal_out = 0;
  }
  else
  { // copy the return value and print the result
    retVal_out = retVal_in;
    printResultOfPath(OUT, retVal_in, T);
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
      if (isnan(endPoint->point->coord[i].r) || isnan(endPoint->point->coord[i].i) || isinf(endPoint->point->coord[i].r) || isinf(endPoint->point->coord[i].i))
        isNumber = 0;

    if (!isNumber)
      retVal_out = retVal_NAN;
  }

  if (retVal_out)
  { // print the path number, error message, time and point to FAIL
    printFailureMsg_d(FAIL, endPoint, dehomP, path_num, retVal_out, isNumber, 0, trackCount, T);
  }

  // print the path footer to OUT
  printPathFooterOut_d(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, dehomP, T, CD->Prog, 1, 0);

  return retVal_out;
}

int printDimbyDimFooter_mp(codim_t *CD, int codim_index, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp dehomP, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for path level_num & path_num        *
\***************************************************************/
{
  int i, retVal_out = 0, isNumber = 1, codim = CD->codim[codim_index].codim;

  if (retVal_in && d_abs_mp(endPoint->time) < T->minTrackT)
  { // display a warning message
    if (T->screenOut)
      printf("WARNING: Path %d for codim %d reached the minimum value of T (%e < %e) with retVal %d.\n", path_num, codim, mpf_get_d(endPoint->time->r), T->minTrackT, retVal_in);
    fprintf(OUT, "WARNING: Path %d for codim %d reached the minimum value of T (%e < %e) with retVal %d.\n", path_num, codim, mpf_get_d(endPoint->time->r), T->minTrackT, retVal_in);

    // initially consider it a success
    retVal_out = 0;
  }
  else
  { // copy the return value and print the result
    retVal_out = retVal_in;
    printResultOfPath(OUT, retVal_in, T);
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
    if (!(mpfr_number_p(endPoint->point->coord[i].r) && mpfr_number_p(endPoint->point->coord[i].i)))
        isNumber = 0;

    if (!isNumber)
      retVal_out = retVal_NAN;
  }

  if (retVal_out)
  { // print the path number, error message, time and point to FAIL
    printFailureMsg_mp(FAIL, endPoint, dehomP, path_num, retVal_out, isNumber, 0, trackCount, T);
  }

  // print the path footer to OUT
  printPathFooterOut_mp(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, first_increase, dehomP, T, CD->Prog, 1, 0);

  return retVal_out;
}

void dimbydim_classifyCodim(codim_t *CD, int codim_index, double final_tol, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: completes the classificiation of endpoints for codim   *
\***************************************************************/
{
  int i, j, cont, indexI = 0, indexJ = 0, num_paths = CD->codim[codim_index].num_paths;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_nonsoln = 0;

  if (MPType == 0)
  { // sort in double precision
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
    vec_d tempVec;

    init_vec_d(tempVec, 0);

    for (i = 0; i < num_paths; i++)
    { // setup sortPts
      sortPts[i].path_num = i;
      sortPts[i].norm = infNormVec_d(CD->codim[codim_index].endPts_d[i].endPt);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

    // do the final classificiation
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        j = i + 1;
        while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
        {
          indexJ = sortPts[j].path_num;
          if (CD->codim[codim_index].endPts_d[indexJ].retVal == 0)
          { // find difference if the jth path is successful
            vec_sub_d(tempVec, CD->codim[codim_index].endPts_d[indexI].endPt, CD->codim[codim_index].endPts_d[indexJ].endPt);

            if (infNormVec_d(tempVec) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
          j++;
        }
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // clear the memory
    clear_vec_d(tempVec);
    free(sortPts);
  }
  else if (MPType == 1)
  { // sort in multi precision
    sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));
    vec_mp tempVec;
    mpf_t norm_diff;

    // initialize tempVec & norm_diff 
    init_vec_mp(tempVec, 0);
    mpf_init(norm_diff);

    for (i = 0; i < num_paths; i++)
    { // setup sortPts
      mpf_init(sortPts[i].norm);
      sortPts[i].path_num = i;
      infNormVec_mp2(sortPts[i].norm, CD->codim[codim_index].endPts_mp[i].endPt);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

    // do the final classificiation
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        cont = 1;
        j = i;
        do
        { // increment the counter - start at i + 1
          j++;

          if (j < num_paths)
          { // check norm_diff
            indexJ = sortPts[j].path_num;
            mpf_sub(norm_diff, sortPts[j].norm, sortPts[i].norm);
            if (mpf_get_d(norm_diff) > final_tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we can continue
          if (cont && CD->codim[codim_index].endPts_mp[indexJ].retVal == 0)
          { // find difference if the jth path is successful
            vec_sub_mp(tempVec, CD->codim[codim_index].endPts_mp[indexI].endPt, CD->codim[codim_index].endPts_mp[indexJ].endPt);

            if (infNormVec_mp(tempVec) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
        } while (cont);
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // clear the memory
    for (i = num_paths - 1; i >= 0; i--)
    {
      mpf_clear(sortPts[i].norm);
    }
    free(sortPts);
    mpf_clear(norm_diff);
    clear_vec_mp(tempVec);
  }
  else
  { // sort using AMP
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));
    mpf_t norm_diff_mp;

    // initialize MP
    mpf_init2(norm_diff_mp, 64);

    for (i = 0; i < num_paths; i++)
    { // setup sortPts
      sortPts[i].path_num = i;
      if (CD->codim[codim_index].endPts_amp[i].curr_prec < 64)
        sortPts[i].norm = infNormVec_d(CD->codim[codim_index].endPts_amp[i].endPt_d);
      else
        sortPts[i].norm = infNormVec_mp(CD->codim[codim_index].endPts_amp[i].endPt_mp);
    }

    // sort the structure - use qsort to make comparisons efficient
    qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

    // do the final classificiation
    for (i = 0; i < num_paths; i++)
    {
      indexI = sortPts[i].path_num;
      if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
      { // compare against successful paths to see if it is equal to any other path
        cont = 1;
        j = i;
        do
        { // increment the counter - start at i + 1
          j++;

          if (j < num_paths)
          { // chech norm_diff
            indexJ = sortPts[j].path_num;
            if (sortPts[i].norm - sortPts[j].norm > final_tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we can continue
          if (cont && CD->codim[codim_index].endPts_amp[indexJ].retVal == 0)
          { // find the difference if the jth path is successful
            findDiff_endpoint_data(norm_diff_mp, CD->codim[codim_index].endPts_amp[indexI], CD->codim[codim_index].endPts_amp[indexJ]);

            if (mpf_get_d(norm_diff_mp) < final_tol)
            { // i & j are the same
              CD->codim[codim_index].endPt_types[indexI] = CD->codim[codim_index].endPt_types[indexJ] = SINGULAR;
            }
          }
        } while (cont);
      }

      // add to count
      if (CD->codim[codim_index].endPt_types[indexI] == retVal_going_to_infinity || CD->codim[codim_index].endPt_types[indexI] == retVal_security_max)
        num_inf++;
      else if (CD->codim[codim_index].endPt_types[indexI] == NON_SINGULAR)
        num_nonsing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == SINGULAR)
        num_sing++;
      else if (CD->codim[codim_index].endPt_types[indexI] == retVal_Bertini_Junk)
        num_nonsoln++;
      else
        num_bad++;
    }

    // clear the memory
    free(sortPts);
    mpf_clear(norm_diff_mp);
  }
 
  // store the counts
  CD->codim[codim_index].num_superset = num_sing + num_nonsing;
  CD->codim[codim_index].num_sing = num_sing;
  CD->codim[codim_index].num_nonsing = num_nonsing;
  CD->codim[codim_index].num_nonsolns = num_nonsoln;
  CD->codim[codim_index].num_bad = num_bad;
  CD->codim[codim_index].num_inf = num_inf;

  return;
}

int store_dimbydim_endPoint(endgame_data_t *endPt, int corank, double smallest, double largest, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, void const *CD_d, void const *CD_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: store endPt                                            *
\***************************************************************/
{
  codim_t *CD = (codim_t *)CD_d;
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
    dimbydimFindOrigVarsDehom_d(orig_vars, dehom_d, endPt->PD_d.point, CD, codim_index);

    // print the footer for the point
    CD->codim[codim_index].endPts_d[path_num].retVal = printDimbyDimFooter_d(CD, codim_index, path_num, &endPt->PD_d, orig_vars, dehom_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

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
    dimbydimFindOrigVarsDehom_mp(orig_vars, dehom_mp, endPt->PD_mp.point, CD, codim_index);

    // print the footer for the point
    CD->codim[codim_index].endPts_mp[path_num].retVal = printDimbyDimFooter_mp(CD, codim_index, path_num, &endPt->PD_mp, orig_vars, dehom_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

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
      dimbydimFindOrigVarsDehom_d(orig_vars, dehom_d, endPt->PD_d.point, CD, codim_index);

      // print the footer for the point
      CD->codim[codim_index].endPts_amp[path_num].retVal = printDimbyDimFooter_d(CD, codim_index, path_num, &endPt->PD_d, orig_vars, dehom_d, endPt->condition_number, endPt->function_residual_d, endPt->latest_newton_residual_d, endPt->t_val_at_latest_sample_point_d, endPt->error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

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
      dimbydimFindOrigVarsDehom_mp(orig_vars, dehom_mp, endPt->PD_mp.point, CD, codim_index);

      // print the footer for the point
      CD->codim[codim_index].endPts_amp[path_num].retVal = printDimbyDimFooter_mp(CD, codim_index, path_num, &endPt->PD_mp, orig_vars, dehom_mp, endPt->condition_number, endPt->first_increase, endPt->function_residual_mp, endPt->latest_newton_residual_mp, endPt->t_val_at_latest_sample_point_mp, endPt->error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, endPt->retVal, T, trackCount);

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

void findDiff_endpoint_data(mpf_t norm_diff, endpoint_data_amp endPt1, endpoint_data_amp endPt2)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the norm of the difference of endPt1 and endPt2  *
\***************************************************************/
{
  findDiff_point(norm_diff, endPt1.endPt_d, endPt1.endPt_mp, endPt1.curr_prec, endPt2.endPt_d, endPt2.endPt_mp, endPt2.curr_prec);

  return;
}

void dimbydimSortEndpoint(int *rankDef, int *finite, int *soln, double condNum, codim_t *CD, int codim_index, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoint - rankDef is already known          *
\***************************************************************/
{
  double residTol = 0;

  // initialize
  *finite = *soln = 1;

  if (endPt_prec >= 64 && T->MPType == 2)
  { // set the precision correctly
    initMP(endPt_prec);
    change_prec(CD, endPt_prec);
  }

  // check to see if finite
  if (endPt_prec < 64)
    *finite = determineDimbyDimFinite_d(T->finiteThreshold, endPt_d->point, CD, codim_index);
  else
    *finite = determineDimbyDimFinite_mp(T->finiteThreshold, endPt_mp->point, endPt_prec, CD, codim_index);

  if (*finite)
  { // check to see if it a solution
    if (endPt_prec < 64)
    {
      if (last_approx_prec > 52)
      { // convert to _d
        point_mp_to_d(last_approx_d, last_approx_mp);
      }
      residTol = MAX(T->funcResTol, 1e-15);
      *soln = determineDimbyDimSoln_d(residTol, T->ratioTol, endPt_d->point, last_approx_d, endPt_d->time, CD, codim_index);
    }
    else
    {
      if (last_approx_prec < 64)
      { // convert to _mp
        setprec_point_mp(last_approx_mp, endPt_prec);
        point_d_to_mp(last_approx_mp, last_approx_d);
      }
      residTol = pow(10, -prec_to_digits(endPt_prec) + 1);
      residTol = MAX(T->funcResTol, residTol);
      *soln = determineDimbyDimSoln_mp(residTol, T->ratioTol, endPt_mp->point, last_approx_mp, endPt_mp->time, endPt_prec, CD, codim_index);
    }
  }

  // we already know rankDef, so we can print the info
  fprintf(OUT, "Path number: %d Finite: %d Soln: %d Rank Def: %d CN: %e\n", path_num, *finite, *soln, *rankDef, condNum);

  return;
}

int dimbydim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  codim_t *CD = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    CD = (codim_t *)ED_d;
    vec_d orig_d;
    init_vec_d(orig_d, 0);

    dimbydimFindOrigVarsDehom_d(orig_d, out_d, in_d, CD, CD->curr_codim_index);

    clear_vec_d(orig_d);
  }
  else
  { // compute out_mp
    CD = (codim_t *)ED_mp;
    vec_mp orig_mp;
    init_vec_mp2(orig_mp, 0, *out_prec);
    setprec_vec_mp(out_mp, *out_prec);

    dimbydimFindOrigVarsDehom_mp(orig_mp, out_mp, in_mp, CD, CD->curr_codim_index);

    clear_vec_mp(orig_mp);
  }

  CD = NULL;

  return 0;
}


