// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "eqbyeq.h"
#include "parallel.h"

// provides the functions to track the system using the equation-by-equation method of Sommese, Verschelde and Wampler
void setup_eqbyeq_omp_d(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, FILE ***NONSOLN_copy, FILE *NONSOLN, char *nonName, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_d **BED_copy, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp);
void setup_eqbyeq_omp_mp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, FILE ***NONSOLN_copy, FILE *NONSOLN, char *nonName, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_mp **BED_copy, basic_eval_data_mp *ED);

void clear_eqbyeq_omp_d(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, char *nonName, tracker_config_t **T_copy, basic_eval_data_d **BED_copy);
void clear_eqbyeq_omp_mp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, char *nonName, tracker_config_t **T_copy, basic_eval_data_mp **BED_copy);

void eqbyeqWitnessTrack_d(int max_threads, int pathMod, basic_eval_data_d ED_d[], tracker_config_t T[], FILE **OUT, FILE **MIDOUT, FILE **RAWOUT, FILE **FAIL, FILE **NONSOLN, int subsystem_num, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void eqbyeqWitnessTrack_mp(int max_threads, int pathMod, basic_eval_data_mp ED[], tracker_config_t T[], FILE **OUT, FILE **MIDOUT, FILE **RAWOUT, FILE **FAIL, FILE **NONSOLN, int subsystem_num, trackingStats trackCount[], int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void eqbyeqWitnessTrackPath_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int subsystem_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void eqbyeqWitnessTrackPath_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int subsystem_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void sortWitnessEndpoints_d(int max_threads, int pathMod, point_data_d PD_d[], basic_eval_data_d ED_copy_d[], basic_eval_data_d *ED_d, tracker_config_t T[], FILE **OUT, int subsystem, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void sortWitnessEndpoints_mp(int max_threads, int pathMod, point_data_mp PD[], basic_eval_data_mp ED_copy[], basic_eval_data_mp *ED, tracker_config_t T[], FILE **OUT, int subsystem, double final_tol, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void eqbyeqStageTrack_d(int max_threads, int pathMod, basic_eval_data_d ED_d[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, FILE **NONSOLN, int stage_num, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void eqbyeqStageTrack_mp(int max_threads, int pathMod, basic_eval_data_mp ED[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, FILE **NONSOLN, int stage_num, trackingStats trackCount[], int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void eqbyeqStageTrackPath_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int stage_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void eqbyeqStageTrackPath_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int stage_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void sortStageEndpoints_d(int max_threads, int pathMod, point_data_d PD_d[], basic_eval_data_d ED_copy_d[], basic_eval_data_d *ED_d, tracker_config_t T[], FILE **OUT, int stage, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void sortStageEndpoints_mp(int max_threads, int pathMod, point_data_mp PD[], basic_eval_data_mp ED_copy[], basic_eval_data_mp *ED, tracker_config_t T[], FILE **OUT, int stage, double final_tol, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int determineEqbyEqFinite(double maxNorm, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int num, int isStage);

void eqbyeq_track_d(FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking via eq-by-eq of SVW     *
\***************************************************************/
{
  int i, num_paths, num_vars, num_crossings, max = max_threads();
  point_data_d *tempPoint = (point_data_d *)bmalloc(max * sizeof(point_data_d));
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  char outName[] = "output", midName[] = "midout", rawName[] = "rawout", failName[] = "fail", nonName[] = "nonsolutions";
  FILE *MIDOUT = fopen(midFile, "w"), *NONSOLN = NULL;
  if (MIDOUT == NULL)
  {
    printf("ERROR: '%s' is not a valid name for a file!\n", midFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  if (!ED_d->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen(nonName, "w");
    fprintf(NONSOLN, "                                    \n\n");
  }  

  // pointers for OpenMP tracking
  tracker_config_t *T_copy = NULL;
  basic_eval_data_d *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, **NONSOLN_copy = NULL;

  // setup the structures
  setup_eqbyeq_omp_d(max, &trackCount_copy, trackCount, &OUT_copy, OUT, outName, &RAWOUT_copy, RAWOUT, rawName, &MIDOUT_copy, MIDOUT, midName, &FAIL_copy, FAIL, failName, &NONSOLN_copy, NONSOLN, nonName, &T_copy, T, &BED_copy, ED_d, ED_mp);
  for (i = 0; i < max; i++)
  {
    init_point_data_d(&tempPoint[i], 0);
  }

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // generate the witness sets for each of the subsystems
  ptr_to_eval_d = &witness_eqbyeq_eval_d;
  ptr_to_eval_mp = &witness_eqbyeq_eval_mp;
  change_prec = &change_eqbyeq_eval_prec;

  for (i = 0; i < ED_d->EqD->num_subsystems; i++)
  { // find the number of paths for this subsystem
    num_paths = ED_d->EqD->witnessData_d[i].num_paths;
    num_vars = ED_d->EqD->witnessData_d[i].depth;

    // track each path for this subsystem
    eqbyeqWitnessTrack_d(max, pathMod, BED_copy, T_copy, OUT_copy, MIDOUT_copy, RAWOUT_copy, FAIL_copy, NONSOLN_copy, i, trackCount_copy, ptr_to_eval_d, ptr_to_eval_mp, change_prec, eqbyeq_witness_dehom);

    // combine the files
    combine_omp_file(OUT, &OUT_copy, max);
    combine_omp_file(MIDOUT, &MIDOUT_copy, max);
    combine_omp_file(RAWOUT, &RAWOUT_copy, max);
    combine_omp_file(FAIL, &FAIL_copy, max);

    // make sure that no path crossing happened when finding this subsystem, if this is not the only subsystem
    if (ED_d->EqD->num_subsystems > 1)
    { // close MIDOUT and check for path crossings
      fclose(MIDOUT);
      num_crossings = 0;
      midpoint_checker(num_paths, num_vars, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // setup OUT_copy for sorting
    setup_omp_file(&OUT_copy, OUT, outName, max);

    // sort the endpoints
    sortWitnessEndpoints_d(max, pathMod, tempPoint, BED_copy, ED_d, T_copy, OUT_copy, i, target_tol, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

    // combine OUT
    combine_omp_file(OUT, &OUT_copy, max);

    // if this is not the last subsystem, reopen the files
    if (i + 1 < ED_d->EqD->num_subsystems)
    { // reopen the files
      MIDOUT = fopen(midFile, "w");
      setup_omp_file(&MIDOUT_copy, MIDOUT, midName, max);
      setup_omp_file(&OUT_copy, OUT, outName, max);
      setup_omp_file(&RAWOUT_copy, RAWOUT, rawName, max);
      setup_omp_file(&FAIL_copy, FAIL, failName, max);
    }
  }
  // setup the first stage
  setupEqbyEqFirstStage_d(ED_d->EqD, T->MPType);
  setup_omp_eqbyeq_first_stage_d(max, BED_copy, ED_d->EqD, T->MPType);
  // clear the first witness data
  clearEqbyEqFirstWitnessData_d(ED_d->EqD, T->MPType);
  clear_omp_eqbyeq_witness_d(max, BED_copy, 0, T->MPType);

  for (i = 1; i < ED_d->EqD->num_subsystems; i++)
  { // setup stage i
    setupEqbyEqNextStage_d(ED_d, i, T->MPType, T->AMP_max_prec);
    setup_omp_eqbyeq_next_stage_d(max, BED_copy, ED_d->EqD, i, T->MPType);
    // clear the stage date from stage i - 1 since it is no longer needed
    clearEqbyEqStageData_d(ED_d->EqD, i - 1, T->MPType);
    clear_omp_eqbyeq_stage_d(max, BED_copy, i - 1, T->MPType);
    // clear the witness data from subsystem i since it is no longer needed
    clearEqbyEqWitnessData_d(ED_d->EqD, i, T->MPType);
    clear_omp_eqbyeq_witness_d(max, BED_copy, i, T->MPType);

    // setup the files to track stage i
    MIDOUT = fopen(midFile, "w");
    setup_omp_file(&MIDOUT_copy, MIDOUT, midName, max);
    setup_omp_file(&OUT_copy, OUT, outName, max);
    setup_omp_file(&RAWOUT_copy, RAWOUT, rawName, max);
    setup_omp_file(&FAIL_copy, FAIL, failName, max);

    // find the number of paths for this next stage
    num_paths = ED_d->EqD->stageData_d[i].num_paths;
    if (ED_d->EqD->stageData_d[i].useIntrinsicSlice)
      num_vars = ED_d->EqD->stageData_d[i].depth_x + ED_d->EqD->stageData_d[i].depth_y;
    else
      num_vars = 2 * ED_d->EqD->num_vars;

    // setup the evaluators for stage tracking
    ptr_to_eval_d = &standard_eqbyeq_eval_d;
    ptr_to_eval_mp = &standard_eqbyeq_eval_mp;

    // track each path for this stage
    eqbyeqStageTrack_d(max, pathMod, BED_copy, T_copy, OUT_copy, RAWOUT_copy, MIDOUT_copy, FAIL_copy, NONSOLN_copy, i, trackCount_copy, ptr_to_eval_d, ptr_to_eval_mp, change_prec, eqbyeq_stage_dehom);

    // combine all of the files
    combine_omp_file(OUT, &OUT_copy, max);
    combine_omp_file(MIDOUT, &MIDOUT_copy, max);
    combine_omp_file(RAWOUT, &RAWOUT_copy, max);
    combine_omp_file(FAIL, &FAIL_copy, max);
 
    // make sure that no path crossing happened when finding this stage, if this is not the last stage
    if (i + 1 < ED_d->EqD->num_subsystems)
    { // close MIDOUT and check for path crossings
      fclose(MIDOUT);
      num_crossings = 0;
      midpoint_checker(num_paths, num_vars, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }
 
    // setup OUT_copy for sorting
    setup_omp_file(&OUT_copy, OUT, outName, max);

    // setup the evaluators for stage sorting
    ptr_to_eval_d = &stage_sort_eqbyeq_eval_d;
    ptr_to_eval_mp = &stage_sort_eqbyeq_eval_mp;

    // sort the endpoints
    sortStageEndpoints_d(max, pathMod, tempPoint, BED_copy, ED_d, T_copy, OUT_copy, i, target_tol, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

    // combine OUT
    combine_omp_file(OUT, &OUT_copy, max);
  }
  // close MIDOUT
  fclose(MIDOUT);

  if (!ED_d->squareSystem.noChanges)
  { // complete NONSOLN
    combine_omp_file(NONSOLN, &NONSOLN_copy, max);
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }  

  // combine data and clear copies - just need to delete the extra files
  clear_eqbyeq_omp_d(max, &trackCount_copy, trackCount, outName, rawName, midName, failName, nonName, &T_copy, &BED_copy);

  // set the total number tracked on last stage to trackCount
  trackCount->numPoints = ED_d->EqD->stageData_d[ED_d->EqD->num_subsystems - 1].num_paths;

  // free tempPoint
  for (i = 0; i < max; i++)
  {
    clear_point_data_d(&tempPoint[i]);
  }
  free(tempPoint);

  return;
}

void setup_eqbyeq_omp_d(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, FILE ***NONSOLN_copy, FILE *NONSOLN, char *nonName, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_d **BED_copy, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup everything needed to do eq-by-eq tracking        *
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
  setup_omp_file(NONSOLN_copy, NONSOLN, nonName, max_threads);

  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED_d;
    (*BED_copy)->BED_mp = ED_mp; // make sure that this is pointed to inside of ED_d
  }
  else // max_threads > 1
  { // allocate memory
    *trackCount_copy = (trackingStats *)bmalloc(max_threads * sizeof(trackingStats));
    *T_copy = (tracker_config_t *)bmalloc(max_threads * sizeof(tracker_config_t));
    *BED_copy = (basic_eval_data_d *)bmalloc(max_threads * sizeof(basic_eval_data_d));

    // copy T, ED_d, ED_mp, & trackCount
    for (i = 0; i < max_threads; i++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[i], T);
      // copy ED_d & ED_mp
      cp_basic_eval_data_d(&(*BED_copy)[i], ED_d, ED_mp, T->MPType);
      // setup EqD
      setup_omp_eqbyeq_d(&(*BED_copy)[i], ED_d->EqD, T->MPType);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[i]);
      (*trackCount_copy)[i].numPoints = trackCount->numPoints;
    }
  }

  return;
}

void clear_eqbyeq_omp_d(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, char *nonName, tracker_config_t **T_copy, basic_eval_data_d **BED_copy)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear everything needed to do eq-by-eq tracking        *
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
    *BED_copy = NULL;
  }
  else if (max_threads > 1)
  { // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);

    // clear the copies of T & ED_d
    for (i = max_threads - 1; i >= 0; i--)
    { // clear BED_copy - 0 since we want to target exactly want needs cleared in the eq-by-eq data
      basic_eval_clear_d(&(*BED_copy)[i], 0, (*T_copy)[i].MPType);
      // clear EqD
      clear_omp_eqData_d(&(*BED_copy)[i], (*T_copy)[i].MPType);
      // clear T_copy
      tracker_config_clear(&(*T_copy)[i]);
    }

    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);

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

      // delete nonsolutions
      size = 1 + snprintf(NULL, 0, "%s_%d", nonName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", nonName, i);
      remove(str);
    }
    free(str);
  }

  return;
}

void eqbyeqWitnessTrack_d(int max_threads, int pathMod, basic_eval_data_d ED_d[], tracker_config_t T[], FILE **OUT, FILE **MIDOUT, FILE **RAWOUT, FILE **FAIL, FILE **NONSOLN, int subsystem_num, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths to find witness set for this subsystem*
* Uses OpenMP if available                                      *
\***************************************************************/
{
  int i, oid, num_paths = ED_d[0].EqD->witnessData_d[subsystem_num].num_paths, num_subs = ED_d[0].EqD->num_subsystems;

  // display messages
  printf("\nFinding witness points for subsystem %d of %d: %d path%s to track.\n", subsystem_num, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Finding witness points for subsystem %d.\n", subsystem_num);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this subsystem
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid) schedule(runtime) if (num_paths > 2 * max_threads)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // track the path
    eqbyeqWitnessTrackPath_d(&ED_d[oid], ED_d[oid].BED_mp, &T[oid], OUT[oid], MIDOUT[oid], RAWOUT[oid], FAIL[oid], NONSOLN[oid], subsystem_num, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  return;
}

void eqbyeqWitnessTrackPath_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int subsystem_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in subsystem_num & path_num of ED_d    *
\***************************************************************/
{
  int num_subs = ED_d->EqD->num_subsystems;

  endgame_data_t endPt;
  point_data_d startPt;
  point_d dehom_d, orig_last_d;
  point_mp dehom_mp, orig_last_mp, tempPoint_mp;

  // initialize
  init_endgame_data(&endPt, T->Precision);
  init_point_data_d(&startPt, 0);
  init_point_d(dehom_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last_mp, 0); init_point_mp(tempPoint_mp, 0);

  // seutp for tracking the path
  point_cp_d(startPt.point, ED_d->EqD->witnessData_d[subsystem_num].startPts[path_num]);
  set_double_d(startPt.time, 1, 0);
  ED_d->EqD->curr_stage_num = subsystem_num;
  ED_d->EqD->curr_path_num = path_num;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 2)
  { // initialize for AMP tracking
    ED_mp->EqD->curr_stage_num = subsystem_num;
    ED_mp->EqD->curr_path_num = path_num;
  }

  // print the header for the point
  printPathHeader_d(OUT, &startPt, T, path_num, ED_d, ptr_to_eval_d);

  zero_dim_track_path_d(path_num, &endPt, &startPt, OUT, MIDOUT, T, ED_d, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

  // check to see if it should be sharpened - only when this is the only subsystem
  if (num_subs == 1 && endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, ED_d, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // store the condition number
  ED_d->EqD->witnessData_d[subsystem_num].condition_nums[path_num] = endPt.condition_number;

  // store if higher dimenaional
  ED_d->EqD->witnessData_d[subsystem_num].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, ED_d, ED_mp, subsystem_num, 0);

  // setup orig_last
  if (endPt.last_approx_prec < 64)
  { // use _d
    intrinsicToExtrinsic_d(orig_last_d, endPt.last_approx_d, ED_d->EqD->witnessData_d[subsystem_num].B, ED_d->EqD->witnessData_d[subsystem_num].p);

    if (endPt.prec > 52)
    { // convert to _mp
      setprec_point_mp(orig_last_mp, endPt.prec);
      point_d_to_mp(orig_last_mp, orig_last_d);
    }
  }
  else
  { // use _mp
    setprec_point_mp(orig_last_mp, endPt.last_approx_prec);

    intrinsicToExtrinsic_mp(orig_last_mp, endPt.last_approx_mp, ED_mp->EqD->witnessData_mp[subsystem_num].B, ED_mp->EqD->witnessData_mp[subsystem_num].p);

    if (endPt.prec < 64)
    { // convert _d
      point_mp_to_d(orig_last_d, orig_last_mp);
    }
  }

  if (endPt.prec < 64)
  { // copy the data to the appropriate spot
    point_cp_d(ED_d->EqD->witnessData_d[subsystem_num].endPts_in[path_num], endPt.PD_d.point);
    set_d(ED_d->EqD->witnessData_d[subsystem_num].finalTs[path_num], endPt.PD_d.time);
  
    // convert to the original coordinates
    intrinsicToExtrinsic_d(ED_d->EqD->witnessData_d[subsystem_num].endPts[path_num], endPt.PD_d.point, ED_d->EqD->witnessData_d[subsystem_num].B, ED_d->EqD->witnessData_d[subsystem_num].p);

    // find dehom_d
    getDehomPoint_d(dehom_d, ED_d->EqD->witnessData_d[subsystem_num].endPts[path_num], ED_d->EqD->witnessData_d[subsystem_num].endPts[path_num]->size, &ED_d->preProcData);

    // print the footer to OUT for the point
    ED_d->EqD->witnessData_d[subsystem_num].endPt_retVals[path_num] = printWitnessFooter_d(ED_d, subsystem_num, path_num, &endPt.PD_d, ED_d->EqD->witnessData_d[subsystem_num].endPts[path_num], orig_last_d, dehom_d, endPt.condition_number, endPt.function_residual_d, endPt.latest_newton_residual_d, endPt.t_val_at_latest_sample_point_d, endPt.error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);
  }
  else
  { // copy over to the appropriate spot - converting to double precision     
    setprec_point_mp(dehom_mp, endPt.prec);
    setprec_point_mp(tempPoint_mp, endPt.prec);

    point_mp_to_d(ED_d->EqD->witnessData_d[subsystem_num].endPts_in[path_num], endPt.PD_mp.point);
    mp_to_d(ED_d->EqD->witnessData_d[subsystem_num].finalTs[path_num], endPt.PD_mp.time);

    // convert to the original coordinates
    intrinsicToExtrinsic_mp(tempPoint_mp, endPt.PD_mp.point, ED_mp->EqD->witnessData_mp[subsystem_num].B, ED_mp->EqD->witnessData_mp[subsystem_num].p);
    intrinsicToExtrinsic_mp(tempPoint_mp, endPt.PD_mp.point, ED_mp->EqD->witnessData_mp[subsystem_num].B, ED_mp->EqD->witnessData_mp[subsystem_num].p);

    // copy tempPoint_mp to endPts
    point_mp_to_d(ED_d->EqD->witnessData_d[subsystem_num].endPts[path_num], tempPoint_mp);

    // find dehom_mp
    getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

    // print the footer to OUT for the point & find the condition number
    ED_d->EqD->witnessData_d[subsystem_num].endPt_retVals[path_num] = printWitnessFooter_mp(ED_mp, subsystem_num, path_num, &endPt.PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);
  }

  // clear endPt
  clear_endgame_data(&endPt);
  clear_point_data_d(&startPt);
  clear_point_d(dehom_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last_mp); clear_point_mp(tempPoint_mp);

  return;
}

void sortWitnessEndpoints_d(int max_threads, int pathMod, point_data_d PD_d[], basic_eval_data_d ED_copy_d[], basic_eval_data_d *ED_d, tracker_config_t T[], FILE **OUT, int subsystem, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points found for current subsystem   *
\***************************************************************/
// Still doing everything as intrinsic linears
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, oid, rankDef, finite, indexI = 0, indexJ = 0, num_paths = ED_d->EqD->witnessData_d[subsystem].num_paths, num_subs = ED_d->EqD->num_subsystems;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  vec_d tempVec;
  sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));

  init_vec_d(tempVec, 0);

  // print header for level
  printf("\nSorting witness points for subsystem %d of %d: %d path%s to sort.\n", subsystem, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Sorting witness points for subsystem %d.\n", subsystem);
  fprintf(OUT[0], "*****************************************************\n");

  // determine if each of the paths is rank deficient
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, oid, rankDef, finite) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sorting %d of %d\n", i, num_paths);

    // only check the ones that were successful
    if (ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_retVals[i] == 0)
    { // setup for evaluation
      ED_copy_d[oid].EqD->curr_stage_num = subsystem;
      ED_copy_d[oid].EqD->curr_path_num = i;

      if (T[oid].MPType == 2)
      {
        ED_copy_d[oid].BED_mp->EqD->curr_stage_num = subsystem;
        ED_copy_d[oid].BED_mp->EqD->curr_path_num = i;
      }

      // setup PD_d[oid]
      point_cp_d(PD_d[oid].point, ED_copy_d[oid].EqD->witnessData_d[subsystem].endPts_in[i]);
      set_zero_d(PD_d[oid].time);

      fprintf(OUT[oid], "Path number: %d\n", i);
      // determine if it is rank deficient, finite 
      eqbyeqWitnessSortEndpoint_d(&rankDef, &finite, ED_copy_d[oid].EqD->witnessData_d[subsystem].higherDim[i], &ED_copy_d[oid], ED_copy_d[oid].BED_mp, subsystem, i, &ED_copy_d[oid].EqD->witnessData_d[subsystem].condition_nums[i], &T[oid], OUT[oid], &PD_d[oid], NULL, 52, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

      // copy back to endPts_in - determineRankDef will update PD_d[oid] if it is improved
      point_cp_d(ED_copy_d[oid].EqD->witnessData_d[subsystem].endPts_in[i], PD_d[oid].point);

      // check for success
      if (T[oid].regen_remove_inf && !finite)
      { // dehom point is infinite
        ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_types[i] = retVal_going_to_infinity;
      } 
      else if (T[oid].regen_higher_dim_check && ED_copy_d[oid].EqD->witnessData_d[subsystem].higherDim[i])
      { // it lies on a higher dimensional component
        ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_types[i] = retVal_higher_dim;
      }
      else if (!rankDef)
      { // classify as non-singular
        ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_types[i] = MOVE_TO_NEXT;
      }
      else
      { // classify as singular
        ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_types[i] = DO_NOT_MOVE_TO_NEXT;
      }
    }
    else
    { // path was not a success - copy over error code
      ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_types[i] = ED_copy_d[oid].EqD->witnessData_d[subsystem].endPt_retVals[i];
    }

    // setup sortPts
    sortPts[i].path_num = i;
    sortPts[i].norm = infNormVec_d(ED_copy_d[oid].EqD->witnessData_d[subsystem].endPts_in[i]);
  }

  // sort the structure - use qsort to make comparisons efficient
  qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

  // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      j = i + 1;
      while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
      {
        indexJ = sortPts[j].path_num;
        if (ED_d->EqD->witnessData_d[subsystem].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          vec_sub_d(tempVec, ED_d->EqD->witnessData_d[subsystem].endPts_in[indexI], ED_d->EqD->witnessData_d[subsystem].endPts_in[indexJ]);
          if (infNormVec_d(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] = ED_d->EqD->witnessData_d[subsystem].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
        j++;
      }
    }

    // add to count
    if (ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == retVal_going_to_infinity || ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (ED_d->EqD->witnessData_d[subsystem].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  ED_d->EqD->witnessData_d[subsystem].num_sing = num_sing;
  ED_d->EqD->witnessData_d[subsystem].num_nonsing = num_nonsing;
  ED_d->EqD->witnessData_d[subsystem].num_bad = num_bad;
  ED_d->EqD->witnessData_d[subsystem].num_higher_dim = num_higher_dim;
  ED_d->EqD->witnessData_d[subsystem].num_inf = num_inf;

  // clear the memory
  free(sortPts);
  clear_point_d(tempVec);

  return;
}

void eqbyeqWitnessSortEndpoint_d(int *rankDef, int *finite, int higherDim, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int subsystem, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoint - rankDef, finite & higherDim       *
\***************************************************************/
{ // initialize
  *rankDef = 0; // we already know higherDim!
  *finite = 1;

  // check to see if finite
  if (T->regen_remove_inf)
  { // determine if the endpoint is finite
    *finite = determineEqbyEqFinite(T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, ED, ED_mp, subsystem, 0);
  }

  // determine if we need to do the rank deficient test (which also refines the endpoint)
  if ((T->regen_higher_dim_check && higherDim == 0) && *finite == 1)
  { // determine if it is rank deficient 
    *rankDef = determineRankDef(condNum, T->final_tolerance, NULL, 52, endPt_d, endPt_mp, endPt_prec, T, OUT, ED, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

    // check to see if finite
    if (T->regen_remove_inf)
    { // determine if the endpoint is finite
      *finite = determineEqbyEqFinite(T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, ED, ED_mp, subsystem, 0);
    }

    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d Rank Def: %d", higherDim, *rankDef);
    }
    else
    {
      fprintf(OUT, "Rank Def: %d", *rankDef);
    }
    if (T->regen_remove_inf)
      fprintf(OUT, " Finite: %d", *finite);
    fprintf(OUT, " CN: %e\n", *condNum);
  }
  else
  { // we know that we are removing this end point
    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d", higherDim);
      if (T->regen_remove_inf)
        fprintf(OUT, " Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
    else
    {
      if (T->regen_remove_inf)
        fprintf(OUT, "Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
  }

  return;
}

void eqbyeqStageTrack_d(int max_threads, int pathMod, basic_eval_data_d ED_d[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, FILE **NONSOLN, int stage_num, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths to find the witness points for stage  *
* Uses OpenMP if available                                      *
\***************************************************************/
{
  int i, oid, num_paths = ED_d[0].EqD->stageData_d[stage_num].num_paths, num_stages = ED_d[0].EqD->num_subsystems;

  // display messages
  printf("\nTracking points for stage %d of %d: %d path%s to track.\n", stage_num, num_stages, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Tracking points for stage %d.\n", stage_num);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this stage
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // track the path
    eqbyeqStageTrackPath_d(&ED_d[oid], ED_d[oid].BED_mp, &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], NONSOLN[oid], stage_num, i, &trackCount[oid], ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);
  }

  return;
}

void eqbyeqStageTrackPath_d(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int stage_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in stage_num & path_num of ED_d        *
\***************************************************************/
{
  int num_subs = ED_d->EqD->num_subsystems;

  endgame_data_t endPt;
  point_data_d startPt;
  point_d dehom_d, orig_last_d;
  point_mp dehom_mp, orig_last_mp, tempPoint_mp;

  // initialize
  init_endgame_data(&endPt, T->Precision);
  init_point_data_d(&startPt, 0);
  init_point_d(dehom_d, 0); init_point_d(orig_last_d, 0);
  init_point_mp(dehom_mp, 0); init_point_mp(orig_last_mp, 0); init_point_mp(tempPoint_mp, 0);

  // seutp for tracking the path
  point_cp_d(startPt.point, ED_d->EqD->stageData_d[stage_num].startPts[path_num]);

  set_double_d(startPt.time, 1, 0);
  ED_d->EqD->curr_stage_num = stage_num;
  ED_d->EqD->curr_path_num = path_num;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 2)
  { // initialize for AMP tracking
    init_endgame_data(&endPt, 64);

    ED_mp->EqD->curr_stage_num = stage_num;
    ED_mp->EqD->curr_path_num = path_num;
  }

  // print the header for the point
  printPathHeader_d(OUT, &startPt, T, path_num, ED_d, ptr_to_eval_d);

  // track the path
  zero_dim_track_path_d(path_num, &endPt, &startPt, OUT, MIDOUT, T, ED_d, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

  // check to see if it should be sharpened - only when this is the last stage
  if (stage_num + 1 == num_subs && endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, ED_d, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
  }

  // store the condition number
  ED_d->EqD->stageData_d[stage_num].condition_nums[path_num] = endPt.condition_number;

  ED_d->EqD->stageData_d[stage_num].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, ED_d, ED_mp, stage_num, 1);

  // find orig_last
  if (endPt.last_approx_prec < 64)
  { 
    if (ED_d->EqD->stageData_d[stage_num].useIntrinsicSlice)
    { // convert to extrinsic coordinates
      intrinsicToExtrinsic_d(orig_last_d, endPt.last_approx_d, ED_d->EqD->stageData_d[stage_num].B0, ED_d->EqD->stageData_d[stage_num].p0);
    }
    else
    { // copy over
      point_cp_d(orig_last_d, endPt.last_approx_d);
    }
    orig_last_d->size /= 2; // remove the bottom half of the coordinates

    if (endPt.prec > 52)
    { // convert to _mp
      setprec_point_mp(orig_last_mp, endPt.prec);
      point_d_to_mp(orig_last_mp, orig_last_d);
    }
  }
  else
  { 
    setprec_point_mp(orig_last_mp, endPt.last_approx_prec);
    if (ED_d->EqD->stageData_d[stage_num].useIntrinsicSlice)
    { // convert to extrinsic coordinates
      intrinsicToExtrinsic_mp(orig_last_mp, endPt.last_approx_mp, ED_d->EqD->stageData_mp[stage_num].B0, ED_d->EqD->stageData_mp[stage_num].p0);
    }
    else
    { // copy over
      point_cp_mp(orig_last_mp, endPt.last_approx_mp);
    }
    orig_last_d->size /= 2; // remove the bottom half of the coordinates

    if (endPt.prec < 64)
    { // convert to _d
      point_mp_to_d(orig_last_d, orig_last_mp);
    }
  }

  if (endPt.prec < 64)
  { // copy over to the appropriate spot
    set_d(ED_d->EqD->stageData_d[stage_num].finalTs[path_num], endPt.PD_d.time);

    if (ED_d->EqD->stageData_d[stage_num].useIntrinsicSlice)
    { // store the intrinsic endpoint and its corresponding point in the original variables
      intrinsicToExtrinsic_d(dehom_d, endPt.PD_d.point, ED_d->EqD->stageData_d[stage_num].B0, ED_d->EqD->stageData_d[stage_num].p0);
      dehom_d->size /= 2; // remove the bottom half of the coordinates

      point_cp_d(ED_d->EqD->stageData_d[stage_num].endPts[path_num], dehom_d);

      // convert to intrinsic coordinates
      mat_d B_transpose;
      init_mat_d(B_transpose, 0, 0);

      transpose_d(B_transpose, ED_d->EqD->stageData_d[stage_num].B);
      extrinsicToIntrinsic_d(ED_d->EqD->stageData_d[stage_num].endPts_in[path_num], dehom_d, B_transpose, ED_d->EqD->stageData_d[stage_num].p);

      clear_mat_d(B_transpose);
    }
    else
    { // use the top set of coordinates so that we can store the original coordinates
      point_cp_d(ED_d->EqD->stageData_d[stage_num].endPts[path_num], endPt.PD_d.point);
      ED_d->EqD->stageData_d[stage_num].endPts[path_num]->size /= 2;
    }

    // find dehom_d
    getDehomPoint_d(dehom_d, ED_d->EqD->stageData_d[stage_num].endPts[path_num], ED_d->EqD->stageData_d[stage_num].endPts[path_num]->size, &ED_d->preProcData);

    // print the footer for the point
    ED_d->EqD->stageData_d[stage_num].endPt_retVals[path_num] = printStageFooter_d(ED_d, stage_num, path_num, &endPt.PD_d, ED_d->EqD->stageData_d[stage_num].endPts[path_num], orig_last_d, dehom_d, endPt.condition_number, endPt.function_residual_d, endPt.latest_newton_residual_d, endPt.t_val_at_latest_sample_point_d, endPt.error_at_latest_sample_point_d, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);
  }
  else
  { // set precision dehom_mp & tempPoint_mp
    setprec_point_mp(dehom_mp, endPt.prec);
    setprec_point_mp(tempPoint_mp, endPt.prec);

    // copy over to the appropriate spot - converting to double precision
    mp_to_d(ED_d->EqD->stageData_d[stage_num].finalTs[path_num], endPt.PD_mp.time);

    if (ED_d->EqD->stageData_d[stage_num].useIntrinsicSlice)
    { // store the intrinsic endpoint and its corresponding point in the original variables
      intrinsicToExtrinsic_mp(tempPoint_mp, endPt.PD_mp.point, ED_d->EqD->stageData_mp[stage_num].B0, ED_d->EqD->stageData_mp[stage_num].p0);
      tempPoint_mp->size /= 2; // remove the bottom half of the coordinates

      point_mp_to_d(ED_d->EqD->stageData_d[stage_num].endPts[path_num], tempPoint_mp);

      // convert to the original coordinates
      mat_mp B_transpose;
      init_mat_mp2(B_transpose, 0, 0, endPt.prec);
      transpose_mp(B_transpose, ED_d->EqD->stageData_mp[stage_num].B);

      extrinsicToIntrinsic_mp(dehom_mp, tempPoint_mp, B_transpose, ED_d->EqD->stageData_mp[stage_num].p);

      // copy to endPts_in
      point_mp_to_d(ED_d->EqD->stageData_d[stage_num].endPts_in[path_num], dehom_mp);

      clear_mat_mp(B_transpose);
    }
    else
    { // use the top set of coordinates so that we can store the original coordinates
      point_cp_mp(tempPoint_mp, endPt.PD_mp.point);
      tempPoint_mp->size /= 2;

      point_mp_to_d(ED_d->EqD->stageData_d[stage_num].endPts[path_num], tempPoint_mp);
    }

    // find dehom_mp
    getDehomPoint_mp(dehom_mp, tempPoint_mp, tempPoint_mp->size, &ED_mp->preProcData);

    // print the footer for the point
    ED_d->EqD->stageData_d[stage_num].endPt_retVals[path_num] = printStageFooter_mp(ED_mp, stage_num, path_num, &endPt.PD_mp, tempPoint_mp, orig_last_mp, dehom_mp, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);
  }

  // clear
  clear_endgame_data(&endPt);
  clear_point_data_d(&startPt);
  clear_point_d(dehom_d); clear_point_d(orig_last_d);
  clear_point_mp(dehom_mp); clear_point_mp(orig_last_mp); clear_point_mp(tempPoint_mp);

  return;
}

void sortStageEndpoints_d(int max_threads, int pathMod, point_data_d PD_d[], basic_eval_data_d ED_copy_d[], basic_eval_data_d *ED_d, tracker_config_t T[], FILE **OUT, int stage, double final_tol, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points found for current stage       *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, oid, rankDef, finite, indexI = 0, indexJ = 0, num_paths = ED_d->EqD->stageData_d[stage].num_paths, num_subs = ED_d->EqD->num_subsystems;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  vec_d tempVec;
  sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_paths * sizeof(sortStruct_d));

  init_vec_d(tempVec, 0);
 
  // print header for level
  printf("\nSorting points for stage %d of %d: %d path%s to sort.\n", stage, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Sorting points for stage %d.\n", stage);
  fprintf(OUT[0], "*****************************************************\n");
 
  // determine if each of the paths is rank deficient
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, finite, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();
 
    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sorting %d of %d\n", i, num_paths);

    // only check the ones that were successful
    if (ED_copy_d[oid].EqD->stageData_d[stage].endPt_retVals[i] == 0)
    { // setup for evaluation
      ED_copy_d[oid].EqD->curr_stage_num = stage;
      ED_copy_d[oid].EqD->curr_path_num = i;

      if (T[oid].MPType == 2)
      {
        ED_copy_d[oid].BED_mp->EqD->curr_stage_num = stage;
        ED_copy_d[oid].BED_mp->EqD->curr_path_num = i;
      }

      // setup PD_d[oid]
      if (ED_copy_d[oid].EqD->stageData_d[stage].useIntrinsicSlice)
      { // sort using the intrinsic coordinates
        point_cp_d(PD_d[oid].point, ED_copy_d[oid].EqD->stageData_d[stage].endPts_in[i]);
      }
      else
      { // sort using original coordinates
        point_cp_d(PD_d[oid].point, ED_copy_d[oid].EqD->stageData_d[stage].endPts[i]);
      }
      set_zero_d(PD_d[oid].time);

      fprintf(OUT[oid], "Path number: %d\n", i);

      // determine if it is rank deficient, finite
      eqbyeqStageSortEndpoint_d(&rankDef, &finite, ED_copy_d[oid].EqD->stageData_d[stage].higherDim[i], &ED_copy_d[oid], ED_copy_d[oid].BED_mp, stage, i, &ED_copy_d[oid].EqD->stageData_d[stage].condition_nums[i], &T[oid], OUT[oid], &PD_d[oid], NULL, 52, ptr_to_eval_d, ptr_to_eval_mp, change_prec); 

      // copy back 
      if (ED_copy_d[oid].EqD->stageData_d[stage].useIntrinsicSlice)
      { // store to the intrinsic coordinates
        point_cp_d(ED_copy_d[oid].EqD->stageData_d[stage].endPts_in[i], PD_d[oid].point);
      }
      else
      { // store using original coordinates
        point_cp_d(ED_copy_d[oid].EqD->stageData_d[stage].endPts[i], PD_d[oid].point);
      }

      // check for success
      if (T[oid].regen_remove_inf && !finite)
      { // dehom point is infinite
        ED_copy_d[oid].EqD->stageData_d[stage].endPt_types[i] = retVal_going_to_infinity;
      }
      else if (T[oid].regen_higher_dim_check && ED_copy_d[oid].EqD->stageData_d[stage].higherDim[i])
      { // it lies on a higher dimensional component
        ED_copy_d[oid].EqD->stageData_d[stage].endPt_types[i] = retVal_higher_dim;
      }
      else if (!rankDef)
      { // classify as non-singular
        ED_copy_d[oid].EqD->stageData_d[stage].endPt_types[i] = MOVE_TO_NEXT;
      }
      else
      { // classify as singular
        ED_copy_d[oid].EqD->stageData_d[stage].endPt_types[i] = DO_NOT_MOVE_TO_NEXT;
      }
    }
    else
    { // path was not a success - copy over error code
      ED_copy_d[oid].EqD->stageData_d[stage].endPt_types[i] = ED_copy_d[oid].EqD->stageData_d[stage].endPt_retVals[i];
    }

    // setup sortPts
    sortPts[i].path_num = i;
    if (ED_copy_d[oid].EqD->stageData_d[stage].useIntrinsicSlice)
    { // setup using the intrinsic coordinates
      sortPts[i].norm = infNormVec_d(ED_copy_d[oid].EqD->stageData_d[stage].endPts_in[i]);
    }
    else
    { // setup using original coordinates
      sortPts[i].norm = infNormVec_d(ED_copy_d[oid].EqD->stageData_d[stage].endPts[i]);
    }
  }

  // sort the structure - use qsort to make comparisons efficient
  qsort(sortPts, num_paths, sizeof(sortStruct_d), sort_order_d);

  // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (ED_d->EqD->stageData_d[stage].endPt_types[indexI] == MOVE_TO_NEXT)
    { // compare against successful paths to see if it is equal to any other path
      j = i + 1;
      while ((j < num_paths) && (sortPts[j].norm - sortPts[i].norm < final_tol))
      {
        indexJ = sortPts[j].path_num;
        if (ED_d->EqD->stageData_d[stage].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          if (ED_d->EqD->stageData_d[stage].useIntrinsicSlice)
          { // subtract using intrinsic coordinates
            vec_sub_d(tempVec, ED_d->EqD->stageData_d[stage].endPts_in[indexI], ED_d->EqD->stageData_d[stage].endPts_in[indexJ]);
          }
          else
          { // subtract using extrinsic coordinates
            vec_sub_d(tempVec, ED_d->EqD->stageData_d[stage].endPts[indexI], ED_d->EqD->stageData_d[stage].endPts[indexJ]);
          }

          if (infNormVec_d(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            ED_d->EqD->stageData_d[stage].endPt_types[indexI] = ED_d->EqD->stageData_d[stage].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
        j++;
      }
    }

    // add to count
    if (ED_d->EqD->stageData_d[stage].endPt_types[indexI] == retVal_going_to_infinity || ED_d->EqD->stageData_d[stage].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (ED_d->EqD->stageData_d[stage].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (ED_d->EqD->stageData_d[stage].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (ED_d->EqD->stageData_d[stage].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  ED_d->EqD->stageData_d[stage].num_sing = num_sing;
  ED_d->EqD->stageData_d[stage].num_nonsing = num_nonsing;
  ED_d->EqD->stageData_d[stage].num_bad = num_bad;
  ED_d->EqD->stageData_d[stage].num_inf = num_inf;
  ED_d->EqD->stageData_d[stage].num_higher_dim = num_higher_dim;

  // clear the memory
  free(sortPts);
  clear_vec_d(tempVec);

  return;
}

void eqbyeqStageSortEndpoint_d(int *rankDef, int *finite, int higherDim, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int stage_num, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoint - rankDef, finite & higherDim       *
\***************************************************************/
{ // initialize
  *rankDef = 0; // we already know higherDim
  *finite = 1;

  // check to see if finite
  if (T->regen_remove_inf)
  { // determine if the endpoint is finite
    *finite = determineEqbyEqFinite(T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, ED, ED_mp, stage_num, 1);
  }

  // determine if we need to do the rank deficient test (which also refines the endpoint)
  if ((T->regen_higher_dim_check && higherDim == 0) && *finite == 1)
  { // determine if it is rank deficient 
    *rankDef = determineRankDef(condNum, T->final_tolerance, NULL, 52, endPt_d, endPt_mp, endPt_prec, T, OUT, ED, ED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec);

    // check to see if finite
    if (T->regen_remove_inf)
    { // determine if the endpoint is finite
      *finite = determineEqbyEqFinite(T->finiteThreshold, endPt_d, endPt_mp, endPt_prec, ED, ED_mp, stage_num, 1);
    }

    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d Rank Def: %d", higherDim, *rankDef);
    }
    else
    {
      fprintf(OUT, "Rank Def: %d", *rankDef);
    }
    if (T->regen_remove_inf)
      fprintf(OUT, " Finite: %d", *finite);
    fprintf(OUT, " CN: %e\n", *condNum);
  }
  else
  { // we know that we are removing this end point
    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d", higherDim);
      if (T->regen_remove_inf)
        fprintf(OUT, " Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
    else
    {
      if (T->regen_remove_inf)
        fprintf(OUT, "Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
  }

  return;
}

int determineEqbyEqFinite(double maxNorm, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int num, int isStage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether it is finite or not                    *
* NOTES:                                                        *
\***************************************************************/
{
  int finite = 0;

  if (Pt_prec < 64)
  { // check using _d
    vec_d orig_d, dehom_d;
    init_vec_d(orig_d, 0); init_vec_d(dehom_d, 0);

    if (isStage)
    { // use stageData_d
      if (ED->EqD->stageData_d[num].useIntrinsicSlice)
      { // convert to original variables
        intrinsicToExtrinsic_d(orig_d, Pt_d->point, ED->EqD->stageData_d[num].B, ED->EqD->stageData_d[num].p);
      }
      else
      { // already ssetup in original variables
        point_cp_d(orig_d, Pt_d->point);
      }
    }
    else
    { // use witnessData_d - convert to original variables
      intrinsicToExtrinsic_d(orig_d, Pt_d->point, ED->EqD->witnessData_d[num].B, ED->EqD->witnessData_d[num].p);
    }
    // find dehom_d
    getDehomPoint_d(dehom_d, orig_d, orig_d->size, &ED->preProcData);

    if (infNormVec_d(dehom_d) < maxNorm)
      finite = 1;
    else 
      finite = 0;

    clear_vec_d(orig_d); clear_vec_d(dehom_d);
  }
  else
  { // check using _mp
    vec_mp orig_mp, dehom_mp;
    init_vec_mp2(orig_mp, 0, Pt_prec); init_vec_mp2(dehom_mp, 0, Pt_prec);

    if (isStage)
    { // use stageData_mp
      if (ED_mp->EqD->stageData_mp[num].useIntrinsicSlice)
      { // convert to original variables
        intrinsicToExtrinsic_mp(orig_mp, Pt_mp->point, ED_mp->EqD->stageData_mp[num].B, ED_mp->EqD->stageData_mp[num].p);
      }
      else
      { // already ssetup in original variables
        point_cp_mp(orig_mp, Pt_mp->point);
      }
    }
    else
    { // use witnessData_mp - convert to original variables
      intrinsicToExtrinsic_mp(orig_mp, Pt_mp->point, ED_mp->EqD->witnessData_mp[num].B, ED_mp->EqD->witnessData_mp[num].p);
    }
    // find dehom_mp
    getDehomPoint_mp(dehom_mp, orig_mp, orig_mp->size, &ED_mp->preProcData);

    if (infNormVec_mp(dehom_mp) < maxNorm)
      finite = 1;
    else
      finite = 0;

    clear_vec_mp(orig_mp); clear_vec_mp(dehom_mp);
  }

  return finite;
}

int determineEqbyEqHigherDim(double tol, double ratio, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int num, int isStage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether it satisfies later functions or not    *
* NOTES:                                                        *
\***************************************************************/
{
  int i, higherDim = 0, startFunc, endFunc, numFuncs = (Pt_prec < 64 ? ED->EqD->num_funcs : ED_mp->EqD->num_funcs);

  // setup startFunc & endFunc
  if (isStage)
  { // use stageData
    startFunc = 0;
    if (Pt_prec < 64)
      endFunc = ED->EqD->stageData_d[num].depth_x + ED->EqD->stageData_d[num].depth_y;
    else
      endFunc = ED_mp->EqD->stageData_mp[num].depth_x + ED_mp->EqD->stageData_mp[num].depth_y;
  }
  else
  { // use witnessData
    if (Pt_prec < 64)
    {
      startFunc = ED->EqD->witnessData_d[num].startFunction;
      endFunc = ED->EqD->witnessData_d[num].startFunction + ED->EqD->witnessData_d[num].depth;
    }
    else
    {
      startFunc = ED_mp->EqD->witnessData_mp[num].startFunction;
      endFunc = ED_mp->EqD->witnessData_mp[num].startFunction + ED_mp->EqD->witnessData_mp[num].depth;
    }
  }

  // see if there is a possibility of being on a higher dim component
  if (startFunc == 0 && endFunc == numFuncs)
  { // no 'extra' functions
    higherDim = 0;
  }
  else 
  { // there are 'extra' functions to test
    int *isZero = (int *)bmalloc(numFuncs * sizeof(int));

    if (Pt_prec < 64)
    { // check using _d
      vec_d orig_d, orig_last_d;
      eval_struct_d e;
      init_vec_d(orig_d, 0);
      init_vec_d(orig_last_d, 0);
      init_eval_struct_d(e, 0, 0, 0);

      if (isStage)
      { // use StageData
        if (ED->EqD->stageData_d[num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_d(orig_d, Pt_d->point, ED->EqD->stageData_d[num].B, ED->EqD->stageData_d[num].p);
        }
        else
        { // copy
          point_cp_d(orig_d, Pt_d->point);
        }

        // evaluate the 'square' function - only need from endFunc to num_funcs since startFunc == 0
        eqbyeq_square_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_d, Pt_d->time, ED, endFunc, numFuncs);

        // setup orig_last
        if (last_approx_prec < 64)
        { // compute orig_last_d
          if (ED->EqD->stageData_d[num].useIntrinsicSlice)
          { // convert to original coordinates first
            intrinsicToExtrinsic_d(orig_last_d, last_approx_d, ED->EqD->stageData_d[num].B, ED->EqD->stageData_d[num].p);
          }
          else
          { // copy
            point_cp_d(orig_last_d, last_approx_d);
          }
        }
        else
        { // convert to _d
          point_mp_to_d(orig_last_d, last_approx_mp);

          if (ED->EqD->stageData_d[num].useIntrinsicSlice)
          { // convert to original coordinates first
            intrinsicToExtrinsic_d(orig_last_d, orig_last_d, ED->EqD->stageData_d[num].B, ED->EqD->stageData_d[num].p);
          }
        }

        // evaluate the 'square' function - only need from endFunc to num_funcs since startFunc == 0
        eqbyeq_square_eval_d(orig_d, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_d, Pt_d->time, ED, endFunc, numFuncs);

        // compare
        nonsolutions_check_compare_d(isZero, e.funcVals, orig_d, endFunc, numFuncs, tol, ratio);

        // check
        for (i = endFunc; i < numFuncs; i++)
          if (isZero[i])
            higherDim = 1;
      }
      else
      { // use witnessData_d - convert to original coordinates first
        intrinsicToExtrinsic_d(orig_d, Pt_d->point, ED->EqD->witnessData_d[num].B, ED->EqD->witnessData_d[num].p);

        // evaluate the 'square' function - might as well find all functions since startFunc may not be 0
        eqbyeq_square_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_d, Pt_d->time, ED, 0, numFuncs);

        if (last_approx_prec < 64)
        { // use _d
          intrinsicToExtrinsic_d(orig_last_d, last_approx_d, ED->EqD->witnessData_d[num].B, ED->EqD->witnessData_d[num].p);
        }
        else
        { // convert to _d
          point_mp_to_d(orig_last_d, last_approx_mp);
          intrinsicToExtrinsic_d(orig_last_d, orig_last_d, ED->EqD->witnessData_d[num].B, ED->EqD->witnessData_d[num].p);
        }

        // evaluate the 'square' function - might as well find all functions since startFunc may not be 0
        eqbyeq_square_eval_d(orig_d, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_d, Pt_d->time, ED, 0, numFuncs);

        // compare
        nonsolutions_check_compare_d(isZero, e.funcVals, orig_d, 0, numFuncs, tol, ratio);

        // check the function above
        for (i = 0; i < startFunc; i++)
          if (isZero[i])
            higherDim = 1;

        // check the functions below
        for (i = endFunc; i < numFuncs; i++)
          if (isZero[i])
            higherDim = 1;
      }

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

      if (isStage)
      { // use StageData
        if (ED_mp->EqD->stageData_mp[num].useIntrinsicSlice)
        { // convert to original coordinates first
          intrinsicToExtrinsic_mp(orig_mp, Pt_mp->point, ED_mp->EqD->stageData_mp[num].B, ED_mp->EqD->stageData_mp[num].p);
        }
        else
        { // copy
          point_cp_mp(orig_mp, Pt_mp->point);
        }

        // evaluate the 'square' function - only need from endFunc to num_funcs since startFunc == 0
        eqbyeq_square_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_mp, Pt_mp->time, ED_mp, endFunc, numFuncs);

        // setup orig_last
        if (last_approx_prec < 64)
        { // convert to _mp
          point_d_to_mp(orig_last_mp, last_approx_d);

          if (ED_mp->EqD->stageData_mp[num].useIntrinsicSlice)
          { // convert to original coordinates first
            intrinsicToExtrinsic_mp(orig_last_mp, orig_last_mp, ED_mp->EqD->stageData_mp[num].B, ED_mp->EqD->stageData_mp[num].p);
          }
        }
        else
        { // compute orig_last_mp
          if (ED_mp->EqD->stageData_mp[num].useIntrinsicSlice)
          { // convert to original coordinates first
            intrinsicToExtrinsic_mp(orig_last_mp, last_approx_mp, ED_mp->EqD->stageData_mp[num].B, ED_mp->EqD->stageData_mp[num].p);
          }
          else
          { // copy
            point_cp_mp(orig_last_mp, last_approx_mp);
          }
        }

        // evaluate the 'square' function - only need from endFunc to num_funcs since startFunc == 0
        eqbyeq_square_eval_mp(orig_mp, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_mp, Pt_mp->time, ED_mp, endFunc, numFuncs); 

        // compare
        nonsolutions_check_compare_mp(isZero, e.funcVals, orig_mp, endFunc, numFuncs, tol, ratio);

        // check
        for (i = endFunc; i < numFuncs; i++)
          if (isZero[i])
            higherDim = 1;
      }
      else
      { // use witnessData_mp - convert to original coordinates first
        intrinsicToExtrinsic_mp(orig_mp, Pt_mp->point, ED_mp->EqD->witnessData_mp[num].B, ED_mp->EqD->witnessData_mp[num].p);

        // evaluate the 'square' function - might as well find all functions since startFunc may not be 0
        eqbyeq_square_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, orig_mp, Pt_mp->time, ED_mp, 0, numFuncs);

        if (last_approx_prec < 64)
        { // convert to _mp
          point_d_to_mp(orig_last_mp, last_approx_d);
          intrinsicToExtrinsic_mp(orig_last_mp, orig_last_mp, ED_mp->EqD->witnessData_mp[num].B, ED_mp->EqD->witnessData_mp[num].p);
        }
        else
        { // use _mp
          intrinsicToExtrinsic_mp(orig_last_mp, last_approx_mp, ED_mp->EqD->witnessData_mp[num].B, ED_mp->EqD->witnessData_mp[num].p);
        }

        // evaluate the 'square' function - might as well find all functions since startFunc may not be 0
        eqbyeq_square_eval_mp(orig_mp, e.parVals, e.parDer, e.Jv, e.Jp, orig_last_mp, Pt_mp->time, ED_mp, 0, numFuncs);

        // compare
        nonsolutions_check_compare_mp(isZero, e.funcVals, orig_mp, 0, numFuncs, tol, ratio);

        // check the function above
        for (i = 0; i < startFunc; i++)
          if (isZero[i])
            higherDim = 1;

        // check the functions below
        for (i = endFunc; i < numFuncs; i++)
          if (isZero[i])
            higherDim = 1;
      }

      clear_point_mp(orig_mp);
      clear_point_mp(orig_last_mp);
      clear_eval_struct_mp(e);
    }

    free(isZero);
  }
 
  return higherDim;
}

int printStageFooter_d(basic_eval_data_d *ED, int stage_num, int path_num, point_data_d *endPoint, point_d orig_vars, point_d orig_last_d, point_d dehom_d, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for path stage_num & path_num        *
\***************************************************************/
{
  int i, true_failure = 0, retVal_out = 0, isNumber = 1, isSoln = 0, dehomInf = 0;

  if (retVal_in)
  {
    if (T->regen_remove_inf && infNormVec_d(dehom_d) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else if (retVal_in == retVal_sharpening_singular_endpoint)
    { // sharpening determined it was singular
      retVal_out = retVal_sharpening_singular_endpoint;
      true_failure = 1;
    }
    else if (endPoint->time->r < T->minTrackT)
    { // display a warning message if it is a finite path
      if (retVal_in == retVal_reached_minTrackT)
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d for stage %d reached the minimum value of T (%e < %e).\n", path_num, stage_num, endPoint->time->r, T->minTrackT);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d for stage %d reached the minimum value of T (%e < %e).\n", path_num, stage_num, endPoint->time->r, T->minTrackT);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }
      else
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d for stage %d had final T of %e but had retVal %d.\n", path_num, stage_num, endPoint->time->r, retVal_in);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d for stage %d had final T of %e but had retVal %d.\n", path_num, stage_num, endPoint->time->r, retVal_in);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }

      // consider it a success
      retVal_out = 0;
    }
    else
    { // the path had an error code
      retVal_out = retVal_in;
      true_failure = 1; // this path truly failed
      printResultOfPath(OUT, retVal_in, T);
    }
  }
  else // no error code
  { // check to make sure that dehom is not at infinity
    if (T->regen_remove_inf && infNormVec_d(dehom_d) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else
    {
      retVal_out = retVal_in;
      printResultOfPath(OUT, retVal_in, T);
    }
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
      if (isnan(endPoint->point->coord[i].r) || isnan(endPoint->point->coord[i].i) || isinf(endPoint->point->coord[i].r) || isinf(endPoint->point->coord[i].i))
        isNumber = 0;

    if (!isNumber)
    {
      retVal_out = retVal_NAN;
      true_failure = 1; // this path truly failed
    }
    else if (stage_num + 1 == ED->EqD->num_subsystems)
    { // check that the point of original variables satisfies the original system if it is a number and we are at the top level
      isSoln = nonsolutions_check_d(ED->squareSystem.size_f, ED->squareSystem.size_r, orig_vars, orig_last_d, endPoint->time, T->funcResTol, T->ratioTol, ED->squareSystem.Prog);

      if (!isSoln)
      {
        retVal_out = retVal_Bertini_Junk;
        true_failure = 1; // this path truly failed
      }
    }
  }

  // print the path footer to OUT
  printPathFooterOut_d(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, dehom_d, T, ED->squareSystem.Prog, ED->preProcData.num_var_gp, 0);

  // determine if we need to print to RAWOUT & FAIL
  if (stage_num + 1 == ED->EqD->num_subsystems)
  { // this is the last subsystem, so we need to either print to FAIL or RAWOUT
    if (retVal_out && !dehomInf)
    { // path was not successful, so update the number of failures and print the info to FAIL
      if (true_failure)
        trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_d(FAIL, endPoint, dehom_d, path_num, retVal_out, isNumber, !isSoln, trackCount, T);

      if (retVal_out == retVal_Bertini_Junk)
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
    { // path was successful, so update the number of successes and print the info to RAWOUT
      trackCount->successes++;

      if (retVal_in == retVal_reached_minTrackT) // success but convergence may not be quite accurate enough
        printSuccess_d(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, func_residual, newton_error, t_val_sample, error_sample, -1);
      else
        printSuccess_d(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, func_residual, newton_error, t_val_sample, error_sample, 1);
    }
  }

  return retVal_out;
}

int printWitnessFooter_d(basic_eval_data_d *ED, int subsystem_num, int path_num, point_data_d *endPoint, point_d orig_vars, point_d orig_last_d, point_d dehom, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for path subsystem_num & path_num    *
\***************************************************************/
{
  int i, true_failure = 0, retVal_out = 0, isNumber = 1, isSoln = 0, dehomInf = 0;

  if (retVal_in)
  {
    if (T->regen_remove_inf && infNormVec_d(dehom) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else if (retVal_in == retVal_sharpening_singular_endpoint)
    { // sharpening determined it was singular
      retVal_out = retVal_sharpening_singular_endpoint;
      true_failure = 1;
    }
    else if (endPoint->time->r < T->minTrackT)
    { // display a warning message if it is a finite path
      if (retVal_in == retVal_reached_minTrackT)
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on subsystem %d reached the minimum value of T (%e < %e).\n", path_num, subsystem_num, endPoint->time->r, T->minTrackT);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on subsystem %d reached the minimum value of T (%e < %e).\n", path_num, subsystem_num, endPoint->time->r, T->minTrackT);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }
      else
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on subsystem %d had final T of %e but had retVal %d.\n", path_num, subsystem_num, endPoint->time->r, retVal_in);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on subsystem %d had final T of %e but had retVal %d.\n", path_num, subsystem_num, endPoint->time->r, retVal_in);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }

      // consider it a success
      retVal_out = 0;
    }
    else
    { // the path had an error code
      retVal_out = retVal_in;
      true_failure = 1; // this path truly failed
      printResultOfPath(OUT, retVal_in, T);
    }
  }
  else // no error code
  { // check to make sure that dehom is not at infinity
    if (T->regen_remove_inf && infNormVec_d(dehom) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 0;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else
    {
      retVal_out = retVal_in;
      printResultOfPath(OUT, retVal_in, T);
    }
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
      if (isnan(endPoint->point->coord[i].r) || isnan(endPoint->point->coord[i].i) || isinf(endPoint->point->coord[i].r) || isinf(endPoint->point->coord[i].i))
        isNumber = 0;

    if (!isNumber)
    {
      retVal_out = retVal_NAN;
      true_failure = 1; // this path truly failed
    }
    else if (ED->EqD->num_subsystems == 1)
    { // this is the only subsystem, so we need to verify that the original variables satisfies the original function
      isSoln = nonsolutions_check_d(ED->squareSystem.size_f, ED->squareSystem.size_r, orig_vars, orig_last_d, endPoint->time, T->funcResTol, T->ratioTol, ED->squareSystem.Prog);
   
      if (!isSoln)
      {
        retVal_out = retVal_Bertini_Junk;
        true_failure = 1; // this path truly failed 
      }
    }
  }

  // print the path footer to OUT
  printPathFooterOut_d(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, dehom, T, ED->squareSystem.Prog, ED->preProcData.num_var_gp, 0);

  // determine if we need to print to RAWOUT & FAIL
  if (ED->EqD->num_subsystems == 1)
  { // this is the only subsystem, so we need to either print to FAIL or RAWOUT
    if (retVal_out && !dehomInf)
    { // path was not successful, so update the number of failures and print the info to FAIL
      if (true_failure)
        trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_d(FAIL, endPoint, dehom, path_num, retVal_out, isNumber, !isSoln, trackCount, T);

      if (retVal_out == retVal_Bertini_Junk)
      { // print to NONSOLN
        for (i = 0; i < dehom->size; i++)
        {
          print_d(NONSOLN, 0, &dehom->coord[i]);
          fprintf(NONSOLN, "\n");
        }
        fprintf(NONSOLN, "\n");
      }
    }
    else
    { // path was successful, so update the number of successes and print the info to RAWOUT
      trackCount->successes++;

      if (retVal_in == retVal_reached_minTrackT) // success but convergence may not be quite accurate enough
        printSuccess_d(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, func_residual, newton_error, t_val_sample, error_sample, -1);
      else
        printSuccess_d(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, func_residual, newton_error, t_val_sample, error_sample, 1);
    }
  }

  return retVal_out;
}

///////// MULTI PRECISION ////////////////

void eqbyeq_track_mp(FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_mp *ED, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does zero dimensional tracking via eq-by-eq of SVW     *
\***************************************************************/
{
  int i, num_paths, num_vars, num_crossings, max = max_threads();
  point_data_mp *tempPoint = (point_data_mp *)bmalloc(max * sizeof(point_data_mp));
  int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  char outName[] = "output", midName[] = "midout", rawName[] = "rawout", failName[] = "fail", nonName[] = "nonsolutions";
  FILE *MIDOUT = fopen(midFile, "w"), *NONSOLN = NULL;
  if (MIDOUT == NULL)
  {
    printf("ERROR: '%s' is not a valid name for a file!\n", midFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  if (!ED->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen(nonName, "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // pointers for OpenMP tracking
  tracker_config_t *T_copy = NULL;
  basic_eval_data_mp *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, **NONSOLN_copy = NULL;

  // setup the structures
  setup_eqbyeq_omp_mp(max, &trackCount_copy, trackCount, &OUT_copy, OUT, outName, &RAWOUT_copy, RAWOUT, rawName, &MIDOUT_copy, MIDOUT, midName, &FAIL_copy, FAIL, failName, &NONSOLN_copy, NONSOLN, nonName, &T_copy, T, &BED_copy, ED);
  for (i = 0; i < max; i++)
  {
    init_point_data_mp(&tempPoint[i], 0);
  }

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // generate the witness sets for each of the subsystems
  ptr_to_eval = &witness_eqbyeq_eval_mp;

  for (i = 0; i < ED->EqD->num_subsystems; i++)
  { // find the number of paths for this subsystem
    num_paths = ED->EqD->witnessData_mp[i].num_paths;
    num_vars = ED->EqD->witnessData_mp[i].depth;

    // track each path for this subsystem
    eqbyeqWitnessTrack_mp(max, pathMod, BED_copy, T_copy, OUT_copy, MIDOUT_copy, RAWOUT_copy, FAIL_copy, NONSOLN_copy, i, trackCount_copy, ptr_to_eval, eqbyeq_witness_dehom);

    // combine the files
    combine_omp_file(OUT, &OUT_copy, max);
    combine_omp_file(MIDOUT, &MIDOUT_copy, max);
    combine_omp_file(RAWOUT, &RAWOUT_copy, max);
    combine_omp_file(FAIL, &FAIL_copy, max);

    // make sure that no path crossing happened when finding this subsystem, if this is not the only subsystem
    if (ED->EqD->num_subsystems > 1)
    { // close MIDOUT and check for path crossings
      fclose(MIDOUT);
      num_crossings = 0;
      midpoint_checker(num_paths, num_vars, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // setup OUT_copy for sorting
    setup_omp_file(&OUT_copy, OUT, outName, max);

    // sort the endpoints
    sortWitnessEndpoints_mp(max, pathMod, tempPoint, BED_copy, ED, T_copy, OUT_copy, i, target_tol, ptr_to_eval);

    // combine OUT
    combine_omp_file(OUT, &OUT_copy, max);

    // if this is not the last subsystem, reopen the files
    if (i + 1 < ED->EqD->num_subsystems)
    { // reopen the files
      MIDOUT = fopen(midFile, "w");
      setup_omp_file(&MIDOUT_copy, MIDOUT, midName, max);
      setup_omp_file(&OUT_copy, OUT, outName, max);
      setup_omp_file(&RAWOUT_copy, RAWOUT, rawName, max);
      setup_omp_file(&FAIL_copy, FAIL, failName, max);
    }
  }

  // setup the first stage
  setupEqbyEqFirstStage_mp(ED->EqD);
  setup_omp_eqbyeq_first_stage_mp(max, BED_copy, ED->EqD);

  // clear the first witness data
  clearEqbyEqFirstWitnessData_mp(ED->EqD);
  clear_omp_eqbyeq_witness_mp(max, BED_copy, 0);

  for (i = 1; i < ED->EqD->num_subsystems; i++)
  { // setup stage i
    setupEqbyEqNextStage_mp(ED, i);
    setup_omp_eqbyeq_next_stage_mp(max, BED_copy, ED->EqD, i);

    // clear the stage date from stage i - 1 since it is no longer needed
    clearEqbyEqStageData_mp(ED->EqD, i - 1);
    clear_omp_eqbyeq_stage_mp(max, BED_copy, i - 1);

    // clear the witness data from subsystem i since it is no longer needed
    clearEqbyEqWitnessData_mp(ED->EqD, i);
    clear_omp_eqbyeq_witness_mp(max, BED_copy, i);

    // setup the files to track stage i
    MIDOUT = fopen(midFile, "w");
    setup_omp_file(&MIDOUT_copy, MIDOUT, midName, max);
    setup_omp_file(&OUT_copy, OUT, outName, max);
    setup_omp_file(&RAWOUT_copy, RAWOUT, rawName, max);
    setup_omp_file(&FAIL_copy, FAIL, failName, max);

    // find the number of paths for this next stage
    num_paths = ED->EqD->stageData_mp[i].num_paths;
    if (ED->EqD->stageData_mp[i].useIntrinsicSlice)
      num_vars = ED->EqD->stageData_mp[i].depth_x + ED->EqD->stageData_mp[i].depth_y;
    else
      num_vars = 2 * ED->EqD->num_vars;

    // setup the evaluators for stage tracking
    ptr_to_eval = &standard_eqbyeq_eval_mp;

    // track each path for this stage
    eqbyeqStageTrack_mp(max, pathMod, BED_copy, T_copy, OUT_copy, RAWOUT_copy, MIDOUT_copy, FAIL_copy, NONSOLN_copy, i, trackCount_copy, ptr_to_eval, eqbyeq_stage_dehom);

    // combine all of the files
    combine_omp_file(OUT, &OUT_copy, max);
    combine_omp_file(MIDOUT, &MIDOUT_copy, max);
    combine_omp_file(RAWOUT, &RAWOUT_copy, max);
    combine_omp_file(FAIL, &FAIL_copy, max);

    // make sure that no path crossing happened when finding this stage, if this is not the last stage
    if (i + 1 < ED->EqD->num_subsystems)
    { // close MIDOUT and check for path crossings
      fclose(MIDOUT);
      num_crossings = 0;
      midpoint_checker(num_paths, num_vars, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
    }

    // setup OUT_copy for sorting
    setup_omp_file(&OUT_copy, OUT, outName, max);

    // setup the evaluators for stage sorting
    ptr_to_eval = &stage_sort_eqbyeq_eval_mp;

    // sort the endpoints
    sortStageEndpoints_mp(max, pathMod, tempPoint, BED_copy, ED, T_copy, OUT_copy, i, target_tol, ptr_to_eval);

    // combine OUT
    combine_omp_file(OUT, &OUT_copy, max);
  }
  // close MIDOUT
  fclose(MIDOUT);

  if (!ED->squareSystem.noChanges)
  { // complete NONSOLN
    combine_omp_file(NONSOLN, &NONSOLN_copy, max);
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // combine data and clear copies - just need to delete the extra files
  clear_eqbyeq_omp_mp(max, &trackCount_copy, trackCount, outName, rawName, midName, failName, nonName, &T_copy, &BED_copy);

  // set the total number tracked on last stage to trackCount
  trackCount->numPoints = ED->EqD->stageData_mp[ED->EqD->num_subsystems - 1].num_paths;

  // free tempPoint
  for (i = 0; i < max; i++)
  {
    clear_point_data_mp(&tempPoint[i]);
  }
  free(tempPoint);

  return;
}

void setup_eqbyeq_omp_mp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, FILE ***NONSOLN_copy, FILE *NONSOLN, char *nonName, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_mp **BED_copy, basic_eval_data_mp *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup everything needed to do eq-by-eq tracking        *
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
  setup_omp_file(NONSOLN_copy, NONSOLN, nonName, max_threads);

  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED;
  }
  else // max_threads > 1
  { // allocate memory
    *trackCount_copy = (trackingStats *)bmalloc(max_threads * sizeof(trackingStats));
    *T_copy = (tracker_config_t *)bmalloc(max_threads * sizeof(tracker_config_t));
    *BED_copy = (basic_eval_data_mp *)bmalloc(max_threads * sizeof(basic_eval_data_mp));

    // copy T, ED, & trackCount
    for (i = 0; i < max_threads; i++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[i], T);
      // copy ED
      cp_basic_eval_data_mp(&(*BED_copy)[i], ED);
      // setup EqD
      setup_omp_eqbyeq_mp(&(*BED_copy)[i], ED->EqD);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[i]);
      (*trackCount_copy)[i].numPoints = trackCount->numPoints;
    }
  }

  return;
}

void clear_eqbyeq_omp_mp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, char *nonName, tracker_config_t **T_copy, basic_eval_data_mp **BED_copy)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear everything needed to do eq-by-eq tracking        *
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
    *BED_copy = NULL;
  }
  else if (max_threads > 1)
  { // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);

    // clear the copies of T & ED_d
    for (i = max_threads - 1; i >= 0; i--)
    { // clear BED_copy - 0 since we want to target exactly want needs cleared in the eq-by-eq data
      basic_eval_clear_mp(&(*BED_copy)[i], 0, 1);
      // clear EqD
      clear_omp_eqData_mp(&(*BED_copy)[i]);
      // clear T_copy
      tracker_config_clear(&(*T_copy)[i]);
    }

    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);

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

      // delete nonsolutions
      size = 1 + snprintf(NULL, 0, "%s_%d", nonName, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", nonName, i);
      remove(str);
    }
    free(str);
  }

  return;
}

void eqbyeqWitnessTrack_mp(int max_threads, int pathMod, basic_eval_data_mp ED[], tracker_config_t T[], FILE **OUT, FILE **MIDOUT, FILE **RAWOUT, FILE **FAIL, FILE **NONSOLN, int subsystem_num, trackingStats trackCount[], int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths to find witness set for this subsystem*
* Uses OpenMP if available                                      *
\***************************************************************/
{
  int i, oid, num_paths = ED[0].EqD->witnessData_mp[subsystem_num].num_paths, num_subs = ED[0].EqD->num_subsystems;

  // display messages
  printf("\nFinding witness points for subsystem %d of %d: %d path%s to track.\n", subsystem_num, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Finding witness points for subsystem %d.\n", subsystem_num);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this subsystem
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // track the path
    eqbyeqWitnessTrackPath_mp(&ED[oid], &T[oid], OUT[oid], MIDOUT[oid], RAWOUT[oid], FAIL[oid], NONSOLN[oid], subsystem_num, i, &trackCount[oid], ptr_to_eval, find_dehom);
  }

  return;
}

void eqbyeqWitnessTrackPath_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int subsystem_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in subsystem_num & path_num of ED      *
\***************************************************************/
{
  int num_subs = ED->EqD->num_subsystems;

  endgame_data_t endPt;
  point_data_mp startPt;
  point_mp dehom, orig_last;

  // initialize MP
  init_endgame_data(&endPt, T->Precision);
  init_point_data_mp(&startPt, 0);
  init_point_mp(dehom, 0); init_point_mp(orig_last, 0);

  // seutp for tracking the path
  point_cp_mp(startPt.point, ED->EqD->witnessData_mp[subsystem_num].startPts[path_num]);
  mpf_set_ui(startPt.time->r, 1); mpf_set_ui(startPt.time->i, 0);
  ED->EqD->curr_stage_num = subsystem_num;
  ED->EqD->curr_path_num = path_num;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  // print the header for the point
  printPathHeader_mp(OUT, &startPt, T, path_num, ED, ptr_to_eval);

  zero_dim_track_path_mp(path_num, &endPt, &startPt, OUT, MIDOUT, T, ED, ptr_to_eval, find_dehom);

  // check to see if it should be sharpened - only when this is the only subsystem
  if (num_subs == 1 && endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, NULL, ED, NULL, ptr_to_eval, NULL);
  }

  // store the condition number
  ED->EqD->witnessData_mp[subsystem_num].condition_nums[path_num] = endPt.condition_number;

  // store if higher dimenaional
  ED->EqD->witnessData_mp[subsystem_num].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, NULL, ED, subsystem_num, 0);

  // copy over to the appropriate spot
  point_cp_mp(ED->EqD->witnessData_mp[subsystem_num].endPts_in[path_num], endPt.PD_mp.point);
  set_mp(ED->EqD->witnessData_mp[subsystem_num].finalTs[path_num], endPt.PD_mp.time);

  // convert to the original coordinates
  intrinsicToExtrinsic_mp(orig_last, endPt.last_approx_mp, ED->EqD->witnessData_mp[subsystem_num].B, ED->EqD->witnessData_mp[subsystem_num].p);
  intrinsicToExtrinsic_mp(ED->EqD->witnessData_mp[subsystem_num].endPts[path_num], endPt.PD_mp.point, ED->EqD->witnessData_mp[subsystem_num].B, ED->EqD->witnessData_mp[subsystem_num].p);

  // find dehom_mp
  getDehomPoint_mp(dehom, ED->EqD->witnessData_mp[subsystem_num].endPts[path_num], ED->EqD->witnessData_mp[subsystem_num].endPts[path_num]->size, &ED->preProcData);

  // print the footer to OUT for the point & find the condition number
  ED->EqD->witnessData_mp[subsystem_num].endPt_retVals[path_num] = printWitnessFooter_mp(ED, subsystem_num, path_num, &endPt.PD_mp, ED->EqD->witnessData_mp[subsystem_num].endPts[path_num], orig_last, dehom, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);

  // clear MP
  clear_endgame_data(&endPt)
  clear_point_data_mp(&startPt);
  clear_point_mp(dehom); clear_point_mp(orig_last);

  return;
}

void sortWitnessEndpoints_mp(int max_threads, int pathMod, point_data_mp PD[], basic_eval_data_mp ED_copy[], basic_eval_data_mp *ED, tracker_config_t T[], FILE **OUT, int subsystem, double final_tol, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points found for current subsystem   *
\***************************************************************/
// Still doing everything as intrinsic linears
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, cont, oid, rankDef, finite, indexI = 0, indexJ = 0, num_paths = ED->EqD->witnessData_mp[subsystem].num_paths, num_subs = ED->EqD->num_subsystems;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  mpf_t norm_diff;
  vec_mp tempVec;
  sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));

  mpf_init(norm_diff);
  init_vec_mp(tempVec, 0);

  // print header for level
  printf("\nSorting witness points for subsystem %d of %d: %d path%s to sort.\n", subsystem, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Sorting witness points for subsystem %d.\n", subsystem);
  fprintf(OUT[0], "*****************************************************\n");

  // determine if each of the paths is rank deficient
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, finite, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sorting %d of %d\n", i, num_paths);

    // only check the ones that were successful
    if (ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_retVals[i] == 0)
    { // setup for evaluation
      ED_copy[oid].EqD->curr_stage_num = subsystem;
      ED_copy[oid].EqD->curr_path_num = i;

      // setup PD[oid]
      point_cp_mp(PD[oid].point, ED_copy[oid].EqD->witnessData_mp[subsystem].endPts_in[i]);
      set_zero_mp(PD[oid].time);

      fprintf(OUT[oid], "Path number: %d\n", i);
      // determine if it is rank deficient, finite
      eqbyeqWitnessSortEndpoint_mp(&rankDef, &finite, ED_copy[oid].EqD->witnessData_mp[subsystem].higherDim[i], &ED_copy[oid], subsystem, i, &ED_copy[oid].EqD->witnessData_mp[subsystem].condition_nums[i], &T[oid], OUT[oid], &PD[oid], ptr_to_eval);

      // copy back to endPts_in - determineRankDef will update PD[oid] if it is improved
      point_cp_mp(ED_copy[oid].EqD->witnessData_mp[subsystem].endPts_in[i], PD[oid].point);

      // check for success
      if (T[oid].regen_remove_inf && !finite)
      { // dehom point is infinite
        ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_types[i] = retVal_going_to_infinity;
      }
      else if (T[oid].regen_higher_dim_check && ED_copy[oid].EqD->witnessData_mp[subsystem].higherDim[i])
      { // it lies on a higher dimensional component
        ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_types[i] = retVal_higher_dim;
      }
      else if (!rankDef)
      { // classify as non-singular
        ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_types[i] = MOVE_TO_NEXT;
      }
      else
      { // classify as singular
        ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_types[i] = DO_NOT_MOVE_TO_NEXT;
      }
    }
    else
    { // path was not a success - copy over error code
      ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_types[i] = ED_copy[oid].EqD->witnessData_mp[subsystem].endPt_retVals[i];
    }

    // setup sortPts
    mpf_init(sortPts[i].norm);
    sortPts[i].path_num = i;
    infNormVec_mp2(sortPts[i].norm, ED_copy[oid].EqD->witnessData_mp[subsystem].endPts_in[i]);
  }

  // sort the structure - use qsort to make comparisons efficient
  qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

  // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == MOVE_TO_NEXT)
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
        if (cont && ED->EqD->witnessData_mp[subsystem].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          vec_sub_mp(tempVec, ED->EqD->witnessData_mp[subsystem].endPts_in[indexI], ED->EqD->witnessData_mp[subsystem].endPts_in[indexJ]);
          if (infNormVec_mp(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] = ED->EqD->witnessData_mp[subsystem].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
      } while (cont);
    }

    // add to count
    if (ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == retVal_going_to_infinity || ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (ED->EqD->witnessData_mp[subsystem].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  ED->EqD->witnessData_mp[subsystem].num_sing = num_sing;
  ED->EqD->witnessData_mp[subsystem].num_nonsing = num_nonsing;
  ED->EqD->witnessData_mp[subsystem].num_bad = num_bad;
  ED->EqD->witnessData_mp[subsystem].num_higher_dim = num_higher_dim;
  ED->EqD->witnessData_mp[subsystem].num_inf = num_inf;

  // clear the memory
  for (i = num_paths - 1; i >= 0; i--)
  {
    mpf_clear(sortPts[i].norm);
  }
  free(sortPts);

  mpf_clear(norm_diff);
  clear_vec_mp(tempVec);

  return;
}

void eqbyeqWitnessSortEndpoint_mp(int *rankDef, int *finite, int higherDim, basic_eval_data_mp *ED, int subsystem, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_mp *endPt_mp, int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoint - rankDef, finite & higherDim       *
\***************************************************************/
{ // initialize
  *rankDef = 0; // we already know higherDim!
  *finite = 1;

  // check to see if finite
  if (T->regen_remove_inf)
  { // determine if the endpoint is finite
    *finite = determineEqbyEqFinite(T->finiteThreshold, NULL, endPt_mp, T->Precision, NULL, ED, subsystem, 0);
  }

  // determine if we need to do the rank deficient test (which also refines the endpoint)
  if ((T->regen_higher_dim_check && higherDim == 0) && *finite == 1)
  { // determine if it is rank deficient 
    *rankDef = determineRankDef(condNum, T->final_tolerance, NULL, 52, NULL, endPt_mp, T->Precision, T, OUT, NULL, ED, NULL, ptr_to_eval_mp, NULL);

    // check to see if finite
    if (T->regen_remove_inf)
    { // determine if the endpoint is finite
      *finite = determineEqbyEqFinite(T->finiteThreshold, NULL, endPt_mp, T->Precision, NULL, ED, subsystem, 0);
    }

    if (T->regen_higher_dim_check)
    {  
      fprintf(OUT, "Higher Dim'l: %d Rank Def: %d", higherDim, *rankDef);
    }
    else
    {
      fprintf(OUT, "Rank Def: %d", *rankDef);
    }
    if (T->regen_remove_inf)
      fprintf(OUT, " Finite: %d", *finite);
    fprintf(OUT, " CN: %e\n", *condNum);
  }
  else
  { // we know that we are removing this end point
    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d", higherDim);
      if (T->regen_remove_inf)
        fprintf(OUT, " Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
    else
    {
      if (T->regen_remove_inf)
        fprintf(OUT, "Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
  }

  return;
}

void eqbyeqStageTrack_mp(int max_threads, int pathMod, basic_eval_data_mp ED[], tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, FILE **NONSOLN, int stage_num, trackingStats trackCount[], int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks all paths to find the witness points for stage  *
* Uses OpenMP if available                                      *
\***************************************************************/
{
  int i, oid, num_paths = ED[0].EqD->stageData_mp[stage_num].num_paths, num_stages = ED[0].EqD->num_subsystems;

  // display messages
  printf("\nTracking points for stage %d of %d: %d path%s to track.\n", stage_num, num_stages, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Tracking points for stage %d.\n", stage_num);
  fprintf(OUT[0], "*****************************************************\n");

  // track each path for this stage
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, num_paths);

    // track the path
    eqbyeqStageTrackPath_mp(&ED[oid], &T[oid], OUT[oid], RAWOUT[oid], MIDOUT[oid], FAIL[oid], NONSOLN[oid], stage_num, i, &trackCount[oid], ptr_to_eval, find_dehom);
  }

  return;
}

void eqbyeqStageTrackPath_mp(basic_eval_data_mp *ED, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int stage_num, int path_num, trackingStats *trackCount, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in stage_num & path_num of ED          *
\***************************************************************/
{
  int num_subs = ED->EqD->num_subsystems;

  endgame_data_t endPt;
  point_data_mp startPt;
  point_mp dehom, orig_last;

  // initialize MP
  init_endgame_data(&endPt, T->Precision);
  init_point_data_mp(&startPt, 0);
  init_point_mp(dehom, 0); init_point_mp(orig_last, 0);

  // seutp for tracking the path
  point_cp_mp(startPt.point, ED->EqD->stageData_mp[stage_num].startPts[path_num]);
  mpf_set_ui(startPt.time->r, 1); mpf_set_ui(startPt.time->i, 0);
  ED->EqD->curr_stage_num = stage_num;
  ED->EqD->curr_path_num = path_num;
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  // print the header for the point
  printPathHeader_mp(OUT, &startPt, T, path_num, ED, ptr_to_eval);

  zero_dim_track_path_mp(path_num, &endPt, &startPt, OUT, MIDOUT, T, ED, ptr_to_eval, find_dehom);

  // check to see if it should be sharpened - only when this is the last stage
  if (stage_num + 1 == num_subs && endPt.retVal == 0 && T->sharpenDigits > 0)
  { // use the sharpener for after an endgame
    sharpen_endpoint_endgame(&endPt, T, OUT, NULL, ED, NULL, ptr_to_eval, NULL);
  }

  // store the condition number
  ED->EqD->stageData_mp[stage_num].condition_nums[path_num] = endPt.condition_number;

  // store if higher dimenaional
  ED->EqD->stageData_mp[stage_num].higherDim[path_num] = determineEqbyEqHigherDim(T->funcResTol, T->ratioTol, &endPt.PD_d, &endPt.PD_mp, endPt.prec, endPt.last_approx_d, endPt.last_approx_mp, endPt.last_approx_prec, NULL, ED, stage_num, 1);

  // copy over to the appropriate spot
  set_mp(ED->EqD->stageData_mp[stage_num].finalTs[path_num], endPt.PD_mp.time);

  if (ED->EqD->stageData_mp[stage_num].useIntrinsicSlice)
  { // store the intrinsic endpoint and its corresponding point in the original variables
    intrinsicToExtrinsic_mp(orig_last, endPt.last_approx_mp, ED->EqD->stageData_mp[stage_num].B0, ED->EqD->stageData_mp[stage_num].p0);
    orig_last->size /= 2; // remove the bottom half of the coordinates
    intrinsicToExtrinsic_mp(dehom, endPt.PD_mp.point, ED->EqD->stageData_mp[stage_num].B0, ED->EqD->stageData_mp[stage_num].p0);
    dehom->size /= 2; // remove the bottom half of the coordinates

    point_cp_mp(ED->EqD->stageData_mp[stage_num].endPts[path_num], dehom);

    // convert to intrinsic coordinates
    mat_mp B_transpose;
    init_mat_mp(B_transpose, 0, 0);

    transpose_mp(B_transpose, ED->EqD->stageData_mp[stage_num].B);
    extrinsicToIntrinsic_mp(ED->EqD->stageData_mp[stage_num].endPts_in[path_num], dehom, B_transpose, ED->EqD->stageData_mp[stage_num].p);
  
    clear_mat_mp(B_transpose);
  }
  else
  { // use the top set of coordinates so that we can store the original coordinates
    point_cp_mp(orig_last, endPt.last_approx_mp);
    orig_last->size /= 2; // remove the bottom half of the coordinates
    point_cp_mp(ED->EqD->stageData_mp[stage_num].endPts[path_num], endPt.PD_mp.point);
    ED->EqD->stageData_mp[stage_num].endPts[path_num]->size /= 2;
  }

  // find dehom
  getDehomPoint_mp(dehom, ED->EqD->stageData_mp[stage_num].endPts[path_num], ED->EqD->stageData_mp[stage_num].endPts[path_num]->size, &ED->preProcData);

  // print the footer for the point & find the condition number
  ED->EqD->stageData_mp[stage_num].endPt_retVals[path_num] = printStageFooter_mp(ED, stage_num, path_num, &endPt.PD_mp, ED->EqD->stageData_mp[stage_num].endPts[path_num], orig_last, dehom, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, RAWOUT, FAIL, NONSOLN, endPt.retVal, T, trackCount);

  // clear MP
  clear_endgame_data(&endPt)
  clear_point_data_mp(&startPt);
  clear_point_mp(dehom); clear_point_mp(orig_last);

  return;
}

void sortStageEndpoints_mp(int max_threads, int pathMod, point_data_mp PD[], basic_eval_data_mp ED_copy[], basic_eval_data_mp *ED, tracker_config_t T[], FILE **OUT, int stage, double final_tol, int (*ptr_to_eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the witness points found for current stage       *
\***************************************************************/
/*
0  -> UNCLASSIFIED
10 -> MOVE TO NEXT LEVEL
15 -> DO NOT MOVE TO NEXT LEVEL
<0 -> bad path (code = retVal returned from tracking the path)
*/
{
  int i, j, cont, oid, rankDef, finite, indexI = 0, indexJ = 0, num_paths = ED->EqD->stageData_mp[stage].num_paths, num_subs = ED->EqD->num_subsystems;
  int num_sing = 0, num_nonsing = 0, num_bad = 0, num_inf = 0, num_higher_dim = 0;
  mpf_t norm_diff;
  vec_mp tempVec;
  sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_paths * sizeof(sortStruct_mp));

  mpf_init(norm_diff);
  init_vec_mp(tempVec, 0);

  // print header for level
  printf("\nSorting points for stage %d of %d: %d path%s to sort.\n", stage, num_subs, num_paths, num_paths == 1 ? "" : "s");
  fprintf(OUT[0], "\n*****************************************************\n");
  fprintf(OUT[0], "Sorting points for stage %d.\n", stage);
  fprintf(OUT[0], "*****************************************************\n");

  // determine if each of the paths is rank deficient
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, rankDef, finite, oid) schedule(runtime)
#endif
  for (i = 0; i < num_paths; i++)
  { // get the current thread number
    oid = thread_num();

    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sorting %d of %d\n", i, num_paths);

    // only check the ones that were successful
    if (ED_copy[oid].EqD->stageData_mp[stage].endPt_retVals[i] == 0)
    { // setup for evaluation
      ED_copy[oid].EqD->curr_stage_num = stage;
      ED_copy[oid].EqD->curr_path_num = i;

      // setup PD[oid]
      if (ED_copy[oid].EqD->stageData_mp[stage].useIntrinsicSlice)
      { // sort using the intrinsic coordinates
        point_cp_mp(PD[oid].point, ED_copy[oid].EqD->stageData_mp[stage].endPts_in[i]);
      }
      else
      { // sort using original coordinates
        point_cp_mp(PD[oid].point, ED_copy[oid].EqD->stageData_mp[stage].endPts[i]);
      }
      set_zero_mp(PD[oid].time);

      fprintf(OUT[oid], "Path number: %d\n", i);
      // determine if it is rank deficient, finite
      eqbyeqStageSortEndpoint_mp(&rankDef, &finite, ED_copy[oid].EqD->stageData_mp[stage].higherDim[i], &ED_copy[oid], stage, i, &ED_copy[oid].EqD->stageData_mp[stage].condition_nums[i], &T[oid], OUT[oid], &PD[oid], ptr_to_eval);
      
      // copy back 
      if (ED_copy[oid].EqD->stageData_mp[stage].useIntrinsicSlice)
      { // store to the intrinsic coordinates
        point_cp_mp(ED_copy[oid].EqD->stageData_mp[stage].endPts_in[i], PD[oid].point);
      }
      else
      { // store using original coordinates
        point_cp_mp(ED_copy[oid].EqD->stageData_mp[stage].endPts[i], PD[oid].point); 
      }

      // check for success
      if (T[oid].regen_remove_inf && !finite)
      { // dehom point is infinite
        ED_copy[oid].EqD->stageData_mp[stage].endPt_types[i] = retVal_going_to_infinity;
      }
      else if (T[oid].regen_higher_dim_check && ED_copy[oid].EqD->stageData_mp[stage].higherDim[i])
      { // it lies on a higher dimensional component
        ED_copy[oid].EqD->stageData_mp[stage].endPt_types[i] = retVal_higher_dim;
      }
      else if (!rankDef)
      { // classify as non-singular
        ED_copy[oid].EqD->stageData_mp[stage].endPt_types[i] = MOVE_TO_NEXT;
      }
      else
      { // classify as singular
        ED_copy[oid].EqD->stageData_mp[stage].endPt_types[i] = DO_NOT_MOVE_TO_NEXT;
      }
    }
    else
    { // path was not a success - copy over error code
      ED_copy[oid].EqD->stageData_mp[stage].endPt_types[i] = ED_copy[oid].EqD->stageData_mp[stage].endPt_retVals[i];
    }

    // setup sortPts
    mpf_init(sortPts[i].norm);
    sortPts[i].path_num = i;
    if (ED_copy[oid].EqD->stageData_mp[stage].useIntrinsicSlice)
      infNormVec_mp2(sortPts[i].norm, ED_copy[oid].EqD->stageData_mp[stage].endPts_in[i]);
    else 
      infNormVec_mp2(sortPts[i].norm, ED_copy[oid].EqD->stageData_mp[stage].endPts[i]);
  }

  // sort the structure - use qsort to make comparisons efficient
  qsort(sortPts, num_paths, sizeof(sortStruct_mp), sort_order_mp);

  // do the final classificiation - not using OpenMP since we could possibly change the same item at the same time in the while loop
  for (i = 0; i < num_paths; i++)
  {
    indexI = sortPts[i].path_num;
    if (ED->EqD->stageData_mp[stage].endPt_types[indexI] == MOVE_TO_NEXT)
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
        if (cont && ED->EqD->stageData_mp[stage].endPt_retVals[indexJ] == 0)
        { // find difference if the jth path is successful
          if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
          { // subtract using intrinsic coordinates
            vec_sub_mp(tempVec, ED->EqD->stageData_mp[stage].endPts_in[indexI], ED->EqD->stageData_mp[stage].endPts_in[indexJ]);
          }
          else 
          { // subtract using extrinsic coordinates
            vec_sub_mp(tempVec, ED->EqD->stageData_mp[stage].endPts[indexI], ED->EqD->stageData_mp[stage].endPts[indexJ]);
          }

          if (infNormVec_mp(tempVec) < final_tol)
          { // i & j are the same - do not move to next level!
            ED->EqD->stageData_mp[stage].endPt_types[indexI] = ED->EqD->stageData_mp[stage].endPt_types[indexJ] = DO_NOT_MOVE_TO_NEXT;
          }
        }
      } while (cont);
    }

    // add to count
    if (ED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_going_to_infinity || ED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_security_max)
      num_inf++;
    else if (ED->EqD->stageData_mp[stage].endPt_types[indexI] == retVal_higher_dim)
      num_higher_dim++;
    else if (ED->EqD->stageData_mp[stage].endPt_types[indexI] == MOVE_TO_NEXT)
      num_nonsing++;
    else if (ED->EqD->stageData_mp[stage].endPt_types[indexI] == DO_NOT_MOVE_TO_NEXT)
      num_sing++;
    else
      num_bad++;
  }

  // store the counts
  ED->EqD->stageData_mp[stage].num_sing = num_sing;
  ED->EqD->stageData_mp[stage].num_nonsing = num_nonsing;
  ED->EqD->stageData_mp[stage].num_bad = num_bad;
  ED->EqD->stageData_mp[stage].num_inf = num_inf;
  ED->EqD->stageData_mp[stage].num_higher_dim = num_higher_dim;

  // clear the memory
  for (i = num_paths - 1; i >= 0; i--)
  {
    mpf_clear(sortPts[i].norm);
  }
  free(sortPts);

  mpf_clear(norm_diff);
  clear_vec_mp(tempVec);

  return;
}

void eqbyeqStageSortEndpoint_mp(int *rankDef, int *finite, int higherDim, basic_eval_data_mp *ED, int stage_num, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_mp *endPt_mp, int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sorts the endpoint - rankDef, finite & higherDim       *
\***************************************************************/
{
  // initialize
  *rankDef = 0; // we already know higherDim
  *finite = 1;

  // check to see if finite
  if (T->regen_remove_inf)
  { // determine if the endpoint is finite
    *finite = determineEqbyEqFinite(T->finiteThreshold, NULL, endPt_mp, T->Precision, NULL, ED, stage_num, 1);
  }

  // determine if we need to do the rank deficient test (which also refines the endpoint)
  if ((T->regen_higher_dim_check && higherDim == 0) && *finite == 1)
  { // determine if it is rank deficient 
    *rankDef = determineRankDef(condNum, T->final_tolerance, NULL, 52, NULL, endPt_mp, T->Precision, T, OUT, NULL, ED, NULL, ptr_to_eval_mp, NULL);

    // check to see if finite
    if (T->regen_remove_inf)
    { // determine if the endpoint is finite
      *finite = determineEqbyEqFinite(T->finiteThreshold, NULL, endPt_mp, T->Precision, NULL, ED, stage_num, 1);
    }

    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d Rank Def: %d", higherDim, *rankDef);
    }
    else
    {
      fprintf(OUT, "Rank Def: %d", *rankDef);
    }
    if (T->regen_remove_inf)
      fprintf(OUT, " Finite: %d", *finite);
    fprintf(OUT, " CN: %e\n", *condNum);
  }
  else
  { // we know that we are removing this end point
    if (T->regen_higher_dim_check)
    {
      fprintf(OUT, "Higher Dim'l: %d", higherDim);
      if (T->regen_remove_inf)
        fprintf(OUT, " Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
    else
    {
      if (T->regen_remove_inf)
        fprintf(OUT, "Finite: %d\n", *finite);
      else
        fprintf(OUT, "\n");
    }
  }

  return;
}

int printStageFooter_mp(basic_eval_data_mp *ED, int stage_num, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp orig_last, point_mp dehom_mp, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for path stage_num & path_num        *
\***************************************************************/
{
  int i, true_failure = 0, retVal_out = 0, isNumber = 1, isSoln = 0, dehomInf = 0;

  if (retVal_in)
  {
    if (T->regen_remove_inf && infNormVec_mp(dehom_mp) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else if (retVal_in == retVal_sharpening_singular_endpoint)
    { // sharpening determined it was singular
      retVal_out = retVal_sharpening_singular_endpoint;
      true_failure = 1;
    }
    else if (mpf_get_d(endPoint->time->r) < T->minTrackT)
    { // display a warning message if it is a finite path

      if (retVal_in == retVal_reached_minTrackT)
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on stage %d reached the minimum value of T (%e < %e).\n", path_num, stage_num, mpf_get_d(endPoint->time->r), T->minTrackT);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on stage %d reached the minimum value of T (%e < %e).\n", path_num, stage_num, mpf_get_d(endPoint->time->r), T->minTrackT);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }
      else
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on stage %d had final T of %e but had retVal %d.\n", path_num, stage_num, mpf_get_d(endPoint->time->r), retVal_in);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on stage %d had final T of %e but had retVal %d.\n", path_num, stage_num, mpf_get_d(endPoint->time->r), retVal_in);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }

      // consider it a success
      retVal_out = 0;
    }
    else
    { // the path had an error code
      retVal_out = retVal_in;
      true_failure = 1; // this path truly failed
      printResultOfPath(OUT, retVal_in, T);
    }
  }
  else // no error code
  { // check to make sure that dehom is not at infinity
    if (T->regen_remove_inf && infNormVec_mp(dehom_mp) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else
    {
      retVal_out = retVal_in;
      printResultOfPath(OUT, retVal_in, T);
    }
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
      if (!(mpfr_number_p(endPoint->point->coord[i].r) && mpfr_number_p(endPoint->point->coord[i].i)))
        isNumber = 0;

    if (!isNumber)
    {
      retVal_out = retVal_NAN;
      true_failure = 1; // this path truly failed
    }
    else if (stage_num + 1 == ED->EqD->num_subsystems)
    { // verify that the original variables satisfies the original functions
      isSoln = nonsolutions_check_mp(ED->squareSystem.size_f, ED->squareSystem.size_r, orig_vars, orig_last, endPoint->time, T->funcResTol, T->ratioTol, ED->squareSystem.Prog);

      if (!isSoln)
      {
        retVal_out = retVal_Bertini_Junk;
        true_failure = 1; // this path truly failed
      }
    }
  }

  // print the path footer to OUT
  printPathFooterOut_mp(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, first_increase, dehom_mp, T, ED->squareSystem.Prog, ED->preProcData.num_var_gp, 0);

  // determine if we need to print to RAWOUT & FAIL
  if (stage_num + 1 == ED->EqD->num_subsystems)
  { // this is the last stage, so we need to either print to FAIL or RAWOUT
    if (retVal_out && !dehomInf)
    { // path was not successful, so update the number of failures and print the info to FAIL
      if (true_failure)
        trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_mp(FAIL, endPoint, dehom_mp, path_num, retVal_out, isNumber, !isSoln, trackCount, T);

      if (retVal_out == retVal_Bertini_Junk)
      { // print to NONSOLN
        for (i = 0; i < dehom_mp->size; i++)
        {
          print_mp(NONSOLN, 0, &dehom_mp->coord[i]);
          fprintf(NONSOLN, "\n");
        }
        fprintf(NONSOLN, "\n");
      }
    }
    else
    { // path was successful, so update the number of successes and print the info to RAWOUT
      trackCount->successes++;

      if (retVal_in == retVal_reached_minTrackT)  // success but convergence may not be quite accurate enough
        printSuccess_mp(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, first_increase, func_residual, newton_error, t_val_sample, error_sample, -1);
      else
        printSuccess_mp(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, first_increase, func_residual, newton_error, t_val_sample, error_sample, 1);
    }
  }

  return retVal_out;
}

int printWitnessFooter_mp(basic_eval_data_mp *ED, int subsystem_num, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp orig_last, point_mp dehomP, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints the footer for path subsystem_num & path_num    *
\***************************************************************/
{
  int i, true_failure = 0, retVal_out = 0, isNumber = 1, isSoln = 0, dehomInf = 0;

  if (retVal_in)
  {
    if (T->regen_remove_inf && infNormVec_mp(dehomP) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else if (retVal_in == retVal_sharpening_singular_endpoint)
    { // sharpening determined it was singular
      retVal_out = retVal_sharpening_singular_endpoint;
      true_failure = 1;
    }
    else if (mpf_get_d(endPoint->time->r) < T->minTrackT)
    { // display a warning message if it is a finite path

      if (retVal_in == retVal_reached_minTrackT)
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on subsystem %d reached the minimum value of T (%e < %e).\n", path_num, subsystem_num, mpf_get_d(endPoint->time->r), T->minTrackT);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on subsystem %d reached the minimum value of T (%e < %e).\n", path_num, subsystem_num, mpf_get_d(endPoint->time->r), T->minTrackT);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }
      else
      {
        if (T->screenOut)
        {
          printf("WARNING: Path %d on subsystem %d had final T of %e but had retVal %d.\n", path_num, subsystem_num, mpf_get_d(endPoint->time->r), retVal_in);
          printf("         It is not at infinity and will initially be considered a success.\n");
        }
        fprintf(OUT, "WARNING: Path %d on subsystem %d had final T of %e but had retVal %d.\n", path_num, subsystem_num, mpf_get_d(endPoint->time->r), retVal_in);
        fprintf(OUT, "         It is not at infinity and will initially be considered a success.\n");
      }

      // consider it a success
      retVal_out = 0;
    }
    else
    { // the path had an error code
      retVal_out = retVal_in;
      true_failure = 1; // this path truly failed
      printResultOfPath(OUT, retVal_in, T);
    }
  }
  else // no error code
  { // check to make sure that dehom is not at infinity
    if (T->regen_remove_inf && infNormVec_mp(dehomP) > T->finiteThreshold)
    { // de-hom point is at infinity
      dehomInf = 1;
      retVal_out = retVal_going_to_infinity;
      if (T->screenOut)
        printf("De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
      fprintf(OUT, "De-hom point is at infinity (infinity-norm of de-hom point approximation exceeded %e).\n", T->finiteThreshold);
    }
    else
    {
      retVal_out = retVal_in;
      printResultOfPath(OUT, retVal_in, T);
    }
  }

  // if it looks like a successful path, check that output is a number
  if (!retVal_out)
  { // make sure that the output value is a number
    for (i = 0; i < endPoint->point->size && isNumber; i++)
      if (!(mpfr_number_p(endPoint->point->coord[i].r) && mpfr_number_p(endPoint->point->coord[i].i)))
        isNumber = 0;

    if (!isNumber)
    {
      retVal_out = retVal_NAN;
      true_failure = 1; // this path truly failed
    }
    else if (ED->EqD->num_subsystems == 1)
    { // this is the only subsystem, so we need to verify that the original variables satisfies the original function
      isSoln = nonsolutions_check_mp(ED->squareSystem.size_f, ED->squareSystem.size_r, orig_vars, orig_last, endPoint->time, T->funcResTol, T->ratioTol, ED->squareSystem.Prog);
   
      if (!isSoln)
      {
        retVal_out = retVal_Bertini_Junk;
        true_failure = 1; // this path truly failed
      }
    }
  }

  // print the path footer to OUT
  printPathFooterOut_mp(OUT, RAWOUT, 0, path_num, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, first_increase, dehomP, T, ED->squareSystem.Prog, ED->preProcData.num_var_gp, 0);

  // determine if we need to print to RAWOUT & FAIL
  if (ED->EqD->num_subsystems == 1)
  { // this is the only subsystem, so we need to either print to FAIL or RAWOUT
    if (retVal_out && !dehomInf)
    { // path was not successful, so update the number of failures and print the info to FAIL
      if (true_failure)
        trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_mp(FAIL, endPoint, dehomP, path_num, retVal_out, isNumber, !isSoln, trackCount, T);

      if (retVal_out == retVal_Bertini_Junk)
      { // print to NONSOLN
        for (i = 0; i < dehomP->size; i++)
        {
          print_mp(NONSOLN, 0, &dehomP->coord[i]);
          fprintf(NONSOLN, "\n");
        }
        fprintf(NONSOLN, "\n");
      }
    }
    else
    { // path was successful, so update the number of successes and print the info to RAWOUT
      trackCount->successes++;

      if (retVal_in == retVal_reached_minTrackT) // success but convergence may not be quite accurate enough
        printSuccess_mp(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, first_increase, func_residual, newton_error, t_val_sample, error_sample, -1);
      else
        printSuccess_mp(RAWOUT, orig_vars, endPoint->cycle_num, path_num, cond_num, first_increase, func_residual, newton_error, t_val_sample, error_sample, 1);
    }
  }

  return retVal_out;
}

int eqbyeq_witness_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  int level_num;
  eqData_t *EqD = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    basic_eval_data_d *BED = (basic_eval_data_d *)ED_d;
    EqD = BED->EqD;
    level_num = EqD->curr_stage_num;

    // convert to the original coordinates
    intrinsicToExtrinsic_d(out_d, in_d, EqD->witnessData_d[level_num].B, EqD->witnessData_d[level_num].p);
    getDehomPoint_d(out_d, out_d, out_d->size, &BED->preProcData);

    BED = NULL;
  }
  else
  { // compute out_mp
    basic_eval_data_mp *BED = (basic_eval_data_mp *)ED_mp;
    EqD = BED->EqD;
    level_num = EqD->curr_stage_num;

    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);

    // convert to the original coordinates
    intrinsicToExtrinsic_mp(out_mp, in_mp, EqD->witnessData_mp[level_num].B, EqD->witnessData_mp[level_num].p);
    getDehomPoint_mp(out_mp, out_mp, out_mp->size, &BED->preProcData);

    BED = NULL;
  }

  EqD = NULL;

  return 0;
}

int eqbyeq_stage_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  int level_num;
  eqData_t *EqD = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    basic_eval_data_d *BED = (basic_eval_data_d *)ED_d;
    EqD = BED->EqD;
    level_num = EqD->curr_stage_num;

    if (EqD->stageData_d[level_num].useIntrinsicSlice)
    { // convert to extrainsic coordinates
      intrinsicToExtrinsic_d(out_d, in_d, EqD->stageData_d[level_num].B0, EqD->stageData_d[level_num].p0);
      out_d->size /= 2; // remove the bottom half of the coordinates
    }
    else
    { // use the top set of coordinates 
      point_cp_d(out_d, in_d);
      out_d->size /= 2; // remove the bottom half of the coordinates
    }

    // convert to the original coordinates
    getDehomPoint_d(out_d, out_d, out_d->size, &BED->preProcData);

    BED = NULL;
  }
  else
  { // compute out_mp
    basic_eval_data_mp *BED = (basic_eval_data_mp *)ED_mp;
    EqD = BED->EqD;
    level_num = EqD->curr_stage_num;

    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);

    if (EqD->stageData_mp[level_num].useIntrinsicSlice)
    { // convert to extrainsic coordinates
      intrinsicToExtrinsic_mp(out_mp, in_mp, EqD->stageData_mp[level_num].B0, EqD->stageData_mp[level_num].p0);
      out_mp->size /= 2; // remove the bottom half of the coordinates
    }
    else
    { // use the top set of coordinates
      point_cp_mp(out_mp, in_mp);
      out_mp->size /= 2; // remove the bottom half of the coordinates
    }

    // convert to the original coordinates
    getDehomPoint_mp(out_mp, out_mp, out_mp->size, &BED->preProcData);

    BED = NULL;
  }

  EqD = NULL;

  return 0;
}

