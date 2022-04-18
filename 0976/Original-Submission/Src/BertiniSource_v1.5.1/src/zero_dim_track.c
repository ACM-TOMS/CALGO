// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

void zero_dim_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does standard zero dimensional tracking                *
*  in either double precision or adaptive precision             *
\***************************************************************/
{
  int i, oid, eachStartPt, startPointIndex, max = max_threads();
  point_data_d *startPts = NULL;
  endgame_data_t *EG = NULL;
  tracker_config_t *T_copy = NULL;
  basic_eval_data_d *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, *NONSOLN = NULL, **NONSOLN_copy = NULL;

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
  // Find the number of start points
  fscanf(START, "%d", &trackCount->numPoints);
  scanRestOfLine(START);

  // setup NONSOLN 
  if (!ED_d->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // verify positive number
  if (trackCount->numPoints <= 0)
  {
    printf("\n\nERROR: The number of startpoints must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (max > 1 || trackCount->numPoints < 1000) // serial processing can read in the paths at the beginning for small problems
    eachStartPt = 0; // read in all of the start points at the beginning
  else
    eachStartPt = 1; // read in each start point when it is needed

  // setup startPts
  if (eachStartPt)
  { // only need 1 start point - setup on each iteration of the loop
    startPts = (point_data_d *)bmalloc(1 * sizeof(point_data_d));
    init_point_data_d(&startPts[0], 0);
  }
  else
  { // setup the start points
    startPts = (point_data_d *)bmalloc(trackCount->numPoints * sizeof(point_data_d));
    for (i = 0; i < trackCount->numPoints; i++)
    { // setup startPts[i]
      init_point_data_d(&startPts[i], 0);
      setupStart_d(T, &startPts[i], START);
    }
  }

  // setup the rest of the structures
  setup_zero_dim_omp_d(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, T, &BED_copy, ED_d, ED_mp);

  // track each of the start points
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid, startPointIndex) schedule(runtime)
#endif
  for (i = 0; i < trackCount->numPoints; i++)
  { // get current thread number
    oid = thread_num();

    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, trackCount->numPoints);

    if (eachStartPt)
    { // setup the next start point
      startPointIndex = 0;  
      setupStart_d(&T_copy[oid], &startPts[startPointIndex], START);
    }
    else
    { // next start point is setup at index i
      startPointIndex = i;
    }

    // print the header of the path to OUT
    printPathHeader_d(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], i, &BED_copy[oid], eval_func_d);

    // track the path
    zero_dim_track_path_d(i, &EG[oid], &startPts[startPointIndex], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    // check to see if it should be sharpened
    if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, eval_func_d, eval_func_mp, change_prec);
    }

    if (EG[oid].prec < 64)
    { // print footer in double precision
      printPathFooter_d(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], &BED_copy[oid]);
    }
    else
    { // print footer in multi precision
      printPathFooter_mp(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], BED_copy[oid].BED_mp); 
    }
  }

  // clear startPts
  if (eachStartPt)
  { // only need to clear 1 start point
    startPointIndex = 1;
  }
  else
  { // clear all of the start points
    startPointIndex = trackCount->numPoints;
  }
  for (i = startPointIndex - 1; i >= 0; i--)
  { // clear startPts[i]
    clear_point_data_d(&startPts[i]);
  }
  free(startPts);

  // clear the structures
  clear_zero_dim_omp_d(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);

  if (!ED_d->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  return;
}

void zero_dim_track_path_d(int pathNum, endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: actually does the zero-dimensional tracking and sets   *
*  up EG_out                                                    *
\***************************************************************/
{
  EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;
  if (T->MPType == 2)
  { // track using AMP
    EG_out->prec = EG_out->last_approx_prec = 52;

    EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    if (EG_out->prec == 52)
    { // copy over values in double precision
      EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
      EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
      EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
      findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
    }
    else
    { // make sure that the other MP things are set to the correct precision
      mpf_clear(EG_out->function_residual_mp);
      mpf_init2(EG_out->function_residual_mp, EG_out->prec);

      mpf_clear(EG_out->latest_newton_residual_mp);
      mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);

      mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
      mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);

      mpf_clear(EG_out->error_at_latest_sample_point_mp);
      mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);

      // copy over the values
      mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
      mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
      mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
      findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
    }
  }
  else
  { // track using double precision
    EG_out->prec = EG_out->last_approx_prec = 52;
    EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);

    EG_out->first_increase = 0;
    // copy over values in double precision
    EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
    EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
    EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
    findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
  }

  return;
}

void zero_dim_track_path_rank_d(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: actually does the zero-dimensional tracking and sets   *
*  up EG_out                                                    *
\***************************************************************/
{
  EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;
  if (T->MPType == 2)
  { // track using AMP
    EG_out->prec = EG_out->last_approx_prec = 52;

    // track using AMP
    EG_out->retVal = endgame_rank_amp(T->endgameNumber, EG_out->pathNum, &EG_out->condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    if (EG_out->prec == 52)
    { // copy over values in double precision
      EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
      EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
      EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
      findFunctionResidual_d(&EG_out->function_residual_d, &EG_out->PD_d, ED_d, eval_func_d);
    }
    else
    { // make sure that the other MP things are set to the correct precision
      mpf_clear(EG_out->function_residual_mp);
      mpf_init2(EG_out->function_residual_mp, EG_out->prec);

      mpf_clear(EG_out->latest_newton_residual_mp);
      mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);

      mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
      mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);

      mpf_clear(EG_out->error_at_latest_sample_point_mp);
      mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);

      // copy over the values
      mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
      mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
      mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
      findFunctionResidual_mp(EG_out->function_residual_mp, &EG_out->PD_mp, ED_mp, eval_func_mp);
    }
  }
  else
  { // track using double precision
    EG_out->retVal = endgame_rank_d(T->endgameNumber, EG_out->pathNum, &EG_out->condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);

    EG_out->prec = EG_out->last_approx_prec = 52;
    EG_out->first_increase = 0;
    // copy over values in double precision
    EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
    EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
    EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
    findFunctionResidual_d(&EG_out->function_residual_d, &EG_out->PD_d, ED_d, eval_func_d);
  }

  return;
}

void setup_zero_dim_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_d **BED_copy, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup everything needed to do zero dimensional tracking*
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

  // allocate space for EG
  *EG = (endgame_data_t *)bmalloc(max_threads * sizeof(endgame_data_t));
  // initialize 
  for (i = 0; i < max_threads; i++)
    if (T->MPType == 2) 
    { // initialize for AMP tracking
      init_endgame_data(&(*EG)[i], 64);
    }
    else 
    { // initialize for double precision tracking
      init_endgame_data(&(*EG)[i], 52);
    }

  // allocate space to hold pointers to the files
  *OUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *MIDOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *RAWOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *FAIL_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *NONSOLN_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));

  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED_d;
    (*BED_copy)->BED_mp = ED_mp; // make sure that this is pointed to inside of ED_d

    (*OUT_copy)[0] = OUT;
    (*RAWOUT_copy)[0] = RAWOUT;
    (*MIDOUT_copy)[0] = MIDOUT;
    (*FAIL_copy)[0] = FAIL;
    (*NONSOLN_copy)[0] = NONSOLN;
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
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[i]);
      (*trackCount_copy)[i].numPoints = trackCount->numPoints;
    }

    // setup the files
    char *str = NULL;
    int size;
    for (i = 0; i < max_threads; i++)
    {
      size = 1 + snprintf(NULL, 0, "output_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", i);
      (*OUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", i);
      (*MIDOUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "rawout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", i);
      (*RAWOUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "fail_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", i);
      (*FAIL_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", i);
      (*NONSOLN_copy)[i] = fopen(str, "w+");
    }
    free(str);
  }

  return;
}

void clear_zero_dim_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, basic_eval_data_d **BED_copy)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy the relevant data back to the standard spot and   *
*  clear the allocated data that was used by OpenMP             *
\***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int i;

  // clear EG
  for (i = max_threads - 1; i >= 0; i--)
  {
    clear_endgame_data(&(*EG)[i]);
  }
  free(*EG);

  if (max_threads == 1)
  { // set the pointers to NULL since they just pointed to the actual values
    *trackCount_copy = NULL;
    *T_copy = NULL;
    *BED_copy = NULL;

    *OUT_copy[0] = NULL;
    *RAWOUT_copy[0] = NULL;
    *MIDOUT_copy[0] = NULL;
    *FAIL_copy[0] = NULL;

    // free the memory of the file pointers
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
  else if (max_threads > 1)
  {
    // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);

    // clear the copies T, ED_d & ED_mp
    for (i = max_threads - 1; i >= 0; i--)
    { // clear BED_copy - 0 since not using regeneration
      basic_eval_clear_d(&(*BED_copy)[i], 0, (*T_copy)[i].MPType);
       // clear T_copy
      tracker_config_clear(&(*T_copy)[i]);
    }

    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);

    // copy all of the files to the appropriate place
    char ch, *str = NULL;
    int size;
    for (i = 0; i < max_threads; i++)
    {
      // rewind to beginning for reading
      rewind((*OUT_copy)[i]);
      // copy over to OUT
      ch = fgetc((*OUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(OUT, "%c", ch);
        ch = fgetc((*OUT_copy)[i]);
      }
      // close file & delete
      fclose((*OUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "output_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*MIDOUT_copy)[i]);
      // copy over to MIDOUT
      ch = fgetc((*MIDOUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(MIDOUT, "%c", ch);
        ch = fgetc((*MIDOUT_copy)[i]);
      }
      // close file & delete
      fclose((*MIDOUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*RAWOUT_copy)[i]);
      // copy over to RAWOUT
      ch = fgetc((*RAWOUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(RAWOUT, "%c", ch);
        ch = fgetc((*RAWOUT_copy)[i]);
      }
      // close file & delete
      fclose((*RAWOUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "rawout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*FAIL_copy)[i]);
      // copy over to FAIL
      ch = fgetc((*FAIL_copy)[i]);
      while (ch != EOF)
      {
        fprintf(FAIL, "%c", ch);
        ch = fgetc((*FAIL_copy)[i]);
      }
      // close file & delete
      fclose((*FAIL_copy)[i]);
      size = 1 + snprintf(NULL, 0, "fail_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*NONSOLN_copy)[i]);
      // copy over to NONSOLN
      ch = fgetc((*NONSOLN_copy)[i]);
      while (ch != EOF)
      {
        fprintf(NONSOLN, "%c", ch);
        ch = fgetc((*NONSOLN_copy)[i]);
      }
      // close file & delete
      fclose((*NONSOLN_copy)[i]);
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", i);
      remove(str);
    }
    free(str);
    // free file memory
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }

  return;
}

void printPathHeader_d(FILE *OUT, point_data_d *PD, tracker_config_t *T, int pathNum, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header for the path to OUT                  *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int timeDigits = 10, pointDigits = 15;
    eval_struct_d e;
    init_eval_struct_d(e, 0, 0, 0);

    // print the information to OUT
    fprintf(OUT, "\nPath number %d:\nStarting with:\n", pathNum);
    fprintf(OUT, "time: <"); print_d(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_d(OUT, pointDigits, PD->point);

    // evaluate the function at the start point
    eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func_d);

    // print the value of the function to OUT
    fprintf(OUT, "H(Point, time):\n     ");
    printVec_d(OUT, pointDigits, e.funcVals);
    fprintf(OUT, "\n\n");

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nPath number %d:\nStarting with:\n", pathNum);
      printf("time: <"); print_d(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_d(stdout, pointDigits, PD->point);
      printf("H(Point, time):\n     ");
      printVec_d(stdout, pointDigits, e.funcVals);
      printf("\n\n");
    }

    clear_eval_struct_d(e);
  }

  return;
}

void printPathFooter_d(trackingStats *trackCount, endgame_data_t *EG, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines if the path was a success or not and then   *
*prints the footer for the path to OUT & other applicable places*
\***************************************************************/
{
  int i, isNumber = 1, isSoln = 1;
  point_d dehomP;
  init_point_d(dehomP, 0);

  // cast ED as basic_eval_data_d
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;

  printResultOfPath(OUT, EG->retVal, T);

  // make sure that the output value is a number
  for (i = 0; i < T->numVars && isNumber; i++)
    if (isnan(EG->PD_d.point->coord[i].r) || isnan(EG->PD_d.point->coord[i].i) || isinf(EG->PD_d.point->coord[i].r) || isinf(EG->PD_d.point->coord[i].i))
      isNumber = 0;

  // run bertini junk checker if it is a number that looks like a successful path
  if (isNumber)
    if ((EG->PD_d.time->r <= T->minTrackT) || (EG->retVal == 0))
    {
      if (EG->last_approx_prec > 52)
      { // convert to _d
        point_mp_to_d(EG->last_approx_d, EG->last_approx_mp);
      }

      isSoln = nonsolutions_check_d(BED->squareSystem.size_f, BED->squareSystem.size_r, EG->PD_d.point, EG->last_approx_d, EG->PD_d.time, T->funcResTol, T->ratioTol, BED->squareSystem.Prog);
    }

  // find the dehomogenized point
  getDehomPoint_d(dehomP, EG->PD_d.point, EG->PD_d.point->size, &BED->preProcData);

  // if (retVal != 0 && t is too large) || is not a number || is junk, then it is a failure
  if ((EG->retVal != 0 && EG->PD_d.time->r > T->minTrackT) || !isNumber || !isSoln)
  { // update the number of failures
    trackCount->failures++;

    // print the path number, error message, time and point to FAIL
    printFailureMsg_d(FAIL, &EG->PD_d, dehomP, EG->pathNum, EG->retVal, isNumber, !isSoln, trackCount, T);

    // print the footer for OUT
    printPathFooterOut_d(OUT, RAWOUT, 0, EG->pathNum, &EG->PD_d, EG->condition_number, EG->function_residual_d, EG->latest_newton_residual_d, EG->t_val_at_latest_sample_point_d, EG->error_at_latest_sample_point_d, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);

    if (!isSoln)
    { // print to NONSOLN
      for (i = 0; i < dehomP->size; i++)
      {
        print_d(NONSOLN, 0, &dehomP->coord[i]);
        fprintf(NONSOLN, "\n");
      }
      fprintf(NONSOLN, "\n");
    }
  }
  else
  { // update the number of successes
    trackCount->successes++;

    // print the footer for OUT and print to RAWOUT
    if (EG->retVal == 0)
    { // we have convergence
      printPathFooterOut_d(OUT, RAWOUT, 1, EG->pathNum, &EG->PD_d, EG->condition_number, EG->function_residual_d, EG->latest_newton_residual_d, EG->t_val_at_latest_sample_point_d, EG->error_at_latest_sample_point_d, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else if (EG->retVal == retVal_sharpening_singular_endpoint)
    { // sharpening failure due to be singular
      printPathFooterOut_d(OUT, RAWOUT, retVal_sharpening_singular_endpoint, EG->pathNum, &EG->PD_d, EG->condition_number, EG->function_residual_d, EG->latest_newton_residual_d, EG->t_val_at_latest_sample_point_d, EG->error_at_latest_sample_point_d, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else if (EG->retVal == retVal_sharpening_failed)
    { // only failure when doing sharpening
      printPathFooterOut_d(OUT, RAWOUT, retVal_sharpening_failed, EG->pathNum, &EG->PD_d, EG->condition_number, EG->function_residual_d, EG->latest_newton_residual_d, EG->t_val_at_latest_sample_point_d, EG->error_at_latest_sample_point_d, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else
    { // some other convergence failure
      printPathFooterOut_d(OUT, RAWOUT, -1, EG->pathNum, &EG->PD_d, EG->condition_number, EG->function_residual_d, EG->latest_newton_residual_d, EG->t_val_at_latest_sample_point_d, EG->error_at_latest_sample_point_d, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
  }

  clear_point_d(dehomP);

  return;
}

void printPathFooterOut_d(FILE *OUT, FILE *RAWOUT, int success, int pathNum, point_data_d *PD, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, point_d dehomP, tracker_config_t *T, prog_t *Prog, int print_dehom, int eval_orig_f)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the footer to OUT                               *
* There are 5 levels of 'success':                              *
*   singular_endpoint - sharpening module found it was singular *
*   sharpening_failed - failure only in the sharpening module   *
*   -1 - no convergence before minTrackT                        *
*    0 - failure                                                *
*    1 - convergence                                            *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int timeDigits = 10, pointDigits = 15;

    // print the information to OUT
    fprintf(OUT, "\nEnding with\n");
    fprintf(OUT, "time: <"); print_d(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_d(OUT, pointDigits, PD->point);

    // print the dehomogenized point, if needed
    if (print_dehom)  
    {
      fprintf(OUT, "De-hom point: ");
      printVec_d(OUT, pointDigits, dehomP);
    }

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nEnding with:\n");
      printf("time: <"); print_d(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_d(stdout, pointDigits, PD->point);

      // print the dehomogenized point, if needed
      if (print_dehom)
      {
        printf("De-hom point: ");
        printVec_d(stdout, pointDigits, dehomP);
      }
    }

    if (eval_orig_f)
    { // evaluate the original function
      eval_struct_d e;
      init_eval_struct_d(e, 0, 0, 0);

      evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, Prog);
      fprintf(OUT, "Orig f(Point, time):\n     ");
      printVec_d(OUT, pointDigits, e.funcVals);

      if (T->screenOut)
      { // print all the same info to the screen
        printf("Orig f(Point, time):\n     ");
        printVec_d(stdout, pointDigits, e.funcVals);
      }

      clear_eval_struct_d(e);
    }

    // print the value of the function to OUT
    fprintf(OUT, "||H(Point, time)||: %.15e\n", func_residual);
    fprintf(OUT, "\n_____________________________________________________\n");

    if (T->screenOut)
    {
      printf("||H(Point, time)||: %.15e\n", func_residual);
      printf("\n_____________________________________________________\n");
    }
  }

  if (success)
  { // print the successful path to RAWOUT
    printSuccess_d(RAWOUT, PD->point, PD->cycle_num, pathNum, cond_num, func_residual, newton_error, t_val_sample, error_sample, success);
  }

  return;
}

void printSuccess_d(FILE *RAWOUT, point_d orig_vars, int cycle_num, int path_num, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, int success_flag)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the info to RAWOUT                              *
\***************************************************************/
{
  int i;

  fprintf(RAWOUT, "%d\n%d\n", path_num, 52); // 52 == prec for double

  for (i = 0; i < orig_vars->size; i++)
    fprintf(RAWOUT, "%.15e %.15e\n", orig_vars->coord[i].r, orig_vars->coord[i].i);

  fprintf(RAWOUT, "%.15e\n%.15e\n%.15e\n%.15e\n%.15e\n", func_residual, cond_num, newton_error, t_val_sample, error_sample);
  fprintf(RAWOUT, "%.15e\n%d\n%d\n", 0.0, cycle_num, success_flag); // since no increase in precision, its first increase did not occur

  return;
}

void printFailureMsg_d(FILE *FAIL, point_data_d *PD, point_d dehomP, int pathNum, int retVal, int isNumber, int isJunk, trackingStats *trackCount, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints error message to FAIL and updates trackCount    *
\***************************************************************/
{
  int j;

  fprintf(FAIL, "-------------------------\nPath Number %d\n", pathNum);

  if (isJunk)
  {
    fprintf(FAIL, "The point does not satisfy the original system.\n");
    trackCount->junkCount++;
  }
  else if (!isNumber)
  {
    fprintf(FAIL, "The point has atleast one coordinate that is not a valid number.\n");
    trackCount->nanCount++;
  }
  else if (retVal == retVal_going_to_infinity)
  {
    fprintf(FAIL, "Path went to infinity (inf-norm of path approximation exceeded %e).\n", T->goingToInfinity);
    trackCount->infCount++;
  }
  else if (retVal == retVal_step_size_too_small)
  {
    fprintf(FAIL, "Step size dropped below min.\n");
    trackCount->sizeCount++;
  }
  else if (retVal == retVal_PSEG_failed)
  {
    fprintf(FAIL, "Error when tracking in power series endgame.\n");
    trackCount->PSEGCount++;
  }
  else if (retVal == retVal_max_prec_reached)
  {
    fprintf(FAIL, "The maximum precision (%d) was reached.\n", T->AMP_max_prec);
    trackCount->precCount++;
  }
  else if (retVal == retVal_cycle_num_too_high)
  {
    fprintf(FAIL, "The cycle num was too high\n");
    trackCount->cycleCount++;
  }
  else if (retVal == retVal_too_many_steps)
  {
    fprintf(FAIL, "The maximum number of steps was reached\n");
    trackCount->stepCount++;
  }
  else if (retVal == retVal_refining_failed)
  {
    fprintf(FAIL, "Refining failed\n");
    trackCount->refineCount++;
  }
  else if (retVal == retVal_sharpening_singular_endpoint)
  {
    fprintf(FAIL, "Sharpening failed - singular endpoint\n");
    trackCount->refineCount++;
  }
  else if (retVal == retVal_sharpening_failed)
  {
    fprintf(FAIL, "Sharpening failed\n");
    trackCount->refineCount++;
  }
  else if (retVal == -1)
  {
    fprintf(FAIL, "Linear solving has failed\n");
    trackCount->otherCount++;
  }
  else if (retVal == retVal_security_max)
  {
    fprintf(FAIL, "Path was truncated for having 2 consecutive endpoint approximations exceed %e.\n", T->securityMaxNorm);
    trackCount->securityCount++;
  }
  else
  {
    fprintf(FAIL, "retVal: %d\n", retVal);
    trackCount->otherCount++;
  }

  fprintf(FAIL, "time: %.15e\n", PD->time->r);
  for (j = 0; j < dehomP->size; j++)
  {
    print_d(FAIL, 16, &dehomP->coord[j]);
    fprintf(FAIL, "\n");
  }

  return;
}

void printResultOfPath(FILE *OUT, int retVal, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints result of the path to OUT                       *
\***************************************************************/
{
  if (retVal == retVal_going_to_infinity)
  {
    if (T->screenOut)
      printf("Path went to infinity (inf-norm of path approximation exceeded %e).\n", T->goingToInfinity);
    fprintf(OUT, "Path went to infinity (inf-norm of path approximation exceeded %e).\n", T->goingToInfinity);
  }
  else if (retVal == retVal_step_size_too_small)
  {
    if (T->screenOut)
      printf("Step size dropped below min.\n");
    fprintf(OUT, "Step size dropped below min.\n");
  }
  else if (retVal == retVal_PSEG_failed)
  {
    if (T->screenOut)
      printf("Error when tracking in power series endgame.\n");
    fprintf(OUT, "Error when tracking in power series endgame.\n");
  }
  else if (retVal == retVal_EG_failed_to_converge)
  {
    if (T->screenOut)
      printf("The endgame failed to converge before reaching NBHDRADIUS.\n");
    fprintf(OUT, "The endgame failed to converge before reaching NBHDRADIUS.\n");
  }
  else if (retVal == retVal_max_prec_reached)
  {
    if (T->screenOut)
      printf("The maximum precision (%d) was reached.\n", T->AMP_max_prec);
    fprintf(OUT, "The maximum precision (%d) was reached.\n", T->AMP_max_prec);
  }
  else if (retVal == retVal_cycle_num_too_high)
  {
    if (T->screenOut)
      printf("The cycle num was too high\n");
    fprintf(OUT, "The cycle num was too high\n");
  }
  else if (retVal == retVal_too_many_steps)
  {
    if (T->screenOut)
      printf("The maximum number of steps was reached\n");
    fprintf(OUT, "The maximum number of steps was reached\n");
  }
  else if (retVal == retVal_refining_failed)
  {
    if (T->screenOut)
      printf("Refining failed\n");
    fprintf(OUT, "Refining failed\n");
  }
  else if (retVal == retVal_Bertini_Junk)
  {
    if (T->screenOut)
      printf("The point does not satisfy the original system.\n");
    fprintf(OUT, "The point does not satisfy the original system.\n");
  }
  if (retVal == retVal_security_max)
  {
    if (T->screenOut)
      printf("Path was truncated for approaching infinity (inf-norm of path approximation exceeded %e).\n", T->securityMaxNorm);
    fprintf(OUT, "Path was truncated for approaching infinity (inf-norm of path approximation exceeded %e).\n", T->securityMaxNorm);
  }

  return;
}

void printBasicFooter_d(FILE *OUT, point_data_d *PD, tracker_config_t *T, double func_residual)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the basic footer to OUT                         *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int timeDigits = 10, pointDigits = 15;

    // print the information to OUT
    fprintf(OUT, "\nEnding with\n");
    fprintf(OUT, "time: <"); print_d(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_d(OUT, pointDigits, PD->point);

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nEnding with:\n");
      printf("time: <"); print_d(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_d(stdout, pointDigits, PD->point);
    }

    // print the value of the function to OUT
    fprintf(OUT, "||H(Point, time)||: %.15e\n", func_residual);
    fprintf(OUT, "\n_____________________________________________________\n");

    if (T->screenOut)
    {
      printf("||H(Point, time)||: %.15e\n", func_residual);
      printf("\n_____________________________________________________\n");
    }
  }

  return;
}

/////////// MP VERSIONS ///////////////////

void zero_dim_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does standard zero dimensional tracking                *
\***************************************************************/
{
  int i, oid, eachStartPt, startPointIndex, max = max_threads();
  point_data_mp *startPts = NULL;
  endgame_data_t *EG = NULL;
  tracker_config_t *T_copy = NULL;
  basic_eval_data_mp *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, *NONSOLN = NULL, **NONSOLN_copy = NULL;

  if (max == 1)
    eachStartPt = 1; // read in each start point when it is needed
  else
    eachStartPt = 0; // read in all of the start points at the beginning

  // top of RAWOUT - number of varialbes and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
  // Find the number of start points
  fscanf(START, "%d", &trackCount->numPoints);
  scanRestOfLine(START);

  if (!ED->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // verify positive number
  if (trackCount->numPoints <= 0)
  {
    printf("\n\nERROR: The number of startpoints must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup startPts
  if (eachStartPt)
  { // only need 1 start point - setup on each iteration of the loop
    startPts = (point_data_mp *)bmalloc(1 * sizeof(point_data_mp));
    init_point_data_mp2(&startPts[0], 0, T->Precision);
  }
  else
  { // setup the start points
    startPts = (point_data_mp *)bmalloc(trackCount->numPoints * sizeof(point_data_mp));
    for (i = 0; i < trackCount->numPoints; i++)
    { // setup startPts[i]
      init_point_data_mp2(&startPts[i], 0, T->Precision);
      setupStart_mp(T, &startPts[i], START);     
    }
  }

  // setup the rest of the structures
  setup_zero_dim_omp_mp(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, T, &BED_copy, ED); 

  // track each of the start points
#ifdef _OPENMP
  #pragma omp parallel for private(i, oid, startPointIndex) schedule(runtime)
#endif
  for (i = 0; i < trackCount->numPoints; i++)
  { // get current thread number
    oid = thread_num();

    // print the path number if needed
    if (pathMod > 0)
      if (!(i % pathMod))
        printf("Tracking path %d of %d\n", i, trackCount->numPoints);

    if (eachStartPt)
    { // setup the next start point
      startPointIndex = 0;
      setupStart_mp(&T_copy[oid], &startPts[startPointIndex], START);
    }
    else
    { // next start point is setup at index i
      startPointIndex = i;
    }

    // print the header of the path to OUT
    printPathHeader_mp(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], i, &BED_copy[oid], eval_func);

    // track the path
    zero_dim_track_path_mp(i, &EG[oid], &startPts[startPointIndex], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], eval_func, find_dehom);

    // check to see if it should be sharpened
    if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], NULL, &BED_copy[oid], NULL, eval_func, NULL);
    }

    // print the footer of the path to OUT and print to where else applicable
    printPathFooter_mp(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], &BED_copy[oid]);
  }

  // clear startPts
  if (eachStartPt)
  { // only need to clear 1 start point
    startPointIndex = 1;
  }
  else
  { // clear all of the start points
    startPointIndex = trackCount->numPoints;
  }
  for (i = startPointIndex - 1; i >= 0; i--)
  { // clear startPts[i]
    clear_point_data_mp(&startPts[i]);
  }
  free(startPts);

  // clear the structures
  clear_zero_dim_omp_mp(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);

  if (!ED->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  return;
}

void zero_dim_track_path_mp(int pathNum, endgame_data_t *EG_out, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: actually does the zero-dimensional tracking and sets   *
*  up EG_out                                                    *
\***************************************************************/
{
  EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;

  // track using MP
  EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED, eval_func_mp, find_dehom);

  EG_out->prec = EG_out->last_approx_prec = T->Precision;
  EG_out->first_increase = 0;

  // copy over the values
  mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
  mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
  mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
  findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED, eval_func_mp);

  return;
}

void zero_dim_track_path_rank_mp(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, endgame_data_t *EG_out, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: actually does the zero-dimensional tracking and sets   *
*  up EG_out                                                    *
\***************************************************************/
{
  EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;

  // track using MP
  EG_out->retVal = endgame_rank_mp(T->endgameNumber, EG_out->pathNum, &EG_out->condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, &EG_out->PD_mp, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED, eval_func_mp, find_dehom);

  EG_out->prec = EG_out->last_approx_prec = T->Precision;
  EG_out->first_increase = 0;

  // copy over the values
  mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
  mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
  mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
  findFunctionResidual_mp(EG_out->function_residual_mp, &EG_out->PD_mp, ED, eval_func_mp);

  return;
}

void setup_zero_dim_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_mp **BED_copy, basic_eval_data_mp *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup everything needed to do zero dimensional tracking*
*  using OpenMP                                                 *
\***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int i;

  // error checking
  if (max_threads <= 0)
  {
    printf("\n\nERROR: The number of threads needs to be positive when setting up for tracking!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // allocate space for EG and initialize to the fixed precision
  *EG = (endgame_data_t *)bmalloc(max_threads * sizeof(endgame_data_t));
  for (i = 0; i < max_threads; i++)
  {
    init_endgame_data(&(*EG)[i], T->Precision);
  }

  // allocate space to hold pointers to the files
  *OUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *MIDOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *RAWOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *FAIL_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *NONSOLN_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));

  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED;

    (*OUT_copy)[0] = OUT;
    (*RAWOUT_copy)[0] = RAWOUT;
    (*MIDOUT_copy)[0] = MIDOUT;
    (*FAIL_copy)[0] = FAIL;
    (*NONSOLN_copy)[0] = NONSOLN;
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
      // copy ED_d & ED_mp
      cp_basic_eval_data_mp(&(*BED_copy)[i], ED);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[i]);
      (*trackCount_copy)[i].numPoints = trackCount->numPoints;
    }

    // setup the files
    char *str = NULL;
    int size;
    for (i = 0; i < max_threads; i++)
    {
      size = 1 + snprintf(NULL, 0, "output_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", i);
      (*OUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", i);
      (*MIDOUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "rawout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", i);
      (*RAWOUT_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "fail_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", i);
      (*FAIL_copy)[i] = fopen(str, "w+");

      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", i);
      (*NONSOLN_copy)[i] = fopen(str, "w+");
    }
    free(str);
  }

  return;
}

void clear_zero_dim_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, basic_eval_data_mp **BED_copy)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy the relevant data back to the standard spot and   *
*  clear the allocated data that was used by OpenMP             *
\***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int i;

  // clear EG
  for (i = max_threads - 1; i >= 0; i--)
  {
    clear_endgame_data(&(*EG)[i]);
  }
  free(*EG);

  if (max_threads == 1)
  { // set the pointers to NULL since they just pointed to the actual values
    *trackCount_copy = NULL;
    *T_copy = NULL;
    *BED_copy = NULL;

    *OUT_copy[0] = NULL;
    *RAWOUT_copy[0] = NULL;
    *MIDOUT_copy[0] = NULL;
    *FAIL_copy[0] = NULL;
    *NONSOLN_copy[0] = NULL;

    // free the memory of the file pointers
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
  else if (max_threads > 1)
  {
    // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);

    // clear the copies T & ED
    for (i = max_threads - 1; i >= 0; i--)
    { // clear BED_copy - 0 since not using regeneration, 1 - for clearing Prog
      basic_eval_clear_mp(&(*BED_copy)[i], 0, 1);
       // clear T_copy
      tracker_config_clear(&(*T_copy)[i]);
    }

    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);

    // copy all of the files to the appropriate place
    char ch, *str = NULL;
    int size;
    for (i = 0; i < max_threads; i++)
    {
      // rewind to beginning for reading
      rewind((*OUT_copy)[i]);
      // copy over to OUT
      ch = fgetc((*OUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(OUT, "%c", ch);
        ch = fgetc((*OUT_copy)[i]);
      }
      // close file & delete
      fclose((*OUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "output_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*MIDOUT_copy)[i]);
      // copy over to MIDOUT
      ch = fgetc((*MIDOUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(MIDOUT, "%c", ch);
        ch = fgetc((*MIDOUT_copy)[i]);
      }
      // close file & delete
      fclose((*MIDOUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "midout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*RAWOUT_copy)[i]);
      // copy over to RAWOUT
      ch = fgetc((*RAWOUT_copy)[i]);
      while (ch != EOF)
      {
        fprintf(RAWOUT, "%c", ch);
        ch = fgetc((*RAWOUT_copy)[i]);
      }
      // close file & delete
      fclose((*RAWOUT_copy)[i]);
      size = 1 + snprintf(NULL, 0, "rawout_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*FAIL_copy)[i]);
      // copy over to FAIL
      ch = fgetc((*FAIL_copy)[i]);
      while (ch != EOF)
      {
        fprintf(FAIL, "%c", ch);
        ch = fgetc((*FAIL_copy)[i]);
      }
      // close file & delete
      fclose((*FAIL_copy)[i]);
      size = 1 + snprintf(NULL, 0, "fail_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", i);
      remove(str);

      // rewind to beginning for reading
      rewind((*NONSOLN_copy)[i]);
      // copy over to FAIL
      ch = fgetc((*NONSOLN_copy)[i]);
      while (ch != EOF)
      {
        fprintf(NONSOLN, "%c", ch);
        ch = fgetc((*NONSOLN_copy)[i]);
      }
      // close file & delete
      fclose((*NONSOLN_copy)[i]);
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", i);
      remove(str);
    }
    free(str);
    // free file memory
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }

  return;
}

void printPathHeader_mp(FILE *OUT, point_data_mp *PD, tracker_config_t *T, int pathNum, void const *ED, int (*eval_f)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header for the path to OUT                  *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int timeDigits = 10, pointDigits = 15;
    eval_struct_mp e;
    init_eval_struct_mp(e, 0, 0, 0);

    // print the information to OUT
    fprintf(OUT, "\nPath number %d:\nStarting with:\n", pathNum);
    fprintf(OUT, "time: <"); print_mp(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_mp(OUT, pointDigits, PD->point);

    // evaluate the function at the start point
    eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_f);

    // print the value of the function to OUT
    fprintf(OUT, "H(Point, time):\n     ");
    printVec_mp(OUT, pointDigits, e.funcVals);
    fprintf(OUT, "\n\n");

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nPath number %d:\nStarting with:\n", pathNum);
      printf("time: <"); print_mp(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_mp(stdout, pointDigits, PD->point);
      printf("H(Point, time):\n     ");
      printVec_mp(stdout, pointDigits, e.funcVals);
      printf("\n\n");
    }

    // clear
    clear_eval_struct_mp(e);
  }

  return;
}

void printPathFooter_mp(trackingStats *trackCount, endgame_data_t *EG, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines if the path was a success or not and then   *
*prints the footer for the path to OUT & other applicable places*
\***************************************************************/
{
  int i, isNumber = 1, isSoln = 1;
  point_mp dehomP;
  comp_d tmpDouble;

  init_point_mp(dehomP, 0);

  // cast ED as basic_eval_data_d
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;

  printResultOfPath(OUT, EG->retVal, T);

  mp_to_d(tmpDouble, EG->PD_mp.time);

  // make sure that the output value is a number
  for (i = 0; i < T->numVars && isNumber; i++)
    if (!(mpfr_number_p(EG->PD_mp.point->coord[i].r) && mpfr_number_p(EG->PD_mp.point->coord[i].i)))
      isNumber = 0;

  // run bertini junk checker if it is a number that looks like a successful path
  if (isNumber)
    if ((tmpDouble->r <= T->minTrackT) || (EG->retVal == 0))
    {
      if (EG->last_approx_prec < 64)
      { // copy to _mp
        point_d_to_mp(EG->last_approx_mp, EG->last_approx_d);
      }

      isSoln = nonsolutions_check_mp(BED->squareSystem.size_f, BED->squareSystem.size_r, EG->PD_mp.point, EG->last_approx_mp, EG->PD_mp.time, T->funcResTol, T->ratioTol, BED->squareSystem.Prog);
    }

  // find the dehomogenized point
  getDehomPoint_mp(dehomP, EG->PD_mp.point, EG->PD_mp.point->size, &BED->preProcData);

  // if (retVal != 0 && t is too large) || is not a number || is junk, then it is a failure
  if ((EG->retVal != 0 && tmpDouble->r > T->minTrackT) || !isNumber || !isSoln)
  { // update the number of failures
    trackCount->failures++;

    // print the path number, error message, time and point to FAIL
    printFailureMsg_mp(FAIL, &EG->PD_mp, dehomP, EG->pathNum, EG->retVal, isNumber, !isSoln, trackCount, T);

    // print the footer for OUT
    printPathFooterOut_mp(OUT, RAWOUT, 0, EG->pathNum, &EG->PD_mp, EG->condition_number, EG->function_residual_mp, EG->latest_newton_residual_mp, EG->t_val_at_latest_sample_point_mp, EG->error_at_latest_sample_point_mp, EG->first_increase, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);

    if (!isSoln)
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
  { // update the number of successes
    trackCount->successes++;

    // print the footer for OUT and print to RAWOUT
    if (EG->retVal == 0)
    { // convergence was good
      printPathFooterOut_mp(OUT, RAWOUT, 1, EG->pathNum, &EG->PD_mp, EG->condition_number, EG->function_residual_mp, EG->latest_newton_residual_mp, EG->t_val_at_latest_sample_point_mp, EG->error_at_latest_sample_point_mp, EG->first_increase, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else if (EG->retVal == retVal_sharpening_singular_endpoint)
    { // only failure when doing sharpening
      printPathFooterOut_mp(OUT, RAWOUT, retVal_sharpening_singular_endpoint, EG->pathNum, &EG->PD_mp, EG->condition_number, EG->function_residual_mp, EG->latest_newton_residual_mp, EG->t_val_at_latest_sample_point_mp, EG->error_at_latest_sample_point_mp, EG->first_increase, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else if (EG->retVal == retVal_sharpening_failed)
    { // only failure when doing sharpening
      printPathFooterOut_mp(OUT, RAWOUT, retVal_sharpening_failed, EG->pathNum, &EG->PD_mp, EG->condition_number, EG->function_residual_mp, EG->latest_newton_residual_mp, EG->t_val_at_latest_sample_point_mp, EG->error_at_latest_sample_point_mp, EG->first_increase, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
    else
    { // some other convergence failure
      printPathFooterOut_mp(OUT, RAWOUT, -1, EG->pathNum, &EG->PD_mp, EG->condition_number, EG->function_residual_mp, EG->latest_newton_residual_mp, EG->t_val_at_latest_sample_point_mp, EG->error_at_latest_sample_point_mp, EG->first_increase, dehomP, T, BED->squareSystem.Prog, BED->preProcData.num_var_gp, 1);
    }
  }

  // clear
  clear_point_mp(dehomP);

  return;
}

void printPathFooterOut_mp(FILE *OUT, FILE *RAWOUT, int success, int pathNum, point_data_mp *PD, double cond_num, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, double first_increase, point_mp dehomP, tracker_config_t *T, prog_t *Prog, int print_dehom, int eval_orig_f)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the footer to OUT                               *
*  5 levels of success                                          *
*   singular_endpoint - sharpening module found it was singular *
*   sharpening_failed - failure only in the sharpening module   *
*   -1 - no convergence before minTrackT                        *
*    0 - failure                                                *
*    1 - convergence                                            *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int base = 10, timeDigits = 10, pointDigits = 15;

    // print the information to OUT
    fprintf(OUT, "\nEnding with\n");
    fprintf(OUT, "time: <"); print_mp(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_mp(OUT, pointDigits, PD->point);

    // print the dehomogenized point, if needed
    if (print_dehom)
    {
      fprintf(OUT, "De-hom point: ");
      printVec_mp(OUT, pointDigits, dehomP);
    }

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nEnding with:\n");
      printf("time: <"); print_mp(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_mp(stdout, pointDigits, PD->point);

      // print the dehomogenized point, if we added variables
      if (print_dehom)
      {
        printf("De-hom point: ");
        printVec_mp(stdout, pointDigits, dehomP);
      }
    }

    if (eval_orig_f)
    { // evaluate the original function
      eval_struct_mp e;
      // initialize
      init_eval_struct_mp(e, 0, 0, 0);

      evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, Prog);
      fprintf(OUT, "Orig f(Point, time):\n     ");
      printVec_mp(OUT, pointDigits, e.funcVals);

      if (T->screenOut)
      {
        printf("Orig f(Point, time):\n     ");
        printVec_mp(stdout, pointDigits, e.funcVals);
      }
      // clear
      clear_eval_struct_mp(e);
    }

    // print the value of the function to OUT
    fprintf(OUT, "||H(Point, time)||: ");
    mpf_out_str(OUT, base, pointDigits, func_residual);
    fprintf(OUT, "\n\n_____________________________________________________\n");

    if (T->screenOut)
    {
      printf("||H(Point, time)||: ");
      mpf_out_str(stdout, base, pointDigits, func_residual);
      printf("\n\n_____________________________________________________\n");
    }
  }

  if (success)
  { // print to RAWOUT
    printSuccess_mp(RAWOUT, PD->point, PD->cycle_num, pathNum, cond_num, first_increase, func_residual, newton_error, t_val_sample, error_sample, success);
  }

  return;
}

void printSuccess_mp(FILE *RAWOUT, point_mp orig_vars, int cycle_num, int path_num, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, int success_flag)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the info to RAWOUT                              *
\***************************************************************/
{
  int i;
  long e1, e2;
  char *ch1, *ch2;

  fprintf(RAWOUT, "%d\n%d\n", path_num, (int) mpf_get_prec(orig_vars->coord[0].r));

  for (i = 0; i < orig_vars->size; i++)
  {
    ch1 = mpf_get_str(NULL, &e1, 10, 0, orig_vars->coord[i].r);
    ch2 = mpf_get_str(NULL, &e2, 10, 0, orig_vars->coord[i].i);

    if (ch1[0] != '-')
      if (ch2[0] != '-')
        fprintf(RAWOUT, "0.%se%ld 0.%se%ld\n", ch1, e1, ch2, e2);
      else
        fprintf(RAWOUT, "0.%se%ld -0.%se%ld\n", ch1, e1, &ch2[1], e2);
    else
      if (ch2[0] != '-')
        fprintf(RAWOUT, "-0.%se%ld 0.%se%ld\n", &ch1[1], e1, ch2, e2);
      else
        fprintf(RAWOUT, "-0.%se%ld -0.%se%ld\n", &ch1[1], e1, &ch2[1], e2);

    mpfr_free_str(ch1);
    mpfr_free_str(ch2);
  }

  // print the data to RAWOUT
  mpf_out_str(RAWOUT, 10, 15, func_residual);
  fprintf(RAWOUT, "\n%.15e\n", cond_num);
  mpf_out_str(RAWOUT, 10, 15, newton_error);
  fprintf(RAWOUT, "\n");
  mpf_out_str(RAWOUT, 10, 15, t_val_sample);
  fprintf(RAWOUT, "\n");
  mpf_out_str(RAWOUT, 10, 15, error_sample);
  fprintf(RAWOUT, "\n%.15e\n%d\n%d\n", first_increase, cycle_num, success_flag);

  return;
}

void printFailureMsg_mp(FILE *FAIL, point_data_mp *PD, point_mp dehomP, int pathNum, int retVal, int isNumber, int isJunk, trackingStats *trackCount, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints error message to FAIL and updates trackCount    *
\***************************************************************/
{
  int j;

  fprintf(FAIL, "-------------------------\nPath Number %d\n", pathNum);

  if (isJunk)
  {
    fprintf(FAIL, "The point does not satisfy the original system.\n");
    trackCount->junkCount++;
  }
  else if (!isNumber)
  {
    fprintf(FAIL, "The point has atleast one coordinate that is not a valid number.\n");
    trackCount->nanCount++;
  }
  else if (retVal == retVal_going_to_infinity)
  {
    fprintf(FAIL, "Path went to infinity (inf-norm of path approximation exceeded %e).\n", T->goingToInfinity);
    trackCount->infCount++;
  }
  else if (retVal == retVal_step_size_too_small)
  {
    fprintf(FAIL, "Step size dropped below min.\n");
    trackCount->sizeCount++;
  }
  else if (retVal == retVal_PSEG_failed)
  {
    fprintf(FAIL, "Error when tracking in power series endgame.\n");
    trackCount->PSEGCount++;
  }
  else if (retVal == retVal_max_prec_reached)
  {
    fprintf(FAIL, "The maximum precision (%d) was reached.\n", T->AMP_max_prec);
    trackCount->precCount++;
  }
  else if (retVal == retVal_cycle_num_too_high)
  {
    fprintf(FAIL, "The cycle num was too high\n");
    trackCount->cycleCount++;
  }
  else if (retVal == retVal_too_many_steps)
  {
    fprintf(FAIL, "The maximum number of steps was reached\n");
    trackCount->stepCount++;
  }
  else if (retVal == retVal_refining_failed)
  {
    fprintf(FAIL, "Refining failed\n");
    trackCount->refineCount++;
  }
  else if (retVal == retVal_sharpening_singular_endpoint)
  {
    fprintf(FAIL, "Sharpening failed - singular endpoint\n");
    trackCount->refineCount++;
  }
  else if (retVal == retVal_sharpening_failed)
  {
    fprintf(FAIL, "Sharpening failed\n");
    trackCount->refineCount++;
  }
  else if (retVal == retVal_security_max)
  {
    fprintf(FAIL, "Path was truncated for having 2 consecutive endpoint approximations exceed %e.\n", T->securityMaxNorm);
    trackCount->securityCount++;
  }
  else
  {
    fprintf(FAIL, "retVal: %d\n", retVal);
    trackCount->otherCount++;
  }

  fprintf(FAIL, "time: %.15e\n", mpf_get_d(PD->time->r));
  for (j = 0; j < dehomP->size; j++)
  {
    print_mp(FAIL, 0, &dehomP->coord[j]);
    fprintf(FAIL, "\n");
  }

  return;
}

void printBasicFooter_mp(FILE *OUT, point_data_mp *PD, tracker_config_t *T, mpf_t func_residual)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the basic footer to OUT                         *
\***************************************************************/
{
  if (T->outputLevel >= 0)
  {
    int base = 10, timeDigits = 10, pointDigits = 15;

    // print the information to OUT
    fprintf(OUT, "\nEnding with\n");
    fprintf(OUT, "time: <"); print_mp(OUT, timeDigits, PD->time); fprintf(OUT, ">\n");
    fprintf(OUT, "Point: ");
    printVec_mp(OUT, pointDigits, PD->point);

    if (T->screenOut)
    { // print all the same info to the screen
      printf("\nEnding with:\n");
      printf("time: <"); print_mp(stdout, timeDigits, PD->time); printf(">\n");
      printf("Point: ");
      printVec_mp(stdout, pointDigits, PD->point);
    }

    // print the value of the function to OUT
    fprintf(OUT, "||H(Point, time)||: ");
    mpf_out_str(OUT, base, pointDigits, func_residual);
    fprintf(OUT, "\n\n_____________________________________________________\n");

    if (T->screenOut)
    {
      printf("||H(Point, time)||: ");
      mpf_out_str(stdout, base, pointDigits, func_residual);
      printf("\n\n_____________________________________________________\n");
    }
  }

  return;
}

int zero_dim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the dehom point                                *
\***************************************************************/
{
  basic_eval_data_d *BED_d = NULL;
  basic_eval_data_mp *BED_mp = NULL;

  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    BED_d = (basic_eval_data_d *)ED_d;
    getDehomPoint_d(out_d, in_d, in_d->size, &BED_d->preProcData);
  }
  else
  { // compute out_mp
    BED_mp = (basic_eval_data_mp *)ED_mp;
    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
    getDehomPoint_mp(out_mp, in_mp, in_mp->size, &BED_mp->preProcData);
  }

  BED_d = NULL;
  BED_mp = NULL;

  return 0;
}

int zero_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy in to out                                         *
\***************************************************************/
{
  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    out_d->size = 0;
  }
  else
  { // compute out_mp
    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
    out_mp->size = 0;
  }

  return 0;
}



