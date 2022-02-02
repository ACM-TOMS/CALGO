// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

void getDehomPoint_d(point_d dehomPoint, point_d inPoint, int num_vars, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the dehom point                             *
\***************************************************************/
{
  int i, j, var = PPD->num_var_gp, count = 0, num_gps = PPD->num_var_gp + PPD->num_hom_var_gp, size = 0;
  comp_d recip; 

  // setup size
  size = inPoint->size - var;

  // make sure dehomPoint is large enough
  increase_size_point_d(dehomPoint, size);
  dehomPoint->size = size;

  // check whether using user-defined homotopy
  if (num_gps > 0)
  { // loop over each variable group
    for (i = 0; i < num_gps; i++)  
    {
      if (PPD->type[i] == 0)
      { // already homoginized - copy over
        for (j = 0; j < PPD->size[i]; j++)
        {
          set_d(&dehomPoint->coord[var - PPD->num_var_gp], &inPoint->coord[var]);
          var++;
        }
      }
      else
      { // make sure homogenous coordinate is a valid number
        if (isnan(inPoint->coord[count].r) || isnan(inPoint->coord[count].i) || isinf(inPoint->coord[count].r) || isinf(inPoint->coord[count].i))
        { // not a valid number
          for (j = 0; j < PPD->size[i]; j++)
          {
            set_double_d(&dehomPoint->coord[var - PPD->num_var_gp], -1e199, -1e199);
            var++;
          }
        }
        else
        { // we have a number - determine if it is 0
          if (inPoint->coord[count].r == 0 && inPoint->coord[count].i == 0)
          { // generate a random perturbation so that we can divide
            get_comp_rand_d(recip);
            mul_rdouble_d(recip, recip, 1e-16);
            recip_d(recip, recip);
          }
          else
          { // reciprocate the homogenous coordinate
            recip_d(recip, &inPoint->coord[count]);
          }
          // find the dehom coordinates
          for (j = 0; j < PPD->size[i]; j++)  
          { // multiply by recip to produce de-hom coord.
            mul_d(&dehomPoint->coord[var - PPD->num_var_gp], &inPoint->coord[var], recip);
            var++;
          }
        }
        count++;
      }
    }
  }
  else
  { // simply copy value over when using user-defined homotopy 
    point_cp_d(dehomPoint, inPoint);
  }

  return;
}

void nonsolutions_check_compare_d(int *isZero, point_d f1, point_d f2, int startFunc, int endFunc, double zero_threshold, double max_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  double n1, n2, tol;

  // setup threshold based on given threshold and precision
  tol = MAX(zero_threshold, 1e-15);

  for (i = startFunc; i < endFunc; i++)
  { // test ith function
    isZero[i] = 1;
    n1 = d_abs_d(&f1->coord[i]);
    n2 = d_abs_d(&f2->coord[i]);

    if (tol <= n1 && n1 <= n2)
    { // compare ratio
      if (n1 > max_ratio * n2)
        isZero[i] = 0;
    }
    else if (tol <= n2 && n2 <= n1)
    { // compare ratio
      if (n2 > max_ratio * n1)
        isZero[i] = 0;
    }
  }

  return;
}

int nonsolutions_check_d(int size_f, int size_r, point_d x, point_d last_x, comp_d time, double zero_threshold, double max_ratio, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int isSoln = 1;

  if (size_f > size_r)
  { // we had to square up the system so we need to check to make sure the answer is an actual solution
    int i;
    double n1, n2, tol;
    point_d f;
    eval_struct_d e;

    // setup threshold based on given threshold and precision
    tol = MAX(zero_threshold, 1e-15);

    init_point_d(f, Prog->numFuncs);
    init_eval_struct_d(e, Prog->numFuncs, Prog->numVars, Prog->numPars);

    evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, x, time, Prog);
    evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, last_x, time, Prog);

    // compare the function values
    isSoln = 1;
    for (i = 0; i < Prog->numFuncs && isSoln; i++)
    {
      n1 = d_abs_d(&e.funcVals->coord[i]);
      n2 = d_abs_d(&f->coord[i]);
      if (tol <= n1 && n1 <= n2)
      { // compare ratio
        if (n1 > max_ratio * n2)
          isSoln = 0;
      }
      else if (tol <= n2 && n2 <= n1)
      { // compare ratio
        if (n2 > max_ratio * n1)
          isSoln = 0;
      }
    }

    // clear
    clear_point_d(f);
    clear_eval_struct_d(e);
  }

  return isSoln;
}

int zero_dim_main_d(int MPType, double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *OUT = NULL, *StartPts = NULL, *FAIL = fopen("failed_paths", "w"), *midOUT = NULL, *rawOUT = fopen("raw_data", "w");
  tracker_config_t T;
  prog_t dummyProg;
  bclock_t time1, time2;
  int i, num_variables = 0, userHom = 0, paramHom = 0, pathMod = 0, convergence_failures = 0, sharpening_failures = 0, sharpening_singular = 0, num_crossings = 0, num_sols = 0;
  int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, usedEq = 0;
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  basic_eval_data_d ED;
  trackingStats trackCount;
  char inputName[] = "this_input";
  double midpoint_tol, track_time, intrinsicCutoffMultiplier;

  bclock(&time1); // initialize the clock.
  init_trackingStats(&trackCount); // initialize trackCount to all 0

  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType); 

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

#ifdef _OPENMP
  #pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }

  // call the setup function
  if (userHom == 1)
  { // setup for user defined homotopy
    num_variables = userHom_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &ptr_to_eval_d, &ptr_to_eval_mp);

    // error checking
    if (num_variables > ED.squareSystem.Prog->numFuncs)
    { // underdetermined
      printf("ERROR: The user-defined homotopy is underdetermined (%d variable%s, %d function%s)!\n       Please add additional functions or remove variables.\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "");
      bexit(ERROR_INPUT_SYSTEM);
    }
    else if (num_variables < ED.squareSystem.Prog->numFuncs)
    { // overdetermined
      printf("WARNING: The user-defined homotopy is overdetermined (%d variable%s, %d function%s).\n         Path tracking will use a least squares approach.\n\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "");
    }
  }
  else if (paramHom == 2 || userHom == 2)
  { // setup for parameter homotopy or user-defined homotopy with patches
    num_variables = paramHom_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &ptr_to_eval_d, &ptr_to_eval_mp, "preproc_data", 1, startName);

    if (num_variables > ED.squareSystem.Prog->numFuncs + ED.patch.num_patches)
    { // underdetermined
      if (paramHom == 2)
        printf("ERROR: The parameter homotopy is underdetermined (%d variable%s, %d function%s, %d patch%s)!\n       Please add additional functions or remove variables.\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "", ED.patch.num_patches, ED.patch.num_patches > 1 ? "es" : "");
      else
        printf("ERROR: The user-defined homotopy is underdetermined (%d variable%s, %d function%s, %d patch%s)!\n       Please add additional functions or remove variables.\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "", ED.patch.num_patches, ED.patch.num_patches > 1 ? "es" : "");
      bexit(ERROR_INPUT_SYSTEM);
    }
    else if (num_variables < ED.squareSystem.Prog->numFuncs + ED.patch.num_patches)
    { // overdetermined
      if (paramHom == 2)
        printf("WARNING: The parameter homotopy is overdetermined (%d variable%s, %d function%s, %d patch%s).\n         Path tracking will use a least squares approach.\n\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "", ED.patch.num_patches, ED.patch.num_patches > 1 ? "es" : "");
      else
        printf("WARNING: The user-defined homotopy is overdetermined (%d variable%s, %d function%s, %d patch%s).\n         Path tracking will use a least squares approach.\n\n", num_variables, num_variables > 1 ? "s" : "", ED.squareSystem.Prog->numFuncs, ED.squareSystem.Prog->numFuncs > 1 ? "s" : "", ED.patch.num_patches, ED.patch.num_patches > 1 ? "es" : "");
    }
  }
  else if (userHom == -59)
  { // setup for equation by equation
    num_variables = eqbyeq_setup_d(&OUT, "output", &T, &ED, &dummyProg, &ptr_to_eval_d, &ptr_to_eval_mp, "preproc_data", "deg.out", "nonhom_start", "start", "depth", intrinsicCutoffMultiplier);
  }
  else
  { // setup for standard tracking - 'useRegen' is used to determine whether or not to setup 'start'
    num_variables = zero_dim_basic_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow, &ptr_to_eval_d, &ptr_to_eval_mp, "preproc_data", "deg.out", !useRegen, "nonhom_start", "start");
  }

  // error checking
  if (userHom <= 0 && paramHom != 2)
  { // no pathvariables or parameters allowed!
    if (dummyProg.numPathVars > 0)
    { // path variable present
      printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
    if (dummyProg.numPars > 0)
    { // parameter present
      printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
  }

  if (T.MPType == 2)  //If we are doing adaptive precision path-tracking, we must set up AMP_eps, AMP_Phi, AMP_Psi based on config settings.
  {
    T.AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T.AMP_Phi = T.AMP_bound_on_degree*(T.AMP_bound_on_degree-1.0)*T.AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
    T.AMP_Psi = T.AMP_bound_on_degree*T.AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
    // initialize latest_newton_residual_mp to the maximum precision
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);
  }

  if (userHom > 0 || paramHom == 2 || (userHom == 0 && useRegen == 0))
  { // we are using either a user-defined homotopy or a standard homotopy (not regeneration/eq-by-eq)
  
    // setup the file containing the start points
    if (paramHom == 2 || userHom == 2)
    { // parameter homotopy or m-hom user-defined homotopy - start points have been moved to patch
      StartPts = fopen("start_param_hom", "r");
      if (StartPts == NULL)
      {
        printf("\n\nERROR: 'start_param_hom' does not exist!!!\n\n");
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }
    else if (userHom == 1)
    { // affine user-defined homotopy - start points in given file
      StartPts = fopen(startName, "r");
      if (StartPts == NULL)
      {
        printf("\n\nERROR: '%s' does not exist!!!\n\n", startName);
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }
    else
    { // internally constructed homotopy
      StartPts = fopen("start", "r");  
      if (StartPts == NULL)
      {
        printf("\n\nERROR: 'start' does not exist!!!\n\n");
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }

    // do the actual tracking now that everything is setup
    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      worker_info sendType;
      if (userHom == 1)
        sendType.dataType = USER_HOM_D;
      else if (userHom == 2 || paramHom == 2)
        sendType.dataType = PARAM_HOM_D;
      else
        sendType.dataType = ZERO_DIM_D;
      bcast_worker_info(&sendType, my_id, headnode);

      if (T.endgameNumber == 3)
        head_zero_dim_endgame(T.endgameNumber == 3, T.minCycleTrackBack, &trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, my_id, num_processes, headnode);
      else
        head_zero_dim_track_d(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, my_id, num_processes, headnode);
#endif
    }
    else
    { // use sequential (which includes OpenMP)
      if (T.endgameNumber == 3)
      { // use the track-back endgame
        zero_dim_trackBack_d(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_basic_eval_prec, zero_dim_dehom);
      }
      else
      { // use regular endgame
        zero_dim_track_d(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_basic_eval_prec, zero_dim_dehom);
      }
    }

    fclose(StartPts);
    fclose(midOUT);
    // delete the file containing the start points
    if (paramHom == 2 || userHom == 2)
      remove("start_param_hom");

    // finish the output to rawOUT
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT

    // check for path crossings
    midpoint_checker(trackCount.numPoints, num_variables, midpoint_tol, &num_crossings);

    // setup num_sols
    num_sols = trackCount.successes;
  }
  else if (userHom == -59)
  { // do the eq-by-eq tracking now that everything is setup
    usedEq = 1;

    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      worker_info sendType;
      sendType.dataType = EQBYEQ_D;
      bcast_worker_info(&sendType, my_id, headnode);

      head_eqbyeq_track_d(&trackCount, OUT, rawOUT, FAIL, "midpath_data", pathMod, &T, midpoint_tol, T.final_tol_times_mult, &ED, ED.BED_mp, my_id, num_processes, headnode);
#endif
    }
    else
    { // use sequential
      eqbyeq_track_d(OUT, rawOUT, FAIL, "midpath_data", pathMod, &T, midpoint_tol, T.final_tol_times_mult, &ED, ED.BED_mp, &trackCount);
    }

    // finish the output to rawOUT
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT

    // check for path crossings - need to run the checker with the proper number of variables - '-1' since we only have 1 variable group
    if (ED.EqD->num_subsystems == 1) // tracking during witness generation
      midpoint_checker(trackCount.numPoints, num_variables - 1, midpoint_tol, &num_crossings);
    else if (ED.EqD->stageData_d[ED.EqD->num_subsystems - 1].useIntrinsicSlice)
      midpoint_checker(trackCount.numPoints, num_variables - 1, midpoint_tol, &num_crossings);
    else
      midpoint_checker(trackCount.numPoints, 2 * num_variables, midpoint_tol, &num_crossings);

    // setup num_sols
    num_sols = trackCount.successes;
  }
  else
  { // do the regeneration tracking
    usedEq = 1;

    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      regen_t regen;
      worker_info sendType;
      sendType.dataType = ZERO_DIM_REGEN;
      bcast_worker_info(&sendType, my_id, headnode);

      // setup regen
      setup_regen_from_zero_dim_seq(1, &regen, "startRegen", regenStartLevel, intrinsicCutoffMultiplier, "depth", "deg.out", &T, &ED, ED.BED_mp, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

      // do the actual tracking
      head_regen_track_zero_dim(regenStartLevel, &regen, &T, pathMod, midpoint_tol, "startRegen", OUT, "midpath_data", FAIL, rawOUT, &trackCount, &ED, ED.BED_mp, my_id, num_processes, headnode);

      // check for path crossings - need to run the checker with the proper number of variables
      if (regen.level[regen.num_levels - 1].useIntrinsicSlice)
        midpoint_checker(regen.level[regen.num_levels - 1].num_paths, regen.level[regen.num_levels - 1].level + regen.level[regen.num_levels - 1].depth, midpoint_tol, &num_crossings);
      else
        midpoint_checker(regen.level[regen.num_levels - 1].num_paths, num_variables, midpoint_tol, &num_crossings);

      // setup num_sols
      num_sols = regen.level[regen.num_levels - 1].num_nonsing;

      // clear regen
      clearRegenRandom_zero_dim(1, &regen, T.MPType);
#endif
    }
    else
    { // use sequential - setup the regeneration structure based on the number of threads
      int max = max_threads();
      regen_t *regen = (regen_t *)bmalloc(max * sizeof(regen_t));

      // setup regen 
      setup_regen_from_zero_dim_seq(max, regen, "startRegen", regenStartLevel, intrinsicCutoffMultiplier, "depth", "deg.out", &T, &ED, ED.BED_mp, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

      // do the actual tracking
      if (T.MPType == 0)
        regen_track_seq(regenStartLevel, regen, &T, pathMod, midpoint_tol, "startRegen", OUT, "midpath_data", FAIL, rawOUT, &trackCount, ED.patch.patchCoeff, NULL, NULL, &dummyProg);
      else
        regen_track_seq(regenStartLevel, regen, &T, pathMod, midpoint_tol, "startRegen", OUT, "midpath_data", FAIL, rawOUT, &trackCount, ED.patch.patchCoeff, ED.BED_mp->patch.patchCoeff, ED.BED_mp->patch.patchCoeff_rat, &dummyProg);

      // check for path crossings - need to run the checker with the proper number of variables
      if (regen[0].level[regen[0].num_levels - 1].useIntrinsicSlice)
        midpoint_checker(regen[0].level[regen[0].num_levels - 1].num_paths, regen[0].level[regen[0].num_levels - 1].level + regen[0].level[regen[0].num_levels - 1].depth, midpoint_tol, &num_crossings);
      else
        midpoint_checker(regen[0].level[regen[0].num_levels - 1].num_paths, regen[0].num_variables, midpoint_tol, &num_crossings);

      // setup num_sols
      num_sols = regen[0].level[regen[0].num_levels - 1].num_nonsing;

      // clear regen
      clearRegenRandom_zero_dim(max, regen, T.MPType);
      free(regen);
    }
  }

  // we report how we did with all paths:
  bclock(&time2);
  totalTime(&track_time, time1, time2);
  if (T.screenOut)
  {
    printf("Number of failures:  %d\n", trackCount.failures);
    printf("Number of successes:  %d\n", trackCount.successes);
    printf("Number of paths:  %d\n", trackCount.numPoints);
    printf("Parse Time = %fs\n", parse_time);
    printf("Track Time = %fs\n", track_time);
  }
  fprintf(OUT, "Number of failures:  %d\n", trackCount.failures);
  fprintf(OUT, "Number of successes:  %d\n", trackCount.successes);
  fprintf(OUT, "Number of paths:  %d\n", trackCount.numPoints);
  fprintf(OUT, "Parse Time = %fs\n", parse_time);
  fprintf(OUT, "Track Time = %fs\n", track_time);

  // print the system to rawOUT 
  printZeroDimRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);

  // close all of the files
  fclose(OUT);
  fclose(rawOUT);
  fprintf(FAIL, "\n");
  fclose(FAIL);

  // reproduce the input file needed for this run
  reproduceInputFile(inputName, "func_input", &T, 0, 0, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // print the output
  if (userHom == -59)
  { // print the eq-by-eq output chart to the screen
    eqbyeqOutputChart_d(ED.EqD, stdout, T.regen_remove_inf);
  }

  // do the standard post-processing
  sort_points(num_crossings, &convergence_failures, &sharpening_failures, &sharpening_singular, inputName, num_sols, num_variables, midpoint_tol, T.final_tol_times_mult, &T, &ED.preProcData, useRegen == 1 && userHom == 0, userHom == -59);

  if (useRegen == 1 && userHom == 0)
  {
    printf("\nSummary for the last regeneration level:");
  }
  else if (userHom == -59)
  {
    printf("\nSummary for the last stage of equation-by-equation:");
  }

  // print the failure summary
  printFailureSummary(&trackCount, convergence_failures, sharpening_failures, sharpening_singular);

  // clear memory
  if (userHom == 0 && paramHom == 0)
  {
    if (ED.squareSystem.Prog->numSubfuncs > 0)
    { // clear subFuncsBelow
      for (i = ED.squareSystem.Prog->numFuncs - 1; i >= 0; i--)
        free(subFuncsBelow[i]);
      free(subFuncsBelow);
    }
    free(startSub);
    free(endSub);
    free(startFunc);
    free(endFunc);
    free(startJvsub);
    free(endJvsub);
    free(startJv);
    free(endJv);
  }
  basic_eval_clear_d(&ED, userHom, T.MPType);
  tracker_config_clear(&T);
  clearMP();

  return 0;
}

void getDehomPoint_mp(point_mp dehomPoint, point_mp inPoint, int num_vars, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the dehom point                             *
\***************************************************************/
{
  int i, j, var = PPD->num_var_gp, count = 0, num_gps = PPD->num_var_gp + PPD->num_hom_var_gp, num_digits = prec_to_digits(inPoint->curr_prec), size = 0;
  mpf_t epsilon;
  comp_mp recip;

  mpf_init(epsilon);
  init_mp(recip);

  // setup size
  size = inPoint->size - var;
  
  // make sure dehomPoint is large enough
  increase_size_point_mp(dehomPoint, size);
  // setup dehomPoint->size
  dehomPoint->size = size;

  if (num_gps > 0)
  { // not using userHom
    // loop over the variable groups
    for (i = 0; i < num_gps; i++)  
    {
      if (PPD->type[i] == 0) 
      { // already homogenized - just copy over
        for (j = 0; j < PPD->size[i]; j++)
        {
          set_mp(&dehomPoint->coord[var - PPD->num_var_gp], &inPoint->coord[var]);
          var++;
        }
      }
      else
      { // make sure homogenous coordinate is a valid number
        if (!mpfr_number_p(inPoint->coord[count].r) || !mpfr_number_p(inPoint->coord[count].i))
        { // not a valid number
          for (j = 0; j < PPD->size[i]; j++)
          {
            set_double_mp(&dehomPoint->coord[var - PPD->num_var_gp], -1e199, -1e199);
            var++;
          }
        }
        else
        { // we have a number - determine if it is 0
          if (mpfr_zero_p(inPoint->coord[count].r) && mpfr_zero_p(inPoint->coord[count].i)) // make sure that we do not have zero in the homogenous coordinate
          { // generate a random perturbation so that we can divide
            get_comp_rand_mp(recip);
            mpfr_ui_pow_ui(epsilon, 10, num_digits, __gmp_default_rounding_mode);
            mpf_ui_div(epsilon, 1, epsilon);
            mul_rmpf_mp(recip, recip, epsilon);
            recip_mp(recip, recip);
          }
          else
          { // reciprocate the homogeneous coordinate
            recip_mp(recip, &inPoint->coord[count]);
          }
          // find the dhom coordinates
          for (j = 0; j < PPD->size[i]; j++)
          { // mulitply by recip to produce de-hom coord
            mul_mp(&dehomPoint->coord[var - PPD->num_var_gp], &inPoint->coord[var], recip);
            var++;
          }
        }
        count++;
      }
    }
  }
  else
  { // user-defined homotopy - simply copy over
    point_cp_mp(dehomPoint, inPoint);
  }

  mpf_clear(epsilon);
  clear_mp(recip);

  return;
}

void nonsolutions_check_compare_mp(int *isZero, point_mp f1, point_mp f2, int startFunc, int endFunc, double zero_threshold, double max_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, num_digits = prec_to_digits((int) mpf_get_default_prec());
  double tol;
  mpf_t n1, n2, zero_thresh, max_rat;
  mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
  mpf_set_d(max_rat, max_ratio);

  // setup threshold based on given threshold and precision
  if (num_digits > 300)
    num_digits = 300;
  num_digits -= 2;
  tol = MAX(zero_threshold, pow(10,-num_digits));
  mpf_set_d(zero_thresh, tol);

  for (i = startFunc; i < endFunc; i++)
  { // test ith function
    isZero[i] = 1;
    mpf_abs_mp(n1, &f1->coord[i]);
    mpf_abs_mp(n2, &f2->coord[i]);

    if (mpf_cmp(zero_thresh, n1) <= 0 && mpf_cmp(n1, n2) <= 0)
    { // compare ratio
      mpf_mul(n2, max_rat, n2);
      if (mpf_cmp(n1, n2) > 0)
        isZero[i] = 0;
    }
    else if (mpf_cmp(zero_thresh, n2) <= 0 && mpf_cmp(n2, n1) <= 0)
    { // compare ratio
      mpf_mul(n1, max_rat, n1);
      if (mpf_cmp(n2, n1) > 0)
        isZero[i] = 0;
    }
  }

  mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);

  return;
}

int nonsolutions_check_mp(int size_f, int size_r, point_mp x, point_mp last_x, comp_mp time, double zero_threshold, double max_ratio, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int isSoln = 1;
  
  if (size_f > size_r)
  { // we had to square up the system so we need to check to make sure the answer is an actual solution
    int i, num_digits = prec_to_digits((int) mpf_get_default_prec());
    double tol;
    mpf_t n1, n2, zero_thresh, max_rat;
    point_mp f;
    eval_struct_mp e;

    mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
    init_point_mp(f, Prog->numFuncs);
    init_eval_struct_mp(e, Prog->numFuncs, Prog->numVars, Prog->numPars);

    mpf_set_d(max_rat, max_ratio);

    // setup threshold based on given threshold and precision
    if (num_digits > 300)
      num_digits = 300;
    num_digits -= 2;
    tol = MAX(zero_threshold, pow(10,-num_digits));
    mpf_set_d(zero_thresh, tol);

    evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, x, time, Prog);
    evalProg_mp(f, e.parVals, e.parDer, e.Jv, e.Jp, last_x, time, Prog);

    // compare the function values
    isSoln = 1;
    for (i = 0; i < Prog->numFuncs && isSoln; i++)
    {
      mpf_abs_mp(n1, &e.funcVals->coord[i]);
      mpf_abs_mp(n2, &f->coord[i]);

      if (mpf_cmp(zero_thresh, n1) <= 0 && mpf_cmp(n1, n2) <= 0)
      { // compare ratio
        mpf_mul(n2, max_rat, n2);
        if (mpf_cmp(n1, n2) > 0)
          isSoln = 0;
      }
      else if (mpf_cmp(zero_thresh, n2) <= 0 && mpf_cmp(n2, n1) <= 0)
      { // compare ratio
        mpf_mul(n1, max_rat, n1);
        if (mpf_cmp(n2, n1) > 0)
          isSoln = 0;
      }
    }

    // clear e
    mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
    clear_point_mp(f);
    clear_eval_struct_mp(e);
  }

  return isSoln;
}

int zero_dim_main_mp(double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *OUT = NULL, *StartPts = NULL, *FAIL = fopen("failed_paths", "w"), *midOUT = NULL, *rawOUT = fopen("raw_data", "w");
  tracker_config_t T;
  int i, num_variables = 0, userHom = 0, paramHom = 0, pathMod = 0, convergence_failures = 0, sharpening_failures = 0, sharpening_singular = 0, num_crossings = 0, num_sols = 0;
  int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, usedEq = 0;
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  prog_t dummyProg;
  bclock_t time1, time2;
  int (*ptr_to_eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  basic_eval_data_mp ED;
  trackingStats trackCount;
  char inputName[] = "this_input";
  double midpoint_tol, track_time, intrinsicCutoffMultiplier;

  bclock(&time1); // initialize the clock.
  init_trackingStats(&trackCount); // initialize trackCount to all 0

  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, 1); //Set up T, the tracker_config_t, using the config file.

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

#ifdef _OPENMP
  #pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }

  // initialize latest_newton_residual_mp
  mpf_init(T.latest_newton_residual_mp);

  // call the setup function based on these parameters
  if (userHom == 1)
  { // setup for user defined homotopy
    num_variables = userHom_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &ptr_to_eval_func);
  }
  else if (paramHom == 2 || userHom == 2)
  { // setup for parameter homotopy
    num_variables = paramHom_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &ptr_to_eval_func, "preproc_data", 1, startName);
  }
  else if (userHom == -59)
  { // setup for eq-by-eq
    num_variables = eqbyeq_setup_mp(&OUT, "output", &T, &ED, &dummyProg, &ptr_to_eval_func, "preproc_data", "deg.out", "nonhom_start", "start", "depth", intrinsicCutoffMultiplier);
  }
  else
  { // setup for standard tracking - 'useRegen' is used to determine whether or not to setup 'start'
    num_variables = zero_dim_basic_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED, &dummyProg, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow, &ptr_to_eval_func, "preproc_data", "deg.out", !useRegen, "nonhom_start", "start");
  }

  // error checking
  if (userHom <= 0 && paramHom != 2)
  { // no pathvariables or parameters allowed!
    if (dummyProg.numPathVars > 0)
    { // path variable present
      printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
    if (dummyProg.numPars > 0)
    { // parameter present
      printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
  }

  if (userHom > 0 || paramHom == 2 || (userHom == 0 && useRegen == 0))
  { // we are using either a user-defined homotopy or a standard homotopy (not regeneration/eq-by-eq)

    // setup the file containing the start points
    if (paramHom == 2 || userHom == 2)
    { // parameter homotopy or m-hom user-defined homotopy - start points have been moved to patch
      StartPts = fopen("start_param_hom", "r");
      if (StartPts == NULL)
      {
        printf("\n\nERROR: 'start_param_hom' does not exist!!!\n\n");
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }
    else if (userHom == 1)
    { // affine user-defined homotopy - start points in given file
      StartPts = fopen(startName, "r");
      if (StartPts == NULL)
      {
        printf("\n\nERROR: '%s' does not exist!!!\n\n", startName);
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }
    else
    { // internally constructed homotopy
      StartPts = fopen("start", "r");
      if (StartPts == NULL)
      {
        printf("\n\nERROR: 'start' does not exist!!!\n\n");
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }

    // do the actual tracking now that everything is setup
    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      worker_info sendType;
      if (userHom == 1)
        sendType.dataType = USER_HOM_MP;
      else if (userHom == 2 || paramHom == 2)
        sendType.dataType = PARAM_HOM_MP;
      else
        sendType.dataType = ZERO_DIM_MP;
      bcast_worker_info(&sendType, my_id, headnode);

      if (T.endgameNumber == 3)
        head_zero_dim_endgame(T.endgameNumber == 3, T.minCycleTrackBack, &trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, NULL, &ED, my_id, num_processes, headnode);
      else
        head_zero_dim_track_mp(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ptr_to_eval_func, my_id, num_processes, headnode);
#endif
    }
    else
    { // use sequential (which includes OpenMP)
      if (T.endgameNumber == 3)
      { // use the track-back endagme
        zero_dim_trackBack_mp(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ptr_to_eval_func, zero_dim_dehom);
      }
      else
      { // use regular endgame
        zero_dim_track_mp(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ptr_to_eval_func, zero_dim_dehom);
      }
    }

    fclose(StartPts);
    fclose(midOUT);
    // delete the file containing the start points
    if (paramHom == 2 || userHom == 2)
      remove("start_param_hom");

    // finish the output to rawOUT
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT

    // check for path crossings
    midpoint_checker(trackCount.numPoints, num_variables, midpoint_tol, &num_crossings);

    // setup num_sols
    num_sols = trackCount.successes;
  }
  else if (userHom == -59)
  { // do the eq-by-eq tracking now that everything is setup
    usedEq = 1;

    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      worker_info sendType;
      sendType.dataType = EQBYEQ_MP;
      bcast_worker_info(&sendType, my_id, headnode);

      head_eqbyeq_track_mp(&trackCount, OUT, rawOUT, FAIL, "midpath_data", pathMod, &T, midpoint_tol, T.final_tol_times_mult, &ED, my_id, num_processes, headnode);
#endif
    }
    else
    { // use sequential
      eqbyeq_track_mp(OUT, rawOUT, FAIL, "midpath_data", pathMod, &T, midpoint_tol, T.final_tol_times_mult, &ED, &trackCount);
    }

    // finish the output to rawOUT
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT

    // check for path crossings - need to run the checker with the proper number of variables - '-1' since we only have 1 variable group
    if (ED.EqD->num_subsystems == 1) // tracking during witness generation
      midpoint_checker(trackCount.numPoints, num_variables - 1, midpoint_tol, &num_crossings);
    else if (ED.EqD->stageData_mp[ED.EqD->num_subsystems - 1].useIntrinsicSlice)
      midpoint_checker(trackCount.numPoints, num_variables - 1, midpoint_tol, &num_crossings);
    else
      midpoint_checker(trackCount.numPoints, 2 * num_variables, midpoint_tol, &num_crossings);

    // setup num_sols
    num_sols = trackCount.successes;
  }
  else
  { // do regeneration tracking
    usedEq = 1;

    if (num_processes > 1)
    { // using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
      regen_t regen;
      worker_info sendType;
      sendType.dataType = ZERO_DIM_REGEN;
      bcast_worker_info(&sendType, my_id, headnode);

      // setup regen
      setup_regen_from_zero_dim_seq(1, &regen, "startRegen", regenStartLevel, intrinsicCutoffMultiplier, "depth", "deg.out", &T, NULL, &ED, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

      // do the actual tracking
      head_regen_track_zero_dim(regenStartLevel, &regen, &T, pathMod, midpoint_tol, "startRegen", OUT, "midpath_data", FAIL, rawOUT, &trackCount, NULL, &ED, my_id, num_processes, headnode);

      // check for path crossings - need to run the checker with the proper number of variables
      if (regen.level[regen.num_levels - 1].useIntrinsicSlice)
        midpoint_checker(regen.level[regen.num_levels - 1].num_paths, regen.level[regen.num_levels - 1].level + regen.level[regen.num_levels - 1].depth, midpoint_tol, &num_crossings);
      else
        midpoint_checker(regen.level[regen.num_levels - 1].num_paths, num_variables, midpoint_tol, &num_crossings);

      // setup num_sols
      num_sols = regen.level[regen.num_levels - 1].num_nonsing;

      // clear regen
      clearRegenRandom_zero_dim(1, &regen, T.MPType);
#endif
    }
    else
    { // use sequential - setup the regeneration structure based on the number of threads
      int max = max_threads();
      regen_t *regen = (regen_t *)bmalloc(max * sizeof(regen_t));

      // setup regen
      setup_regen_from_zero_dim_seq(max, regen, "startRegen", regenStartLevel, intrinsicCutoffMultiplier, "depth", "deg.out", &T, NULL, &ED, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv , endJv, subFuncsBelow);

      // do the actual tracking
      regen_track_seq(regenStartLevel, regen, &T, pathMod, midpoint_tol, "startRegen", OUT, "midpath_data", FAIL, rawOUT, &trackCount, NULL, ED.patch.patchCoeff, ED.patch.patchCoeff_rat, &dummyProg);

      // check for path crossings - need to run the checker with the proper number of variables
      if (regen[0].level[regen[0].num_levels - 1].useIntrinsicSlice)
        midpoint_checker(regen[0].level[regen[0].num_levels - 1].num_paths, regen[0].level[regen[0].num_levels - 1].level + regen[0].level[regen[0].num_levels - 1].depth, midpoint_tol, &num_crossings);
      else
        midpoint_checker(regen[0].level[regen[0].num_levels - 1].num_paths, regen[0].num_variables, midpoint_tol, &num_crossings);

      // setup num_sols
      num_sols = regen[0].level[regen[0].num_levels - 1].num_nonsing;

      // clear regen
      clearRegenRandom_zero_dim(max, regen, T.MPType);
      free(regen);
    }
  }

  // report how this run did
  bclock(&time2);
  totalTime(&track_time, time1, time2);
  if (T.screenOut)
  {
    printf("Number of failures:  %d\n", trackCount.failures);
    printf("Number of successes:  %d\n", trackCount.successes);
    printf("Number of paths:  %d\n", trackCount.numPoints);
    printf("Parse Time = %fs\n", parse_time);
    printf("Track Time = %fs\n", track_time);
  }
  fprintf(OUT, "Number of failures:  %d\n", trackCount.failures);
  fprintf(OUT, "Number of successes:  %d\n", trackCount.successes);
  fprintf(OUT, "Number of paths:  %d\n", trackCount.numPoints);
  fprintf(OUT, "Parse Time = %fs\n", parse_time);
  fprintf(OUT, "Track Time = %fs\n", track_time);

  // print the system to rawOUT
  printZeroDimRelevantData(NULL, &ED, T.MPType, usedEq, rawOUT);

  // close all of the files
  fclose(OUT);
  fclose(rawOUT);
  fprintf(FAIL, "\n");
  fclose(FAIL);

  // reproduce the input file needed for this run
  reproduceInputFile(inputName, "func_input", &T, 0, 0, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // print the output
  if (userHom == -59)
  { // print the eq-by-eq output chart to the screen
    eqbyeqOutputChart_mp(ED.EqD, stdout, T.regen_remove_inf);
  }

  // do the standard post-processing
  sort_points(num_crossings, &convergence_failures, &sharpening_failures, &sharpening_singular, inputName, num_sols, num_variables, midpoint_tol, T.final_tol_times_mult, &T, &ED.preProcData, useRegen == 1 && userHom == 0, userHom == -59);

  if (useRegen == 1 && userHom == 0)
  {
    printf("\nSummary for the last regeneration level:");
  }
  else if (userHom == -59)
  {
    printf("\nSummary for the last stage of equation-by-equation:");
  }

  // print the failure summary
  printFailureSummary(&trackCount, convergence_failures, sharpening_failures, sharpening_singular);

  // clear the memory
  if (userHom == 0 && paramHom == 0)
  {
    if (ED.squareSystem.Prog->numSubfuncs > 0)
    { // clear subFuncsBelow
      for (i = ED.squareSystem.Prog->numFuncs - 1; i >= 0; i--)
        free(subFuncsBelow[i]);
      free(subFuncsBelow);
    }
    free(startSub);
    free(endSub);
    free(startFunc);
    free(endFunc);
    free(startJvsub);
    free(endJvsub);
    free(startJv);
    free(endJv);
  }
  basic_eval_clear_mp(&ED, userHom, 1);
  clearMP();
  tracker_config_clear(&T);

  return 0;
}

void clearProg(prog_t *P, int MPType, int clearEvalProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears P                                               *
\***************************************************************/
{
  if (clearEvalProg)
  { // free the structures to evaluate the program
    freeEvalProg(MPType);
  }

  // clear nums
  clearNums(&P->nums, P->numNums); 

  free(P->prog);
  free(P->var_gp_sizes);
  
  return;
}

void clearNums(num_t **nums, int numNums)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears nums                                            *
\***************************************************************/
{
  int i;

  for (i = numNums - 1; i >= 0; i--)
  { // clear the nums
    mpf_clear((*nums)[i].real);
    mpq_clear((*nums)[i].rat);
  }
  free(*nums);

  return;
}

void setupStart_d(tracker_config_t *T, point_data_d *PD, FILE *StartPts)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  double a, b;

  // make sure PD->point is large enough
  increase_size_point_d(PD->point, T->numVars);
  PD->point->size = T->numVars;

  for (i = 0; i < T->numVars; i++)
  {
    fscanf(StartPts, "%lf%lf", &a, &b);
    scanRestOfLine(StartPts);
    set_double_d(&PD->point->coord[i], a, b);
  }

  // finish setup
  set_one_d(PD->time);
  T->endgameOnly = 0;

  return;
}

void setupStart_mp(tracker_config_t *T, point_data_mp *PD, FILE *StartPts)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  // make sure PD->point is large enough
  increase_size_point_mp(PD->point, T->numVars);
  PD->point->size = T->numVars;

  for (i = 0; i < T->numVars; i++)
  {
    mpf_inp_str(PD->point->coord[i].r, StartPts, 10);
    mpf_inp_str(PD->point->coord[i].i, StartPts, 10);
    scanRestOfLine(StartPts);
  }

  set_one_mp(PD->time);

  return;
}

void remove_output_files(int trackType, int sharpenOnly, int removeRawData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: removes useless output files                           *
\***************************************************************/
{
  remove("isosingular_summary");
  if (sharpenOnly)
  { // the sharpending module can create every zero dimensional file but failed paths
    remove("nonsolutions");
    remove("failed_paths");
    remove("raw_data_old");

    remove("main_cascade_data");
    remove("cascade_output");
    remove("raw_cascade_data");
    if (trackType == 0)
    {
      remove("witness_data");
      remove("witness_superset");
    }
  }
  else
  { // doing a standard run
    remove("nonsolutions");
    remove("real_solutions");
    remove("real_finite_solutions");
    remove("finite_solutions");
    remove("infinite_solutions");
    remove("nonsingular_solutions");
    remove("singular_solutions");
    remove("main_data");
    remove("output");
    remove("real_project");
    remove("raw_solutions");
    if (removeRawData)
    {
      remove("raw_data");
      remove("failed_paths");
    }

    remove("midpath_data");

    remove("main_cascade_data");
    remove("cascade_output");
    remove("raw_cascade_data");
    if (trackType < 2)
    {
      remove("witness_data");
      remove("witness_superset");
    }
  }

  // remove old output_'id' files if they exist
  FILE *testFILE = NULL;
  char *str = NULL;
  int size, count = 1;

  // try to open first output file
  size = 1 + snprintf(NULL, 0, "output_%d", count);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "output_%d", count);
  count++;

  testFILE = fopen(str, "r");
  while (testFILE != NULL)
  {
    // close the file & delete it
    fclose(testFILE);
    remove(str);

    // try to open next output file
    size = 1 + snprintf(NULL, 0, "output_%d", count);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "output_%d", count);
    count++;

    testFILE = fopen(str, "r");
  }

  // free the memory
  free(str);

  return; 
}

void remove_temp_files()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  remove("arr.out");
  remove("compareI");
  remove("compareR");
  remove("config");
  remove("const.out");
  remove("eval.out");
  remove("eval2.out");
  remove("finalFile.out");
  remove("gmon.out");
  remove("jacP.out");
  remove("jacV.out");
  remove("num.out");
  remove("paramDerivs.out");
  remove("par.out");
  remove("func_input");  
  remove("nonhom_start");
  remove("preproc_data");
  remove("deg.out");
  remove("names.out");
  remove("witness_data_old"); 
 
  return;
} 

void zero_dim_main(int MPType, double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  // call the corresponding function based on MPType
  if ((MPType == 0) || (MPType == 2)) // double or AMP 
    zero_dim_main_d(MPType, parse_time, currentSeed, startName, my_id, num_processes, headnode);
  else // fixed MP
    zero_dim_main_mp(parse_time, currentSeed, startName, my_id, num_processes, headnode);

  return;
}

// *** initialize_mpi must be called before and MPI routines can be used ***
void initialize_mpi(int argc, char *args[], int *num_of_processes, int *my_id)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes MPI                                        *
\***************************************************************/
{
#ifdef _HAVE_MPI
  // initialize MPI
  MPI_Init(&argc, &args);

  // determine total number of processes and the processor's id number
  MPI_Comm_size(MPI_COMM_WORLD, num_of_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, my_id);
#else
  *num_of_processes = 1;
  *my_id = 0;
#endif

  return;
}

// *** finalize_mpi must be called before the end of the program, but after all MPI routines have been called ***
void finalize_mpi(void)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finalizes MPI                                          *
\***************************************************************/
{
#ifdef _HAVE_MPI
  MPI_Finalize();
#endif

  return;
}


