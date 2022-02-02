// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

int function_eval_main(int printJacobian, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *FuncOut = NULL, *JvOut = NULL, *JpOut = NULL, *StartPts = NULL, *TimeIn = NULL;
  tracker_config_t T;
  preproc_data PPD;
  prog_t dummyProg;
  int i, j, k, num_variables = 0, num_var_gps = 0, userHom = 0, pathMod = 0, paramHom = 0;
  int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, timeIn = 0, useParameters = 0;
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  trackingStats trackCount;
  char start_time[] = "start_time";
  double midpoint_tol, intrinsicCutoffMultiplier;
  point_data_d startPt_d;
  point_data_mp startPt_mp;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  init_trackingStats(&trackCount); // initialize trackCount to all 0

  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType);

  // setup useParameters
  if (userHom || paramHom == 2)
    useParameters = 1;

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // setup a SLP
  T.numVars = setupProg_count(&dummyProg, T.Precision, T.MPType, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow);

  // setup preProcData
  setupPreProcData("preproc_data", &PPD);
  num_var_gps = PPD.num_var_gp;

  init_point_data_d(&startPt_d, T.numVars);
  init_point_data_mp(&startPt_mp, T.numVars);
  startPt_d.point->size = startPt_mp.point->size = T.numVars;

  if (T.MPType == 2)  //If we are doing adaptive precision path-tracking, we must set up AMP_eps, AMP_Phi, AMP_Psi based on config settings.
  { 
    T.AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T.AMP_Phi = T.AMP_bound_on_degree*(T.AMP_bound_on_degree-1.0)*T.AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
    T.AMP_Psi = T.AMP_bound_on_degree*T.AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
    // initialize latest_newton_residual_mp to the maximum precision
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);
  }
  else if (T.MPType == 1)
  { // initialize latest_newton_residual_mp
    mpf_init(T.latest_newton_residual_mp);
  }

#ifdef _HAVE_MPI
  if (num_processes > 1)
  { // using MPI - tell the workers what they will be doing
    worker_info sendType;
    sendType.dataType = STOPCODE;
    bcast_worker_info(&sendType, my_id, headnode);
  }
#endif

  // open StartPts
  StartPts = fopen(startName, "r");
  if (StartPts == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n", startName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // open FuncOut & JvOut
  FuncOut = fopen("function", "w");
  if (printJacobian)
    JvOut = fopen("Jv", "w");

  if (useParameters)
  { // setup JpOut, if needed
    if (printJacobian)
      JpOut = fopen("Jp", "w");

    // see if the user has t-values
    TimeIn = fopen(start_time, "r");
    timeIn = !(TimeIn == NULL);

    if (timeIn)
    { // read in the number of t-values
      fscanf(TimeIn, "%d", &j);
    }
  }

  // read in the number of start points
  fscanf(StartPts, "%d", &trackCount.numPoints);

  // verify positive number
  if (trackCount.numPoints <= 0)
  {
    printf("\n\nERROR: The number of startpoints must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // verify number of points if using user-hom with t-values
  if (useParameters && timeIn)
  { // compare
    if (trackCount.numPoints > j)
    { // not enough t-values!
      printf("\n\nERROR: '%s' needs to contain at least %d time value%s!\n", start_time, trackCount.numPoints, trackCount.numPoints == 1 ? "" : "s");
      bexit(ERROR_INVALID_SIZE);
    }
  }

  // print the number of start points
  fprintf(FuncOut, "%d\n\n", trackCount.numPoints);
  if (printJacobian)
  {
    fprintf(JvOut, "%d\n\n", trackCount.numPoints);
    if (useParameters)
      fprintf(JpOut, "%d\n\n", trackCount.numPoints);
  }

  for (i = 0; i < trackCount.numPoints; i++)
  { // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Evaluating point %d of %d\n", i, trackCount.numPoints);

    // read in the point, evaluate it and print the data
    if (T.MPType == 1)
    { // use fixed precision
      for (j = 0; j < T.numVars; j++)
        if (j < num_var_gps)
        { // set hom variable to 1
          set_one_mp(&startPt_mp.point->coord[j]);
        }
        else
        { // read in from file
          mpf_inp_str(startPt_mp.point->coord[j].r, StartPts, 10);
          mpf_inp_str(startPt_mp.point->coord[j].i, StartPts, 10);
          scanRestOfLine(StartPts);
        }	

      // setup time
      if (useParameters && timeIn)
      { // read in time
        mpf_inp_str(startPt_mp.time->r, TimeIn, 10);
        mpf_inp_str(startPt_mp.time->i, TimeIn, 10);
        scanRestOfLine(TimeIn);
      }
      else
      { // set time to zero
        set_zero_mp(startPt_mp.time);
      }

      // evaluate
      evalProg_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, startPt_mp.point, startPt_mp.time, &dummyProg);

      // print function to FuncOut
      for (j = 0; j < e_mp.funcVals->size; j++)
      {
        print_mp(FuncOut, 0, &e_mp.funcVals->coord[j]);
        fprintf(FuncOut, "\n");
      }
 
      if (printJacobian)
      { // print Jacobian to JvOut
        for (k = 0; k < e_mp.Jv->rows; k++)
        { // print kth row
          for (j = 0; j < T.numVars; j++)
            if (j >= num_var_gps)
            { // print jth column
              print_mp(JvOut, 0, &e_mp.Jv->entry[k][j]);
              fprintf(JvOut, "\n");
            }
        }

        // print Jacobian w.r.t. parameters to JpOut
        if (useParameters)
        {
          for (k = 0; k < e_mp.Jp->rows; k++)
          { // print kth row
            for (j = 0; j < e_mp.Jp->cols; j++)
            { // print jth column
              print_mp(JpOut, 0, &e_mp.Jp->entry[k][j]);
              fprintf(JpOut, "\n");
            }
          }
        }
      }
    }
    else
    { // use double precision
      for (j = 0; j < T.numVars; j++)
        if (j < num_var_gps)
        { // set hom variable to 1
          set_one_d(&startPt_d.point->coord[j]);
        }
        else
        { // read in from file
          fscanf(StartPts, "%lf%lf", &startPt_d.point->coord[j].r, &startPt_d.point->coord[j].i);
          scanRestOfLine(StartPts);
        }

      // setup time
      if (useParameters && timeIn)
      { // read in time
        fscanf(TimeIn, "%lf%lf", &startPt_d.time->r, &startPt_d.time->i);
        scanRestOfLine(TimeIn);
      }
      else
      { // set time to zero
        set_zero_d(startPt_d.time);
      }

      // evaluate
      evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, startPt_d.point, startPt_d.time, &dummyProg);

      // print function to FuncOut
      for (j = 0; j < e_d.funcVals->size; j++)
      {
        print_d(FuncOut, 0, &e_d.funcVals->coord[j]);
        fprintf(FuncOut, "\n");
      }

      if (printJacobian)
      { // print Jacobian to JvOut
        for (k = 0; k < e_d.Jv->rows; k++)
        { // print kth row
          for (j = 0; j < T.numVars; j++)
            if (j >= num_var_gps)
            { // print jth column
              print_d(JvOut, 0, &e_d.Jv->entry[k][j]);
              fprintf(JvOut, "\n");
            }
        }

        // print Jacobian w.r.t. parameters to JpOut
        if (useParameters)
        {
          for (k = 0; k < e_d.Jp->rows; k++)
          { // print kth row
            for (j = 0; j < e_d.Jp->cols; j++)
            { // print jth column
              print_d(JpOut, 0, &e_d.Jp->entry[k][j]);
              fprintf(JpOut, "\n");
            }
          } 
        }
      }
    }
    fprintf(FuncOut, "\n");

    if (printJacobian)
    {
      fprintf(JvOut, "\n");
      if (useParameters)
        fprintf(JpOut, "\n");
    }
  }

  printf("\n------------------------------------------------------------------------------------------\n");
  printf("The following files may be of interest to you:\n\n");
  printf("function:  The function values evaluated at the given points.\n");
  if (printJacobian)
  {
    printf("Jv:        The Jacobian w.r.t. variables evaluated at the given points.\n");

    if (useParameters)
    {
      if (timeIn)
        printf("Jp:        The Jacobian w.r.t. parameters evaluated at the given points and time.\n");
      else
        printf("Jp:        The Jacobian w.r.t. parameters evaluated at the given points at time 0.\n");
    }
  }
  printf("------------------------------------------------------------------------------------------\n");

  // close files
  fclose(StartPts);
  fclose(FuncOut);
  if (JvOut != NULL)
    fclose(JvOut);
  if (TimeIn != NULL)
    fclose(TimeIn);
  if (JpOut != NULL)
    fclose(JpOut);

  // clear memory
  clear_point_data_d(&startPt_d);
  clear_point_data_mp(&startPt_mp);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  if (!userHom)
  {
    if (dummyProg.numSubfuncs > 0)
    { // clear subFuncsBelow
      for (i = dummyProg.numFuncs - 1; i >= 0; i--)
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

  clearProg(&dummyProg, T.MPType, 1);

  tracker_config_clear(&T);
  clearMP();
 
  return 0;
}



