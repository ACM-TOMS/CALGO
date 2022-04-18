// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

int proj_newton_residual_d(vec_d dX, double *residual, int return_norm, double *CN, vec_d P, comp_d time, int **degrees, preproc_data *PPD, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int proj_newton_residual_mp(vec_mp dX, mpf_t residual, int return_norm, double *CN, vec_mp P, comp_mp time, int **degrees, preproc_data *PPD, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int newton_eval_main(int computeCN, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *RVOut = NULL, *NewOut = NULL, *ResOut = NULL, *CNOut = NULL, *StartPts = NULL, *TimeIn = NULL;
  tracker_config_t T;
  preproc_data PPD;
  prog_t dummyProg;
  int i, j, rV, num_variables = 0, num_var_gps = 0, userHom = 0, pathMod = 0, paramHom = 0, useParameters = 0;
  int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, timeIn = 0;
  int **degrees = NULL;
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  trackingStats trackCount;
  char start_time[] = "start_time";
  double midpoint_tol, intrinsicCutoffMultiplier, norm_J, norm_J_inv, CN;
  point_data_d startPt_d;
  point_data_mp startPt_mp;
  vec_d orig_d, dX_d;
  vec_mp orig_mp, dX_mp;
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

  init_vec_d(dX_d, 0);
  init_vec_mp(dX_mp, 0);

  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // setup a SLP
  T.numVars = setupProg_count(&dummyProg, T.Precision, T.MPType, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow);

  // setup preProcData
  setupPreProcData("preproc_data", &PPD);
  num_var_gps = PPD.num_var_gp;

  // read in degrees
  TimeIn = fopen("deg.out", "r");
  if (TimeIn == NULL)
  {
    printf("\nERROR: The file 'deg.out' does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  degrees = (int **)bmalloc(dummyProg.numFuncs * sizeof(int *));
  rV = PPD.num_var_gp + PPD.num_hom_var_gp;
  for (i = 0; i < dummyProg.numFuncs; i++)
  { // allocate
    degrees[i] = (int *)bmalloc(rV * sizeof(int));
    for (j = 0; j < rV; j++)
      fscanf(TimeIn, "%d\n", &degrees[i][j]);
  }
  fclose(TimeIn);
  TimeIn = NULL;

  // verify square system or overdetermined system
  if (T.numVars > dummyProg.numFuncs + PPD.num_var_gp + PPD.num_hom_var_gp)
  {
    printf("\nERROR: This operation assumes that the system is square or overdetermined!\n");
    bexit(ERROR_CONFIGURATION);
  }

  init_vec_d(orig_d, T.numVars - num_var_gps);
  init_vec_mp(orig_mp, T.numVars - num_var_gps);
  orig_d->size = orig_mp->size = T.numVars - num_var_gps;

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

  // open RVOut, NewOut, ResOut & CNOut
  RVOut = fopen("successFlag", "w");
  NewOut = fopen("newPoints", "w");
  ResOut = fopen("newtonResidual", "w");
  if (computeCN)
    CNOut = fopen("CN", "w");

  if (useParameters)
  { // see if the user has t-values
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
  fprintf(RVOut, "%d\n\n", trackCount.numPoints);
  fprintf(NewOut, "%d\n\n", trackCount.numPoints);
  fprintf(ResOut, "%d\n\n", trackCount.numPoints);
  if (computeCN)
    fprintf(CNOut, "%d\n\n", trackCount.numPoints);

  for (i = 0; i < trackCount.numPoints; i++)
  { // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Evaluating point %d of %d\n", i, trackCount.numPoints);

    // read in the point, compute newton iteration, and print the data
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
          set_mp(&orig_mp->coord[j - num_var_gps], &startPt_mp.point->coord[j]);
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

      if (userHom == 1)
      { // compute usual newton iteration
        rV = newton_residual_mp(dX_mp, T.latest_newton_residual_mp, computeCN, &norm_J, &norm_J_inv, startPt_mp.point, startPt_mp.time, &e_mp, &dummyProg, &evalProg_mp_void);

        // print rV to RVOut
        fprintf(RVOut, "%d\n", rV);

        if (!rV)
        { // successful
          if (computeCN)
            CN = norm_J * norm_J_inv;

          // print residual to ResOut
          mpf_out_str(ResOut, 10, 0, T.latest_newton_residual_mp);
          fprintf(ResOut, "\n");

          // update start point
          vec_sub_mp(startPt_mp.point, startPt_mp.point, dX_mp);
        }
        else
        { // not successful
          CN = -1;

          // print residual to ResOut
          fprintf(ResOut, "0\n");
        }

        if (computeCN)
        { // print CN
          if (CN > 0)
            fprintf(CNOut, "%.15e\n", CN);
          else
            fprintf(CNOut, "%d\n", (int) CN);
        }

        // print point to NewOut
        for (j = 0; j < T.numVars; j++)
          if (j >= num_var_gps)
          { // print jth column
            print_mp(NewOut, 0, &startPt_mp.point->coord[j]);
            fprintf(NewOut, "\n");
          }
      }
      else
      { // compute m-hom newton iteratio
        rV = proj_newton_residual_mp(dX_mp, T.latest_newton_residual_mp, computeCN, &CN, startPt_mp.point, startPt_mp.time, degrees, &PPD, &e_mp, &dummyProg, &evalProg_mp_void);

        // print rV to RVOut
        fprintf(RVOut, "%d\n", rV);

        if (!rV)
        { // successful - compute new point
          vec_sub_mp(e_mp.parVals, startPt_mp.point, dX_mp);

          // dehomogenize
          getDehomPoint_mp(e_mp.funcVals, e_mp.parVals, T.numVars, &PPD);

          // print point to NewOut
          for (j = num_var_gps; j < T.numVars; j++)
          { // print jth column
            print_mp(NewOut, 0, &e_mp.funcVals->coord[j - num_var_gps]);
            fprintf(NewOut, "\n");
          }

          // compute residual
          vec_sub_mp(e_mp.parVals, e_mp.funcVals, orig_mp);

          // compute residual in 2-norm
          twoNormVec_mp2(e_mp.parVals, T.latest_newton_residual_mp);

          // print residual to ResOut
          mpf_out_str(ResOut, 10, 0, T.latest_newton_residual_mp);
          fprintf(ResOut, "\n");
        }
        else
        { // failure
          CN = -1;

          // print point to NewOut
          for (j = num_var_gps; j < T.numVars; j++)
          { // print jth column
            print_mp(NewOut, 0, &orig_mp->coord[j - num_var_gps]);
            fprintf(NewOut, "\n");
          }

          // print residual to ResOut
          fprintf(ResOut, "0\n");
        }

        if (computeCN)
        { // print CN
          if (CN > 0)
            fprintf(CNOut, "%.15e\n", CN);
          else
            fprintf(CNOut, "%d\n", (int) CN);
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
          set_d(&orig_d->coord[j - num_var_gps], &startPt_d.point->coord[j]);
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

      if (userHom == 1)
      { // compute usual newton iteration
        rV = newton_residual_d(dX_d, &T.latest_newton_residual_d, computeCN, &norm_J, &norm_J_inv, startPt_d.point, startPt_d.time, &e_d, &dummyProg, &evalProg_d_void);

        // print rV to RVOut
        fprintf(RVOut, "%d\n", rV);

        if (!rV)
        { // successful
          if (computeCN)
            CN = norm_J * norm_J_inv;

          // print residual to ResOut
          fprintf(ResOut, "%.15e\n", T.latest_newton_residual_d);

          // update start point
          vec_sub_d(startPt_d.point, startPt_d.point, dX_d);
        }
        else
        { // not successful
          CN = -1;

          // print residual to ResOut
          fprintf(ResOut, "0\n");
        }

        if (computeCN)
        { // print CN
          if (CN > 0)
            fprintf(CNOut, "%.15e\n", CN);
          else
            fprintf(CNOut, "%d\n", (int) CN);
        }

        // print point to NewOut
        for (j = num_var_gps; j < T.numVars; j++)
        { // print jth column
          print_d(NewOut, 0, &startPt_d.point->coord[j]);
          fprintf(NewOut, "\n");
        }
      }
      else
      { // compute m-hom newton iteratio
        rV = proj_newton_residual_d(dX_d, &T.latest_newton_residual_d, computeCN, &CN, startPt_d.point, startPt_d.time, degrees, &PPD, &e_d, &dummyProg, &evalProg_d_void);

        // print rV to RVOut
        fprintf(RVOut, "%d\n", rV);

        if (!rV)
        { // successful - compute new point
          vec_sub_d(e_d.parVals, startPt_d.point, dX_d);

          // dehomogenize
          getDehomPoint_d(e_d.funcVals, e_d.parVals, T.numVars, &PPD);

          // print point to NewOut
          for (j = num_var_gps; j < T.numVars; j++)
          { // print jth column
            print_d(NewOut, 0, &e_d.funcVals->coord[j - num_var_gps]);
            fprintf(NewOut, "\n");
          }

          // compute residual
          vec_sub_d(e_d.parVals, e_d.funcVals, orig_d);

          // compute residual in 2-norm
          twoNormVec_d(e_d.parVals, &T.latest_newton_residual_d); 

          // print residual to ResOut
          fprintf(ResOut, "%.15e\n", T.latest_newton_residual_d);
        }
        else
        { // failure
          CN = -1;

          // print point to NewOut
          for (j = num_var_gps; j < T.numVars; j++)
          { // print jth column
            print_d(NewOut, 0, &orig_d->coord[j - num_var_gps]);
            fprintf(NewOut, "\n");
          }

          // print residual to ResOut
          fprintf(ResOut, "0\n");
        }

        if (computeCN)
        { // print CN
          if (CN > 0)
            fprintf(CNOut, "%.15e\n", CN);
          else
            fprintf(CNOut, "%d\n", (int) CN);
        }
      }
    }

    fprintf(RVOut, "\n");
    fprintf(NewOut, "\n");
    fprintf(ResOut, "\n");
    if (computeCN)
      fprintf(CNOut, "\n");
  }

  printf("\n------------------------------------------------------------------------------------------\n");
  printf("The following files may be of interest to you:\n\n");
  printf("successFlag:    Success or failure of Newton iterations.\n");
  printf("newPoints:      The point resulting from the Newton iteration.\n");
  printf("newtonResidual: The Newton residual for each point.\n");
  if (computeCN)
    printf("CN:             Approximation of the condition number at each point.\n");
  printf("------------------------------------------------------------------------------------------\n");

  // close files
  fclose(StartPts);
  if (TimeIn != NULL)
    fclose(TimeIn);
  fclose(RVOut);
  fclose(NewOut);
  fclose(ResOut);
  if (computeCN)
    fclose(CNOut);

  // clear memory
  clear_point_data_d(&startPt_d);
  clear_point_data_mp(&startPt_mp);
  clear_vec_d(orig_d);
  clear_vec_mp(orig_mp);
  clear_vec_d(dX_d);
  clear_vec_mp(dX_mp);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  for (i = 0; i < dummyProg.numFuncs; i++)
    free(degrees[i]);
  free(degrees);

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

///////////////// PROJECTIVE NEWTON ITERATION WITH CONDITION NUMBER /////////////////

int proj_newton_residual_d(vec_d dX, double *residual, int return_norm, double *CN, vec_d P, comp_d time, int **degrees, preproc_data *PPD, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find CN                  *
* NOTES: finds the newton residual dX but does not update P     *
\***************************************************************/
{
  int i, j, k, retVal, rows, cols, num_gps = PPD->num_hom_var_gp + PPD->num_var_gp, var = PPD->num_var_gp, count = 0;
  double norm_J, norm_J_inv;

  // evaluate the function
  evalProg_d_void(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED);

  // verify size
  if (e->Jv->rows + num_gps < e->Jv->cols)
  {
    printf("ERROR: The system must either be square or overdetermined when performing projective Newton iterations!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // increase the size of funcVals & Jv
  rows = e->Jv->rows;
  cols = e->Jv->cols;
  increase_rows_mat_d(e->Jv, rows + num_gps);
  increase_size_vec_d(e->funcVals, rows + num_gps);
  e->Jv->rows = e->funcVals->size = rows + num_gps;

  // input correct data into extra locations
  for (i = 0; i < num_gps; i++)
  { // set func[k] to 0
    k = i + rows;
    set_zero_d(&e->funcVals->coord[k]);

    // initialize new rows in Jv to 0
    for (j = 0; j < cols; j++)
      set_zero_d(&e->Jv->entry[k][j]);

    // set Jv[i + rows][:] to conjugate of ith group
    if (PPD->type[i] != 0)
    { // added homogeneous variable
      conjugate_d(&e->Jv->entry[k][count], &P->coord[count]);
      count++;
    }
    // other variables in group
    for (j = 0; j < PPD->size[i]; j++)
    {
      conjugate_d(&e->Jv->entry[k][var], &P->coord[var]);
      var++;
    }
  }

  if (return_norm)
  { // compute Newton residual and its condition number
    retVal = matrixSolve_cond_num_norms_d(dX, e->Jv, e->funcVals, CN, &norm_J, &norm_J_inv);
  }
  else
  { // just compute newton residual
    retVal = matrixSolve_d(dX, e->Jv, e->funcVals);
  }

  // check for matrixSolve failure
  *residual = 0;
  if (!retVal)
  { // find the size of dX
    *residual = infNormVec_d(dX);
  }

  return retVal;
}

int proj_newton_residual_mp(vec_mp dX, mpf_t residual, int return_norm, double *CN, vec_mp P, comp_mp time, int **degrees, preproc_data *PPD, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if return_norm != 0 - find CN                  *
* NOTES: finds the newton residual dX but does not update P     *
\***************************************************************/
{
  int i, j, k, retVal, rows, cols, num_gps = PPD->num_hom_var_gp + PPD->num_var_gp, var = PPD->num_var_gp, count = 0;
  double norm_J, norm_J_inv;

  // evaluate the function
  evalProg_mp_void(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P, time, ED);

  // verify size
  if (e->Jv->rows + num_gps < e->Jv->cols)
  {
    printf("ERROR: The system must either be square or overdetermined when performing projective Newton iterations!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // increase the size of funcVals & Jv
  rows = e->Jv->rows;
  cols = e->Jv->cols;
  increase_rows_mat_mp(e->Jv, rows + num_gps);
  increase_size_vec_mp(e->funcVals, rows + num_gps);
  e->Jv->rows = e->funcVals->size = rows + num_gps;

  // input correct data into extra locations
  for (i = 0; i < num_gps; i++)
  { // set func[k] to 0
    k = i + rows;
    set_zero_mp(&e->funcVals->coord[k]);

    // initialize new rows in Jv to 0
    for (j = 0; j < cols; j++)
      set_zero_mp(&e->Jv->entry[k][j]);

    // set Jv[i + rows][:] to conjugate of ith group
    if (PPD->type[i] != 0)
    { // added homogeneous variable
      conjugate_mp(&e->Jv->entry[k][count], &P->coord[count]);
      count++;
    }
    // other variables in group
    for (j = 0; j < PPD->size[i]; j++)
    {
      conjugate_mp(&e->Jv->entry[k][var], &P->coord[var]);
      var++;
    }
  }

  if (return_norm)
  { // compute Newton residual and its condition number
    retVal = matrixSolve_cond_num_norms_mp(dX, e->Jv, e->funcVals, CN, &norm_J, &norm_J_inv);
  }
  else
  { // just compute newton residual
    retVal = matrixSolve_mp(dX, e->Jv, e->funcVals);
  }

  // check for matrixSolve failure
  mpf_set_ui(residual, 0);
  if (!retVal)
  { // find the size of dX
    infNormVec_mp2(residual, dX);
  }

  return retVal;
}



