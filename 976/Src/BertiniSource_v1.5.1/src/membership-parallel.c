// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

void setupMembershipRandomSlicePoint(witness_t *W, int pure_codim_index, point_d inputPt_d, point_mp inputPt_mp, int inputPt_prec, int MPType, int max_prec, point_d endPt_d, point_mp endPt_mp, int *endPt_prec);
void printIncidenceMatrix(int ***incidence, int num_points, int num_codim, witnessCodim_t *codim);

////////// membership tests /////////////////

void membershipTest(witness_t *W, tracker_config_t *T, int pathMod, char *testPoints);
int pureDimMembershipTest_seq(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, double midpoint_tol, tracker_config_t *T, FILE *OUT, char *midName);
void setupMembershipRandomSlice(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, int MPType, int max_prec, point_d endPt_d, point_mp endPt_mp, int *endPt_prec);
int membershipTrack(int *isMember, witness_t *W, int pathNum_startPt, int pure_codim_index, point_d endPt_d, point_mp endPt_mp, int endPt_prec, double tol, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void setupMembershipEndpoint_d(point_d endPt, point_d origEndpoint, int num_var_gps, vec_d patch, vec_d H, comp_d homVarConst);
void setupMembershipEndpoint_mp(point_mp endPt, point_mp origEndpoint, int num_var_gps, vec_mp patch, vec_mp H, comp_mp homVarConst);
void membershipTesting(witness_t *W, int num_points, point_d *testPts_d, point_mp *testPts_mp, int pt_prec, int pathMod, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);


void membershipMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: performs membership test                               *
\***************************************************************/
{
  int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  tracker_config_t T;
  witness_t witnessSet;

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

  // quit if we are using MPI 
  if (num_processes > 1)
  { // exit since parallel sampling is not implemented
#ifdef _HAVE_MPI
    printf("ERROR: Parallel membership test is not implemented. Please use sequential version!\n");
    bexit(ERROR_OTHER);
#endif
  }

  // setup witnessSet
  setupWitnessDataFromFile("witness_data", "witness_data_old", "preproc_data", "deg.out", &witnessSet, &T, 1);

  // error checking
  if (witnessSet.Prog->numPathVars > 0)
  { // path variable present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (witnessSet.Prog->numPars > 0)
  { // parameter present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup the number of variables
  T.numVars = witnessSet.Prog->numVars;

  // setup the rest of T, if needed
  if (T.MPType == 1)
  { // initialize latest_newton_residual_mp
    mpf_init(T.latest_newton_residual_mp);
  }
  else if (T.MPType == 2)
  { // initialize latest_newton_residual_mp
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);

    // setup eps, Phi & Psi
    T.AMP_eps = (double) T.numVars * T.numVars;
    T.AMP_Phi = T.AMP_bound_on_degree * (T.AMP_bound_on_degree - 1) * T.AMP_bound_on_abs_vals_of_coeffs;
    T.AMP_Psi = T.AMP_bound_on_degree * T.AMP_bound_on_abs_vals_of_coeffs;
  }

  // create output files now that everything is setup
  numIrredDecompOutput(&witnessSet, &T, 3, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 3

  // loop through the test points to test to see if and which component they are on
  membershipTest(&witnessSet, &T, pathMod, "member_points");

  // display decomposition chart
  numIrredDecompChart(&witnessSet, stdout, T.MPType, reducedOnly);

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void membershipTest(witness_t *W, tracker_config_t *T, int pathMod, char *testPoints)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: read in the test points and check which component it is*
\***************************************************************/
{
  int i, j, rV, num_points = 0, prec = T->Precision, base = 10;
  point_d *testPoints_d = NULL;
  point_mp *testPoints_mp = NULL;
  FILE *OUT = fopen("output_membership", "w"), *MIDOUT = fopen("midpath_data", "w"), *IN = fopen(testPoints, "r");

  // check that IN exists
  if (IN == NULL)
  { // file does not exist
    printf("\n\nERROR: '%s' does not exist!!!\n\n", testPoints);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in the number of points that are to be tested
  fscanf(IN, "%d", &num_points);
  scanRestOfLine(IN);

  // setup the precision to use
  if (T->MPType == 0)
  { // use D
    prec = 52;
  }
  else if (T->MPType == 1)
  { // use MP
    prec = T->Precision;
  }
  else
  { // use enough precision based on the final tolerance
    prec = -floor(log10(T->final_tolerance) - 0.5);
    prec = digits_to_prec(prec);
  }

  if (prec < 64)
  { // setup testPoints_d
    testPoints_d = (point_d *)bmalloc(num_points * sizeof(point_d));

    // loop through to read in the test points to _d
    for (i = 0; i < num_points; i++)
    { // setup the size
      init_point_d(testPoints_d[i], W->orig_variables);
      testPoints_d[i]->size = W->orig_variables;

      if (W->PPD.num_var_gp)
      { // Bertini homogenized the variable group       
        set_one_d(&testPoints_d[i]->coord[0]); // set hom coord == 1
        for (j = 1; j < W->orig_variables; j++)
        {
          rV = fscanf(IN, "%lf%lf", &testPoints_d[i]->coord[j].r, &testPoints_d[i]->coord[j].i);
          if (rV < 0)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", testPoints);
            bexit(ERROR_INVALID_SIZE);
          }
          // scan in rest of line
          rV = scanRestOfLine(IN);
          if (rV && j + 1 < W->orig_variables)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", testPoints);
            bexit(ERROR_INVALID_SIZE);
          }
        }
      }
      else
      { // already homogenized variable group
        for (j = 0; j < W->orig_variables; j++)
        {
          rV = fscanf(IN, "%lf%lf", &testPoints_d[i]->coord[j].r, &testPoints_d[i]->coord[j].i);
          if (rV < 0)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", testPoints);
            bexit(ERROR_INVALID_SIZE);
          }
          // scan in rest of line
          rV = scanRestOfLine(IN);
          if (rV && j + 1 < W->orig_variables)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", testPoints);
            bexit(ERROR_INVALID_SIZE);
          }
        }
      }
    }
  }
  else 
  { // setup testPoints_mp
    testPoints_mp = (point_mp *)bmalloc(num_points * sizeof(point_mp));

    // loop through to read in the test points to _mp
    for (i = 0; i < num_points; i++)
    { // initialize the point
      init_point_mp2(testPoints_mp[i], W->orig_variables, prec);
      // setup the size
      testPoints_mp[i]->size = W->orig_variables;

      if (W->PPD.num_var_gp)
      { // Bertini homogenized the variable group       
        set_one_mp(&testPoints_mp[i]->coord[0]); // set hom coord == 1
        for (j = 1; j < W->orig_variables; j++)
        {
          mpf_inp_str(testPoints_mp[i]->coord[j].r, IN, base);
          mpf_inp_str(testPoints_mp[i]->coord[j].i, IN, base);
          // scan rest of line
          scanRestOfLine(IN);
        }
      }
      else
      { // already homogenized variable group
        for (j = 0; j < W->orig_variables; j++)
        {
          mpf_inp_str(testPoints_mp[i]->coord[j].r, IN, base);
          mpf_inp_str(testPoints_mp[i]->coord[j].i, IN, base);
          // scan rest of line
          scanRestOfLine(IN);
        }
      }
    }
  }

  // do the membership testing
  membershipTesting(W, num_points, testPoints_d, testPoints_mp, prec, pathMod, T, OUT, MIDOUT);

  // close files
  fclose(OUT);
  fclose(MIDOUT);

  // clear testPoints
  if (prec < 64)
  { // clear _d
    for (i = num_points - 1; i >= 0; i--)
      clear_point_d(testPoints_d[i]);
    free(testPoints_d);
  }
  else
  { // clear testPoints_mp
    for (i = num_points - 1; i >= 0; i--)
      clear_point_mp(testPoints_mp[i]);
    free(testPoints_mp);
  }

  return;
}

void membershipTesting(witness_t *W, int num_points, point_d *testPts_d, point_mp *testPts_mp, int pt_prec, int pathMod, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do the actual membership testing for each test point   *
\***************************************************************/
{
  int i, j, k, retVal, setupGamma, isMemberOfSomething, isMemberOfComp, badCount, maxBadCount = 20, init_prec = MAX(64, pt_prec), num_codim = W->num_codim;
  double funcRes_d;
  mpf_t funcRes_mp;
  mat_d tempMat_d;
  mat_mp tempMat_mp;
  point_data_d tempPD_d;
  point_data_mp tempPD_mp;
  int **fullRankProgInfo = NULL; // == 1 if fullRankProgs needs to be cleared, otherwise it will be NULLed out
  membership_slice_moving_t *sliceMover = NULL;
  prog_t ***fullRankProgs = NULL;
  endpoint_data_d **endPts_d = NULL;
  endpoint_data_mp **endPts_mp = NULL;
  endpoint_data_amp **endPts_amp = NULL;

  // allocate & initialize the incidence matrix
  int ***incidence = (int ***)bmalloc(num_points * sizeof(int **));
  for (i = 0; i < num_points; i++)
  {
    incidence[i] = (int **)bmalloc(num_codim * sizeof(int *));
    for (j = 0; j < num_codim; j++)
    {
      incidence[i][j] = (int *)bmalloc(W->codim[j].num_components * sizeof(int));
      for (k = 0; k < W->codim[j].num_components; k++)
        incidence[i][j][k] = 0;
    }
  }
  
  mpf_init2(funcRes_mp, init_prec);
  init_mat_d(tempMat_d, 0, 0);
  init_mat_mp2(tempMat_mp, 0, 0, init_prec);
  init_point_data_d(&tempPD_d, 0);
  init_point_data_mp2(&tempPD_mp, 0, init_prec);

  // allocate for each codimension
  sliceMover = (membership_slice_moving_t *)bmalloc(num_codim * sizeof(membership_slice_moving_t));
  fullRankProgs = (prog_t ***)bmalloc(num_codim * sizeof(prog_t **));
  fullRankProgInfo = (int **)bmalloc(num_codim * sizeof(int *));
  if (T->MPType == 0)
    endPts_d = (endpoint_data_d **)bmalloc(num_codim * sizeof(endpoint_data_d *));
  else if (T->MPType == 1)
    endPts_mp = (endpoint_data_mp **)bmalloc(num_codim * sizeof(endpoint_data_mp *));
  else
    endPts_amp = (endpoint_data_amp **)bmalloc(num_codim * sizeof(endpoint_data_amp *));

  // setup sliceMover, fullRankProgs & endPts for each codimension
  for (i = 0; i < num_codim; i++)
  {
    basic_setup_slice_moving(&sliceMover[i], W, i, T->MPType, T->AMP_max_prec);
    if (T->MPType == 0)
      deflate_for_junkRemoval(&fullRankProgs[i], &fullRankProgInfo[i], &endPts_d[i], NULL, NULL, &sliceMover[i], W, i, -1, T, OUT);
    else if (T->MPType == 1)
      deflate_for_junkRemoval(&fullRankProgs[i], &fullRankProgInfo[i], NULL, &endPts_mp[i], NULL, &sliceMover[i], W, i, -1, T, OUT);
    else
      deflate_for_junkRemoval(&fullRankProgs[i], &fullRankProgInfo[i], NULL, NULL, &endPts_amp[i], &sliceMover[i], W, i, -1, T, OUT);
  }

  // display message
  printf("\nTesting membership: %d point%s to test.\n", num_points, num_points == 1 ? "" : "s");

  for (i = 0; i < num_points; i++)
  { // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Testing %d of %d\n", i, num_points);

    fprintf(OUT, "\n*****************************************************\n");
    fprintf(OUT, "Testing %d\n", i);
    fprintf(OUT, "*****************************************************\n");

    // verify that this point actually satisfies the original system
    retVal = 0;
    if (pt_prec < 64)
    { // use testPts_d
      point_cp_d(tempPD_d.point, testPts_d[i]);
      set_zero_d(tempPD_d.time);
      findFunctionResidual_d(&funcRes_d, &tempPD_d, W->Prog, &evalProg_d_void);

      // verify that the function residual is small enough
      if (funcRes_d > T->final_tol_times_mult)
      { // print error message
        printf("\nIt appears that point %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", i, funcRes_d, T->final_tol_times_mult);
        retVal = 1;
      }
    }
    else
    { // use testPts_mp
      if (T->MPType == 2)
      { // set the precision properly - all other things have this precision already set
        initMP(pt_prec);
        change_witness_prec(W, pt_prec);
      }
      point_cp_mp(tempPD_mp.point, testPts_mp[i]);
      set_zero_mp(tempPD_mp.time);
      findFunctionResidual_mp(funcRes_mp, &tempPD_mp, W->Prog, &evalProg_mp_void);

      if (mpf_get_d(funcRes_mp) > T->final_tol_times_mult)
      { // print error message
        printf("\nIt appears that point %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", i, mpf_get_d(funcRes_mp), T->final_tol_times_mult);
        retVal = 1;
      }
    }

    // if this point satisfies the original system, perform the membership test
    if (!retVal)
    { // perform the membership test by looping over the different codimensions
      printf("\n");
      isMemberOfSomething = 0;
      for (j = 0; j < num_codim; j++)
      { // perform membership test with the pure-dim witness set

        // move the point to the proper patch and setup a random slice through the point
        if (pt_prec < 64)
        { // setup _d
          change_size_mat_d(tempMat_d, 1, sliceMover[j].p_d->size);
          tempMat_d->rows = 1;
          tempMat_d->cols = sliceMover[j].p_d->size;
          for (k = 0; k < tempMat_d->cols; k++)
          {
            set_d(&tempMat_d->entry[0][k], &sliceMover[j].p_d->coord[k]);
          }
          // move the point to the patch
          move_to_patch_mat_d(testPts_d[i], testPts_d[i], tempMat_d, &W->PPD);

          // generate a random slice through the test point
          setup_slice_moving_slice(&sliceMover[j], testPts_d[i], NULL, pt_prec, T->MPType, T->AMP_max_prec);
        }
        else
        { // setup _mp
          change_size_mat_mp(tempMat_mp, 1, sliceMover[j].p_mp->size);
          tempMat_mp->rows = 1;
          tempMat_mp->cols = sliceMover[j].p_mp->size;
          for (k = 0; k < tempMat_mp->cols; k++)
          {
            set_mp(&tempMat_mp->entry[0][k], &sliceMover[j].p_mp->coord[k]);
          }
          // move the point to the patch
          move_to_patch_mat_mp(testPts_mp[i], testPts_mp[i], tempMat_mp, &W->PPD);

          // generate a random slice through the test point
          setup_slice_moving_slice(&sliceMover[j], NULL, testPts_mp[i], pt_prec, T->MPType, T->AMP_max_prec);
        }

        // move all of the points for this codim to the random slice and see if any hit endPt
        do
        { // we need to setup gamma
          setupGamma = 1;
          // initialize badCount & retVal
          badCount = retVal = 0;

          // loop through the witness points and move the slice to see if any of the endpoints hit the testPt
          // continue until we have checked all of them successfully OR we have an error
          for (k = 0; k < W->codim[j].num_set && !retVal; k++)
          { // initialize
            retVal = isMemberOfComp = 0;

            // verify that deflation worked
            if (fullRankProgInfo[j][k] != -1)
            { // finish the setup for sliceMover for this endpoint
              final_setup_slice_moving(&sliceMover[j], fullRankProgs[j][k], T->MPType, T->AMP_max_prec, setupGamma);

              // we only setup gamma once for the whole set of witness points to avoid monodromy
              if (setupGamma)
                setupGamma = 0;

              // track the path and check for equality
              if (T->MPType == 0)
              { // compare using _d
                retVal = membership_slice_moving_track(&isMemberOfComp, testPts_d[i], NULL, pt_prec, &sliceMover[j], endPts_d[j][k].endPt, NULL, 52, k, T->final_tol_times_mult, T, OUT, MIDOUT);
              }
              else if (T->MPType == 1)
              { // compare using _mp
                retVal = membership_slice_moving_track(&isMemberOfComp, NULL, testPts_mp[i], pt_prec, &sliceMover[j], NULL, endPts_mp[j][k].endPt, T->Precision, k, T->final_tol_times_mult, T, OUT, MIDOUT);
              }
              else
              { // compare using _amp
                if (pt_prec < 64)
                  retVal = membership_slice_moving_track(&isMemberOfComp, testPts_d[i], NULL, pt_prec, &sliceMover[j], endPts_amp[j][k].endPt_d, endPts_amp[j][k].endPt_mp, endPts_amp[j][k].curr_prec, k, T->final_tol_times_mult, T, OUT, MIDOUT);
                else
                  retVal = membership_slice_moving_track(&isMemberOfComp, NULL, testPts_mp[i], pt_prec, &sliceMover[j], endPts_amp[j][k].endPt_d, endPts_amp[j][k].endPt_mp, endPts_amp[j][k].curr_prec, k, T->final_tol_times_mult, T, OUT, MIDOUT);
              }
            }

            if (retVal)
            { // we had bad random numbers - increment the number of bad loops and try again
              badCount++;
  
              if (badCount >= maxBadCount)
              { // exit since we have failure - this should never happen
                printf("\nERROR: Failed to move slice (after trying %d different choices of random numbers).\n", maxBadCount);
                bexit(ERROR_LOOP_FAILURE);
              }
            }
            else if (isMemberOfComp)
            { // print the membership
              printf("It appears that point %d lies on component %d of dimension %d.\n", i, W->codim[j].component_nums[k], W->orig_variables - W->codim[j].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp);
              isMemberOfSomething = 1;
              
              // update incidence matrix
              incidence[i][j][W->codim[j].component_nums[k]] = 1;
            }
          }
        } while (retVal);
      }
      // check to see if we have found which component this point is a member of 
      if (!isMemberOfSomething)  
        printf("\nThe membership test for point %d has failed.  There are several possible reasons:\n  (1) the point actually lies on a component that was not found by Bertini, or\n  (2) the tracking tolerances are too restrictive.\nTo correct for (2), please try looser tracking tolerances.\n", i);
    }
  }
  printf("\n");

  // construct incidence matrix file
  printIncidenceMatrix(incidence, num_points, num_codim, W->codim);

  // clear incidence matrix
  for (i = num_points - 1; i >= 0; i--)
  {
    for (j = num_codim - 1; j >= 0; j--)
      free(incidence[i][j]);
    free(incidence[i]);
  }
  free(incidence);

  // clear memory
  mpf_clear(funcRes_mp);
  clear_mat_d(tempMat_d);
  clear_mat_mp(tempMat_mp);
  clear_point_data_d(&tempPD_d);
  clear_point_data_mp(&tempPD_mp); 

  // clear the deflation information
  clear_sliceMover_fullRankProgs(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, W, T->MPType);

  return;
}

void printIncidenceMatrix(int ***incidence, int num_points, int num_codim, witnessCodim_t *codim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                * 
* NOTES: prints the incidence matrix to a file                  *
\***************************************************************/
{
  int i, j, k;
  FILE *FP = fopen("incidence_matrix","w");

  // print the number of codimensions
  fprintf(FP, "%d\n", num_codim);
  // print the codim and number of components in each codimension
  for (i = 0; i < num_codim; i++)
    fprintf(FP, "%d %d\n", codim[i].codim, codim[i].num_components);

  // print the number of points
  fprintf(FP, "\n%d\n\n", num_points);

  // loop over the points printing the incidence vector for each point
  for (i = 0; i < num_points; i++)
  {
    for (j = 0; j < num_codim; j++)
      for (k = 0; k < codim[j].num_components; k++)
        fprintf(FP, "%d ", incidence[i][j][k]);
    fprintf(FP, "\n"); 
  }
  fprintf(FP, "\n");

  fclose(FP);

  return;
}

int pureDimMembershipTest(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, double midpoint_tol, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - the endpoint given by pathNum is a member  *
* of the pure_codim_index pure-dim witness set, 0 - otherwise   *
* NOTES: determines if the endpoint is a member of algebraic set*
\***************************************************************/
{
  int isMember = 1;

  if (num_processes > 1)
  {
    printf("ERROR: Parallel membership test is not implemented. Please use sequential version!\n");
    bexit(ERROR_OTHER);
  }
  else
  {
    isMember = pureDimMembershipTest_seq(W, pure_codim_index, pathNum, pathNum_codim_index, midpoint_tol, T, OUT, midName);
  }

  return isMember;
}

int pureDimMembershipTest_seq(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, double midpoint_tol, tracker_config_t *T, FILE *OUT, char *midName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - the endpoint given by pathNum is a member  *
* of the pure_codim_index pure-dim witness set, 0 - otherwise   *
* NOTES: perform membership test in sequential fashion          *
\***************************************************************/
{
  int i, j, endPt_prec, retVal = 0, badCount = 0, isMember = 0, num_paths = W->codim[pure_codim_index].num_nonsing;
  int *pathNum_startPts = (int *)bmalloc(num_paths * sizeof(int));
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_prec)(void const *, int);
  FILE *MIDOUT = fopen(midName, "w");
  point_d endPt_d;
  point_mp endPt_mp;

  // initialize endPt
  init_point_d(endPt_d, 0);
  init_point_mp(endPt_mp, 0);

  // setup the evaluators
  ptr_to_eval_d = &membership_slice_moving_eval_d;
  ptr_to_eval_mp = &membership_slice_moving_eval_mp;
  change_prec = &change_witness_prec;

  // setup pathNum_startPts
  j = 0;
  for (i = 0; i < num_paths; i++)
  { // find the next non-singular endpoint
    while (W->codim[pure_codim_index].witnessPt_types[j] != NON_SINGULAR)
      j++;

    // store the endpoint
    pathNum_startPts[i] = j;

    // move past
    j++;
  }

  // generate a random slice through the 'pathNum' endpoint with dimension that of pure_codim_index
  setupMembershipRandomSlice(W, pure_codim_index, pathNum, pathNum_codim_index, T->MPType, T->AMP_max_prec, endPt_d, endPt_mp, &endPt_prec);
  
  // loop through the nonsingular witness points and move them to the random slice to see if any of the endpoints are equal to the 'pathNum' endpoint
  for (i = 0; i < num_paths && !isMember; i++)
  { // reset the bad counter
    badCount = 0;

    do 
    { // track the path and check for equality
      retVal = membershipTrack(&isMember, W, pathNum_startPts[i], pure_codim_index, endPt_d, endPt_mp, endPt_prec, T->final_tol_times_mult, T, OUT, MIDOUT, ptr_to_eval_d, ptr_to_eval_mp, change_prec); 

      if (retVal)
      { // we had bad random numbers
        badCount++;

        if (badCount < 20)
        { // generate new random slice
          setupMembershipRandomSlice(W, pure_codim_index, pathNum, pathNum_codim_index, T->MPType, T->AMP_max_prec, endPt_d, endPt_mp, &endPt_prec);
        }
        else
        { // exit since we have failure - this should never happen
          printf("\nERROR: Failed to eliminate junk correctly (after trying 20 different choices of random numbers).\n");
          bexit(ERROR_LOOP_FAILURE);
        }
      }
    } while (retVal);
  }

  // close MIDOUT
  fclose(MIDOUT);

  // free pathNum_startPts
  free(pathNum_startPts);

  // clear endPt
  clear_point_d(endPt_d);
  clear_point_mp(endPt_mp);

  return isMember; 
}

int membershipTrack(int *isMember, witness_t *W, int pathNum_startPt, int pure_codim_index, point_d endPt_d, point_mp endPt_mp, int endPt_prec, double tol, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the path in codim_index & path_num of CD        *
\***************************************************************/
{
  int retVal;
  endgame_data_t endPt;
  init_endgame_data(&endPt, T->Precision);
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = &zero_dehom;

  // setup for tracking
  T->first_step_of_path = 1;
  T->endgameOnly = 0;
  W->curr_codim_index = pure_codim_index;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;
    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    point_cp_d(startPt.point, W->codim[pure_codim_index].witnessPts_d[pathNum_startPt].endPt);
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, W, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(pathNum_startPt, &endPt, &startPt, OUT, MIDOUT, T, W, W, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // print the footer for the path and find retVal
    retVal = printMembershipFooter_d(&endPt.PD_d, endPt.condition_number, endPt.function_residual_d, endPt.latest_newton_residual_d, endPt.t_val_at_latest_sample_point_d, endPt.error_at_latest_sample_point_d, OUT, endPt.retVal, T);

    if (!retVal)
    { // find the difference between endpoint of the path and the given endPt
      vec_sub_d(endPt.PD_d.point, endPt.PD_d.point, endPt_d);
      if (infNormVec_d(endPt.PD_d.point) < tol)
      { // they are the same!
        *isMember = 1;
      }
    }

    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;
    init_point_data_mp(&startPt, 0);

    // setup for tracking the path
    point_cp_mp(startPt.point, W->codim[pure_codim_index].witnessPts_mp[pathNum_startPt].endPt);
    set_one_mp(startPt.time);

    // print the header for the path
    printPathHeader_mp(OUT, &startPt, T, pathNum_startPt, W, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_mp(pathNum_startPt, &endPt, &startPt, OUT, MIDOUT, T, W, ptr_to_eval_mp, find_dehom);

    // print the footer in multi precision and find retVal
    retVal = printMembershipFooter_mp(&endPt.PD_mp, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, endPt.retVal, T);

    if (!retVal)
    { // find the difference between endpoint of the path and the given endPt
      vec_sub_mp(endPt.PD_mp.point, endPt.PD_mp.point, endPt_mp);
      if (infNormVec_mp(endPt.PD_mp.point) < tol)
      { // they are the same!
        *isMember = 1;
      }
    }

    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    mpf_t norm;
    point_data_d startPt;

    mpf_init(norm);
    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    if (W->codim[pure_codim_index].witnessPts_amp[pathNum_startPt].curr_prec < 64)
    { // copy to startPt.point in double precision
      point_cp_d(startPt.point, W->codim[pure_codim_index].witnessPts_amp[pathNum_startPt].endPt_d);
    }
    else
    { // covert to double precision and copy to startPt.point
      point_mp_to_d(startPt.point, W->codim[pure_codim_index].witnessPts_amp[pathNum_startPt].endPt_mp);
    }
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, W, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_d(pathNum_startPt, &endPt, &startPt, OUT, MIDOUT, T, W, W, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // print the footer for the path and find retVal
    if (endPt.prec < 64)
    { // print the footer in double precision and find retVal
      retVal = printMembershipFooter_d(&endPt.PD_d, endPt.condition_number, endPt.function_residual_d, endPt.latest_newton_residual_d, endPt.t_val_at_latest_sample_point_d, endPt.error_at_latest_sample_point_d, OUT, endPt.retVal, T);
    }
    else
    { // print the footer in multi precision and find retVal
      retVal = printMembershipFooter_mp(&endPt.PD_mp, endPt.condition_number, endPt.first_increase, endPt.function_residual_mp, endPt.latest_newton_residual_mp, endPt.t_val_at_latest_sample_point_mp, endPt.error_at_latest_sample_point_mp, OUT, endPt.retVal, T);
    }

    if (!retVal)
    { // find the difference between endpoint of the path and the given endPt

      findDiff_point(norm, endPt.PD_d.point, endPt.PD_mp.point, endPt.prec, endPt_d, endPt_mp, endPt_prec);
      if (mpf_get_d(norm) < tol)
      { // they are the same!
        *isMember = 1;
      }
    }

    mpf_clear(norm);
    clear_point_data_d(&startPt);
  }

  clear_endgame_data(&endPt);

  return 0;
}

int printMembershipFooter_d(point_data_d *endPoint, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, int retVal_in, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints a basic footer for membership test              *
\***************************************************************/
{
  int i, retVal_out = 0, isNumber = 1;

  if (retVal_in && d_abs_d(endPoint->time) < T->minTrackT)
  { // we have come close enough, consider it a success

    // display a warning message
    if (T->screenOut)
      printf("WARNING: Path reached the minimum value of T (%e < %e) with retVal %d.\n", endPoint->time->r, T->minTrackT, retVal_in);
    fprintf(OUT, "WARNING: Path reached the minimum value of T (%e < %e) with retVal %d.\n", endPoint->time->r, T->minTrackT, retVal_in);

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

  // print the path footer to OUT
  printPathFooterOut_d(OUT, NULL, 0, 0, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, NULL, T, NULL, 0, 0);

  return retVal_out;
}

int printMembershipFooter_mp(point_data_mp *endPoint, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, int retVal_in, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: correct return value for the path              *
* NOTES: prints a basic footer for membership test              *
\***************************************************************/
{
  int i, retVal_out = 0, isNumber = 1;

  if (retVal_in && d_abs_mp(endPoint->time) < T->minTrackT)
  { // we have come close enough, consider it a success

    // display a warning message
    if (T->screenOut)
      printf("WARNING: Path reached the minimum value of T (%e < %e) with retVal %d.\n", mpf_get_d(endPoint->time->r), T->minTrackT, retVal_in);
    fprintf(OUT, "WARNING: Path reached the minimum value of T (%e < %e) with retVal %d.\n", mpf_get_d(endPoint->time->r), T->minTrackT, retVal_in);

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

  // print the path footer to OUT
  printPathFooterOut_mp(OUT, NULL, 0, 0, endPoint, cond_num, func_residual, newton_error, t_val_sample, error_sample, first_increase, NULL, T, NULL, 0, 0);

  return retVal_out;
}

void setupMembershipRandomSlice(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, int MPType, int max_prec, point_d endPt_d, point_mp endPt_mp, int *endPt_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the target slice in W based on the info          *
\***************************************************************/
{
  if (MPType == 0)
    setupMembershipRandomSlicePoint(W, pure_codim_index, W->codim[pathNum_codim_index].witnessPts_d[pathNum].endPt, NULL, 52, MPType, max_prec, endPt_d, endPt_mp, endPt_prec);
  else if (MPType == 1)
    setupMembershipRandomSlicePoint(W, pure_codim_index, NULL, W->codim[pathNum_codim_index].witnessPts_mp[pathNum].endPt, W->curr_precision, MPType, max_prec, endPt_d, endPt_mp, endPt_prec);
  else
    setupMembershipRandomSlicePoint(W, pure_codim_index, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_d, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_mp, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].curr_prec, MPType, max_prec, endPt_d, endPt_mp, endPt_prec);

  return;
}

void setupMembershipEndpoint_d(point_d endPt, point_d origEndpoint, int num_var_gps, vec_d patch, vec_d H, comp_d homVarConst)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the endPt on the correct patch                    *
\***************************************************************/
{
  int i, size;
  comp_d tempComp, tempComp2;

  size = origEndpoint->size;
  increase_size_point_d(endPt, size);
  endPt->size = size;

  if (num_var_gps)
  { // we want endPt * patch = 1
    set_zero_d(tempComp);
    for (i = 0; i < size; i++)
    {
      sum_mul_d(tempComp, &origEndpoint->coord[i], &patch->coord[i]);
    }
    // reciprocate
    recip_d(tempComp, tempComp);
    // normalize
    for (i = 0; i < size; i++)
    {
      mul_d(&endPt->coord[i], &origEndpoint->coord[i], tempComp);
    }
  }
  else
  { // we want endPt * (patch - H) = homVarConst
    set_zero_d(tempComp);
    for (i = 0; i < size; i++)
    {
      sub_d(tempComp2, &patch->coord[i], &H->coord[i]);
      sum_mul_d(tempComp, &origEndpoint->coord[i], tempComp2);
    }
    // divide
    div_d(tempComp, homVarConst, tempComp);
    // normalize
    for (i = 0; i < size; i++)
    {
      mul_d(&endPt->coord[i], &origEndpoint->coord[i], tempComp);
    }
  }

  return;
}

void setupMembershipEndpoint_mp(point_mp endPt, point_mp origEndpoint, int num_var_gps, vec_mp patch, vec_mp H, comp_mp homVarConst)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the endPt on the correct patch                    *
\***************************************************************/
{
  int i, size;
  comp_mp tempComp, tempComp2;
  init_mp(tempComp); init_mp(tempComp2);

  size = origEndpoint->size;
  increase_size_point_mp(endPt, size);
  endPt->size = size;

  if (num_var_gps)
  { // we want endPt * patch = 1
    set_zero_mp(tempComp);
    for (i = 0; i < size; i++)
    {
      sum_mul_mp(tempComp, &origEndpoint->coord[i], &patch->coord[i]);
    }
    // reciprocate
    recip_mp(tempComp, tempComp);
    // normalize
    for (i = 0; i < size; i++)
    {
      mul_mp(&endPt->coord[i], &origEndpoint->coord[i], tempComp);
    }
  }
  else
  { // we want endPt * (patch - H) = homVarConst
    set_zero_mp(tempComp);
    for (i = 0; i < size; i++)
    {
      sub_mp(tempComp2, &patch->coord[i], &H->coord[i]);
      sum_mul_mp(tempComp, &origEndpoint->coord[i], tempComp2);
    }
    // divide
    div_mp(tempComp, homVarConst, tempComp);
    // normalize
    for (i = 0; i < size; i++)
    {
      mul_mp(&endPt->coord[i], &origEndpoint->coord[i], tempComp);
    }
  }

  clear_mp(tempComp); clear_mp(tempComp2);

  return;
}

void setupMembershipRandomSlicePoint(witness_t *W, int pure_codim_index, point_d inputPt_d, point_mp inputPt_mp, int inputPt_prec, int MPType, int max_prec, point_d endPt_d, point_mp endPt_mp, int *endPt_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: endPt - point that we are trying to hit        *
* NOTES: use the given point to setup the target slice in W     *
\***************************************************************/
{
  int i, j, rows, cols;

  if (MPType == 0)
  { // setup extrinsic slice in double precision
    comp_d tempComp, hom_y;
    point_d tempPoint;
    init_point_d(tempPoint, 0);

    // setup tempPoint
    if (inputPt_prec < 64)
    { // copy inputPt_d to tempPoint
      point_cp_d(tempPoint, inputPt_d);
    }
    else
    { // copy inputPt_mp to tempPoint
      point_mp_to_d(tempPoint, inputPt_mp);
    }

    // setup endPt_d - this is what we want to try and hit (move tempPoint to the correct patch)
    *endPt_prec = 52;
    setupMembershipEndpoint_d(endPt_d, tempPoint, W->PPD.num_var_gp, W->codim[pure_codim_index].p_d, W->codim[pure_codim_index].H_d, W->codim[pure_codim_index].homVarConst_d);

    // find hom_y = H * endPt + homVarConst
    set_d(hom_y, W->codim[pure_codim_index].homVarConst_d);
    for (i = 0; i < endPt_d->size; i++)
    {
      sum_mul_d(hom_y, &endPt_d->coord[i], &W->codim[pure_codim_index].H_d->coord[i]);
    }

    // see if targetSliceMat & targetSliceVec need initialized
    if (!W->targetSliceInit)
    { // initialize
      init_mat_d(W->targetSliceMat_d, 0, 0);
      init_vec_d(W->targetSliceVec_d, 0);
      W->targetSliceInit = 1;
    }

    // generate a random matrix
    rows = W->codim[pure_codim_index].B_d->rows;
    cols = W->codim[pure_codim_index].B_d->cols;
    make_matrix_random_d(W->targetSliceMat_d, rows, cols);

    // find matrix * endPt
    mul_mat_vec_d(W->targetSliceVec_d, W->targetSliceMat_d, endPt_d);

    // update matrix = hom_y * matrix - (matrix * endPt)*H
    // update vector = vector * homVarConst
    for (i = 0; i < rows; i++)
    {
      for (j = 0; j < cols; j++)
      {
        mul_d(tempComp, &W->targetSliceVec_d->coord[i], &W->codim[pure_codim_index].H_d->coord[j]);
        mul_d(&W->targetSliceMat_d->entry[i][j], &W->targetSliceMat_d->entry[i][j], hom_y);
        sub_d(&W->targetSliceMat_d->entry[i][j], &W->targetSliceMat_d->entry[i][j], tempComp);
      }

      // update vector = vector * homVarConst
      mul_d(&W->targetSliceVec_d->coord[i], &W->targetSliceVec_d->coord[i], W->codim[pure_codim_index].homVarConst_d);
    }

    clear_point_d(tempPoint);
  }
  else if (MPType == 1)
  { // setup extrinsic slice in fixed multi precision
    comp_mp tempComp, hom_y;
    point_mp tempPoint;

    init_mp(tempComp); init_mp(hom_y);
    init_point_mp(tempPoint, 0);

    // setup tempPoint
    if (inputPt_prec < 64)
    { // copy inputPt_d to tempPoint
      point_d_to_mp(tempPoint, inputPt_d);
    }
    else
    { // copy inputPt_mp to tempPoint
      point_cp_mp(tempPoint, inputPt_mp);
    }

    // setup endPt_mp - this is what we want to try and hit (move tempPoint to the correct patch)
    *endPt_prec = W->curr_precision;
    setprec_point_mp(endPt_mp, *endPt_prec);
    setupMembershipEndpoint_mp(endPt_mp, tempPoint, W->PPD.num_var_gp, W->codim[pure_codim_index].p_mp, W->codim[pure_codim_index].H_mp, W->codim[pure_codim_index].homVarConst_mp);

    // find hom_y = H * endPt + homVarConst
    set_mp(hom_y, W->codim[pure_codim_index].homVarConst_mp);
    for (i = 0; i < endPt_mp->size; i++)
    {
      sum_mul_mp(hom_y, &endPt_mp->coord[i], &W->codim[pure_codim_index].H_mp->coord[i]);
    }

    // see if targetSliceMat & targetSliceVec need initialized
    if (!W->targetSliceInit)
    { // initialize
      init_mat_mp(W->targetSliceMat_mp, 0, 0);
      init_vec_mp(W->targetSliceVec_mp, 0);
      W->targetSliceInit = 1;
    }

    // generate a random matrix
    rows = W->codim[pure_codim_index].B_mp->rows;
    cols = W->codim[pure_codim_index].B_mp->cols;
    make_matrix_random_mp(W->targetSliceMat_mp, rows, cols, W->curr_precision);

    // find matrix * endPt
    mul_mat_vec_mp(W->targetSliceVec_mp, W->targetSliceMat_mp, endPt_mp);

    // update matrix = hom_y * matrix - (matrix*endPt)*H
    // update vector = vector * homVarConst
    for (i = 0; i < rows; i++)
    {
      for (j = 0; j < cols; j++)
      {
        mul_mp(tempComp, &W->targetSliceVec_mp->coord[i], &W->codim[pure_codim_index].H_mp->coord[j]);
        mul_mp(&W->targetSliceMat_mp->entry[i][j], &W->targetSliceMat_mp->entry[i][j], hom_y);
        sub_mp(&W->targetSliceMat_mp->entry[i][j], &W->targetSliceMat_mp->entry[i][j], tempComp);
      }

      // update vector = vector * homVarConst
      mul_mp(&W->targetSliceVec_mp->coord[i], &W->targetSliceVec_mp->coord[i], W->codim[pure_codim_index].homVarConst_mp);
    }

    clear_mp(tempComp); clear_mp(hom_y);
    clear_point_mp(tempPoint);
  }
  else
  { // setup extrinsic slice in adaptive precision
    int point_prec = MAX(64, inputPt_prec);
    mpq_t tempMPQ[2];
    comp_mp tempComp, hom_y;
    vec_mp tempVec_mp;

    mpq_init(tempMPQ[0]); mpq_init(tempMPQ[1]);
    init_mp2(tempComp, point_prec); init_mp2(hom_y, point_prec);
    init_vec_mp2(tempVec_mp, 0, point_prec);

    // setup endPt - this is what we want to try and hit (move inputPt to the correct patch)
    *endPt_prec = inputPt_prec;
    if (*endPt_prec < 64)
    { // setup endPt_d & tempVec_mp using inputPt_d
      setupMembershipEndpoint_d(endPt_d, inputPt_d, W->PPD.num_var_gp, W->codim[pure_codim_index].p_d, W->codim[pure_codim_index].H_d, W->codim[pure_codim_index].homVarConst_d);

      point_d_to_mp(tempVec_mp, endPt_d);
    }
    else
    { // setup endPt_mp & tempVec_mp using inputPt_mp
      setprec_point_mp(endPt_mp, *endPt_prec);
      setupMembershipEndpoint_mp(endPt_mp, inputPt_mp, W->PPD.num_var_gp, W->codim[pure_codim_index].p_mp, W->codim[pure_codim_index].H_mp, W->codim[pure_codim_index].homVarConst_mp);

      point_cp_mp(tempVec_mp, endPt_mp);
    }

    // find hom_y = H * tempVec + homVarConst
    mpf_set_q(hom_y->r, W->codim[pure_codim_index].homVarConst_rat[0]);
    mpf_set_q(hom_y->i, W->codim[pure_codim_index].homVarConst_rat[1]);
    for (i = 0; i < tempVec_mp->size; i++)
    {
      mpf_set_q(tempComp->r, W->codim[pure_codim_index].H_rat[i][0]);
      mpf_set_q(tempComp->i, W->codim[pure_codim_index].H_rat[i][1]);

      sum_mul_mp(hom_y, &tempVec_mp->coord[i], tempComp);
    }

    // find number of rows and columns for random matrix
    rows = W->codim[pure_codim_index].B_d->rows;
    cols = W->codim[pure_codim_index].B_d->cols;

    // if target slices are already intialized, clear them out
    if (W->targetSliceInit)
    { // clear Mat
      clear_mat(W->targetSliceMat_d, W->targetSliceMat_mp, W->targetSliceMat_rat, MPType);
      // clear Vec
      clear_vec(W->targetSliceVec_d, W->targetSliceVec_mp, W->targetSliceVec_rat, MPType);
    }

    // initialize Mat correctly
    init_mat_d(W->targetSliceMat_d, rows, cols);
    init_mat_mp2(W->targetSliceMat_mp, rows, cols, W->curr_precision);
    init_mat_rat(W->targetSliceMat_rat, rows, cols);

    // initialize Vec correctly
    init_vec_d(W->targetSliceVec_d, rows); 
    init_vec_mp2(W->targetSliceVec_mp, rows, W->curr_precision);
    init_vec_rat(W->targetSliceVec_rat, rows);

    W->targetSliceInit = 1;

    // generate a random matrix
    make_matrix_random_rat(W->targetSliceMat_d, W->targetSliceMat_mp, W->targetSliceMat_rat, rows, cols, W->curr_precision, max_prec, 0, 0);

    // find matrix * tempVec
    W->targetSliceVec_d->size = W->targetSliceVec_mp->size = rows;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&W->targetSliceVec_mp->coord[i]);
      for (j = 0; j < cols; j++)
      {
        mpf_set_q(tempComp->r, W->targetSliceMat_rat[i][j][0]);
        mpf_set_q(tempComp->i, W->targetSliceMat_rat[i][j][1]);

        sum_mul_mp(&W->targetSliceVec_mp->coord[i], tempComp, &tempVec_mp->coord[j]);
      }
      // convert to a rational number
      mpf_t_to_rat(W->targetSliceVec_rat[i][0], W->targetSliceVec_mp->coord[i].r);
      mpf_t_to_rat(W->targetSliceVec_rat[i][1], W->targetSliceVec_mp->coord[i].i);
    }

    // update matrix = hom_y * matrix - (matrix*tempVec)*H
    // update vector = vector * homVarConst
    for (i = 0; i < rows; i++)
    {
      for (j = 0; j < cols; j++)
      {
        mpf_set_q(tempComp->r, W->codim[pure_codim_index].H_rat[j][0]);
        mpf_set_q(tempComp->i, W->codim[pure_codim_index].H_rat[j][1]);

        mul_mp(tempComp, &W->targetSliceVec_mp->coord[i], tempComp);
        mul_mp(&W->targetSliceMat_mp->entry[i][j], &W->targetSliceMat_mp->entry[i][j], hom_y);
        sub_mp(&W->targetSliceMat_mp->entry[i][j], &W->targetSliceMat_mp->entry[i][j], tempComp);

        // convert to a rational number
        mpf_t_to_rat(W->targetSliceMat_rat[i][j][0], W->targetSliceMat_mp->entry[i][j].r);
        mpf_t_to_rat(W->targetSliceMat_rat[i][j][1], W->targetSliceMat_mp->entry[i][j].i);

        // convert _mp & _d
        mpf_set_q(W->targetSliceMat_mp->entry[i][j].r, W->targetSliceMat_rat[i][j][0]);
        mpf_set_q(W->targetSliceMat_mp->entry[i][j].i, W->targetSliceMat_rat[i][j][1]);
        W->targetSliceMat_d->entry[i][j].r = mpq_get_d(W->targetSliceMat_rat[i][j][0]);
        W->targetSliceMat_d->entry[i][j].i = mpq_get_d(W->targetSliceMat_rat[i][j][1]);
      }

      // update vector = vector * homVarConst
      mpq_mul(tempMPQ[0], W->targetSliceVec_rat[i][0], W->codim[pure_codim_index].homVarConst_rat[0]);
      mpq_mul(tempMPQ[1], W->targetSliceVec_rat[i][1], W->codim[pure_codim_index].homVarConst_rat[1]);
      mpq_sub(tempMPQ[0], tempMPQ[0], tempMPQ[1]);
      mpq_mul(tempMPQ[1], W->targetSliceVec_rat[i][0], W->codim[pure_codim_index].homVarConst_rat[1]);
      mpq_mul(W->targetSliceVec_rat[i][1], W->targetSliceVec_rat[i][1], W->codim[pure_codim_index].homVarConst_rat[0]);
      mpq_add(W->targetSliceVec_rat[i][1], W->targetSliceVec_rat[i][1], tempMPQ[1]);
      mpq_set(W->targetSliceVec_rat[i][0], tempMPQ[0]);

      // convert _mp & _d
      mpf_set_q(W->targetSliceVec_mp->coord[i].r, W->targetSliceVec_rat[i][0]);
      mpf_set_q(W->targetSliceVec_mp->coord[i].i, W->targetSliceVec_rat[i][1]);
      W->targetSliceVec_d->coord[i].r = mpq_get_d(W->targetSliceVec_rat[i][0]);
      W->targetSliceVec_d->coord[i].i = mpq_get_d(W->targetSliceVec_rat[i][1]);
    }

    mpq_clear(tempMPQ[0]); mpq_clear(tempMPQ[1]);
    clear_mp(tempComp); clear_mp(hom_y);
    clear_vec_mp(tempVec_mp);
  }

  return;
}

