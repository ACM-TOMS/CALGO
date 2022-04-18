// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"

// This file contains the sharpening module.

//////////////////////////// REDO THE POST PROCESSING WITH SHARPENING THE ENDPOINTS //////////////////////////////////////

void sharpen_process_zero_dim_main(int MPType, unsigned int currentSeed);

void sharpening_menu(int pathMod, int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void sharpen_all_endpoints(int pathMod, int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void sharpen_endpoints_file(int pathMod, int userHom, int useRegen, char *fileName, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void sharpen_endpoints_manually(int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int sharpen_endpoint(post_process_t *endPoint, tracker_config_t *T, FILE *OUT, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int sharpen_endpoint_using_endgame(post_process_t *endPoint, int pathNum, tracker_config_t *T, FILE *OUT, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void countFiniteSuccess(int *finiteSuccess, int *infiniteSuccess, int *nonfiniteSuccess, int numStartPts, post_process_t *endPoints, int num_path_nums, int *path_nums);


void sharpen_process_main(int MPType, int trackType, unsigned int currentSeed, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: main sharpening function when called directly          *
\***************************************************************/
{
  if (trackType == 0)
  { // use zero dimensional sharpener
    sharpen_process_zero_dim_main(MPType, currentSeed);
  }
  else
  { // use the sampling menu to find 'sharp' sample points on the selected components
    sampleComponent(currentSeed, MPType, 1, my_id, num_processes, headnode);
  }

  return;
}

void sharpen_process_zero_dim_main(int MPType, unsigned int currentSeed)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: main sharpening function                               *
\***************************************************************/
{
  FILE *OUT = NULL, *midIN = NULL, *midOUT = NULL, *rawIN = NULL, *rawOUT = NULL;
  tracker_config_t T;
  prog_t dummyProg;
  int i, j, size, count, max_prec, retVal = 0, num_variables = 0, userHom = 0, paramHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, eqbyeqMethod = 0;
  int midOutExists = 0;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;  
  basic_eval_data_d  ED_d;
  basic_eval_data_mp ED_mp;
  double midpoint_tol, intrinsicCutoffMultiplier;
  post_process_t *endPoints = NULL;

  // check for existence of 'raw_data'
  rawIN = fopen("raw_data", "r");
  if (rawIN == NULL)
  {
    printf("\nERROR: 'raw_data' does not exist!  This file is needed to sharpen the endpoints!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  // now that we have its exitence, move it!
  fclose(rawIN);
  rename("raw_data", "raw_data_old");
  rawIN = fopen("raw_data_old", "r");
  rawOUT = fopen("raw_data", "w");

  // check for existence of 'midpath_data'
  midIN = fopen("midpath_data", "r");
  if (midIN == NULL)
  {
    midOutExists = 0;
  }
  else
  { // its exitence, move it!
    midOutExists = 1;
    fclose(midIN);
    rename("midpath_data", "midpath_data_old");
    midIN = fopen("midpath_data_old", "r");
  }

  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType); // Set up T using the config file.
  reproduceInputFile("this_input", "func_input", &T, 0, 0, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  if (userHom == -59 || useRegen)
    eqbyeqMethod = 1;

  // call the setup function based on these parameters
  if (userHom == 1)
  { // setup for user defined homotopy
    if (T.MPType == 0 || T.MPType == 2)
      num_variables = userHom_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED_d, &dummyProg, &ptr_to_eval_d, &ptr_to_eval_mp);
    else
      num_variables = userHom_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED_mp, &dummyProg, &ptr_to_eval_mp);
  }
  else if (paramHom == 2 || userHom == 2)
  { // setup for parameter homotopy
    if (T.MPType == 0 || T.MPType == 2)
      num_variables = paramHom_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED_d, &dummyProg, &ptr_to_eval_d, &ptr_to_eval_mp, "preproc_data", 0, "start"); // do not setup start points
    else
      num_variables = paramHom_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED_mp, &dummyProg, &ptr_to_eval_mp, "preproc_data", 0, "start"); // do not setup start points
  }
  else if (userHom == 0 || userHom == -59)
  { // setup for standard tracking
    int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
    if (T.MPType == 0 || T.MPType == 2)
    {
      num_variables = zero_dim_basic_setup_d(&OUT, "output", &midOUT, "midpath_data", &T, &ED_d, &dummyProg, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow, &ptr_to_eval_d, &ptr_to_eval_mp, "preproc_data", "deg.out", 0, NULL, NULL);

      // clear subFuncsBelow
      if (ED_d.squareSystem.Prog->numSubfuncs > 0)
      { // clear subFuncsBelow
        for (i = ED_d.squareSystem.Prog->numFuncs - 1; i >= 0; i--)
          free(subFuncsBelow[i]);
        free(subFuncsBelow);
      }
    }
    else
    {
      num_variables = zero_dim_basic_setup_mp(&OUT, "output", &midOUT, "midpath_data", &T, &ED_mp, &dummyProg, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow, &ptr_to_eval_mp, "preproc_data", "deg.out", 0, NULL, NULL);

      // clear subFuncsBelow
      if (ED_mp.squareSystem.Prog->numSubfuncs > 0)
      { // clear subFuncsBelow
        for (i = ED_mp.squareSystem.Prog->numFuncs - 1; i >= 0; i--)
          free(subFuncsBelow[i]);
        free(subFuncsBelow);
      }
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
  else // SHOULD NEVER HAPPEN!
  { // put raw_data_old back at raw_data, midpath_data_old back at midpath_data, and exit!
    fclose(rawIN);
    fclose(rawOUT);
    fclose(midIN);
    rename("raw_data_old", "raw_data");
    if (midOutExists)
      rename("midpath_data_old", "midpath_data");
    printf("ERROR: The sharpening module needs to have userHom (%d) == 0, 1, or 2!\n", userHom);
    bexit(ERROR_CONFIGURATION);
  }

  // setup change_prec - regeneration structures are ignored since we just use the original functions
  change_prec = &change_basic_eval_prec;

  if (T.MPType == 2)  //If we are doing adaptive precision path-tracking, we must set up AMP_eps, AMP_Phi, AMP_Psi based on config settings.
  {
    T.AMP_eps = (double) num_variables * num_variables;  // According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T.AMP_Phi = T.AMP_bound_on_degree * (T.AMP_bound_on_degree - 1.0) * T.AMP_bound_on_abs_vals_of_coeffs;  // Phi from the AMP paper.
    T.AMP_Psi = T.AMP_bound_on_degree * T.AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
    // initialize latest_newton_residual_mp to the maximum precision
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);
  }
  else if (T.MPType == 1)
  { // initialize latest_newton_residual_mp to the fixed precision
    mpf_init2(T.latest_newton_residual_mp, T.Precision);
  }

  // now that everything is setup, we need to read in the data from rawIN
  count = size = -1;
  fscanf(rawIN, "%d\n%d\n", &count, &size); // read in the number of variables and dimension

  // error checking
  if (count != num_variables)
  { // put raw_data_old back at raw_data, midpath_data_old back at midpath_data, and exit!
    fclose(rawIN);
    fclose(rawOUT);
    fclose(midIN);
    fclose(midOUT);
    rename("raw_data_old", "raw_data");
    if (midOutExists)
      rename("midpath_data_old", "midpath_data");
    printf("ERROR: The number of variables described in 'raw_data' is not correct!\n");
    bexit(ERROR_CONFIGURATION);
  } 
  if (size != 0)
  { // put raw_data_old back at raw_data, midpath_data_old back at midpath_data, and exit!
    fclose(rawIN);
    fclose(rawOUT);
    fclose(midIN);
    fclose(midOUT);
    rename("raw_data_old", "raw_data");
    if (midOutExists)
      rename("midpath_data_old", "midpath_data");
    printf("ERROR: The dimension described in 'raw_data' is not correct!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  count = retVal = 0;
  size = 1;
  endPoints = (post_process_t *)brealloc(endPoints, size * sizeof(post_process_t));

  // loop until we have reached the bottom of rawIN
  max_prec = 0;
  while (!retVal)
  { // make sure we have room to store the next one
    if (size <= count)
    {
      endPoints = (post_process_t *)brealloc(endPoints, 2 * size * sizeof(post_process_t));
      size *= 2;
    }
    retVal = setupPostProcess(&j, rawIN, &endPoints[count], num_variables, T.MPType);
    
    if (!retVal)
    { // check to see if this is the maximum precision used
      if (j > max_prec)
        max_prec = j;
      // update the number read in
      count++;
    }
  }
  // store the number of endpoints that was read in
  size = count;

  // check for errors
  if ((T.MPType == 0 && max_prec > 52) || (T.MPType == 1 && max_prec > T.Precision))
  { 
    printf("NOTE: The original run used higher precision than currently utilized and so the data will be truncated to the current precision.\n");
  }

  // check to make sure that there is something in endPoints
  if (size == 0)
  { // no endpoints!!!! - put raw_data_old back at raw_data, midpath_data_old back at midpath_data, and exit!
    fclose(rawIN);
    fclose(rawOUT);
    fclose(midIN);
    rename("raw_data_old", "raw_data");
    if (midOutExists)
      rename("midpath_data_old", "midpath_data");
    printf("NOTE: There are no endpoints in 'raw_data' to sharpen!\n");
  }
  else
  { // setup the patch with the old values stored in rawIN
    if (T.MPType == 0)
      updateZeroDimRelevantData(&ED_d, NULL, T.MPType, eqbyeqMethod, rawIN);
    else if (T.MPType == 1) 
      updateZeroDimRelevantData(NULL, &ED_mp, T.MPType, eqbyeqMethod, rawIN);
    else
      updateZeroDimRelevantData(&ED_d, ED_d.BED_mp, T.MPType, eqbyeqMethod, rawIN);

     // go to the menu to see what the user wants to do
    if (T.MPType == 0 || T.MPType == 2)
      sharpening_menu(pathMod, userHom, useRegen, endPoints, size, &T, OUT, rawOUT, midOutExists, midIN, &ED_d, ED_d.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    else
      sharpening_menu(pathMod, userHom, useRegen, endPoints, size, &T, OUT, rawOUT, midOutExists, midIN, NULL, &ED_mp, NULL, ptr_to_eval_mp, NULL);

    // close rawIN, rawOUT
    fclose(rawIN);
    fclose(rawOUT);
  }
  // close other files
  fclose(midOUT);
  fclose(OUT);

  // free the data
  for (i = size - 1; i >= 0; i--)
  {
    if (endPoints[i].sol_prec >= 64)
    {
      mpf_clear(endPoints[i].function_resid_mp);
      mpf_clear(endPoints[i].newton_resid_mp);
      for (j = 0; j < num_variables; j++)
      {
        clear_mp(endPoints[i].sol_mp[j]);
      }
    }
    free(endPoints[i].sol_mp);
    free(endPoints[i].sol_d);
  }
  free(endPoints);

  // clear the memory - 0 for no regeneration - only sharpening the whole system
  if (T.MPType == 0 || T.MPType == 2)
    basic_eval_clear_d(&ED_d, 0, T.MPType);
  else
   basic_eval_clear_mp(&ED_mp, 0, T.MPType);

  clearMP();
  tracker_config_clear(&T);

  // remove temporary files
  remove("raw_data_old");
  if (midOutExists)
    remove("midpath_data_old");
  else
    remove("midpath_data");

  return;
}

int pos_dim_sharpening_menu(tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: displays the menu for positive dim sharpending module  *
\***************************************************************/
{
  int i, rV, option, selection_found = 0;
  char ch;

  do
  { // find the sharpening option
    printf("\nPositive-Dimensional Sharpening options:\n");
    printf("  1. Component sampling\n");
    printf("  2. Change the number of sharpening digits (currently %d digits)\n", T->sharpenDigits);
    printf("  3. Quit\n");
    printf("Please enter your option: ");
    rV = scanf("%d", &option);

    if (rV < 0)
    { // at EOF - so we need to quit
      option = 3;
      selection_found = 1;
    }
    else // we are not at EOF
    { // flush the buffer
      do
      { 
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");
        selection_found = 0;
      }
      else if (option < 1 || option > 3)
      { // invalid option
        printf("\n%d is not a valid option!\n", option);
        selection_found = 0;
      }
      else if (option == 2)
      { // change the number of sharpening digits
        do
        {
          printf("\nPlease input the number of sharpening digits (currently %d digits): ", T->sharpenDigits);
          rV = scanf("%d", &i);

          if (rV < 0)
          { // at EOF - so we need to quit
            option = 3;
            selection_found = 1;
            i = 1;
          }
          else // we are not at EOF
          { // flush the buffer
            do
            {
              ch = getchar();
            } while (ch != EOF && ch != '\n');

            if (rV == 0)
            { // input was not read in correctly - so quit
              printf("\nThe input was not read in correctly so the number of sharpening digits will not be changed!\n");
              i = 1;
            }
            else // scanf was successful
            { // make sure the current precision can handle this many digits
              if (T->MPType == 0)
                rV = 15;
              else if (T->MPType == 1)
                rV = prec_to_digits(T->Precision);
              else
                rV = INT_MAX;

              if (i > rV)
              { // too many digits for this precision
                T->sharpenDigits = rV - 1;
                printf("\nThe current precision cannot handle %d digits reliably so the number of sharpening digits has been set to %d digits!\n", i, T->sharpenDigits);
              }
              else if (i > 0) // && i <= rV
              { // input is valid so store to sharpenDigits
                T->sharpenDigits = i;
                printf("\nThe number of sharpening digits has been changed to %d digits!\n", T->sharpenDigits);
              }
              else // i <= 0
              { // input not valid positive number
                printf("\nThe number of sharpening digits must be positive!\n");
              }
            }
            selection_found = 0;
          }
        } while (i <= 0);
      }
      else // valid selection
        selection_found = 1;
    }
  } while (!selection_found);

  if (option == 1)
  { // continue on to sampling menu
    rV = 1; 
  }
  else 
  { // do not continue on to sampling menu
    rV = 0;
  }

  return rV;
}

void sharpening_menu(int pathMod, int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: displays the menu for sharpending module               *
\***************************************************************/
{
  int i, rV, selection_found = 0, option = 0, size_of_string = 255, sharpening_method = 0, usedEq = (userHom == -59 || useRegen) ? 1 : 0;
  char ch, *inputFile = NULL, *tempStr = NULL, sharpening_string[2][20] = {"Newton's method", "an endgame approach"};
  FILE *LIST = NULL;

  do
  { // find the sharpening option
    printf("\nSharpening options:\n");
    printf("  1. Sharpen all endpoints\n");
    printf("  2. Sharpen endpoints listed in a file\n");
    printf("  3. Manually input endpoints to sharpen\n");
    printf("  4. Recreate output (i.e., run the post-processor)\n");
    printf("  5. Change the number of sharpening digits (currently %d digits)\n", T->sharpenDigits);
    printf("  6. Change the sharpening method (currently %s)\n", sharpening_string[sharpening_method]);
    printf("  7. Quit\n");
    printf("Please enter your option: ");
    rV = scanf("%d", &option);

    if (rV < 0)
    { // at EOF - so we need to quit
      option = 7;
      selection_found = 1;
    }
    else // we are not at EOF
    {
      // flush the buffer
      do
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");
        selection_found = 0;
      }
      else if (option < 1 || option > 7)
      {
        printf("\n%d is not a valid option!\n", option);
        selection_found = 0;
      }
      else if (option == 5)
      { // change the number of sharpening digits
        do
        {
          printf("\nPlease input the number of sharpening digits (currently %d digits): ", T->sharpenDigits);
          rV = scanf("%d", &i);

          if (rV < 0)
          { // at EOF - so we need to quit
            option = 7;
            selection_found = 1;
            i = 1;
          }
          else // we are not at EOF
          { 
            // flush the buffer
            do
            {
              ch = getchar();
            } while (ch != EOF && ch != '\n');

            if (rV == 0)
            { // input was not read in correctly - so quit
              printf("\nThe input was not read in correctly so the number of sharpening digits will not be changed!\n");
              i = 1;
            }
            else // scanf was successful
            { // make sure the current precision can handle this many digits

              if (T->MPType == 0)
                rV = prec_to_digits(52); 
              else if (T->MPType == 1)
                rV = prec_to_digits(T->Precision) - 1;
              else
                rV = INT_MAX;

              if (i > rV)
              { // too many digits for this precision
                T->sharpenDigits = rV - 1;
                printf("\nThe current precision cannot handle %d digits reliably so the number of sharpening digits has been set to %d digits!\n", i, T->sharpenDigits);
              }
              else if (i > 0) // && i <= rV
              { // input is valid so store to sharpenDigits
                T->sharpenDigits = i;
                printf("\nThe number of sharpening digits has been changed to %d digits!\n", T->sharpenDigits);
              }
              else // i <= 0
              { // input not valid positive number
                printf("\nThe number of sharpening digits must be positive!\n");
              }
            }
            selection_found = 0;
          }
        } while (i <= 0);
      }
      else if (option == 6)
      { // change the sharpening method
        if (userHom == -59)
        { // not possible to change the method
          printf("\nSince a diagonal equation-by-equation method was used, only Newton's method is available.\n");
        }
        else if (useRegen)
        { // not possible to change the method
          printf("\nSince regeneration was used, only Newton's method is available.\n");
        }
        else if (midOutExists == 0)
        { // not possible to change the method
          printf("\nSince 'midpath_data' does not exist, only Newton's method is available.\n");
        }
        else
        { // decide on either Newton's method or endgame
          printf("\nSharpening method:\n");
          printf("  1. Newton's method\n");
          printf("  2. Endgame\n");
          printf("Please enter your option: ");
          rV = scanf("%d", &i);

          if (rV < 0)
          { // at EOF - so we need to quit
            option = 7;
            selection_found = 1;
            i = 1;
          }
          else // we are not at EOF
          {
            // flush the buffer
            do
            {
              ch = getchar();
            } while (ch != EOF && ch != '\n');

            if (rV == 0)
            { // input was not read in correctly - so quit
              printf("\nThe input was not read in correctly so the sharpening method will not be changed!\n");
              i = 1;
            }
            else // scanf was successful
            { 
              if (i < 1 || i > 2)
              {
                printf("\nSince the option (%d) was invalid, the sharpening method will not be changed!\n", i);
              }
              else
              {
                sharpening_method = i - 1;
              }
            }
            selection_found = 0;
          }
        } 
      }
      else // valid selection
        selection_found = 1;
    }
  } while (!selection_found);
  printf("\n");

  if (option == 1)
  { // sharpen all successful endpoints

    // test for sharpening method and number of digits
    if (sharpening_method && T->sharpenDigits > 300)
    { // print message
      printf("\nSince an endgame is being utilized, sharpening is limited to 300 digits.\n\n");
      T->sharpenDigits = 300;
    }

    sharpen_all_endpoints(pathMod, userHom, useRegen, endPoints, num_endPoints, T, OUT, rawOUT, sharpening_method, midOutExists, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }
  else if (option == 2)
  { // sharpen endpoints listed in a file
    tempStr = (char *)bmalloc(((int) log10(size_of_string) + 10) * sizeof(char));
    snprintf(tempStr, size_of_string + 10, "%%%ds", size_of_string);
    inputFile = (char *)bmalloc((size_of_string + 1) * sizeof(char));
    for (i = 0; i <= size_of_string; i++)
      inputFile[i] = '\0';
    do
    { // find the name of the file
      printf("Please input the name of the file that lists the\n endpoints to sharpen or type 'quit' or 'exit' (max of %d characters): ", size_of_string);
      rV = scanf(tempStr, inputFile);

      if (rV <= 0)
      { // at EOF so we need to quit
        option = 7; 
        selection_found = 1;
      }
      else
      {
        // flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (!strcmp(inputFile, "'quit'") || !strcmp(inputFile, "quit") || !strcmp(inputFile, "'exit'") || !strcmp(inputFile, "exit"))
        { // user wants to quit
          option = 7;
          selection_found = 1;
        }
        else
        { // check for existence of this file
          LIST = fopen(inputFile, "r");
          if (LIST == NULL)
          { // not exist
            printf("\nA file named \"%s\" does not exist!\n\n", inputFile);
            selection_found = 0;
          }
          else
          { // exists!
            fclose(LIST);
            selection_found = 1;
          }
        }
      }
    } while (!selection_found);
    printf("\n");

    // make sure the user still wants to continue
    if (option == 2)
    { // sharpen endpoints listed in a file

      // test for sharpening method and number of digits
      if (sharpening_method && T->sharpenDigits > 300)
      { // print message
        printf("\nSince an endgame is being utilized, sharpening is limited to 300 digits.\n\n");
        T->sharpenDigits = 300;
      }

      sharpen_endpoints_file(pathMod, userHom, useRegen, inputFile, endPoints, num_endPoints, T, OUT, rawOUT, sharpening_method, midOutExists, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }

    // free the memory
    free(tempStr);
    free(inputFile);
  }
  else if (option == 3)
  { // sharpen endpoints the user inputs manually

    // test for sharpening method and number of digits
    if (sharpening_method && T->sharpenDigits > 300)
    { // print message
      printf("\nSince an endgame is being utilized, sharpening is limited to 300 digits.\n\n");
      T->sharpenDigits = 300;
    }

    sharpen_endpoints_manually(userHom, useRegen, endPoints, num_endPoints, T, OUT, rawOUT, sharpening_method, midOutExists, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }
  else if (option == 4)
  { // run the post-processor
    FILE *mainOUT;

    // setup raw_data for this run
    create_raw_data_from_endPoints(endPoints, num_endPoints, T->numVars, rawOUT);
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
    printZeroDimRelevantData(ED_d, ED_mp, T->MPType, usedEq, rawOUT);
    
    // move midpath_data
    if (midOutExists)
    {
      fclose(midIN);
      rename("midpath_data_old", "midpath_data");
    }

    // do the rest of the output
    mainOUT = fopen("main_data", "w");
    if (T->MPType == 0 || T->MPType == 2)
      zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_d->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
    else
      zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_mp->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
    fclose(mainOUT);
  }

  // this is broken off of the if-else tree for a reason - since the user can change options inside the tree!
  if (option == 7)
  { // quit - so just setup raw_data
    create_raw_data_from_endPoints(endPoints, num_endPoints, T->numVars, rawOUT);
    fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
    printZeroDimRelevantData(ED_d, ED_mp, T->MPType, usedEq, rawOUT);

    // move midpath_data
    if (midOutExists)
    {
      fclose(midIN);
      rename("midpath_data_old", "midpath_data");
    }
  }

  remove("this_input");

  return;
}

void sharpen_all_endpoints(int pathMod, int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens every sucessful endpoint in endPoints         *
\***************************************************************/
{
  int i, retVal, failures, singular;
  int finiteSuccess = 0, infiniteSuccess = 0, nonfiniteSuccess = 0, usedEq = (userHom == -59 || useRegen) ? 1 : 0;
  FILE *mainOUT;

  // initialize counters
  failures = singular = 0;

  for (i = 0; i < num_endPoints; i++)
  { 
    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sharpening endpoint %d of %d\n", i, num_endPoints);

    // do the sharpening
    if (sharpening_option == 0)
      retVal = sharpen_endpoint(&endPoints[i], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    else // use endgame
      retVal = sharpen_endpoint_using_endgame(&endPoints[i], endPoints[i].path_num, T, OUT, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

    // check the outcome
    if (retVal == retVal_sharpening_singular_endpoint)
      singular++;
    else if (retVal == retVal_sharpening_failed)
      failures++;
  }

  // setup raw_data for this run
  create_raw_data_from_endPoints(endPoints, num_endPoints, T->numVars, rawOUT);
  fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
  printZeroDimRelevantData(ED_d, ED_mp, T->MPType, usedEq, rawOUT);

  // setup midpath_data
  if (midOutExists)
  {
    fclose(midIN);
    rename("midpath_data_old", "midpath_data");
  }

  // do the rest of the output
  mainOUT = fopen("main_data", "w");
  if (T->MPType == 0 || T->MPType == 2)
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_d->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  else
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_mp->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  fclose(mainOUT);

  // count the number of successes
  countFiniteSuccess(&finiteSuccess, &infiniteSuccess, &nonfiniteSuccess, num_endPoints, endPoints, -1, NULL);

  // print the summary
  printSharpeningFailureSummary(num_endPoints, failures, singular, finiteSuccess, infiniteSuccess, nonfiniteSuccess);

  return;
}

void countFiniteSuccess(int *finiteSuccess, int *infiniteSuccess, int *nonfiniteSuccess, int numStartPts, post_process_t *endPoints, int num_path_nums, int *path_nums)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: counts the number of success that are finite, infinite *
*   and unable to determine finite from infinite                *
* num_path_nums: -1 - all endpoints, else number of path_nums   *
\***************************************************************/
{
  int i, j, count;

  // initialize all to 0
  *finiteSuccess = *infiniteSuccess = *nonfiniteSuccess = 0;

  if (num_path_nums == -1)
  { // look at all of the endpoints
    for (i = 0; i < numStartPts; i++)
      if (endPoints[i].success == 1)
      {

        if (endPoints[i].isFinite == 1)
          *finiteSuccess = *finiteSuccess + 1;
        else if (endPoints[i].isFinite == 0)
          *infiniteSuccess = *infiniteSuccess + 1;
        else // do not know what finite is - user homotopy or already homogenized
          *nonfiniteSuccess = *nonfiniteSuccess + 1;
      }
  }
  else
  { // look through path_nums
    for (i = 0; i < num_path_nums; i++)
    { // check to see if path_num[i] is one of the endPoints
      count = -1;
      for (j = 0; j < numStartPts; j++)
        if (endPoints[j].path_num == path_nums[i])
        { // found the array number
          count = j;
          j = numStartPts;
        }

      if (count >= 0) 
      { // path_nums[i] exists
        if (endPoints[count].success == 1)
        {
          if (endPoints[count].isFinite == 1)
            *finiteSuccess = *finiteSuccess + 1;
          else if (endPoints[count].isFinite == 0)
            *infiniteSuccess = *infiniteSuccess + 1;
          else // do not know what finite is - user homotopy or already homogenized
            *nonfiniteSuccess = *nonfiniteSuccess + 1;
        }
      }
    }
  }

  return;
}

void sharpen_endpoints_file(int pathMod, int userHom, int useRegen, char *fileName, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens the endpoints listed in fileName              *
\***************************************************************/
{
  int i, j, count = 0, size = 1, cont = 1, retVal, attempted, failures, singular, usedEq = (userHom == -59 || useRegen) ? 1 : 0;
  int finiteSuccess = 0, infiniteSuccess = 0, nonfiniteSuccess = 0;
  int *path_nums = (int *)bmalloc(size * sizeof(int));
  FILE *mainOUT, *IN = fopen(fileName, "r");

  // initialize counters
  count = attempted = failures = singular = 0;

  // loop to read the selected endpoints in from IN
  cont = 1;
  while (cont)
  { // make sure we have room to store the next one
    if (size <= count)
    {
      path_nums = (int *)brealloc(path_nums, 2 * size * sizeof(int));
      size *= 2;
    }

    // scan in an integer
    i = -1; // used to see the end of the file
    fscanf(IN, "%d", &i);
    if (i < 0)
    { // invalid integer scanned in so we get out of the loop
      cont = 0;
    }
    else
    { // non-negative integer scanned in so we add it to the list
      path_nums[count] = i;
      count++;
      cont = 1;
    }
  }
  // store the number of endpoints that was read in & close IN
  size = count;
  fclose(IN);

  // loop through the endpoints and sharpen the ones that we can 
  for (i = 0; i < size; i++)
  {
    // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Sharpening endpoint %d of %d\n", i, size);

    // check to see if path_num[i] is one of the endPoints
    count = -1;
    for (j = 0; j < num_endPoints; j++)
      if (endPoints[j].path_num == path_nums[i])
      { // found the array number
        count = j;
        j = num_endPoints;
      }

    // check to see we can be sharpen the path
    if (count >= 0)
    { // check to see if we can sharpen it
      attempted++;

      // do the sharpening
      if (sharpening_option == 0)
        retVal = sharpen_endpoint(&endPoints[count], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      else // use endgame
        retVal = sharpen_endpoint_using_endgame(&endPoints[count], endPoints[count].path_num, T, OUT, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

      // check the outcome
      if (retVal == retVal_sharpening_singular_endpoint)
        singular++;
      else if (retVal == retVal_sharpening_failed)
        failures++;
    }
    else
    { // path number does not exist
      printf("NOTE: The path number %d does not exist.\n", path_nums[i]);
    }
  }

  // setup raw_data for this run
  create_raw_data_from_endPoints(endPoints, num_endPoints, T->numVars, rawOUT);
  fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
  printZeroDimRelevantData(ED_d, ED_mp, T->MPType, usedEq, rawOUT);

  // setup midpath_data
  if (midOutExists)
  {
    fclose(midIN);
    rename("midpath_data_old", "midpath_data");
  }

  // do the rest of the output
  mainOUT = fopen("main_data", "w");
  if (T->MPType == 0 || T->MPType == 2)
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_d->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  else
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_mp->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  fclose(mainOUT);


  // count the number of successes
  countFiniteSuccess(&finiteSuccess, &infiniteSuccess, &nonfiniteSuccess, num_endPoints, endPoints, size, path_nums);

  // print the summary
  printSharpeningFailureSummary(attempted, failures, singular, finiteSuccess, infiniteSuccess, nonfiniteSuccess);

  // free the memory
  free(path_nums);

  return;
}

void sharpen_endpoints_manually(int userHom, int useRegen, post_process_t *endPoints, int num_endPoints, tracker_config_t *T, FILE *OUT, FILE *rawOUT, int sharpening_option, int midOutExists, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens the endpoints listed in fileName              *
\***************************************************************/
{
  int i, retVal, path_num, array_num, min_path_num, max_path_num, usedEq = (userHom == -59 || useRegen) ? 1 : 0;
  char ch;
  FILE *mainOUT;

  // find the statistics on the endPoints
  min_path_num = INT_MAX;
  max_path_num = 0;
  for (i = 0; i < num_endPoints; i++)
  {
    if (min_path_num > endPoints[i].path_num)
      min_path_num = endPoints[i].path_num;
    if (max_path_num < endPoints[i].path_num)
      max_path_num = endPoints[i].path_num;
  }

  // initialize path_num to something other than -9
  path_num = -2;
  do 
  { // read in an endpoint
    if (path_num != -9)
      printf("The path numbers range between %d and %d.\n", min_path_num, max_path_num);
    printf("Please enter a path number to sharpen (-1 to quit and -9 to display the list of path numbers): ");
    retVal = scanf("%d", &path_num);

    if (retVal < 0)
    { // at EOF - so we need to quit
      path_num = -1;
    }
    else
    {
      // flush the buffer
      do
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (retVal == 0)
      { // invalid entry
        printf("\nThe input was not read in correctly!\n");
        path_num = -2;
      }
      else if (path_num == -9)
      { // display all path numbers
        printf("\nThe list of path numbers:\n");
        for (i = 0; i < num_endPoints; i++)
          printf("%d\n", endPoints[i].path_num);
        printf("\n");
      }
      else if (path_num != -1)
      { // check to see if the path number is valid
        if (path_num < min_path_num || path_num > max_path_num)
        { // entry not in the valid range
          printf("\n%d is not a valid path number!\n", path_num);
        }
        else
        { // do a loop to find the existence of the path number
          array_num = -1; // initialize to not found
          for (i = 0; i < num_endPoints; i++)
            if (endPoints[i].path_num == path_num)
            { // found it!
              array_num = i;
              i = num_endPoints;
            }

          if (array_num == -1)
          { // path_num is not valid
            printf("\n%d is not a valid path number!\n", path_num);
          }
          else
          { // we have a valid endpoint so we can sharpen it
            printf("\nSharpening path number %d\n", path_num);

            // do the sharpening
            if (sharpening_option == 0)
              retVal = sharpen_endpoint(&endPoints[array_num], T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
            else // use endgame
              retVal = sharpen_endpoint_using_endgame(&endPoints[array_num], endPoints[array_num].path_num, T, OUT, midIN, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

            // check the outcome
            if (retVal == 0)
              printf("Sharpening was successful.\n");
            else if (retVal == retVal_sharpening_singular_endpoint)
              printf("The endpoint failed the singularity test.\n");
            else if (retVal == retVal_sharpening_failed)
              printf("Sharpening failed to converge - try adjusting the number of sharpening digits\n or consider using adaptive precision.\n");
          }
        }
      }
    }
    if (path_num != -1 && path_num != -9)
      printf("\n");
  } while (path_num != -1);
  printf("\n");  

  // setup raw_data for this run
  create_raw_data_from_endPoints(endPoints, num_endPoints, T->numVars, rawOUT);
  fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
  printZeroDimRelevantData(ED_d, ED_mp, T->MPType, usedEq, rawOUT);

  // setup midpath_data
  if (midOutExists)
  {
    fclose(midIN);
    rename("midpath_data_old", "midpath_data");
  }

  // do the rest of the output
  mainOUT = fopen("main_data", "w");
  if (T->MPType == 0 || T->MPType == 2)
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_d->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  else
    zeroDimPostProcess(mainOUT, endPoints, num_endPoints, T->numVars, T->final_tol_times_mult, T, &ED_mp->preProcData, 0, 0, "this_input", useRegen == 1 && userHom == 0, userHom == -59);
  fclose(mainOUT);

  return;
}

int sharpen_endpoint_using_endgame(post_process_t *endPoint, int pathNum, tracker_config_t *T, FILE *OUT, FILE *midIN, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: retVal of sharpening module                    *
* NOTES: sharpens the endPoint                                  *
\***************************************************************/
{
  int i, numVars = T->numVars, retVal = 0, foundPathNum = 0, inputPathNum = -1;

  // locate the midpoint 
  rewind(midIN);
  while (foundPathNum == 0)
  { // read in path number
    retVal = fscanf(midIN, "%d\n", &inputPathNum);
    if (retVal <= 0)
    { // error reading in number
      printf("\nSince a point on the path was not found, Newton's method is being used!\n\n");
      foundPathNum = 1;
      retVal = sharpen_endpoint(endPoint, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }
    else if (inputPathNum == pathNum)
    { // found the proper path number - do the actual endgame sharpening
      foundPathNum = 1;

      // setup the final tolerance based on sharpen digits
      T->final_tolerance = pow(10, -T->sharpenDigits);

      // setup the other items
      T->endgameSwitch = 1;
      T->currentNewtonTol = T->endgameNewtonTol;
      T->currentStepSize = MIN(T->maxStepSize, T->endgameBoundary);
      T->minStepSize = T->minStepSizeDuringEndGame;
      T->error_at_latest_sample_point = 0;

      // read in the endpoint and run the endgame
      if (T->MPType == 0)
      { // use double
        point_data_d Final;
        point_d startPoint, lastApprox;

        init_point_data_d(&Final, numVars);
        init_point_d(startPoint, numVars);
        init_point_d(lastApprox, numVars);

        // setup startPoint
        startPoint->size = numVars;
        for (i = 0; i < numVars; i++)
          fscanf(midIN, "%lf%lf", &startPoint->coord[i].r, &startPoint->coord[i].i);

        // run the endgame
        if (T->endgameNumber == 1)
        { // run PS EG
          PSEG_samples_struct_d PSEG_samples;
          init_PSEG_samples_struct_d(&PSEG_samples, T->num_PSEG_sample_points);
          PSEG_samples.num_samples = 1;
          // setup start point & time
          point_cp_d(PSEG_samples.samples[0].point, startPoint);
          set_double_d(PSEG_samples.samples[0].time, T->endgameBoundary, 0);
          PSEG_samples.samples[0].cycle_num = 0;

          // run the endgame
          retVal = PSEG_struct_d(&Final, lastApprox, &PSEG_samples, T, OUT, ED_d, eval_func_d, zero_dim_dehom);

          // clear PSEG_samples
          clear_PSEG_samples_struct_d(&PSEG_samples);
        }
        else
        { // run Cauchy EG
          int cycle = 0, samples_per_loop = 0;
          comp_d finalTime;
          point_data_d *endSamples = NULL;

          // setup start point & time
          point_cp_d(Final.point, startPoint);
          set_double_d(Final.time, T->endgameBoundary, 0);
          Final.cycle_num = 0;

          // setup ending time
          set_double_d(finalTime, T->targetT, 0);

          // rund the endgame
          retVal = CauchyEG_main_d2(&Final, lastApprox, &endSamples, &cycle, &samples_per_loop, finalTime, &Final, T, OUT, ED_d, eval_func_d, zero_dim_dehom);

          // clear endSamples since is it not being used
          if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
          { // clear
            for (i = samples_per_loop * cycle - 1; i >= 0; i--)
              clear_point_data_d(&endSamples[i]);
            free(endSamples);
          }
        }

        // check the outcome
        if (retVal == 0)
        { // success
          endPoint->success = 1;
        }
        else
        { // failure
          retVal = endPoint->success = retVal_sharpening_failed;
        }

        if (endPoint->success == 1)
        { // save data
          if (endPoint->sol_prec > 52)
          { // clear _mp and save in _d
            endPoint->sol_d = (comp_d *)bmalloc(Final.point->size * sizeof(comp_d));
            for (i = 0; i < endPoint->size_sol; i++)
            {
              clear_mp(endPoint->sol_mp[i]);
            }
            mpf_clear(endPoint->function_resid_mp);
            mpf_clear(endPoint->newton_resid_mp);
            free(endPoint->sol_mp);
            endPoint->sol_mp = NULL;
          }
          else if (endPoint->size_sol < Final.point->size)
          { // this should never happen!!!
            endPoint->sol_d = (comp_d *)brealloc(endPoint->sol_d, Final.point->size * sizeof(comp_d));
          }

          endPoint->sol_prec = 52;
          endPoint->size_sol = Final.point->size;
          for (i = 0; i < endPoint->size_sol; i++)
          {
            set_d(endPoint->sol_d[i], &Final.point->coord[i]);
          }
          endPoint->newton_resid_d = T->latest_newton_residual_d;
          endPoint->final_t = T->t_val_at_latest_sample_point; 
          endPoint->accuracy_estimate = T->error_at_latest_sample_point;
          endPoint->cycle_num = Final.cycle_num;
          endPoint->first_increase = 0;
          // find function residual and condition number
          findFunctionResidual_conditionNumber_d(&endPoint->function_resid_d, &endPoint->cond_est, &Final, ED_d, eval_func_d);
        }

        // clear memory
        clear_point_data_d(&Final);
        clear_point_d(startPoint);
        clear_point_d(lastApprox);
      }
      else if (T->MPType == 1)
      { // use multi precision
        point_data_mp Final;
        point_mp startPoint, lastApprox;

        init_point_data_mp(&Final, numVars);
        init_point_mp(startPoint, numVars);
        init_point_mp(lastApprox, numVars);

        // setup startPoint
        startPoint->size = numVars;
        for (i = 0; i < numVars; i++)
        {
          mpf_inp_str(startPoint->coord[i].r, midIN, 10);
          mpf_inp_str(startPoint->coord[i].i, midIN, 10);
        }

        // run the endgame
        if (T->endgameNumber == 1)
        { // run PS EG
          PSEG_samples_struct_mp PSEG_samples;
          init_PSEG_samples_struct_mp(&PSEG_samples, T->num_PSEG_sample_points, T->Precision);
          PSEG_samples.num_samples = 1;
          // setup start point & time
          point_cp_mp(PSEG_samples.samples[0].point, startPoint);
          set_double_mp(PSEG_samples.samples[0].time, T->endgameBoundary, 0);
          PSEG_samples.samples[0].cycle_num = 0;

          // run the endgame
          retVal = PSEG_struct_mp(&Final, lastApprox, &PSEG_samples, T, OUT, ED_mp, eval_func_mp, zero_dim_dehom);

          // clear PSEG_samples
          clear_PSEG_samples_struct_mp(&PSEG_samples);
        }
        else
        { // run Cauchy EG
          int cycle = 0, samples_per_loop = 0;
          comp_mp finalTime;
          point_data_mp *endSamples = NULL;

          init_mp(finalTime);

          // setup start point & time
          point_cp_mp(Final.point, startPoint);
          set_double_mp(Final.time, T->endgameBoundary, 0);
          Final.cycle_num = 0;

          // setup ending time
          set_double_mp(finalTime, T->targetT, 0);

          // rund the endgame
          retVal = CauchyEG_main_mp2(&Final, lastApprox, &endSamples, &cycle, &samples_per_loop, finalTime, &Final, T, OUT, ED_mp, eval_func_mp, zero_dim_dehom);

          // clear endSamples since is it not being used
          if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
          { // clear
            for (i = samples_per_loop * cycle - 1; i >= 0; i--)
              clear_point_data_mp(&endSamples[i]);
            free(endSamples);
          }

          clear_mp(finalTime);
        }

        // check the outcome
        if (retVal == 0)
        { // success
          endPoint->success = 1;
        }
        else
        { // failure
          retVal = endPoint->success = retVal_sharpening_failed;
        }

        if (endPoint->success == 1)
        { // save data
          if (endPoint->sol_prec < 64)
          { // clear _d and save in _mp
            free(endPoint->sol_d);
            endPoint->sol_d = NULL;

            endPoint->sol_mp = (comp_mp *)bmalloc(Final.point->size * sizeof(comp_mp));
            for (i = 0; i < Final.point->size; i++)
            {
              init_mp(endPoint->sol_mp[i]);
            }
            mpf_init(endPoint->function_resid_mp);
            mpf_init(endPoint->newton_resid_mp);
          }
          else if (endPoint->size_sol < Final.point->size)
          { // this should never happen!!!
            for (i = 0; i < endPoint->size_sol; i++)
            {
              clear_mp(endPoint->sol_mp[i]);
            }
            endPoint->sol_mp = (comp_mp *)brealloc(endPoint->sol_mp, Final.point->size * sizeof(comp_mp));
            for (i = 0; i < Final.point->size; i++)
            {
              init_mp(endPoint->sol_mp[i]);
            }
          }
          endPoint->sol_prec = T->Precision;
          endPoint->size_sol = Final.point->size;
          for (i = 0; i < endPoint->size_sol; i++)
          {
            set_mp(endPoint->sol_mp[i], &Final.point->coord[i]);
          }
          mpf_set(endPoint->newton_resid_mp, T->latest_newton_residual_mp);
          endPoint->final_t = T->t_val_at_latest_sample_point;
          endPoint->accuracy_estimate = T->error_at_latest_sample_point;
          endPoint->cycle_num = Final.cycle_num;
          endPoint->first_increase = 0;
          // find function residual and condition number
          findFunctionResidual_conditionNumber_mp(endPoint->function_resid_mp, &endPoint->cond_est, &Final, ED_mp, eval_func_mp);
        }

        // clear memory
        clear_point_data_mp(&Final);
        clear_point_mp(startPoint);
        clear_point_mp(lastApprox);
      }
      else
      { // use adaptive precision
        int final_prec, lastApprox_prec;

        point_data_d Final_d;
        point_data_mp Final_mp;
        point_d startPoint_d, lastApprox_d;
        point_mp startPoint_mp, lastApprox_mp;
       
        init_point_data_d(&Final_d, numVars);
        init_point_d(startPoint_d, numVars);
        init_point_d(lastApprox_d, numVars);
        init_point_data_mp2(&Final_mp, numVars, MAX(endPoint->sol_prec, T->Precision));
        init_point_mp2(startPoint_mp, numVars, MAX(endPoint->sol_prec, T->Precision));
        init_point_mp2(lastApprox_mp, numVars, MAX(endPoint->sol_prec, T->Precision));

        // setup startPoint -- use same precision as solution
        if (endPoint->sol_prec < 64)
        { // setup using _d
          startPoint_d->size = numVars;
          for (i = 0; i < numVars; i++)
            fscanf(midIN, "%lf%lf", &startPoint_d->coord[i].r, &startPoint_d->coord[i].i);
        }
        else
        { // setup using _mp
          startPoint_mp->size = numVars;
          for (i = 0; i < numVars; i++)
          {
            mpf_inp_str(startPoint_mp->coord[i].r, midIN, 10);
            mpf_inp_str(startPoint_mp->coord[i].i, midIN, 10);
          }
        }

        // run the endgame
        if (T->endgameNumber == 1)
        { // run PS EG
          PSEG_samples_struct_amp PSEG_samples;
          init_PSEG_samples_struct_amp(&PSEG_samples, T->num_PSEG_sample_points, endPoint->sol_prec);
          PSEG_samples.num_samples = 1;

          PSEG_samples.max_prec[0] = endPoint->sol_prec;
          if (endPoint->sol_prec < 64)
          { // setup start point & time
            point_cp_d(PSEG_samples.samples_d[0].point, startPoint_d);
            set_double_d(PSEG_samples.samples_d[0].time, T->endgameBoundary, 0);
            PSEG_samples.samples_d[0].cycle_num = 0;
          }
          else
          { // setup start point & time
            point_cp_mp(PSEG_samples.samples_mp[0].point, startPoint_mp);
            set_double_mp(PSEG_samples.samples_mp[0].time, T->endgameBoundary, 0);
            PSEG_samples.samples_mp[0].cycle_num = 0;
          }

          // run the endgame
          retVal = PSEG_struct_amp(&final_prec, &endPoint->first_increase, &Final_d, &Final_mp, lastApprox_d, lastApprox_mp, &lastApprox_prec, &PSEG_samples, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, zero_dim_dehom);

          // clear PSEG_samples
          clear_PSEG_samples_struct_amp(&PSEG_samples);
        }
        else
        { // run Cauchy EG
          int cycle = 0, samples_per_loop = 0;
          comp_d finalTime_d;
          comp_mp finalTime_mp;
          point_data_d *endSamples_d = NULL;
          point_data_mp *endSamples_mp = NULL;

          init_mp2(finalTime_mp, MAX(endPoint->sol_prec, T->Precision));

          // setup start point & time
          if (endPoint->sol_prec < 64)
          { // setup start point & time
            point_cp_d(Final_d.point, startPoint_d);
            set_double_d(Final_d.time, T->endgameBoundary, 0);
            Final_d.cycle_num = 0;
          }
          else
          { // setup start point & time
            point_cp_mp(Final_mp.point, startPoint_mp);
            set_double_mp(Final_mp.time, T->endgameBoundary, 0);
            Final_mp.cycle_num = 0;
          }

          // setup ending time
          set_double_d(finalTime_d, T->targetT, 0);
          set_double_mp(finalTime_mp, T->targetT, 0);

          // rund the endgame
          retVal = CauchyEG_main_amp2(&final_prec, &Final_d, &Final_mp, lastApprox_d, lastApprox_mp, &lastApprox_prec, &endSamples_d, &endSamples_mp, &cycle, &samples_per_loop, &endPoint->first_increase, finalTime_d, finalTime_mp, &Final_d, &Final_mp, endPoint->sol_prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, zero_dim_dehom);

          // clear endSamples since is it not being used
          if (retVal == 0 || retVal == retVal_EG_failed_to_converge)
          { // clear
            if (final_prec < 64) 
            {
              for (i = samples_per_loop * cycle - 1; i >= 0; i--)
                clear_point_data_d(&endSamples_d[i]);
              free(endSamples_d);
            }
            else
            {
              for (i = samples_per_loop * cycle - 1; i >= 0; i--)
                clear_point_data_mp(&endSamples_mp[i]);
              free(endSamples_mp);
            }
          }
        }

        // check the outcome
        if (retVal == 0)
        { // success
          endPoint->success = 1;
        }
        else
        { // failure
          retVal = endPoint->success = retVal_sharpening_failed;
        }

        if (endPoint->success == 1)
        { // save data
          if (final_prec < 64)
          { // save to _d
            if (endPoint->sol_prec > 52)
            { // clear _mp and save in _d
              endPoint->sol_d = (comp_d *)bmalloc(Final_d.point->size * sizeof(comp_d));
              for (i = 0; i < endPoint->size_sol; i++)
              {
                clear_mp(endPoint->sol_mp[i]);
              }
              mpf_clear(endPoint->function_resid_mp);
              mpf_clear(endPoint->newton_resid_mp);
              free(endPoint->sol_mp);
              endPoint->sol_mp = NULL;
            }
            else if (endPoint->size_sol < Final_d.point->size)
            { // this should never happen!!!
              endPoint->sol_d = (comp_d *)brealloc(endPoint->sol_d, Final_d.point->size * sizeof(comp_d));
            }

            endPoint->sol_prec = 52;
            endPoint->size_sol = Final_d.point->size;
            for (i = 0; i < endPoint->size_sol; i++)
            {
              set_d(endPoint->sol_d[i], &Final_d.point->coord[i]);
            }
            endPoint->newton_resid_d = T->latest_newton_residual_d;
            endPoint->final_t = T->t_val_at_latest_sample_point;
            endPoint->accuracy_estimate = T->error_at_latest_sample_point;
            endPoint->cycle_num = Final_d.cycle_num;
            // find function residual and condition number
            findFunctionResidual_conditionNumber_d(&endPoint->function_resid_d, &endPoint->cond_est, &Final_d, ED_d, eval_func_d);
          }
          else
          { // save to _mp
            if (endPoint->sol_prec < 64)
            { // clear _d and save in _mp
              free(endPoint->sol_d);
              endPoint->sol_d = NULL;

              endPoint->sol_mp = (comp_mp *)bmalloc(Final_mp.point->size * sizeof(comp_mp));
              for (i = 0; i < Final_mp.point->size; i++)
              {
                init_mp2(endPoint->sol_mp[i], final_prec);
              }
              mpf_init2(endPoint->function_resid_mp, final_prec);
              mpf_init2(endPoint->newton_resid_mp, final_prec);
            }
            else if (endPoint->size_sol < Final_mp.point->size)
            { // this should never happen!!!
              for (i = 0; i < endPoint->size_sol; i++)
              {
                clear_mp(endPoint->sol_mp[i]);
              }
              endPoint->sol_mp = (comp_mp *)brealloc(endPoint->sol_mp, Final_mp.point->size * sizeof(comp_mp));
              for (i = 0; i < Final_mp.point->size; i++)
              {
                init_mp2(endPoint->sol_mp[i], final_prec);
              }
            }
            endPoint->sol_prec = final_prec; 
            endPoint->size_sol = Final_mp.point->size;
            for (i = 0; i < endPoint->size_sol; i++)
            {
              set_mp(endPoint->sol_mp[i], &Final_mp.point->coord[i]);
            }
            mpf_set(endPoint->newton_resid_mp, T->latest_newton_residual_mp);
            endPoint->final_t = T->t_val_at_latest_sample_point;
            endPoint->accuracy_estimate = T->error_at_latest_sample_point;
            endPoint->cycle_num = Final_mp.cycle_num;
            endPoint->first_increase = 0;
            // find function residual and condition number
            findFunctionResidual_conditionNumber_mp(endPoint->function_resid_mp, &endPoint->cond_est, &Final_mp, ED_mp, eval_func_mp);  
          }
        }

        // clear memory
        clear_point_data_d(&Final_d);
        clear_point_d(startPoint_d);
        clear_point_d(lastApprox_d);
        clear_point_data_mp(&Final_mp);
        clear_point_mp(startPoint_mp);
        clear_point_mp(lastApprox_mp);
      }
    }
    else
    { // move past this point
      for (i = 0; i < numVars; i++)
        scanRestOfLine(midIN);
    }
  }
  
  return retVal; 
}

int sharpen_endpoint(post_process_t *endPoint, tracker_config_t *T, FILE *OUT, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: retVal of sharpening module                    *
* NOTES: sharpens the endPoint                                  *
\***************************************************************/
{
  int j, prec_out, retVal, prec_in = MAX(endPoint->sol_prec, 64);
  double current_tol_d = 0;
  mpf_t  current_tol_mp;
  point_data_d  in_d,  out_d;
  point_data_mp in_mp, out_mp;

  init_point_data_d(&in_d, 0); init_point_data_d(&out_d, 0);
  init_point_data_mp2(&in_mp, 0, prec_in); init_point_data_mp2(&out_mp, 0, prec_in);
  mpf_init2(current_tol_mp, prec_in);

  if (endPoint->sol_prec == 52)
  { // copy to in_d
    prec_in = 52;
    increase_size_point_d(in_d.point, endPoint->size_sol);
    in_d.point->size = endPoint->size_sol;
    for (j = 0; j < endPoint->size_sol; j++)
    {
      set_d(&in_d.point->coord[j], endPoint->sol_d[j]);
    }
    set_double_d(in_d.time, T->targetT, 0);
    // estimate the current accuracy
    current_tol_d = endPoint->newton_resid_d;
    if (current_tol_d <= 0)
      current_tol_d = 1e-14; // based on double precision
  }
  else
  { // copy to in_mp
    prec_in = endPoint->sol_prec;
    increase_size_point_mp(in_mp.point, endPoint->size_sol);
    in_mp.point->size = endPoint->size_sol;
    for (j = 0; j < endPoint->size_sol; j++)
    {
      set_mp(&in_mp.point->coord[j], endPoint->sol_mp[j]);
    }
    set_double_mp(in_mp.time, T->targetT, 0);

    // estimate the current accuracy
    mpf_set(current_tol_mp, endPoint->newton_resid_mp);
    if (mpf_cmp_ui(current_tol_mp, 0) <= 0)
    { // set based on precision
      retVal = (int) floor(prec_in * log10(2) - 2.5);
      int size = 1 + snprintf(NULL, 0, "1e%d", -retVal);
      char *str = (char *)bmalloc(size * sizeof(char));
      sprintf(str, "1e%d", -retVal);
      mpf_set_str(current_tol_mp, str, 10);

      free(str);
    }
  }

  // do the sharpening
  retVal = sharpen_zero_main(current_tol_d, current_tol_mp, prec_in, &out_d, &out_mp, &prec_out, &in_d, &in_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // check the outcome
  if (retVal == 0)
  { // success
    endPoint->success = 1;
  }
  else
  { // failure
    endPoint->success = retVal;
  }

  // check to see if sharpening did anything
  if (retVal != retVal_sharpening_singular_endpoint)
  { // sharpening was attempted so we need to copy the results back to endPoint
    if (prec_out == 52)
    { // copy to sol_d
      if (endPoint->sol_prec > 52)
      { // setup sol_d
        endPoint->sol_d = (comp_d *)bmalloc(out_d.point->size * sizeof(comp_d));
        // clear sol_mp & other MP
        for (j = 0; j < endPoint->size_sol; j++)
        {
          clear_mp(endPoint->sol_mp[j]);
        }
        mpf_clear(endPoint->function_resid_mp);
        mpf_clear(endPoint->newton_resid_mp);
        free(endPoint->sol_mp);
        endPoint->sol_mp = NULL;
      }
      else if (endPoint->size_sol < out_d.point->size)
      { // this should never happen!!!
        endPoint->sol_d = (comp_d *)brealloc(endPoint->sol_d, out_d.point->size * sizeof(comp_d));
      }
      // copy to sol_d
      for (j = 0; j < out_d.point->size; j++)
      {
        set_d(endPoint->sol_d[j], &out_d.point->coord[j]);
      }
      endPoint->size_sol = out_d.point->size;
      endPoint->sol_prec = prec_out;
      endPoint->newton_resid_d = T->latest_newton_residual_d;
      // find function residual & condition estimate
      findFunctionResidual_conditionNumber_d(&endPoint->function_resid_d, &endPoint->cond_est, &out_d, ED_d, eval_func_d);
    }
    else
    { // copy to sol_mp
      if (endPoint->sol_prec == 52)
      { // clear sol_d
        free(endPoint->sol_d);
        endPoint->sol_d = NULL;
        // setup sol_mp & other MP
        endPoint->sol_mp = (comp_mp *)bmalloc(out_mp.point->size * sizeof(comp_mp));
        for (j = 0; j < out_mp.point->size; j++)
        {
          init_mp2(endPoint->sol_mp[j], prec_out);
        }
        mpf_init2(endPoint->function_resid_mp, prec_out);
        mpf_init2(endPoint->newton_resid_mp, prec_out);
      }
      else if (endPoint->size_sol < out_mp.point->size)
      { // this should never happen!!!
        for (j = 0; j < endPoint->size_sol; j++)
        {
          clear_mp(endPoint->sol_mp[j]);
        }
        endPoint->sol_mp = (comp_mp *)brealloc(endPoint->sol_mp, out_mp.point->size * sizeof(comp_mp));
        for (j = 0; j < out_mp.point->size; j++)
        {
          init_mp2(endPoint->sol_mp[j], prec_out);
        }
      }
      else if (endPoint->sol_prec != prec_out)
      { // set to the corerct precision
        for (j = 0; j < endPoint->size_sol; j++)
        {
          change_prec_mp(endPoint->sol_mp[j], prec_out);
        }
        mpf_set_prec(endPoint->function_resid_mp, prec_out);
        mpf_set_prec(endPoint->newton_resid_mp, prec_out);
      }
      // copy to sol_mp
      for (j = 0; j < out_mp.point->size; j++)
      {
        set_mp(endPoint->sol_mp[j], &out_mp.point->coord[j]);
      }
      endPoint->size_sol = out_mp.point->size;
      endPoint->sol_prec = prec_out;
      mpf_set(endPoint->newton_resid_mp, T->latest_newton_residual_mp);
      // find function residual & condition estimate
      findFunctionResidual_conditionNumber_mp(endPoint->function_resid_mp, &endPoint->cond_est, &out_mp, ED_mp, eval_func_mp);
    }
  }

  // clear
  clear_point_data_d(&in_d); clear_point_data_d(&out_d);
  clear_point_data_mp(&in_mp); clear_point_data_mp(&out_mp);
  mpf_clear(current_tol_mp);

  return retVal;
}

void create_raw_data_from_endPoints(post_process_t *endPoints, int num_endPoints, int num_variables, FILE *rawOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the raw_data file to rawOUT                    *
\***************************************************************/
{
  int i, j;

  // top of rawOUT - number of variables and that we are doing zero dimensional
  fprintf(rawOUT, "%d\n%d\n", num_variables, 0);

  for (i = 0; i < num_endPoints; i++)
  { // print each endpoint to rawOUT

    // print path number & precision
    fprintf(rawOUT, "%d\n%d\n", endPoints[i].path_num, endPoints[i].sol_prec);

    if (endPoints[i].sol_prec == 52)
    { // print using double precision
      for (j = 0; j < num_variables; j++)
        fprintf(rawOUT, "%.15e %.15e\n", endPoints[i].sol_d[j]->r, endPoints[i].sol_d[j]->i);

      fprintf(rawOUT, "%.15e\n%.15e\n%.15e\n", endPoints[i].function_resid_d, endPoints[i].cond_est, endPoints[i].newton_resid_d);
    }
    else
    { // print using multiprecision
      long e1, e2;
      char *ch1, *ch2;

      for (j = 0; j < num_variables; j++)
      {
        if (mpfr_number_p(endPoints[i].sol_mp[j]->r) && mpfr_number_p(endPoints[i].sol_mp[j]->i))
        {
          ch1 = mpf_get_str(NULL, &e1, 10, 0, endPoints[i].sol_mp[j]->r);
          ch2 = mpf_get_str(NULL, &e2, 10, 0, endPoints[i].sol_mp[j]->i);

          if (mpf_sgn(endPoints[i].sol_mp[j]->r) >= 0)
            if (mpf_sgn(endPoints[i].sol_mp[j]->i) >= 0)
              fprintf(rawOUT, "0.%se%ld 0.%se%ld\n", ch1, e1, ch2, e2);
            else
              fprintf(rawOUT, "0.%se%ld -0.%se%ld\n", ch1, e1, &ch2[1], e2);
          else
            if (mpf_sgn(endPoints[i].sol_mp[j]->i) >= 0)
              fprintf(rawOUT, "-0.%se%ld 0.%se%ld\n", &ch1[1], e1, ch2, e2);
            else
              fprintf(rawOUT, "-0.%se%ld -0.%se%ld\n", &ch1[1], e1, &ch2[1], e2);

          mpfr_free_str(ch1);
          mpfr_free_str(ch2);
        }
        else
          fprintf(rawOUT, "NaN NaN\n");
      }

      mpf_out_str(rawOUT, 10, 15, endPoints[i].function_resid_mp);
      fprintf(rawOUT, "\n%.15e\n", endPoints[i].cond_est);
      mpf_out_str(rawOUT, 10, 15, endPoints[i].newton_resid_mp);
      fprintf(rawOUT, "\n");
    }
    // print the rest of the data
    fprintf(rawOUT, "%.15e\n%.15e\n%.15e\n%d\n%d\n", endPoints[i].final_t, endPoints[i].accuracy_estimate, endPoints[i].first_increase, endPoints[i].cycle_num, endPoints[i].success);
  }

  return;
}

void sharpen_endpoint_endgame(endgame_data_t *endPoint, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: retVal inside of endPoint is updated           *
* NOTES: sharpens endpoint and updates the endgame_data_t struct*
\***************************************************************/
{
  int prec_before = endPoint->prec;

  // try to sharpen the endpoint
  endPoint->retVal = sharpen_zero_main(T->final_tol_times_mult, NULL, 52, &endPoint->PD_d, &endPoint->PD_mp, &endPoint->prec, &endPoint->PD_d, &endPoint->PD_mp, endPoint->prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // update the values inside of endPoint
  if (endPoint->prec == 52)
  { // store to double precision locations
    endPoint->latest_newton_residual_d = T->latest_newton_residual_d;

    // find function residual and condition number
    findFunctionResidual_conditionNumber_d(&endPoint->function_residual_d, &endPoint->condition_number, &endPoint->PD_d, ED_d, eval_func_d);

    if (prec_before > 52)
    { // convert _mp to _d
      endPoint->t_val_at_latest_sample_point_d = mpf_get_d(endPoint->t_val_at_latest_sample_point_mp);
      endPoint->error_at_latest_sample_point_d = mpf_get_d(endPoint->error_at_latest_sample_point_mp);
    }
  }
  else
  { // set precision correctly for the values that need updated & store to mp locations
    mpf_clear(endPoint->function_residual_mp);
    mpf_init2(endPoint->function_residual_mp, endPoint->prec);

    mpf_clear(endPoint->latest_newton_residual_mp);
    mpf_init2(endPoint->latest_newton_residual_mp, endPoint->prec);

    // update values
    mpf_set(endPoint->latest_newton_residual_mp, T->latest_newton_residual_mp);
    findFunctionResidual_conditionNumber_mp(endPoint->function_residual_mp, &endPoint->condition_number, &endPoint->PD_mp, ED_mp, eval_func_mp);

    if (prec_before == 52)
    { // convert _d to _mp
      mpf_clear(endPoint->t_val_at_latest_sample_point_mp);
      mpf_init2(endPoint->t_val_at_latest_sample_point_mp, endPoint->prec);

      mpf_clear(endPoint->error_at_latest_sample_point_mp);
      mpf_init2(endPoint->error_at_latest_sample_point_mp, endPoint->prec);

      mpf_set_d(endPoint->t_val_at_latest_sample_point_mp, endPoint->t_val_at_latest_sample_point_d);
      mpf_set_d(endPoint->error_at_latest_sample_point_mp, endPoint->error_at_latest_sample_point_d);
    }
  }

  return;
}

/////////////////////// ZERO DIMENSIONAL SHARPENING ////////////////////////////////

int sharpen_d(int outputLevel, int sharpenDigits, double *latest_newton_residual, point_data_d *out, point_data_d *in, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int sharpen_mp(int outputLevel, int sharpenDigits, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, point_data_mp *out, point_data_mp *in, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int sharpen_amp(int outputLevel, int sharpenDigits, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int sharpen_zero_main(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: USED FOR ZERO DIMENSIONAL SHARPENING                   *
\***************************************************************/
{
  int retVal;
  double CN;

  retVal = determineRankDef(&CN, current_tol_d, current_tol_mp, tol_prec, in_d, in_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  if (retVal)
  { // Jacobian is singular so we just copy in to out (of the same precision) and return
    *prec_out = prec_in;

    if (prec_in == 52)
    { // copy using double precision
      point_data_cp_d(out_d, in_d);
    } 
    else
    { // verify precision
      point_data_mp tempPt;
      init_point_data_mp2(&tempPt, in_mp->point->size, prec_in);

      // copy in_mp to tempPt
      point_data_cp_mp(&tempPt, in_mp);

      // setup out_mp
      change_prec_point_data_mp(out_mp, prec_in);
      point_data_cp_mp(out_mp, in_mp);

      // clear tempPt
      clear_point_data_mp(&tempPt);
    }

    return retVal_sharpening_singular_endpoint;
  }
  // so we can assume that the Jacobian is non-singular

  // determine if we are using fixed precision or adaptive precision
  if (T->MPType == 0)
  { // using fixed double precision
    if (prec_in != 52)
    { // input was in higher precision so convert to double precision
      point_data_d tempPD_d;
      init_point_data_d(&tempPD_d, in_mp->point->size);
      convert_point_data_mp_to_d(&tempPD_d, in_mp);

      // try to sharpen using double precision
      retVal = sharpen_d(T->outputLevel, T->sharpenDigits, &T->latest_newton_residual_d, out_d, &tempPD_d, OUT, ED_d, eval_func_d);

      clear_point_data_d(&tempPD_d);
    }
    else
    { // try to sharpen using double precision
      retVal = sharpen_d(T->outputLevel, T->sharpenDigits, &T->latest_newton_residual_d, out_d, in_d, OUT, ED_d, eval_func_d);
    }
    *prec_out = 52;
  }
  else if (T->MPType == 1)
  { // using fixed multiprecision
    if (prec_in == 52)
    { // input was in double precision so convert to multiprecision
      point_data_mp tempPD_mp;
      init_point_data_mp2(&tempPD_mp, in_d->point->size, T->Precision);
      convert_point_data_d_to_mp(&tempPD_mp, in_d);

      // verify precision matches
      change_prec_point_data_mp(out_mp, T->Precision);

      // try to sharpen using fixed multiprecision
      retVal = sharpen_mp(T->outputLevel, T->sharpenDigits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, out_mp, &tempPD_mp, OUT, ED_mp, eval_func_mp);

      clear_point_data_mp(&tempPD_mp);
    }
    else
    { // input was in multiprecision

      // verify precision matches
      change_prec_point_data_mp(out_mp, T->Precision);

      // try to sharpen using fixed multiprecision
      retVal = sharpen_mp(T->outputLevel, T->sharpenDigits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, out_mp, in_mp, OUT, ED_mp, eval_func_mp);
    }
    *prec_out = T->Precision;
  }
  else if (T->MPType == 2)
  { // using adaptive precision
    retVal = sharpen_amp(T->outputLevel, T->sharpenDigits, &T->latest_newton_residual_d, T->latest_newton_residual_mp, current_tol_d, current_tol_mp, tol_prec, out_d, out_mp, prec_out, in_d, in_mp, prec_in, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }
  else
  {
    printf("ERROR: Invalid MPType (%d)!!\n", T->MPType);
    bexit(ERROR_CONFIGURATION);
  } 

  return retVal;
}

int sharpen_d(int outputLevel, int sharpenDigits, double *latest_newton_residual, point_data_d *out, point_data_d *in, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens using double precision                        *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, correct_digits = 0, maxIts = log10(sharpenDigits) / log10(2.0) + 5;
  double residual_old = 1e300;
  eval_struct_d e;
  init_eval_struct_d(e, 0, 0, 0);

  do // We want to keep refining until we hit the desired residual or the residual levels off.
  { 
    // perform a newton iteration
    retVal = newton_iteration_d(latest_newton_residual, 0, NULL, NULL, in->point, in->time, &e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in sharpen_d, but Bertini can still continue.\n");
      retVal = retVal_sharpening_failed;
      cont = 0;
    }
    else
    { 
      if (outputLevel > 1)
      {
        fprintf(OUT, "sharpening residual = %e\n", *latest_newton_residual);
      }

      // find the number of digits correct
      if (*latest_newton_residual == 0)
      { // no error - so all digits are correct
        correct_digits = 16;
      }
      else
      { // do the calculation
        correct_digits = (int) floor(-log10(*latest_newton_residual));
      }

      // check for convergence
      if (correct_digits >= sharpenDigits)
      { // converged so quit
        cont = 0;
        retVal = 0;
      }
      else if (*latest_newton_residual > residual_old)
      { // no better than the previous so quit
        cont = 0;
        retVal = retVal_sharpening_failed;
      }
      else
      { // otherwise, update the size of the residual and the number of iterations
        residual_old = *latest_newton_residual;
        its++;
  
        // check to see if we have done too many iterations
        if (its > maxIts)
        {
          cont = 0;
          retVal = retVal_sharpening_failed;
        }
      }
    }
  } while (cont);

  // copy the best approximation so far
  point_data_cp_d(out, in);

  clear_eval_struct_d(e);

  // determine if sharpening was successful 
  return retVal;
}

int sharpen_mp(int outputLevel, int sharpenDigits, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, point_data_mp *out, point_data_mp *in, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens using multi precision                         *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, correct_digits = 0, curr_prec = mpf_get_prec(in->point->coord[0].r), maxIts = log10(sharpenDigits) / log10(2.0) + 5;
  mpf_t residual_old, tempMPF;
  eval_struct_mp e;

  // initialize residual_old
  mpf_init_set_d(residual_old, 1e300);
  mpf_init(tempMPF);
  init_eval_struct_mp(e, 0, 0, 0);

  do // We want to keep refining until we hit the desired residual or the residual levels off.
  {
    // perform a newton iteration
    retVal = newton_iteration_mp(latest_newton_residual_mp, 0, NULL, NULL, in->point, in->time, &e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in sharpen_mp, but Bertini can still continue.\n");
      cont = 0;
      retVal = retVal_sharpening_failed;
    }
    else
    { // convert residual to double
      *latest_newton_residual_d = mpf_get_d(latest_newton_residual_mp);

      if (outputLevel > 1)
      {
        fprintf(OUT, "sharpening residual = "); mpf_out_str(OUT, 10, 6, latest_newton_residual_mp); fprintf(OUT, "\n");
      }

      // find the number of correct digits
      if (mpfr_zero_p(latest_newton_residual_mp))
      { // no error so all digits are correct
        correct_digits = (int) floor(curr_prec * log10(2) - 0.5);
      }
      else
      { // do the calculation
        mpfr_log10(tempMPF, latest_newton_residual_mp, __gmp_default_rounding_mode);
        mpf_neg(tempMPF, tempMPF);
        mpfr_floor(tempMPF, tempMPF);
        correct_digits = (int) mpf_get_si(tempMPF);
      }

      // check for convergence
      if (correct_digits >= sharpenDigits)
      { // converged so quit
        cont = 0;
        retVal = 0;
      }
      else if (mpfr_less_p(residual_old, latest_newton_residual_mp))
      { // no better than the previous so quit
        cont = 0;
        retVal = retVal_sharpening_failed;
      }
      else
      { // otherwise, update the size of the residual and the number of iterations
        mpf_set(residual_old, latest_newton_residual_mp);
        its++;

        // check to see if we have done too many iterations
        if (its > maxIts)
        {
          cont = 0;
          retVal = retVal_sharpening_failed;
        }
      }
    }
  } while (cont);

  // copy the best approximation so far
  point_data_cp_mp(out, in);

  // clear MP
  mpf_clear(residual_old); 
  mpf_clear(tempMPF);
  clear_eval_struct_mp(e);

  return retVal;
}

int sharpen_amp(int outputLevel, int sharpenDigits, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens using adpative precision                      *
\***************************************************************/
{
  int retVal = 0, its = 0, cont = 0, correct_digits = 0, prec_digits = 0, curr_prec = 0, maxIts = log10(sharpenDigits) / log10(2.0) + 5;
  int latest_newton_prec = mpf_get_prec(latest_newton_residual_mp);
  double residual_old_d  = 1e300;
  comp_mp residual_old_mp;
  point_data_mp tempPoint_mp;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  // MP things will be initialized to this, and so it needs to be atleast 64 bits
  curr_prec = MAX(prec_in, 64);

  // initialize e_d
  init_eval_struct_d(e_d, 0, 0, 0);

  if (tol_prec == 52)
  { // use current_tol_d to get an idea as to how many correct digits there should be currently
    if (current_tol_d <= 0)
    { // the tolerance is 0 so we set the number of correct digits based on the current precision
      if (prec_in == 52)
        correct_digits = 14;
      else
        correct_digits = (int) floor(prec_in * log10(2) - 1.5);
    }
    else
      correct_digits = (int) floor(-log10(current_tol_d));
  }
  else
  { // use current_tol_mp to get an idea as to how many correct digits there should be currently
    if (mpf_cmp_ui(current_tol_mp, 0) <= 0)
    { // the tolerance is 0 so we set the number of correct digits based on the current precision
      if (prec_in == 52)
        correct_digits = 14;
      else
        correct_digits = (int) floor(prec_in * log10(2) - 1.5);
    }
    else
    {
      mpf_t tempMPF;
      mpf_init2(tempMPF, tol_prec);
 
      mpfr_log10(tempMPF, current_tol_mp, __gmp_default_rounding_mode);
      mpf_neg(tempMPF, tempMPF);
      mpfr_floor(tempMPF, tempMPF);
      correct_digits = (int) mpf_get_si(tempMPF);

      mpf_clear(tempMPF);
    }
  }
  if (correct_digits < 1) // fail safe checking
    correct_digits = 1;

  if (prec_in == 52 && correct_digits < 7)
  { // loop using double precision until it no longer suffices

    // double precision can try to handle 14 digits
    prec_digits = 14;

    cont = 1;
    while (cont == 1)
    {
      // perform a newton iteration
      retVal = newton_iteration_d(latest_newton_residual_d, 0, NULL, NULL, in_d->point, in_d->time, &e_d, ED_d, eval_func_d);

      // check for newton success (i.e. if matrixSolve failed)
      if (retVal)
      {
        fprintf(OUT, "NOTE: matrixSolve failed in sharpen_amp, but Bertini can still continue.\n");
        cont = -1; // need to continue in MP
      }
      else
      {
        if (outputLevel > 1)
        {
          fprintf(OUT, "sharpening residual = %e\n", *latest_newton_residual_d);
        }

        // find the number of digits correct
        if (*latest_newton_residual_d == 0)
        { // no error - so all digits are correct
          correct_digits = prec_digits;
        }
        else
        { // do the calculation
          correct_digits = (int) floor(-log10(*latest_newton_residual_d));
          if (correct_digits < 1) // fail safe checking
            correct_digits = 1;

          if (correct_digits > prec_digits) // only double precision
            correct_digits = prec_digits;
        }

        // check for convergence
        if (correct_digits >= sharpenDigits)
        {
          cont = 0; // convergence so quit
          retVal = 0;
        } 
        else if (*latest_newton_residual_d > residual_old_d)
        { // no better than the previous so move to MP
          cont = -1; // continue on in MP
        }
        else
        { // otherwise, update the size of the residual and the number of iterations
          residual_old_d = *latest_newton_residual_d;
          its++;

          // check to see if we should continue
          if (its > maxIts)
          {
            cont = 0; // to many iterations so we need to quit
            retVal = retVal_sharpening_failed;
          }
          else if (2 * correct_digits > prec_digits)
          { // double precision no longer suffices so continue on in MP
            cont = -1;
          }
        }
      }
    }

    if (cont == -1)
    { // setup to continue on using MP in 64-bit precision
      
      // initialize residual_old_mp
      init_mp2(residual_old_mp, curr_prec);
      set_zero_mp(residual_old_mp);
      mpf_set_d(residual_old_mp->r, 1e300);

      // initialize tempPoint_mp
      init_point_data_mp2(&tempPoint_mp, in_d->point->size, curr_prec);
      convert_point_data_d_to_mp(&tempPoint_mp, in_d);

      // intialize e
      init_eval_struct_mp(e_mp, 0, 0, 0);
    }
    else
    { // we are finished so copy data over 
      *prec_out = 52;
      point_data_cp_d(out_d, in_d);
    }
  }
  else
  { // setup to use MP
    cont = -1; 

    // initialize residual_old_mp
    init_mp2(residual_old_mp, curr_prec);
    set_zero_mp(residual_old_mp);
    mpf_set_d(residual_old_mp->r, 1e300);

    // initialize tempPoint_mp
    init_point_data_mp2(&tempPoint_mp, 0, curr_prec);

    // initialize e
    init_eval_struct_mp(e_mp, 0, 0, 0);

    // copy to tempPoint_mp
    if (prec_in == 52)
    {
      convert_point_data_d_to_mp(&tempPoint_mp, in_d);
    }
    else
    {
      point_data_cp_mp(&tempPoint_mp, in_mp);
    }
  }

  if (cont == -1)
  { // loop over MP

    do
    { // find the precision needed for the next iteration - first calculate the '32-bit' multiplier
      prec_digits = (int) ceil((correct_digits + 1.5) / (32 * log10(2)));
      if (prec_digits < 1)
        prec_digits = 64;  // need precision atleast of 64 bits
      else
        prec_digits *= 4 * 32;
      // find the precision we should use
      if (prec_digits > 3 * curr_prec + 32)
        prec_digits = 3 * curr_prec + 32;

      // see if we have enough precision
      if (curr_prec < prec_digits)
      { // precision needs increased
        curr_prec = prec_digits;

        // set everything to this precision
        initMP(curr_prec);
        change_prec(ED_mp, curr_prec);
        mpf_set_prec(latest_newton_residual_mp, curr_prec);
 
        change_prec_mp2(residual_old_mp, curr_prec);
        set_zero_mp(residual_old_mp);
        mpf_set_d(residual_old_mp->r, 1e300); // since precision was increased, the old residual is meaningless 

        change_prec_point_data_mp(&tempPoint_mp, curr_prec);

        setprec_eval_struct_mp(e_mp, curr_prec);
      }

      // perform a newton iteration
      retVal = newton_iteration_mp(latest_newton_residual_mp, 0, NULL, NULL, tempPoint_mp.point, tempPoint_mp.time, &e_mp, ED_mp, eval_func_mp);

      // check for newton success (i.e. if matrixSolve failed)
      if (retVal)
      {
        fprintf(OUT, "NOTE: matrixSolve failed in sharpen_amp, but Bertini can still continue.\n");
        retVal = retVal_sharpening_failed;
        cont = 0; // need to quit
      }
      else
      { // convert residual to double
        *latest_newton_residual_d = mpf_get_d(latest_newton_residual_mp);

        if (outputLevel > 1)
        {
          fprintf(OUT, "sharpening residual = "); mpf_out_str(OUT, 10, 6, latest_newton_residual_mp); fprintf(OUT, "\n");
        }

        // find the number of correct digits
        if (mpfr_zero_p(latest_newton_residual_mp))
        { // no error so all digits are correct
          correct_digits = (int) floor(curr_prec * log10(2) - 0.5);
        }
        else
        { // do the calculation
          mpfr_log10(residual_old_mp->i, latest_newton_residual_mp, __gmp_default_rounding_mode);
          mpf_neg(residual_old_mp->i, residual_old_mp->i);
          mpfr_floor(residual_old_mp->i, residual_old_mp->i);
          correct_digits = (int) mpf_get_si(residual_old_mp->i);
          if (correct_digits < 1) // fail safe checking
            correct_digits = 1;

          if (correct_digits > prec_to_digits(curr_prec)) // resticted to current precision
            correct_digits = prec_to_digits(curr_prec); 
        }

        // check for convergence
        if (correct_digits >= sharpenDigits)
        {
          cont = 0; // convergence so quit
          retVal = 0;
        }
        else if (mpfr_less_p(residual_old_mp->r, latest_newton_residual_mp))
        { // no better than the previous so quit
          cont = 0; 
          retVal = retVal_sharpening_failed;
        }
        else
        { // otherwise, update the size of the residual and the number of iterations
          mpf_set(residual_old_mp->r, latest_newton_residual_mp);
          mpf_mul_ui(residual_old_mp->r, residual_old_mp->r, 10);
          its++;

          // check to see if we should continue
          if (its > maxIts)
          {
            cont = 0; // to many iterations so we need to quit
            retVal = retVal_sharpening_failed;
          }
        }
      }
    } while (cont);

    // calculate the exact precision required based on sharpenDigits
    prec_digits = (int) ceil((sharpenDigits + 0.5)/ (32 * log10(2)));
    if (prec_digits < 2) 
    { // set to use 64-bit
      prec_digits = 64;
    }
    else
      prec_digits *= 32;

    // record this precision
    *prec_out = prec_digits;

    // set everything to this precision
    initMP(*prec_out);
    change_prec(ED_mp, *prec_out);
    change_prec_point_mp(out_mp->point, *prec_out);
    change_prec_mp(out_mp->time, *prec_out);

    // set latest_newton_residual_mp back to its original precision
    mpf_set(residual_old_mp->r, latest_newton_residual_mp);
    mpf_set_prec(latest_newton_residual_mp, latest_newton_prec);
    mpf_set(latest_newton_residual_mp, residual_old_mp->r);

    // setup out_mp
    point_data_cp_mp(out_mp, &tempPoint_mp);

    // clear MP
    clear_mp(residual_old_mp);
    clear_point_data_mp(&tempPoint_mp);
    clear_eval_struct_mp(e_mp);
  }

  clear_eval_struct_d(e_d);

  return retVal;
}

