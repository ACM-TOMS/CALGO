// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

void deflate_for_sampling(prog_t **fullRankProg, int *fullRankProgInfo, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT);

void pos_dim_main(int trackType, int genType, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: main control function for positive dimensional         *
\***************************************************************/
{
  if (trackType == 1)
  { // numerical irreducible decomposition
    // genType: 0 - cascade, 1 - dim-by-dim, 2 - regen cascade
    numericalIrredDecomp(currentSeed, MPType, genType, my_id, num_processes, headnode);
  }
  else if (trackType == 2)
  { // sampling a component
    sampleComponent(currentSeed, MPType, 0, my_id, num_processes, headnode);
  }
  else if (trackType == 3)
  { // membership test
    membershipMain(currentSeed, MPType, my_id, num_processes, headnode);
  }
  else if (trackType == 4)
  { // print witness set
    printWitnessMain(currentSeed, MPType, my_id, num_processes, headnode);
  }
  else if (trackType == 5)
  { // witness projection
    witnessProjectionMain(currentSeed, MPType, my_id, num_processes, headnode);
  }
  else if (trackType == 6)
  { // isosingular stabilization test & witness generation
    witnessGeneration(currentSeed, MPType, startName, my_id, num_processes, headnode);
  }
  else if (trackType == 7)
  { // regeneration extension
    regenExtendMain(currentSeed, MPType, my_id, num_processes, headnode);
  }
  else
  {
    printf("ERROR: TrackType must be between 1 & 7!\n");
    bexit(ERROR_CONFIGURATION);
  }

  return;
}

////////// Sample a component ////////////////

void sampleComponent(unsigned int currentSeed, int MPType, int useSharpeningMenu, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds sample points on a given component               *
\***************************************************************/
{
  int rV, userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0;
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
    printf("ERROR: Parallel sampling is not implemented. Please use sequential version!\n");
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
  numIrredDecompOutput(&witnessSet, &T, 2, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 2

  if (useSharpeningMenu)
  { // display the menu that can change the number of digits
    rV = pos_dim_sharpening_menu(&T);
  }
  else
  { // setup to go to the sampling menu
    rV = 1;
  }

  if (rV)
  { // ask the user which component to sample and do the actual sampling
    sampleComponentMenu(&witnessSet, &T, pathMod);
  }

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void sampleComponentMenu(witness_t *W, tracker_config_t *T, int pathMod)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: displays a menu of components to sample and does the   *
* actual sampling                                               *
\***************************************************************/
{
  int i, j, codim_index, min_deg, max_deg, count, rV, dim_number, component_number, num_samples, outputType, selection_made, size_of_string = 255;
  int *degrees = NULL, *dim = (int *)bmalloc(W->num_codim * sizeof(int)), *codim_good = (int *)bmalloc(W->num_codim * sizeof(int));
  char ch, *tempStr = NULL, *outputFile = NULL;

  // find the dimension and good dimensions, and make sure one exists
  rV = 0;
  for (codim_index = 0; codim_index < W->num_codim; codim_index++)
  { // see if there are classified components for this codim
    codim_good[codim_index] = (W->codim[codim_index].num_components > 0 ? 1 : 0);

    if (codim_good[codim_index])
      rV = 1;

    // determine the dimension of this codim
    dim[codim_index] = W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
  }

  if (!rV)
  { // there are no classified components!!
    printf("\nThere are no classified components to sample!\n\n");
    free(dim);
    free(codim_good);
    return;
  }

  // so we have atleast one classified component
  do 
  { // initialize
    selection_made = 0;

    // print title
    printf("\n\n*************** Components to Sample ****************\n\n");

    // display a catalog of the available components in each codim
    for (codim_index = 0; codim_index < W->num_codim; codim_index++)
    { // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      // find the minimum and maximum degree
      max_deg = 0;
      min_deg = W->codim[codim_index].num_set + 1;
      for (i = 0; i < W->codim[codim_index].num_components; i++)
      { // find maximum degree
        if (max_deg < degrees[i])
          max_deg = degrees[i];

        // find minimum degree
        if (min_deg > degrees[i])
          min_deg = degrees[i];
      }

      if (W->codim[codim_index].num_components > 0)
      {
        printf("Dimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
        if (W->codim[codim_index].num_components == 1)
          printf("\n");
        else
          printf("s\n");
        printf("-----------------------------------------------------\n");

        // display the summary
        for (i = min_deg; i <= max_deg; i++)
        { // count the number that have degree == i
          count = 0;
          for (j = 0; j < W->codim[codim_index].num_components; j++)
            if (degrees[j] == i)
              count++;

          if (count > 0)
          { // display the number
            printf("   degree %d: %d component", i, count);
            if (count == 1)
              printf("\n");
            else
              printf("s\n");
          }
        }
        printf("\n");
      }
    }

    printf("\nPlease select a dimension to sample (-1 to quit): ");
    rV = scanf("%d", &dim_number);

    if (rV < 0)
    { // at EOF - so we need to quit
      dim_number = -1;
      selection_made = 1;
    }
    else
    { // we are not at EOF - flush the buffer
      do
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");
        selection_made = 0;
      }
      else if (dim_number == -1)
      { // dim_number is valid
        selection_made = 1;
      }
      else
      { // verify that the dimension is one of the good ones
        selection_made = 0;
        for (j = 0; j < W->num_codim; j++)
          if (codim_good[j] && dim[j] == dim_number)
          { // selection is valid
            selection_made = 1;
            break;
          }
          
        if (!selection_made) 
        { // dim_number is not valid
          printf("\nThe dimension %d is not valid!\n", dim_number);
        }
      }
    }
  } while (!selection_made);

  // so, either dim_number == -1 OR dim_number is one of the ones that have a classified component
  if (dim_number != -1)
  { // find the codim_index for the selected dimension
    codim_index = 0;
    for (i = 0; i < W->num_codim; i++)
    { // see if the dimension agrees
      if (dim[i] == dim_number)
      { // store the index and exit loop
        codim_index = i;
        break;
      }
    }

    do
    { // initialize
      selection_made = 0;

      // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      printf("\nDimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
      if (W->codim[codim_index].num_components == 1)
        printf("\n");
      else
        printf("s\n");
      printf("-----------------------------------------------------\n");

      // display the summary of components
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        printf("   component %d has degree %d\n", i, degrees[i]);
      printf("\n");

      printf("\nPlease select a component to sample (-1 to quit): ");
      rV = scanf("%d", &component_number);

      if (rV < 0)
      { // at EOF - so we need to quit
        component_number = -1;
        selection_made = 1;
      }
      else
      { // we are not at EOF - flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (rV == 0)
        { // invalid input
          printf("\nThe input was not read in correctly!\n");
          selection_made = 0;
        }
        else if (component_number != -1 && (component_number < 0 || component_number >= W->codim[codim_index].num_components))
        { // component_number is not valid
          printf("\nThe component %d is not valid!\n", component_number);
          selection_made = 0;
        }
        else
        { // component_number is valid
          selection_made = 1;
        }
      }
    } while (!selection_made);

    // so, either component_number == -1 OR 0 <= component_number < num_components
    if (0 <= component_number && component_number < W->codim[codim_index].num_components)
    { // determine how many samples to generate on the given component in the given dimension

      do
      { // initialize
        selection_made = 0;

        printf("How many points would you like to find on component %d of dimension %d (-1 to quit)?  ", component_number, dim_number);
        rV = scanf("%d", &num_samples);

        if (rV < 0)
        { // at EOF - so we need to quit
          num_samples = -1;
          selection_made = 1;
        }
        else
        { // we are not at EOF - flush the buffer
          do
          {
            ch = getchar();
          } while (ch != EOF && ch != '\n');

          if (rV == 0)
          { // invalid input
            printf("\nThe input was not read in correctly!\n");
            selection_made = 0;
          }
          else if (num_samples < -1)
          { // num_samples is not valid
            printf("\nThe number of samples %d is not valid!\n", num_samples);
            selection_made = 0;
          }
          else
          { // num_samples is valid
            selection_made = 1;
          }
        }
      } while (!selection_made);

      // so, either num_samples == -1, num_samples == 0 OR num_samples >= 1
      if (num_samples >= 1)
      { // determine what kind of output the user wants
        do
        { // initialize
          selection_made = 0;

          printf("Enter 0 to write the sample points to a file or enter 1 to display the sample points to the screen: ");
          rV = scanf("%d", &outputType);

          if (rV < 0)
          { // at EOF - so we need to quit
            outputType = 0;
            selection_made = 1;
          }
          else
          { // we are not at EOF - flush the buffer
            do
            {
              ch = getchar();
            } while (ch != EOF && ch != '\n');

            if (rV == 0)
            { // invalid input
              printf("\nThe input was not read in correctly!\n");
              selection_made = 0;
            }
            else if (outputType != 0 && outputType != 1)
            { // outputType is not valid
              printf("\nThe output type of %d is not valid!\n", outputType);
              selection_made = 0;
            }
            else
            { // outputType is valid
              selection_made = 1;
            }
          }
        } while (!selection_made);

        // determine if we need to setup a file name
        if (outputType == 0)
        { // read in a file name
          tempStr = (char *)bmalloc(((int) log10(size_of_string) + 10) * sizeof(char));
          snprintf(tempStr, size_of_string + 10, "%%%ds", size_of_string);
          outputFile = (char *)bmalloc((size_of_string + 1) * sizeof(char));
          for (i = 0; i <= size_of_string; i++)
            outputFile[i] = '\0';
          do 
          { // initialize
            selection_made = 0;

            printf("Enter the name of the file where to write the sample points (max of %d characters): ", size_of_string);
            rV = scanf(tempStr, outputFile);

            if (rV < 0)
            { // at EOF - setup to be 'sample_points'
              printf("\nThe file will be 'sample_points'.\n");
              sprintf(outputFile, "sample_points");
              selection_made = 1;
            }
            else
            { // we are not at EOF - flush the buffer
              do
              {
                ch = getchar();
              } while (ch != EOF && ch != '\n');

              // check to see if it is a valid file name
              FILE *TEMP = fopen(outputFile, "w");
              if (TEMP == NULL)
              { // not valid
                printf("\nThe name \"%s\" is not valid!\n\n", outputFile);
                selection_made = 0;
              }
              else
              { // valid - so close
                fclose(TEMP);
                selection_made = 1;
              }
            }
          } while (!selection_made);
        }

        // generate the samples
        int *samples_prec = (int *)bmalloc(num_samples * sizeof(int));
        point_d *samples_d = NULL;
        point_mp *samples_mp = NULL;

        if (T->MPType == 0 || T->MPType == 2)
        { // allocate samples_d
          samples_d = (point_d *)bmalloc(num_samples * sizeof(point_d));
        }
        if (T->MPType == 1 || T->MPType == 2)
        { // allocate samples_mp
          samples_mp = (point_mp *)bmalloc(num_samples * sizeof(point_mp));
        }

        num_samples = generateSamplePoints(samples_d, samples_mp, samples_prec, W, T, codim_index, component_number, num_samples, pathMod);

        if (outputType && num_samples > 0)
        { // print samples to screen using matlab notation
          if (num_samples == 1)
            printf("\nHere is the sample point:\n");
          else
            printf("\nHere are the sample points:\n");

          printSamplePoints(samples_d, samples_mp, samples_prec, num_samples, stdout, 1);
        }
        else if (num_samples > 0)
        { // print samples to 'sample_points' using 'real imag' notation
          FILE *OUT = fopen(outputFile, "w");
          printf("\nWriting the sample points to '%s'.\n\n", outputFile);
          printSamplePoints(samples_d, samples_mp, samples_prec, num_samples, OUT, 0);
          fclose(OUT);
        }

        // release memory
        for (i = num_samples - 1; i >= 0; i--)
          if (samples_prec[i] >= 64)
          { // clear samples_mp
            clear_point_mp(samples_mp[i]);
          }
        free(samples_prec);
        if (samples_d != NULL)
          free(samples_d);
        if (samples_mp != NULL)
          free(samples_mp);
      }
    }
  }

  free(tempStr);
  free(outputFile);
  free(degrees);
  free(dim);
  free(codim_good);

  return;
}

void printSamplePoints(point_d *samples_d, point_mp *samples_mp, int *samples_prec, int num_samples, FILE *fp, int useMatlab)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the samples that were generated to fp           *
\***************************************************************/
{
  int i, j;

  if (!useMatlab)
    fprintf(fp, "%d\n", num_samples);
  fprintf(fp, "\n");

  for (i = 0; i < num_samples; i++)
    if (useMatlab)
    { // print using Matlab format
      if (samples_prec[i] < 64)
        printVec_Matlab_d(fp, 0, samples_d[i]); 
      else
        printVec_Matlab_mp(fp, 0, samples_mp[i]); 
    }
    else
    { // print using 'real imag' format
      if (samples_prec[i] < 64)
      { // print using samples_d
        for (j = 0; j < samples_d[i]->size; j++)
        {
          print_d(fp, 0, &samples_d[i]->coord[j]);
          fprintf(fp, "\n");
        }
      }
      else
      { // print using samples_mp
        for (j = 0; j < samples_mp[i]->size; j++)
        {
          print_mp(fp, 0, &samples_mp[i]->coord[j]);
          fprintf(fp, "\n");
        }
      }
      fprintf(fp, "\n");
    }

  return;
}

int generateSamplePoints(point_d *samples_d, point_mp *samples_mp, int *samples_prec, witness_t *W, tracker_config_t *T, int codim_index, int component_number, int num_samples, int pathMod)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of samples actually generated           *
* NOTES: generates num_samples of points on the component       *
* component_number in dimension given by codim_index            *
\***************************************************************/
{
  int i, startPt, retVal, last_found, fullRankProgInfo, num_found = 0, num_attempts = 0, max_num_attempts_multiplier = 10, dim = W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
  endpoint_data_d endPt_d;
  endpoint_data_mp endPt_mp;
  endpoint_data_amp endPt_amp;
  prog_t *fullRankProg = (prog_t *)bmalloc(1 * sizeof(prog_t));
  membership_slice_moving_t sliceMover;
  FILE *OUT = fopen("output_sampling", "w"), *MIDOUT = fopen("midpath_data", "w");

  // first, find a point that is on the given component
  startPt = -1;
  for (i = 0; i < W->codim[codim_index].num_set; i++)
    if (W->codim[codim_index].component_nums[i] == component_number)
    { // store the startPt that is on the component and exit loop
      startPt = i;
      break;
    }

  // verify that we have a start point
  if (startPt < 0)
  { // no start point exist!
    printf("ERROR: No point on the given component exists!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup sliceMover for this codim
  basic_setup_slice_moving(&sliceMover, W, codim_index, T->MPType, T->AMP_max_prec);

  // now that we have a start point on the component, we can find a full rank SLP to generate the sample points
  deflate_for_sampling(&fullRankProg, &fullRankProgInfo, &endPt_d, &endPt_mp, &endPt_amp, &sliceMover, W, codim_index, startPt, T, OUT);

  // check for a successful deflation
  if (fullRankProgInfo == -1)
  { // deflation was not successful
    printf("\nBertini was unable to deflate component %d of dimension %d.\n", component_number, dim);
  }
  else
  { // deflation was successful
    printf("\nGenerating samples points on component %d of dimension %d: %d points to generate.\n", component_number, dim, num_samples);

    // generate random points on the component
    last_found = num_found - 1;
    while (num_found < num_samples && num_attempts < max_num_attempts_multiplier * num_samples)
    { // print the path number if needed
      if (pathMod > 0 && !(num_found % pathMod) && (last_found != num_found)) // last_found is used so that we do not keep constantly printing this message upon failure
      {
        printf("Generating %d of %d\n", num_found, num_samples);
        last_found = num_found;
      }

      // setup a random target vector
      setup_random_slice_moving(&sliceMover, T->MPType, T->AMP_max_prec);

      // setup a random gamma and finish setting up sliceMover
      final_setup_slice_moving(&sliceMover, fullRankProg, T->MPType, T->AMP_max_prec, 1); 

      // try to generate the next sample
      if (T->MPType == 0)
      { // track to the slice
        retVal = sampleTrack(samples_d[num_found], NULL, &samples_prec[num_found], &sliceMover, endPt_d.endPt, NULL, 52, num_found, T, OUT, MIDOUT);
      }
      else if (T->MPType == 1)
      { // track to the slice
        retVal = sampleTrack(NULL, samples_mp[num_found], &samples_prec[num_found], &sliceMover, NULL, endPt_mp.endPt, T->Precision, num_found, T, OUT, MIDOUT);
      }
      else
      { // track to the slice
        retVal = sampleTrack(samples_d[num_found], samples_mp[num_found], &samples_prec[num_found], &sliceMover, endPt_amp.endPt_d, endPt_amp.endPt_mp, endPt_amp.curr_prec, num_found, T, OUT, MIDOUT);
      }

      // check for success
      if (retVal == 0)
      { // this one is good - convert to dehomogeneous coordinates
        if (samples_prec[num_found] < 64)
        { // convert samples_d
          witnessFindDehom_d(samples_d[num_found], samples_d[num_found], W, codim_index);
        }
        else
        { // convert samples_mp
          witnessFindDehom_mp(samples_mp[num_found], samples_mp[num_found], W, codim_index, samples_prec[num_found]);
        }

        // increment the number found
        num_found++;
      }

      // increment the number of attempts
      num_attempts++;
    }

    // see if we have too many failures
    if (num_found < num_samples)
    {
      printf("\nWARNING: Using %d paths, Bertini was only able to generate %d of the requested %d samples.\n\n", num_attempts, num_found, num_samples);
    }
  }

  // close files
  fclose(OUT);
  fclose(MIDOUT);

  // clear deflation information
  clear_fullRankProg_endPt(&fullRankProg, fullRankProgInfo, &endPt_d, &endPt_mp, &endPt_amp, T->MPType); 
  free(fullRankProg);

  return num_found;
}

void setup_random_slice_moving(membership_slice_moving_t *sliceMover, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the sliceMover - random                          *
\***************************************************************/
{
  int i, size;

  if (MPType == 0)
  { // find the size
    size = sliceMover->B_d->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // seutp startSliceVec_d == 0
    sliceMover->startSliceVec_d->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
    }

    // setup targetSliceVec_d
    make_vec_random_d(sliceMover->targetSliceVec_d, size);
  }
  else if (MPType == 1)
  { // find the size
    size = sliceMover->B_mp->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // seutp startSliceVec_mp == 0
    sliceMover->startSliceVec_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
    }

    // setup targetSliceVec_mp
    make_vec_random_mp(sliceMover->targetSliceVec_mp, size);
  }
  else 
  { // find the size
    size = sliceMover->B_d->rows;

    // initialize slice vectors
    initialize_slice_moving_sliceVec(sliceMover, size, MPType);

    // seutp startSliceVec_d,_mp,_rat == 0
    sliceMover->startSliceVec_d->size = sliceMover->startSliceVec_mp->size = size;
    for (i = 0; i < size; i++)
    {
      set_zero_d(&sliceMover->startSliceVec_d->coord[i]);
      set_zero_mp(&sliceMover->startSliceVec_mp->coord[i]);
      set_zero_rat(sliceMover->startSliceVec_rat[i]);
    }

    // setup targetSliceVec_d, _mp, _rat
    make_vec_random_rat(sliceMover->targetSliceVec_d, sliceMover->targetSliceVec_mp, sliceMover->targetSliceVec_rat, size, sliceMover->curr_precision, max_prec, 0, 0);
  }

  return;
}

int sampleTrack(point_d endPt_d, point_mp endPt_mp, int *endPt_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - success, otherwise failure                 *
*   sample point - in 'dehom' coordinates                       *
* NOTES: tracks from the current path number to the random slice*
\***************************************************************/
{
  int retVal;
  endgame_data_t endPt;

  // initialize endPt
  if (T->MPType == 1)
  {
    init_endgame_data(&endPt, T->Precision);
  }
  else
  {
    init_endgame_data(&endPt, 64);
  }

  // move the slice - look to sharpen if needed
  retVal = slice_moving_track(&endPt, sliceMover, startPt_d, startPt_mp, startPt_prec, pathNum, 1, T, OUT, MIDOUT);

  // check for success
  if (retVal == 0)
  { // setup endPt
    *endPt_prec = endPt.prec;
    
    if (*endPt_prec < 64)
    { // setup endPt_d
      init_point_d(endPt_d, sliceMover->orig_variables);
      endPt.PD_d.point->size = sliceMover->orig_variables;
      point_cp_d(endPt_d, endPt.PD_d.point);
    }
    else
    { // setup endPt_mp
      init_point_mp2(endPt_mp, sliceMover->orig_variables, endPt.prec);
      endPt.PD_mp.point->size = sliceMover->orig_variables;
      point_cp_mp(endPt_mp, endPt.PD_mp.point);
    }
  }

  // clear endPt
  clear_endgame_data(&endPt);

  return retVal;
}

void deflate_for_sampling(prog_t **fullRankProg, int *fullRankProgInfo, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup fullRankProg and endPt to be reduced             *
\***************************************************************/
{
  int rV, input_prec, output_prec, newRandomizedProg, randomizedProgUsed = 0;
  point_data_d inputPD_d, outputPD_d;
  point_data_mp inputPD_mp, outputPD_mp;
  prog_t *randomizedProg = (prog_t *)bmalloc(1 * sizeof(prog_t));

  init_point_data_d(&inputPD_d, 0); init_point_data_d(&outputPD_d, 0);
  init_point_data_mp(&inputPD_mp, 0); init_point_data_mp(&outputPD_mp, 0);

  // setup the randomized SLP - does [I A]f
  if (W->codim[codim_index].A_rows > 0 && W->codim[codim_index].A_cols > 0)
  { // need to randomize
    randomize_SLP(randomizedProg, W->Prog, W->codim[codim_index].A_rat, W->codim[codim_index].A_rows, W->codim[codim_index].A_cols);
    // a new randomized prog was made
    newRandomizedProg = 1;
  }
  else
  { // point to original SLP
    randomizedProg = W->Prog;
    // a new randomized prog was not made
    newRandomizedProg = 0;
  }

  if (T->MPType == 0)
  { // copy all of the endPt data
    init_endpoint_data_d(endPt_d);
    endpoint_data_cp_d(endPt_d, &W->codim[codim_index].witnessPts_d[pathNum]);

    // see if we need to deflate
    if (W->codim[codim_index].witnessPt_types[pathNum] == NON_SINGULAR)
    { // simply point to the randomized prog
      *fullRankProg = randomizedProg;
      // setup Info
      if (newRandomizedProg && !randomizedProgUsed)
      { // we have a new SLP than has not been used before - so it will be cleared
        *fullRankProgInfo = 1;  
      }
      else
      { // this does not need cleared out
        *fullRankProgInfo = 0;
      }
      // we have used the randomized SLP
      randomizedProgUsed = 1;

      // setup the number of deflations needed == 0
      W->codim[codim_index].deflations_needed[pathNum] = 0;
    }
    else
    { // deflate
      prog_t *tempProg = (prog_t *)bmalloc(1 * sizeof(prog_t));

      // setup inputPD_d
      point_cp_d(inputPD_d.point, W->codim[codim_index].witnessPts_d[pathNum].endPt);
      set_d(inputPD_d.time, W->codim[codim_index].witnessPts_d[pathNum].finalT);
      input_prec = 52;

      // do the deflation
      rV = deflation(&W->codim[codim_index].deflations_needed[pathNum], tempProg, &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, endPt_d->corank, endPt_d->smallest_nonzero_SV, endPt_d->largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, W->codim[codim_index].witnessPts_d[pathNum].last_approx, NULL, 52, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[pathNum] - 1);

      // setup fullRankProg & endPt
      *fullRankProg = tempProg;
      point_cp_d(endPt_d->endPt, outputPD_d.point);

      // check for deflation success
      if (rV == 0)
        *fullRankProgInfo = 1; // deflated properly
      else
        *fullRankProgInfo = -1; // did not deflate properly

      tempProg = NULL; // pointed to by fullRankProg
    }
  }
  else if (T->MPType == 1)
  {
    init_endpoint_data_mp(endPt_mp);
    endpoint_data_cp_mp(endPt_mp, &W->codim[codim_index].witnessPts_mp[pathNum]);

    // see if we need to deflate
    if (W->codim[codim_index].witnessPt_types[pathNum] == NON_SINGULAR)
    { // simply point to the randomized prog
      *fullRankProg = randomizedProg;
      // setup Info
      if (newRandomizedProg && !randomizedProgUsed)
      { // we have a new SLP than has not been used before - so it will be cleared
        *fullRankProgInfo = 1;
      }
      else
      { // this does not need cleared out
        *fullRankProgInfo = 0;
      }
      // we have used the randomized SLP
      randomizedProgUsed = 1;

      // setup the number of deflations needed == 0
      W->codim[codim_index].deflations_needed[pathNum] = 0;
    }
    else
    { // deflate
      prog_t *tempProg = (prog_t *)bmalloc(1 * sizeof(prog_t));

      // setup inputPD_mp
      point_cp_mp(inputPD_mp.point, W->codim[codim_index].witnessPts_mp[pathNum].endPt);
      set_mp(inputPD_mp.time, W->codim[codim_index].witnessPts_mp[pathNum].finalT);
      input_prec = T->Precision;

      // do the deflation
      rV = deflation(&W->codim[codim_index].deflations_needed[pathNum], tempProg, &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, endPt_mp->corank, endPt_mp->smallest_nonzero_SV, endPt_mp->largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, NULL, W->codim[codim_index].witnessPts_mp[pathNum].last_approx, T->Precision, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[pathNum] - 1);

      // setup fullRankProg & endPt
      *fullRankProg = tempProg;
      point_cp_mp(endPt_mp->endPt, outputPD_mp.point);

      // check for deflation success
      if (rV == 0)
        *fullRankProgInfo = 1; // deflated properly
      else
        *fullRankProgInfo = -1; // did not deflate properly

      tempProg = NULL; // pointed to by fullRankProg
    }
  }
  else
  {
    init_endpoint_data_amp(endPt_amp, W->codim[codim_index].witnessPts_amp[pathNum].curr_prec, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_prec);
    endpoint_data_cp_amp(endPt_amp, &W->codim[codim_index].witnessPts_amp[pathNum]);

    // see if we need to deflate
    if (W->codim[codim_index].witnessPt_types[pathNum] == NON_SINGULAR)
    { // simply point to the randomized prog
      *fullRankProg = randomizedProg;
      // setup Info
      if (newRandomizedProg && !randomizedProgUsed)
      { // we have a new SLP than has not been used before - so it will be cleared
        *fullRankProgInfo = 1;
      }
      else
      { // this does not need cleared out
        *fullRankProgInfo = 0;
      }
      // we have used the randomized SLP
      randomizedProgUsed = 1;

      // setup the number of deflations needed == 0
      W->codim[codim_index].deflations_needed[pathNum] = 0;
    }
    else
    { // deflate
      prog_t *tempProg = (prog_t *)bmalloc(1 * sizeof(prog_t));

      // setup inputPD
      input_prec = W->codim[codim_index].witnessPts_amp[pathNum].curr_prec;
      if (input_prec < 64)
      {
        point_cp_d(inputPD_d.point, W->codim[codim_index].witnessPts_amp[pathNum].endPt_d);
        set_d(inputPD_d.time, W->codim[codim_index].witnessPts_amp[pathNum].finalT_d);
      }
      else
      {
        setprec_point_mp(inputPD_mp.point, input_prec);
        point_cp_mp(inputPD_mp.point, W->codim[codim_index].witnessPts_amp[pathNum].endPt_mp);

        setprec_mp(inputPD_mp.time, input_prec);
        set_mp(inputPD_mp.time, W->codim[codim_index].witnessPts_amp[pathNum].finalT_mp);
      }

      // do the deflation
      rV = deflation(&W->codim[codim_index].deflations_needed[pathNum], tempProg, &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, endPt_amp->corank, endPt_amp->smallest_nonzero_SV, endPt_amp->largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_d, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_mp, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_prec, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[pathNum] - 1);

      // setup fullRankProg & endPt
      *fullRankProg = tempProg;
      if (output_prec < 64)
      {
        point_cp_d(endPt_amp->endPt_d, outputPD_d.point);
      }
      else
      { // set precision correctly
        if (endPt_amp->curr_prec != output_prec)
        {
          setprec_point_mp(endPt_amp->endPt_mp, output_prec);
        }
        point_cp_mp(endPt_amp->endPt_mp, outputPD_mp.point);
      }

      // check for deflation success
      if (rV == 0)
        *fullRankProgInfo = 1; // deflated properly
      else
        *fullRankProgInfo = -1; // did not deflate properly

      tempProg = NULL; // pointed to by fullRankProg
    }
  }

  // clear data
  clear_point_data_d(&inputPD_d); clear_point_data_d(&outputPD_d);
  clear_point_data_mp(&inputPD_mp); clear_point_data_mp(&outputPD_mp);

  // see if we can clear randomizedProg
  if (newRandomizedProg && !randomizedProgUsed)
  { // clear out randomizedProg
    clearProg(randomizedProg, T->MPType, 0);
    free(randomizedProg);
  }
  else
  { // set to NULL
    randomizedProg = NULL;
  }

  return;
}

//////// copy dim-by-dim structure to witnessSuperset //////////

void dimbydim_copyWitness_clear(witness_t *witnessSuperset, codim_t *CD, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy CD to witnessSuperset and clear CD                *
\***************************************************************/
{
  int i, count;

  // target slices have not been initialized
  witnessSuperset->targetSliceInit = 0;

  // Prog
  witnessSuperset->Prog = CD->Prog;
  CD->Prog = NULL;

  // PPD
  cp_preproc_data(&witnessSuperset->PPD, &CD->PPD);
  preproc_data_clear(&CD->PPD);

  // orig_degrees, new_degrees, P
  witnessSuperset->orig_degrees = CD->orig_degrees;
  witnessSuperset->new_degrees = CD->new_degrees;
  witnessSuperset->P = CD->P;
  CD->orig_degrees = CD->new_degrees = CD->P = NULL;

  // system_rank, orig_variables, new_variables, curr_precision, num_funcs
  witnessSuperset->system_rank = CD->system_rank;
  witnessSuperset->orig_variables = CD->orig_variables;
  witnessSuperset->new_variables = CD->new_variables;
  witnessSuperset->curr_precision = CD->curr_precision;
  witnessSuperset->num_funcs = CD->num_funcs;

  // C (if needed)
  if (CD->orig_variables != CD->new_variables)
  {
    if (MPType == 0 || MPType == 2)
    { // copy and clear C_d
      init_mat_d(witnessSuperset->C_d, CD->C_d->rows, CD->C_d->cols);
      mat_cp_d(witnessSuperset->C_d, CD->C_d);
      clear_mat_d(CD->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // copy and clear C_mp
      init_mat_mp2(witnessSuperset->C_mp, CD->C_mp->rows, CD->C_mp->cols, witnessSuperset->curr_precision);
      mat_cp_mp(witnessSuperset->C_mp, CD->C_mp);
      clear_mat_mp(CD->C_mp);
    }

    if (MPType == 2)
    { // copy and clear C_rat
      witnessSuperset->C_rat = CD->C_rat;
      CD->C_rat = NULL;
    }
  }

  // gamma
  if (MPType == 0 || MPType == 2)
  { // copy and clear gamma_d
    set_d(witnessSuperset->gamma_d, CD->gamma_d);
    clear_d(CD->gamma_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy and clear gamma_mp
    init_mp2(witnessSuperset->gamma_mp, witnessSuperset->curr_precision);
    set_mp(witnessSuperset->gamma_mp, CD->gamma_mp);
    clear_mp(CD->gamma_mp);
  }

  if (MPType == 2)
  { // copy and clear gamma_rat
    mpq_init(witnessSuperset->gamma_rat[0]); mpq_init(witnessSuperset->gamma_rat[1]);
    mpq_set(witnessSuperset->gamma_rat[0], CD->gamma_rat[0]); mpq_set(witnessSuperset->gamma_rat[1], CD->gamma_rat[1]);
    mpq_clear(CD->gamma_rat[0]); mpq_clear(CD->gamma_rat[1]);
  }

  // count the number of codim that have witness points
  count = 0;
  for (i = 0; i < CD->num_codim; i++)
    if (CD->codim[i].num_superset > 0)
      count++;

  // set num_codim to count
  witnessSuperset->num_codim = count;

  // allocate memory for the 'num_codim' witness sets
  witnessSuperset->codim = (witnessCodim_t *)bmalloc(witnessSuperset->num_codim * sizeof(witnessCodim_t));

  // codim
  count = 0;
  for (i = 0; i < CD->num_codim; i++)
    if (CD->codim[i].num_superset > 0)
    { // copy and clear the ith codim
      dimbydim_copyWitness_clear_codim(&witnessSuperset->codim[count], &CD->codim[i], MPType, witnessSuperset->curr_precision, max_prec);
      count++;
    }
    else
    { // just clear this codim since there are no witness points
      clearCodimData(&CD->codim[i], MPType);
    }

  // free CD->codim
  free(CD->codim);

  return;
}

void dimbydim_copyWitness_clear_codim(witnessCodim_t *witCodim, codimData_t *CD, int MPType, int curr_prec, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy CD to witCodim and clear CD                       *
\***************************************************************/
{
  int i, j, num_paths = CD->num_paths;

  // codim, num_set, num_nonsing, num_sing
  witCodim->codim = CD->codim;
  witCodim->num_set = CD->num_superset;
  witCodim->num_nonsing = CD->num_nonsing;
  witCodim->num_sing = CD->num_sing;

  // W
  witCodim->W = CD->W;
  CD->W = NULL;

  // allocate witnessPt_types
  witCodim->witnessPt_types = (int *)bmalloc(witCodim->num_set * sizeof(int));

  // copy and clear A, H, homVarConst, B, & p, setup witnessPts & witnessPt_types, clear startPts, endPts & endPt_types
  if (MPType == 0)
  { // copy and clear A_d
    witCodim->A_rows = CD->A_d->rows;
    witCodim->A_cols = CD->A_d->cols;

    init_mat_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);
    mat_cp_d(witCodim->A_d, CD->A_d);
    clear_mat_d(CD->A_d);

    // setup A_rat
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpq_set_d(witCodim->A_rat[i][j][0], witCodim->A_d->entry[i][j].r);
        mpq_set_d(witCodim->A_rat[i][j][1], witCodim->A_d->entry[i][j].i);
      }

    // copy and clear H_d
    init_vec_d(witCodim->H_d, CD->H_d->size);
    vec_cp_d(witCodim->H_d, CD->H_d);
    clear_vec_d(CD->H_d);

    // copy and clear homVarConst_d
    set_d(witCodim->homVarConst_d, CD->homVarConst_d);
    clear_d(CD->homVarConst_d);

    // setup B_d & p_d
    if (CD->useIntrinsicSlice)
    { // we need to setup B_d & p_d to be the extrinsic version of B_d & p_d
      init_mat_d(witCodim->B_d, 0, 0);
      init_vec_d(witCodim->p_d, 0);
      intrinsicToExtrinsicSlices_d(witCodim->B_d, witCodim->p_d, CD->B_d, CD->p_d);
    }
    else
    { // copy B_d & p_d
      init_mat_d(witCodim->B_d, CD->B_d->rows, CD->B_d->cols);
      mat_cp_d(witCodim->B_d, CD->B_d);
      init_vec_d(witCodim->p_d, CD->p_d->size);
      vec_cp_d(witCodim->p_d, CD->p_d);
    }

    // allocate witnessPts_d
    witCodim->witnessPts_d = (endpoint_data_d *)bmalloc(witCodim->num_set * sizeof(endpoint_data_d));

    // setup witnessPts_d & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->endPt_types[j] == NON_SINGULAR || CD->endPt_types[j] == SINGULAR)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_d(&witCodim->witnessPts_d[i]);
        endpoint_data_cp_d(&witCodim->witnessPts_d[i], &CD->endPts_d[j]);

        if (CD->useIntrinsicSlice)
        { // convert endPts[j].endPt to extrinsic variables
          intrinsicToExtrinsic_d(witCodim->witnessPts_d[i].endPt, witCodim->witnessPts_d[i].endPt, CD->B_d, CD->p_d);
          // convert endPts[j].last_approx to extrinsic variables
          intrinsicToExtrinsic_d(witCodim->witnessPts_d[i].last_approx, witCodim->witnessPts_d[i].last_approx, CD->B_d, CD->p_d);
        }

        // copy type
        witCodim->witnessPt_types[i] = CD->endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_d[j] & endPts_d[j]
      if (CD->startPts_d != NULL)
        clear_point_d(CD->startPts_d[j]);
      clear_endpoint_data_d(&CD->endPts_d[j]);
    }

    // clear B_d & p_d
    clear_mat_d(CD->B_d);
    clear_vec_d(CD->p_d);

    // clear startPts_d, endPts_d & endPt_types
    if (CD->startPts_d != NULL)
      free(CD->startPts_d);
    free(CD->endPts_d);
    free(CD->endPt_types);
  }
  else if (MPType == 1)
  { // copy and clear A_mp
    init_mat_mp2(witCodim->A_mp, CD->A_mp->rows, CD->A_mp->cols, curr_prec);
    mat_cp_mp(witCodim->A_mp, CD->A_mp);
    clear_mat_mp(CD->A_mp);

    // setup A_rat
    witCodim->A_rows = witCodim->A_mp->rows;
    witCodim->A_cols = witCodim->A_mp->cols;
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpf_t_to_rat(witCodim->A_rat[i][j][0], witCodim->A_mp->entry[i][j].r);
        mpf_t_to_rat(witCodim->A_rat[i][j][1], witCodim->A_mp->entry[i][j].i);
      }

    // copy and clear H_mp
    init_vec_mp2(witCodim->H_mp, CD->H_mp->size, curr_prec);
    vec_cp_mp(witCodim->H_mp, CD->H_mp);
    clear_vec_mp(CD->H_mp);

    // copy and clear homVarConst_mp
    init_mp2(witCodim->homVarConst_mp, curr_prec);
    set_mp(witCodim->homVarConst_mp, CD->homVarConst_mp);
    clear_mp(CD->homVarConst_mp);

    // setup B_mp & p_mp
    if (CD->useIntrinsicSlice)
    { // we need to setup B_mp & p_mp to be the extrinsic version of B_mp & p_mp
      init_mat_mp2(witCodim->B_mp, 0, 0, curr_prec);
      init_vec_mp2(witCodim->p_mp, 0, curr_prec);
      intrinsicToExtrinsicSlices_mp(witCodim->B_mp, witCodim->p_mp, CD->B_mp, CD->p_mp);
    }
    else
    { // copy B_mp & p_mp
      init_mat_mp2(witCodim->B_mp, CD->B_mp->rows, CD->B_mp->cols, curr_prec);
      mat_cp_mp(witCodim->B_mp, CD->B_mp);
      init_vec_mp2(witCodim->p_mp, CD->p_mp->size, curr_prec);
      vec_cp_mp(witCodim->p_mp, CD->p_mp);
    }

    // allocate witnessPts_mp
    witCodim->witnessPts_mp = (endpoint_data_mp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_mp));

    // setup witnessPts_mp & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->endPt_types[j] == NON_SINGULAR || CD->endPt_types[j] == SINGULAR)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_mp(&witCodim->witnessPts_mp[i]);
        endpoint_data_cp_mp(&witCodim->witnessPts_mp[i], &CD->endPts_mp[j]);

        if (CD->useIntrinsicSlice)
        { // convert endPts[j].endPt to extrinsic variables
          intrinsicToExtrinsic_mp(witCodim->witnessPts_mp[i].endPt, witCodim->witnessPts_mp[i].endPt, CD->B_mp, CD->p_mp);
          // convert endPts[j].last_approx to extrinsic variables
          intrinsicToExtrinsic_mp(witCodim->witnessPts_mp[i].last_approx, witCodim->witnessPts_mp[i].last_approx, CD->B_mp, CD->p_mp);
        }
 
        // copy type
        witCodim->witnessPt_types[i] = CD->endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_mp[j] & endPts_mp[j]
      if (CD->startPts_mp != NULL)
        clear_point_mp(CD->startPts_mp[j]);
      clear_endpoint_data_mp(&CD->endPts_mp[j]);
    }

    // clear B_mp & p_mp
    clear_mat_mp(CD->B_mp);
    clear_vec_mp(CD->p_mp);

    // clear startPts_mp, endPts_mp & endPt_types
    if (CD->startPts_mp != NULL)
      free(CD->startPts_mp);
    free(CD->endPts_mp);
    free(CD->endPt_types);
  }
  else
  { // copy and clear A_d, A_mp, A_rat
    init_mat_d(witCodim->A_d, CD->A_d->rows, CD->A_d->cols);
    init_mat_mp2(witCodim->A_mp, CD->A_mp->rows, CD->A_mp->cols, curr_prec);
    mat_cp_d(witCodim->A_d, CD->A_d);
    mat_cp_mp(witCodim->A_mp, CD->A_mp);
    witCodim->A_rows = witCodim->A_d->rows;
    witCodim->A_cols = witCodim->A_d->cols;
    witCodim->A_rat = CD->A_rat;
    clear_mat_d(CD->A_d);
    clear_mat_mp(CD->A_mp);
    CD->A_rat = NULL;

    // copy and clear H_d, H_mp, H_rat
    init_vec_d(witCodim->H_d, CD->H_d->size);
    init_vec_mp2(witCodim->H_mp, CD->H_mp->size, curr_prec);
    vec_cp_d(witCodim->H_d, CD->H_d);
    vec_cp_mp(witCodim->H_mp, CD->H_mp);
    witCodim->H_rat = CD->H_rat;
    clear_vec_d(CD->H_d);
    clear_vec_mp(CD->H_mp);
    CD->H_rat = NULL;

    // copy and clear homVarConst_d, homVarConst_mp, homVarConst_rat
    init_d(witCodim->homVarConst_d);
    init_mp2(witCodim->homVarConst_mp, curr_prec);
    init_rat(witCodim->homVarConst_rat);
    set_d(witCodim->homVarConst_d, CD->homVarConst_d);
    set_mp(witCodim->homVarConst_mp, CD->homVarConst_mp);
    set_rat(witCodim->homVarConst_rat, CD->homVarConst_rat);
    clear_d(CD->homVarConst_d);
    clear_mp(CD->homVarConst_mp);
    clear_rat(CD->homVarConst_rat);

    // setup B_rat & p_rat
    init_mat_d(witCodim->B_d, 0, 0);
    init_mat_mp2(witCodim->B_mp, 0, 0, curr_prec);
    init_vec_d(witCodim->p_d, 0);
    init_vec_mp2(witCodim->p_mp, 0, curr_prec);
    if (CD->useIntrinsicSlice)
    { // we need to setup B_rat & p_rat to be the extrinsic version of B_rat & p_rat
      int B_rows, B_cols, p_size;

      intrinsicToExtrinsicSlices_rat(&witCodim->B_rat, &B_rows, &B_cols, &witCodim->p_rat, &p_size, CD->B_rat, CD->p_rat, CD->B_d->rows, CD->B_d->cols, CD->p_d->size, curr_prec, max_prec);

      // set sizes
      change_size_mat_d(witCodim->B_d, B_rows, B_cols);
      change_size_mat_mp(witCodim->B_mp, B_rows, B_cols);
      witCodim->B_d->rows = witCodim->B_mp->rows = B_rows;
      witCodim->B_d->cols = witCodim->B_mp->cols = B_cols;
      change_size_vec_d(witCodim->p_d, p_size);
      change_size_vec_mp(witCodim->p_mp, p_size);
      witCodim->p_d->size = witCodim->p_mp->size = p_size;
    }
    else
    { // point to B_rat & p_rat
      witCodim->B_rat = CD->B_rat;
      witCodim->p_rat = CD->p_rat;

      // set sizes
      change_size_mat_d(witCodim->B_d, CD->B_d->rows, CD->B_d->cols);
      change_size_mat_mp(witCodim->B_mp, CD->B_mp->rows, CD->B_mp->cols);
      witCodim->B_d->rows = witCodim->B_mp->rows = CD->B_d->rows;
      witCodim->B_d->cols = witCodim->B_mp->cols = CD->B_d->cols;
      change_size_vec_d(witCodim->p_d, CD->p_d->size);
      change_size_vec_mp(witCodim->p_mp, CD->p_mp->size);
      witCodim->p_d->size = witCodim->p_mp->size = CD->p_d->size;
    }

    // copy B_rat to B_d & B_mp
    for (i = 0; i < witCodim->B_d->rows; i++)
      for (j = 0; j < witCodim->B_d->cols; j++)
      {
        mpf_set_q(witCodim->B_mp->entry[i][j].r, witCodim->B_rat[i][j][0]);
        mpf_set_q(witCodim->B_mp->entry[i][j].i, witCodim->B_rat[i][j][1]);
        witCodim->B_d->entry[i][j].r = mpq_get_d(witCodim->B_rat[i][j][0]);
        witCodim->B_d->entry[i][j].i = mpq_get_d(witCodim->B_rat[i][j][1]);
      }

    // copy p_rat to p_d & p_mp
    for (i = 0; i < witCodim->p_d->size; i++)
    {
      mpf_set_q(witCodim->p_mp->coord[i].r, witCodim->p_rat[i][0]);
      mpf_set_q(witCodim->p_mp->coord[i].i, witCodim->p_rat[i][1]);
      witCodim->p_d->coord[i].r = mpq_get_d(witCodim->p_rat[i][0]);
      witCodim->p_d->coord[i].i = mpq_get_d(witCodim->p_rat[i][1]);
    }

    // allocate witnessPts_amp
    witCodim->witnessPts_amp = (endpoint_data_amp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_amp));

    // setup witnessPts_amp & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->endPt_types[j] == NON_SINGULAR || CD->endPt_types[j] == SINGULAR)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_amp(&witCodim->witnessPts_amp[i], CD->endPts_amp[j].curr_prec, CD->endPts_amp[j].last_approx_prec);
        endpoint_data_cp_amp(&witCodim->witnessPts_amp[i], &CD->endPts_amp[j]);

        if (CD->useIntrinsicSlice)
        { // convert endPts[j].endPt to extrinsic variables
          if (witCodim->witnessPts_amp[i].curr_prec < 64)
          {
            intrinsicToExtrinsic_d(witCodim->witnessPts_amp[i].endPt_d, witCodim->witnessPts_amp[i].endPt_d, CD->B_d, CD->p_d);
          }
          else
          {
            intrinsicToExtrinsic_mp(witCodim->witnessPts_amp[i].endPt_mp, witCodim->witnessPts_amp[i].endPt_mp, CD->B_mp, CD->p_mp);
          }

          // convert endPts[j].last_approx to extrinsic variables
          if (witCodim->witnessPts_amp[i].last_approx_prec < 64)
          {
            intrinsicToExtrinsic_d(witCodim->witnessPts_amp[i].last_approx_d, witCodim->witnessPts_amp[i].last_approx_d, CD->B_d, CD->p_d);
          }
          else
          {
            intrinsicToExtrinsic_mp(witCodim->witnessPts_amp[i].last_approx_mp, witCodim->witnessPts_amp[i].last_approx_mp, CD->B_mp, CD->p_mp);
          }
        }

        // copy type
        witCodim->witnessPt_types[i] = CD->endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_d[j] & endPts_amp[j]
      if (CD->startPts_d != NULL)
        clear_point_d(CD->startPts_d[j]);
      clear_endpoint_data_amp(&CD->endPts_amp[j]);
    }

    // clear B_d, B_mp & B_rat, p_d, p_mp, & p_rat
    if (CD->useIntrinsicSlice)
    { // B_rat & p_rat needs cleared
      clear_mat(CD->B_d, CD->B_mp, CD->B_rat, MPType);
      clear_vec(CD->p_d, CD->p_mp, CD->p_rat, MPType);
    }
    else
    { // B_rat & p_rat need set to NULL
      clear_mat_d(CD->B_d);
      clear_mat_mp(CD->B_mp);
      CD->B_rat = NULL;

      clear_vec_d(CD->p_d);
      clear_vec_mp(CD->p_mp);
      CD->p_rat = NULL;
    }

    // clear startPts_d, endPts_amp & endPt_types
    if (CD->startPts_d != NULL)
      free(CD->startPts_d);
    free(CD->endPts_amp);
    free(CD->endPt_types);
  }

  // intialize other values in witCodim
  witCodim->num_components = 0;
  witCodim->multiplicities = witCodim->component_nums = witCodim->deflations_needed = NULL;

  // NULL out the pointers
  CD->startPts_d = NULL;
  CD->startPts_mp = NULL;
  CD->endPts_d = NULL;
  CD->endPts_mp = NULL;
  CD->endPts_amp = NULL;

  return;
}

//////// copy casacde structure to witnessSuperset //////////

int cascade_copyWitness_clear(witness_t *witnessSuperset, cascade_t *CD, int MPType, int max_prec, int specificCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy CD to witnessSuperset and clear CD                *
\***************************************************************/
{
  int i, count, retVal = 0;

  // target slices have not been initialized
  witnessSuperset->targetSliceInit = 0;

  // Prog
  witnessSuperset->Prog = CD->Prog;
  CD->Prog = NULL;

  // PPD
  cp_preproc_data(&witnessSuperset->PPD, &CD->PPD);
  preproc_data_clear(&CD->PPD);

  // orig_degrees, new_degrees, P
  witnessSuperset->orig_degrees = CD->orig_degrees;
  witnessSuperset->new_degrees = CD->new_degrees;
  witnessSuperset->P = CD->P;

  // system_rank, orig_variables, new_variables, curr_precision, num_funcs
  witnessSuperset->system_rank = CD->system_rank;
  witnessSuperset->orig_variables = CD->orig_variables;
  witnessSuperset->new_variables = CD->new_variables;
  witnessSuperset->curr_precision = CD->curr_precision;
  witnessSuperset->num_funcs = CD->num_funcs;

  // C (if needed)
  if (CD->orig_variables != CD->new_variables)
  {
    if (MPType == 0 || MPType == 2)
    { // copy and clear C_d
      init_mat_d(witnessSuperset->C_d, CD->C_d->rows, CD->C_d->cols);
      mat_cp_d(witnessSuperset->C_d, CD->C_d);
      clear_mat_d(CD->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // copy and clear C_mp
      init_mat_mp2(witnessSuperset->C_mp, CD->C_mp->rows, CD->C_mp->cols, witnessSuperset->curr_precision);
      mat_cp_mp(witnessSuperset->C_mp, CD->C_mp);
      clear_mat_mp(CD->C_mp);
    }

    if (MPType == 2)
    { // copy and clear C_rat
      witnessSuperset->C_rat = CD->C_rat;
      CD->C_rat = NULL;
    }
  }

  // gamma
  if (MPType == 0 || MPType == 2)
  { // copy and clear gamma_d
    set_d(witnessSuperset->gamma_d, CD->gamma_d);
    clear_d(CD->gamma_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy and clear gamma_mp
    init_mp2(witnessSuperset->gamma_mp, witnessSuperset->curr_precision);
    set_mp(witnessSuperset->gamma_mp, CD->gamma_mp);
    clear_mp(CD->gamma_mp);
  }

  if (MPType == 2)
  { // copy and clear gamma_rat
    init_rat(witnessSuperset->gamma_rat);
    set_rat(witnessSuperset->gamma_rat, CD->gamma_rat);
    clear_rat(CD->gamma_rat);
  }

  // count the number of codim that have witness points that we care about
  count = 0;
  retVal = 1;
  for (i = 0; i < CD->num_codim; i++)
    if (CD->codim[i].num_superset > 0)
    {
      if (specificCodim == 0 || CD->codim[i].codim == specificCodim)
        count++;
      else if (specificCodim > 0)
        retVal = 0;
    }

  // set num_codim to count
  witnessSuperset->num_codim = count;

  // allocate memory for the 'num_codim' witness sets
  witnessSuperset->codim = (witnessCodim_t *)bmalloc(witnessSuperset->num_codim * sizeof(witnessCodim_t));

  // setup codim
  count = 0;
  for (i = 0; i < CD->num_codim; i++)
    if (CD->codim[i].num_superset > 0 && (specificCodim == 0 || CD->codim[i].codim == specificCodim))
    { // copy and clear the ith codim
      cascade_copyWitness_clear_codim(&witnessSuperset->codim[count], CD, i, MPType, witnessSuperset->curr_precision, max_prec);
      count++;
    }
    else
    { // just clear this codim since we are not interested in it 
      clearCascadeCodim(CD, i, MPType);
    }

  // clear the rest of the structures (orig_degrees, new_degrees, P, A, W, H, homVarConst, R, T, B, p, W_prime)

  // clear orig_degrees, new_degrees, P
  CD->orig_degrees = CD->new_degrees = CD->P = NULL;

  // clear A, if needed
  if (CD->num_funcs != CD->system_rank)
  { // clear A
    clear_mat(CD->A_d, CD->A_mp, CD->A_rat, MPType);
  }

  // clear W, if needed
  if (CD->num_funcs != CD->system_rank)
  { // clear W
    for (i = CD->system_rank - 1; i >= 0; i--)
      free(CD->W[i]);
    free(CD->W);
  }

  // clear H
  clear_vec(CD->H_d, CD->H_mp, CD->H_rat, MPType);

  // clear homVarConst
  clear_d_mp_rat(CD->homVarConst_d, CD->homVarConst_mp, CD->homVarConst_rat, MPType);

  // clear R
  clear_mat(CD->R_d, CD->R_mp, CD->R_rat, MPType);

  // clear T
  if (MPType == 0 || MPType == 2)
  { // clear T_d
    clear_vec_d(CD->T_d);
  }
  if (MPType == 1 || MPType == 2)
  { // clear T_mp
    clear_vec_mp(CD->T_mp);
  }

  // clear B
  clear_mat(CD->B_d, CD->B_mp, CD->B_rat, MPType);

  // clear p
  clear_vec(CD->p_d, CD->p_mp, CD->p_rat, MPType);

  // free W_prime
  free(CD->W_prime);

  // free CD->codim
  free(CD->codim);

  return retVal;
}

void cascade_copyWitness_clear_codim(witnessCodim_t *witCodim, cascade_t *CD, int codim_index, int MPType, int curr_prec, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy CD to witCodim and clear CD[codim_index]          *
\***************************************************************/
{
  int i, j, num_paths = CD->codim[codim_index].num_paths, codim = CD->codim[codim_index].codim, num_funcs = CD->num_funcs;

  // codim, num_set, num_nonsing, num_sing
  witCodim->codim = CD->codim[codim_index].codim;
  witCodim->num_set = CD->codim[codim_index].num_superset;
  witCodim->num_nonsing = CD->codim[codim_index].num_nonsing;
  witCodim->num_sing = CD->codim[codim_index].num_sing;

  // setup W
  witCodim->W = (int **)bmalloc(codim * sizeof(int *));
  for (i = 0; i < codim; i++)
  {
    witCodim->W[i] = (int *)bmalloc((num_funcs - codim) * sizeof(int));
    for (j = 0; j < num_funcs - codim; j++)
      witCodim->W[i][j] = CD->new_degrees[i] - CD->new_degrees[j + codim];
  }

  // setup A
  witCodim->A_rows = codim;
  witCodim->A_cols = num_funcs - codim;
  if (MPType == 0)
  { // setup A_d
    init_mat_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);
    make_matrix_random_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);

    // setup A_rat
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpq_set_d(witCodim->A_rat[i][j][0], witCodim->A_d->entry[i][j].r);
        mpq_set_d(witCodim->A_rat[i][j][1], witCodim->A_d->entry[i][j].i);
      }
  }
  else if (MPType == 1)
  { // setup A_mp
    init_mat_mp2(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);
    make_matrix_random_mp(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);

    // setup A_rat
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpf_t_to_rat(witCodim->A_rat[i][j][0], witCodim->A_mp->entry[i][j].r);
        mpf_t_to_rat(witCodim->A_rat[i][j][1], witCodim->A_mp->entry[i][j].i);
      }
  }
  else
  { // setup A_d, A_mp, A_rat

    // allocate for A_d, A_mp, A_rat
    init_mat_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);
    init_mat_mp2(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    // setup A_d, A_mp, A_rat
    make_matrix_random_rat(witCodim->A_d, witCodim->A_mp, witCodim->A_rat, witCodim->A_rows, witCodim->A_cols, curr_prec, max_prec, 0, 0);
  }

  // copy over H & homVarConst
  if (MPType == 0 || MPType == 2)
  { // setup H_d & homVarConst_d
    init_vec_d(witCodim->H_d, CD->H_d->size);
    vec_cp_d(witCodim->H_d, CD->H_d);
    set_d(witCodim->homVarConst_d, CD->homVarConst_d);
  }
  if (MPType == 1 || MPType == 2)
  { // setup H_mp & homVarConst_mp
    init_vec_mp2(witCodim->H_mp, CD->H_mp->size, curr_prec);
    init_mp2(witCodim->homVarConst_mp, curr_prec);
    vec_cp_mp(witCodim->H_mp, CD->H_mp);
    set_mp(witCodim->homVarConst_mp, CD->homVarConst_mp);
  }
  if (MPType == 2)
  { // setup H_rat & homVarConst_rat

    // allocate, initialize and set H_rat
    init_vec_rat(witCodim->H_rat, witCodim->H_mp->size);
    for (i = 0; i < witCodim->H_mp->size; i++)
    {
      set_rat(witCodim->H_rat[i], CD->H_rat[i]);
    }

    // initialize and set homVarConst_rat
    init_rat(witCodim->homVarConst_rat);
    set_rat(witCodim->homVarConst_rat, CD->homVarConst_rat);
  }

  // now the extrinsic slices B & p
  if (MPType == 0 || MPType == 2)
  { // setup B_d & p_d
    init_mat_d(witCodim->B_d, CD->new_variables - codim - 1, CD->new_variables);
    init_vec_d(witCodim->p_d, CD->new_variables);
    witCodim->B_d->rows = CD->new_variables - codim - 1;
    witCodim->B_d->cols = witCodim->p_d->size = CD->new_variables;
    for (j = 0; j < CD->new_variables; j++)
    { // copy p_d[j]
      set_d(&witCodim->p_d->coord[j], &CD->p_d->coord[j]);
      for (i = 0; i < witCodim->B_d->rows; i++)
      { // copy B_d[i][j]
        set_d(&witCodim->B_d->entry[i][j], &CD->B_d->entry[i][j]);
      }
    }
  }
  if (MPType == 1 || MPType == 2)
  { // setup B_mp & p_mp
    init_mat_mp2(witCodim->B_mp, CD->new_variables - codim - 1, CD->new_variables, curr_prec);
    init_vec_mp2(witCodim->p_mp, CD->new_variables, curr_prec);
    witCodim->B_mp->rows = CD->new_variables - codim - 1;
    witCodim->B_mp->cols = witCodim->p_mp->size = CD->new_variables;
    for (j = 0; j < CD->new_variables; j++)
    { // copy p_mp[j]
      set_mp(&witCodim->p_mp->coord[j], &CD->p_mp->coord[j]);
      for (i = 0; i < witCodim->B_mp->rows; i++)
      { // copy B_mp[i][j]
        set_mp(&witCodim->B_mp->entry[i][j], &CD->B_mp->entry[i][j]);
      }
    }
  }
  if (MPType == 2)
  { // allocate & initialize B_rat and p_rat
    init_mat_rat(witCodim->B_rat, witCodim->B_mp->rows, witCodim->B_mp->cols);
    init_vec_rat(witCodim->p_rat, witCodim->p_mp->size);

    // setup p_rat & B_rat
    for (j = 0; j < CD->new_variables; j++)
    { // setup p_rat[j]
      mpq_set(witCodim->p_rat[j][0], CD->p_rat[j][0]);
      mpq_set(witCodim->p_rat[j][1], CD->p_rat[j][1]);

      for (i = 0; i < witCodim->B_mp->rows; i++)
      { // copy B_rat[i][j]
        mpq_set(witCodim->B_rat[i][j][0], CD->B_rat[i][j][0]);
        mpq_set(witCodim->B_rat[i][j][1], CD->B_rat[i][j][1]);
      }
    }
  }

  // allocate witnessPt_types
  witCodim->witnessPt_types = (int *)bmalloc(witCodim->num_set * sizeof(int));

  // setup witnessPts & witnessPt_types and clear startPts, endPts & endPt_types
  if (MPType == 0)
  { // allocate witnessPts_d
    witCodim->witnessPts_d = (endpoint_data_d *)bmalloc(witCodim->num_set * sizeof(endpoint_data_d));

    // setup witnessPts_d & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_NONSING || CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_SING)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_d(&witCodim->witnessPts_d[i]);
        endpoint_data_cp_d(&witCodim->witnessPts_d[i], &CD->codim[codim_index].endPts_d[j]);

        // copy type
        witCodim->witnessPt_types[i] = CD->codim[codim_index].endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_d[j] & endPts_d[j]
      if (CD->codim[codim_index].startPts_d != NULL)
        clear_point_d(CD->codim[codim_index].startPts_d[j]);
      clear_endpoint_data_d(&CD->codim[codim_index].endPts_d[j]);
    }

    // clear startPts_d, endPts_d & endPt_types
    if (CD->codim[codim_index].startPts_d != NULL)
      free(CD->codim[codim_index].startPts_d);
    free(CD->codim[codim_index].endPts_d);
    free(CD->codim[codim_index].endPt_types);
  }
  else if (MPType == 1)
  { // allocate witnessPts_mp
    witCodim->witnessPts_mp = (endpoint_data_mp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_mp));

    // setup witnessPts_mp & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_NONSING || CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_SING)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_mp(&witCodim->witnessPts_mp[i]);
        endpoint_data_cp_mp(&witCodim->witnessPts_mp[i], &CD->codim[codim_index].endPts_mp[j]);

        // copy type
        witCodim->witnessPt_types[i] = CD->codim[codim_index].endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_mp[j] & endPts_mp[j]
      if (CD->codim[codim_index].startPts_mp != NULL)
        clear_point_mp(CD->codim[codim_index].startPts_mp[j]);
      clear_endpoint_data_mp(&CD->codim[codim_index].endPts_mp[j]);
    }

    // clear startPts_mp, endPts_mp & endPt_types
    if (CD->codim[codim_index].startPts_mp != NULL)
      free(CD->codim[codim_index].startPts_mp);
    free(CD->codim[codim_index].endPts_mp);
    free(CD->codim[codim_index].endPt_types);
  }
  else
  { // allocate witnessPts_amp
    witCodim->witnessPts_amp = (endpoint_data_amp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_amp));

    // setup witnessPts_amp & witnessPt_types
    i = 0;
    for (j = 0; j < num_paths; j++)
    { // find the next endpoint that is in the witness set
      if (CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_NONSING || CD->codim[codim_index].endPt_types[j] == SOLUTION_AND_SING)
      { // copy endPts[j] to witnessPts[i]
        init_endpoint_data_amp(&witCodim->witnessPts_amp[i], CD->codim[codim_index].endPts_amp[j].curr_prec, CD->codim[codim_index].endPts_amp[j].last_approx_prec);
        endpoint_data_cp_amp(&witCodim->witnessPts_amp[i], &CD->codim[codim_index].endPts_amp[j]);

        // copy type
        witCodim->witnessPt_types[i] = CD->codim[codim_index].endPt_types[j];

        // increment i
        i++;
      }

      // clear startPts_d[j] & endPts_amp[j]
      if (CD->codim[codim_index].startPts_d != NULL)
        clear_point_d(CD->codim[codim_index].startPts_d[j]);
      clear_endpoint_data_amp(&CD->codim[codim_index].endPts_amp[j]);
    }

    // clear startPts_d, endPts_amp & endPt_types
    if (CD->codim[codim_index].startPts_d != NULL)
      free(CD->codim[codim_index].startPts_d);
    free(CD->codim[codim_index].endPts_amp);
    free(CD->codim[codim_index].endPt_types);
  }

  // intialize other values in witCodim
  witCodim->num_components = 0;
  witCodim->multiplicities = witCodim->component_nums = witCodim->deflations_needed = NULL;

  // NULL out the pointers
  CD->codim[codim_index].startPts_d = NULL;
  CD->codim[codim_index].startPts_mp = NULL;
  CD->codim[codim_index].endPts_d = NULL;
  CD->codim[codim_index].endPts_mp = NULL;
  CD->codim[codim_index].endPts_amp = NULL;

  return;
}

///////////// change precision //////////////

int change_witness_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for witness_t                         *
\***************************************************************/
{
  int i, j, k;

  // cast ED as W
  witness_t *W = (witness_t *)ED;

  // set the SLP to the correct precision
  W->Prog->precision = prec;

  if (prec != W->curr_precision)
  { // need to change the precision
    W->curr_precision = prec;

    // change the precision for gamma
    setprec_mp(W->gamma_mp, prec);
    mpf_set_q(W->gamma_mp->r, W->gamma_rat[0]);
    mpf_set_q(W->gamma_mp->i, W->gamma_rat[1]);

    // change the precision for C, if needed
    if (W->new_variables != W->orig_variables)
    {
      for (i = 0; i < W->C_mp->rows; i++)
        for (j = 0; j < W->C_mp->cols; j++)
        {
          setprec_mp(&W->C_mp->entry[i][j], prec);
          mpf_set_q(W->C_mp->entry[i][j].r, W->C_rat[i][j][0]);
          mpf_set_q(W->C_mp->entry[i][j].i, W->C_rat[i][j][1]);
        }
    }

    // change the precision for targetSliceMat & Vec, if needed
    if (W->targetSliceInit)
    {
      for (i = 0; i < W->targetSliceMat_mp->rows; i++)
        for (j = 0; j < W->targetSliceMat_mp->cols; j++)
        {
          setprec_mp(&W->targetSliceMat_mp->entry[i][j], prec);
          mpf_set_q(W->targetSliceMat_mp->entry[i][j].r, W->targetSliceMat_rat[i][j][0]);
          mpf_set_q(W->targetSliceMat_mp->entry[i][j].i, W->targetSliceMat_rat[i][j][1]);
        }

      for (i = 0; i < W->targetSliceVec_mp->size; i++)
      {
        setprec_mp(&W->targetSliceVec_mp->coord[i], prec);
        mpf_set_q(W->targetSliceVec_mp->coord[i].r, W->targetSliceVec_rat[i][0]);
        mpf_set_q(W->targetSliceVec_mp->coord[i].i, W->targetSliceVec_rat[i][1]);
      }
    }

    // change the precision for all of A_mp, B_mp & p_mp in codim
    for (k = 0; k < W->num_codim; k++)
    { // change the precision for H
      for (i = 0; i < W->codim[k].H_mp->size; i++)
      {
        setprec_mp(&W->codim[k].H_mp->coord[i], prec);
        mpf_set_q(W->codim[k].H_mp->coord[i].r, W->codim[k].H_rat[i][0]);
        mpf_set_q(W->codim[k].H_mp->coord[i].i, W->codim[k].H_rat[i][1]);
      }

      // change the precision for homVarConst
      setprec_mp(W->codim[k].homVarConst_mp, prec);
      mpf_set_q(W->codim[k].homVarConst_mp->r, W->codim[k].homVarConst_rat[0]);
      mpf_set_q(W->codim[k].homVarConst_mp->i, W->codim[k].homVarConst_rat[1]);

      // change the precision for A
      for (i = 0; i < W->codim[k].A_mp->rows; i++)
        for (j = 0; j < W->codim[k].A_mp->cols; j++)
        {
          setprec_mp(&W->codim[k].A_mp->entry[i][j], prec);
          mpf_set_q(W->codim[k].A_mp->entry[i][j].r, W->codim[k].A_rat[i][j][0]);
          mpf_set_q(W->codim[k].A_mp->entry[i][j].i, W->codim[k].A_rat[i][j][1]);
        }

      // change the precision for B
      for (i = 0; i < W->codim[k].B_mp->rows; i++)
        for (j = 0; j < W->codim[k].B_mp->cols; j++)
        {
          setprec_mp(&W->codim[k].B_mp->entry[i][j], prec);
          mpf_set_q(W->codim[k].B_mp->entry[i][j].r, W->codim[k].B_rat[i][j][0]);
          mpf_set_q(W->codim[k].B_mp->entry[i][j].i, W->codim[k].B_rat[i][j][1]);
        }

      // change the precision for p
      for (i = 0; i < W->codim[k].p_mp->size; i++)
      {
        setprec_mp(&W->codim[k].p_mp->coord[i], prec);
        mpf_set_q(W->codim[k].p_mp->coord[i].r, W->codim[k].p_rat[i][0]);
        mpf_set_q(W->codim[k].p_mp->coord[i].i, W->codim[k].p_rat[i][1]);
      }
    }
  }

  return 0;
}

/////////// main clearing function ///////////

void witness_clear(witness_t *W, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears W                                               *
\***************************************************************/
{
  int i;

  // first clear all of the codimensions
  for (i = W->num_codim - 1; i >= 0; i--)
    witness_clear_codim(&W->codim[i], MPType); 
  free(W->codim);

  // now clear the other structures

  // clear Prog
  if (W->Prog != NULL)
  {
    clearProg(W->Prog, MPType, 1); // clear evaluation structures since we are done
    free(W->Prog);
  }

  // clear PPD
  preproc_data_clear(&W->PPD);

  // clear orig_degrees
  if (W->orig_degrees != NULL)
    free(W->orig_degrees);

  // clear new_degrees
  if (W->new_degrees != NULL)
    free(W->new_degrees);

  // clear P
  if (W->P != NULL)
    free(W->P);

  // clear C, if needed
  if (W->orig_variables != W->new_variables)
  { // clear C
    clear_mat(W->C_d, W->C_mp, W->C_rat, MPType);
  }

  // clear gamma
  clear_d_mp_rat(W->gamma_d, W->gamma_mp, W->gamma_rat, MPType);

  // clear targetSliceMat & targetSliceVec, if needed
  if (W->targetSliceInit)
  { // clear Mat
    clear_mat(W->targetSliceMat_d, W->targetSliceMat_mp, W->targetSliceMat_rat, MPType);
    // clear Vec
    clear_vec(W->targetSliceVec_d, W->targetSliceVec_mp, W->targetSliceVec_rat, MPType);
  }

  return;
}

void witness_clear_codim(witnessCodim_t *WC, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears WC                                              *
\***************************************************************/
{
  int i, rows = MPType == 1 ? WC->A_mp->rows : WC->A_d->rows;

  // clear W
  for (i = rows - 1; i >= 0; i--)
    free(WC->W[i]);
  free(WC->W);

  // clear A
  clear_mat(WC->A_d, WC->A_mp, WC->A_rat, MPType);

  // clear H
  clear_vec(WC->H_d, WC->H_mp, WC->H_rat, MPType);

  // clear homVarConst
  clear_d_mp_rat(WC->homVarConst_d, WC->homVarConst_mp, WC->homVarConst_rat, MPType);

  // clear B
  clear_mat(WC->B_d, WC->B_mp, WC->B_rat, MPType);

  // clear p
  clear_vec(WC->p_d, WC->p_mp, WC->p_rat, MPType);

  // clear witnessPt_types
  if (WC->witnessPt_types != NULL)
    free(WC->witnessPt_types);

  // clear multiplicities
  if (WC->multiplicities != NULL)
    free(WC->multiplicities);

  // clear deflations_needed
  if (WC->deflations_needed != NULL)
    free(WC->deflations_needed);

  // clear component_nums
  if (WC->component_nums != NULL)
    free(WC->component_nums);

  // clear witnessPts
  if (MPType == 0)
  { // clear witnessPts_d
    if (WC->witnessPts_d != NULL)
    { // clear D memory
      for (i = WC->num_set - 1; i >= 0; i--)
      {
        clear_endpoint_data_d(&WC->witnessPts_d[i]);
      }
      // free memory
      free(WC->witnessPts_d);
    }
  }
  else if (MPType == 1)
  { // clear witnessPts_mp
    if (WC->witnessPts_mp != NULL)
    { // clear MP memory
      for (i = WC->num_set - 1; i >= 0; i--)
      {
        clear_endpoint_data_mp(&WC->witnessPts_mp[i]);
      }
      // free memory
      free(WC->witnessPts_mp);
    }
  }
  else
  { // clear witnessPts_amp
    if (WC->witnessPts_amp != NULL)
    { // clear memory
      for (i = WC->num_set - 1; i >= 0; i--)
      {
        clear_endpoint_data_amp(&WC->witnessPts_amp[i]);
      }
      // free memory
      free(WC->witnessPts_amp);
    }
  }

  return;
}

/////////////// Output ///////////////////////

void printCodimWitnessStructures(FILE *OUT, witnessCodim_t *WC, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the witness structures for the codim to OUT      *
\***************************************************************/
{
  int i, j, rows, cols;

  // print A
  printRawMat(OUT, WC->A_d, WC->A_mp, WC->A_rat, MPType);

  // print W
  if (MPType == 0 || MPType == 2)
  {
    rows = WC->A_d->rows;
    cols = WC->A_d->cols;
  }
  else
  {
    rows = WC->A_mp->rows;
    cols = WC->A_mp->cols;
  }
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      fprintf(OUT, "%d\n", WC->W[i][j]);

  // print H
  printRawVec(OUT, WC->H_d, WC->H_mp, WC->H_rat, MPType);

  // print homVarConst_d
  printRawComp(OUT, WC->homVarConst_d, WC->homVarConst_mp, WC->homVarConst_rat, MPType);

  // print B
  printRawMat(OUT, WC->B_d, WC->B_mp, WC->B_rat, MPType);

  // print p
  printRawVec(OUT, WC->p_d, WC->p_mp, WC->p_rat, MPType);

  return;
}

void printCodimWitnessSet(FILE *OUT, witnessCodim_t *WC, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the witness data for the codim to OUT            *
\***************************************************************/
{
  int i, j, num_set = WC->num_set, base = 10, digits = 0;

  // print the codimension and the number of points in this witness set
  fprintf(OUT, "%d\n%d\n", WC->codim, num_set);

  // loop through the witness points and print the relevant data
  for (i = 0; i < num_set; i++)
  { // print the precision-dependent data
    if (MPType == 0)
    { // print the precision for the endpoint (52)
      fprintf(OUT, "%d\n", 52);
      // print the endpoint
      for (j = 0; j < WC->witnessPts_d[i].endPt->size; j++)
        fprintf(OUT, "%.15e %.15e\n", WC->witnessPts_d[i].endPt->coord[j].r, WC->witnessPts_d[i].endPt->coord[j].i);
      // print the precision for the last approximation (52)
      fprintf(OUT, "%d\n", 52);
      // print the last approximation
      for (j = 0; j < WC->witnessPts_d[i].last_approx->size; j++)
        fprintf(OUT, "%.15e %.15e\n", WC->witnessPts_d[i].last_approx->coord[j].r, WC->witnessPts_d[i].last_approx->coord[j].i);

      // print the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fprintf(OUT, "%.15e\n%d\n%.15e\n%.15e\n", WC->witnessPts_d[i].cond_num, WC->witnessPts_d[i].corank, WC->witnessPts_d[i].smallest_nonzero_SV, WC->witnessPts_d[i].largest_zero_SV);
    }
    else if (MPType == 1)
    { // print the precision for the end point
      fprintf(OUT, "%d\n", (int) mpf_get_prec(WC->witnessPts_mp[i].endPt->coord[0].r));
      // print the end point
      for (j = 0; j < WC->witnessPts_mp[i].endPt->size; j++)
      {
       if (mpfr_number_p(WC->witnessPts_mp[i].endPt->coord[j].r) && mpfr_number_p(WC->witnessPts_mp[i].endPt->coord[j].i))
       {
         mpf_out_str(OUT, base, digits, WC->witnessPts_mp[i].endPt->coord[j].r);
         fprintf(OUT, " ");
         mpf_out_str(OUT, base, digits, WC->witnessPts_mp[i].endPt->coord[j].i);
         fprintf(OUT, "\n");
       }
       else
         fprintf(OUT, "NaN NaN\n");
      }

      // print the precision for the last approximation 
      fprintf(OUT, "%d\n", (int) mpf_get_prec(WC->witnessPts_mp[i].last_approx->coord[0].r)); 
      // print the point 
      for (j = 0; j < WC->witnessPts_mp[i].last_approx->size; j++) 
      { 
       if (mpfr_number_p(WC->witnessPts_mp[i].last_approx->coord[j].r) && mpfr_number_p(WC->witnessPts_mp[i].last_approx->coord[j].i))
       { 
         mpf_out_str(OUT, base, digits, WC->witnessPts_mp[i].last_approx->coord[j].r); 
         fprintf(OUT, " "); 
         mpf_out_str(OUT, base, digits, WC->witnessPts_mp[i].last_approx->coord[j].i); 
         fprintf(OUT, "\n"); 
       } 
       else 
         fprintf(OUT, "NaN NaN\n");
      }

      // print the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fprintf(OUT, "%.15e\n%d\n%.15e\n%.15e\n", WC->witnessPts_mp[i].cond_num, WC->witnessPts_mp[i].corank, WC->witnessPts_mp[i].smallest_nonzero_SV, WC->witnessPts_mp[i].largest_zero_SV);
    }
    else
    { // print the precision for the end point 
      fprintf(OUT, "%d\n", WC->witnessPts_amp[i].curr_prec);
      // print the end point
      if (WC->witnessPts_amp[i].curr_prec < 64)
      { // print using endPt_d
        for (j = 0; j < WC->witnessPts_amp[i].endPt_d->size; j++)
          fprintf(OUT, "%.15e %.15e\n", WC->witnessPts_amp[i].endPt_d->coord[j].r, WC->witnessPts_amp[i].endPt_d->coord[j].i);
      }
      else
      { // print using endPt_mp
        for (j = 0; j < WC->witnessPts_amp[i].endPt_mp->size; j++)
        {
         if (mpfr_number_p(WC->witnessPts_amp[i].endPt_mp->coord[j].r) && mpfr_number_p(WC->witnessPts_amp[i].endPt_mp->coord[j].i))
         {
           mpf_out_str(OUT, base, digits, WC->witnessPts_amp[i].endPt_mp->coord[j].r);
           fprintf(OUT, " ");
           mpf_out_str(OUT, base, digits, WC->witnessPts_amp[i].endPt_mp->coord[j].i);
           fprintf(OUT, "\n");
         }
         else
           fprintf(OUT, "NaN NaN\n");
        }
      }

      // print the precision for the last approximation
      fprintf(OUT, "%d\n", WC->witnessPts_amp[i].last_approx_prec);
      // print the last approximation
      if (WC->witnessPts_amp[i].last_approx_prec < 64)
      { // print using last_approx_d
        for (j = 0; j < WC->witnessPts_amp[i].last_approx_d->size; j++)
          fprintf(OUT, "%.15e %.15e\n", WC->witnessPts_amp[i].last_approx_d->coord[j].r, WC->witnessPts_amp[i].last_approx_d->coord[j].i);
      }
      else
      { // print using last_approx_mp
        for (j = 0; j < WC->witnessPts_amp[i].last_approx_mp->size; j++)
        {
         if (mpfr_number_p(WC->witnessPts_amp[i].last_approx_mp->coord[j].r) && mpfr_number_p(WC->witnessPts_amp[i].last_approx_mp->coord[j].i))
         {
           mpf_out_str(OUT, base, digits, WC->witnessPts_amp[i].last_approx_mp->coord[j].r);
           fprintf(OUT, " ");
           mpf_out_str(OUT, base, digits, WC->witnessPts_amp[i].last_approx_mp->coord[j].i);
           fprintf(OUT, "\n");
         }
         else
           fprintf(OUT, "NaN NaN\n");
        }
      }

      // print the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fprintf(OUT, "%.15e\n%d\n%.15e\n%.15e\n", WC->witnessPts_amp[i].cond_num, WC->witnessPts_amp[i].corank, WC->witnessPts_amp[i].smallest_nonzero_SV, WC->witnessPts_amp[i].largest_zero_SV);
    }

    // print other data - type, multiplicity, component number, deflations needed
    fprintf(OUT, "%d\n%d\n%d\n%d\n", WC->witnessPt_types[i], WC->multiplicities[i], WC->component_nums[i], WC->deflations_needed[i]);
  }

  return;
}

//////////////// Read in from file ///////////////

void setupCodimWitnessStructuresFromFile(FILE *IN, int inputType, witnessCodim_t *WC, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the witness structures for the codim using IN    *
\***************************************************************/
{
  int i, j, rows, cols;

  // setup A
  setupRawMat(IN, WC->A_d, WC->A_mp, &WC->A_rat, MPType, inputType);
  // setup A_rat, A_rows & A_cols
  if (MPType == 0)
  { // setup using A_d
    rows = WC->A_rows = WC->A_d->rows;
    cols = WC->A_cols = WC->A_d->cols;
    init_mat_rat(WC->A_rat, WC->A_rows, WC->A_cols);
    for (i = 0; i < rows; i ++)
      for (j = 0; j < cols; j++)
      {
        mpq_set_d(WC->A_rat[i][j][0], WC->A_d->entry[i][j].r);
        mpq_set_d(WC->A_rat[i][j][1], WC->A_d->entry[i][j].i);
      }
  }
  else if (MPType == 1)
  { // setup using A_mp
    rows = WC->A_rows = WC->A_mp->rows;
    cols = WC->A_cols = WC->A_mp->cols;
    init_mat_rat(WC->A_rat, WC->A_rows, WC->A_cols);
    for (i = 0; i < rows; i ++)
      for (j = 0; j < cols; j++)
      {
        mpf_t_to_rat(WC->A_rat[i][j][0], WC->A_mp->entry[i][j].r);
        mpf_t_to_rat(WC->A_rat[i][j][1], WC->A_mp->entry[i][j].i);
      }
  }
  else
  { // A_rat is already setup so we just need A_rows & A_cols
    WC->A_rows = WC->A_d->rows;
    WC->A_cols = WC->A_d->cols;
  }

  // setup W
  if (MPType == 0 || MPType == 2)
  {
    rows = WC->A_d->rows;
    cols = WC->A_d->cols;
  }
  else
  {
    rows = WC->A_mp->rows;
    cols = WC->A_mp->cols;
  }
  WC->W = (int **)bmalloc(rows * sizeof(int *));
  for (i = 0; i < rows; i++)
  {
    WC->W[i] = (int *)bmalloc(cols * sizeof(int));
    for (j = 0; j < cols; j++)
      fscanf(IN, "%d\n", &WC->W[i][j]);
  }

  // setup H
  setupRawVec(IN, WC->H_d, WC->H_mp, &WC->H_rat, MPType, inputType);

  // setup homVarConst_d
  setupRawComp(IN, WC->homVarConst_d, WC->homVarConst_mp, WC->homVarConst_rat, MPType, inputType, 1, 1);

  // setup B
  setupRawMat(IN, WC->B_d, WC->B_mp, &WC->B_rat, MPType, inputType);

  // setup p
  setupRawVec(IN, WC->p_d, WC->p_mp, &WC->p_rat, MPType, inputType);

  return;
}

void setupCodimWitnessSetFromFile(FILE *IN, witnessCodim_t *WC, int MPType, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the witness data for the codim using IN          *
\***************************************************************/
{
  int i, j, base = 10, pt_prec, approx_prec, count, *components_seen = NULL;
  vec_d tempPt_d, tempApprox_d;
  vec_mp tempPt_mp, tempApprox_mp;

  // initialize
  init_vec_d(tempPt_d, size); init_vec_d(tempApprox_d, size);
  init_vec_mp(tempPt_mp, size); init_vec_mp(tempApprox_mp, size);
  tempPt_d->size = tempPt_mp->size = tempApprox_d->size = tempApprox_mp->size = size;

  // read in the codimension and the number of points for this witness set
  fscanf(IN, "%d\n%d\n", &WC->codim, &WC->num_set);

  // allocate the witness set structures
  WC->witnessPt_types = (int *)bmalloc(WC->num_set * sizeof(int));
  WC->component_nums = (int *)bmalloc(WC->num_set * sizeof(int));
  WC->multiplicities = (int *)bmalloc(WC->num_set * sizeof(int));
  WC->deflations_needed = (int *)bmalloc(WC->num_set * sizeof(int));
  if (MPType == 0)
    WC->witnessPts_d = (endpoint_data_d *)bmalloc(WC->num_set * sizeof(endpoint_data_d));
  else if (MPType == 1)
    WC->witnessPts_mp = (endpoint_data_mp *)bmalloc(WC->num_set * sizeof(endpoint_data_mp));
  else
    WC->witnessPts_amp = (endpoint_data_amp *)bmalloc(WC->num_set * sizeof(endpoint_data_amp));

  // loop through to setup the witness points
  for (i = 0; i < WC->num_set; i++)
  { // read in the precision for the end point
    fscanf(IN, "%d\n", &pt_prec);

    // read in the end point
    if (pt_prec < 64)
    { // read in the point in double precision
      for (j = 0; j < size; j++)
      { // read in the jth coordinate
        fscanf(IN, "%lf %lf\n", &tempPt_d->coord[j].r, &tempPt_d->coord[j].i);
      }
    }
    else
    { // set the precision correctly
      setprec_vec_mp(tempPt_mp, pt_prec);

      // read in the end point in the current precision
      for (j = 0; j < size; j++)
      {
        mpf_inp_str(tempPt_mp->coord[j].r, IN, base);
        mpf_inp_str(tempPt_mp->coord[j].i, IN, base);
        // scan rest of line
        fscanf(IN, "\n");
      }
    }

    // read in the precision for the approximation
    fscanf(IN, "%d\n", &approx_prec);

    // read in the approximation
    if (approx_prec < 64)
    { // read in the point in double precision
      for (j = 0; j < size; j++)
      { // read in the jth coordinate
        fscanf(IN, "%lf %lf\n", &tempApprox_d->coord[j].r, &tempApprox_d->coord[j].i);
      }
    }
    else
    { // set the precision correctly
      setprec_vec_mp(tempApprox_mp, pt_prec);

      // read in the end point in the current precision
      for (j = 0; j < size; j++)
      {
        mpf_inp_str(tempApprox_mp->coord[j].r, IN, base);
        mpf_inp_str(tempApprox_mp->coord[j].i, IN, base);
        // scan rest of line
        fscanf(IN, "\n");
      }
    }

    // setup witnessPts
    if (MPType == 0)
    { // copy to witnessPts_d
      init_endpoint_data_d(&WC->witnessPts_d[i]);

      // setup endPt
      if (pt_prec < 64)
      { // copy _d
        point_cp_d(WC->witnessPts_d[i].endPt, tempPt_d);
      }
      else
      { // copy _mp
        point_mp_to_d(WC->witnessPts_d[i].endPt, tempPt_mp);
      }
      // setup last_approx
      if (approx_prec < 64)
      { // copy _d
        point_cp_d(WC->witnessPts_d[i].last_approx, tempApprox_d);
      }
      else
      { // copy _mp
        point_mp_to_d(WC->witnessPts_d[i].last_approx, tempApprox_mp);
      }
      // read in the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fscanf(IN, "%lf\n%d\n%lf\n%lf\n", &WC->witnessPts_d[i].cond_num, &WC->witnessPts_d[i].corank, &WC->witnessPts_d[i].smallest_nonzero_SV, &WC->witnessPts_d[i].largest_zero_SV);
      // setup other values
      set_zero_d(WC->witnessPts_d[i].finalT);
      WC->witnessPts_d[i].retVal = 0;
    }
    else if (MPType == 1)
    { // copy to witnessPts_mp
      init_endpoint_data_mp(&WC->witnessPts_mp[i]);

      // setup endPt
      if (pt_prec < 64)
      { // copy _d
        point_d_to_mp(WC->witnessPts_mp[i].endPt, tempPt_d);
      }
      else
      { // copy _mp
        point_cp_mp(WC->witnessPts_mp[i].endPt, tempPt_mp);
      }
      // setup last_approx
      if (approx_prec < 64)
      { // copy _d
        point_d_to_mp(WC->witnessPts_mp[i].last_approx, tempApprox_d);
      }
      else
      { // copy _mp
        point_cp_mp(WC->witnessPts_mp[i].last_approx, tempApprox_mp);
      }
      // read in the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fscanf(IN, "%lf\n%d\n%lf\n%lf\n", &WC->witnessPts_mp[i].cond_num, &WC->witnessPts_mp[i].corank, &WC->witnessPts_mp[i].smallest_nonzero_SV, &WC->witnessPts_mp[i].largest_zero_SV);
      // setup other values
      set_zero_mp(WC->witnessPts_mp[i].finalT);
      WC->witnessPts_mp[i].retVal = 0;
    }
    else
    { // setup witnessPts_amp
      init_endpoint_data_amp(&WC->witnessPts_amp[i], pt_prec, approx_prec);
      // setup endPt
      if (pt_prec < 64)
      { // copy _d
        point_cp_d(WC->witnessPts_amp[i].endPt_d, tempPt_d);
      }
      else
      { // copy _mp
        point_cp_mp(WC->witnessPts_amp[i].endPt_mp, tempPt_mp);
      }
      // setup last_approx
      if (approx_prec < 64)
      { // copy _d
        point_cp_d(WC->witnessPts_amp[i].last_approx_d, tempApprox_d);
      }
      else
      { // copy _mp
        point_cp_mp(WC->witnessPts_amp[i].last_approx_mp, tempApprox_mp);
      }
      // read in the condition number, corank, smallest_nonzero_SV, largest_zero_SV
      fscanf(IN, "%lf\n%d\n%lf\n%lf\n", &WC->witnessPts_amp[i].cond_num, &WC->witnessPts_amp[i].corank, &WC->witnessPts_amp[i].smallest_nonzero_SV, &WC->witnessPts_amp[i].largest_zero_SV);
      // setup other values
      if (pt_prec < 64)
      {
        set_zero_d(WC->witnessPts_amp[i].finalT_d);
      }
      else
      {
        set_zero_mp(WC->witnessPts_amp[i].finalT_mp);
      }
      WC->witnessPts_amp[i].retVal = 0;
    }

    // read in other data - type, multiplicity, component number, deflations needed
    fscanf(IN, "%d\n%d\n%d\n%d\n", &WC->witnessPt_types[i], &WC->multiplicities[i], &WC->component_nums[i], &WC->deflations_needed[i]);
  }

  // setup the counts
  WC->num_nonsing = WC->num_sing = WC->num_components = count = 0;
  for (i = 0; i < WC->num_set; i++)
  { // see if sing or non-singular
    if (WC->witnessPt_types[i] == NON_SINGULAR)
      WC->num_nonsing++;
    else
      WC->num_sing++;

    // see if this is a valid component
    if (WC->component_nums[i] >= 0)
    { // see if we have seen the component number
      pt_prec = 0;
      for (j = 0; j < count && !pt_prec; j++)
        if (components_seen[j] == WC->component_nums[i])
          pt_prec = 1;

      if (!pt_prec)
      { // we have not seen this component number
        components_seen = (int *)brealloc(components_seen, (count + 1) * sizeof(int));
        components_seen[count] = WC->component_nums[i];
        count++;
        WC->num_components++;
      }
    }
  }

  clear_vec_d(tempPt_d); clear_vec_d(tempApprox_d);
  clear_vec_mp(tempPt_mp); clear_vec_mp(tempApprox_mp);
  free(components_seen);

  return;
}

void setupWitnessDataFromFile(char *witnessName, char *newWitnessName, char *preprocFile, char *degreeFile, witness_t *W, tracker_config_t *T, int needToSetupProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup W using the file described by witnessName        *
\***************************************************************/
{
  int codim_index, count, inputType;
  FILE *IN = NULL;

  // setup curr_precision
  W->curr_precision = T->Precision;

  // setup W->PPD
  setupPreProcData(preprocFile, &W->PPD);

  // verify that we are using only 1 homogenous variable group
  if (W->PPD.num_hom_var_gp + W->PPD.num_var_gp != 1)
  { // exit immediately
    printf("ERROR: Positive dimensional setup is implemented for systems with exactly one variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (needToSetupProg)
  { // setup the slp
    W->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    W->orig_variables = W->new_variables = setupProg(W->Prog, T->Precision, T->MPType);
  }
  else
  { // do no setup slp
    W->Prog = NULL;
    // the following assumes 1 variable group
    W->orig_variables = W->new_variables = W->PPD.type[0] + W->PPD.size[0];
  }

  // setup the number of functions
  W->num_funcs = W->PPD.num_funcs;

  if (needToSetupProg)
  { // compute the rank
    if (T->MPType == 0 || T->MPType == 2)
      W->system_rank = rank_finder_d(&W->PPD, W->Prog, T, W->orig_variables);
    else
      W->system_rank = rank_finder_mp(&W->PPD, W->Prog, T, W->orig_variables);

    // setup orig_degrees, new_degrees & P
    setupDegrees_orig_new_perm(&W->orig_degrees, &W->new_degrees, &W->P, W->num_funcs, W->PPD.num_var_gp + W->PPD.num_hom_var_gp, degreeFile);
  }
  else
  { // assume each function is independent
    W->system_rank = MAX(W->orig_variables, W->PPD.num_funcs);

    W->orig_degrees = W->new_degrees = W->P = NULL;
  }

  // now that all of the general data is setup, check for existence of the witness data
  IN = fopen(witnessName, "r");
  if (IN == NULL)
  { // file does not exist
    printf("\n\nERROR: '%s' does not exist!!!\n\n", witnessName);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  // now that we have its existence, move it!
  fclose(IN);
  rename(witnessName, newWitnessName);
  IN = fopen(newWitnessName, "r");

  // read in the number of variables
  fscanf(IN, "%d\n", &count);

  // make sure that the number of variables is correct
  if (count != W->orig_variables)
  { // put witness_data_old back at witnessName and exit!
    fclose(IN);
    rename(newWitnessName, witnessName);
    printf("ERROR: The number of variables described in '%s' is not correct!\n", witnessName);
    bexit(ERROR_INVALID_SIZE);
  }

  // find the number of codimensions
  fscanf(IN, "%d\n", &W->num_codim);

  // allocate the codim
  W->codim = (witnessCodim_t *)bmalloc(W->num_codim * sizeof(witnessCodim_t));

  // loop through to setup the codim
  for (codim_index = 0; codim_index < W->num_codim; codim_index++)
  {
    setupCodimWitnessSetFromFile(IN, &W->codim[codim_index], T->MPType, W->orig_variables);
  }

  // read in '-1' at the end of the witness set and the MPType for the codim structures
  fscanf(IN, "%d\n\n%d\n", &count, &inputType);

  // loop through to setup the structures for each of the codimensions
  for (codim_index = 0; codim_index < W->num_codim; codim_index++)
  {
    setupCodimWitnessStructuresFromFile(IN, inputType, &W->codim[codim_index], T->MPType);
  }

  // setup gamma
  if (T->MPType == 0)
  { // only setup gamma_d
    get_comp_rand_d(W->gamma_d);
  }
  else if (T->MPType == 1)
  { // only setup gamma_mp
    init_mp(W->gamma_mp);
    get_comp_rand_mp(W->gamma_mp);
  }
  else
  { // setup gamma_rat, gamma_mp & gamma_d
    get_comp_rand_rat(W->gamma_d, W->gamma_mp, W->gamma_rat, W->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // target slice info is not setup
  W->targetSliceInit = 0;

  fclose(IN);

  return;
}

void setupDegrees_orig_new_perm(int **orig_degrees, int **new_degrees, int **perm, int num_funcs, int num_var_gps, char *degreeFile)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads in the degrees from degreeFile and sets up arrays*
\***************************************************************/
{
  int i, j, max, max_loc, *tempInts = NULL;
  FILE *tempFile = fopen(degreeFile, "r");
  if (tempFile == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // make sure num_var_gps == 1
  if (num_var_gps > 1)
  { // exit immediately
    printf("ERROR: Bertini can only setup the degrees when there is only 1 variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup original degrees
  *orig_degrees = (int *)bmalloc(num_funcs * sizeof(int));
  for (i = 0; i < num_funcs; i++)
    fscanf(tempFile, "%d\n\n", &(*orig_degrees)[i]);

  // close the file containing the degrees
  fclose(tempFile);

  // allocate perm
  *perm = (int *)bmalloc(num_funcs * sizeof(int));

  // order the functions so that the ith largest degree is perm[i]
  tempInts = (int *)bmalloc(num_funcs * sizeof(int));
  for (i = 0; i < num_funcs; i++)
    tempInts[i] = 1; // will be = 0 when the function has been used

  for (i = 0; i < num_funcs; i++)
  { // find the largest degreee still available
    max = max_loc = -1;
    for (j = 0; j < num_funcs; j++)
    { // see if this function is still available and larger than the current max
      if (tempInts[j] == 1 && max < (*orig_degrees)[j])
      { // update max & max_loc
        max = (*orig_degrees)[j];
        max_loc = j;
      }
    }
    (*perm)[i] = max_loc;
    tempInts[max_loc] = 0;
  }

  // clear tempInts
  free(tempInts);

  // setup new_degrees
  *new_degrees = (int *)bmalloc(num_funcs * sizeof(int));
  for (i = 0; i < num_funcs; i++)
    (*new_degrees)[i] = (*orig_degrees)[(*perm)[i]];

  return;
}

