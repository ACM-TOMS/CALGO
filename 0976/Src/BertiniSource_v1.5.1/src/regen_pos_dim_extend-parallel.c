// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"
#include "parallel.h"
#include "ppParse.h"

void regenExtendMenu(int *codim_index_selected, int *component_selected, witness_t *W);
void regenExtend(witness_t *W, tracker_config_t *T, int codim_index, int component_num);
int  regenExtendEvaluation(prog_t *SLP, witness_t *W, int codim_index, int point_num, tracker_config_t *T);
void regenExtendVerification(int numComponents, int *codim_index, int *component_number, witness_t *witnessSets, tracker_config_t *T, prog_t *SLP, int pathMod);
void regenExtendSetup_points(point_d extra_d, point_mp extra_mp, mat_d Points_d, mat_mp Points_mp, int numComponents, int *codim_index, int *component_number, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP);
int  regenExtendSetup_points_finish(int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, regen_pos_dim_t *RPD);
void regenExtendSetup_system(int *maxCodim, int *specificCodim, double intrinsicCutoffMultiplier, int system_rank, point_d extra_d, point_mp extra_mp, int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP, preproc_data *PPD, regen_pos_dim_t *RPD);
void regenExtendSetup_system_coeff(point_d extra_d, point_mp extra_mp, int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, regen_pos_dim_t *RPD);
int  regenExtendSetup_point(int totalPts, mat_d Point_d, mat_mp Points_mp, int numComponents, int *codim_index, int *point_index, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP, FILE *CODIM, FILE *PTS, point_d extra_d, point_mp extra_mp);
void setupDegrees_extend(int **orig_degrees, int **new_degrees, int **perm, int currCodim, int top_funcs, int bottom_funcs, int num_var_gps, char *degreeFile);
void regenExtendSetup_witness(char *origFile, char *newFile, tracker_config_t *T, regen_pos_dim_t *RPD);
int  regenExtendSetup_start(char *origFile, char *newFile, int curr_codim_index, tracker_config_t *T, regen_pos_dim_t *RPD);
void regenExtend_run(double midpoint_tol, int maxCodim, int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, int curr_codim_index, int my_id, int num_processes, int headnode);
void head_regenExtend_run(double midpoint_tol, int maxCodim, int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, int curr_codim_index, int my_id, int num_processes, int headnode);
void regenExtend_reclassify(tracker_config_t *T, regen_pos_dim_t *RPD, int maxCodim, int specificCodim);
void regenExtend_reclassify_codim(tracker_config_t *T, regen_pos_dim_t *RPD, int codim_index);

void regenExtendMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: extend a witness set using regenerative cascade        *
\***************************************************************/
{
  int i, j, userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0, numComponents = 0;
  int rV, size_of_string = 255, quit = 0, selection_found = 0, system_rank = 0, curr_codim_index, numRandom;
  int *codim_index = NULL, *component_number = NULL;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  char ch, *tempStr = NULL, *inputFile = NULL;
  tracker_config_t T;
  witness_t witnessSet, *subWitnessSets = NULL;
  FILE *IN = NULL, *OUT = NULL;
  preProcessArray ppArray;
  variablegroupArray vargpArray;
  prog_t Prog;
  preproc_data PPD;
  point_d extra_d;
  point_mp extra_mp;
  mat_d Pts_d;
  mat_mp Pts_mp;
  regen_pos_dim_t RPD;

  init_mat_d(Pts_d, 0, 0);
  init_mat_mp(Pts_mp, 0, 0);

  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType);

  // since using regenerative extension -- only generically reduced
  reducedOnly = 1;

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

#ifdef _OPENMP
  #pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }

  // read in the number of nontrivial components
  printf("\n\n*************** Regeneration Extension ****************\n\n");
  do
  { // initialize to some value
    numComponents = -2;

    // find the number
    printf("Please enter the number of nontrivial components (-1 to quit): ");
    rV = scanf("%d", &numComponents);

    if (rV < 0)
    { // at EOF - so we quit
      numComponents = -1;
    }
    else
    { // not at EOF - flush the buffer
      do 
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");  
      }
      else if (numComponents < -1)
      { // verify >= -1
        printf("\nThe number of components must be nonnegative (or -1 to quit)!\n");
      }    
    }
  } while (numComponents < -1);

  if (numComponents < 0)
  { // clear MP and exit
#ifdef _HAVE_MPI
    worker_info sendType;
    sendType.dataType = STOPCODE;
    bcast_worker_info(&sendType, my_id, headnode);
#endif

    // clear MP
    clearMP();
      
    return;
  }
  else if (numComponents == 0)
  { // use standard regen cascade approach
    printf("\n\nSince no components are known, Bertini will use the Regenerative Cascade algorithm.\n");
    
    // call NID using regen cascade
    numericalIrredDecomp(currentSeed, MPType, 2, my_id, num_processes, headnode);

    return;
  }

  // We indeed have some components! Print message about components and what is computed
  printf("\nNOTE: Regeneration extension is only implemented for generically reduced components (both input and output)!\n");

  if (numComponents > 1)
  { // more than one component -- message about independence of witness sets 
    printf("\nNOTE: Regeneration extension assumes the witness sets for the %d components are independent!\n", numComponents);
  }

  // Load the witness sets for the subsystems
  subWitnessSets = (witness_t *)bmalloc(numComponents * sizeof(witness_t));
  codim_index = (int *)bmalloc(numComponents * sizeof(int));
  component_number = (int *)bmalloc(numComponents * sizeof(int));

  // setup the components
  for (j = 0; j < numComponents; j++)
  { // read in the name of the old input file
    tempStr = (char *)bmalloc(((int) log10(size_of_string) + 10) * sizeof(char));
    snprintf(tempStr, log10(size_of_string) + 10, "%%%ds", size_of_string);
    inputFile = (char *)bmalloc((size_of_string + 1) * sizeof(char));
    for (i = 0; i <= size_of_string; i++)
      inputFile[i] = '\0';
    do
    { // find the name of the file
      printf("\nSetup for component %d of %d.\n", j+1, numComponents);
      printf("Please enter the name of the corresponding input file or type\n'quit' or 'exit' (max of %d characters): ", size_of_string);
      rV = scanf(tempStr, inputFile);

      if (rV <= 0)
      { // at EOF so we need to quit
        quit = 1;
        selection_found = 1;
      }
      else
      { // flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (!strcmp(inputFile, "'quit'") || !strcmp(inputFile, "quit") || !strcmp(inputFile, "'exit'") || !strcmp(inputFile, "exit"))
        { // user wants to quit
          quit = 1;
          selection_found = 1;
        }
        else
        { // check for existence of this file
          IN = fopen(inputFile, "r");
          if (IN == NULL)
          { // not exist
            printf("\nA file named \"%s\" does not exist!\n\n", inputFile);
            selection_found = 0;
          }
          else
          { // exists!
            fclose(IN);
            selection_found = 1;
          }
        }
      }
    } while (!selection_found);
    printf("\n");

    if (quit)
    { // clear MP and exit
#ifdef _HAVE_MPI
      worker_info sendType;
      sendType.dataType = STOPCODE;
      bcast_worker_info(&sendType, my_id, headnode);
#endif

      // free memory
      free(tempStr);
      free(inputFile);

      for (i = 0; i < j - 1; i++)
        witness_clear(&subWitnessSets[i], T.MPType);
      free(subWitnessSets);
      free(codim_index);
      free(component_number);

      // clear MP
      clearMP();

      return;
    }

    // we need to do a basic parsing of the old input file
    IN = fopen(inputFile, "r");
    splitParse(IN, "func_input_rand", "config_old");
    fclose(IN);

    // check for 'random' and 'random_real' -- not permited
    IN = fopen("func_input_rand", "r");
    OUT = fopen("func_input_old", "w");
    numRandom = setupRandomValues(OUT, IN, 1, T.AMP_max_prec);
    fclose(IN);
    fclose(OUT);
    remove("func_input_rand");

    // check number of random
    if (numRandom > 0)
    { // random not allowed in regeneration extension!
      printf("\nERROR: Regeneration extension assumes that 'random' and 'random_real' are not used!\n       Please replace with 'constant' and define with actual values.\n");
      bexit(ERROR_CONFIGURATION);
    }

    // preprocess func_input_old
    IN = fopen("func_input_old", "r");
    preProcessParse(&ppArray, IN, paramHom, 1); 
    fclose(IN);

    // add the subfunction number information
    setup_subFuncData_numbers(&ppArray);

    // seutp variable group information
    setup_variablegroupArray(&ppArray, &vargpArray, userHom == 1);
 
    if (userHom == 1)
    { // verify information for user defined homotopy
      verify_userdefined_homotopy(&ppArray, userHom);
    }
    else
    { // verify information for standard tracking
      verify_standard_homotopy(&ppArray, &vargpArray, 1);
    }

    // error checking
    if (vargpArray.numHomGps + vargpArray.numVarGps > 1)
    { // we assume only 1 variable group!
      printf("ERROR: Bertini excepts only one variable group!\n");
      remove("func_input_old");
      remove("config_old");
      bexit(ERROR_INPUT_SYSTEM);
    }
    if (vargpArray.numHomGps > 0)
    { // we assume we are in affine coordinates!
      printf("ERROR: Bertini does not expect a homogeneous variable group!\n");
      remove("func_input_old");
      remove("config_old");
      bexit(ERROR_INPUT_SYSTEM);
    }

    // print preproc_data
    IN = fopen("preproc_data_old", "w");
    print_preproc_data(IN, ppArray.types[FUNCTIONTYPE].numType, vargpArray.numVarGps, vargpArray.numHomGps, vargpArray.types, vargpArray.sizes);
    fclose(IN);

    // clear parsing information
    clear_variablegroupArray(&vargpArray);
    clear_preProcessArray(&ppArray);
 
    // read in the name of the witness data file
    for (i = 0; i <= size_of_string; i++)
      inputFile[i] = '\0';
    do
    { // find the name of the file
      printf("Please enter the name of the corresponding witness_data file or type\n'quit' or 'exit' (max of %d characters): ", size_of_string);
      rV = scanf(tempStr, inputFile);

      if (rV <= 0)
      { // at EOF so we need to quit
        quit = 1;
        selection_found = 1;
      }
      else
      { // flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (!strcmp(inputFile, "'quit'") || !strcmp(inputFile, "quit") || !strcmp(inputFile, "'exit'") || !strcmp(inputFile, "exit"))
        { // user wants to quit
          quit = 1;
          selection_found = 1;
        }
        else
        { // check for existence of this file
          IN = fopen(inputFile, "r");
          if (IN == NULL)
          { // not exist
            printf("\nA file named \"%s\" does not exist!\n\n", inputFile);
            selection_found = 0;
          }
          else
          { // exists!
            fclose(IN);
            selection_found = 1;
          }
        }
      }
    } while (!selection_found);
    printf("\n");

    if (quit)
    { // clear MP and exit
#ifdef _HAVE_MPI
      worker_info sendType;
      sendType.dataType = STOPCODE;
      bcast_worker_info(&sendType, my_id, headnode);
#endif

      // free memory
      free(tempStr);
      free(inputFile);

      for (i = 0; i < j - 1; i++)
        witness_clear(&subWitnessSets[i], T.MPType);
      free(subWitnessSets);
      free(codim_index);
      free(component_number);

      // clear MP
      clearMP();

      return;
    }

    // setup subWitnessSet[j] -- do not setup SLP or degrees
    setupWitnessDataFromFile(inputFile, inputFile, "preproc_data_old", "", &subWitnessSets[j], &T, 0);

    // remove old data
    remove("preproc_data_old");
    remove("func_input_old");
    remove("config_old");

    // select dimension and component number for the irreducible component
    regenExtendMenu(&codim_index[j], &component_number[j], &subWitnessSets[j]);

    if (codim_index[j] == -1 || component_number[j] == -1)
    { // clear MP and exit
#ifdef _HAVE_MPI
      worker_info sendType;
      sendType.dataType = STOPCODE;
      bcast_worker_info(&sendType, my_id, headnode);
#endif

      // free memory
      free(tempStr);
      free(inputFile);

      for (i = 0; i < j - 1; i++)
        witness_clear(&subWitnessSets[i], T.MPType);
      free(subWitnessSets);
      free(codim_index);
      free(component_number);

      // clear MP
      clearMP();

      return;
    }
  }

  // setup SLP for whole system 
  setupProg(&Prog, T.Precision, T.MPType);  

  // setup the rest of T, if needed
  T.numVars = Prog.numVars;
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

  // setup PPD
  setupPreProcData("preproc_data", &PPD);
 
  // verify the setup
  regenExtendVerification(numComponents, codim_index, component_number, subWitnessSets, &T, &Prog, pathMod);

  // compute the system rank
  if (T.MPType == 0 || T.MPType == 2)
    system_rank = rank_finder_d(&PPD, &Prog, &T, Prog.numVars);
  else
    system_rank = rank_finder_mp(&PPD, &Prog, &T, Prog.numVars);

  // setup starting points & system
  init_point_d(extra_d, 0);
  init_point_mp(extra_mp, 0);
  regenExtendSetup_points(extra_d, extra_mp, Pts_d, Pts_mp, numComponents, codim_index, component_number, subWitnessSets, &T, &Prog); 
  regenExtendSetup_system(&maxCodim, &specificCodim, intrinsicCutoffMultiplier, system_rank, extra_d, extra_mp, numComponents, codim_index, subWitnessSets, &T, &Prog, &PPD, &RPD);

  // Use the system to finish setting up the points
  curr_codim_index = regenExtendSetup_points_finish(numComponents, codim_index, subWitnessSets, &T, &RPD);

  // perform through the extension
  regenExtend_run(midpoint_tol, maxCodim, pathMod, &RPD, &T, curr_codim_index, my_id, num_processes, headnode);

  // print initial output chart
  regen_pos_dim_OutputChart(&RPD, stdout, maxCodim);

  // recompute the data without using the special randomization and remove singular points
  regenExtend_reclassify(&T, &RPD, maxCodim, specificCodim);

  // initially setup witnessSet & clear RPD
  regen_pos_dim_copyWitness_clear(&witnessSet, &RPD, "witness_superset", T.MPType, T.AMP_max_prec, specificCodim);

  // remove all of the extra multiple points and setup multplicity to be 1 for nonsingular --> special randomization can cause issues 
  for (i = 0; i < witnessSet.num_codim; i++)
  {
    multiplicity_witness(&witnessSet, i, T.MPType, T.final_tol_times_mult);

    // set the multiplicity on all nonsingular points to be 1
    for (j = 0; j < witnessSet.codim[i].num_set; j++)
      if (witnessSet.codim[i].witnessPt_types[j] == NON_SINGULAR)
        witnessSet.codim[i].multiplicities[j] = 1;
  }

  // setup to do everything extrinsic (use original variables)
  setupWitnessTotallyExtrinisic(&witnessSet, T.MPType, T.AMP_max_prec);

  // print a witness superset file
  witnessSupersetOutput(&witnessSet, T.MPType);

  // if we are using MPI - we need to tell the workers if they are doing decompositions
  if (num_processes > 1)
  { // broadcast numerical irreducible decomposition
#ifdef _HAVE_MPI
    worker_info sendType;
    if (supersetOnly)
      sendType.dataType = STOPCODE;
    else
      sendType.dataType = IRREDDECOMP;
    bcast_worker_info(&sendType, my_id, headnode);
#endif
  }

  // determine if we should do the decomposition
  if (!supersetOnly)
  { 
    int **fullRankProgInfo = NULL; // 0 - just needs NULLed out, 1 - deflated properly and needs to be cleared, -1 - did not deflate properly but needs to be cleared
    prog_t ***fullRankProgs = NULL;
    endpoint_data_d **endPts_d = NULL;
    endpoint_data_mp **endPts_mp = NULL;
    endpoint_data_amp **endPts_amp = NULL;
    membership_slice_moving_t *sliceMover = NULL;

    // setup info and remove any singular points
    junkRemoval(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, &T, pathMod, midpoint_tol, reducedOnly, specificCodim, 0, my_id, num_processes, headnode); // topDimension is unknown

    // do the break-up into irreducible components
    pureDecomp(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, &witnessSet, &T, pathMod, my_id, num_processes, headnode);

    // display decomposition chart
    numIrredDecompChart(&witnessSet, stdout, T.MPType, reducedOnly);

    // create output files
    numIrredDecompOutput(&witnessSet, &T, 7, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 7, genType == 2

    // display deflation errors
    displayDeflationSummary(fullRankProgInfo, &witnessSet);

    // clear the deflation information
    clear_sliceMover_fullRankProgs(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, T.MPType);
  }

  // clear memory
  free(tempStr);
  free(inputFile);
  clearProg(&Prog, T.MPType, 0);
  preproc_data_clear(&PPD);

  for (i = 0; i < numComponents; i++)
    witness_clear(&subWitnessSets[i], T.MPType);
  free(subWitnessSets);
  free(codim_index);
  free(component_number);

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear Pts
  clear_point_d(extra_d);
  clear_point_mp(extra_mp);
  clear_mat_d(Pts_d);
  clear_mat_mp(Pts_mp);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void regenExtend_run(double midpoint_tol, int maxCodim, int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, int curr_codim_index, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: run the actual extension now that everything is setup  *
\***************************************************************/
{
  if (num_processes > 1)
  { // do parallel tracking using MPI - tell the workers what they will be doing
#ifdef _HAVE_MPI
    worker_info sendType;
    sendType.dataType = REGEN_EXTEND;
    bcast_worker_info(&sendType, my_id, headnode);

    // do parallel preparation
    head_regenExtend_run(midpoint_tol, maxCodim, pathMod, RPD, T, curr_codim_index, my_id, num_processes, headnode);
#endif
  }
  else
  { // prepare the first codimension
    if (curr_codim_index + 1 < maxCodim)
    { // need to run
      int size, num_paths, num_crossings;
      char *str = NULL;
      char witName[] = "witness_superset", startName[] = "startRPD", outName[] = "output_extend", midName[] = "midpath_data", failName[] = "fail_extend", rawPrepareFile[] = "rawout_prepare";
      FILE *MIDOUT, *START, *RAWOUT, *FAIL, *OUT;
      trackingStats trackCount;
      init_trackingStats(&trackCount);
  
      // open OUT & FAIL
      OUT = fopen(outName, "w");
      FAIL = fopen(failName, "w");
      MIDOUT = fopen(midName, "w");
      START = fopen("startRPD_test", "r");
    
      // open RAWOUT for preparing
      size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, RPD->codim[curr_codim_index + 1].codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawPrepareFile, RPD->codim[curr_codim_index + 1].codim);
      RAWOUT = fopen(str, "w");

      // setup the name of the file to contain the next set of start points
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[curr_codim_index + 1].codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, RPD->codim[curr_codim_index + 1].codim);

      // create the next set of start points for the next codim
      num_paths = regen_pos_dim_PrepareNextCodim(pathMod, RPD, T, OUT, RAWOUT, MIDOUT, FAIL, START, str, curr_codim_index, maxCodim, &trackCount, regen_pos_dim_moving_linear_eval_d, regen_pos_dim_moving_linear_eval_mp, change_regen_pos_dim_prec);

      // close files 
      fclose(RAWOUT);
      fclose(MIDOUT);
      fclose(START);
      fclose(OUT);
      fclose(FAIL);

      // check for path crossings
      num_crossings = 0;

      // check to see if using intrinsic slice
      if (RPD->codim[curr_codim_index].useIntrinsicSlice)
        midpoint_checker(num_paths, RPD->codim[curr_codim_index].codim + 1, midpoint_tol, &num_crossings);
      else
        midpoint_checker(num_paths, RPD->new_variables, midpoint_tol, &num_crossings);

      // print message to screen about path crossing
      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
  
      free(str);

      // now that the next codimension is prepared, track through the rest of them!
      regen_pos_dim_seq_track(curr_codim_index + 1, maxCodim, &trackCount, pathMod, midpoint_tol, T, RPD, startName, witName);
    }
  }

  return;
}

void regenExtendMenu(int *codim_index_selected, int *component_selected, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: menu for selecting dim/components to extend            *
\***************************************************************/
{
  int i, j, codim_index, min_deg, max_deg, count, rV, dim_number, component_number, selection_made;
  int *degrees = NULL, *dim = (int *)bmalloc(W->num_codim * sizeof(int)), *codim_good = (int *)bmalloc(W->num_codim * sizeof(int));
  int *genReduced = NULL;
  char ch;

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
    printf("\nThere are no classified components to regenerate!\n\n");
    *codim_index_selected = *component_selected = -1;
    free(dim);
    free(codim_good);
    return;
  }

  // so we have atleast one classified component
  do
  { // initialize
    selection_made = 0;

    // print title
    printf("\n\n*************** Components to Regenerate ****************\n\n");

    // display a catalog of the available components in each codim
    for (codim_index = 0; codim_index < W->num_codim; codim_index++)
    { // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      genReduced = (int *)brealloc(genReduced, W->codim[codim_index].num_components * sizeof(int));

      for (i = 0; i < W->codim[codim_index].num_components; i++)
      {      
        degrees[i] = 0;
        genReduced[i] = 1;
      }

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
        {
          degrees[W->codim[codim_index].component_nums[i]]++;
          if (genReduced[W->codim[codim_index].component_nums[i]] && (W->codim[codim_index].multiplicities[i] != 1 || W->codim[codim_index].deflations_needed[i] != 0))
            genReduced[W->codim[codim_index].component_nums[i]] = 0;
        }
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

    printf("\nPlease select a dimension to regenerate (-1 to quit): ");
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
        codim_index = *codim_index_selected = i;
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
        printf("   component %d has degree %d (gen. reduced: %s)\n", i, degrees[i], genReduced[i] ? "Yes" : "No");
      printf("\n");

      printf("\nPlease select a component to regenerate (-1 to quit, -2 to regenerate all): ");
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
        else if (component_number < -2 || component_number >= W->codim[codim_index].num_components)
        { // component_number is not valid
          printf("\nThe component %d is not valid!\n", component_number);
          selection_made = 0;
        }
        else
        { // component_number is valid
          selection_made = 1;
          *component_selected = component_number;
        }
      }
    } while (!selection_made);
  }
  else
  {
    *codim_index_selected = *component_selected = -1;
  }

  free(degrees);
  free(genReduced);
  free(dim);
  free(codim_good);

  return;
}

int regenExtendSetup_point(int totalPoints, mat_d Pts_d, mat_mp Pts_mp, int numComponents, int *codim_index, int *point_index, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP, FILE *CODIM, FILE *PTS, point_d extra_d, point_mp extra_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - printed to PTS, 0 - printed to CODIM       *
* NOTES: setup the specific point for extending regeneration    *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, j, rV = 0, currVar = 0;

  if (T->MPType == 0)
  { // use _d
    eval_struct_d e;
    point_d pt, approx;
    comp_d time, homValue, homValue_approx;

    // initialize
    init_eval_struct_d(e, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_d(pt, SLP->numVars); init_point_d(approx, SLP->numVars);
    pt->size = approx->size = SLP->numVars;
    init_d(time); init_d(homValue); init_d(homValue_approx);
    set_zero_d(time);

    // increase Pts_d
    increase_size_mat_d(Pts_d, totalPoints, SLP->numVars);
    Pts_d->rows = totalPoints;
    Pts_d->cols = SLP->numVars;

    // pick a homogenizing coordinate at random
    get_comp_rand_d(&pt->coord[0]);
    set_d(&approx->coord[0], &pt->coord[0]);

    // setup rest of pt
    currVar = 1;
    for (i = 0; i < numComponents; i++)
    { // compute normalizing value and then setup coordinates
      div_d(homValue, &pt->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_d[point_index[i]].endPt->coord[0]);
      div_d(homValue_approx, &approx->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_d[point_index[i]].last_approx->coord[0]);
      for (j = 1; j < subWitnessSets[i].orig_variables; j++)
      {
        mul_d(&pt->coord[currVar], homValue, &subWitnessSets[i].codim[codim_index[i]].witnessPts_d[point_index[i]].endPt->coord[j]);
        mul_d(&approx->coord[currVar], homValue_approx, &subWitnessSets[i].codim[codim_index[i]].witnessPts_d[point_index[i]].last_approx->coord[j]);
        currVar++;
      }
    }

    // setup extra coordinates
    for (i = currVar; i < SLP->numVars; i++)
    {
      mul_d(&pt->coord[i], &extra_d->coord[i - currVar], &pt->coord[0]);
      mul_d(&approx->coord[i], &extra_d->coord[i - currVar], &approx->coord[0]);
    }

    // normalize the point & copy to Pts_d
    homValue->r = infNormVec_d(pt);
    homValue_approx->r = infNormVec_d(approx);
    for (i = 0; i < SLP->numVars; i++)
    {
      mul_rdouble_d(&pt->coord[i], &pt->coord[i], homValue->r);
      set_d(&Pts_d->entry[totalPoints - 1][i], &pt->coord[i]);
      mul_rdouble_d(&approx->coord[i], &approx->coord[i], homValue_approx->r);
    }

    // evaluate
    evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, pt, time, SLP);

    // see if vanishes sufficiently
    if (infNormVec_d(e.funcVals) < T->final_tol_times_mult)
    { // print to CODIM
      rV = 0;

      fprintf(CODIM, "52\n");
      for (i = 0; i < SLP->numVars; i++)
        fprintf(CODIM, "%.15e %.15e\n", pt->coord[i].r, pt->coord[i].i);
      fprintf(CODIM, "\n52\n");
      for (i = 0; i < SLP->numVars; i++)
        fprintf(CODIM, "%.15e %.15e\n", approx->coord[i].r, approx->coord[i].i);
      fprintf(CODIM, "\n");
    }
    else
    { // print to PTS
      rV = 1;
      
      for (i = 0; i < SLP->numVars; i++)
        fprintf(PTS, "%.15e %.15e\n", pt->coord[i].r, pt->coord[i].i);
      fprintf(PTS, "\n");
    }

    // clear
    clear_eval_struct_d(e);
    clear_point_d(pt); clear_point_d(approx);
    clear_d(time); clear_d(homValue); clear_d(homValue_approx);
  }
  else if (T->MPType == 1)
  { // use _mp
    eval_struct_mp e;
    point_mp pt, approx;
    comp_mp time, homValue, homValue_approx;

    // initialize
    init_eval_struct_mp(e, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_mp(pt, SLP->numVars); init_point_mp(approx, SLP->numVars);
    pt->size = approx->size = SLP->numVars;
    init_mp(time); init_mp(homValue); init_mp(homValue_approx);
    set_zero_mp(time);

    // increase Pts_mp_
    increase_size_mat_mp(Pts_mp, totalPoints, SLP->numVars);
    Pts_mp->rows = totalPoints;
    Pts_mp->cols = SLP->numVars;

    // pick a homogenizing coordinate at random
    get_comp_rand_mp(&pt->coord[0]);
    set_mp(&approx->coord[0], &pt->coord[0]);

    // setup rest of pt
    currVar = 1;
    for (i = 0; i < numComponents; i++)
    { // compute normalizing value and then setup coordinates
      div_mp(homValue, &pt->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_mp[point_index[i]].endPt->coord[0]);
      div_mp(homValue_approx, &approx->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_mp[point_index[i]].last_approx->coord[0]);
      for (j = 1; j < subWitnessSets[i].orig_variables; j++)
      {
        mul_mp(&pt->coord[currVar], homValue, &subWitnessSets[i].codim[codim_index[i]].witnessPts_mp[point_index[i]].endPt->coord[j]);
        mul_mp(&approx->coord[currVar], homValue_approx, &subWitnessSets[i].codim[codim_index[i]].witnessPts_mp[point_index[i]].last_approx->coord[j]);
        currVar++;
      }
    }

    // setup extra coordinates
    for (i = currVar; i < SLP->numVars; i++)
    {
      mul_mp(&pt->coord[i], &extra_mp->coord[i - currVar], &pt->coord[0]);
      mul_mp(&approx->coord[i], &extra_mp->coord[i - currVar], &approx->coord[0]);
    }

    // normalize the point & copy to Pts_mp
    infNormVec_mp2(homValue->r, pt);
    infNormVec_mp2(homValue_approx->r, approx);
    for (i = 0; i < SLP->numVars; i++)
    {
      mul_rmpf_mp(&pt->coord[i], &pt->coord[i], homValue->r);
      set_mp(&Pts_mp->entry[totalPoints - 1][i], &pt->coord[i]);
      mul_rmpf_mp(&approx->coord[i], &approx->coord[i], homValue_approx->r);
    }

    // evaluate
    evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, pt, time, SLP);

    // see if vanishes sufficiently
    if (infNormVec_mp(e.funcVals) < T->final_tol_times_mult)
    { // print to CODIM
      rV = 0;

      fprintf(CODIM, "%d\n", T->Precision);
      for (i = 0; i < SLP->numVars; i++)
      {
        print_mp(CODIM, 0, &pt->coord[i]);
        fprintf(CODIM, "\n");
      }
      fprintf(CODIM, "\n%d\n", T->Precision);
      for (i = 0; i < SLP->numVars; i++)
      {
        print_mp(CODIM, 0, &approx->coord[i]);
        fprintf(CODIM, "\n");
      }
      fprintf(CODIM, "\n");
    }
    else
    { // print to PTS
      rV = 1;
      
      for (i = 0; i < SLP->numVars; i++)
      {
        print_mp(PTS, 0, &pt->coord[i]);
        fprintf(PTS, "\n");
      }
      fprintf(PTS, "\n");
    }

    // clear
    clear_eval_struct_mp(e);
    clear_point_mp(pt); clear_point_mp(approx);
    clear_mp(time); clear_mp(homValue); clear_mp(homValue_approx);
  }
  else
  { // use either _d or _mp
    int maxPrec = MAX(52, digits_to_prec(ceil(-log10(T->final_tol_times_mult) + 3)));
    for (i = 0; i < numComponents; i++)
    {
      maxPrec = MAX(maxPrec, subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].curr_prec);
      maxPrec = MAX(maxPrec, subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_prec);
    }

    if (maxPrec < 64)
    { // use _d
      eval_struct_d e;
      point_d pt, approx;
      comp_d time, homValue, homValue_approx;

      // initialize
      init_eval_struct_d(e, SLP->numFuncs, SLP->numVars, SLP->numPars);
      init_point_d(pt, SLP->numVars); init_point_d(approx, SLP->numVars);
      pt->size = approx->size = SLP->numVars;
      init_d(time); init_d(homValue); init_d(homValue_approx);
      set_zero_d(time);

      // increase Pts_d & Pts_mp
      increase_size_mat_d(Pts_d, totalPoints, SLP->numVars);
      increase_size_mat_mp(Pts_mp, totalPoints, SLP->numVars);
      Pts_d->rows = Pts_mp->rows = totalPoints;
      Pts_d->cols = Pts_mp->cols = SLP->numVars;

      // pick a homogenizing coordinate at random
      get_comp_rand_d(&pt->coord[0]);
      set_d(&approx->coord[0], &pt->coord[0]);

      // setup rest of pt
      currVar = 1;
      for (i = 0; i < numComponents; i++)
      { // compute normalizing value and then setup coordinates
        div_d(homValue, &pt->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_d->coord[0]);
        div_d(homValue_approx, &approx->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_d->coord[0]);
        for (j = 1; j < subWitnessSets[i].orig_variables; j++)
        {
          mul_d(&pt->coord[currVar], homValue, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_d->coord[j]);
          mul_d(&approx->coord[currVar], homValue_approx, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_d->coord[j]);
          currVar++;
        }
      }

      // setup extra coordinates
      for (i = currVar; i < SLP->numVars; i++)
      {
        mul_d(&pt->coord[i], &extra_d->coord[i - currVar], &pt->coord[0]);
        mul_d(&approx->coord[i], &extra_d->coord[i - currVar], &approx->coord[0]);
      }

      // normalize the point & copy to Pts_d
      homValue->r = infNormVec_d(pt);
      homValue_approx->r = infNormVec_d(approx);
      for (i = 0; i < SLP->numVars; i++)
      {
        mul_rdouble_d(&pt->coord[i], &pt->coord[i], homValue->r);
        set_d(&Pts_d->entry[totalPoints - 1][i], &pt->coord[i]);
        d_to_mp(&Pts_mp->entry[totalPoints - 1][i], &Pts_d->entry[totalPoints - 1][i]);
        mul_rdouble_d(&approx->coord[i], &approx->coord[i], homValue_approx->r);
      }

      // evaluate
      evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, pt, time, SLP);

      // see if vanishes sufficiently
      if (infNormVec_d(e.funcVals) < T->final_tol_times_mult)
      { // print to CODIM
        rV = 0;

        fprintf(CODIM, "52\n");
        for (i = 0; i < SLP->numVars; i++)
          fprintf(CODIM, "%.15e %.15e\n", pt->coord[i].r, pt->coord[i].i);
        fprintf(CODIM, "\n52\n");
        for (i = 0; i < SLP->numVars; i++)
          fprintf(CODIM, "%.15e %.15e\n", approx->coord[i].r, approx->coord[i].i);
        fprintf(CODIM, "\n");
      }
      else
      { // print to PTS
        rV = 1;
        
        for (i = 0; i < SLP->numVars; i++)
          fprintf(PTS, "%.15e %.15e\n", pt->coord[i].r, pt->coord[i].i);
        fprintf(PTS, "\n");
      }

      // clear
      clear_eval_struct_d(e);
      clear_point_d(pt); clear_point_d(approx);
      clear_d(time); clear_d(homValue); clear_d(homValue_approx);
    }
    else
    { // use _mp
      initMP(maxPrec);

      eval_struct_mp e;
      point_mp pt, approx;
      comp_mp time, homValue, homValue_approx, tempComp;

      // initialize
      init_eval_struct_mp(e, SLP->numFuncs, SLP->numVars, SLP->numPars);
      init_point_mp(pt, SLP->numVars); init_point_mp(approx, SLP->numVars);
      pt->size = approx->size = SLP->numVars;
      init_mp(time); init_mp(homValue); init_mp(homValue_approx); init_mp(tempComp);
      set_zero_mp(time);

      // increase Pts_mp_
      increase_size_mat_d(Pts_d, totalPoints, SLP->numVars);
      increase_size_mat_mp(Pts_mp, totalPoints, SLP->numVars);
      Pts_d->rows = Pts_mp->rows = totalPoints;
      Pts_d->cols = Pts_mp->cols = SLP->numVars;

      // pick a homogenizing coordinate at random
      get_comp_rand_mp(&pt->coord[0]);
      set_mp(&approx->coord[0], &pt->coord[0]);

      // setup rest of pt
      currVar = 1;
      for (i = 0; i < numComponents; i++)
      { // compute normalizing value and then setup coordinates
        if (subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].curr_prec < 64)
        {
          d_to_mp(tempComp, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_d->coord[0]);
          div_mp(homValue, &pt->coord[0], tempComp);

          for (j = 1; j < subWitnessSets[i].orig_variables; j++)
          {
            d_to_mp(tempComp, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_d->coord[j]);
            mul_mp(&pt->coord[currVar + j - 1], homValue, tempComp);
          }
        }
        else
        {
          div_mp(homValue, &pt->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_mp->coord[0]);

          for (j = 1; j < subWitnessSets[i].orig_variables; j++)
          {
            mul_mp(&pt->coord[currVar + j - 1], homValue, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].endPt_mp->coord[j]);
          }
        }

        if (subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_prec < 64)
        {
          d_to_mp(tempComp, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_d->coord[0]);
          div_mp(homValue_approx, &approx->coord[0], tempComp);

          for (j = 1; j < subWitnessSets[i].orig_variables; j++)
          {
            d_to_mp(tempComp, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_d->coord[j]);
            mul_mp(&approx->coord[currVar], homValue_approx, tempComp); 
            currVar++;
          }
        }
        else
        {
          div_mp(homValue_approx, &approx->coord[0], &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_mp->coord[0]);

          for (j = 1; j < subWitnessSets[i].orig_variables; j++)
          {
            mul_mp(&approx->coord[currVar], homValue_approx, &subWitnessSets[i].codim[codim_index[i]].witnessPts_amp[point_index[i]].last_approx_mp->coord[j]);
            currVar++;
          }
        }
      }

      // setup extra coordinates
      for (i = currVar; i < SLP->numVars; i++)
      {
        mul_mp(&pt->coord[i], &extra_mp->coord[i - currVar], &pt->coord[0]);
        mul_mp(&approx->coord[i], &extra_mp->coord[i - currVar], &approx->coord[0]);
      }

      // normalize the point & copy to Pts_mp
      infNormVec_mp2(homValue->r, pt);
      infNormVec_mp2(homValue_approx->r, approx);
      for (i = 0; i < SLP->numVars; i++)
      {
        mul_rmpf_mp(&pt->coord[i], &pt->coord[i], homValue->r);
        set_mp(&Pts_mp->entry[totalPoints - 1][i], &pt->coord[i]);
        mp_to_d(&Pts_d->entry[totalPoints - 1][i], &Pts_mp->entry[totalPoints - 1][i]);
        mul_rmpf_mp(&approx->coord[i], &approx->coord[i], homValue_approx->r);
      }

      // evaluate
      evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, pt, time, SLP);
  
      // see if vanishes sufficiently
      if (infNormVec_mp(e.funcVals) < T->final_tol_times_mult)
      { // print to CODIM
        rV = 0;

        fprintf(CODIM, "%d\n", maxPrec);
        for (i = 0; i < SLP->numVars; i++)
        {
          print_mp(CODIM, 0, &pt->coord[i]);
          fprintf(CODIM, "\n");
        }
        fprintf(CODIM, "\n%d\n", maxPrec);
        for (i = 0; i < SLP->numVars; i++)
        {
          print_mp(CODIM, 0, &approx->coord[i]);
          fprintf(CODIM, "\n");
        }
        fprintf(CODIM, "\n");
      }
      else
      { // print to PTS
        rV = 1;
      
        for (i = 0; i < SLP->numVars; i++)
        {
          print_mp(PTS, 0, &pt->coord[i]);
          fprintf(PTS, "\n");
        }
        fprintf(PTS, "\n");
      }

      // clear
      clear_eval_struct_mp(e);
      clear_point_mp(pt); clear_point_mp(approx);
      clear_mp(time); clear_mp(homValue); clear_mp(homValue_approx); clear_mp(tempComp);
    }
  }

  return rV;
}

void regenExtendSetup_system(int *maxCodim, int *specificCodim, double intrinsicCutoffMultiplier, int system_rank, point_d extra_d, point_mp extra_mp, int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP, preproc_data *PPD, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup RPD for extending regeneration                   *
\***************************************************************/
{ 
  int i, j, k, intrinsicCutoff = 0, currCodim = 0, origFuncs = 0, newFuncs = 0;

  // print message about what is happening
  printf("\nConstructing slices.\n");

  // compute the number of original functions and new functions, and currCodim
  for (i = 0; i < numComponents; i++)
  {
    currCodim += subWitnessSets[i].codim[codim_index[i]].codim;
    origFuncs += subWitnessSets[i].num_funcs;
  }
  newFuncs = PPD->num_funcs - origFuncs;

  // setup precision
  RPD->curr_precision = T->Precision;

  // setup the SLP
  RPD->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(RPD->Prog, SLP);
  RPD->orig_variables = SLP->numVars;

  // error checking
  if (RPD->Prog->numPathVars > 0)
  { // path variable present
    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (RPD->Prog->numPars > 0)
  { // parameter present
    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  
  // setup PPD
  cp_preproc_data(&RPD->PPD, PPD);

  // verify that we are using only 1 homogenous variable group
  if (RPD->PPD.num_hom_var_gp + RPD->PPD.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: Positive dimensional regeneration is implemented for systems with only one variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }
  
  // setup rank and number of variables for tracking
  RPD->system_rank = system_rank;
  T->numVars = RPD->new_variables = RPD->system_rank + RPD->PPD.num_var_gp + RPD->PPD.num_hom_var_gp;

  // setup the number of functions & codimension
  RPD->num_funcs = RPD->PPD.num_funcs;
  RPD->num_codim = RPD->system_rank;

  // error checking on specific codimension
  if (*specificCodim > 0 && *specificCodim > MIN(currCodim + newFuncs, RPD->num_codim))
  {
    printf("NOTE: Based on the setup, there will be no components of codimension %d!\n", *specificCodim);
    bexit(ERROR_INPUT_SYSTEM);
  }
  
  // determine where to switch from intrinsic to extrinsic
  intrinsicCutoff = floor(intrinsicCutoffMultiplier * RPD->new_variables);

  // setup degree based on the special structure
  setupDegrees_extend(&RPD->orig_degrees, &RPD->new_degrees, &RPD->P, currCodim, origFuncs, newFuncs, RPD->PPD.num_var_gp + RPD->PPD.num_hom_var_gp, "deg.out");

  // setup W
  RPD->W = (int ***)bmalloc(RPD->num_codim * sizeof(int **));
  for (i = 0; i < RPD->num_codim; i++)
  { // setup for codimension i + 1
    RPD->W[i] = (int **)bmalloc((i + 1) * sizeof(int *));
    for (j = 0; j <= i; j++)
    {
      RPD->W[i][j] = (int *)bmalloc((RPD->num_funcs - j - 1) * sizeof(int));
      for (k = 0; k < RPD->num_funcs - j - 1; k++)
      {
        RPD->W[i][j][k] = MAX(RPD->new_degrees[j] - RPD->new_degrees[j+ k + 1], 0);
      }
    }
  }

  // setup gamma
  if (T->MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(RPD->gamma_d);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    init_mp(RPD->gamma_mp);
    get_comp_rand_mp(RPD->gamma_mp);
  }
  else
  { // setup gamma_d, gamma_mp & gamma_rat
    get_comp_rand_rat(RPD->gamma_d, RPD->gamma_mp, RPD->gamma_rat, RPD->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // allocate codim
  RPD->codim = (regenCodim_t *)bmalloc(RPD->num_codim * sizeof(regenCodim_t));

  // initialize
  RPD->sameA = 0;

  // setup C (convert between old and new variables using intrinsic slicing) and coeff (slices)
  regenExtendSetup_system_coeff(extra_d, extra_mp, numComponents, codim_index, subWitnessSets, T, RPD);

  // setup H, homVarConst, patchCoeff, A
  if (T->MPType == 0)
  { // setup H_d & homVarConst_d
    init_vec_d(RPD->H_d, RPD->new_variables);
    RPD->H_d->size = RPD->new_variables;
    if (RPD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (RPD->orig_variables != RPD->new_variables)
      { // H_d is first row of C
        for (i = 0; i < RPD->new_variables; i++)
        {
          set_d(&RPD->H_d->coord[i], &RPD->C_d->entry[0][i]);
        }
      }
      else
      { // H_d = [1,0..0]
        set_one_d(&RPD->H_d->coord[0]);
        for (i = 1; i < RPD->new_variables; i++)
        {
          set_zero_d(&RPD->H_d->coord[i]);
        }
      }

      // setup homVarConst_d to be 0
      set_zero_d(RPD->homVarConst_d);
    }
    else
    { // using a homogeneous variable group

      // setup H_d to be random
      make_vec_random_d(RPD->H_d, RPD->new_variables);

      // setup homVarConst_d to be random
      get_comp_rand_d(RPD->homVarConst_d);
    }

    // setup patchCoeff_d
    init_vec_d(RPD->patchCoeff_d, RPD->new_variables);
    make_vec_random_d(RPD->patchCoeff_d, RPD->new_variables);

    // setup A_d
    RPD->sameA = 1;
    RPD->A_d = (mat_d *)bmalloc(RPD->num_codim * sizeof(mat_d));
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // setup A_d[i]
      init_mat_d(RPD->A_d[i], i + 1, RPD->num_funcs);
      RPD->A_d[i]->rows = i + 1;
      RPD->A_d[i]->cols = RPD->num_funcs;

      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_d(RPD->A_d[i], i + 1, RPD->num_funcs);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_d(&RPD->A_d[i]->entry[j][k]);
            }
            else
            { // set to 0
              set_zero_d(&RPD->A_d[i]->entry[j][k]);
            }

        // clear out extra parts that should be zero
        for (j = 0; j < currCodim; j++) // first set does not depend on new functions
          for (k = j+1; k < currCodim + newFuncs; k++)
          {
            set_zero_d(&RPD->A_d[i]->entry[j][k]);
          }
        for (j = currCodim; j <= i; j++) // last set does not depend on old functions
          for (k = currCodim + newFuncs; k < RPD->num_funcs; k++)
          {
            set_zero_d(&RPD->A_d[i]->entry[j][k]);
          }
      }
      else
      { // copy the top of A_d[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
          }
      }
    }
  }
  else if (T->MPType == 1)
  { // setup H_mp & homVarConst_mp
    init_vec_mp(RPD->H_mp, RPD->new_variables);
    RPD->H_mp->size = RPD->new_variables;
    init_mp(RPD->homVarConst_mp);
    if (RPD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (RPD->orig_variables != RPD->new_variables)
      { // H_mp is first row of C
        for (i = 0; i < RPD->new_variables; i++)
        {
          set_mp(&RPD->H_mp->coord[i], &RPD->C_mp->entry[0][i]);
        }
      }
      else
      { // H_d = [1,0..0]
        set_one_mp(&RPD->H_mp->coord[0]);
        for (i = 1; i < RPD->new_variables; i++)
        {
          set_zero_mp(&RPD->H_mp->coord[i]);
        }
      }

      // setup homVarConst_mp to be 0
      set_zero_mp(RPD->homVarConst_mp);
    }
    else
    { // using a homogeneous variable group

      // setup H_mp to be random
      make_vec_random_mp(RPD->H_mp, RPD->new_variables);

      // setup homVarConst_mp to be random
      get_comp_rand_mp(RPD->homVarConst_mp);
    }

    // setup patchCoeff_mp
    init_vec_mp(RPD->patchCoeff_mp, RPD->new_variables);
    make_vec_random_mp(RPD->patchCoeff_mp, RPD->new_variables);

    // setup A_mp
    RPD->sameA = 1;
    RPD->A_mp = (mat_mp *)bmalloc(RPD->num_codim * sizeof(mat_mp));
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // setup A_mp[i]
      init_mat_mp(RPD->A_mp[i], i + 1, RPD->num_funcs);
      RPD->A_mp[i]->rows = i + 1;
      RPD->A_mp[i]->cols = RPD->num_funcs;

      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_mp(RPD->A_mp[i], i + 1, RPD->num_funcs, T->Precision);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_mp(&RPD->A_mp[i]->entry[j][k]);
            }
            else
            { // set to 0
              set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
            }

        // clear out extra parts that should be zero
        for (j = 0; j < currCodim; j++) // first set does not depend on new functions
          for (k = j+1; k < currCodim + newFuncs; k++)
          {
            set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
          }
        for (j = currCodim; j < currCodim + newFuncs; j++) // last set does not depend on old functions
          for (k = currCodim + newFuncs; k < RPD->num_funcs; k++)
          {
            set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
          }
      }
      else
      { // copy the top of A_mp[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
          }
      }
    }
  }
  else
  { // setup H_d, H_mp, H_rat & homVarConst_d, homVarConst_mp, homVarConst_rat
    init_vec_d(RPD->H_d, RPD->new_variables);
    init_vec_mp2(RPD->H_mp, RPD->new_variables, RPD->curr_precision);
    init_vec_rat(RPD->H_rat, RPD->new_variables);
    RPD->H_d->size = RPD->H_mp->size = RPD->new_variables;
    init_mp2(RPD->homVarConst_mp, RPD->curr_precision);
    init_rat(RPD->homVarConst_rat);

    if (RPD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (RPD->orig_variables != RPD->new_variables)
      { // H is first row of C
        for (i = 0; i < RPD->new_variables; i++)
        {
          set_d(&RPD->H_d->coord[i], &RPD->C_d->entry[0][i]);
          set_mp(&RPD->H_mp->coord[i], &RPD->C_mp->entry[0][i]);
          set_rat(RPD->H_rat[i], RPD->C_rat[0][i]);
        }
      }
      else
      { // H_d = [1,0..0]
        set_one_d(&RPD->H_d->coord[0]);
        set_one_mp(&RPD->H_mp->coord[0]);
        set_one_rat(RPD->H_rat[0]);
        for (i = 1; i < RPD->new_variables; i++)
        {
          set_zero_d(&RPD->H_d->coord[i]);
          set_zero_mp(&RPD->H_mp->coord[i]);
          set_zero_rat(RPD->H_rat[i]);
        }
      }

      // setup homVarConst to be 0
      set_zero_d(RPD->homVarConst_d);
      set_zero_mp(RPD->homVarConst_mp);
      set_zero_rat(RPD->homVarConst_rat);
    }
    else
    { // using a homogeneous variable group

      // setup H to be random
      make_vec_random_rat(RPD->H_d, RPD->H_mp, RPD->H_rat, RPD->new_variables, RPD->curr_precision, T->AMP_max_prec, 0, 0);

      // setup homVarConst to be random
      get_comp_rand_rat(RPD->homVarConst_d, RPD->homVarConst_mp, RPD->homVarConst_rat, RPD->curr_precision, T->AMP_max_prec, 0, 0);
    }

    // setup patchCoeff
    init_vec_d(RPD->patchCoeff_d, RPD->new_variables);
    init_vec_mp2(RPD->patchCoeff_mp, RPD->new_variables, RPD->curr_precision);
    init_vec_rat(RPD->patchCoeff_rat, RPD->new_variables);
    make_vec_random_rat(RPD->patchCoeff_d, RPD->patchCoeff_mp, RPD->patchCoeff_rat, RPD->new_variables, RPD->curr_precision, T->AMP_max_prec, 0, 0);

    // setup A
    RPD->sameA = 1;
    RPD->A_d = (mat_d *)bmalloc(RPD->num_codim * sizeof(mat_d));
    RPD->A_mp = (mat_mp *)bmalloc(RPD->num_codim * sizeof(mat_mp));
    RPD->A_rat = (mpq_t ****)bmalloc(RPD->num_codim * sizeof(mpq_t ***));
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // setup A[i]
      init_mat_d(RPD->A_d[i], i + 1, RPD->num_funcs);
      init_mat_mp2(RPD->A_mp[i], i + 1, RPD->num_funcs, RPD->curr_precision);
      init_mat_rat(RPD->A_rat[i], i + 1, RPD->num_funcs);
      RPD->A_d[i]->rows = RPD->A_mp[i]->rows = i + 1;
      RPD->A_d[i]->cols = RPD->A_mp[i]->cols = RPD->num_funcs;

      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_rat(RPD->A_d[i], RPD->A_mp[i], RPD->A_rat[i], i + 1, RPD->num_funcs, RPD->curr_precision, T->AMP_max_prec, 0, 0);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_d(&RPD->A_d[i]->entry[j][k]);
              set_one_mp(&RPD->A_mp[i]->entry[j][k]);
              set_one_rat(RPD->A_rat[i][j][k]);
            }
            else
            { // set to 0
              set_zero_d(&RPD->A_d[i]->entry[j][k]);
              set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
              set_zero_rat(RPD->A_rat[i][j][k]);
            }

        // clear out extra parts that should be zero
        for (j = 0; j < currCodim; j++) // first set does not depend on new functions
          for (k = j+1; k < currCodim + newFuncs; k++)
          {
            set_zero_d(&RPD->A_d[i]->entry[j][k]);
            set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
            set_zero_rat(RPD->A_rat[i][j][k]);
          }
        for (j = currCodim; j <= i; j++) // last set does not depend on old functions
          for (k = currCodim + newFuncs; k < RPD->num_funcs; k++)
          {
            set_zero_d(&RPD->A_d[i]->entry[j][k]);
            set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
            set_zero_rat(RPD->A_rat[i][j][k]);
          }
      }
      else
      { // copy the top of A_d[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
            set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
            set_rat(RPD->A_rat[i][j][k], RPD->A_rat[i+1][j][k]);
          }
      }
    }
  }

  // initialize the codimensions and construct intrinsic slicing, if needed
  for (i = 0; i < RPD->num_codim; i++)
  { // setup the ith codimData
    setupRegenCodimData(RPD, i, i + 1, T->MPType, T->AMP_max_prec, intrinsicCutoff);
  }

  // print message about codimensions
  currCodim = MIN(RPD->num_codim, currCodim + newFuncs);
  if (*specificCodim > 0)
  { // print a message
    printf("\nNOTE: You have requested to compute only codimension %d.\n", *specificCodim);
  }
  else if (*maxCodim > 0 && *maxCodim < currCodim)
  { // print a message
    printf("\nNOTE: You have requested a maximum codimension of %d.\n", *maxCodim);
  }
  // setup max codimension based on current codim and number of new functions
  *maxCodim = *maxCodim > 0 ? MIN(*maxCodim, currCodim) : currCodim;

  return;
}

int regenExtendSetup_points_finish(int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: curr_codim_index                               *
* NOTES: finish setuping the points for extending regeneration  *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, curr_codim = 0, size, numPoints;
  char *str = NULL;
  FILE *CODIM = NULL;

  // print message about what is happening
  printf("\nConstructing start points.\n");

  // setup
  for (i = 0; i < numComponents; i++)
    curr_codim += subWitnessSets[i].codim[codim_index[i]].codim;

  // setup witness_superset for current codimension 
  size = 1 + snprintf(NULL, 0, "witness_superset_%d", curr_codim);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "witness_superset_%d", curr_codim);
  CODIM = fopen(str, "r");
  if (CODIM == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", str);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in number of points
  numPoints = 0;
  fscanf(CODIM, "%d", &numPoints);
  fclose(CODIM);

  if (numPoints <= 0)
  { // create empty file
    CODIM = fopen(str, "w");
    fclose(CODIM);
  }
  else
  { // setup the witness superset properly
    rename(str, "witness_superset_test");

    RPD->curr_codim = curr_codim - 1;

    // setup wtiness superset from points
    regenExtendSetup_witness("witness_superset_test", str, T, RPD);

    // delete temp file
    remove("witness_superset_test");

    // make sure RPD knows that this codim has nonsingular witness points
    RPD->codim[curr_codim - 1].num_superset = RPD->codim[curr_codim - 1].num_nonsing = numPoints;
  }

  // setup startRPD for preparing the next codimension
  size = 1 + snprintf(NULL, 0, "startRPD_%d", curr_codim+1);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "startRPD_%d", curr_codim+1);
  RPD->codim[curr_codim - 1].num_nonsolns = regenExtendSetup_start(str, "startRPD_test", curr_codim - 1, T, RPD);

  free(str);

  return (curr_codim - 1); 
}

int regenExtendSetup_start(char *origFile, char *newFile, int curr_codim_index, tracker_config_t *T, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of points                               *
* NOTES: setup the next start points                            *
\***************************************************************/
{
  int i, j, numPoints = 0;
  point_d pt_d;
  point_mp pt_mp;
  mat_d tr_d, B_transpose_d;
  mat_mp tr_mp, B_transpose_mp;
  FILE *IN = fopen(origFile, "r"), *OUT = fopen(newFile, "w");

  init_point_d(pt_d, RPD->orig_variables);
  init_point_mp(pt_mp, RPD->orig_variables);
  init_mat_d(tr_d, 0, 0);
  init_mat_d(B_transpose_d, 0, 0);
  init_mat_mp(tr_mp, 0, 0);
  init_mat_mp(B_transpose_mp, 0, 0);

  // loop over the points -- print to OUT
  fscanf(IN, "%d", &numPoints);

  if (T->MPType == 0 || T->MPType == 2)
  { // setup using _d
    transpose_d(tr_d, RPD->C_d);

    if (RPD->codim[curr_codim_index].useIntrinsicSlice)
    { // convert to intrinsic coords
      transpose_d(B_transpose_d, RPD->codim[curr_codim_index].B_d);
    }

    for (i = 0; i < numPoints; i++)
    { // read in point 
      increase_size_point_d(pt_d, RPD->orig_variables);
      pt_d->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
        fscanf(IN, "%lf%lf", &pt_d->coord[j].r, &pt_d->coord[j].i);

      // convert to new variables, if needed
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_d(pt_d, tr_d, pt_d);
      }

      move_to_patch_vec_d(pt_d, pt_d, RPD->patchCoeff_d);

      // convert to the next set of coordinates
      if (RPD->codim[curr_codim_index].useIntrinsicSlice)
      { // convert to intrinsic coords
        extrinsicToIntrinsic_d(pt_d, pt_d, B_transpose_d, RPD->codim[curr_codim_index].p_d);
      }

      // print data to OUT
      fprintf(OUT, "%d\n%d\n", i, 52); 
      for (j = 0; j < pt_d->size; j++)
        fprintf(OUT, "%.15e %.15e\n", pt_d->coord[j].r, pt_d->coord[j].i);
      fprintf(OUT, "\n");
    }
  }
  else
  { // setup using _mp
    transpose_mp(tr_mp, RPD->C_mp);

    if (RPD->codim[curr_codim_index].useIntrinsicSlice)
    { // convert to intrinsic coords
      transpose_mp(B_transpose_mp, RPD->codim[curr_codim_index].B_mp);
    }

    for (i = 0; i < numPoints; i++)
    { // read in point
      increase_size_point_mp(pt_mp, RPD->orig_variables);
      pt_mp->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
      {
        mpf_inp_str(pt_mp->coord[j].r, IN, 10);
        mpf_inp_str(pt_mp->coord[j].i, IN, 10);
        scanRestOfLine(IN);
      }

      // convert to new variables, if needed
      if (RPD->new_variables != RPD->orig_variables)
      {
        mul_mat_vec_mp(pt_mp, tr_mp, pt_mp);
      }

      move_to_patch_vec_mp(pt_mp, pt_mp, RPD->patchCoeff_mp);

      // convert to the next set of coordinates
      if (RPD->codim[curr_codim_index].useIntrinsicSlice)
      { // convert to intrinsic coords
        extrinsicToIntrinsic_mp(pt_mp, pt_mp, B_transpose_mp, RPD->codim[curr_codim_index].p_mp);
      }

      // print data to OUT
      fprintf(OUT, "%d\n%d\n", i, T->Precision); 
      for (j = 0; j < pt_mp->size; j++)
      {
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].r);
        fprintf(OUT, " ");
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].i);
        fprintf(OUT, "\n");
      }
      fprintf(OUT, "\n");
    }
  }

  fclose(IN);
  fclose(OUT);

  clear_point_d(pt_d);
  clear_point_mp(pt_mp);
  clear_mat_d(tr_d);
  clear_mat_d(B_transpose_d);
  clear_mat_mp(tr_mp);
  clear_mat_mp(B_transpose_mp);

  return numPoints;
}

void regenExtendSetup_witness(char *origFile, char *newFile, tracker_config_t *T, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the witness superset                             *
\***************************************************************/
{
  int i, j, numPoints = 0, pt_prec, last_prec, corank, codim_index = RPD->curr_codim;
  double CN, smallest_nonzero_SV, largest_zero_SV;
  comp_d final_d;
  comp_mp final_mp;
  point_d pt_d, last_d, int_d, int_last_d;
  point_mp pt_mp, last_mp, int_mp, int_last_mp;
  mat_d tr_d;
  mat_mp tr_mp;
  eval_struct_d e_d;
  eval_struct_mp e_mp;
  FILE *IN = fopen(origFile, "r"), *OUT = fopen(newFile, "w"), *TEMP = fopen("temp_file", "w");

  // initialize
  init_d(final_d);
  init_mp(final_mp);
  init_point_d(pt_d, RPD->orig_variables);
  init_point_d(last_d, RPD->orig_variables);
  init_point_mp(pt_mp, RPD->orig_variables);
  init_point_mp(last_mp, RPD->orig_variables);
  init_point_d(int_d, 0);
  init_point_d(int_last_d, 0);
  init_point_mp(int_mp, 0);
  init_point_mp(int_last_mp, 0);
  init_mat_d(tr_d, 0, 0);
  init_mat_mp(tr_mp, 0, 0);
  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);
  pt_d->size = last_d->size = pt_mp->size = last_mp->size = RPD->orig_variables;
  set_zero_d(final_d);
  set_zero_mp(final_mp);

  // loop over the points -- print to OUT
  fscanf(IN, "%d", &numPoints);
  for (i = 0; i < numPoints; i++)
  { // read in the precision for the point
    fscanf(IN, "%d", &pt_prec);

    if (pt_prec < 64)
    { // setup _d
      increase_size_point_d(pt_d, RPD->orig_variables);
      pt_d->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
        fscanf(IN, "%lf%lf", &pt_d->coord[j].r, &pt_d->coord[j].i);
    }
    else 
    { // setup _mp
      setprec_point_mp(pt_mp, pt_prec);
      increase_size_point_mp(pt_mp, RPD->orig_variables);
      pt_mp->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
      {
        mpf_inp_str(pt_mp->coord[j].r, IN, 10);
        mpf_inp_str(pt_mp->coord[j].i, IN, 10);
        scanRestOfLine(IN);
      }
    }

    // read in the precision for the last approximation
    fscanf(IN, "%d", &last_prec);

    if (last_prec < 64)
    { // setup _d
      increase_size_point_d(last_d, RPD->orig_variables);
      last_d->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
        fscanf(IN, "%lf%lf", &last_d->coord[j].r, &last_d->coord[j].i);
    }
    else
    { // setup _mp
      setprec_point_mp(last_mp, pt_prec);
      increase_size_point_mp(last_mp, RPD->orig_variables);
      last_mp->size = RPD->orig_variables;
      for (j = 0; j < RPD->orig_variables; j++)
      {
        mpf_inp_str(last_mp->coord[j].r, IN, 10);
        mpf_inp_str(last_mp->coord[j].i, IN, 10);
        scanRestOfLine(IN);
      }
    } 

    // convert to new variables, if needed
    if (RPD->new_variables != RPD->orig_variables)
    {
      if (T->MPType == 0)
      { // use _d
        transpose_d(tr_d, RPD->C_d);
        mul_mat_vec_d(pt_d, tr_d, pt_d);
        mul_mat_vec_d(last_d, tr_d, last_d);
      }
      else if (T->MPType == 1)
      { // use _mp
        transpose_mp(tr_mp, RPD->C_mp);
        mul_mat_vec_mp(pt_mp, tr_mp, pt_mp);
        mul_mat_vec_mp(last_mp, tr_mp, last_mp);
      }
      else
      { // AMP
        if (pt_prec < 64)
        { // use _d
          transpose_d(tr_d, RPD->C_d);
          mul_mat_vec_d(pt_d, tr_d, pt_d);
        }
        else
        { // use _mp
          change_regen_pos_dim_prec(RPD, pt_prec);
          setprec_mat_mp(tr_mp, pt_prec);
          transpose_mp(tr_mp, RPD->C_mp);
          mul_mat_vec_mp(pt_mp, tr_mp, pt_mp);
        }

        if (last_prec < 64)
        { // use _d
          transpose_d(tr_d, RPD->C_d);
          mul_mat_vec_d(last_d, tr_d, last_d);
        }
        else
        { // use _mp
          change_regen_pos_dim_prec(RPD, last_prec);
          setprec_mat_mp(tr_mp, last_prec);
          transpose_mp(tr_mp, RPD->C_mp);
          mul_mat_vec_mp(last_mp, tr_mp, last_mp);
        }
      }
    }

    // move to patch
    if (T->MPType == 0)
    {
      move_to_patch_vec_d(pt_d, pt_d, RPD->patchCoeff_d);
      move_to_patch_vec_d(last_d, last_d, RPD->patchCoeff_d);
    }
    else if (T->MPType == 1)
    {
      move_to_patch_vec_mp(pt_mp, pt_mp, RPD->patchCoeff_mp);
      move_to_patch_vec_mp(last_mp, last_mp, RPD->patchCoeff_mp);
    }
    else
    {
      if (pt_prec < 64)
      {
        move_to_patch_vec_d(pt_d, pt_d, RPD->patchCoeff_d);
      }
      else
      {
        change_regen_pos_dim_prec(RPD, pt_prec);
        move_to_patch_vec_mp(pt_mp, pt_mp, RPD->patchCoeff_mp);
      }

      if (last_prec < 64)
      {
        move_to_patch_vec_d(last_d, last_d, RPD->patchCoeff_d);
      }
      else
      {
        change_regen_pos_dim_prec(RPD, last_prec);
        move_to_patch_vec_mp(last_mp, last_mp, RPD->patchCoeff_mp);
      }
    }

    // convert to intrinsic coordinates, if needed
    if (RPD->codim[codim_index].useIntrinsicSlice)
    {
      if (pt_prec < 64)
      {
        transpose_d(tr_d, RPD->codim[codim_index].B_d);
        extrinsicToIntrinsic_d(int_d, pt_d, tr_d, RPD->codim[codim_index].p_d);
      }
      else
      {
        transpose_mp(tr_mp, RPD->codim[codim_index].B_mp);
        extrinsicToIntrinsic_mp(int_mp, pt_mp, tr_mp, RPD->codim[codim_index].p_mp);
      }

      if (last_prec < 64)
      {
        transpose_d(tr_d, RPD->codim[codim_index].B_d);
        extrinsicToIntrinsic_d(int_last_d, last_d, tr_d, RPD->codim[codim_index].p_d);
      }
      else
      {
        transpose_mp(tr_mp, RPD->codim[codim_index].B_mp);
        extrinsicToIntrinsic_mp(int_last_mp, last_mp, tr_mp, RPD->codim[codim_index].p_mp);
      }
    }
    else
    {
      if (pt_prec < 64)
      {
        point_cp_d(int_d, pt_d);
      }
      else
      {
        point_cp_mp(int_mp, pt_mp);
      }

      if (last_prec < 64)
      {
        point_cp_d(int_last_d, last_d);
      }
      else
      {
        point_cp_mp(int_last_mp, last_mp);
      }
    }

    // compute the data
    if (T->MPType == 0)
    {
      corank = Cauchy_corank_d(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_d, int_last_d, final_d, T, TEMP, &e_d, RPD, regen_pos_dim_eval_d);
    }
    else if (T->MPType == 1)
    {
      corank = Cauchy_corank_mp(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_mp, int_last_mp, final_mp, T, TEMP, &e_mp, RPD, regen_pos_dim_eval_mp);
    }
    else
    {
      corank = Cauchy_corank_amp(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_d, int_mp, pt_prec, int_last_d, int_last_mp, last_prec, final_d, final_mp, T, TEMP, &e_d, &e_mp, RPD, RPD, regen_pos_dim_eval_d, regen_pos_dim_eval_mp, change_regen_pos_dim_prec);
    }

    // print data to OUT
    fprintf(OUT, "%d\n", pt_prec);
    if (pt_prec < 64)
    { // print pt_d
      fprintf(OUT, "%d\n", pt_d->size);
      for (j = 0; j < pt_d->size; j++)
        fprintf(OUT, "%.15e %.15e\n", pt_d->coord[j].r, pt_d->coord[j].i);
    }
    else
    { // print pt_mp
      fprintf(OUT, "%d\n", pt_mp->size);
      for (j = 0; j < pt_mp->size; j++)
      {
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].r);
        fprintf(OUT, " ");
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].i);
        fprintf(OUT, "\n");
      }
    }
    // print prec for last_approx
    fprintf(OUT, "%d\n", last_prec);
    if (last_prec < 64)
    { // print last_d
      for (j = 0; j < last_d->size; j++)
        fprintf(OUT, "%.15e %.15e\n", last_d->coord[j].r, last_d->coord[j].i);
    }
    else
    { // print last_mp
      for (j = 0; j < last_mp->size; j++)
      {
        mpf_out_str(OUT, 10, 0, last_mp->coord[j].r);
        fprintf(OUT, " ");
        mpf_out_str(OUT, 10, 0, last_mp->coord[j].i);
        fprintf(OUT, "\n");
      }
    }

    if (pt_prec < 64)
    { // print final_d
      fprintf(OUT, "%.15e %.15e\n", final_d->r, final_d->i);
    }
    else
    { // print final_mp
      mpf_out_str(OUT, 10, 0, final_mp->r);
      fprintf(OUT, " ");
      mpf_out_str(OUT, 10, 0, final_mp->i);
      fprintf(OUT, "\n");
    }

    // print other data
    fprintf(OUT, "%.15e\n%d\n%.15e\n%.15e\n%d\n", CN, corank, smallest_nonzero_SV, largest_zero_SV, 0);
  }

  fclose(IN);
  fclose(OUT);
  fclose(TEMP);
  remove("temp_file");

  clear_d(final_d);
  clear_mp(final_mp);
  clear_point_d(pt_d);
  clear_point_d(last_d);
  clear_point_mp(pt_mp); 
  clear_point_mp(last_mp);
  clear_point_d(int_d);
  clear_point_d(int_last_d);
  clear_point_mp(int_mp);
  clear_point_mp(int_last_mp);
  clear_mat_d(tr_d);
  clear_mat_mp(tr_mp);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  return;
}

void regenExtendSetup_points(point_d extra_d, point_mp extra_mp, mat_d Points_d, mat_mp Points_mp, int numComponents, int *codim_index, int *component_number, witness_t *subWitnessSets, tracker_config_t *T, prog_t *SLP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the points for extending regeneration            *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, isGood, extra_start = 0, extra_end = 0, curr_codim = 0, size = 0, numPts = 0, numStartPts = 0;
  int *curr_index = (int *)bmalloc(numComponents * sizeof(int)), *max_index = (int *)bmalloc(numComponents * sizeof(int));
  char *str = NULL;
  FILE *CODIM = NULL, *PTS = NULL;

  // print message about what is happening
  printf("\nSorting points.\n");

  // setup 
  extra_start = 1;
  for (i = 0; i < numComponents; i++)
  {
    curr_codim += subWitnessSets[i].codim[codim_index[i]].codim;
    extra_start += subWitnessSets[i].orig_variables - 1;
    curr_index[i] = 0;
    max_index[i] = subWitnessSets[i].codim[codim_index[i]].num_set;
  }
  extra_end = SLP->numVars;

  // fix the extra coordinates randomly
  make_vec_random_mp(extra_mp, extra_end - extra_start);
  point_mp_to_d(extra_d, extra_mp);

  // setup witness_superset_* for lower codimensions -- empty file
  for (i = 1; i < curr_codim; i++)
  {
    size = 1 + snprintf(NULL, 0, "witness_superset_%d", i);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "witness_superset_%d", i);
    CODIM = fopen(str, "w");
    fclose(CODIM);
  }

  // setup witness_superset for current codimension and points to track for the regenerative cascade
  size = 1 + snprintf(NULL, 0, "witness_superset_%d", curr_codim);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "witness_superset_%d", curr_codim);
  CODIM = fopen(str, "w");
  fprintf(CODIM, "                              \n\n");

  size = 1 + snprintf(NULL, 0, "startRPD_%d", curr_codim+1);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "startRPD_%d", curr_codim+1);
  PTS = fopen(str, "w+");
  fprintf(PTS, "                              \n\n");

  // loop over the points -- after one is found, test it and print appropriately
  while (1)
  { // check the current point -> verify each is on correct component
    isGood = 1;
    for (i = 0; i < numComponents && isGood; i++)
      if (!(curr_index[i] < max_index[i] && (component_number[i] == -2 || subWitnessSets[i].codim[codim_index[i]].component_nums[curr_index[i]] == component_number[i])))
        isGood = 0;

    if (isGood)
    { // this point is good
      numPts++;

      // print to CODIM or PTS
      numStartPts += regenExtendSetup_point(numPts, Points_d, Points_mp, numComponents, codim_index, curr_index, subWitnessSets, T, SLP, CODIM, PTS, extra_d, extra_mp);
    }

    // update curr_index
    for (i = numComponents - 1; i >= 0; i--)
    {
      if (curr_index[i] >= max_index[i] - 1)
        curr_index[i] = 0;
      else
      { // increment and leave for loop
        curr_index[i]++;
        break;
      }
    }

    // check if complete -- exit while loop
    if (i < 0)
      break;
  }
  
  // finish up
  rewind(CODIM);
  fprintf(CODIM, "%d", numPts - numStartPts);
  rewind(PTS);
  fprintf(PTS, "%d", numStartPts);

  // close & clear
  fclose(PTS);
  fclose(CODIM);

  return;
}

void regenExtendVerification_Jacobian(int numComponents, witness_t *witnessSets, tracker_config_t *T, prog_t *SLP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify the Jacobian for extending regeneration         *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, j, k, start_func, end_func, start_vars, end_vars;

  // compute a jacobian matrix and make sure the functions have the proper dependencies
  if (T->MPType == 0 || T->MPType == 2)
  { // perform using _d
    eval_struct_d e_d;
    point_d pt_d;
    comp_d time_d;

    init_eval_struct_d(e_d, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_d(pt_d, SLP->numVars);
    init_d(time_d);

    // pt = random, time = 0
    make_vec_random_d(pt_d, SLP->numVars);
    set_zero_d(time_d);

    // evaluate
    evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, pt_d, time_d, SLP);

    // analyze Jacobian based on the functions
    start_func = 0;
    start_vars = 1; // ignore homogenizing variable
    for (j = 0; j < numComponents; j++)
    { // look at the next section of functions and make sure only depend on correct variables
      end_func = start_func + witnessSets[j].num_funcs;
      end_vars = start_vars + witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP

      // verify within bounds
      if (end_func > SLP->numFuncs)
      {
        printf("\n\nERROR: Too many functions prescribed by the known components!\n\n");
        bexit(ERROR_CONFIGURATION);
      }
      if (end_vars > SLP->numVars) 
      {
        printf("\n\nERROR: Too many variables prescribed by the known components!\n\n");
        bexit(ERROR_CONFIGURATION);
      }

      for (i = start_func; i < end_func; i++)
        for (k = 1; k < SLP->numVars; k++) // ignore homogenizing coordinate
          if (k < start_vars || k >= end_vars)
          { // this entry should be exactly zero!
            if (e_d.Jv->entry[i][k].r != 0 || e_d.Jv->entry[i][k].i != 0)
            {

              printf("\n\nERROR: The functions corresponding to system %d depend on improper variables!\n\n", j+1);
              bexit(ERROR_CONFIGURATION);
            }
          } 

      // update
      start_func = end_func;
      start_vars = end_vars;
    }

    // clear
    clear_eval_struct_d(e_d);
    clear_point_d(pt_d);
    clear_d(time_d);
  }
  else
  { // perform using _mp 
    eval_struct_mp e_mp;
    point_mp pt_mp;
    comp_mp time_mp;

    init_eval_struct_mp(e_mp, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_mp(pt_mp, SLP->numVars);
    init_mp(time_mp);

    // pt = random, time = 0
    make_vec_random_mp(pt_mp, SLP->numVars);
    set_zero_mp(time_mp);

    // evaluate
    evalProg_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, pt_mp, time_mp, SLP);

    // analyze Jacobian based on the functions
    start_func = 0;
    start_vars = 1; // ignore homogenizing variable
    for (j = 0; j < numComponents; j++)
    { // look at the next section of functions and make sure only depend on correct variables
      end_func = start_func + witnessSets[j].num_funcs;
      end_vars = start_vars + witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP

      // verify within bounds
      if (end_func > SLP->numFuncs)
      {
        printf("\n\nERROR: Too many functions prescribed by the known components!\n\n");
        bexit(ERROR_CONFIGURATION);
      }
      if (end_vars > SLP->numVars) // -1 for homogeneous coordinate
      {
        printf("\n\nERROR: Too many variables prescribed by the known components!\n\n");
        bexit(ERROR_CONFIGURATION);
      }

      for (i = start_func; i < end_func; i++)
        for (k = 1; k < SLP->numVars; k++) // ignore homogenizing coordinate
          if (k < start_vars || k >= end_vars)
          { // this entry should be exactly zero!
            if (!(mpfr_zero_p(e_mp.Jv->entry[i][k].r) && mpfr_zero_p(e_mp.Jv->entry[i][k].i)))
            {

              printf("\n\nERROR: The functions corresponding to system %d depend on improper variables!\n\n", j+1);
              bexit(ERROR_CONFIGURATION);
            }
          } 

      // update
      start_func = end_func;
      start_vars = end_vars;
    }

    // clear
    clear_eval_struct_mp(e_mp);
    clear_point_mp(pt_mp);
    clear_mp(time_mp);
  }

  return;
}

void regenExtendVerification_point(int currNumber, int currCodimIndex, int currPoint, int numComponents, witness_t *witnessSets, tracker_config_t *T, prog_t *SLP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify the point for extending regeneration            *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, j, start_func, end_func, start_vars, end_vars;

  if (T->MPType == 0)
  { // use _d
    eval_struct_d e_d;
    point_d pt_d;
    comp_d time_d;

    init_eval_struct_d(e_d, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_d(pt_d, SLP->numVars);
    init_d(time_d);

    // set time = 0,  pt = initial set to random numbers and then adjust the specific coordinates
    set_zero_d(time_d);
    make_vec_random_d(pt_d, SLP->numVars);
    // set homogenizing coordinate
    set_d(&pt_d->coord[0], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_d[currPoint].endPt->coord[0]);

    start_func = 0;
    start_vars = 1; // ignore homogenizing variable
    for (j = 0; j < currNumber; j++)
    { // look at the next section of functions and make sure only depend on correct variables
      start_func += witnessSets[j].num_funcs;
      start_vars += witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP
    }
    end_func = start_func + witnessSets[currNumber].num_funcs;
    end_vars = start_vars + witnessSets[currNumber].orig_variables - 1;

    // update the coordinates
    for (i = start_vars; i < end_vars; i++)
    {
      set_d(&pt_d->coord[i], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_d[currPoint].endPt->coord[i - start_vars + 1]);
    }

    // evaluate 
    evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, pt_d, time_d, SLP);

    // check the functions
    for (i = start_func; i < end_func; i++)
      if (d_abs_d(&e_d.funcVals->coord[i]) > T->final_tol_times_mult)
      {
        printf("\nIt appears that a point for system %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", currNumber + 1, d_abs_d(&e_d.funcVals->coord[i]), T->final_tol_times_mult);
        bexit(ERROR_CONFIGURATION);
      }

    // clear
    clear_eval_struct_d(e_d);
    clear_point_d(pt_d);
    clear_d(time_d);
  }
  else if (T->MPType == 1)
  { // use _mp
    eval_struct_mp e_mp;
    point_mp pt_mp;
    comp_mp time_mp;

    init_eval_struct_mp(e_mp, SLP->numFuncs, SLP->numVars, SLP->numPars);
    init_point_mp(pt_mp, SLP->numVars);
    init_mp(time_mp);

    // set time = 0,  pt = initial set to random numbers and then adjust the specific coordinates
    set_zero_mp(time_mp);
    make_vec_random_mp(pt_mp, SLP->numVars);
    // set homogenizing coordinate
    set_mp(&pt_mp->coord[0], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_mp[currPoint].endPt->coord[0]);

    start_func = 0;
    start_vars = 1; // ignore homogenizing variable
    for (j = 0; j < currNumber; j++)
    { // look at the next section of functions and make sure only depend on correct variables
      start_func += witnessSets[j].num_funcs;
      start_vars += witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP
    }
    end_func = start_func + witnessSets[currNumber].num_funcs;
    end_vars = start_vars + witnessSets[currNumber].orig_variables - 1;

    // update the coordinates
    for (i = start_vars; i < end_vars; i++)
    {
      set_mp(&pt_mp->coord[i], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_mp[currPoint].endPt->coord[i - start_vars + 1]);
    }

    // evaluate 
    evalProg_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, pt_mp, time_mp, SLP);

    // check the functions
    for (i = start_func; i < end_func; i++)
      if (d_abs_mp(&e_mp.funcVals->coord[i]) > T->final_tol_times_mult)
      {
        printf("\nIt appears that a point for system %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", currNumber + 1, d_abs_mp(&e_mp.funcVals->coord[i]), T->final_tol_times_mult);
        bexit(ERROR_CONFIGURATION);
      }

    // clear
    clear_eval_struct_mp(e_mp);
    clear_point_mp(pt_mp);
    clear_mp(time_mp);
  }
  else 
  { // use appropriate selection based on tolerance and precision of point
    int pt_prec = witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].curr_prec;
    int tol_prec = digits_to_prec(ceil(-log10(T->final_tol_times_mult) + 3));
    int prec = MAX(pt_prec,tol_prec);

    if (prec < 64)
    { // use _d
      eval_struct_d e_d;
      point_d pt_d;
      comp_d time_d;

      init_eval_struct_d(e_d, SLP->numFuncs, SLP->numVars, SLP->numPars);
      init_point_d(pt_d, SLP->numVars);
      init_d(time_d);

      // set time = 0,  pt = initial set to random numbers and then adjust the specific coordinates
      set_zero_d(time_d);
      make_vec_random_d(pt_d, SLP->numVars);
      // set homogenizing coordinatei -- all in double precision
      set_d(&pt_d->coord[0], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_d->coord[0]);

      start_func = 0;
      start_vars = 1; // ignore homogenizing variable
      for (j = 0; j < currNumber; j++)
      { // look at the next section of functions and make sure only depend on correct variables
        start_func += witnessSets[j].num_funcs;
        start_vars += witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP
      }
      end_func = start_func + witnessSets[currNumber].num_funcs;
      end_vars = start_vars + witnessSets[currNumber].orig_variables - 1;

      // update the coordinates
      for (i = start_vars; i < end_vars; i++)
      {
        set_d(&pt_d->coord[i], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_d->coord[i - start_vars + 1]);
      }

      // evaluate 
      evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, pt_d, time_d, SLP);

      // check the functions
      for (i = start_func; i < end_func; i++)
        if (d_abs_d(&e_d.funcVals->coord[i]) > T->final_tol_times_mult)
        {
          printf("\nIt appears that a point for system %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", currNumber + 1, d_abs_d(&e_d.funcVals->coord[i]), T->final_tol_times_mult);
          bexit(ERROR_CONFIGURATION);
        }

      // clear
      clear_eval_struct_d(e_d);
      clear_point_d(pt_d);
      clear_d(time_d);
    }
    else
    { // use _mp
      initMP(prec);

      eval_struct_mp e_mp;
      point_mp pt_mp;
      comp_mp time_mp;

      init_eval_struct_mp(e_mp, SLP->numFuncs, SLP->numVars, SLP->numPars);
      init_point_mp(pt_mp, SLP->numVars);
      init_mp(time_mp);

      // set time = 0,  pt = initial set to random numbers and then adjust the specific coordinates
      set_zero_mp(time_mp);
      make_vec_random_mp(pt_mp, SLP->numVars);
      // set homogenizing coordinate
      if (pt_prec < 64)
      {
        d_to_mp(&pt_mp->coord[0], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_d->coord[0]);
      }
      else  
      {
        set_mp(&pt_mp->coord[0], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_mp->coord[0]);
      }

      start_func = 0;
      start_vars = 1; // ignore homogenizing variable
      for (j = 0; j < currNumber; j++)
      { // look at the next section of functions and make sure only depend on correct variables
        start_func += witnessSets[j].num_funcs;
        start_vars += witnessSets[j].orig_variables - 1; // ASSUME EACH SUBSYSTEM HAS ONE AFFINE VARIABLE GROUP
      }
      end_func = start_func + witnessSets[currNumber].num_funcs;
      end_vars = start_vars + witnessSets[currNumber].orig_variables - 1;

      // update the coordinates
      if (pt_prec < 64)
      {
        for (i = start_vars; i < end_vars; i++)
        {
          d_to_mp(&pt_mp->coord[i], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_d->coord[i - start_vars + 1]);
        }
      }
      else  
      {
        for (i = start_vars; i < end_vars; i++)
        {
          set_mp(&pt_mp->coord[i], &witnessSets[currNumber].codim[currCodimIndex].witnessPts_amp[currPoint].endPt_mp->coord[i - start_vars + 1]);
        }
      }

      // evaluate 
      evalProg_mp(e_mp.funcVals, e_mp.parVals, e_mp.parDer, e_mp.Jv, e_mp.Jp, pt_mp, time_mp, SLP);

      // check the functions
      for (i = start_func; i < end_func; i++)
        if (d_abs_mp(&e_mp.funcVals->coord[i]) > T->final_tol_times_mult)
        {
          printf("\nIt appears that a point for system %d does not sufficiently satisfy the original system (residual: %e > tolerance: %e).\n", currNumber + 1, d_abs_mp(&e_mp.funcVals->coord[i]), T->final_tol_times_mult);
          bexit(ERROR_CONFIGURATION);
        }

      // clear
      clear_eval_struct_mp(e_mp);
      clear_point_mp(pt_mp);
      clear_mp(time_mp);
    }
  }

  return;
}

void regenExtendVerification(int numComponents, int *codim_index, int *component_number, witness_t *witnessSets, tracker_config_t *T, prog_t *SLP, int pathMod)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify the setup for extending regeneration            *
\***************************************************************/
{ // ASSUME ONE AFFINE VARIABLE GROUP
  int i, j, num_points;

  // display message about what is happening
  if (numComponents == 1)
    printf("\nVerifying the component.\n");
  else
    printf("\nVerifying the %d components.\n", numComponents);
 
  // compute a jacobian matrix and make sure the functions have the proper dependencies
  regenExtendVerification_Jacobian(numComponents, witnessSets, T, SLP);

  // loop over the components: verify mult 1 & sufficnet vanishing of proper polynomials 
  for (j = 0; j < numComponents; j++)
  {
    if (pathMod > 0 && !(j % pathMod))
      printf("Verifying component %d of %d\n", j, numComponents);

    num_points = 0;
    for (i = 0; i < witnessSets[j].codim[codim_index[j]].num_set; i++)
      if (component_number[j] == -2 || witnessSets[j].codim[codim_index[j]].component_nums[i] == component_number[j])
      { // verify mult = 1 & deflations needed = 0
        if (witnessSets[j].codim[codim_index[j]].multiplicities[i] != 1 || witnessSets[j].codim[codim_index[j]].deflations_needed[i] != 0)
        { // error
          printf("\n\nERROR: The implemenation of regeneration extension assumes all components are generically reduced!\n");
          bexit(ERROR_INPUT_SYSTEM);
        }

        // verify the poit sufficiently satisfies the system
        regenExtendVerification_point(j, codim_index[j], i, numComponents, witnessSets, T, SLP);

        // show the existence of points
        num_points++;
      }

    if (num_points == 0)
    {
      printf("\n\nERROR: Invalid component number (%d) for system %d.\n", component_number[j], j+1);
      bexit(ERROR_CONFIGURATION);
    }
  }

  return;
}

void setupDegrees_extend(int **orig_degrees, int **new_degrees, int **perm, int currCodim, int top_funcs, int bottom_funcs, int num_var_gps, char *degreeFile)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup degrees for extending regeneration               *
\***************************************************************/
{
  int i, j, max, max_loc, *tempInts = NULL, num_funcs = top_funcs + bottom_funcs, currLoc = 0;
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
    fscanf(tempFile, "%d", &(*orig_degrees)[i]);
  
  // close the file containing the degrees
  fclose(tempFile);

  // allocate perm
  *perm = (int *)bmalloc(num_funcs * sizeof(int));
  for (i = 0; i < num_funcs; i++)
    (*perm)[i] = i;

  // order the top and bottom functions independently in order
  tempInts = (int *)bmalloc(num_funcs * sizeof(int));
  for (i = 0; i < num_funcs; i++)
    tempInts[i] = 1; // will be = 0 when the function has been used

  // order up to currCodim using top functions
  currLoc = 0;
  for (i = 0; i < currCodim; i++)
  { // find the largest degreee still available
    max = max_loc = -1;
    for (j = 0; j < top_funcs; j++)
    { // see if this function is still available and larger than the current max
      if (tempInts[j] == 1 && max < (*orig_degrees)[j])
      { // update max & max_loc
        max = (*orig_degrees)[j];
        max_loc = j;
      }
    }
    (*perm)[currLoc] = max_loc;
    currLoc++;
    tempInts[max_loc] = 0;
  }

  // order the bottom functions based on degree
  for (i = top_funcs; i < num_funcs; i++)
  { // find the largest degreee still available
    max = max_loc = -1;
    for (j = top_funcs; j < num_funcs; j++)
    { // see if this function is still available and larger than the current max
      if (tempInts[j] == 1 && max < (*orig_degrees)[j])
      { // update max & max_loc
        max = (*orig_degrees)[j];
        max_loc = j;
      }
    }
    (*perm)[currLoc] = max_loc;
    currLoc++;
    tempInts[max_loc] = 0;
  }

  // order the remaining functions 
  for (i = currCodim; i < top_funcs; i++)
  { // find the largest degreee still available
    max = max_loc = -1;
    for (j = 0; j < top_funcs; j++)
    { // see if this function is still available and larger than the current max
      if (tempInts[j] == 1 && max < (*orig_degrees)[j])
      { // update max & max_loc
        max = (*orig_degrees)[j];
        max_loc = j;
      }
    }
    (*perm)[currLoc] = max_loc;
    currLoc++;
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

void regenExtendSetup_system_coeff(point_d extra_d, point_mp extra_mp, int numComponents, int *codim_index, witness_t *subWitnessSets, tracker_config_t *T, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup C and coeff for system -- use original slices    *
\***************************************************************/
{ // ASSUME 1 AFFINE VARIABLE GROUP
  int i, j, k, currSlice, currVar, numSlices = 0, origVars = 0, newVars = 0, currCodim = 0;

  for (i = 0; i < numComponents; i++)
  {
    origVars += subWitnessSets[i].orig_variables - 1; // -1 due to homogenizing variable that was added
    currCodim += subWitnessSets[i].codim[codim_index[i]].codim;
  }
  newVars = (RPD->orig_variables - 1) - origVars;

  // setup the system of original slices
  if (T->MPType == 0)
  { // use _d
    mat_d Slices_d, rand_d;
    init_mat_d(rand_d, 0, 0);

    for (i = 0; i < numComponents; i++)
      numSlices += subWitnessSets[i].codim[codim_index[i]].B_d->rows;
    numSlices += newVars;

    // initialize to 0
    init_mat_d(Slices_d, numSlices, RPD->orig_variables);
    Slices_d->rows = numSlices;
    Slices_d->cols = RPD->orig_variables;
    for (j = 0; j < numSlices; j++)
      for (k = 0; k < RPD->orig_variables; k++)
      {
        set_zero_d(&Slices_d->entry[j][k]);
      }

    // setup Slices
    currSlice = 0;
    currVar = 1; // skip homogenizing
    for (i = 0; i < numComponents; i++)
    {
      for (j = subWitnessSets[i].codim[codim_index[i]].B_d->rows - 1; j >= 0; j--)
      { // set homogenizing variable
        set_d(&Slices_d->entry[currSlice][0], &subWitnessSets[i].codim[codim_index[i]].B_d->entry[j][0]);
        // set others
        for (k = subWitnessSets[i].codim[codim_index[i]].B_d->cols - 1; k > 0; k--)
        {
          set_d(&Slices_d->entry[currSlice][currVar + k - 1], &subWitnessSets[i].codim[codim_index[i]].B_d->entry[j][k]);
        }
        currSlice++;
      }
      currVar += subWitnessSets[i].codim[codim_index[i]].B_d->cols - 1;
    }
    for (i = 0; i < newVars; i++)
    {
      neg_d(&Slices_d->entry[currSlice][0], &extra_d->coord[i]);
      set_one_d(&Slices_d->entry[currSlice][currVar + i]);
      currSlice++;
    }

    // randomize the slices
    make_matrix_random_d(rand_d, numSlices, numSlices);
    mat_mul_d(Slices_d, rand_d, Slices_d);

    // setup C, if needed
    if (RPD->orig_variables != RPD->new_variables)
    { // only setup C_d -- go intrinsic on 'extra' slices
      mat_d Q, R, P, tempMat_d;
      init_mat_d(Q, 0, 0);
      init_mat_d(R, 0, 0);
      init_mat_d(P, 0, 0);
      init_mat_d(tempMat_d, RPD->orig_variables - RPD->new_variables, RPD->orig_variables); 
      tempMat_d->rows = RPD->orig_variables - RPD->new_variables;
      tempMat_d->cols = RPD->orig_variables;

      // compute basis for null space via QR decomposition
      k = numSlices - tempMat_d->rows;
      for (i = 0; i < tempMat_d->rows; i++)
        for (j = 0; j < tempMat_d->cols; j++)
        {
          set_d(&tempMat_d->entry[i][j], &Slices_d->entry[i + k][j]);
        }
      transpose_d(tempMat_d, tempMat_d);
      QR_d_prec(Q, R, P, tempMat_d, 0);

      // setup C_d -- last columns of Q
      init_mat_d(RPD->C_d, RPD->orig_variables, RPD->new_variables);
      RPD->C_d->rows = RPD->orig_variables;
      RPD->C_d->cols = RPD->new_variables;
      k = Q->cols - 1;
      for (i = 0; i < RPD->orig_variables; i++)
        for (j = 0; j < RPD->new_variables; j++)
        {
          set_d(&RPD->C_d->entry[i][j], &Q->entry[i][k - j]);
        } 

      clear_mat_d(Q);
      clear_mat_d(R);
      clear_mat_d(P);
      clear_mat_d(tempMat_d);
    } 
    else
    { // set to 0
      RPD->C_d->rows = RPD->C_d->cols = RPD->C_mp->rows = RPD->C_mp->cols = 0; 
    }

    // setup coeff -- start off using random and then adjust based on Slices
    RPD->coeff_mp = NULL;
    RPD->coeff_rat = NULL;

    // allocate coeff_d make each one random
    RPD->coeff_d = (comp_d ***)bmalloc(RPD->num_codim * sizeof(comp_d **));
    for (i = 0; i < RPD->num_codim; i++)
    { // allocate for degree
      RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          get_comp_rand_d(RPD->coeff_d[i][j][k]);
        }
      }
    }

    // setup Slices in new variables
    if (RPD->orig_variables != RPD->new_variables)
    { // remove the 'extra' slices
      Slices_d->rows = Slices_d->rows - (RPD->orig_variables - RPD->new_variables);
      mat_mul_d(Slices_d, Slices_d, RPD->C_d);
    }

    // copy over to coeff in degree 0
    for (i = 0; i < RPD->num_codim - currCodim; i++)
      for (k = 0; k < RPD->new_variables; k++)
      {
        set_d(RPD->coeff_d[i + currCodim][0][k], &Slices_d->entry[i][k]);
      }

    clear_mat_d(Slices_d);
    clear_mat_d(rand_d);
  }
  else if (T->MPType == 1)
  { // use _mp
    mat_mp Slices_mp, rand_mp;
    init_mat_mp(rand_mp, 0, 0);

    for (i = 0; i < numComponents; i++)
      numSlices += subWitnessSets[i].codim[codim_index[i]].B_mp->rows;
    numSlices += newVars;

    // initialize to 0
    init_mat_mp(Slices_mp, numSlices, RPD->orig_variables);
    Slices_mp->rows = numSlices;
    Slices_mp->cols = RPD->orig_variables;
    for (j = 0; j < numSlices; j++)
      for (k = 0; k < RPD->orig_variables; k++)
      {
        set_zero_mp(&Slices_mp->entry[j][k]);
      }

    // setup Slices
    currSlice = 0;
    currVar = 1; // skip homogenizing
    for (i = 0; i < numComponents; i++)
    {
      for (j = subWitnessSets[i].codim[codim_index[i]].B_mp->rows - 1; j >= 0; j--)
      { // set homogenizing variable
        set_mp(&Slices_mp->entry[currSlice][0], &subWitnessSets[i].codim[codim_index[i]].B_mp->entry[j][0]);
        // set others
        for (k = subWitnessSets[i].codim[codim_index[i]].B_mp->cols - 1; k > 0; k--)
        {
          set_mp(&Slices_mp->entry[currSlice][currVar + k - 1], &subWitnessSets[i].codim[codim_index[i]].B_mp->entry[j][k]);
        }
        currSlice++;
      }
      currVar += subWitnessSets[i].codim[codim_index[i]].B_mp->cols - 1;
    }
    for (i = 0; i < newVars; i++)
    {
      neg_mp(&Slices_mp->entry[currSlice][0], &extra_mp->coord[i]);
      set_one_mp(&Slices_mp->entry[currSlice][currVar + i]);
      currSlice++;
    }

    // randomize the slices
    make_matrix_random_mp(rand_mp, numSlices, numSlices, T->Precision);
    mat_mul_mp(Slices_mp, rand_mp, Slices_mp);

    // setup C, if needed
    if (RPD->orig_variables != RPD->new_variables)
    { // only setup C_mp -- go intrinsic on 'extra' slices
      mat_mp Q, R, P, tempMat_mp;
      init_mat_mp(Q, 0, 0);
      init_mat_mp(R, 0, 0);
      init_mat_mp(P, 0, 0);
      init_mat_mp(tempMat_mp, RPD->orig_variables - RPD->new_variables, RPD->orig_variables);
      tempMat_mp->rows = RPD->orig_variables - RPD->new_variables;
      tempMat_mp->cols = RPD->orig_variables;

      // compute basis for null space extra slices via QR decomposition
      k = numSlices - tempMat_mp->rows;
      for (i = 0; i < tempMat_mp->rows; i++)
        for (j = 0; j < tempMat_mp->cols; j++)
        {
          set_mp(&tempMat_mp->entry[i][j], &Slices_mp->entry[i + k][j]);
        }
      transpose_mp(tempMat_mp, tempMat_mp);
      QR_mp_prec(Q, R, P, tempMat_mp, 0, T->Precision);

      // setup C_mp -- last columns of Q
      init_mat_mp(RPD->C_mp, RPD->orig_variables, RPD->new_variables);
      RPD->C_mp->rows = RPD->orig_variables;
      RPD->C_mp->cols = RPD->new_variables;
      k = Q->cols - 1;
      for (i = 0; i < RPD->orig_variables; i++)
        for (j = 0; j < RPD->new_variables; j++)
        {
          set_mp(&RPD->C_mp->entry[i][j], &Q->entry[i][k - j]);
        }

      clear_mat_mp(Q);
      clear_mat_mp(R);
      clear_mat_mp(P);
      clear_mat_mp(tempMat_mp);
    }
    else
    { // set to 0
      RPD->C_d->rows = RPD->C_d->cols = RPD->C_mp->rows = RPD->C_mp->cols = 0;
    }

    // setup coeff -- start off using random and then adjust based on Slices
    RPD->coeff_d = NULL;
    RPD->coeff_rat = NULL;

    // allocate coeff_mp make each one random
    RPD->coeff_mp = (comp_mp ***)bmalloc(RPD->num_codim * sizeof(comp_mp **));
    for (i = 0; i < RPD->num_codim; i++)
    { // allocate for degree
      RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          init_mp(RPD->coeff_mp[i][j][k]);
          get_comp_rand_mp(RPD->coeff_mp[i][j][k]);
        }
      }
    }

    // setup Slices in new variables
    if (RPD->orig_variables != RPD->new_variables)
    { // remove the 'extra' slices
      Slices_mp->rows = Slices_mp->rows - (RPD->orig_variables - RPD->new_variables);
      mat_mul_mp(Slices_mp, Slices_mp, RPD->C_mp);
    }

    // copy over to coeff in degree 0
    for (i = 0; i < RPD->num_codim - currCodim; i++)
      for (k = 0; k < RPD->new_variables; k++)
      {
        set_mp(RPD->coeff_mp[i + currCodim][0][k], &Slices_mp->entry[i][k]);
      }

    clear_mat_mp(Slices_mp);
    clear_mat_mp(rand_mp);
  }
  else
  { // _d, _mp, and _rat
    mat_mp Slices_mp, rand_mp;
    init_mat_mp2(rand_mp, 0, 0, T->AMP_max_prec);

    for (i = 0; i < numComponents; i++)
      numSlices += subWitnessSets[i].codim[codim_index[i]].B_mp->rows;
    numSlices += newVars;

    // initialize to 0
    init_mat_mp2(Slices_mp, numSlices, RPD->orig_variables, T->AMP_max_prec);
    Slices_mp->rows = numSlices;
    Slices_mp->cols = RPD->orig_variables;
    for (j = 0; j < numSlices; j++)
      for (k = 0; k < RPD->orig_variables; k++)
      {
        set_zero_mp(&Slices_mp->entry[j][k]);
      }

    // setup Slices
    currSlice = 0;
    currVar = 1; // skip homogenizing
    for (i = 0; i < numComponents; i++)
    {
      for (j = subWitnessSets[i].codim[codim_index[i]].B_mp->rows - 1; j >= 0; j--)
      { // set homogenizing variable
        rat_to_mp(&Slices_mp->entry[currSlice][0], subWitnessSets[i].codim[codim_index[i]].B_rat[j][0]);
        // set others
        for (k = subWitnessSets[i].codim[codim_index[i]].B_mp->cols - 1; k > 0; k--)
        {
          rat_to_mp(&Slices_mp->entry[currSlice][currVar + k - 1], subWitnessSets[i].codim[codim_index[i]].B_rat[j][k]);
        }
        currSlice++;
      }
      currVar += subWitnessSets[i].codim[codim_index[i]].B_mp->cols - 1;
    }
    for (i = 0; i < newVars; i++)
    {
      neg_mp(&Slices_mp->entry[currSlice][0], &extra_mp->coord[i]);
      set_one_mp(&Slices_mp->entry[currSlice][currVar + i]);
      currSlice++;
    }

    // randomize the slices
    make_matrix_random_mp(rand_mp, numSlices, numSlices, T->AMP_max_prec);
    mat_mul_mp(Slices_mp, rand_mp, Slices_mp);

    // setup C, if needed
    if (RPD->orig_variables != RPD->new_variables)
    { // only setup C_mp -- go intrinsic on 'extra' slices
      mat_mp Q, R, P, tempMat_mp;

      initMP(T->AMP_max_prec);

      init_mat_mp(Q, 0, 0);
      init_mat_mp(R, 0, 0);
      init_mat_mp(P, 0, 0);
      init_mat_mp(tempMat_mp, RPD->orig_variables - RPD->new_variables, RPD->orig_variables);
      tempMat_mp->rows = RPD->orig_variables - RPD->new_variables;
      tempMat_mp->cols = RPD->orig_variables;

      // compute basis for null space extra slices via QR decomposition
      k = numSlices - tempMat_mp->rows;
      for (i = 0; i < tempMat_mp->rows; i++)
        for (j = 0; j < tempMat_mp->cols; j++)
        {
          set_mp(&tempMat_mp->entry[i][j], &Slices_mp->entry[i + k][j]);
        }
      transpose_mp(tempMat_mp, tempMat_mp);
      QR_mp_prec(Q, R, P, tempMat_mp, 0, T->AMP_max_prec);

      // setup C -- last columns of Q
      init_mat_d(RPD->C_d, RPD->orig_variables, RPD->new_variables);
      init_mat_mp2(RPD->C_mp, RPD->orig_variables, RPD->new_variables, RPD->curr_precision);
      init_mat_rat(RPD->C_rat, RPD->orig_variables, RPD->new_variables);
      RPD->C_d->rows = RPD->C_mp->rows = RPD->orig_variables;
      RPD->C_d->cols = RPD->C_mp->cols = RPD->new_variables;
      k = Q->cols - 1;
      for (i = 0; i < RPD->orig_variables; i++)
        for (j = 0; j < RPD->new_variables; j++)
        {
          mp_to_rat(RPD->C_rat[i][j], &Q->entry[i][k - j]);
          set_mp(&RPD->C_mp->entry[i][j], &Q->entry[i][k - j]);
          mp_to_d(&RPD->C_d->entry[i][j], &RPD->C_mp->entry[i][j]);
        }

      initMP(T->Precision);

      clear_mat_mp(Q);
      clear_mat_mp(R);
      clear_mat_mp(P);
      clear_mat_mp(tempMat_mp);
    }
    else
    { // set to 0
      RPD->C_d->rows = RPD->C_d->cols = RPD->C_mp->rows = RPD->C_mp->cols = 0;
    }

    // setup coeff -- start off using random and then adjust based on Slices

    // allocate coeff_mp make each one random
    RPD->coeff_d = (comp_d ***)bmalloc(RPD->num_codim * sizeof(comp_d **));
    RPD->coeff_mp = (comp_mp ***)bmalloc(RPD->num_codim * sizeof(comp_mp **));
    RPD->coeff_rat = (mpq_t ****)bmalloc(RPD->num_codim * sizeof(mpq_t ***));
    for (i = 0; i < RPD->num_codim; i++)
    { // allocate for degree
      RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
      RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
      RPD->coeff_rat[i] = (mpq_t ***)bmalloc(RPD->new_degrees[i] * sizeof(mpq_t **));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
        RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
        RPD->coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          RPD->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
          init_rat(RPD->coeff_rat[i][j][k]);
          get_comp_rand_rat(RPD->coeff_d[i][j][k], RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k], RPD->curr_precision, T->AMP_max_prec, 0, 0);
        }
      }
    }

    // setup Slices in new variables
    if (RPD->orig_variables != RPD->new_variables)
    { // remove the 'extra' slices
      Slices_mp->rows = Slices_mp->rows - (RPD->orig_variables - RPD->new_variables);
      mat_mul_mp(Slices_mp, Slices_mp, RPD->C_mp);
    }

    // copy over to coeff in degree 0
    for (i = 0; i < RPD->num_codim - currCodim; i++)
      for (k = 0; k < RPD->new_variables; k++)
      {
        mp_to_rat(RPD->coeff_rat[i + currCodim][0][k], &Slices_mp->entry[i][k]);
        set_mp(RPD->coeff_mp[i + currCodim][0][k], &Slices_mp->entry[i][k]);
        mp_to_d(RPD->coeff_d[i + currCodim][0][k], RPD->coeff_mp[i + currCodim][0][k]);
      }

    clear_mat_mp(Slices_mp);
    clear_mat_mp(rand_mp);
  }

  return;
}

void regenExtend_reclassify(tracker_config_t *T, regen_pos_dim_t *RPD, int maxCodim, int specificCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reclassify the points in the witness superset using    *
* a fully randomized system                                     *
\***************************************************************/
{
  int i, j, k, *coeff_degrees = RPD->new_degrees;

  // print message about what is happening
  printf("\nReclassifying points.\n");

  // reconstruct the randomization used in RPD --> change sorting by degrees etc.
  free(RPD->orig_degrees);
  free(RPD->P);
  RPD->new_degrees = NULL;
  setupDegrees_orig_new_perm(&RPD->orig_degrees, &RPD->new_degrees, &RPD->P, RPD->num_funcs, RPD->PPD.num_var_gp + RPD->PPD.num_hom_var_gp, "deg.out");

  // recompute W
  for (i = 0; i < RPD->num_codim; i++)
  { // setup for codimension i + 1
    for (j = 0; j <= i; j++)
    {
      for (k = 0; k < RPD->num_funcs - j - 1; k++)
      {
        RPD->W[i][j][k] = RPD->new_degrees[j] - RPD->new_degrees[k + j + 1];
      }
    }
  }

  // randomly generate A and adjust coeff
  if (T->MPType == 0)
  {  // setup A_d
    RPD->sameA = 1;
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // reconstruct A_d[i]
      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_d(RPD->A_d[i], i + 1, RPD->num_funcs);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_d(&RPD->A_d[i]->entry[j][k]);
            }
            else
            { // set to 0
              set_zero_d(&RPD->A_d[i]->entry[j][k]);
            }
      }
      else
      { // copy the top of A_d[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
          }
      }
    }

    // save the linear slices
    mat_d coeff_d;
    init_mat_d(coeff_d, RPD->num_codim, RPD->new_variables);
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_d(&coeff_d->entry[i][j], RPD->coeff_d[i][0][j]);
      }

    // reallocate coeff_d and setup randomly 
    for (i = 0; i < RPD->num_codim; i++)
    { // clear based on old degree
      for (j = 0; j < coeff_degrees[i]; j++)
      { // clear 
        free(RPD->coeff_d[i][j]);
      }
    
      // reallocate for new degrees 
      RPD->coeff_d[i] = (comp_d **)brealloc(RPD->coeff_d[i], RPD->new_degrees[i] * sizeof(comp_d *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          get_comp_rand_d(RPD->coeff_d[i][j][k]);
        }
      }
    }

    // copy coeff_d
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_d(RPD->coeff_d[i][0][j], &coeff_d->entry[i][j]);
      }

    clear_mat_d(coeff_d);
  }
  else if (T->MPType == 1)
  { // setup A_mp
    RPD->sameA = 1;
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // reconstruct A_mp[i]
      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_mp(RPD->A_mp[i], i + 1, RPD->num_funcs, T->Precision);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_mp(&RPD->A_mp[i]->entry[j][k]);
            }
            else
            { // set to 0
              set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
            }
      }
      else
      { // copy the top of A_mp[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
          }
      }
    }

    // save the linear slices
    mat_mp coeff_mp;
    init_mat_mp(coeff_mp, RPD->num_codim, RPD->new_variables);
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_mp(&coeff_mp->entry[i][j], RPD->coeff_mp[i][0][j]);
      }

    // reallocate coeff_mp and setup randomly
    for (i = 0; i < RPD->num_codim; i++)
    { // clear based on old degree
      for (j = 0; j < coeff_degrees[i]; j++)
      { // clear
        for (k = 0; k < RPD->new_variables; k++)
        {
          clear_mp(RPD->coeff_mp[i][j][k]);
        }
        free(RPD->coeff_mp[i][j]);
      }

      // reallocate for new degrees
      RPD->coeff_mp[i] = (comp_mp **)brealloc(RPD->coeff_mp[i], RPD->new_degrees[i] * sizeof(comp_mp *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          init_mp(RPD->coeff_mp[i][j][k]);
          get_comp_rand_mp(RPD->coeff_mp[i][j][k]);
        }
      }
    }

    // copy coeff_mp
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_mp(RPD->coeff_mp[i][0][j], &coeff_mp->entry[i][j]);
      }

    clear_mat_mp(coeff_mp);
  }
  else
  { // setup A
    RPD->sameA = 1;
    // make them in reverse order
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { // reconstruct A[i]
      if (i == RPD->num_codim - 1)
      { // we generate the main matrix
        make_matrix_random_rat(RPD->A_d[i], RPD->A_mp[i], RPD->A_rat[i], i + 1, RPD->num_funcs, RPD->curr_precision, T->AMP_max_prec, 0, 0);
        // make upper triangular with 1's on diagonal
        for (j = 0; j <= i; j++)
          for (k = 0; k <= j; k++)
            if (j == k)
            { // set to 1
              set_one_d(&RPD->A_d[i]->entry[j][k]);
              set_one_mp(&RPD->A_mp[i]->entry[j][k]);
              set_one_rat(RPD->A_rat[i][j][k]);
            }
            else
            { // set to 0
              set_zero_d(&RPD->A_d[i]->entry[j][k]);
              set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
              set_zero_rat(RPD->A_rat[i][j][k]);
            }
      }
      else
      { // copy the top of A_d[i+1]
        for (j = 0; j <= i; j++)
          for (k = 0; k < RPD->num_funcs; k++)
          {
            set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
            set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
            set_rat(RPD->A_rat[i][j][k], RPD->A_rat[i+1][j][k]);
          }
      }
    }

    // save the linear slices
    mpq_t ***coeff_rat = NULL;
    init_mat_rat(coeff_rat, RPD->num_codim, RPD->new_variables);
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_rat(coeff_rat[i][j], RPD->coeff_rat[i][0][j]);
      }

    // reallocate coeff_d, coeff_mp, and coeff_rat, and setup randomly
    for (i = 0; i < RPD->num_codim; i++)
    { // clear based on old degree
      for (j = 0; j < coeff_degrees[i]; j++)
      { // clear
        for (k = 0; k < RPD->new_variables; k++)
        {
          clear_mp(RPD->coeff_mp[i][j][k]);
          clear_rat(RPD->coeff_rat[i][j][k]);
          free(RPD->coeff_rat[i][j][k]);
        }
        free(RPD->coeff_d[i][j]);
        free(RPD->coeff_mp[i][j]);
        free(RPD->coeff_rat[i][j]);
      }

      // reallocate for new degrees
      RPD->coeff_d[i] = (comp_d **)brealloc(RPD->coeff_d[i], RPD->new_degrees[i] * sizeof(comp_d *));
      RPD->coeff_mp[i] = (comp_mp **)brealloc(RPD->coeff_mp[i], RPD->new_degrees[i] * sizeof(comp_mp *));
      RPD->coeff_rat[i] = (mpq_t ***)brealloc(RPD->coeff_rat[i], RPD->new_degrees[i] * sizeof(mpq_t **));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      { // allocate for variables
        RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
        RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
        RPD->coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
        for (k = 0; k < RPD->new_variables; k++)
        { // make random
          RPD->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
          init_rat(RPD->coeff_rat[i][j][k]);
          get_comp_rand_rat(RPD->coeff_d[i][j][k], RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k], RPD->curr_precision, T->AMP_max_prec, 0, 0);
        }
      }
    }

    // copy coeff_rat and setup others
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_rat(RPD->coeff_rat[i][0][j], coeff_rat[i][j]);
        rat_to_mp(RPD->coeff_mp[i][0][j], coeff_rat[i][j]);
        mp_to_d(RPD->coeff_d[i][0][j], RPD->coeff_mp[i][0][j]);
      }

    clear_mat_rat(coeff_rat, RPD->num_codim, RPD->new_variables);
  }

  // loop through the codimensions
  for (i = 0; i < maxCodim; i++)
  {
    if (specificCodim <= 0 || RPD->codim[i].codim == specificCodim)
    {
      printf("Reclassifying codimension %d: %d point%s to consider.\n", RPD->codim[i].codim, RPD->codim[i].num_superset, RPD->codim[i].num_superset == 1 ? "" : "s");

      if (RPD->codim[i].num_superset > 0)
      { // reclassify this codimension
        regenExtend_reclassify_codim(T, RPD, i);
      }
    }
  }

  free(coeff_degrees);

  return;
}

void regenExtend_reclassify_codim(tracker_config_t *T, regen_pos_dim_t *RPD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reclassify the points in this witness superset         *
\***************************************************************/
{
  int i, j, rV, pt_prec, last_prec, corank, size, curr_codim = RPD->codim[codim_index].codim;
  int num_nonsing = 0;
  char *str = NULL;
  FILE *CODIM = NULL, *OUT = NULL, *TEMP = fopen("temp_file", "w");
  double CN, smallest_nonzero_SV, largest_zero_SV;
  comp_d final_d;
  comp_mp final_mp;
  point_d pt_d, last_d, int_d, int_last_d;
  point_mp pt_mp, last_mp, int_mp, int_last_mp;
  mat_d tr_d;
  mat_mp tr_mp;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  // initialize
  RPD->curr_codim = codim_index;
  init_d(final_d);
  init_mp(final_mp);
  init_point_d(pt_d, RPD->new_variables);
  init_point_d(last_d, RPD->new_variables);
  init_point_mp(pt_mp, RPD->new_variables);
  init_point_mp(last_mp, RPD->new_variables);
  init_point_d(int_d, 0);
  init_point_d(int_last_d, 0);
  init_point_mp(int_mp, 0);
  init_point_mp(int_last_mp, 0);
  init_mat_d(tr_d, 0, 0);
  init_mat_mp(tr_mp, 0, 0);
  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);
  pt_d->size = last_d->size = pt_mp->size = last_mp->size = RPD->new_variables;
  set_zero_d(final_d);
  set_zero_mp(final_mp);

  size = 1 + snprintf(NULL, 0, "witness_superset_%d", curr_codim);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "witness_superset_%d", curr_codim);
  CODIM = fopen(str, "r");
  if (CODIM == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", str);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  OUT = fopen("witness_superset_test", "w");

  // compute tr_d
  if (RPD->codim[codim_index].useIntrinsicSlice && T->MPType != 1)
  {  
    transpose_d(tr_d, RPD->codim[codim_index].B_d);
  }

  // loop over the points
  for (i = 0; i < RPD->codim[codim_index].num_superset; i++)
  { // read in point
    fscanf(CODIM, "%d%d", &pt_prec, &size);
    if (pt_prec < 64)
    {
      increase_size_point_d(pt_d, RPD->new_variables);
      pt_d->size = RPD->new_variables;
      for (j = 0; j < RPD->new_variables; j++)
        fscanf(CODIM, "%lf%lf", &pt_d->coord[j].r, &pt_d->coord[j].i);

      // setup int_d
      if (RPD->codim[codim_index].useIntrinsicSlice)
      {
        extrinsicToIntrinsic_d(int_d, pt_d, tr_d, RPD->codim[codim_index].p_d);
      }
      else
      {
        point_cp_d(int_d, pt_d);
      }
    }
    else
    {
      setprec_point_mp(pt_mp, pt_prec);
      increase_size_point_mp(pt_mp, RPD->new_variables);
      pt_mp->size = RPD->new_variables;
      for (j = 0; j < RPD->new_variables; j++)
      {
        mpf_inp_str(pt_mp->coord[j].r, CODIM, 10);
        mpf_inp_str(pt_mp->coord[j].i, CODIM, 10);
        scanRestOfLine(CODIM);
      }

      // setup int_mp
      setprec_point_mp(int_mp, pt_prec);
      if (RPD->codim[codim_index].useIntrinsicSlice)
      {
        change_regen_pos_dim_prec(RPD, pt_prec);
        setprec_mat_mp(tr_mp, pt_prec);
        transpose_mp(tr_mp, RPD->codim[codim_index].B_mp);
        extrinsicToIntrinsic_mp(int_mp, pt_mp, tr_mp, RPD->codim[codim_index].p_mp);
      }
      else
      {
        point_cp_mp(int_mp, pt_mp);
      }
    }

    // read in last approximation
    fscanf(CODIM, "%d", &last_prec);
    if (last_prec < 64)
    {
      increase_size_point_d(last_d, RPD->new_variables);
      last_d->size = RPD->new_variables;
      for (j = 0; j < RPD->new_variables; j++)
        fscanf(CODIM, "%lf%lf", &last_d->coord[j].r, &last_d->coord[j].i);

      // setup int_last_d
      if (RPD->codim[codim_index].useIntrinsicSlice)
      {
        extrinsicToIntrinsic_d(int_last_d, last_d, tr_d, RPD->codim[codim_index].p_d);
      }
      else
      {
        point_cp_d(int_last_d, last_d);
      }
    }
    else
    {
      setprec_point_mp(last_mp, last_prec);
      increase_size_point_mp(last_mp, RPD->new_variables);
      last_mp->size = RPD->new_variables;
      for (j = 0; j < RPD->new_variables; j++)
      {
        mpf_inp_str(last_mp->coord[j].r, CODIM, 10);
        mpf_inp_str(last_mp->coord[j].i, CODIM, 10);
        scanRestOfLine(CODIM);
      }

      // setup int_last_mp
      setprec_point_mp(int_last_mp, pt_prec);
      if (RPD->codim[codim_index].useIntrinsicSlice)
      {
        change_regen_pos_dim_prec(RPD, last_prec);
        setprec_mat_mp(tr_mp, last_prec);
        transpose_mp(tr_mp, RPD->codim[codim_index].B_mp);
        extrinsicToIntrinsic_mp(int_last_mp, last_mp, tr_mp, RPD->codim[codim_index].p_mp);
      }
      else
      {
        point_cp_mp(int_last_mp, pt_mp);
      }
    }

    // read in other data 
    if (pt_prec < 64)
    { // read in final_d
      fscanf(CODIM, "%lf%lf", &final_d->r, &final_d->i);
    }
    else
    { // read in final_mp
      setprec_mp(final_mp, pt_prec);
      mpf_inp_str(final_mp->r, CODIM, 10);
      mpf_inp_str(final_mp->i, CODIM, 10);
    }

    // other data
    fscanf(CODIM, "%lf%d%lf%lf%d\n", &CN, &corank, &smallest_nonzero_SV, &largest_zero_SV, &rV);

    // compute the data
    if (T->MPType == 0)
    {
      corank = Cauchy_corank_d(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_d, int_last_d, final_d, T, TEMP, &e_d, RPD, regen_pos_dim_eval_d);
    }
    else if (T->MPType == 1)
    {
      corank = Cauchy_corank_mp(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_mp, int_last_mp, final_mp, T, TEMP, &e_mp, RPD, regen_pos_dim_eval_mp);
    }
    else
    {
      corank = Cauchy_corank_amp(&CN, &smallest_nonzero_SV, &largest_zero_SV, int_d, int_mp, pt_prec, int_last_d, int_last_mp, last_prec, final_d, final_mp, T, TEMP, &e_d, &e_mp, RPD, RPD, regen_pos_dim_eval_d, regen_pos_dim_eval_mp, change_regen_pos_dim_prec);
    }

    // check for nonsingular
    if (corank == 0)
    { // increment
      num_nonsing++;
    }

    // print data to OUT
    fprintf(OUT, "%d\n", pt_prec);
    if (pt_prec < 64)
    { // print pt_d
      fprintf(OUT, "%d\n", pt_d->size);
      for (j = 0; j < pt_d->size; j++)
        fprintf(OUT, "%.15e %.15e\n", pt_d->coord[j].r, pt_d->coord[j].i);
    }
    else
    { // print pt_mp
      fprintf(OUT, "%d\n", pt_mp->size);
      for (j = 0; j < pt_mp->size; j++)
      {
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].r);
        fprintf(OUT, " ");
        mpf_out_str(OUT, 10, 0, pt_mp->coord[j].i);
        fprintf(OUT, "\n");
      }
    }
    // print prec for last_approx
    fprintf(OUT, "%d\n", last_prec);
    if (last_prec < 64)
    { // print last_d
      for (j = 0; j < last_d->size; j++)
        fprintf(OUT, "%.15e %.15e\n", last_d->coord[j].r, last_d->coord[j].i);
    }
    else
    { // print last_mp
      for (j = 0; j < last_mp->size; j++)
      {
        mpf_out_str(OUT, 10, 0, last_mp->coord[j].r);
        fprintf(OUT, " ");
        mpf_out_str(OUT, 10, 0, last_mp->coord[j].i);
        fprintf(OUT, "\n");
      }
    }

    if (pt_prec < 64)
    { // print final_d
      fprintf(OUT, "%.15e %.15e\n", final_d->r, final_d->i);
    }
    else
    { // print final_mp
      mpf_out_str(OUT, 10, 0, final_mp->r);
      fprintf(OUT, " ");
      mpf_out_str(OUT, 10, 0, final_mp->i);
      fprintf(OUT, "\n");
    }

    // print other data
    fprintf(OUT, "%.15e\n%d\n%.15e\n%.15e\n%d\n", CN, corank, smallest_nonzero_SV, largest_zero_SV, rV);
  }

  // reset the numbers
  RPD->codim[codim_index].num_nonsing = num_nonsing;
  RPD->codim[codim_index].num_sing = RPD->codim[codim_index].num_superset - num_nonsing;

  // close and rename
  fclose(OUT);
  fclose(CODIM);
  rename("witness_superset_test", str);
  fclose(TEMP);
  remove("temp_file");

  // clear
  free(str);
  clear_d(final_d);
  clear_mp(final_mp);
  clear_point_d(pt_d);
  clear_point_d(last_d);
  clear_point_mp(pt_mp);
  clear_point_mp(last_mp);
  clear_point_d(int_d);
  clear_point_d(int_last_d);
  clear_point_mp(int_mp);
  clear_point_mp(int_last_mp);
  clear_mat_d(tr_d);
  clear_mat_mp(tr_mp);
  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  return;
}

// ******************************** PARALLEL FUNCTIONS ***************************************

#ifdef _HAVE_MPI

void head_regenExtend_run(double midpoint_tol, int maxCodim, int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, int curr_codim_index, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: run the extension in parallel -- headnode              * 
\***************************************************************/
{
  int codim_index, minPacketSize = 1, maxPacketSize = 20;

  // setup maxCodim
  MPI_Bcast(&maxCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // send 'curr_codim_index
  MPI_Bcast(&curr_codim_index, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (curr_codim_index + 1 < maxCodim)
  { // send T to the workers
    bcast_tracker_config_t(T, my_id, headnode);

    // send RPD to the workers
    bcast_regen_pos_dim_t(RPD, T->MPType, my_id, headnode);

    // send codimension information to the workers
    for (codim_index = 0; codim_index < RPD->num_codim; codim_index++)
    { // send the codimension data to the workers
      bcast_regenCodim_t(&RPD->codim[codim_index], RPD->curr_precision, T->MPType, my_id, headnode);
    }

    // prepare the first codimension
    int size, num_paths, num_crossings, codim;
    char *str = NULL;
    char witName[] = "witness_superset", startName[] = "startRPD", outName[] = "output", failName[] = "fail", rawPrepareFile[] = "rawout_prepare", midName[] = "midpath_data", midPrepareFile[] = "midout_prepare", rawTrackFile[] = "rawout_track", midTrackFile[] = "midout_track", rawSortFile[] = "rawout_sort";
    FILE *START, *RAWOUT, *FAIL, *OUT, *WITSUPER;
    trackingStats trackCount;
    init_trackingStats(&trackCount);

    // open OUT & FAIL
    OUT = fopen(outName, "w");
    FAIL = fopen(failName, "w");
    START = fopen("startRPD_test", "r");

    // open RAWOUT for preparing
    size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, RPD->codim[curr_codim_index + 1].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", rawPrepareFile, RPD->codim[curr_codim_index + 1].codim);
    RAWOUT = fopen(str, "w");

    // setup the name of the file to contain the next set of start points
    size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[curr_codim_index + 1].codim);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", startName, RPD->codim[curr_codim_index + 1].codim);

    // create the next set of start points for the next codim
    num_paths = head_regen_pos_dim_PrepareNextCodim(minPacketSize, maxPacketSize, &trackCount, curr_codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, str, my_id, headnode, num_processes);

    // close files
    fclose(RAWOUT);
    fclose(START);

    // wait until the files are closed
    MPI_Barrier(MPI_COMM_WORLD);

    // check for path crossings for this codimension - do not delete!
    if (RPD->codim[codim_index].useIntrinsicSlice)
      num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, RPD->codim[curr_codim_index].codim + 1, midpoint_tol, my_id, headnode, num_processes);
    else
      num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, RPD->new_variables, midpoint_tol, my_id, headnode, num_processes);

    if (num_crossings > 0)
      printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

    // now that the next codimension is prepared, track through the rest of them!
    for (codim_index = curr_codim_index + 1; codim_index < maxCodim; codim_index++)
    { // track the paths for this codim
      codim = RPD->codim[codim_index].codim;
      num_paths = RPD->codim[codim_index].num_paths;

      // setup START
      size = 1 + snprintf(NULL, 0, "%s_%d", startName, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", startName, codim);
      START = fopen(str, "r");
      // make sure START exits
      if (START == NULL)
      {
        printf("ERROR: The file to contain the start points, '%s', does not exist!\n", str);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // setup RAWOUT
      size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawTrackFile, codim);
      RAWOUT = fopen(str, "w");

      // read in the number of paths from START
      fscanf(START, "%d", &num_crossings);
      // make sure that we have agreement
      if (num_paths != num_crossings)
      {
        printf("ERROR: The number of paths (%d vs %d) described in '%s' is not correct!\n", num_paths, num_crossings, str);
        bexit(ERROR_INVALID_SIZE);
      }

      // track the paths for this codimension
      head_regen_pos_dim_TrackCodim(minPacketSize, maxPacketSize, &trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, my_id, headnode, num_processes);

      // close START & RAWOUT
      fclose(START);
      fclose(RAWOUT);

      // wait until the files are closed
      MPI_Barrier(MPI_COMM_WORLD);

      // check for path crossings for this codimension - do not delete!
      if (RPD->codim[codim_index].useIntrinsicSlice)
        num_crossings = parallel_midpoint_checking(midName, midTrackFile, 0, num_paths, codim, midpoint_tol, my_id, headnode, num_processes);
      else
        num_crossings = parallel_midpoint_checking(midName, midTrackFile, 0, num_paths, RPD->new_variables, midpoint_tol, my_id, headnode, num_processes);

      if (num_crossings > 0)
        printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);

      // open WITSUPER for sorting - witness superset points
      size = 1 + snprintf(NULL, 0, "%s_%d", witName, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", witName, codim);
      WITSUPER = fopen(str, "w");

      // setup START to be the file containing the endpoints
      size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawTrackFile, codim);
      START = fopen(str, "r");
      // make sure START exits
      if (START == NULL)
      {
        printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // reopen RAWOUT for sorting
      size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawSortFile, codim);
      RAWOUT = fopen(str, "w");

      // sort the paths for this codimension
      head_regen_pos_dim_SortCodim(minPacketSize, maxPacketSize, &trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, WITSUPER, my_id, headnode, num_processes);

      // close START, RAWOUT & WITSUPER
      fclose(START);
      fclose(RAWOUT);
      fclose(WITSUPER);

      // print the regen summary
      RAWOUT = fopen("regenSummary", "w");
      printRPDSummaryData(RPD, codim_index, RAWOUT);
      fclose(RAWOUT);

      // see if we need to prepare the next codim
      if (codim_index + 1 < maxCodim)
      { // setup START to be the file containing the ones that need moved forward
        size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
        str = (char *)brealloc(str, size * sizeof(char));
        sprintf(str, "%s_%d", rawSortFile, codim);
        START = fopen(str, "r");
        // make sure START exits
        if (START == NULL)
        {
          printf("ERROR: The file to contain the end points, '%s', does not exist!\n", str);
          bexit(ERROR_FILE_NOT_EXIST);
        }

        // reopen RAWOUT for preparing
        size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, RPD->codim[codim_index + 1].codim);
        str = (char *)brealloc(str, size * sizeof(char));
        sprintf(str, "%s_%d", rawPrepareFile, RPD->codim[codim_index + 1].codim);
        RAWOUT = fopen(str, "w");

        // setup the name of the file to contain the next set of start points
        size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[codim_index + 1].codim);
        str = (char *)brealloc(str, size * sizeof(char));
        sprintf(str, "%s_%d", startName, RPD->codim[codim_index + 1].codim);

        // create the next set of start points for the next codim
        num_paths = head_regen_pos_dim_PrepareNextCodim(minPacketSize, maxPacketSize, &trackCount, codim_index, maxCodim, pathMod, T, RPD, START, OUT, RAWOUT, FAIL, str, my_id, headnode, num_processes);

        // close RAWOUT & START
        fclose(RAWOUT);
        fclose(START);

        // wait until the files are closed
        MPI_Barrier(MPI_COMM_WORLD);

        // check for path crossings for this codimension - do not delete!
        if (RPD->codim[codim_index].useIntrinsicSlice)
          num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, codim + 1, midpoint_tol, my_id, headnode, num_processes);
        else
          num_crossings = parallel_midpoint_checking(midName, midPrepareFile, 0, num_paths, RPD->new_variables, midpoint_tol, my_id, headnode, num_processes);

        if (num_crossings > 0)
          printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame for this level.\n\n", num_crossings);
      }
      else
      { // print the bottom codim file
        size = 1 + snprintf(NULL, 0, "%s_%d", startName, RPD->codim[codim_index].codim + 1);
        str = (char *)brealloc(str, size * sizeof(char));
        sprintf(str, "%s_%d", startName, RPD->codim[codim_index].codim + 1);
        RAWOUT = fopen(str, "w");
        printRPDRelevantData(RPD, T->MPType, codim_index + 1, RAWOUT);
        fclose(RAWOUT);
      }

      // delete RAWOUT files
      size = 1 + snprintf(NULL, 0, "%s_%d", rawTrackFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawTrackFile, codim);
      remove(str);

      size = 1 + snprintf(NULL, 0, "%s_%d", rawSortFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawSortFile, codim);
      remove(str);

      size = 1 + snprintf(NULL, 0, "%s_%d", rawPrepareFile, codim);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", rawPrepareFile, codim);
      remove(str);
    }

    // wait until all workers have closed the files
    MPI_Barrier(MPI_COMM_WORLD);

    // clear memory
    free(str);
  }

  return;
}

void worker_regenExtend(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: run the extension in parallel -- worker                * 
\***************************************************************/
{
  int size, codim_index, curr_codim_index, maxCodim;
  worker_info recvType;
  char *str = NULL, midPrepareFile[] = "midout_prepare", midTrackFile[] = "midout_track";
  FILE *OUT = NULL, *MIDOUT = NULL, *RAWOUT = NULL, *FAIL = NULL;

  // setup maxCodim
  MPI_Bcast(&maxCodim, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // send 'curr_codim_index
  MPI_Bcast(&curr_codim_index, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (curr_codim_index + 1 < maxCodim)
  { // need to run
    trackingStats trackCount;
    tracker_config_t T;
    regen_pos_dim_t RPD;

    // recv T
    bcast_tracker_config_t(&T, my_id, headnode);

    initMP(T.Precision);

    // recv RPD 
    bcast_regen_pos_dim_t(&RPD, T.MPType, my_id, headnode);

    // send codimension information to the workers
    RPD.codim = (regenCodim_t *)bmalloc(RPD.num_codim * sizeof(regenCodim_t));
    for (codim_index = 0; codim_index < RPD.num_codim; codim_index++)
    { // send the codimension data to the workers
      bcast_regenCodim_t(&RPD.codim[codim_index], RPD.curr_precision, T.MPType, my_id, headnode);
    }

    // setup the local files - OUT & FAIL & RAWOUT
    size = 1 + snprintf(NULL, 0, "output_%d", my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "output_%d", my_id);
    OUT = fopen(str, "w");

    size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "fail_%d", my_id);
    FAIL = fopen(str, "w");

    size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "rawout_%d", my_id);
    RAWOUT = fopen(str, "w+");

    // open MIDOUT for preparing next level
    size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", midPrepareFile, my_id);
    MIDOUT = fopen(str, "w");

    // prepare for the next codim
    worker_regen_pos_dim_PrepareNextCodim(&trackCount, curr_codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

    // close
    fclose(MIDOUT);

    // wait until everybody has closed MIDOUT
    MPI_Barrier(MPI_COMM_WORLD);

    // consider doing midpoint checking in parallel

    // now that the next codimension is prepared, track through the rest of them!
    for (codim_index = curr_codim_index + 1; codim_index < maxCodim; codim_index++)
    { // setup MIDOUT for this codim
      size = 1 + snprintf(NULL, 0, "%s_%d", midTrackFile, my_id);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", midTrackFile, my_id);
      MIDOUT = fopen(str, "w");

      // track this codim
      worker_regen_pos_dim_TrackCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

      // close MIDOUT
      fclose(MIDOUT);

      // wait until everybody has closed MIDOUT
      MPI_Barrier(MPI_COMM_WORLD);

      // consider doing midpoint checking in parallel

      // setup MIDOUT for sorting
      MIDOUT = NULL;

      // sort this codim
      worker_regen_pos_dim_SortCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

      // see if we need to prepare the next codim
      if (codim_index + 1 < maxCodim)
      { // reopen MIDOUT for preparing next level
        size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
        str = (char *)brealloc(str, size * sizeof(char));
        sprintf(str, "%s_%d", midPrepareFile, my_id);
        MIDOUT = fopen(str, "w");

        // prepare for the next codim
        worker_regen_pos_dim_PrepareNextCodim(&trackCount, codim_index, &T, &RPD, OUT, RAWOUT, MIDOUT, FAIL, my_id, headnode, num_processes);

        // close MIDOUT
        fclose(MIDOUT);

        // wait until everybody has closed MIDOUT
        MPI_Barrier(MPI_COMM_WORLD);

        // consider doing midpoint checking in parallel
      }
    }

    // close OUT, RAWOUT & FAIL
    fclose(OUT);
    fclose(RAWOUT);
    fclose(FAIL);

    // delete RAWOUT
    size = 1 + snprintf(NULL, 0, "rawout_%d", my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "rawout_%d", my_id);
    remove(str);

    // delete MIDOUT
    size = 1 + snprintf(NULL, 0, "%s_%d", midTrackFile, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", midTrackFile, my_id);
    remove(str);

    size = 1 + snprintf(NULL, 0, "%s_%d", midPrepareFile, my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", midPrepareFile, my_id);
    remove(str);

    // delete FAIL
    size = 1 + snprintf(NULL, 0, "fail_%d", my_id);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "fail_%d", my_id);
    remove(str);

    // clear memory
    regen_pos_dim_clear(&RPD, T.MPType);
    tracker_config_clear(&T);
    clearMP();
    free(str);

    // wait until all workers have closed the files
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // determine if decomposition is needed
  bcast_worker_info(&recvType, my_id, headnode);

  if (recvType.dataType != STOPCODE)
  { // now do the junk removal & break up into irreducible components
    worker_witness_superset_decomposition(my_id, num_processes, headnode);
  }

  return;
}

#endif
 
