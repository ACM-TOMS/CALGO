// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "ppParse.h"
#include "parallel.h"
#include <limits.h>

// these functions are local to parse_input.c
void addBounds(double degBound, double coeffBound, int printDeg, char *fileName);
int findConfig(int *found_loc, char *ch, int *count, FILE *IN, char **str, int *len, int size);
void errorConfig(char *ch, int count, FILE *IN);
void configErrorChecking(tracker_config_t *T, int *maxCodim, int *specificCodim, int trackType);
void printConfigLine(FILE *OUT, int valueType, void const *value, int *title, char *str);
void loadDefaultConfig(tracker_config_t *T, int MPType, int *genType, int *userHom, int *useRegen, int *regenStartLevel, int *maxCodim, int *specificCodim, int *printMod, double *intrinsicMult, int *reducedOnly, int *constructWitnessSet, int *supersetOnly, int *paramHom);
void parseConfigurations(FILE *IN, int *trackType, int *MPType, int *genType, int *userHom, unsigned int *randomSeed, int *sharpenOnly, int *needToDiff, int *remove_temp, int *paramHom, int *maxPrec, int *useRegen, int *regenStartLevel);


void parse_input(char *inputName, int *trackType, int *MPType, int *genType, int *userHom, unsigned int *randomSeed, int *sharpenOnly, int *needToDiff, int *remove_temp, int useParallelDiff, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int paramHom = 0, maxPrec = 0, useRegen = 0, regenStartLevel = 0, setupRandom = 1, numRandom = 0;
  FILE *IN = fopen(inputName, "r"), *OUT = NULL;
  if (IN == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", inputName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // initialize to defaults
  *trackType = 0;  // zero dimension tracking
  *MPType = 2;     // use AMP
  *genType = 2;    // use regen cascade for generating witness supersets
  *randomSeed = 0; // use random seed
  *userHom = 0;    // add our own start system
  *sharpenOnly = 0;// do not only use sharpening
  *needToDiff = 1; // need to compute the derivative of the system
  *remove_temp = 1;// remove the temporary files

  // split 'IN' into CON & INPUT
  splitParse(IN, "func_input_rand", "config");

  // close the input file
  fclose(IN);

  // detect certain config values
  IN = fopen("config", "r");
  parseConfigurations(IN, trackType, MPType, genType, userHom, randomSeed, sharpenOnly, needToDiff, remove_temp, &paramHom, &maxPrec, &useRegen, &regenStartLevel);
  fclose(IN);

  // generate a random seed if it is needed
  if (*randomSeed == 0)
    *randomSeed = (unsigned int) time(NULL);

  // seed the random number generator
  srand(*randomSeed);

  // setup random values
  IN = fopen("func_input_rand", "r");
  OUT = fopen("func_input", "w");
  setupRandom = !(*sharpenOnly || ((*genType == 2 || useRegen == 1) && regenStartLevel > 0) || (2 <= *trackType && *trackType <= 6)); // not to include regeneration extension
  numRandom = setupRandomValues(OUT, IN, setupRandom, maxPrec);
  fclose(IN);
  fclose(OUT);
  remove("func_input_rand");

  // check number of random
  if (numRandom > 0 && *trackType == 7)
  { // random not allowed in regeneration extension!
    printf("\nERROR: Regeneration extension assumes that 'random' and 'random_real' are not used!\n       Please replace with 'constant' and define with actual values.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // verify paramHom
  if (*trackType > 0)
  { // set paramHom to 0
    paramHom = 0;
  }
  else if (paramHom && *userHom > 0)
  { // error!
    printf("\nERROR: Please choose either a parameter homotopy or a user-defined homotopy!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (paramHom > 0 && *needToDiff == 0)
  { // set needToDiff to 2 if it is 0 (update only the numbers)
    *needToDiff = 2;
  }
  else if (paramHom == 2 && *sharpenOnly == 0)
  { // display message
    printf("\nNOTE: Please make sure that you have also setup the proper start points file for your parameter homotopy\n");
    printf("      with the starting parameters in 'start_parameters' and the target parameters in 'final_parameters'.\n\n");
  }

  if (*userHom > 0)
  { // the user defined a homotopy
    if (*sharpenOnly == 0)
      printf("\nNOTE: Please make sure that you have also setup the proper start points file for your homotopy.\n\n");
  }

#ifdef _HAVE_MPI
  // send data to workers
  useParallelDiff = useParallelDiff && *needToDiff;
  MPI_Bcast(&useParallelDiff, 1, MPI_INT, headnode, MPI_COMM_WORLD);
#endif

  if (*needToDiff == 0)
  { // display message
    printf("\nNOTE: You have requested that Bertini not differentiate the system again!\n");
  }
  else if (*needToDiff == 2)
  { // display message
    printf("\nNOTE: You have requested that Bertini only update the numbers used in the system.\n");
    printf("      This option may yield an invalid system!\n");
  }

  if (*needToDiff)
  { // setup the SLP
    double degBound, coeffBound;
    preProcessArray ppArray; 
    parseArray parArray;
    variablegroupArray vargpArray;
    subFuncData *sfData = NULL;

    // preprocess func_input
    IN = fopen("func_input", "r");
    preProcessParse(&ppArray, IN, paramHom, 1); 
    rewind(IN);

    // add the subfunction number information
    setup_subFuncData_numbers(&ppArray);

    // seutp variable group information
    setup_variablegroupArray(&ppArray, &vargpArray, *userHom == 1);

    if (*userHom > 0)
    { // verify information for user defined homotopy
      verify_userdefined_homotopy(&ppArray, *userHom);
    }
    else if (paramHom)
    { // verify information for a parameter homotopy
      verify_parameter_homotopy(&ppArray);
    }
    else
    { // verify information for standard tracking
      verify_standard_homotopy(&ppArray, &vargpArray, *trackType > 0);
    }

    if (*needToDiff == 1)
    { // print preproc_data
      OUT = fopen("preproc_data", "w");
      print_preproc_data(OUT, ppArray.types[FUNCTIONTYPE].numType, vargpArray.numVarGps, vargpArray.numHomGps, vargpArray.types, vargpArray.sizes);
      fclose(OUT);
      OUT = NULL;
    }

    // setup system if using a parameter homotopy -- reparses system 
    if (paramHom)
    { // setup new file
      OUT = fopen("func_input_param", "w+");
      setupParameterHomotopy(OUT, IN, paramHom, *sharpenOnly || (useRegen == 1 && regenStartLevel > 0), maxPrec, &ppArray, &vargpArray);
      // close IN & OUT
      fclose(IN);
      fclose(OUT);
      // setup IN
      IN = fopen("func_input_param", "r");
    }

    if (*needToDiff == 1)
    { // process the system -- always add all of the numbers
      processSystem(&parArray, &sfData, &ppArray, &vargpArray, *userHom == 1, 1, IN);
      fclose(IN);

      // compute the SLP
      compute_SLP(&degBound, &coeffBound, !(*userHom == 1), &parArray, sfData, &vargpArray, *userHom == 1, *trackType > 0, useParallelDiff, my_id, num_processes, headnode);

      // add degBond & coeffBound to config (put first so user can override them)
      if (!(*userHom == 1))
        addBounds(degBound, coeffBound, 1, "config");
    }
    else
    { // print the new num.out file
      setupSLPNumOut(ppArray.types[NUMBERTYPE].numType, ppArray.types[NUMBERTYPE].name, "num.out");
    }

    if (paramHom)
    { // delete func_input_param
      remove("func_input_param");
    }

    // clear memory
    clear_variablegroupArray(&vargpArray);
    clear_subFuncData(ppArray.types[DEFINEDSUBFUNCTIONTYPE].numType, &sfData);
    clear_preProcessArray(&ppArray);
    if (*needToDiff == 1)
      clear_parseArray(&parArray);
  }

  if (*trackType == 3 && *sharpenOnly == 0)
  { // performing membership test 
    printf("\nNOTE: Please make sure that you have also put the proper test points for the membership test in 'member_points'.\n\n");
  }

  return;
}

void clearConfigurations(char **str, int *len, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the configurations                               *
\***************************************************************/
{
  int i;

  for (i = 0; i < size; i++)
    free(str[i]);
  free(str);
  free(len);

  return;
}

void loadParseConfigurations(char ***str, int **len, int *size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: load the configurations needed in the initial parsing  *
\***************************************************************/
{
  int i;

  *size = 13;
  *len = (int *)bmalloc(*size * sizeof(int));
  *str = (char **)bmalloc(*size * sizeof(char *));
  for (i = 0; i < *size; i++)
    (*str)[i] = (char *)bmalloc(25 * sizeof(char));

  strcpy((*str)[0], "TRACKTYPE:");
  strcpy((*str)[1], "MPTYPE:");
  strcpy((*str)[2], "PRECISION:");
  strcpy((*str)[3], "WITNESSGENTYPE:");
  strcpy((*str)[4], "USERHOMOTOPY:");
  strcpy((*str)[5], "RANDOMSEED:");
  strcpy((*str)[6], "SHARPENONLY:");
  strcpy((*str)[7], "NEEDTODIFF:");
  strcpy((*str)[8], "DELETETEMPFILES:");
  strcpy((*str)[9], "PARAMETERHOMOTOPY:");
  strcpy((*str)[10], "AMPMAXPREC:");
  strcpy((*str)[11], "USEREGENERATION:");
  strcpy((*str)[12], "REGENSTARTLEVEL:");

  for (i = 0; i < *size; i++)
    (*len)[i] = strlen((*str)[i]);

  return;
}

void loadConfigurations(char ***str, int **len, int *size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: load the configurations                                *
\***************************************************************/
{
  int i;

  *size = 76;
  *len = (int *)bmalloc(*size * sizeof(int));
  *str = (char **)bmalloc(*size * sizeof(char *));
  for (i = 0; i < *size; i++)
    (*str)[i] = (char *)bmalloc(25 * sizeof(char));

  strcpy((*str)[0], "TRACKTYPE:");
  strcpy((*str)[1], "MPTYPE:");
  strcpy((*str)[2], "PRECISION:");
  strcpy((*str)[3], "COEFFBOUND:");
  strcpy((*str)[4], "DEGREEBOUND:");
  strcpy((*str)[5], "AMPMAXPREC:");
  strcpy((*str)[6], "AMPSAFETYDIGITS1:");
  strcpy((*str)[7], "AMPSAFETYDIGITS2:");
  strcpy((*str)[8], "TRACKTOLBEFOREEG:");
  strcpy((*str)[9], "TRACKTOLDURINGEG:");
  strcpy((*str)[10], "ODEPREDICTOR:");
  strcpy((*str)[11], "MAXNEWTONITS:");
  strcpy((*str)[12], "MAXSTEPSIZE:");
  strcpy((*str)[13], "MINSTEPSIZEBEFOREEG:");
  strcpy((*str)[14], "MINSTEPSIZEDURINGEG:");
  strcpy((*str)[15], "MAXNUMBERSTEPS:");
  strcpy((*str)[16], "STEPSFORINCREASE:");
  strcpy((*str)[17], "STEPFAILFACTOR:");
  strcpy((*str)[18], "STEPSUCCESSFACTOR:");
  strcpy((*str)[19], "PATHTRUNCATIONTHRESHOLD:");
  strcpy((*str)[20], "FINALTOL:");
  strcpy((*str)[21], "ENDGAMENUM:");
  strcpy((*str)[22], "SAMPLEFACTOR:");
  strcpy((*str)[23], "NUMSAMPLEPOINTS:");
  strcpy((*str)[24], "ENDGAMEBDRY:");
  strcpy((*str)[25], "MINCYCLETRACKBACK:");
  strcpy((*str)[26], "NBHDRADIUS:");
  strcpy((*str)[27], "MAXCYCLENUM:");
  strcpy((*str)[28], "SECURITYLEVEL:");
  strcpy((*str)[29], "SECURITYMAXNORM:");
  strcpy((*str)[30], "TARGETTOLMULTIPLIER:");
  strcpy((*str)[31], "IMAGTHRESHOLD:");
  strcpy((*str)[32], "CONDNUMTHRESHOLD:");
  strcpy((*str)[33], "ENDPOINTFINITETHRESHOLD:");
  strcpy((*str)[34], "SHARPENDIGITS:");
  strcpy((*str)[35], "SHARPENONLY:");
  strcpy((*str)[36], "USERHOMOTOPY:");
  strcpy((*str)[37], "PARAMETERHOMOTOPY:");
  strcpy((*str)[38], "WITNESSGENTYPE:");
  strcpy((*str)[39], "MAXCODIMENSION:");
  strcpy((*str)[40], "JUNKREMOVALTEST:");
  strcpy((*str)[41], "MAXLDTDEPTH:");
  strcpy((*str)[42], "WITNESSSUPERSETONLY:");
  strcpy((*str)[43], "MULTONEONLY:");
  strcpy((*str)[44], "MAXNUMPTSFORTRACE:");
  strcpy((*str)[45], "MAXNUMMONOLOOPS:");
  strcpy((*str)[46], "MAXNUMBADMONOLOOPS:");
  strcpy((*str)[47], "USEREGENERATION:");
  strcpy((*str)[48], "SLICETOLBEFOREEG:");
  strcpy((*str)[49], "SLICETOLDURINGEG:");
  strcpy((*str)[50], "SLICEFINALTOL:");
  strcpy((*str)[51], "REGENSTARTLEVEL:");
  strcpy((*str)[52], "REGENREMOVEINF:");
  strcpy((*str)[53], "USEDIAGONAL:");
  strcpy((*str)[54], "OUTPUTLEVEL:");
  strcpy((*str)[55], "PRINTPATHPROGRESS:");
  strcpy((*str)[56], "SCREENOUT:");
  strcpy((*str)[57], "RANDOMSEED:");
  strcpy((*str)[58], "INTRINSICMULTIPLIER:");
  strcpy((*str)[59], "SINGVALZEROTOL:");
  strcpy((*str)[60], "DELETETEMPFILES:");
  strcpy((*str)[61], "NEEDTODIFF:");
  strcpy((*str)[62], "CYCLETIMECUTOFF:");
  strcpy((*str)[63], "RATIOTIMECUTOFF:");
  strcpy((*str)[64], "REGENHIGHERDIMTEST:");
  strcpy((*str)[65], "FUNCTIONTOLERANCE:");
  strcpy((*str)[66], "RATIOTOLERANCE:");
  strcpy((*str)[67], "MAXSTEPSBEFORENEWTON:");
  strcpy((*str)[68], "SPECIFICCODIMENSION:");
  strcpy((*str)[69], "CONSTRUCTWITNESSSET:");

  // deprecated
  strcpy((*str)[70], "MAXNORM:");
  strcpy((*str)[71], "MAXNUMMONLINEARS:");
  strcpy((*str)[72], "MAXNUMBADLOOPSINMON:");
  strcpy((*str)[73], "PRINTPATHMODULUS:");
  strcpy((*str)[74], "REDUCEDONLY:");
  strcpy((*str)[75], "TARGETTIME:");

  for (i = 0; i < *size; i++)
    (*len)[i] = strlen((*str)[i]);

  return;
}

void setupParseConfiguration(int num, FILE *IN, int *trackType, int *MPType, int *precision, int *genType, int *userHom, unsigned int *randomSeed, int *sharpenOnly, int *needToDiff, int *remove_temp, int *paramHom, int *AMPMaxPrec, int *useRegen, int *regenStartLevel)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the configuration described by num               *
*    DOES NOT PRINT ERROR MESSAGES                              *
\***************************************************************/
{
  int t;
  double d;

  // rules on how to setup configuration
  if (num == 0)
  { // TrackType
    fscanf(IN, "%d", &t);
    if (-4 <= t && t <= 7)
      *trackType = t;
  }
  else if (num == 1)
  { // MPType
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      *MPType = t;
  }
  else if (num == 2)
  { // Precision
    fscanf(IN, "%d", &t);
    if (64 <= t) // error checking on value
      *precision = t;
  }
  else if (num == 3)
  { // WitnessGenType
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      *genType = t;
  }
  else if (num == 4)
  { // UserHomotopy
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      *userHom = t;
  }
  else if (num == 5)
  { // RandomSeed
    fscanf(IN, "%lf", &d);
    if ((unsigned int) d >= 0)
      *randomSeed = (unsigned int) d;
  }
  else if (num == 6)
  { // SharpenOnly
   fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      *sharpenOnly = t;
  }
  else if (num == 7)
  { // NeedToDiff
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)
      *needToDiff = t;
  }
  else if (num == 8)
  { // DeleteTempFiles
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      *remove_temp = t;
  }
  else if (num == 9)
  { // ParameterHomotopy
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2) // error checking on value
      *paramHom = t;
  }
  else if (num == 10)
  { // AMPMaxPrec
    fscanf(IN, "%d", &t);
    if (t >= 64) // error checking on value
      *AMPMaxPrec = t;
  }
  else if (num == 11)
  { // UseRegeneration
   fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
      *useRegen = t;
  }
  else if (num == 12)
  { // RegenStartLevel
    fscanf(IN, "%d", &t);
    if (t >= 0) // error checking on value
      *regenStartLevel = t;
  }

  return;
}

void setupConfiguration(int num, FILE *IN, tracker_config_t *T, int *trackType, int *userHom, int *useRegen, int *regenStartLevel, int *maxCodim, int *specificCodim, int *printMod, double *intrinsicCutoffMultiplier, int *reducedOnly, int *constructWitnessSet, int *supersetOnly, int *paramHom, unsigned int *randomSeed, int *genType, int *remove_temp, int *needToDiff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the configuration described by num               *
\***************************************************************/
{
  int t;
  double d;

  // rules on how to setup configuration
  if (num == 0)
  { // TrackType
    fscanf(IN, "%d", &t); 
    if (-4 <= t && t <= 7)
      *trackType = t;
    else
      printf("WARNING: TrackType needs to be between -4 and 7. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 1)
  { // MPType
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      T->MPType = t;
    else
      printf("WARNING: MPType needs to be either 0, 1, or 2. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 2)
  { // Precision
    fscanf(IN, "%d", &t);
    if (64 <= t) // error checking on value
      T->Precision = t;
    else
      printf("WARNING: Precision needs to be at least 64. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 3)
  { // COeffBound
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->AMP_bound_on_abs_vals_of_coeffs = d;
    else
      printf("WARNING: CoeffBound needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 4)
  { // DegreeBound
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->AMP_bound_on_degree = d;
    else
      printf("WARNING: DegreeBound needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 5)
  { // AMPMaxPrec
    fscanf(IN, "%d", &t);
    if (t >= 64) // error checking on value
      T->AMP_max_prec = t;
    else
      printf("WARNING: AMPMaxPrec needs to be at least 64. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 6)
  { // AMPSafetyDigits1
    fscanf(IN, "%d", &t); // no error checking needed
    T->AMP_safety_digits_1 = t;
  }
  else if (num == 7)
  { // AMPSafetyDigits2
    fscanf(IN, "%d", &t); // no error checking needed
    T->AMP_safety_digits_2 = t;
  }
  else if (num == 8)
  { // TrackTolBeforeEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->basicNewtonTol = d;
    else
      printf("WARNING: TrackTolBeforeEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 9)
  { // TrackTolDuringEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->endgameNewtonTol = d;
    else
      printf("WARNING: TrackTolDuringEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 10)
  { // ODEPredictor
    fscanf(IN, "%d", &t);
    if (-1 <= t && t <= 8) // error check on value
      T->odePredictor = t;
    else
      printf("WARNING: ODEPredictor needs to be between -1 and 8. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 11)
  { // MaxNewtonIts
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->maxNewtonIts = t;
    else
      printf("WARNING: MaxNewtonIts needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 12)
  { // MaxStepSize
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->maxStepSize = d;
    else
      printf("WARNING: MaxStepSize needs to be strictly between 0 & 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 13)
  { // MinStepSizeBeforeEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->minStepSizeBeforeEndGame = d;
    else
      printf("WARNING: MinStepSizeBeforeEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 14)
  { // MinStepSizeDuringEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->minStepSizeDuringEndGame = d;
    else
      printf("WARNING: MinStepSizeDuringEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 15)
  { // MaxNumberSteps
    fscanf(IN, "%d", &t);
    if (t > 0)  // error checking on value
      T->maxNumSteps = t;
    else
      printf("WARNING: MaxNumberSteps needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 16)
  { // StepsForIncrease
    fscanf(IN, "%d", &t);
    if (t >= 1) // error checking on value
      T->cSecInc = t;
    else
      printf("WARNING: StepsForIncrease needs to be >= 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 17)
  { // StepFailFactor
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1) // error checking on value
      T->step_fail_factor = d;
    else
      printf("WARNING: StepFailFactor needs to be strictly between 0 and 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 18)
  { // Step SuccessFactor
    fscanf(IN, "%lf", &d);
    if (d >= 1) // error checking on value
      T->step_success_factor = d;
    else
      printf("WARNING: StepSuccessFactor needs to be >= 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 19)
  { // PathTruncationThreshold
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error check on value
      T->goingToInfinity = d;
    else
      printf("WARNING: PathTruncationThreshold needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 20)
  { // FinalTol
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->final_tolerance = d;
    else
      printf("WARNING: FinalTol needs to be > 0. The input value (%e) has been ignored.\n", d);
  } 
  else if (num == 21)
  { // EndgameNum
    fscanf(IN, "%d", &t);
    if (t == 1 || t == 2 || t == 3)  // error checking on value
      T->endgameNumber = t;
    else
      printf("WARNING: EndgameNum needs to be either 1, 2, or 3. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 22)
  { // SampleFactor
    fscanf(IN, "%lf", &d);
    if ((0 < d) && (d < 1)) // error checking on value
      T->power_series_sample_factor = d;
    else
      printf("WARNING: SampleFactor needs to be > 0 and < 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 23)
  { // NumSamplePoints
    fscanf(IN, "%d", &t);
    if (t >= 2) // error checking on value
      T->num_PSEG_sample_points = t;
    else
      printf("WARNING: NumSamplePoints needs to be at least 2. The input value (%d) has been ignored.\n", t);  
  }
  else if (num == 24)
  { // EndgameBdry
    fscanf(IN, "%lf", &d);
    if ((0 < d) && (d < 1))  // error checking on value
      T->endgameBoundary = d;
    else
      printf("WARNING: EndgameBdry needs to be strictly between 0 & 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 25)
  { // MinCycleTrackBack
    fscanf(IN, "%d", &t);
    if (t > 0)  // error checking on value
      T->minCycleTrackBack = t;
    else
      printf("WARNING: MinCycleTrackback needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 26)
  { // NbhdRadius
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->minTrackT = d;
    else
      printf("WARNING: NbhdRadius needs to be strictly between 0 & 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 27)
  { // MaxCycleNum
    fscanf(IN, "%d", &t);
    if (t > 0)  // error checking on value
      T->cycle_num_max = t;
    else
      printf("WARNING: MaxCycleNum needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 28)
  { // SecurityLevel
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      T->securityLevel = t;
    else
      printf("WARNING: SecurityLevel needs to be 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 29)
  { // SecurityMaxNorm
    fscanf(IN, "%lf", &d);
    if (d > 0) // error checking on value
      T->securityMaxNorm = d;
    else
      printf("WARNING: SecurityMaxNorm needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 30)
  { // TargetTolMultiplier
   fscanf(IN, "%lf", &d);
    if (d >= 1)  // error checking on value
      T->final_tol_multiplier = d;
    else
      printf("WARNING: TargetTolMultiplier needs to be >= 1. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 31)
  { // ImagThreshold
    fscanf(IN, "%lf", &d);
    if (d > 0) // error checking on value
      T->real_threshold = d;
    else
      printf("WARNING: ImagThreshold needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 32)
  { // CondNumThreshold
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->cond_num_threshold = d;
    else
      printf("WARNING: CondNumThreshold needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 33)
  { // EndpointFiniteThreshold
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error check on value
      T->finiteThreshold = d;
    else
      printf("WARNING: EndpointFiniteThreshold needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 34)
  { // SharpenDigits
    fscanf(IN, "%d", &t);
    if (t > 0) // error checking on value
      T->sharpenDigits = t;
    else
      printf("WARNING: SharpenDigits needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 35)
  { // SharpenOnly
   fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      T->sharpenOnly = t;
    else
      printf("WARNING: SharpenOnly needs to be 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 36)
  { // UserHomotopy
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      *userHom = t;
    else
      printf("WARNING: UserHomotopy needs to be either 0, 1, or 2. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 37)
  { // ParameterHomotopy
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2) // error checking on value
      *paramHom = t;
    else
      printf("WARNING: ParameterHomotopy needs to be between either 0, 1, or 2.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 38)
  { // WitnessGenType
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)  // error checking on value
      *genType = t;
    else
      printf("WARNING: WitnessGenType needs to be either 0, 1, or 2. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 39)
  { // MaxCodimension
    fscanf(IN, "%d", &t);
    if (t >= 0) // error checking on value
      *maxCodim = t;
    else
      printf("WARNING: MaxCodimension needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 40)
  { // JunkRemovalTest
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      T->junkRemovalTest = t;
    else
      printf("WARNING: JunkRemovalTest needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 41)
  { // MaxLDTDepth
    fscanf(IN, "%d", &t);
    if (t > 0) // error checking on value
      T->maxDepthLDT = t;
    else
      printf("WARNING: MaxLDTDepth needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 42)
  { // WitnessSupersetOnly
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
      *supersetOnly = t;
    else
      printf("WARNING: WitnessSupersetOnly needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 43)
  { // MultOneOnly
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      *reducedOnly = t;
    else
      printf("WARNING: MultOneOnly needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 44)
  { // MaxNumPtsForTrace
    fscanf(IN, "%d", &t);
    if (t > 0)  // error checking on value
      T->max_num_pts_for_trace = t;
    else
      printf("WARNING: MaxNumPtsForTrace needs to be > 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 45)
  { // MaxNumMonoLoops
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->max_num_mon_linears = t;
    else
      printf("WARNING: MaxNumMonoLoops needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 46)
  { // MaxNumBadMonoLoops
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->max_num_bad_loops_in_mon = t;
    else
      printf("WARNING: MaxNumBadMonoLoops needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 47)
  { // UseRegeneration
   fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
      *useRegen = t;
    else
      printf("WARNING: UseRegeneration needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 48)
  { // SliceTolBeforeEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->sliceBasicNewtonTol = d;
    else
      printf("WARNING: SliceTolBeforeEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 49)
  { // SliceTolDuringEG
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->sliceEndgameNewtonTol = d;
    else
      printf("WARNING: SliceTolDuringEG needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 50)
  { // SliceFinalTol
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->sliceFinalTol = d;
    else
      printf("WARNING: SliceFinalTol needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 51)
  { // RegenStartLevel
    fscanf(IN, "%d", &t);
    if (t >= 0) // error checking on value
      *regenStartLevel = t;
    else
      printf("WARNING: RegenStartLevel needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 52)
  { // RegenRemoveInf
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      T->regen_remove_inf = t;
    else
      printf("WARNING: RegenRemoveInf needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 53)
  { // UseDiagonal
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
    {
      if (t == 1)
        *userHom = -59;
    }
    else
      printf("WARNING: UseDiagonal needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 54)
  { // OutputLevel
    fscanf(IN, "%d", &t);
    if (-1 <= t && t <= 3)  // error checking on value
      T->outputLevel = t;
    else
      printf("WARNING: OutputLevel needs to be between -1 and 3.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 55)
  { // PrintPathProgress
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      *printMod = t;
    else
      printf("WARNING: PrintPathProgress needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 56)
  { // ScreenOut
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)
      T->screenOut = t;
    else
      printf("WARNING: ScreenOut needs to be either 0 or 1 . The input value (%d) has been ignored.\n", t);
  }
  else if (num == 57)
  { // RandomSeed
    fscanf(IN, "%lf", &d);
    if ((unsigned int) d >= 0)
      *randomSeed = (unsigned int) d;
  }
  else if (num == 58)
  { // IntrinsicMultiplier
    fscanf(IN, "%lf", &d);
    if (0 <= d &&  d <= 1) // error checking on value
      *intrinsicCutoffMultiplier = d;
    else
      printf("WARNING: IntrinsicMultiplier needs to be in [0,1]. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 59)
  { // SingValZeroTol
    fscanf(IN, "%lf", &d);
    if (d > 0) // error checking on value
      T->sing_val_zero_tol = d;
    else
      printf("WARNING: SingValZeroTol needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 60)
  { // DeleteTempFiles
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      *remove_temp = t;
    else
      printf("WARNING: DeleteTempFiles needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 61)
  { // NeedToDiff
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1 || t == 2)
      *needToDiff = t;
    else
      printf("WARNING: NeedToDiff needs to be either 0, 1, or 2. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 62)
  { // CycleTimeCutoff
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->cutoffCycleTime = d;
    else
      printf("WARNING: CycleTimeCutoff needs to be strictly between 0 & 1.  The input value (%e) has been ignored.\n", d);
  }
  else if (num == 63)
  { // RatioTimeCutoff
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->cutoffRatioTime = d;
    else
      printf("WARNING: RatioTimeCutoff needs to be strictly between 0 & 1.  The input value (%e) has been ignored.\n", d);
  }
  else if (num == 64)
  { // RegenHigherDimTest 
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
      T->regen_higher_dim_check = t;
    else
      printf("WARNING: RegenHigherDimTest needs to be either 0 or 1.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 65)
  { // FunctionTolerance
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->funcResTol = d;
    else
      printf("WARNING: FunctionTolearnce needs to be strictly between 0 & 1.  The input value (%e) has been ignored.\n", d);
  }
  else if (num == 66)
  { // RatioTolerance
    fscanf(IN, "%lf", &d);
    if (0 < d && d < 1)  // error checking on value
      T->ratioTol = d;
    else
      printf("WARNING: RatioTolearnce needs to be strictly between 0 & 1.  The input value (%e) has been ignored.\n", d);
  }
  else if (num == 67)
  { // MaxStepsBefo8eNewton
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->maxStepsBeforeNewton = t;
    else
      printf("WARNING: MaxStepsBeforeNewton needs to be >= 0.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 68)
  { // SpecificCodimension 
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      *specificCodim = t;
    else
      printf("WARNING: SpecificCodimension needs to be >= 0.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 69)
  { // ConstructWitnessSet
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1)  // error checking on value
      *constructWitnessSet = t;
    else
      printf("WARNING: ConstructWitnessSet needs to be either 0 or 1.  The input value (%d) has been ignored.\n", t);
  }
  else if (num == 70)
  { // MaxNorm
    printf("WARNING: MaxNorm is deprecated.  Please use PathTruncationThreshold and EndpointFiniteThreshold.\n");
    fscanf(IN, "%lf", &d);
    if (d > 0)  // error checking on value
      T->goingToInfinity = T->finiteThreshold = d;
    else
      printf("WARNING: MaxNorm needs to be > 0. The input value (%e) has been ignored.\n", d);
  }
  else if (num == 71)
  { // MaxNumMonLinears
    printf("WARNING: MaxNumMonLinears is deprecated.  Please use MaxNumMonoLoops.\n");
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->max_num_mon_linears = t;
    else
      printf("WARNING: MaxNumMonoLoops needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 72)
  { // MaxNumBadLoopsInMon 
    printf("WARNING: MaxNumBadLoopsInMon is deprecated.  Please use MaxNumBadMonoLoops.\n");
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      T->max_num_bad_loops_in_mon = t;
    else
      printf("WARNING: MaxNumBadMonoLoops needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 73)
  { // PrintPathModulus
    printf("WARNING: PrintPathModulus is deprecated.  Please use PrintPathProgress.\n");
    fscanf(IN, "%d", &t);
    if (t >= 0)  // error checking on value
      *printMod = t;
    else
      printf("WARNING: PrintPathProgress needs to be >= 0. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 74)
  { // ReducedOnly
    printf("WARNING: ReducedOnly is deprecated.  Please use MultOneOnly.\n");
    fscanf(IN, "%d", &t);
    if (t == 0 || t == 1) // error checking on value
      *reducedOnly = t;
    else
      printf("WARNING: MultOneOnly needs to be either 0 or 1. The input value (%d) has been ignored.\n", t);
  }
  else if (num == 75)
  { // TargetTime
    printf("WARNING: TargetTime is deprecated.\n");
    fscanf(IN, "%lf", &d);
    if (0 <= d && d < 1)  // error checking on value
      T->targetT = d;
    else
      printf("WARNING: TargetTime needs to be in [0,1). The input value (%e) has been ignored.\n", d);
  }

  return;
}

void parseConfigurations(FILE *IN, int *trackType, int *MPType, int *genType, int *userHom, unsigned int *randomSeed, int *sharpenOnly, int *needToDiff, int *remove_temp, int *paramHom, int *maxPrec, int *useRegen, int *regenStartLevel)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the configurations needed for initial parsing    *
\***************************************************************/
{
  int rV, size, *len = NULL, prec = 96, AMPMaxPrec = 1024, chSize = 0, loc = 0;
  char ch[25], **str = NULL;

  // load the configurations
  loadParseConfigurations(&str, &len, &size);

  // initialize the data
  *trackType = *userHom = *sharpenOnly = *paramHom = *useRegen = *regenStartLevel = 0;
  *randomSeed = 0;
  *needToDiff = *remove_temp = 1;
  *MPType = *genType = 2;

  // read in the first character of the line
  ch[0] = fgetc(IN);
  chSize = 1;
  while (ch[0] != EOF)
  { // see if we should ignore line
    if (ch[0] != '%')
    { // determine the configuration
      rV = findConfig(&loc, ch, &chSize, IN, str, len, size);
      if (rV)
      { // found a configuration to setup
        setupParseConfiguration(loc, IN, trackType, MPType, &prec, genType, userHom, randomSeed, sharpenOnly, needToDiff, remove_temp, paramHom, &AMPMaxPrec, useRegen, regenStartLevel);
      }

      // scan rest of line
      scanRestOfLine(IN);
    }
    else
    { // ignore test of line
      scanRestOfLine(IN);
    }

    // read in the next line
    ch[0] = fgetc(IN);
    chSize = 1;
  }

  // setup the maximum precision
  *maxPrec = MAX(prec, AMPMaxPrec);

  // clear configurations
  clearConfigurations(str, len, size);

  return;
}

void addBounds(double degBound, double coeffBound, int printDeg, char *fileName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add bound data to config                               *
\***************************************************************/
{
  int configExist = 0, size = 0;
  char ch, *newName = NULL;
  FILE *CONFIG = NULL, *OUT = NULL;

  // verify existence
  CONFIG = fopen(fileName, "r");
  if (CONFIG == NULL)
  { // does not exist
    configExist = 0;
  }
  else
  { // exists
    configExist = 1;
  }
  fclose(CONFIG);

  if (configExist)
  { // move to another location and recreate
    size = snprintf(NULL, 0, "%s_old", fileName);
    newName = (char *)bmalloc((size + 1) * sizeof(char));
    sprintf(newName, "%s_old", fileName);
    rename(fileName, newName);
  }

  // create a new file
  OUT = fopen(fileName, "w");

  // round the value up
  if (degBound - (int) degBound > 1e-15)
    degBound = (int) degBound + 1;
  // print configurations, if needed
  fprintf(OUT, "CoeffBound: %e;\n", coeffBound);
  if (printDeg)
    fprintf(OUT, "DegreeBound: %d;\n", (int) degBound);

  if (configExist)
  { // copy old configurations
    CONFIG = fopen(newName, "r");
    ch = fgetc(CONFIG);
    while (ch != EOF)
    {
      fprintf(OUT, "%c", ch);
      ch = fgetc(CONFIG);
    }

    // close and delete
    fclose(CONFIG);
    remove(newName);
  }

  // close OUT
  fclose(OUT);

  // clear memory
  free(newName);

  return;
}

int setupConfig(tracker_config_t *T, double *midpointTol, int *userHom, int *useRegen, int *regenStartLevel, int *maxCodim, int *specificCodim, int *printMod, double *intrinsicCutoffMultiplier, int *reducedOnly, int *constructWitnessSet, int *supersetOnly, int *paramHom, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *IN = NULL;
  int rV, size, *len = NULL, chSize = 0, loc = 0, trackType, genType, remove_temp, needToDiff;
  unsigned int randomSeed;
  char ch[25], **str = NULL;

  // load the configurations
  loadConfigurations(&str, &len, &size);

  // give a little spacing so that the error messages stand out
  printf("\n");

  IN = fopen("config", "r");
  if (IN == NULL)
  {
    printf("\nconfig does not exist!\n\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // load the default configurations
  loadDefaultConfig(T, MPType, &genType, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, printMod, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // scan to see what configurations the user wants changed
  ch[0] = fgetc(IN);
  chSize = 1;
  while (ch[0] != EOF)
  { // see if we should ignore line
    if (ch[0] != '%')
    { // determine the configuration
      rV = findConfig(&loc, ch, &chSize, IN, str, len, size);
      if (rV)
      { // found a configuration to setup
        setupConfiguration(loc, IN, T, &trackType, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, printMod, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom, &randomSeed, &genType, &remove_temp, &needToDiff);

        // scan rest of line
        scanRestOfLine(IN);
      }
      else 
      { // print error message
        errorConfig(ch, chSize, IN);
      }
    }
    else
    { // ignore test of line
      scanRestOfLine(IN);
    }

    // read in the next line
    ch[0] = fgetc(IN);
    chSize = 1;
  }
  // close the file
  fclose(IN);

  // do error checking
  configErrorChecking(T, maxCodim, specificCodim, trackType);

  if (*userHom == -59 && *useRegen == 1)
    *userHom = 0;

  // set rest of values
  *midpointTol = T->basicNewtonTol;
  T->final_tol_times_mult = T->final_tolerance * T->final_tol_multiplier;

  // give a little spacing after the warning messages
  printf("\n");

  // clear configurations
  clearConfigurations(str, len, size);

  return 0;
}

void configErrorChecking(tracker_config_t *T, int *maxCodim, int *specificCodim, int trackType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does error checking on T                               *
\***************************************************************/
{
  if ((T->endgameNumber == 1 || T->endgameNumber == 0) && T->targetT != 0)
  { // set target time to 0 when using PSEG
    printf("WARNING: TargetTime != 0 when using the either the zero order or power series endgames.  Bertini, by default, will make TargetTime = 0.\n");
    T->targetT = 0;
  }

  if (T->endgameBoundary < T->targetT)
  { // set endgame boundary to be the target time
    printf("WARNING: EndgameBdry < TargetTime.  Bertini, by default, will make them equal.\n");
    T->endgameBoundary = T->targetT;
  }

  if (T->minStepSizeBeforeEndGame > T->maxStepSize)
  { // minimum step size needs to be smaller than the maximum
    printf("WARNING: MinStepSizeBeforeEG > MaxStepSize.  Bertini, by default, will make them equal.\n");
    T->minStepSizeBeforeEndGame = T->maxStepSize;
  }
 
  if (T->minStepSizeDuringEndGame > T->maxStepSize)
  { // minimum step size needs to be smaller than the maximum
    printf("WARNING: MinStepSizeDuringEG > MaxStepSize.  Bertini, by default, will make them equal.\n");
    T->minStepSizeDuringEndGame = T->maxStepSize;
  }

  if (T->sharpenOnly == 0)
  { // regular path tracking
    if (T->sharpenDigits > 0)
    { // using sharpening

      int digits_tol = (int) floor(-log10(T->final_tolerance)), digits_prec;
      if (T->MPType == 0)
        digits_prec = prec_to_digits(52);
      else if (T->MPType == 1)
        digits_prec = prec_to_digits(T->Precision) - 1;
      else
        digits_prec = INT_MAX;

      // compare with sharpenDigits
      if (T->sharpenDigits <= digits_tol)
      { // the number of sharpening digits needs to be more than the number of digits requested by the final tolerance
        printf("WARNING: SharpenDigits <= the number of digits required based on FinalTol for path tracking. Bertini, by default,\n");
        printf("         will ignore SharpenDigits.\n");
        T->sharpenDigits = 0;
      }
      else if (T->sharpenDigits > digits_prec)
      { // the number of sharpening digits needs to be smaller than the digits of precision
        printf("WARNING: SharpenDigits cannot be reliably calculated using the current precision. Bertini, by default, will\n");
        printf("         set SharpenDigits to %d.\n", digits_prec - 1);
        T->sharpenDigits = digits_prec - 1;
      }
    }

    if (T->MPType == 2)
    { // display message when using AMP
      printf("\nNOTE: You have requested to use adaptive path tracking.  Please make sure that you have\nsetup the following tolerances appropriately:\n");
      printf("CoeffBound: %.12e, DegreeBound: %.12e\nAMPSafetyDigits1: %d, AMPSafetyDigits2: %d, AMPMaxPrec: %d\n", T->AMP_bound_on_abs_vals_of_coeffs, T->AMP_bound_on_degree, T->AMP_safety_digits_1, T->AMP_safety_digits_2, T->AMP_max_prec);
    }
    else if (T->endgameNumber == 2)
    { // display message when using Cauchy EG
      printf("\nNOTE: You have requested to use the Cauchy Endgame. Please make sure that you have\nsetup the following tolerances appropriately:\n");
      printf("CoeffBound: %.12e, DegreeBound: %.12e\n", T->AMP_bound_on_abs_vals_of_coeffs, T->AMP_bound_on_degree);
    }

    // setup funcResTol
    if (T->MPType == 0)
    { // setup to 14 digits if still default value
      if (T->funcResTol == 1e-300)
        T->funcResTol = 1e-14;
    }
    else if (T->MPType == 1)
    { // setup based on precision if still default value
      if (T->funcResTol == 1e-300)
      {
        int digits = prec_to_digits(T->Precision) - 2;
        if (digits < 300)
          T->funcResTol = pow(10, -digits);
      }
    }
  }
  else
  { // sharpening module
    if (T->sharpenDigits == 0)
    { // by default, will initialize sharpenDigits to 14
      T->sharpenDigits = 14;
    }
    
    // make sure sharpen digits is acceptable
    int digits_prec;
    if (T->MPType == 0)
      digits_prec = prec_to_digits(52);
    else if (T->MPType == 1)
      digits_prec = prec_to_digits(T->Precision) - 1;
    else 
      digits_prec = INT_MAX;

    // compare with sharpenDigits
    if (T->sharpenDigits > digits_prec)
    { // the number of sharpening digits needs to be smaller than the digits of precision
      printf("\nWARNING: SharpenDigits cannot be reliably calculated using the current precision. Bertini, by default, will\n");
      printf("         set SharpenDigits to %d.\n", digits_prec - 1);
      T->sharpenDigits = digits_prec - 1;
    }
  }

  if (trackType == 1 || trackType == 7)
  { // setup maximum codimension equal to the specific codimension, if needed
    if (*specificCodim > 0)
    {
      *maxCodim = *specificCodim;
      if (T->junkRemovalTest == 0)
      {
        printf("\nNOTE: When computing a specific codimension, the local dimension test will automatically be used!\n");
        T->junkRemovalTest = 1;
      }
    }
  }

  return;
}

int findConfig(int *found_loc, char *ch, int *count, FILE *IN, char **str, int *len, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether it is found or not                     *
* NOTES: Read in configuration name                             *
\***************************************************************/
{
  int i, t, error = 0, found = 0;

  // loop over the possible configurations
  for (i = 0; i < size; i++)
  { // initialize
    error = 0;

    // compare against already there
    for (t = 0; t < *count && t < len[i]; t++)
      if ((t < len[i] - 1 && ch[t] != str[i][t] && ch[t] != str[i][t] + 32) || (t == len[i] - 1 && ch[t] != str[i][t])) // str uses uppercase, + 32 for lower case, end has ':'
        error = 1;

    // read in more if needed
    for (*count = *count; *count < len[i] && !error; (*count)++)
    { // read in next character
      ch[*count] = fgetc(IN);
      if ((*count < len[i] - 1 && ch[*count] != str[i][*count] && ch[*count] != str[i][*count] + 32) || (*count == len[i] - 1 && ch[*count] != str[i][*count]))
        error = 1;
    }

    // see if it agrees
    if (error == 0)
    { // we have a match!
      *found_loc = i;
      found = 1;
    }
  }

  return found;
}

void errorConfig(char *ch, int count, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the error line                                   *
\***************************************************************/
{
  int i;

  printf("WARNING: the following line in the configurations was ignored\n");

  // print the characters already read in
  for (i = 0; i < count; i++)
    printf("%c", ch[i]);

  // see if we need to print more
  if (ch[count - 1] != '\n')
  { // print the rest of the line
    printRestOfLine(stdout, IN);
  }

  return;
}

void printConfigLine(FILE *OUT, int voidType, void const *value, int *title, char *str)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the configuration line to OUT-print title if nec *
\***************************************************************/
{
  if (!(*title))
  {
    fprintf(OUT, "CONFIG\n\n");
    *title = 1;
  }

  if (voidType == -199)
    fprintf(OUT, "%s: %.12e;\n", str, *(double *)value);
  else if (voidType == -11)
    fprintf(OUT, "%s: %u;\n", str, *(unsigned int *)value);
  else
    fprintf(OUT, "%s: %d;\n", str, *(int *)value);

  return;
}

void printConfigValues(FILE *OUT, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print values to OUT that are different than the default*
\***************************************************************/
{
  int title_printed = 0; // determine if title has been printed
  int genType_d, userHom_d, useRegen_d, regenStartLevel_d, maxCodim_d, specificCodim_d, pathMod_d, reducedOnly_d, constructWitnessSet_d, supersetOnly_d, paramHom_d; // hold default values
  double intrinsicCutoffMultiplier_d; // hold default values
  tracker_config_t T_default;

  // load default values to see if they are changed in T
  loadDefaultConfig(&T_default, T->MPType, &genType_d, &userHom_d, &useRegen_d, &regenStartLevel_d, &maxCodim_d, &specificCodim_d, &pathMod_d, &intrinsicCutoffMultiplier_d, &reducedOnly_d, &constructWitnessSet_d, &supersetOnly_d, &paramHom_d);

  // print out configuration lines if the values are not the default values
  if (trackType != 0)
    printConfigLine(OUT, 0, &trackType, &title_printed, "TrackType");
  if (trackType == 1 && genType != genType_d)
    printConfigLine(OUT, 0, &genType, &title_printed, "WitnessGenType");
  if (trackType != 0 && reducedOnly_d != reducedOnly)
    printConfigLine(OUT, 0, &reducedOnly, &title_printed, "MultOneOnly");
  if (trackType != 0 && constructWitnessSet_d != constructWitnessSet)
    printConfigLine(OUT, 0, &constructWitnessSet, &title_printed, "ConstructWitnessSet");
  if (trackType != 0 && supersetOnly_d != supersetOnly)
    printConfigLine(OUT, 0, &supersetOnly, &title_printed, "WitnessSupersetOnly");
  if (pathMod != pathMod_d)
    printConfigLine(OUT, 0, &pathMod, &title_printed, "PrintPathProgress");
  if (userHom == -59)
  { // use diagonal
    userHom = 1;
    printConfigLine(OUT, 0, &userHom, &title_printed, "UseDiagonal");
  }
  else if (userHom != userHom_d)
    printConfigLine(OUT, 0, &userHom, &title_printed, "UserHomotopy");
  if (paramHom != paramHom_d)
    printConfigLine(OUT, 0, &paramHom, &title_printed, "ParameterHomotopy");
  if (intrinsicCutoffMultiplier != intrinsicCutoffMultiplier_d)
    printConfigLine(OUT, -199, &intrinsicCutoffMultiplier, &title_printed, "IntrinsicMultiplier");
  if (useRegen != useRegen_d)
    printConfigLine(OUT, 0, &useRegen, &title_printed, "UseRegeneration");
  if (regenStartLevel != regenStartLevel_d)
    printConfigLine(OUT, 0, &regenStartLevel, &title_printed, "RegenStartLevel");
  if (maxCodim != maxCodim_d)
    printConfigLine(OUT, 0, &maxCodim, &title_printed, "MaxCodimension");
  if (specificCodim != specificCodim_d)
    printConfigLine(OUT, 0, &specificCodim, &title_printed, "SpecificCodimension");

  if (T->MPType != T_default.MPType)
    printConfigLine(OUT, 0, &T->MPType, &title_printed, "MPType");
  if ((T->MPType == 1) && (T->Precision != T_default.Precision)) // only need to print precision if we are dealing with fixed multiprecision
    printConfigLine(OUT, 0, &T->Precision, &title_printed, "Precision");
  if (T->screenOut != T_default.screenOut)
    printConfigLine(OUT, 0, &T->screenOut, &title_printed, "ScreenOut");
  if (T->outputLevel != T_default.outputLevel)
    printConfigLine(OUT, 0, &T->outputLevel, &title_printed, "OutputLevel");
  if (T->cSecInc != T_default.cSecInc)
    printConfigLine(OUT, 0, &T->cSecInc, &title_printed, "StepsForIncrease");
  if (T->maxNewtonIts != T_default.maxNewtonIts)
    printConfigLine(OUT, 0, &T->maxNewtonIts, &title_printed, "MaxNewtonIts");
  if (T->maxStepSize != T_default.maxStepSize)
    printConfigLine(OUT, -199, &T->maxStepSize, &title_printed, "MaxStepSize");
  if (T->minStepSizeBeforeEndGame != T_default.minStepSizeBeforeEndGame)
    printConfigLine(OUT, -199, &T->minStepSizeBeforeEndGame, &title_printed, "MinStepSizeBeforeEG");
  if (T->minStepSizeDuringEndGame != T_default.minStepSizeDuringEndGame)
    printConfigLine(OUT, -199, &T->minStepSizeDuringEndGame, &title_printed, "MinStepSizeDuringEG");
  if (T->minTrackT != T_default.minTrackT)
    printConfigLine(OUT, -199, &T->minTrackT, &title_printed, "NbhdRadius");
  if (T->basicNewtonTol != T_default.basicNewtonTol)
    printConfigLine(OUT, -199, &T->basicNewtonTol, &title_printed, "TrackTolBeforeEG");
  if (T->endgameNewtonTol != T_default.endgameNewtonTol)
    printConfigLine(OUT, -199, &T->endgameNewtonTol, &title_printed, "TrackTolDuringEG");
  if (T->final_tolerance != T_default.final_tolerance)
    printConfigLine(OUT, -199, &T->final_tolerance, &title_printed, "FinalTol");
  if (T->goingToInfinity != T_default.goingToInfinity)
    printConfigLine(OUT, -199, &T->goingToInfinity, &title_printed, "MaxNorm");
  if (T->endgameBoundary != T_default.endgameBoundary)
    printConfigLine(OUT, -199, &T->endgameBoundary, &title_printed, "EndgameBdry");
  if (T->targetT != T_default.targetT)
    printConfigLine(OUT, -199, &T->targetT, &title_printed, "TargetTime");
  if (T->maxNumSteps != T_default.maxNumSteps)
    printConfigLine(OUT, 0, &T->maxNumSteps, &title_printed, "MaxNumberSteps");
  if (T->endgameNumber != T_default.endgameNumber)
    printConfigLine(OUT, 0, &T->endgameNumber, &title_printed, "EndgameNum");
  if (T->power_series_sample_factor != T_default.power_series_sample_factor)
    printConfigLine(OUT, -199, &T->power_series_sample_factor, &title_printed, "SampleFactor");
  if (T->cycle_num_max != T_default.cycle_num_max)
    printConfigLine(OUT, 0, &T->cycle_num_max, &title_printed, "MaxCycleNum");
  if (T->num_PSEG_sample_points != T_default.num_PSEG_sample_points)
    printConfigLine(OUT, 0, &T->num_PSEG_sample_points, &title_printed, "NumSamplePoints");
  if (T->real_threshold != T_default.real_threshold)
    printConfigLine(OUT, -199, &T->real_threshold, &title_printed, "ImagThreshold");
  if (T->AMP_bound_on_abs_vals_of_coeffs != T_default.AMP_bound_on_abs_vals_of_coeffs)
    printConfigLine(OUT, -199, &T->AMP_bound_on_abs_vals_of_coeffs, &title_printed, "CoeffBound");
  if (T->AMP_bound_on_degree != T_default.AMP_bound_on_degree)
    printConfigLine(OUT, -199, &T->AMP_bound_on_degree, &title_printed, "DegreeBound");
  if (T->final_tol_multiplier != T_default.final_tol_multiplier)
    printConfigLine(OUT, -199, &T->final_tol_multiplier, &title_printed, "TargetTolMultiplier");
  if (T->AMP_safety_digits_1 != T_default.AMP_safety_digits_1)
    printConfigLine(OUT, 0, &T->AMP_safety_digits_1, &title_printed, "AMPSafetyDigits1");
  if (T->AMP_safety_digits_2 != T_default.AMP_safety_digits_2)
    printConfigLine(OUT, 0, &T->AMP_safety_digits_2, &title_printed, "AMPSafetyDigits2");
  if (T->AMP_max_prec != T_default.AMP_max_prec)
    printConfigLine(OUT, 0, &T->AMP_max_prec, &title_printed, "AMPMaxPrec");
  if (T->cond_num_threshold != T_default.cond_num_threshold)
    printConfigLine(OUT, -199, &T->cond_num_threshold, &title_printed, "CondNumThreshold");
  if (T->step_fail_factor != T_default.step_fail_factor)
    printConfigLine(OUT, -199, &T->step_fail_factor, &title_printed, "StepFailFactor");
  if (T->step_success_factor != T_default.step_success_factor)
    printConfigLine(OUT, 0, &T->step_success_factor, &title_printed, "StepSuccessFactor");
  if (T->max_num_pts_for_trace != T_default.max_num_pts_for_trace)
    printConfigLine(OUT, 0, &T->max_num_pts_for_trace, &title_printed, "MaxNumPtsForTrace");
  if (T->max_num_mon_linears != T_default.max_num_mon_linears)
    printConfigLine(OUT, 0, &T->max_num_mon_linears, &title_printed, "MaxNumMonoLoops");
  if (T->max_num_bad_loops_in_mon != T_default.max_num_bad_loops_in_mon)
    printConfigLine(OUT, 0, &T->max_num_bad_loops_in_mon, &title_printed, "MaxNumBadMonoLoops");
  if (T->sing_val_zero_tol != T_default.sing_val_zero_tol)
    printConfigLine(OUT, -199, &T->sing_val_zero_tol, &title_printed, "SingValZeroTol");
  if (T->sharpenDigits != T_default.sharpenDigits)
    printConfigLine(OUT, 0, &T->sharpenDigits, &title_printed, "SharpenDigits");
  if (T->sharpenOnly != T_default.sharpenOnly)
    printConfigLine(OUT, 0, &T->sharpenOnly, &title_printed, "SharpenOnly");
  if (T->regen_remove_inf != T_default.regen_remove_inf)
    printConfigLine(OUT, 0, &T->regen_remove_inf, &title_printed, "RegenRemoveInf");
  if (T->sliceBasicNewtonTol != T_default.sliceBasicNewtonTol)
    printConfigLine(OUT, -199, &T->sliceBasicNewtonTol, &title_printed, "SliceTolBeforeEG");
  if (T->sliceEndgameNewtonTol != T_default.sliceEndgameNewtonTol)
    printConfigLine(OUT, -199, &T->sliceEndgameNewtonTol, &title_printed, "SliceTolDuringEG");
  if (T->sliceFinalTol != T_default.sliceFinalTol)
    printConfigLine(OUT, -199, &T->sliceFinalTol, &title_printed, "SliceFinalTol");
  if (T->minCycleTrackBack != T_default.minCycleTrackBack)
    printConfigLine(OUT, 0, &T->minCycleTrackBack, &title_printed, "MinCycleTrackback");
  if (T->junkRemovalTest != T_default.junkRemovalTest)
    printConfigLine(OUT, 0, &T->junkRemovalTest, &title_printed, "JunkRemovalTest");
  if (T->maxDepthLDT != T_default.maxDepthLDT)
    printConfigLine(OUT, 0, &T->maxDepthLDT, &title_printed, "MaxLDTDepth");
  if (T->odePredictor != T_default.odePredictor)
    printConfigLine(OUT, 0, &T->odePredictor, &title_printed, "ODEPredictor");
  if (T->securityLevel != T_default.securityLevel)
    printConfigLine(OUT, 0, &T->securityLevel, &title_printed, "SecurityLevel");
  if (T->securityMaxNorm != T_default.securityMaxNorm)
    printConfigLine(OUT, -199, &T->securityMaxNorm, &title_printed, "SecurityMaxNorm");
  if (T->cutoffCycleTime != T_default.cutoffCycleTime)
    printConfigLine(OUT, -199, &T->cutoffCycleTime, &title_printed, "CycleTimeCutoff");
  if (T->cutoffRatioTime != T_default.cutoffRatioTime)
    printConfigLine(OUT, -199, &T->cutoffRatioTime, &title_printed, "RatioTimeCutoff");
  if (T->finiteThreshold != T_default.finiteThreshold)
    printConfigLine(OUT, -199, &T->finiteThreshold, &title_printed, "EndpointFiniteThreshold");
  if (T->regen_higher_dim_check != T_default.regen_higher_dim_check)
    printConfigLine(OUT, 0, &T->regen_higher_dim_check, &title_printed, "RegenHigherDimTest");
  if (T->funcResTol != T_default.funcResTol)
    printConfigLine(OUT, -199, &T->funcResTol, &title_printed, "FunctionTolerance");
  if (T->ratioTol != T_default.ratioTol)
    printConfigLine(OUT, -199, &T->ratioTol, &title_printed, "RatioTolerance");
  if (T->maxStepsBeforeNewton != T_default.maxStepsBeforeNewton)
    printConfigLine(OUT, 0, &T->maxStepsBeforeNewton, &title_printed, "MaxStepsBeforeNewton");

  printConfigLine(OUT, -11, &randomSeed, &title_printed, "RandomSeed");

  if (title_printed)
    fprintf(OUT, "\nEND;\n");

  return;
}

void loadDefaultConfig(tracker_config_t *T, int MPType, int *genType, int *userHom, int *useRegen, int *regenStartLevel, int *maxCodim, int *specificCodim, int *printMod, double *intrinsicMult, int *reducedOnly, int *constructWitnessSet, int *supersetOnly, int *paramHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: load the default values                                *
\***************************************************************/
{
  *genType = 2;
  *userHom = 0;
  *useRegen = 0;
  *regenStartLevel = 0;
  *maxCodim = 0;
  *specificCodim = 0;
  *printMod = 20;
  *intrinsicMult = 0.75;
  *reducedOnly = 0;
  *constructWitnessSet = 0;
  *supersetOnly = 0;
  *paramHom = 0;

  T->MPType = 2;
  T->Precision = 96;
  T->screenOut = 0;
  T->outputLevel = 0;
  T->cSecInc = 5;
  T->maxNewtonIts = 2;

  // step size settings
  T->maxStepSize = 0.1;
  T->step_fail_factor = 0.5;
  T->step_success_factor = 2.0;
  if (MPType == 2)
  { // load default values for AMP

    // minimum step size is controlled in the endgames - so we need to be not restrictive, but fail-safe
    T->minStepSizeBeforeEndGame = 1e-250;
    T->minStepSizeDuringEndGame = 1e-300;
  }
  else
  {
    T->minStepSizeBeforeEndGame = 1e-14;
    T->minStepSizeDuringEndGame = 1e-15;
  }

  // beyond 1e-150, computations in the endgames using double precision break down
  T->minTrackT = 1e-100;

  // the tolerances
  T->basicNewtonTol = 1e-5;
  T->endgameNewtonTol = 1e-6;
  T->final_tolerance = 1e-11;

  // post processing
  T->final_tol_multiplier = 10.0;

  T->real_threshold = 1e-8;

  T->goingToInfinity = 1e5;
  T->cond_num_threshold = 1e8;

  T->endgameBoundary = 0.1;
  T->targetT = 0.0;
  T->maxNumSteps = 10000;
  T->endgameNumber = 1;

  // power series settings
  T->power_series_sample_factor = 0.25;
  T->cycle_num_max = 6;
  T->num_PSEG_sample_points = 2;

  // AMP settings
  T->AMP_bound_on_abs_vals_of_coeffs = 1e3;
  T->AMP_bound_on_degree = 5;
  T->AMP_safety_digits_1 = 1;
  T->AMP_safety_digits_2 = 1;
  T->AMP_max_prec = 1024;

  T->sing_val_zero_tol = 1e-12;

  T->max_num_pts_for_trace = 10;
  T->max_num_mon_linears = 1000000;
  T->max_num_bad_loops_in_mon = 10;

  T->sharpenDigits = 0;
  T->sharpenOnly = 0;

  T->regen_remove_inf = 1;
  T->regen_higher_dim_check = 1;
  T->sliceBasicNewtonTol = 1e-7;
  T->sliceEndgameNewtonTol = 1e-8;
  T->sliceFinalTol = 1e-11;

  T->minCycleTrackBack = 4;
  T->junkRemovalTest = 1;
  T->maxDepthLDT = 3;
  T->odePredictor = 5;

  T->securityLevel = 0;
  T->securityMaxNorm = 1e4;

  T->cutoffCycleTime = 1e-8;
  T->cutoffRatioTime = 1e-14;
  T->finiteThreshold = 1e5;

  T->maxStepsBeforeNewton = 1;

  T->funcResTol = 1e-300; 
  T->ratioTol = 0.99;

  return;
}

void findNames(char ***names, int *numNames, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the random values                                *
\***************************************************************/
{
  int currSize = 0, currNames = 0;
  char ch;

  do 
  { // add a new name
    currSize = 0;
    *names = (char **)brealloc(*names, (currNames + 1) * sizeof(char *));
    (*names)[currNames] = (char *)bmalloc((currSize + 1) * sizeof(char));
    (*names)[currNames][currSize] = '\0';

    ch = fgetc(IN);
    while (ch != ',' && ch != ';')
    { // add ch
      (*names)[currNames] = (char *)brealloc((*names)[currNames], (currSize + 2) * sizeof(char));
      (*names)[currNames][currSize] = ch;
      currSize++;
      (*names)[currNames][currSize] = '\0';

      // read in next ch
      ch = fgetc(IN);
    }
    currNames++;

  } while (ch != ';'); 

  // setup numNames
  *numNames = currNames;

  return;
}

int setupRandomValues(FILE *OUT, FILE *IN, int createNewValues, int maxPrec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of random values setup                  *
* NOTES: setup the random values                                *
\***************************************************************/
{
  int i, j, size = 2, count = 0, numRand = 0, numNames = 0, isFileOpen = 0, rV = 0;
  int *len = (int *)bmalloc(size * sizeof(int));
  char ch[22] = "";
  char **str = (char **)bmalloc(size * sizeof(char *)), **names = NULL;
  FILE *RAND = NULL;
  comp_d rand_d;
  comp_mp rand_mp;
  mpq_t rand_rat[2];

  for (i = 0; i < size; i++)
    str[i] = (char *)bmalloc(22 * sizeof(char));

  strcpy(str[0], "random ");
  len[0] = strlen(str[0]);
  strcpy(str[1], "random_real ");
  len[1] = strlen(str[1]);

  // initialize MP
  initMP(maxPrec);
  init_mp(rand_mp);
  init_rat(rand_rat);

  if (createNewValues)
  { // need to generate the random numbers - only open file when needed
    while ((ch[0] = getc(IN)) != EOF)
    { // determine if 'r'
      if (ch[0] == 'r')
      { // setup count and see if we have a match
        count = 1;
        if (findConfig(&i, ch, &count, IN, str, len, size))
        { // find the names
          findNames(&names, &numNames, IN);
          scanRestOfLine(IN);
          numRand += numNames;

          // print the constant list
          fprintf(OUT, "constant ");
          for (j = 0; j < numNames; j++)
          { // print the names
            if (j + 1 < numNames)
              fprintf(OUT, "%s,", names[j]);
            else
              fprintf(OUT, "%s;\n", names[j]);
          }

          // print the values
          for (j = 0; j < numNames; j++)
          { // print the name
            fprintf(OUT, "%s=", names[j]);
            // generate a random number
            get_comp_rand_rat(rand_d, rand_mp, rand_rat, maxPrec, maxPrec, 0, 0);

            // print to OUT
            if (i == 0)
            { // print
              mpf_out_str(OUT, 10, 0, rand_mp->r);
              fprintf(OUT, "+I*");
              mpf_out_str(OUT, 10, 0, rand_mp->i);
              fprintf(OUT, ";\n");
            }
            else
            { // make real
              mpf_set_ui(rand_mp->i, 0);    
              mpf_out_str(OUT, 10, 0, rand_mp->r);
              fprintf(OUT, ";\n");
            }

            if (!isFileOpen)
            { // open file
              RAND = fopen("random_values", "w");
              fprintf(RAND, "                            \n\n");
              isFileOpen = 1;
            }

            // print to RAND
            print_mp(RAND, 0, rand_mp);
            fprintf(RAND, "\n");
          }

          // clear names
          for (j = 0; j < numNames; j++)
            free(names[j]);
          free(names);
          names = NULL;
        }
        else
        { // print line to OUT
          for (j = 0; j < count; j++)
            fprintf(OUT, "%c", ch[j]);
          printRestOfLine(OUT, IN);
        }
      }
      else
      { // print line to OUT
        fprintf(OUT, "%c", ch[0]);
        printRestOfLine(OUT, IN);
      }
    }

    // print the number of random values
    if (isFileOpen)
    {
      rewind(RAND);
      fprintf(RAND, "%d", numRand);
      fclose(RAND);
    }

    rV = numRand;
  }
  else
  { // look to read in the random numbers, they should exist - only open file when needed
    while ((ch[0] = getc(IN)) != EOF)
    { // determine if 'r'
      if (ch[0] == 'r')
      { // setup count and see if we have a match
        count = 1;
        if (findConfig(&i, ch, &count, IN, str, len, size))
        { // see if file is open
          if (!isFileOpen)
          { // open the file
            RAND = fopen("random_values", "r");
            if (RAND == NULL)
            { // error
              printf("\nERROR: 'random_values' does not exist!\n");
              bexit(ERROR_FILE_NOT_EXIST);
            }

            // read in the total number of random values
            fscanf(RAND, "%d", &numRand);
            rV = numRand;

            isFileOpen = 1;
          }

          // find the names
          findNames(&names, &numNames, IN);
          scanRestOfLine(IN);
          numRand -= numNames;

          if (numRand < 0)
          { // there are not enough random values!
            printf("\nERROR: There are not enough random values in 'random_values'!\n");
            bexit(ERROR_INVALID_SIZE);
          }

          // print the constant list
          fprintf(OUT, "constant ");
          for (j = 0; j < numNames; j++)
          { // print the names
            if (j + 1 < numNames)
              fprintf(OUT, "%s,", names[j]);
            else
              fprintf(OUT, "%s;\n", names[j]);
          }

          // print the values
          for (j = 0; j < numNames; j++)
          { // print the name
            fprintf(OUT, "%s=", names[j]);
            // read in a random number
            mpf_inp_str(rand_mp->r, RAND, 10);
            mpf_inp_str(rand_mp->i, RAND, 10);
            scanRestOfLine(RAND);

            // print to OUT
            if (i == 0)
            { // print
              mpf_out_str(OUT, 10, 0, rand_mp->r);
              fprintf(OUT, "+I*");
              mpf_out_str(OUT, 10, 0, rand_mp->i);
              fprintf(OUT, ";\n");
            }
            else
            { // make real
              mpf_set_ui(rand_mp->i, 0);
              mpf_out_str(OUT, 10, 0, rand_mp->r);
              fprintf(OUT, ";\n");
            }
          }

          // clear names
          for (j = 0; j < numNames; j++)
            free(names[j]);
          free(names);
          names = NULL;
        }
        else
        { // print line to OUT
          for (j = 0; j < count; j++)
            fprintf(OUT, "%c", ch[j]);
          printRestOfLine(OUT, IN);
        }
      }
      else
      { // print line to OUT
        fprintf(OUT, "%c", ch[0]);
        printRestOfLine(OUT, IN);
      }
    }

    if (isFileOpen)
    { // close RAND
      fclose(RAND);
    }
  }

  // clear memory
  clear_mp(rand_mp);
  clear_rat(rand_rat);
  for (i = size - 1; i >= 0; i--)
    free(str[i]);
  free(str);
  free(len);

  return numRand;
}

